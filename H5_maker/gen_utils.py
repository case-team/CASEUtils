
import ROOT
from ROOT import TLorentzVector, TFile
import numpy as np
import h5py

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import *
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.postprocessing.tools import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.JetSysColl import JetSysColl, JetSysObj
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import eventLoop
from PhysicsTools.NanoAODTools.postprocessing.framework.preskimming import preSkim

def isFinal(statusFlag):
    mask = 1 << 13
    return (statusFlag & mask) != 0

def isFirstCopy(statusFlag):
    mask = 1 << 12
    return (statusFlag & mask) != 0

Z_ID = 23
W_ID = 24
top_ID = 6
MAX_LEP_ID = 19

def prune_genparts(qs, n):
    if(len(qs) <= n): return qs

    max_dist = -9999
    prune_idxs = []
    need_to_prune = len(qs) - n
    for i, (dist, q) in enumerate(qs):
        if(dist > max_dist):
            max_dist = dist
            prune_idxs = [i]
        if(dist == max_dist):
            prune_idxs.append(i)

    while(len(prune_idxs) > need_to_prune):
        max_pt = -9999
        max_pt_idx = None
        for j,idx in enumerate(prune_idxs):
            if(qs[idx][1].pt > max_pt):
                max_pt = qs[idx][1].pt
                max_pt_idx = j
        prune_idxs.pop(max_pt_idx)
    qs = [q for i,q in enumerate(qs) if i not in prune_idxs]
    return prune_genparts(qs, n)






def findMother(coll, part, mother_ids, dist = 0):
    if(part.genPartIdxMother < 0): return None, -1
    mother =  coll[part.genPartIdxMother]
    if(abs(mother.pdgId) in mother_ids and isFirstCopy(mother.statusFlags)): return mother, dist+1
    return findMother(coll, mother, mother_ids, dist+1)


def parse_ZpToTpTp(event):
    GenPartsColl = Collection(event, "GenPart")

    Tp1 = Tp2 = top1 = top2 = Z1 = Z2 = None
    q1s = []
    q2s = []


    for gen_part in GenPartsColl:
        #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)
        #if(isFinal(gen_part.statusFlags)): print(gen_part.pt, gen_part.pdgId, gen_part.genPartIdxMother)
        if(isFirstCopy(gen_part.statusFlags)):
            if(abs(gen_part.pdgId) == Z_ID):
                if(Z1 is None): Z1 = gen_part
                elif(Z2 is None): Z2 = gen_part
                else: print("Extra Z!")
                #print("Z", gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)

            if(abs(gen_part.pdgId) == top_ID):
                if(top1 is None): top1 = gen_part
                elif(top2 is None): top2 = gen_part
                else: print("Extra top!")
                #print("top", gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)

    #Might need to swap order sometimes ? 
    #Z1, Z2 = (Z2,Z1) if (top1.genPartIdxMother == Z2.genPartIdxMother and top2.genPartIdxMother == Z1.genPartIdxMother) else (Z1,Z2)

    #Tp1v1 = ROOT.Math.PtEtaPhiMVector(Z1.pt, Z1.eta, Z1.phi, Z1.mass) + ROOT.Math.PtEtaPhiMVector(top1.pt, top1.eta, top1.phi, top1.mass)
    #Tp2v1 = ROOT.Math.PtEtaPhiMVector(Z2.pt, Z2.eta, Z2.phi, Z2.mass) + ROOT.Math.PtEtaPhiMVector(top2.pt, top2.eta, top2.phi, top2.mass)

    #Zp = Tp1v1 + Tp2v1

    #print("Gen", Tp1v1.mass(), Tp2v1.mass(), Tp1v2.mass(), Tp2v2.mass(), Zp.mass())

    parent_ids = {top_ID, W_ID, Z_ID}

    #avoid low pt garbage from shower
    pt_cut = 3.

    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids
                and gen_part.pt > pt_cut):
            mother, dist = findMother(GenPartsColl, gen_part, {Z_ID, top_ID}, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(mother is top1 or mother is Z1): q1s.append((dist, gen_part))
            elif(mother is top2 or mother is Z2): q2s.append((dist, gen_part))

    #gen matching isn't perfect, do some attempt at cleanup here
    #if(len(q1s) != 5 or len(q2s) != 5):
    #    print("Issue in quark finding!")
    #    print(len(q1s), len(q2s))
    #    print(q1s)
    #    print(q2s)
    if(len(q1s) > 5): 
        q1s = prune_genparts(q1s, 5)
    if(len(q2s) > 5): 
        q2s = prune_genparts(q2s, 5)
    q1_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1_vecs) < 5): q1_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2_vecs) < 5): q2_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1_vecs + q2_vecs




def parse_Wkk(event):
    GenPartsColl = Collection(event, "GenPart")

    Radion = W1 = W2 = W_ISO = None
    WKK_ID = 9000024
    RADION_ID = 9000025
    qs_iso = []
    qs_radion = []


    for i, gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        if(isFirstCopy(gen_part.statusFlags)):
            if(abs(gen_part.pdgId) == W_ID and gen_part.genPartIdxMother >= 0):
                m = GenPartsColl[gen_part.genPartIdxMother]
                if(abs(m.pdgId) == WKK_ID or gen_part.genPartIdxMother==0 ): W_ISO = gen_part
                elif(abs(m.pdgId) == RADION_ID): 
                    if(W1 is None): W1 = gen_part
                    elif(W2 is None): W2 = gen_part
                    else: print("Extra W!")    

    #exit(1)
    parent_ids = {W_ID}


    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):
            mother, dist = findMother(GenPartsColl, gen_part, parent_ids, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(mother is W_ISO): qs_iso.append((dist, gen_part))
            elif(mother is W1 or mother is W2): qs_radion.append((dist, gen_part))

    #if(len(qs_iso) != 2 or len(qs_radion) != 4):
    #    print("Issue in quark finding!")
    #    print(qs_iso)
    #    print(qs_radion)
    #    exit(1)

    #gen matching isn't perfect, do some attempt at cleanup here
    if(len(qs_iso) > 2): qs_iso = prune_genparts(qs_iso, 2)
    if(len(qs_radion) > 4): qs_radion = prune_genparts(qs_radion, 4)

    qs_iso_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in qs_iso ]
    qs_radion_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in qs_radion ]

    #zero pad if we missed some quarks
    while(len(qs_iso_vecs) < 2): qs_iso_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(qs_radion_vecs) < 4): qs_radion_vecs.append([-1.0, 0.0, 0.0, 0])

    return qs_radion_vecs + qs_iso_vecs


#number of gen particles to save and parsing function
gen_dict = {
        'ZpToTpTp' : (5*2, parse_ZpToTpTp),
        'Wkk' : (4+2, parse_Wkk),
         }

