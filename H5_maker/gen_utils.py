
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

def fromHardProcess(statusFlag):
    mask = 1 << 8
    return (statusFlag & mask) != 0


Z_ID = 23
W_ID = 24
top_ID = 6
b_ID = 5
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
    if(part.genPartIdxMother == 0): return None, dist
    mother =  coll[part.genPartIdxMother]
    if(abs(mother.pdgId) in mother_ids and isFirstCopy(mother.statusFlags)): return mother, dist+1
    return findMother(coll, mother, mother_ids, dist+1)


def parse_Wp(event):
    GenPartsColl = Collection(event, "GenPart")

    Z = top = None
    #Wp -> Bp t, Bp -> b Z
    Wp_ID = 6000024
    Bp_ID = 6000007
    q1s = []
    q2s = []

    parent_ids = {Wp_ID, Bp_ID}


    for i,gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        #if(isFinal(gen_part.statusFlags)): print(gen_part.pt, gen_part.pdgId, gen_part.genPartIdxMother)
        if(isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0  and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):

            if(abs(gen_part.pdgId) == Z_ID):
                if(Z is None): Z = gen_part
                else: print("Extra Z!")
                #print("Z", gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)
            if(abs(gen_part.pdgId) == top_ID):
                if(top is None): top = gen_part
                else: print("Extra top!")
                #print("top", gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)

    parent_ids = {Bp_ID, top_ID, Z_ID}


    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 and fromHardProcess(gen_part.statusFlags)):
            if(abs(gen_part.pdgId) == b_ID and (abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == Bp_ID or abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == Wp_ID)):
                q1s.append((0, gen_part))
            else:
                #TODO WIP 
                daughter_cand, daughter_dist = findMother(GenPartsColl, gen_part, {Z_ID, top_ID}, dist=0)
                if(daughter_cand is not None):
                    WP_cand, WP_dist = findMother(GenPartsColl, daughter_cand, {Wp_ID}, dist=0)
                    if(abs(daughter_cand.pdgId) == Z_ID and WP_dist >= 0 and GenPartsColl[gen_part.genPartIdxMother].pdgId == Z_ID): q1s.append((WP_dist + daughter_dist, gen_part))
                    elif(abs(daughter_cand.pdgId) == top_ID and WP_dist >=0 and 
                            (abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == top_ID or abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == W_ID)): 
                        q2s.append((WP_dist + daughter_dist, gen_part))

    #gen matching isn't always perfect, do some attempt at cleanup afterwards
    if(len(q1s) != 3 or len(q2s) != 3):
        print("Issue in quark finding!")
        print(len(q1s), len(q2s))
        print(q1s)
        print(q2s)


    if(len(q1s) > 3): q1s = prune_genparts(q1s, 3)
    if(len(q2s) > 3): q2s = prune_genparts(q2s, 3)

    q1_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1_vecs) < 3): q1_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2_vecs) < 3): q2_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1_vecs + q2_vecs


def parse_YToHH(event):
    GenPartsColl = Collection(event, "GenPart")

    #Y -> HH, H-> tt
    H1 = H2 = None
    H_ID = 25
    Y_ID = 39
    q1s = []
    q2s = []

    for i, gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        if(isFirstCopy(gen_part.statusFlags) and fromHardProcess(gen_part.statusFlags)):
        #if(isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0  and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):
            if(abs(gen_part.pdgId) == H_ID):
                if(H1 is None): H1 = gen_part
                elif(H2 is None): H2 = gen_part
                else: print("Extra H!")
                #print("Z", gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, gen_part.genPartIdxMother)


    parent_ids = {top_ID, W_ID}

    #avoid low pt garbage from shower

    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 and fromHardProcess(gen_part.statusFlags)
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):

            H_cand, H_dist = findMother(GenPartsColl, gen_part, {H_ID, Y_ID}, dist=0)
            if(H1 is not None and H2 is not None):
                if(H_cand is H1 ): q1s.append((H_dist, gen_part))
                elif(H_cand is H2 ): q2s.append((H_dist, gen_part))
            else:
                #No H candidates saved, just do the best we can to split the partons between the two H candidates
                #b quarks usually come first, 2 per H cand
                if(abs(gen_part.pdgId) == b_ID and  abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == top_ID and len(q1s) < 2): q1s.append((len(q1s), gen_part))
                elif(abs(gen_part.pdgId) == b_ID and  abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == top_ID): q2s.append((len(q2s), gen_part))
                elif(abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == W_ID and len(q1s) < 6): q1s.append((len(q1s), gen_part))
                elif(abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) == W_ID): q2s.append((len(q2s), gen_part))


            #mother, dist = findMother(GenPartsColl, gen_part, {H_ID}, dist=0)
            ##print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            #if(mother is top1 or mother is Z1): q1s.append((dist, gen_part))
            #elif(mother is top2 or mother is Z2): q2s.append((dist, gen_part))

    #gen matching isn't always perfect, do some attempt at cleanup here
    if(len(q1s) != 6 or len(q2s) != 6):
        print("Issue in quark finding!")
        print(len(q1s), len(q2s))
        print(q1s)
        print(q2s)
        for i, gen_part in enumerate(GenPartsColl):
            print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        exit(1)

    if(len(q1s) > 6): q1s = prune_genparts(q1s, 6)
    if(len(q2s) > 6): q2s = prune_genparts(q2s, 6)
    q1_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1_vecs) < 6): q1_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2_vecs) < 6): q2_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1_vecs + q2_vecs




def parse_ZpToTpTp(event):
    GenPartsColl = Collection(event, "GenPart")

    #Zp -> Tp Tp, Tp-> t Z
    Tp1 = Tp2 = top1 = top2 = Z1 = Z2 = None
    q1s = []
    q2s = []


    for i, gen_part in enumerate(GenPartsColl):
        if(isFirstCopy(gen_part.statusFlags) and fromHardProcess(gen_part.statusFlags)):
        #if(isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0  and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):
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


    parent_ids = {top_ID, W_ID, Z_ID}

    #avoid low pt garbage from shower

    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 and fromHardProcess(gen_part.statusFlags)
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):
            mother, dist = findMother(GenPartsColl, gen_part, {Z_ID, top_ID}, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(mother is top1 or mother is Z1): q1s.append((dist, gen_part))
            elif(mother is top2 or mother is Z2): q2s.append((dist, gen_part))

    #gen matching isn't always perfect, do some attempt at cleanup here
    if(len(q1s) != 5 or len(q2s) != 5):
        print("Issue in quark finding!")
        print(len(q1s), len(q2s))
        print(q1s)
        print(q2s)

    if(len(q1s) > 5): q1s = prune_genparts(q1s, 5)
    if(len(q2s) > 5): q2s = prune_genparts(q2s, 5)
    q1_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1_vecs) < 5): q1_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2_vecs) < 5): q2_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1_vecs + q2_vecs



def parse_Qstar(event):
    GenPartsColl = Collection(event, "GenPart")

    Radion = W1 = W2 = W_ISO = None
    #decay to W + q
    #qstar IDs 4000000 + q ID
    Qstar_up = 4000001
    Qstar_down = 4000002
    Qstar_strange = 4000003
    Qstar_charm = 4000004
    Qstar_b = 4000004
    W_ID = 24
    q1s = []
    q2s = []

    parent_ids = {W_ID, Qstar_up, Qstar_down, Qstar_strange, Qstar_charm, Qstar_b}


    for i, gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 
                and (abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids) and isFinal(GenPartsColl[gen_part.genPartIdxMother].statusFlags)):
            mother, dist = findMother(GenPartsColl, gen_part, parent_ids, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(abs(mother.pdgId) == W_ID): q1s.append((dist, gen_part))
            elif(abs(mother.pdgId) >= Qstar_up ): q2s.append((dist, gen_part))


    if(len(q1s) != 2 or len(q2s) != 1):
        print("Issue in quark finding!")
        print(len(q1s), len(q2s))
        print(q1s)
        print(q2s)

    #gen matching isn't always perfect, do some attempt at cleanup here
    if(len(q1s) > 2): q1s = prune_genparts(q1s, 2)
    if(len(q2s) > 1): q2s = prune_genparts(q2s, 1)

    q1s_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2s_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1s_vecs) < 2): q1s_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2s_vecs) < 1): q2s_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1s_vecs + q2s_vecs

def parse_XYY(event):
    GenPartsColl = Collection(event, "GenPart")

    Radion = W1 = W2 = W_ISO = None
    #decays use higgs and Z as Y/Y'
    X_ID = 9000001
    Y_ID = 23
    YPrime_ID = 25
    q1s = []
    q2s = []

    parent_ids = {Y_ID, YPrime_ID}


    for i, gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids):
            mother, dist = findMother(GenPartsColl, gen_part, parent_ids, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(mother.pdgId == Y_ID): q1s.append((dist, gen_part))
            elif(mother.pdgId == YPrime_ID): q2s.append((dist, gen_part))

    if(len(q1s) != 2 or len(q2s) != 2):
        print("Issue in quark finding!")
        print(len(q1s), len(q2s))
        print(q1s)
        print(q2s)

    #gen matching isn't perfect, do some attempt at cleanup here
    if(len(q1s) > 2): q1s = prune_genparts(q1s, 2)
    if(len(q2s) > 2): q2s = prune_genparts(q2s, 2)

    q1s_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q1s ]
    q2s_vecs =  [ [gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId] for dist,gen_part in q2s ]

    #zero pad if we missed some quarks
    while(len(q1s_vecs) < 2): q1s_vecs.append([-1.0, 0.0, 0.0, 0])
    while(len(q2s_vecs) < 2): q2s_vecs.append([-1.0, 0.0, 0.0, 0])

    return q1s_vecs + q2s_vecs


def parse_Wkk(event):
    GenPartsColl = Collection(event, "GenPart")

    Radion = W1 = W2 = W_ISO = None
    WKK_ID = 9000024
    RADION_ID = 9000025
    qs_iso = []
    qs_radion = []


    for i, gen_part in enumerate(GenPartsColl):
        #print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass)
        if(isFirstCopy(gen_part.statusFlags) and fromHardProcess(gen_part.statusFlags)):
            if(abs(gen_part.pdgId) == W_ID and gen_part.genPartIdxMother >= 0):
                m = GenPartsColl[gen_part.genPartIdxMother]
                if(abs(m.pdgId) == WKK_ID or gen_part.genPartIdxMother==0 ): W_ISO = gen_part
                elif(abs(m.pdgId) == RADION_ID): 
                    if(W1 is None): W1 = gen_part
                    elif(W2 is None): W2 = gen_part
                    else: print("Extra W!")    

    parent_ids = {W_ID}


    for gen_part in GenPartsColl:
        #if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0):
        if(abs(gen_part.pdgId) < MAX_LEP_ID and isFirstCopy(gen_part.statusFlags) and gen_part.genPartIdxMother > 0 
                and abs(GenPartsColl[gen_part.genPartIdxMother].pdgId) in parent_ids and fromHardProcess(gen_part.statusFlags)):
            mother, dist = findMother(GenPartsColl, gen_part, parent_ids, dist=0)
            #print(gen_part.pt, gen_part.eta, gen_part.phi, gen_part.pdgId, GenPartsColl[gen_part.genPartIdxMother].pdgId, mother.pdgId)
            if(mother is W_ISO): qs_iso.append((dist, gen_part))
            elif(mother is W1 or mother is W2): qs_radion.append((dist, gen_part))

    #gen matching isn't perfect, do some attempt at cleanup here

    if(len(qs_iso) != 2 or len(qs_radion) != 4):
        print("Issue in quark finding!")
        print(qs_iso)
        print(qs_radion)
        for i,gen_part in enumerate(GenPartsColl):
            print(i, gen_part.pdgId, gen_part.genPartIdxMother, gen_part.pt, gen_part.eta, gen_part.phi, gen_part.mass, fromHardProcess(gen_part.statusFlags), isFirstCopy(gen_part.statusFlags))




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
        'YtoHH' : (6*2, parse_YToHH),
        'ZpToTpTp' : (5*2, parse_ZpToTpTp),
        'Wkk' : (4+2, parse_Wkk),
        'XYY' : (2+2, parse_XYY),
        'Qstar' : (2+1, parse_Qstar),
        'Wp' : (3+3, parse_Wp),
         }

