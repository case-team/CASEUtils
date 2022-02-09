# H5 Maker

Make h5's from the extended NanoAOD (pancakes) after applying a preselection. 

## Preselection

The preselection is as follows

Filters: 

```
Flag\_goodVertices, Flag\_globalTightHalo2016Filter, Flag\_eeBadScFilter, Flag\_HBHENoiseFilter, Flag\_HBHENoiseIsoFilter, Flag\_ecalBadCalibFilter, Flag\_EcalDeadCellTriggerPrimitiveFilter, Flag\_BadChargedCandidateFilter
```

Triggers:

2016: ```HLT\_PFHT800, HLT\_PFHT900, HLT\_PFJet450, HLT\_PFJet500, HLT\_AK8PFJet450, HLT\_AK8PFJet500 ```

2017/8: ```HLT\_PFHT1050, HLT\_PFJet500, HLT\_AK8PFJet380\_TrimMass30, HLT\_AK8PFJet400\_TrimMass30 ```

Mjj > 1200.

2 AK8 Jets with tight ID,  pt > 300 , |eta| < 2.5


## Output
The output is an h5 dataset with several different keys
In general the highest two pt jets are the dijet candidate.
The jets are labeled so that 'j1' is the more massive jet in the event and 'j2'
the less massive jet


The data in each keys are:


**preselection\_eff**: Single float for whole file. Efficiency of the preselection requirements on this sample

**truth\_label** :  single int. Labels  the type of event. 0 is QCD, signals are => 1 (depends on dataset), other backgrounds are TBD but will be < 0

**event\_info** : 6 floats. [eventNum, MET, MET\_phi, genWeight, leptonic\_decay, runNum]
leptonic decay is a flag that will check the generator level
decay to see if it is leptonic / semi-leptonic or not (1 for leptonic, 0 for full hadronic)

**jet\_kinematics** :  14 floats. Mjj, delta\_eta (between j1 and j2), followed by the 4 vectors of j1, j2 and j3  in pt,eta,phi,m\_softdrop format (if no 3rd jet passing cuts, zeros)

**jet1(2)\_extraInfo** :  7 floats. 4 nsubjettiness scores (tau\_i), followed by LSF3,
btag score (max deepB score of ak4 subjets) and number of PF constituents for j1 (j2)

**jet1(2)\_PFCands** : 400 floats. 4 vectors of (up to) 100 PFCands of j1 (j2) in  Px, Py, Pz, E  format. Zero
padded


If `--sys` flag is used, additional columns with info necessary for systematics
computation are added

**sys\_weights** 19 floats: See `sys\_weights\_map` dictionary inside H5\_maker.py for map of variable name to index. 
    "nom_weight" is the nominal weight and all others are the weights for given systematic variation.
WARNING: Careful with the overall normalizations, they are not accurate currently (including nom_weight)! (eg some weights may average to 1.3 instead of 1.0 )
        Best to work with relative normalizations/shape changes for now.
        ```[nom_weight, pdf_up, pdf_down, prefire_up, prefire_down, pileup_up, pileup_down, btag_up, btag_down, 
        PS_ISR_up, PS_ISR_down, PS_FSR_up, PS_FSR_down, F_up, F_down, R_up, R_down, RF_up, RF_down] ```

**jet1(2)\_JME\_vars**:  12 floats. See `JME\_weights\_map` dictionary inside H5_maker.py for map of variable name to index.
The standard corrections are applied if the systematics turned on, this stores
variations.
        ```[pt_JES_up m_JES_up pt_JES_down m_JES_down pt_JER_up m_JER_up pt_JER_down m_JER_down m_JMS_up m_JMS_down m_JMR_up m_JMR_down]```



## H5_merge
H5\_merge.py is also useful. It combines different h5 files together like hadd.
Does a weighted average for preselection\_eff
