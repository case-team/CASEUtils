# H5 Maker

Make h5's from the extended NanoAOD (pancakes) after applying a preselection. 

To run locally, try something like `python make_h5_local.py -i input_file.root -o output_file.h5 -y 201X -f 0` 

The option `-f` sets the truth label of the output. Usually signal is 1, QCD is 0, single top -1, ttbar -2, V+jets -3
Use the `--sys` option if you are running on a signal and want to save the event weights for the systematic variation. 
The options are described a bit in `make_h5_local.py`

## Preselection

The preselection is as follows

Filters: 

```
Flag\_goodVertices Flag\_globalSuperTightHalo2016Filter, Flag\_HBHENoiseFilter, Flag\_HBHENoiseIsoFilter, Flag\_EcalDeadCellTriggerPrimitiveFilter, 
Flag\_BadPFMuonFilter, Flag\_BadPFMuonDzFilter, Flag\_eeBadScFilter, Flag\_CSCTightHaloFilter (2016 only), Flag\_ecalBadCalibFilter (2017/8 only)

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
**d\_eta\_eff** : Single float for whole file. Efficiency of |dEta| < 1.3 cut (after preselection)

**truth\_label** :  single int. Labels  the type of event. 0 is QCD, signals are => 1 (depends on dataset), other backgrounds are TBD but will be < 0

**event\_info** : 8 floats. [eventNum, MET, MET\_phi, genWeight, leptonic\_decay, runNum, year, nJets]
leptonic decay is a flag that will check the generator level
decay to see if it is leptonic / semi-leptonic or not (1 for leptonic, 0 for full hadronic)
nJets is the number of Ak8 jets with pt > 50, |eta| < 2.5 and passing tight ID

**jet\_kinematics** :  14 floats. Mjj, delta\_eta (between j1 and j2), followed by the 4 vectors of j1, j2 and j3  in pt,eta,phi,m\_softdrop format (if no 3rd jet passing cuts, zeros)

**jet1(2)\_extraInfo** :  7 floats. 4 nsubjettiness scores (tau\_i), followed by LSF3,
btag score (max deepB score of ak4 subjets) and number of PF constituents for j1 (j2)

**jet1(2)\_PFCands** : 400 floats. 4 vectors of (up to) 100 PFCands of j1 (j2) in  Px, Py, Pz, E  format. Zero
padded


If `--sys` flag is used, additional columns with info necessary for systematics
computation are added

**sys\_weights** 21 floats: See `sys_weights_map` dictionary inside H5\_maker.py for map of variable name to index. 
    "nom_weight" is the nominal weight and all others are the weights for given systematic variation.
        
        ```[nom_weight, pdf_up, pdf_down, prefire_up, prefire_down, pileup_up, pileup_down, btag_up, btag_down, 
        PS_ISR_up, PS_ISR_down, PS_FSR_up, PS_FSR_down, F_up, F_down, R_up, R_down, RF_up, RF_down, top_ptrw_up, top_ptrw_down] ```

**jet1(2)\_JME\_vars**:  12 floats. See `JME_weights_map` dictionary inside H5_maker.py for map of variable name to index.
The nominal corrections are applied if the systematics turned on, this stores
variations.

        ```[pt_JES_up, m_JES_up, pt_JES_down, m_JES_down, pt_JER_up, m_JER_up, pt_JER_down, m_JER_down, m_JMS_up, m_JMS_down, m_JMR_up, m_JMR_down]```

Variations of the **preselection\_eff** for JES and JER systematic variations are also saved (eg `preselection_eff_JES_up`)




## H5_merge
H5\_merge.py is also useful. It combines different h5 files together like hadd.
Does a weighted average for preselection\_eff
