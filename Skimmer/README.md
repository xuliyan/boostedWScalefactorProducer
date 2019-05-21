
### Producing samples

## Running locally
First you need to produce your input files by skimming nanoAOD samples.

For local skimming tests, the syntax is (remember to change infile!):
```
 python doSkim.py
```

## Running with qsub
To submit batch jobs go to directory ```Skimmer/qsub```

Submit with:
```
python submit_qsub_official.py [ TT | ST | VV | Wjets | data ] outfolder
```  

The first time, run with option ```-c``` as second argument to create the filelists from DAS

Also remember to perform VOMS login: ```voms-proxy-init --voms cms --valid 48:00```

## Running with crab (obsolete)
To submit with crab go to directory Skimmer/crab and source crab3 if you have not already done so
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
Submit with:
```
python submit_all.py -c Pset_nanoSkim.py -d DirectoryName -f listOfDatasets.txt
```  

To make JMAR Skims (Loose skim with >=1 AK8 Jet with Pt > 200 GeV and >=1 Lepton ):(Skim these to create the files we need)
```
cd crab/JMARskim
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py .
python  submit_all_uif.py  -c  PSet.py -d June5_nanoskim-JetsAndLepton  -f test94X_DY_madgraph.txt
 
```


## Running outside of CMSSW (only require python2.7 and ROOT)
```
cd PhysicsTools/NanoAODTools/
(JUST ONCE:)
bash standalone/env_standalone.sh build
(EVERY TIME:)
source standalone/env_standalone.sh
```
