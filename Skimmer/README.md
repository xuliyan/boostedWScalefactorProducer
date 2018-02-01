
### Producing samples

## Running locally
First you need to produce your input files by skimming nanoAOD samples.

For local skimming tests, the syntax is (remember to change infile!):
```
 python crab/crab_script_SFs.py 
```

## Running with crab
To submit with crab go to directory Skimmer/crab and source crab3 if you have not already done so
```
source /cvmfs/cms.cern.ch/crab3/crab.sh
```
Submit with:
```
python submit_all.py -c Pset_nanoSkim.py -d DirectoryName -f listOfDatasets.txt
```  
## Running with qsub
To submit with qsyb go to directory Skimmer/qsub

Submit with:
```
python submit_qsub.py [ TT | ST | VV | Wjets | data ] outfolder
```  
(This is submitting process_nanoAOD.py,rather sync with A.)
## Running outside of CMSSW (only require python2.7 and ROOT)
```
cd PhysicsTools/NanoAODTools/
(JUST ONCE:)
bash standalone/env_standalone.sh build
(EVERY TIME:)
source standalone/env_standalone.sh
```
