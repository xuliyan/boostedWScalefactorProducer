
### How to skim nanoAOD samples ###
#########################################

The selections are given in selectionModule.py. Additional cuts and input file is defined in process_nanoAOD.py
```
### Producing samples locally

```
python process_nanoAOD.py
```

### running with crab
Change input sample and SE in crab_cfg.py. Selections and module in crab_script.py
```
cd crab/
voms-proxy-init -voms cms --valid 200:00
source /cvmfs/cms.cern.ch/crab3/crab.sh
python submit_all.py -f listOfDatasets.txt
```