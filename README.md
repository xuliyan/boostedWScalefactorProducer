
### How to run the W-tagging scalefactor code ###
#########################################

## installation instructions
Setup CMSSW and get nanoAOD packages
```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv

git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py .

scram build

cd PhysicsTools/NanoAODTools/

git remote add sal https://github.com/rappoccio/nanoAOD-tools.git
git fetch sal
git checkout -b TTbarResHad remotes/sal/TTbarResHad


```

### getting the code

```
git clone -b nanoAOD git@github.com:BoostedScalefactors/WTopScalefactorProducer.git
cd WTopScalefactorProducer

#For public version:
#git clone https://github.com/${GITUSER}/WTopScalefactorProducer 
#cd WTopScalefactorProducer
#git remote add originalRemote https://github.com/BoostedScalefactors/WTopScalefactorProducer.git
#git fetch originalRemote
#git checkout -b nanoOAD originalRemote/nanoAOD
#git fetch originalRemote

```
### Producing samples

See README in subdirectory skim_nanoAOD/
```
cd skim_nanoAOD/
python process_nanoAOD.py
```

### Working locally (without CMSSW, just python2.7 and ROOT)
```
cd PhysicsTools/NanoAODTools/
(JUST ONCE:)
bash standalone/env_standalone.sh build
(EVERY TIME:)
source standalone/env_standalone.sh
```

### running scalefactor code

```
python Automatic_Setup.py #To compile
python runSF.py -b   #To run
```

The basic script to be run is 

```
python runSF.py
```
It takes as input .root files containing a TTree with a branch for the mass distribution you want to calculate a scale factor for. This branch can contain events after full selection is applied, or new selections can be implemented on the fly in wtagSFfits.py. In addition to a data and the separate background MC files, you need one file called "*pseudodata* which contains all MC added together (with their appropriate weights, using ROOT hadd).

   
   General Options:
```
    -b : To run without X11 windows
    -c : channel you are using(electron,muon or electron+muon added together)
    --HP : HP working point
    --LP : LP working point
    --fitTT : Only do fits to truth matched tt MC
    --fitMC : Only do fits to MC (test fit functions)
    --sample : name of TT MC eg --sample "herwig"
    --doBinned : to do binned simultaneous fit (default is unbinned)
    --76X : Use files with postfix "_76X" (change to postfix of your choice if running on several different samples)
    --useDDT : Uses DDT tagger instead of pruning+softdrop (ops! Requires softdrop variables)
    --usePuppiSD : Uses PUPPI + softdrop and PUPPI n-subjettiness
```
