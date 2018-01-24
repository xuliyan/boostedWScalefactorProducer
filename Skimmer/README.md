### How to run the W-tagging scalefactor code ###
#########################################

## installation instructions
Setup CMSSW and get nanoAOD packages
```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
cmsenv

git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools

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


For public version:
git clone https://github.com/${GITUSER}/WTopScalefactorProducer 
cd WTopScalefactorProducer
git remote add originalRemote https://github.com/BoostedScalefactors/WTopScalefactorProducer.git
git fetch originalRemote
git checkout -b nanoOAD originalRemote/nanoAOD
git fetch originalRemote
cd Skimmer/
ln -s $CMSSW_BASE/src/PhysicsTools/NanoAODTools/scripts/haddnano.py .
cd ..

```



### Producing samples

First you need to produce your input files by skimming nanoAOD samples. For this, see README in subdirectory Skimmer/.
The syntax is: python process_nanoAOD.py <infile> <outdir> <outtreename>. To submit with crab go to Skimmer/crab
```
cd Skimmer/
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

The fitting code is located in Fitter/. For scalefactors from merged W AK8 jet, see README in Fitter/partiallyMerged. For scalefactors from merged top AK8 jet, see README in Fitter/fullyMerged
