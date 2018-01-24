### How to run the W-tagging scalefactor code ###
#########################################

## Installation instructions
Setup CMSSW and get nanoAOD packages
```
cmsrel CMSSW_9_4_2
cd CMSSW_9_4_2/src
git cms-init
git remote add cms-nanoAOD https://github.com/cms-nanoAOD/cmssw.git
git fetch cms-nanoAOD
git checkout -b nanoAOD remotes/cms-nanoAOD/master
git cms-addpkg PhysicsTools/NanoAOD
scram b -j 10
```

### getting the code
First fork your own version of the repository at https://github.com/BoostedScalefactors/WTopScalefactorProducer
```
export GITUSER=`git config user.github`
echo "Your github username has been set to \"$GITUSER\""
git clone -b nanoOAD git@github.com:${GITUSER}/WTopScalefactorProducer.git
cd WTopScalefactorProducer
git remote add originalRemote git@github.com:BoostedScalefactors/WTopScalefactorProducer.git
```



### Producing samples

First you need to produce your input files by skimming nanoAOD samples. For this, see README in subdirectory Skimmer/.
The syntax is: python process_nanoAOD.py <infile> <outdir> <outtreename>
```
cd Skimmer/
python process_nanoAOD.py
```

### running scalefactor code

The fitting code is located in Fitter/. For scalefactors from merged W AK8 jet, see README in Fitter/partiallyMerged. For scalefactors from merged top AK8 jet, see README in Fitter/fullyMerged
