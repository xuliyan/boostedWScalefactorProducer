# Skimmer

Skimmer running on UZH ntuples to create input for W/top scalefactor code

## Getting started

You need ROOT or CMSSW to get started, here CMSSW will be used, and the installation will be done into a new directory:
```
cd Skimmer
source /cvmfs/cms.cern.ch/cmsset_default.sh
cmsrel CMSSW_8_0_20
cd CMSSW_8_0_20/src
cmsenv
cd ../..
```
Get SFrame, initialise submodules and create some necessary repositories:
```
source init.sh
```

Setup and compile SFrame:
```
cd SFrame
source setup.sh
make
cd ..
```

## Compilation

To compile all submodules:
```
source make.sh
```
To ```make distclean``` for all submodules:
```
source makedistclean.sh
```
Mind that SFrame itself and any other directory, which is not a submodule, will not be touched.

## Further documentation

Please read https://git-scm.com/book/en/v2/Git-Tools-Submodules for information on how to work with submodules.
