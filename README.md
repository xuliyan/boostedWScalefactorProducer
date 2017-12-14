
### How to run the W-tagging scalefactor code ###
#########################################

## installation instructions
Source ROOT version 5.34/36 independently (currently no CMSSW version with ROOT 5.34.X with X>18 available), e.g at PSI:
```
source /swshare/ROOT/root_v5.34.32_precompiled/root/bin/thisroot/.sh
```
### getting the code

```
export GITUSER=`git config user.github`
echo "Your github username has been set to \"$GITUSER\""
git clone git@github.com:${GITUSER}/WTopScalefactorProducer.git
cd WTopScalefactorProducer
git remote add originalRemote git@github.com:thaarres/WTopScalefactorProducer.git
git fetch originalRemote
git checkout -b DevelopmentBranch originalRemote/master
```

### running scalefactor code

```
cd boostedWScalefactorProducer
export ROOFITSYS=$ROOTSYS
python Automatic_Setup.py --clean
python Automatic_Setup.py #To compile
python wtagSFfits.py -b   #To run
```

The basic script to be run is 

```
python wtagSFfits.py
```
It takes as input .root files containing a TTree with a branch for the mass distribution you want to calculate a scalefactor for. This branch can contain events after full selection is applied, or new selections can be implemented on the fly in wtagSFfits.py. In addition to a data and the separate background MC files, you need one file called "*pseudodata* wchich contains all MC added together (with their appropriate weights, using ROOT hadd).

   
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
