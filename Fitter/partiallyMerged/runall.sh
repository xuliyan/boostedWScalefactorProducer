#!/bin/bash

# cd /uscms_data/d3/xyan/WTagging/TEST
# ./haddall.sh

cd /uscms/homes/x/xyan/work/CMSSW_10_2_6/src/boostedWScalefactorProducer/Fitter/partiallyMerged
python addWeight.py

# python runSF_nanoAOD.py -b --doBinned --tagger SelectedJet_tau21 --workspace workspace_tau21_0p35_0p75 --HP 0.35 --LP 0.75 --doWS
