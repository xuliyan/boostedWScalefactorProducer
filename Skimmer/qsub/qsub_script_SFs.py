#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

sys.path.append("../python")
from TTSkimmer import *


if len(sys.argv)>1:
	 infile = sys.argv[1]
else:
	infile = "root://cms-xrd-global.cern.ch//store/user/asparker/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/TTToSemiLeptonicTuneCP5PSweights13TeV-powheg-pythia8/180130_175206/0000/80XNanoV0-TTbar_SemiLep_102.root"
if len(sys.argv)>2:
	 outputDir = sys.argv[2]
else:
	outputDir = "TTToSemiLeptonic_TuneCP5_PSweights_13TeV"	
if len(sys.argv)>3:
	 name = sys.argv[3]
else:
	name = "tree"
  
# NOTE : This file is configured to Process MC.
# If you want to process Data then ADD "jsonInput='Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'"
# to the PostProcessor arguement below

p=PostProcessor(outputDir, [infile], 
                 "nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& ( (nElectron > 0 && HLT_Ele115_CaloIdVT_GsfTrkIdT) || (nMuon > 0 && HLT_Mu50))" ,"keep_and_drop.txt",
                  modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False ,
                  haddFileName=name+".root" )

p.run()

print "DONE"
