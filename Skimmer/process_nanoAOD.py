#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from  selectionModule import *

if len(sys.argv)>1:
	 infile = sys.argv[1]
else:
	infile = "root://cms-xrd-global.cern.ch/store/user/thaarres/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv2-PUMoriond1780XmcRun2asymptotic/180116_083632/0000/test80X_NANO_99.root"
if len(sys.argv)>2:
	 outputDir = sys.argv[2]
else:
	outputDir = "TT_TuneCUETP8M2T4_13TeV-powheg-pythia8"	
if len(sys.argv)>3:
	 name = sys.argv[3]
else:
	name = "tree"
p=PostProcessor(outputDir,[infile],"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_msoftdrop<140&&FatJet_eta<2.5&&FatJet_pt>200&&MET_sumEt>40","keep_and_drop.txt",[selectionModule()],postfix=name,provenance=True,haddFileName=name+".root",jsonInput="Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt")
p.run()
