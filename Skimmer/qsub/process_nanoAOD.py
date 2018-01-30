#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from PhysicsTools.NanoAODTools.WTopScalefactorProducer.Skimmer.TTbar_SemiLep import *
if len(sys.argv)>1:
     infileis = sys.argv[1]
else:
     infileis = ["root://cmsxrootd.fnal.gov//store/user/srappocc/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8RunIISummer16MiniAODv2-PUMoriond1780XmcRun2/180112_155912/0000/test80X_NANO_10.root"]
     #"root://cms-xrd-global.cern.ch//store/user/thaarres/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/TTTuneCUETP8M2T413TeV-powheg-pythia8RunIISummer16MiniAODv2-PUMoriond1780XmcRun2asymptotic/180116_083632/0000/test80X_NANO_99.root"
if len(sys.argv)>2:
    outputDir = sys.argv[2]
else:
    outputDir = "TTJetsTuneCUETP8M113TeV-madgraphMLM-pythia8"	
if len(sys.argv)>3:
    name = sys.argv[3]
else:
    name = "tree"
p=PostProcessor(".", infileis  ,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& ( (nElectron > 0 && HLT_Ele115_CaloIdVT_GsfTrkIdT) || (nMuon > 0 && HLT_Mu50))" ,modules=[TTbar_SemiLep()],provenance=True,fwkJobReport=True ,histFileName= 'TTbar_SemiLep_hists.root', histDirName='ttbar_semilep'  , haddFileName=  '80XNanoV0-TTbar_SemiLep.root'  )
p.run()
