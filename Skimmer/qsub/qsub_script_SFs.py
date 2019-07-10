#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from WTopScalefactorProducer.Skimmer.skimmer import Skimmer
if len(sys.argv)>1:
   infile = sys.argv[1].split(',')
else:
  # infile = ["root://cms-xrd-global.cern.ch//store/user/asparker/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/TTToSemiLeptonicTuneCP5PSweights13TeV-powheg-pythia8/180130_175206/0000/80XNanoV0-muSkim_1.root"]
  infile = ["SingleMuon.root"]

if len(sys.argv)>2:
	 outputDir = sys.argv[2]
else:
	outputDir = "TEST"	

if len(sys.argv)>3:
	 name = sys.argv[3]
else:
	name = "test.root"

if len(sys.argv)>4:
	 chunck = sys.argv[4]
else:
	chunck = ""  

# HLT_Mu50&&nMuon>0&&Muon_pt[0]>55.&&Muon_pfRelIso03_chg[0]<0.15&&Muon_highPtId>1&&nFatJet>0&&FatJet_pt>200

if infile[0].find("SingleMuon")!=-1:
  channel = "mu"
  print "Processing a Single Muon dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/work/zucchett/WTagging/CMSSW_10_2_6/src/WTopScalefactorProducer/Skimmer/python/JSON/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt',
                    )
                   
elif infile[0].find("SingleElectron")!=-1:
  channel = "el"
  print "Processing a Single Electron dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
                    # jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')

else:
  print "Processing MC..."
  channel = "mu"
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
p.run()
print "DONE"
