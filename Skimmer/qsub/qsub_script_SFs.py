#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from WTopScalefactorProducer.Skimmer.TTSkimmer import *


if len(sys.argv)>1:
   infile = sys.argv[1].split(',')
else:
	infile = ["root://cms-xrd-global.cern.ch//store/user/asparker/TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8/TTToSemiLeptonicTuneCP5PSweights13TeV-powheg-pythia8/180130_175206/0000/80XNanoV0-TTbar_SemiLep_1.root"]
if len(sys.argv)>2:
	 outputDir = sys.argv[2]
else:
	outputDir = "TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8"	
if len(sys.argv)>3:
	 name = sys.argv[3]
else:
	name = "TTToSemiLeptonic.root"
if len(sys.argv)>4:
	 chunck = sys.argv[4]
else:
	chunck = ""  
# if infile.find("data")==-1:
#   p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& ( (nElectron > 0 && HLT_Ele115_CaloIdVT_GsfTrkIdT) || (nMuon > 0 && HLT_Mu50))" ,"keep_and_drop.txt",
#                     modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False ,
#                     haddFileName=  name+'.root' )
# else:
#   p=PostProcessor(outputDir,infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& ( (nElectron > 0 && HLT_Ele115_CaloIdVT_GsfTrkIdT) || (nMuon > 0 && HLT_Mu50))" ,"keep_and_drop.txt",
#                     modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False ,
#                     haddFileName=  name+'.root',jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt' )

if infile[0].find("SingleMuon")==-1:
  p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& nMuon>0 && HLT_Mu50" ,"keep_and_drop.txt",
                    modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False ,
                    haddFileName= name+chunck+'.root')
                    # histFileName=outputDir+"/genEV_"+name+'.root', histDirName="genEV_"+outputDir )
else:
  p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& nMuon>0 && HLT_Mu50" ,"keep_and_drop.txt",
                    modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False ,
                    haddFileName=  name+chunck+'.root',jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')
                    # histFileName=outputDir+"/genEV_"+name+'.root', histDirName="genEV_"+outputDir )
                    
p.run()

print "DONE"
