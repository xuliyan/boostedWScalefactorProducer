#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from WTopScalefactorProducer.Skimmer.skimmer import Skimmer
if len(sys.argv)>1:
   infile = sys.argv[1].split(',')
else:
  infile = ["root://cms-xrd-global.cern.ch//store/mc/RunIIFall17NanoAOD/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/0686D954-B442-E811-9AB3-0CC47AD98BC2.root"]

if len(sys.argv)>2: outputDir = sys.argv[2]
else: outputDir = "Samples"	

if len(sys.argv)>3: name = sys.argv[3]
else: name = "test.root"

if len(sys.argv)>4: chunck = sys.argv[4]
else: chunck = ""  

if infile[0].find("SingleMuon")!=-1:
  channel = "mu"
  print "Processing a Single Muon dataset file..."
  p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&&nMuon>0&&HLT_Mu50" ,"keep_and_drop.txt",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                   jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')
                   
elif infile[0].find("SingleElectron")!=-1:
  channel = "el"
  print "Processing a Single Electron dataset file..."
  p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&&nElectron>0&&HLT_Ele115_CaloIdVT_GsfTrkIdT" ,"keep_and_drop.txt",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/PromptReco/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt')
                    
else:
  print "Processing MC..."
  channel = "mu"
  p=PostProcessor(outputDir, infile,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&&((nElectron>0&&HLT_Ele115_CaloIdVT_GsfTrkIdT)||(nMuon>0&&HLT_Mu50))" ,"keep_and_drop.txt",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
p.run()
print "DONE"
