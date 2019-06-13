#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from WTopScalefactorProducer.Skimmer.skimmer import Skimmer
if len(sys.argv)>1:
   infile = sys.argv[1].split(',')
else:
  infile = [
      "root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/NANOAOD/Nano14Dec2018_ver2-v1/40000/4BF53ABE-BB2D-B147-8BDE-888B21C3E07C.root"
  ]

if len(sys.argv)>2: outputDir = sys.argv[2]
else: outputDir = "TEST"	

if len(sys.argv)>3: name = sys.argv[3]
else: name = "test.root"

if len(sys.argv)>4: chunck = sys.argv[4]
else: chunck = ""  

#"keep_and_drop.txt"
#FatJet_msoftdrop>30&&MET_sumEt>40&&    ((nElectron>0&&HLT_Ele115_CaloIdVT_GsfTrkIdT)||
if infile[0].find("SingleMuon")!=-1:
  channel = "mu"
  print "Processing a Single Muon dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt')
                   
elif infile[0].find("SingleElectron")!=-1:
  channel = "el"
  print "Processing a Single Electron dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt')
                    
else:
  print "Processing MC..."
  channel = "mu"
  p=PostProcessor(outputDir, infile, "" ,"",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
p.run()
print "DONE"
