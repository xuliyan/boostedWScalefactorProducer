#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from boostedWScalefactorProducer.Skimmer.skimmer import Skimmer
if len(sys.argv)>1:
    infile = ["root://cms-xrd-global.cern.ch/"+sys.argv[1]]
    print infile
else:
    infile = ["F.root"]
    print infile

if len(sys.argv)>2: outputDir = sys.argv[2]
else: outputDir = "TEST"	

if len(sys.argv)>3: name = sys.argv[3]
else: name = "test.root"

if len(sys.argv)>4: chunck = sys.argv[4]
else: chunck = ""  

# HLT_Mu50&&nMuon>0&&Muon_pt[0]>55.&&Muon_pfRelIso03_chg[0]<0.15&&Muon_highPtId>1&&nFatJet>0&&FatJet_pt>200

if infile[0].find("SingleMuon")!=-1:
    channel = "mu"
    print "Processing a Single Muon dataset file..."
    p=PostProcessor(outputDir, infile, None, None, #"HLT_Mu50 && nMuon>0 && Muon_pt[0]>55. && nFatJet>0"
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSONv1.txt',
                    #jsonInput='Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt',
                    )

elif infile[0].find("SingleElectron")!=-1:
    channel = "el"
    print "Processing a Single Electron dataset file..."
    p=PostProcessor(outputDir, infile, None, None, #"(event.HLT_Ele32_WPTight_Gsf || event.HLT_Ele35_WPTight_Gsf || event.HLT_Ele40_WPTight_Gsf || HLT_Ele115_CaloIdVT_GsfTrkIdT) && nElectron>0 && Electron_pt[0]>55. && nFatJet>0"
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSONv1.txt',
                    )

else:
    print "Processing MC..."
    channel = "elmu"
    p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
p.run()
print "DONE"
