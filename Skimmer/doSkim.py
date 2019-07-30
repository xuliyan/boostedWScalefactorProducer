#!/Usr-/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 
from WTopScalefactorProducer.Skimmer.skimmer import Skimmer
from optparse import OptionParser

#sys.stdout = open('./testskimoutput.txt', 'w')

parser = OptionParser()

parser.add_option('--infiles', '-i', action="store", default="", help="Input files")
parser.add_option('--output', '-o', dest="outputDir", action="store",  default="TEST", help="Output dierctory")
parser.add_option('--name', dest="name", action="store", default="test.root", help="Name of output file")
parser.add_option('--chunk', dest="chunk", action="store", default="", help="Number of chunks")


args, other = parser.parse_args()

if len(args.infiles)>1:
   infile = args.infiles.split(',')
else:
  infile = [
#      "root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/NANOAOD/Nano14Dec2018_ver2-v1/40000/4BF53ABE-BB2D-B147-8BDE-888B21C3E07C.root"
#      "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv4/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/120000/BB978C6D-1770-CB42-A45E-AFAA249DAA82.root"
#      "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv4/ZZ_TuneCP5_13TeV-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/110000/56817FDC-2D4C-7E4B-ABA0-9A1F7BA505C4.root"
#      "root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv4/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/60000/29C0425B-F520-D141-BC3E-8E46A9E88E87.root Skimmed_2019_06_15/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8 ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8"
      #"root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18NanoAODv4/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/Nano14Dec2018_102X_upgrade2018_realistic_v16-v1/120000/BB978C6D-1770-CB42-A45E-AFAA249DAA82.root"
       "file:/work/mhuwiler/data/WScaleFactors/LocalTestFiles/BB978C6D-1770-CB42-A45E-AFAA249DAA82.root"
  ]


outputDir = args.outputDir
name = args.name
chunk = args.chunk

print infile, outputDir, name, chunk


#"keep_and_drop.txt"
#FatJet_msoftdrop>30&&MET_sumEt>40&&    ((nElectron>0&&HLT_Ele115_CaloIdVT_GsfTrkIdT)||
if infile[0].find("SingleMuon")!=-1:
  channel = "mu"
  print "Processing a Single Muon dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt')
                   
elif infile[0].find("EGamma")!=-1:
  channel = "el"
  print "Processing a Single Electron dataset file..."
  p=PostProcessor(outputDir, infile, None, None,
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False,
                    jsonInput='/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/PromptReco/Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt')
                  
elif infile[0].find("SemiLeptonic")!=-1:
  channel = "elmu"
  print "Processing a Muon and Electron dataset file..."
  p=PostProcessor(outputDir, infile, "" ,"",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
else:
  print "Processing MC..."
  channel = "elmu"
  p=PostProcessor(outputDir, infile, "" ,"",
                    modules=[Skimmer(channel)],provenance=False,fwkJobReport=False)
p.run()
print "DONE"
