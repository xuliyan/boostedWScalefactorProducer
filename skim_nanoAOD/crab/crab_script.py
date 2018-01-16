#!/usr/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

# from  PhysicsTools.NanoAODTools.postprocessing.examples.mhtProducer import *
# p=PostProcessor(".",inputFiles(),"Jet_pt>200",modules=[mht()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
# p.run()
sys.path.append("../")
from  selectionModule import *
p=PostProcessor(".",inputFiles(),"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_msoftdrop<140&&FatJet_eta<2.5&&FatJet_pt>200&&MET_sumEt>40","keep_and_drop.txt",[selectionModule()],provenance=True,fwkJobReport=True,jsonInput=runsAndLumis())
p.run()

print "DONE"
os.system("ls -lR")





