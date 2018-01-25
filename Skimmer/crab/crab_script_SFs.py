#!/Usr-/bin/env python
import os
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

from  TTbar_SemiLep import *


p=PostProcessor(".", inputFiles()  ,"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_pt>200&&MET_sumEt>40&& ( (nElectron > 0 && HLT_Ele115_CaloIdVT_GsfTrkIdT) || (nMuon > 0 && HLT_Mu50))" ,modules=[TTbar_SemiLep()],provenance=True,fwkJobReport=True,jsonInput='Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON.txt' ,histFileName= 'TTbar_SemiLep_hists.root', histDirName='ttbar_semilep'  , haddFileName=  '80XNanoV0-TTbar_SemiLep.root'  )

p.run()

print "DONE"
os.system("ls -lR")
