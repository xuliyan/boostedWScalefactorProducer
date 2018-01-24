#!/usr/bin/env python
import os,sys
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import * 

#this takes care of converting the input files from CRAB
from PhysicsTools.NanoAODTools.postprocessing.framework.crabhelper import inputFiles,runsAndLumis

# from  PhysicsTools.NanoAODTools.postprocessing.analysis.smp.xs.TTbar_SemiLep import *
# p=PostProcessor(".",inputFiles(),"Jet_pt>200",modules=[TTbar_SemiLep()],provenance=False,fwkJobReport=False,jsonInput=runsAndLumis(),histFileName= requestName() + '_TTbar_SemiLep_hists.root', histDirName='ttbar_semilep'  , postfix='skimmed-TTbar_SemiLep', haddFileName= requestName() + '_TTbar_SemiLep_hadded.root'  )

sys.path.append("../")
from  selectionModule import *
p=PostProcessor(".",inputFiles(),"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_msoftdrop<140&&FatJet_eta<2.5&&FatJet_pt>200&&MET_sumEt>40","../keep_and_drop.txt",[selectionModule()],provenance=True,fwkJobReport=True)

p.run()
