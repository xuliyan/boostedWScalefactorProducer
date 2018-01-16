#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True
from importlib import import_module
from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor

from  selectionModule import *
p=PostProcessor(".",["TT_sal.root"],"nFatJet>0&&FatJet_msoftdrop>30&&FatJet_msoftdrop<140&&FatJet_eta<2.5&&FatJet_pt>200&&MET_sumEt>40","keep_and_drop.txt",[selectionModule()],provenance=True)
p.run()
