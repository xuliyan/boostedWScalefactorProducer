#!/usr/bin/env python
import os, sys
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import PostProcessor
from eventSelector import *
from PhysicsTools.NanoAODTools.postprocessing.modules.btv.btagSFProducer import *
from PhysicsTools.NanoAODTools.postprocessing.modules.jme.jecUncertainties import *

files=["../nanolzma_1.root"]

selection='''
 (Sum$(Electron_pt > 20 && Electron_mvaSpring16GP_WP80) >= 1   ||
 Sum$(Muon_pt > 20 && Muon_tightId) >= 1 ) &&
 Sum$((abs(Jet_eta)<2.5 && Jet_pt > 20)) >= 2
'''

p=PostProcessor(".",files,selection.replace('\n',' '),"keep_and_drop.txt",[eventSel()])

p.run()
