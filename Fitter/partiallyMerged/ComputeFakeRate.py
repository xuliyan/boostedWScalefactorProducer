#!/usr/bin/env python

import ROOT
from rootpy.tree import Cut

ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")


directory = "/eos/home-m/mhuwiler/data/Wtagging/Mergedefinition2017/"

basecut = Cut("SelectedJet_pt>480. && SelectedJet_pt<600. && passedMETfilters && maxAK4CSV>0.8484 && SelectedJet_mass>50 && SelectedJet_mass<130") # && genmatchedAK82017==1 "SelectedJet_pt>200 && SelectedJet_pt<10000 && passedMETfilters && maxAK4CSV>0.8484 && Jet_mass>50 && Jet_mass<130"
cut = Cut("SelectedJet_tau21<0.35")
cutHP = Cut("SelectedJet_tau21<0.35")
cutLP = Cut("SelectedJet_tau21<0.75 && SelectedJet_tau21>0.35")

#cutHP = Cut("SelectedJet_tau21<0.45")
#cutLP = Cut("SelectedJet_tau21<0.75 && SelectedJet_tau21>0.45")

#cutHP = Cut("SelectedJet_tau21_ddt_retune<0.43")
#cutLP = Cut("SelectedJet_tau21_ddt_retune<0.79 && SelectedJet_tau21_ddt_retune>0.43")

#cutHP = Cut("SelectedJet_tau21_ddt_retune<0.50")
#cutLP = Cut("SelectedJet_tau21_ddt_retune<0.80 && SelectedJet_tau21_ddt_retune>0.50")

weight = "eventweightlumi"

variable = "eventweightlumi"


signalfiles = [ "TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"]

backgroundfiles = [ 
	"QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
	"QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8.root", 
]

#backgroundfiles=signalfiles


signal = ROOT.TChain("Events")
background = ROOT.TChain("Events")

for file in backgroundfiles:
	background.Add(directory+file)


overflowmargin = 20.
maxVal = background.GetMaximum(variable)+overflowmargin
minVal = background.GetMinimum(variable)-overflowmargin

passcutHP = ROOT.TH1D("passHP", "passHP", 1, minVal, maxVal)
passcutLP = ROOT.TH1D("passLP", "passLP", 1, minVal, maxVal)
total = ROOT.TH1D("total", "total", 1, minVal, maxVal)
passcutHP.Sumw2()
passcutLP.Sumw2()
total.Sumw2()


background.Draw(variable+">>"+passcutHP.GetName(), weight*(basecut & cutHP))
background.Draw(variable+">>"+passcutLP.GetName(), weight*(basecut & cutLP))
background.Draw(variable+">>"+total.GetName(), weight*(basecut))


NpassHP = passcutHP.Integral()
NpassLP = passcutLP.Integral()
Ntotal = total.Integral()

print NpassHP, Ntotal 

effHP = passcutHP.Divide(total)
effErrHP = passcutLP.Divide(total)

#effLP = ROOT.getEfficiency(NpassLP, Ntotal)
#effErrLP = ROOT.getEfficiencyError(effLP, Ntotal)

print passcutHP.GetBinContent(1), passcutHP.GetBinError(1), passcutLP.GetBinContent(1), passcutLP.GetBinError(1)



