#!/usr/bin/env python

import ROOT
from rootpy.tree import Cut

ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")


directory = "/eos/home-m/mhuwiler/data/Wtagging/Mergedefinition2017/"

basecut = Cut("GenJetAK8_pt>200")
cut = Cut("SelectedJet_tau21<0.35")

weight = "eventweightlumi"

variable = "puweight"


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


signal = ROOT.TChain("Events")
background = ROOT.TChain("Events")

for file in backgroundfiles:
	background.Add(directory+file)


overflowmargin = 20
maxVal = background.GetMaximum(variable)-overflowmargin
minVal = background.GetMaximum(variable)+overflowmargin

passcut = ROOT.TH1D("pass", "pass", 1, minVal, maxVal)
total = ROOT.TH1D("total", "total", 1, minVal, maxVal)


background.Draw(variable+">>"+passcut.GetName(), weight*(basecut & cut))
background.Draw(variable+">>"+total.GetName(), weight*(basecut & -cut))


Npass = passcut.Integral()
Ntotal = total.Integral()

print Npass, Ntotal 

eff = ROOT.getEfficiency(Npass, Ntotal)
effErr = ROOT.getEfficiencyError(eff, Ntotal)

print eff, effErr



