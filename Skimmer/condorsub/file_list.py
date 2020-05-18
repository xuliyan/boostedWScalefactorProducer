#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re
import datetime
import subprocess
import itertools

now = datetime.datetime.now()
timestamp =  now.strftime("%Y_%m_%d")
nfilesperjob = 1

listDirectory = "fileLists_Fall17_nanoAODv6/"

def createLists(dataset, name):
    instance="prod/phys03"
    cmd='das_client --limit=0 --query="file dataset=%s"'%(dataset)
    print "Executing ",cmd
    cmd_out = getoutput( cmd )
    tmpList = cmd_out.split(os.linesep)
    files = []
    for l in tmpList:
        if l.find(".root") != -1:
            files.append(l)
         
    fileName = listDirectory+"/"+name+".txt"
    with open(fileName, "w") as f:
        for l in files:
            f.write("%s\n" % l)
    print "Wrote file list", fileName
    return

if __name__ == "__main__":
    out = "Skimmed_%s/"%timestamp
    createlists = False

    patternsData  = [
        "/SingleMuon/Run2017B-Nano25Oct2019-v1/NANOAOD",
        "/SingleMuon/Run2017C-Nano25Oct2019-v1/NANOAOD",
        "/SingleMuon/Run2017D-Nano25Oct2019-v1/NANOAOD",
        "/SingleMuon/Run2017E-Nano25Oct2019-v1/NANOAOD",
        "/SingleMuon/Run2017F-Nano25Oct2019-v1/NANOAOD",
        "/SingleElectron/Run2017B-Nano25Oct2019-v1/NANOAOD",
        "/SingleElectron/Run2017C-Nano25Oct2019-v1/NANOAOD",
        "/SingleElectron/Run2017D-Nano25Oct2019-v1/NANOAOD",
        "/SingleElectron/Run2017E-Nano25Oct2019-v1/NANOAOD",
        "/SingleElectron/Run2017F-Nano25Oct2019-v1/NANOAOD",
    ]
    patternsTT    = [
        "/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7_ext1-v1/NANOAODSIM",
        "/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_new_pmx_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/TT_TuneCH3_13TeV-powheg-herwig7/RunIIFall17NanoAODv5-PU2017_12Apr2018_Nano1June2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_new_pmx_102X_mc2017_realistic_v7-v1/NANOAODSIM",
    ]
    patternsST    = [
        "/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_new_pmx_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
    ]
    patternsVV    = [
        "/WZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WW_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v2/NANOAODSIM",
    ]
    patternsWJets = [
        "/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
    ]
    patternsQCD = [
        "/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_new_pmx_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
        "/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17NanoAODv6-PU2017_12Apr2018_Nano25Oct2019_102X_mc2017_realistic_v7-v1/NANOAODSIM",
    ]
    
    patterns = []

    if len(sys.argv) > 1:
        if sys.argv[1].find("data")!=-1:  patterns = patternsData  
        if sys.argv[1].find("TT")!=-1:    patterns = patternsTT    
        if sys.argv[1].find("ST")!=-1:    patterns = patternsST    
        if sys.argv[1].find("VV")!=-1:    patterns = patternsVV    
        if sys.argv[1].find("WJets")!=-1: patterns = patternsWJets 
        if sys.argv[1].find("QCD")!=-1:   patterns = patternsQCD
        if sys.argv[1].find("mc")!=-1:    patterns = patternsTT+patternsST+patternsVV+patternsWJets+patternsQCD
        if sys.argv[1].find("all")!=-1:   patterns = patternsTT+patternsST+patternsVV+patternsWJets+patternsQCD+patternsData
        print 'Location of input files',  patterns
    else:
        print "No location given, give folder with files"
        exit(0)
        
    if len(sys.argv) > 2:
        if sys.argv[2].find("-c")!=-1: createlists = True
        else:
            out = sys.argv[2]
            print 'Output goes here: ', out
    else:
        print "Using default output folder: ", out
      
    try: os.stat(out)
    except: os.mkdir(out)
  
    for pattern in patterns:
        name = pattern.split("/")[1].replace("/","") + ("-" + pattern.split("/")[2].split("-")[0] if 'Run201' in pattern else "")
        createLists(pattern, name)
  
