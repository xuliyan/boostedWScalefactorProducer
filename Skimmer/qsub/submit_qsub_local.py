#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re
import datetime
import subprocess

now = datetime.datetime.now()
timestamp =  now.strftime("%Y_%m_%d")

def getFileListDAS(dataset,instance="prod/phys03",run=-1):
	cmd='das_client --limit=0 --query="file dataset=%s instance=%s"'%(dataset,instance)
	print "Executing ",cmd
	cmd_out = getoutput( cmd )
	tmpList = cmd_out.split(os.linesep)
	files = []
	for l in tmpList:
	   if l.find(".root") != -1:
	      files.append(l)
	         
	return files 
   
def createJobs(f, outfolder,name):
    cmd = 'python qsub_script_SFs.py root://t3dcachedb.psi.ch:1094//%s %s %s \n'%(f, outfolder,name)
    print cmd
    jobs.write(cmd)
    return 1

def submitJobs(jobList, nchunks, outfolder, batchSystem):
    print 'Reading joblist'
    jobListName = jobList
    print jobList
#    subCmd = 'qsub -t 1-%s -o logs nafbatch_runner_GEN.sh %s' %(nchunks,jobListName)
    subCmd = 'qsub -t 1-%s -o %s/logs/ %s %s' %(nchunks,outfolder,batchSystem,jobListName)
    print 'Going to submit', nchunks, 'jobs with', subCmd
    os.system(subCmd)

    return 1


if __name__ == "__main__":
  location = '/pnfs/psi.ch/cms/trivcat/store/user/thaarres/'
  outfolder = "%s"%timestamp
  batchSystem = 'psibatch_runner.sh'
  if len(sys.argv) > 1:
    if sys.argv[1].find("TT80X")!=-1: 	patterns = ["TT_TuneCUETP8M2T4_13TeV-powheg-pythia8"]	
    elif sys.argv[1].find("TT")!=-1:	  patterns = ["TTToSemiLeptonic_TuneCP5_PSweights_13TeV-powheg-pythia8"]	
    elif sys.argv[1].find("ST")!=-1:    patterns = ["ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8","ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8","ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8","ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8", "ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8"]
    elif sys.argv[1].find("VV")!=-1:    patterns = ["WW_TuneCP5_13TeV-pythia8","WZ_TuneCP5_13TeV-pythia8","ZZ_TuneCP5_13TeV-pythia8" ]
    elif sys.argv[1].find("Wjets")!=-1: patterns = [""]
    elif sys.argv[1].find("data")!=-1:  patterns = ["SingleMuon/SingleMuon_Run2017B-17Nov2017-v1","SingleMuon/SingleMuon_Run2017C-17Nov2017-v1","SingleMuon/SingleMuon_Run2017D-17Nov2017-v1","SingleMuon/SingleMuon_Run2017E-17Nov2017-v1","SingleMuon/SingleMuon_Run2017F-17Nov2017-v1"]
    
    
    
    print 'Location of input files', patterns
  else:
    print "No location given, give folder with files"
    exit(0)
	
  if len(sys.argv) > 2:
		outfolder = sys.argv[2]
		print 'Output goes here: ', outfolder
  else:
  	print "Using default output folder(s): ", [p for p in patterns]
  
  
  
  for pattern in patterns:
    print "Looking for file in " ,location+'/'+pattern+'/*/*/*/*.root'
    if sys.argv[1].find("data")!=-1: files = glob.glob(location+'/'+pattern+'*/*/*/*.root')
    else: files = glob.glob(location+'/'+pattern+'/*/*/*/*.root')
    print "FILELIST = ", files
    name = pattern
    if sys.argv[1].find("data")!=-1: name = pattern.split("/")[1].replace("/","")
    print "creating job file " ,'joblist%s.txt'%name
    jobList = 'joblist%s.txt'%name
    jobs = open(jobList, 'w')
    nChunks = 0
    outfolder = name
    try: os.stat(outfolder)
    except: os.mkdir(outfolder)
    try: os.stat(outfolder+'/logs/')
    except: os.mkdir(outfolder+'/logs/')
    for f in files:
    	createJobs(f,outfolder,name)
    	nChunks = nChunks+1
    
    jobs.close()
    submit = raw_input("Do you also want to submit the jobs to the batch system? [y/n] ")
    if submit == 'y' or submit=='Y':
    	submitJobs(jobList,nChunks, outfolder, batchSystem)
    else:
    	print "Not submitting jobs"
    
  
  
  