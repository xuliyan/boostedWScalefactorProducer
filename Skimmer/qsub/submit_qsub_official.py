#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re
import datetime
import subprocess
import itertools

now = datetime.datetime.now()
timestamp =  now.strftime("%Y_%m_%d")

def split_seq(iterable, size):
    it = iter(iterable)
    item = list(itertools.islice(it, size))
    while item:
        yield item
        item = list(itertools.islice(it, size))
        
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
   
def createJobs(f, outfolder,name,nchunks):
  infiles = []
  for files in f:
    infiles.append("root://cms-xrd-global.cern.ch/"+files)
  cmd = 'python qsub_script_SFs.py %s %s %s %i \n'%(','.join(infiles), outfolder,name,nchunks)
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
	out = "Skimmed_%s/"%timestamp
	batchSystem = 'psibatch_runner.sh'
	
	patternsTT    = ["/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER"]
	patternsdata  = ["/SingleMuon/arizzi-RunII2017ReReco17Nov17-94X-Nano01-e70630e8aef2c186cd650f6150c31168/USER"]
	patternsST    = ["/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/WW_TuneCP5_13TeV-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER"]
	patternsVV    = ["/WZ_TuneCP5_13TeV-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/ZZ_TuneCP5_13TeV-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/W1JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER"]
	patternsWJets = ["/W1JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/W2JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/W2JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER","/W3JetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/arizzi-RunIIFall17MiniAOD-94X-Nano01Fall17-e273b12d9f89d622a34e4bc98b05ee29/USER"]
	
	patternsTT    = ["/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"]
	patternsST    = ["/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM""/ST_t-channel_antitop_5f_TuneCP5_PSweights_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM""/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"]
	patternsVV    = ["/ZZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM","/WZ_TuneCP5_13TeV-pythia8/RunIIFall17NanoAOD-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/NANOAODSIM"]
	patternsWJets = []
	
	if len(sys.argv) > 1:
	  if sys.argv[1].find("TT")!=-1:		patterns = patternsTT    
	  if sys.argv[1].find("data")!=-1:  patterns = patternsdata  
	  if sys.argv[1].find("ST")!=-1:		patterns = patternsST    
	  if sys.argv[1].find("VV")!=-1:		patterns = patternsVV    
	  if sys.argv[1].find("WJets")!=-1: patterns = patternsWJets 
	  if sys.argv[1].find("ALL")!=-1:   patterns = patternsTT+patternsdata+patternsST+patternsVV+patternsWJets  
	
	  print 'Location of input files',  patterns
	else:
		print "No location given, give folder with files"
		exit(0)
	
	if len(sys.argv) > 2:
		out = sys.argv[2]
		print 'Output goes here: ', out
	else:
		print "Using default output folder: ", out
	
	try: os.stat(out)
	except: os.mkdir(out)
	
	for pattern in patterns:
		files = getFileListDAS(pattern)
		print "FILELIST = ", files
		name = pattern.split("/")[1].replace("/","")
		print "creating job file " ,'joblist%s.txt'%name
		jobList = 'joblist%s.txt'%name
		jobs = open(jobList, 'w')
		nChunks = 0
		outfolder = out+name
		try: os.stat(outfolder)
		except: os.mkdir(outfolder)
		try: os.stat(outfolder+'/logs/')
		except: os.mkdir(outfolder+'/logs/')
		filelists = list(split_seq(files,10))	
		for f in filelists:
			print "FILES = ",f
			createJobs(f,outfolder,name,nChunks)
			nChunks = nChunks+1
		
		jobs.close()
    # submit = raw_input("Do you also want to submit the jobs to the batch system? [y/n] ")
		submit = 'y'
		if submit == 'y' or submit=='Y':
			submitJobs(jobList,nChunks, outfolder, batchSystem)
		else:
			print "Not submitting jobs"
		
		
		
	