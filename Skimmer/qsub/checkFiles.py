#!/usr/bin/env python
import os, glob, sys
from commands import getoutput
import re
import ROOT

pattern = "Wprime"
if len(sys.argv) > 1: pattern = sys.argv[1]
		
filelist = glob.glob(pattern+'*.root')
# filelist = glob.glob('/scratch/thaarres/SUBSTRUCTURE/LOLAoutput/*.root')

for file in filelist:
	f = ROOT.TFile(file, "read")
	if not f.GetListOfKeys().Contains("Runs"):
  # if not f.GetListOfKeys().Contains("tree"):
	    print "FILE IS BUGGY. WILL REMOVE!"
	    cmd = 'rm %s' %file
	    print 'Going to execute: ' , cmd
	    print "Remember to resubmit %s , job number %s" %(pattern,file.split("_")[2])
	    os.system(cmd)
	else:
		continue
		print "FILE IS GOOD, KEEPING IT"	