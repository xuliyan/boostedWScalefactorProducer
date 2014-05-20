#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess
from subprocess import Popen
from optparse import OptionParser

from ROOT import TFile, TString, TH1F, TF1, TStyle, TCanvas, TPad, TGraphErrors, TFormula, TColor, kGreen, TLatex, kYellow, kBlue, kRed, gROOT

##########################
### Parse the options ####
##########################

parser = OptionParser()
parser.add_option('-d', '--inputDirectory', action = "store", type = "string", dest = "inputDirectory", default = "./")
parser.add_option('-c', '--channel',        action = 'store', type = "string", dest = "channel", default = "em")
parser.add_option('-j', '--jetBin',         action = 'store', type = "string", default = "", dest = "jetBin" )

(options, args) = parser.parse_args()

mass   = [600,700,800,900,1000]
cprime = [01,02,03,05,07,10]
brNew  = [00,01,02,03,04,05]

if __name__ == "__main__":

  ## set the working directory where datacards are located
  os.chdir(options.inputDirectory);
     

  for imass in mass:
    for icprime in cprime :
        for ibrnew in brNew :

            fileName = "higgsCombinehwwlvj_pval_exp_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(imass,options.channel,options.jetBin,icprime,ibrnew,imass);
            command  = "hadd -f %s higgsCombinehwwlvj_pval_exp_ggH%03d_%s%s_%02d_%02d_unbin*"%(fileName,imass,options.channel,options.jetBin,icprime,ibrnew);
            print "command ",command ;
            os.system(command);

          
