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
parser.add_option('-d', '--inputDirectory', action = "store",     type = "string", dest = "inputDirectory", default = "./")
parser.add_option('-c', '--channel', action = 'store', type = "string", default = "em", dest = "channel" )
parser.add_option('-j', '--jetBin',  action = 'store', type = "string", default = "", dest = "jetBin" )
parser.add_option('--skipJetSystematics',      action="store", type="int", dest="skipJetSystematics",  default=0)
parser.add_option('--pseudodata',  action='store', type="int",    dest='pseudodata',  default=1, help='pseudodata 0 -> use real data, else use stack of MC backgrounds')

(options, args) = parser.parse_args()

####### main part of the code ##############
if __name__ == "__main__":

 ## set the working directory where datacards are located   
 os.chdir(options.inputDirectory);
 
 if options.channel == "mu":

  os.system("sed -i \"#CMS_hwwlvj_TTbar   lnN   -         -       -      1.070    -       -/d\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_trigger_mu/CMS_trigger_m/g\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_Top_norm_mu/CMS_Top_norm_m/g\" *_mu_*.txt");
  os.system("sed -i \"s/Wjet_Norm_mu/Wjet_Norm_m/g\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_eff_mu/CMS_eff_m/g\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_btag_eff/CMS_eff_b/g\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_scale_l/CMS_scale_m/g\" *_mu_*.txt");
  os.system("sed -i \"s/CMS_res_l/CMS_res_m/g\" *_mu_*.txt");

  if options.skipJetSystematics :
      
   if options.pseudodata :   
    os.system("ls | grep _mu_ | grep -v _ggH_ | grep -v _vbfH_ | grep unbin | awk '{print \"sed -i \"s/Wjet_Norm_m   lnN     -         -    1.044     -      -        -/Wjet_Norm_m   lnN     -         -    1.047     -      -        -\\nCMS_scale_j lnN   1.029     1.044     0.975\\/1.035    1.050\\/0.950   1.051\\/0.949   1.032\\/0.968\\nCMS_res_j lnN   1.008     1.035     0.992\\/1.008    1.004\\/0.996   1.009\\/0.991   1.007\\/0.993/g\" \" $1}' | /bin/sh");
    os.system("ls | grep _mu_ | grep -v _ggH_ | grep -v _vbfH_ | grep unbin | awk '{print \"sed -i \"s/Wjet_Norm_m   lnN     -         -    1.047     -      -        -/Wjet_Norm_m   lnN     -         -    1.047     -      -        -\n CMS_scale_j lnN   1.029     1.044     0.975\\/1.035    1.050\\/0.950   1.051\\/0.949   1.032\\/0.968\\nCMS_res_j lnN   1.008     1.035     0.992\\/1.008    1.004\\/0.996   1.009\\/0.991   1.007\\/0.993/g\" \" $1}' | /bin/sh");

 if options.channel == "el":

  os.system("sed -i \"#CMS_hwwlvj_TTbar   lnN   -         -       -      1.070    -       -/d\" *_el_*.txt");
  os.system("sed -i \"s/CMS_trigger_el/CMS_trigger_e/g\" *_el_*.txt");
  os.system("sed -i \"s/CMS_Top_norm_el/CMS_Top_norm_e/g\" *_el_*.txt");
  os.system("sed -i \"s/Wjet_Norm_el/Wjet_Norm_e/g\" *_el_*.txt");
  os.system("sed -i \"s/CMS_eff_el/CMS_eff_e/g\" *_el_*.txt");
  os.system("sed -i \"s/CMS_btag_eff/CMS_eff_b/g\" *_el_*.txt");
  os.system("sed -i \"s/CMS_scale_l/CMS_scale_e/g\" *_el_*.txt");
  os.system("sed -i \"s/CMS_res_l/CMS_res_e/g\" *_el_*.txt");

  if options.skipJetSystematics :

   if options.pseudodata :   
    os.system("ls | grep _el_ | grep -v _ggH_ | grep -v _vbfH_ | grep unbin | awk '{print \"sed -i \"s/Wjet_Norm_e   lnN     -         -    1.096     -      -        -/Wjet_Norm_e   lnN     -         -    1.090     -      -        -\\nCMS_scale_j lnN   1.034     1.072     0.971\\/1.029    1.052\\/0.948   1.054\\/0.946   1.037\\/0.963\\nCMS_res_j lnN   1.007     1.030     0.985\\/1.015    1.004\\/0.996   1.014\\/0.986   1.004\\/0.996/g\" \" $1}' | /bin/sh");
    
    os.system("ls | grep _el_ | grep -v _ggH_ | grep -v _vbfH_ | grep unbin | awk '{print \"sed -i \"s/Wjet_Norm_e   lnN     -         -    1.090     -      -        -/Wjet_Norm_e   lnN     -         -    1.090     -      -        -\\nCMS_scale_j lnN   1.034     1.072     0.971\\/1.029    1.052\\/0.948   1.054\\/0.946   1.037\\/0.963\\nCMS_res_j lnN   1.007     1.030     0.985\\/1.015    1.004\\/0.996   1.014\\/0.986   1.004\\/0.996/g\" \" $1}' | /bin/sh");

 if options.channel == "em" and options.jetBin == "_2jet" :

  os.system("sed -i \"#CMS_hwwlvj_TTbar    lnN    -          -        -       1.070    -       -       -/d\" *_em_2jet_*.txt");
  os.system("sed -i \"s/CMS_trigger_em/CMS_trigger_l/g\" *_em_2jet_*.txt");
  os.system("sed -i \"s/CMS_Top_norm_em/CMS_Top_norm_l_2jet/g\" *_em_2jet_*.txt");
  os.system("sed -i \"s/Wjet_Norm_em__2jet/Wjet_Norm_l_2jet/g\" *_em_2jet_*.txt");
  os.system("sed -i \"s/CMS_eff_em/CMS_eff_l/g\" *_em_2jet_*.txt");
  os.system("sed -i \"s/CMS_btag_eff/CMS_eff_b/g\" *_em_2jet_*.txt");

  if options.skipJetSystematics :

   if options.pseudodata :   
    os.system("ls | grep _em_2jet_ | grep -v _ggH_ | grep -v _vbfH_ | grep unbin | awk '{print \"sed -i \"s/Wjet_Norm_l_2jet lnN     -           -         1.283      -       -        -        -/Wjet_Norm_l_2jet lnN     -           -         1.283      -       -        -        -\\nCMS_scale_j lnN   1.065     1.043     0.959/1.041    1.071/0.929   1.071/0.929   1.075/0.925    1.046/0.954\\nCMS_res_j lnN   1.070     1.141     0.981/1.019    1.051/0.949   1.084/0.916   1.043/0.957    1.063/0.937/g\" \" $1}' | /bin/sh");
   
   
