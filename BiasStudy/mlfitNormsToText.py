#run with: python mlfitNormsToText.py -d cards_EXP_TIGHT_em/ -o OUTPUT -m 900

#! /usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
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

from ROOT import TFile, TString, TH1F, TF1, TStyle, TCanvas, TPad, TGraphErrors, TFormula, TColor, kGreen, TLatex, kYellow, kBlue, kRed, TMath, gROOT

parser = OptionParser()
parser.add_option('-d', '--inputDirectory', action = "store",     type = "string", dest = "inputDirectory", default = "./")
parser.add_option('-o', '--outputDir',  action = 'store', type = "string" , default = "outputDir" , dest = "outputDir" )
parser.add_option('-m', '--mass', action = 'store', type = "string", default = "600", dest = "mass")

(options, args) = parser.parse_args()

def setPlotStyle():

 ROOT.gStyle.SetPadLeftMargin(0.10);
 ROOT.gStyle.SetPadRightMargin(0.05);
 ROOT.gStyle.SetOptTitle(0);
 ROOT.gStyle.SetOptStat(0);
 ROOT.gStyle.SetStatBorderSize(0);
 ROOT.gStyle.SetStatColor(10);
 ROOT.gStyle.SetStatFont(42);
 ROOT.gStyle.SetStatX(0.94);
 ROOT.gStyle.SetStatY(0.91);
 ROOT.gStyle.cd();
 ROOT.gStyle.SetGridStyle(2);                          

####MAIN PROGRAM
 
if __name__ == "__main__":

 ROOTStyle = os.getenv ("ROOTStyle");
 
 ROOT.gStyle.SetOptFit(111);
 ROOT.gStyle.SetOptStat(1111);
 ROOT.gStyle.SetStatFont(42);
 ROOT.gStyle.SetTitleSize(0.04,"XYZ");
    

###Search for all the files with a certain mass point, and fill the list temp
 
 nameInputDirectory = options.inputDirectory ;
 if not os.path.isdir(nameInputDirectory):
    print " bad input directory --> exit ";
    sys.exit();

 os.system("ls "+nameInputDirectory+" | grep root | grep mlfit | grep "+options.mass+" > list_temp.txt");

 vector_root_file = [];
 os.system("mkdir -p "+options.outputDir);

 with open("list_temp.txt") as input_list:
  for line in input_list:
   for name in line.split():
     name.replace(" ", "");
     vector_root_file.append(TFile(TString(nameInputDirectory+"/"+name).Data(),"READ"));
                              

 if len(argv) == 0: raise RuntimeError, "Usage: mlfitNormsToText.py [ -u ] mlfit.root";

 errors = True

###histogram creation

 histo_pull_bkg_B = ROOT.TH1F ("histo_pull_bkg_B", "histo_pull_bkg_B", 30, -5, 5);
 histo_pull_bkg_SB = ROOT.TH1F ("histo_pull_bkg_SB", "histo_pull_bkg_SB", 30, -5, 5);
 histo_pull_sig = ROOT.TH1F ("histo_pull_sig", "histo_pull_sig", 30, -5, 5);
 histo_residual_bkg_B = ROOT.TH1F ("histo_residual_bkg_B", "histo_residual_bkg_B", 30, -50, 50);
 histo_residual_bkg_SB = ROOT.TH1F ("histo_residual_bkg_SB", "histo_residual_bkg_SB", 30, -50, 50);
 histo_residual_sig = ROOT.TH1F ("histo_residual_sig", "histo_residual_sig", 30, -50, 50); 
 
 histo_sigma_bkg = ROOT.TH1F ("histo_sigma_bkg", "histo_sigma_bkg", 30, 0, 50);
 histo_sigma_sig = ROOT.TH1F ("histo_sigma_sig", "histo_sigma_sig", 30, 0, 50); 

 pull_bkg_B=[];
 pull_bkg_SB=[];
 pull_sig=[];


###start to loop on the files
 
 count=0;
 for ifile in range(len(vector_root_file)):

  print ifile;
    
#if len(argv) > 2 and argv[1] == "-u": 
#    errors = True
#    argv[1] = argv[2];
#  file = ROOT.TFile.Open(vector_root_file[ifile]);
  vector_root_file[ifile].cd();
  print vector_root_file[ifile].GetName();  
  prefit = vector_root_file[ifile].Get("norm_prefit")
  fit_s = vector_root_file[ifile].Get("norm_fit_s")
  fit_b = vector_root_file[ifile].Get("norm_fit_b")
  if prefit == None: stderr.write("Missing fit_s in %s. Did you run MaxLikelihoodFit in a recent-enough version of combine and with --saveNorm?\n" % vector_root_file[ifile]);
  if fit_s  == None:
     print "Missing fit_s in %s. Did you run MaxLikelihoodFit with --saveNorm?" % vector_root_file[ifile];
     continue;
 
  if fit_b  == None:
     print "Missing fit_b in %s. Did you run MaxLikelihoodFit with --saveNorm?" % vector_root_file[ifile];
     continue;
  count = count+1;

  bkg=0;
  bkg_pre=0;
  bkgSB=0;
  sig=0;
  sig_pre=0;
  sigma2_preS=0;
  sigma2_pre=0;

###Loop over the different backgrounds (or signals)
 
  iter = fit_s.createIterator()
  while True:
    norm_s = iter.Next()
    if norm_s == None: break;
    norm_b = fit_b.find(norm_s.GetName())
#    print norm_s.GetName();
    norm_p = prefit.find(norm_s.GetName()) if prefit else None
    m = re.match(r"(\w+)/(\w+)", norm_s.GetName());
    if m == None: m = re.match(r"n_exp_(?:final_)?(?:bin)+(\w+)_proc_(\w+)", norm_s.GetName());
    if m == None: raise RuntimeError, "Non-conforming object name %s" % norm_s.GetName()
    if norm_b == None: raise RuntimeError, "Missing normalization %s for background fit" % norm_s.GetName()    
    if prefit and norm_p and errors:
        print "%-30s %-30s %7.3f +/- %7.3f   %7.3f +/- %7.3f  %7.3f +/- %7.3f" % (m.group(1), m.group(2), norm_p.getVal(), norm_p.getError(), norm_s.getVal(), norm_s.getError(), norm_b.getVal(), norm_b.getError())
    else:
        if errors:
            print "%-30s %-30s %7.3f +/- %7.3f  %7.3f +/- %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_s.getError(), norm_b.getVal(), norm_b.getError())
        else:
            print "%-30s %-30s %7.3f %7.3f" % (m.group(1), m.group(2), norm_s.getVal(), norm_b.getVal())

###Fill the residual and pull vector/histograms

    if norm_s.GetName().find("ggH")==-1 and norm_s.GetName().find("vbfH")==-1:
        bkg = bkg + norm_b.getVal();
        bkgSB = bkgSB + norm_s.getVal();
        bkg_pre = bkg_pre + norm_p.getVal();
        sigma2_pre = sigma2_pre + norm_p.getError()*norm_p.getError();

    else:
        sig = sig + norm_s.getVal();
        sig_pre = sig_pre + norm_p.getVal();
        sigma2_preS = sigma2_preS + norm_p.getError()*norm_p.getError();        
        

  histo_residual_bkg_B.Fill(bkg-bkg_pre);
  histo_residual_bkg_SB.Fill(bkgSB-bkg_pre);
  histo_residual_sig.Fill(sig-sig_pre);
  
  sigma_pre = TMath.Sqrt(sigma2_pre);
  sigma_pre_S = TMath.Sqrt(sigma2_preS);

  pull_bkg_B.append(bkg-bkg_pre);
  pull_bkg_SB.append(bkgSB-bkg_pre);
  pull_sig.append(sig-sig_pre);

  histo_sigma_bkg.Fill(sigma_pre);
  histo_sigma_sig.Fill(sigma_pre_S);

#  histo_sigma.SetDirectory(0);
#  vector_root_file[ifile].Close();

 print count;
# print residual[ifile]
 for i in range(0,count):
    histo_pull_bkg_B.Fill(pull_bkg_B[i]/histo_sigma_bkg.GetMean());
    histo_pull_bkg_SB.Fill(pull_bkg_SB[i]/histo_sigma_bkg.GetMean());
    histo_pull_sig.Fill(pull_sig[i]/histo_sigma_sig.GetMean());
    
 print "Mean: ";
 print histo_residual_bkg_B.GetMean();
# gaus = ROOT.TF1("gaus","gaus",-3,3);
# histo_pull_bkg_B.Fit("gaus");

# setPlotStyle();


###Plot time
  
 c_pull_bkg_B = ROOT.TCanvas("pull_bkg_B","pull_bkg_B");
 c_pull_bkg_B.cd();
 histo_pull_bkg_B.Draw("");
# gaus.Draw("same");
 c_pull_bkg_B.SaveAs(options.outputDir+"/"+"pull_bkg_B_"+options.mass+".root","root");
# c_pull_bkg_B.Close();

 c_pull_bkg_SB = ROOT.TCanvas("pull_bkg_SB","pull_bkg_SB");
 c_pull_bkg_SB.cd();
 histo_pull_bkg_SB.Draw("");
# gaus.Draw("same");
 c_pull_bkg_SB.SaveAs(options.outputDir+"/"+"pull_bkg_SB_"+options.mass+".root","root");
# c_pull_bkg_B.Close();

 c_pull_sig = ROOT.TCanvas("pull_sig","pull_sig");
 c_pull_sig.cd();
 histo_pull_sig.Draw("");
# gaus.Draw("same");
 c_pull_sig.SaveAs(options.outputDir+"/"+"pull_sig_"+options.mass+".root","root");
# c_pull_sig.Close();

 c_residual_bkg_B = ROOT.TCanvas("residual_bkg_B","residual_bkg_B");
 c_residual_bkg_B.cd();
 histo_residual_bkg_B.Draw("");
# gaus.Draw("same");
 c_residual_bkg_B.SaveAs(options.outputDir+"/"+"residual_bkg_B_"+options.mass+".root","root");
# c_residual_bkg_B.Close();

 c_residual_bkg_SB = ROOT.TCanvas("residual_bkg_SB","residual_bkg_SB");
 c_residual_bkg_SB.cd();
 histo_residual_bkg_SB.Draw("");
# gaus.Draw("same");
 c_residual_bkg_SB.SaveAs(options.outputDir+"/"+"residual_bkg_SB_"+options.mass+".root","root");
# c_residual_bkg_B.Close();

 c_residual_sig = ROOT.TCanvas("residual_sig","residual_sig");
 c_residual_sig.cd();
 histo_residual_sig.Draw("");
# gaus.Draw("same");
 c_residual_sig.SaveAs(options.outputDir+"/"+"residual_sig_"+options.mass+".root","root");
# c_residual_sig.Close();



 c_sigma_bkg = ROOT.TCanvas("sigma_bkg","sigma_bkg",800,800);
 c_sigma_bkg.cd();
 histo_sigma_bkg.Draw();
 c_sigma_bkg.SaveAs(options.outputDir+"/"+"sigma_bkg_"+options.mass+".root","root");
# c_sigma_bkg.Close();

 c_sigma_sig = ROOT.TCanvas("sigma_sig","sigma_sig",800,800);
 c_sigma_sig.cd();
 histo_sigma_sig.Draw();
 c_sigma_sig.SaveAs(options.outputDir+"/"+"sigma_sig_"+options.mass+".root","root");
# c_sigma_sig.Close();

 os.system("rm list_temp.txt");
 
