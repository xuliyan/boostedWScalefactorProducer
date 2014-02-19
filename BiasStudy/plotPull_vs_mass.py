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

from ROOT import TFile, TString, TH1F, TF1, TStyle, TCanvas, TPad, TGraphErrors, TFormula

parser = OptionParser()
parser.add_option('-d', '--inputDirectory', action = "store",     type = "string", dest = "inputDirectory", default = "./")
parser.add_option('-b',                     action = 'store_true', dest='noX',     default = False, help = 'no X11 windows')
parser.add_option('-c', '--channel', action = 'store', type = "string", default = "em", dest = "channel" )
parser.add_option('-m', '--isMC',    action = 'store', type = "int" , default = 0 , dest = "isMC" )
parser.add_option('-o', '--outputDir',  action = 'store', type = "string" , default = "" , dest = "outputDir" )


(options, args) = parser.parse_args()

mass = [600,700,800,900,1000]

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
     
## main part of the code
if __name__ == "__main__":

 setPlotStyle();

 nameInputDirectory = options.inputDirectory ;
 if not os.path.isdir(nameInputDirectory):
    print " bad input directory --> exit ";
    sys.exit();


 os.system("ls "+nameInputDirectory+" | grep pull | grep root > list_temp.txt");

 vector_root_file = [];

 with open("list_temp.txt") as input_list:
  for line in input_list:
      for name in line.split():
        name.replace(" ", "");
        if nameInputDirectory == "./":
         vector_root_file.append(TFile(TString(name).Data(),"READ"));
        else:
         vector_root_file.append(TFile(TString(nameInputDirectory+"/"+name).Data(),"READ"));


  histogram_pull_vs_mass_nback = ROOT.TH1F("histogram_pull_vs_mass_nback","",len(mass),0,len(mass));
  gaussian_pull_vs_mass_nback  = ROOT.TH1F("gaussian_pull_vs_mass_nback","",len(mass),0,len(mass));
  
  histogram_pull_vs_mass_nsig  = ROOT.TH1F("histogram_pull_vs_mass_nsig","",len(mass),0,len(mass));
  gaussian_pull_vs_mass_nsig   = ROOT.TH1F("gaussian_pull_vs_mass_nsig","",len(mass),0,len(mass));

  graph = ROOT.TGraphErrors();

  for ifile in range(len(vector_root_file)):

   vector_root_file[ifile].cd();

   if options.isMC == 0: 

    histo_pull_data_wjet    = vector_root_file[ifile].Get("rrv_number_data_sb_lo_fit_%s_mlvj_data_pull"%(options.channel));
#    gaussian_pull_data_wjet = vector_root_file[ifile].Get("Gaussian_rrv_number_data_sb_lo_fit_%s_mlvj_data_pull"%(options.channel)); 
    gaussian_pull_data_wjet = vector_root_file[ifile].Get("Gaussian_pull_1"); 
    histo_pull_signal       = vector_root_file[ifile].Get("rrv_number_signal_region_fit_H_data_pull");
#    gaussian_pull_signal    = vector_root_file[ifile].Get("Gaussian_rrv_number_signal_region_fit_H_data_pull"); 
    gaussian_pull_signal    = vector_root_file[ifile].Get("Gaussian_pull_2"); 
   else:  
    histo_pull_data_wjet    = vector_root_file[ifile].Get("rrv_number_WJets0_sb_lo_fit_%s_mlvj_data_pull"%(options.channel));
    gaussian_pull_data_wjet = vector_root_file[ifile].Get("Gaussian_pull_1"); 
#    gaussian_pull_signal = vector_root_file[ifile].Get("Gaussian_rrv_number_data_sb_lo_fit_%s_mlvj_data_pull"%(options.channel)); 
    histo_pull_signal       = vector_root_file[ifile].Get("rrv_number_signal_region_fit_H_data_pull");
#    gaussian_pull_signal    = vector_root_file[ifile].Get("Gaussian_rrv_number_signal_region_fit_H_data_pull"); 
    gaussian_pull_signal    = vector_root_file[ifile].Get("Gaussian_pull_2"); 

   masspoint = 0 ;
   for imass in range(len(mass)):
      if TString(vector_root_file[ifile].GetName()).Contains("%s"%(mass[imass])):
         masspoint = mass[imass];
         break;
            
   histogram_pull_vs_mass_nback.SetBinContent(imass+1,histo_pull_data_wjet.GetMean());
   histogram_pull_vs_mass_nback.SetBinError(imass+1,histo_pull_data_wjet.GetRMS());
   histogram_pull_vs_mass_nback.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));

   histogram_pull_vs_mass_nsig.SetBinContent(imass+1,histo_pull_signal.GetMean());
   histogram_pull_vs_mass_nsig.SetBinError(imass+1,histo_pull_signal.GetRMS());
   histogram_pull_vs_mass_nsig.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));

   gaussian_pull_vs_mass_nback.SetBinContent(imass+1,gaussian_pull_data_wjet.GetParameter(1));
   gaussian_pull_vs_mass_nback.SetBinError(imass+1,gaussian_pull_data_wjet.GetParameter(2));
   gaussian_pull_vs_mass_nback.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));

   gaussian_pull_vs_mass_nsig.SetBinContent(imass+1,gaussian_pull_signal.GetParameter(1));
   gaussian_pull_vs_mass_nsig.SetBinError(imass+1,gaussian_pull_signal.GetParameter(2));
   gaussian_pull_vs_mass_nsig.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));

   graph.SetPoint(imass+1,imass,0);
   graph.SetPointError(imass+1,0,1);


  graph.SetPoint(len(mass)+1,len(mass),0);
  graph.SetPointError(len(mass)+1,0,1);
           
 canvas_1 = ROOT.TCanvas("histo_bkg", "histo_bkg")
 histogram_pull_vs_mass_nback.GetYaxis().SetTitle("Pull");
 histogram_pull_vs_mass_nback.GetYaxis().SetTitleOffset(1.05);
 histogram_pull_vs_mass_nback.GetXaxis().SetLabelSize(0.065);
 histogram_pull_vs_mass_nback.GetYaxis().SetLabelSize(0.045);
 histogram_pull_vs_mass_nback.GetYaxis().SetTitleSize(0.045);
 histogram_pull_vs_mass_nback.SetTitle("Pull of Background Number of Events Binned")
 histogram_pull_vs_mass_nback.SetMarkerStyle(20);
 histogram_pull_vs_mass_nback.SetMarkerSize(1.0);
 histogram_pull_vs_mass_nback.SetLineColor(1);
 histogram_pull_vs_mass_nback.SetLineWidth(2);
 graph.SetFillColor(3);
 graph.SetFillStyle(3001);


 oneLine = ROOT.TF1("oneLine","0",0,100);
 oneLine.SetLineColor(ROOT.kRed);
 oneLine.SetLineWidth(2);
                     
 histogram_pull_vs_mass_nback.SetMaximum(1.5);
 histogram_pull_vs_mass_nback.SetMinimum(-1.5);
 histogram_pull_vs_mass_nback.Draw("pe");
 graph.Draw("e3same");
 oneLine.Draw("Lsame");
 histogram_pull_vs_mass_nback.Draw("pesame");

 canvas_1.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nback.GetName()+".png","png");
 canvas_1.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nback.GetName()+".pdf","pdf");


 canvas_2 = ROOT.TCanvas("histo_sig", "histo_sig")
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitle("Pull");
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleOffset(1.05);
 gaussian_pull_vs_mass_nsig.GetXaxis().SetLabelSize(0.065);
 gaussian_pull_vs_mass_nsig.GetYaxis().SetLabelSize(0.045);
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleSize(0.045);
 gaussian_pull_vs_mass_nsig.SetTitle("Pull of Signal Number of Events Binned")
 gaussian_pull_vs_mass_nsig.SetMarkerStyle(20);
 gaussian_pull_vs_mass_nsig.SetMarkerSize(1.0);
 gaussian_pull_vs_mass_nsig.SetLineColor(1);
 gaussian_pull_vs_mass_nsig.SetLineWidth(2);
 graph.SetFillColor(3);
 graph.SetFillStyle(3001);

 oneLine = ROOT.TF1("oneLine","0",0,100);
 oneLine.SetLineColor(ROOT.kRed);
 oneLine.SetLineWidth(2);
                     
 gaussian_pull_vs_mass_nsig.SetMaximum(1.5);
 gaussian_pull_vs_mass_nsig.SetMinimum(-1.5);
 gaussian_pull_vs_mass_nsig.Draw("pe");
 graph.Draw("e3same");
 oneLine.Draw("Lsame");
 gaussian_pull_vs_mass_nsig.Draw("pesame");

 canvas_2.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+".png","png");
 canvas_2.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+".pdf","pdf");


 canvas_3 = ROOT.TCanvas("gaus_bkg", "gaus_bkg")
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitle("Pull");
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitleOffset(1.05);
 gaussian_pull_vs_mass_nback.GetXaxis().SetLabelSize(0.065);
 gaussian_pull_vs_mass_nback.GetYaxis().SetLabelSize(0.045);
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitleSize(0.045);
 gaussian_pull_vs_mass_nback.SetTitle("Pull of Background Number of Events Gaussian")
 gaussian_pull_vs_mass_nback.SetMarkerStyle(20);
 gaussian_pull_vs_mass_nback.SetMarkerSize(1.0);
 gaussian_pull_vs_mass_nback.SetLineColor(1);
 gaussian_pull_vs_mass_nback.SetLineWidth(2);
 graph.SetFillColor(3);
 graph.SetFillStyle(3001);


 oneLine = ROOT.TF1("oneLine","0",0,100);
 oneLine.SetLineColor(ROOT.kRed);
 oneLine.SetLineWidth(2);
                     
 gaussian_pull_vs_mass_nback.SetMaximum(1.5);
 gaussian_pull_vs_mass_nback.SetMinimum(-1.5);
 gaussian_pull_vs_mass_nback.Draw("pe");
 graph.Draw("e3same");
 oneLine.Draw("Lsame");
 gaussian_pull_vs_mass_nback.Draw("pesame");

 canvas_3.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nback.GetName()+".png","png");
 canvas_3.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nback.GetName()+".pdf","pdf");


 canvas_4 = ROOT.TCanvas("gaus_sig", "gaus_sig")
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitle("Pull");
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleOffset(1.05);
 gaussian_pull_vs_mass_nsig.GetXaxis().SetLabelSize(0.065);
 gaussian_pull_vs_mass_nsig.GetYaxis().SetLabelSize(0.045);
 gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleSize(0.045);
 gaussian_pull_vs_mass_nsig.SetTitle("Pull of Signal Number of Events Gaussian")
 gaussian_pull_vs_mass_nsig.SetMarkerStyle(20);
 gaussian_pull_vs_mass_nsig.SetMarkerSize(1.0);
 gaussian_pull_vs_mass_nsig.SetLineColor(1);
 gaussian_pull_vs_mass_nsig.SetLineWidth(2);
 graph.SetFillColor(3);
 graph.SetFillStyle(3001);

 gaussian_pull_vs_mass_nsig.SetMaximum(1.5);
 gaussian_pull_vs_mass_nsig.SetMinimum(-1.5);
 gaussian_pull_vs_mass_nsig.Draw("pe");
 graph.Draw("e3same");
 oneLine.Draw("Lsame");
 gaussian_pull_vs_mass_nsig.Draw("pesame");

 canvas_4.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+".png","png");
 canvas_4.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+".pdf","pdf");



 os.system("rm list_temp.txt");
