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
parser.add_option('-b',                     action = 'store_true', dest='noX',     default = False, help = 'no X11 windows')
parser.add_option('-c', '--channel', action = 'store', type = "string", default = "em", dest = "channel" )
parser.add_option('-m', '--isMC',    action = 'store', type = "int" , default = 0 , dest = "isMC" )
parser.add_option('-o', '--outputDir',  action = 'store', type = "string" , default = "" , dest = "outputDir" )
parser.add_option('-g','--fgen',      help='function to generate toys Exp,ExpTail,Pow2,ExpN)', type="string", default="ExpN")
parser.add_option('-r','--fres',      help='function to fit toys (Exp,ExpTail,Pow2,ExpN)',     type="string", default="ExpN")
parser.add_option('-f','--fitjetmass',   help='flag in order to specify if the fit is done on the jet mass', type="int", default=0)
parser.add_option('-t','--ttbarcontrolregion',   help='in order to specify if the study is performed in the ttbar control region', type="int", default=0)
parser.add_option('-j','--mlvjregion',   help='in order to specify if the mWW region where the fit is performed', type="string", default="_sb_lo")
parser.add_option('-k','--onlybackgroundfit',   help='in order to specify if an only background fit has been perfomed', type="int", default=0)

(options, args) = parser.parse_args()

####### mass vector #############
mass = [600]
mass = [600,700,800,900,1000]

####### plot style ##############
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
     
####### main part of the code ##############
if __name__ == "__main__":

 ## take the list of TFile to be used for the plot

 ROOTStyle = os.getenv ("ROOTStyle");
 gROOT.ProcessLine((".x "+ROOTStyle+"/rootLogon.C"));
 gROOT.ProcessLine((".x "+ROOTStyle+"/rootPalette.C"));
 gROOT.ProcessLine((".x "+ROOTStyle+"/rootColors.C"));
 gROOT.ProcessLine((".x "+ROOTStyle+"/setTDRStyle.C"));
                
 ROOT.gStyle.SetOptFit(111);
 ROOT.gStyle.SetOptStat(1111);
 ROOT.gStyle.SetStatFont(42);
 ROOT.gStyle.SetTitleSize(0.04,"XYZ");
 
 nameInputDirectory = options.inputDirectory ;
 if not os.path.isdir(nameInputDirectory):
    print " bad input directory --> exit ";
    sys.exit();

 if options.fitjetmass :
   spectrum = "_mj";
   options.mlvjregion = "";
 else:
   spectrum = "_mlvj";
   
 if options.fitjetmass and options.ttbarcontrolregion:
  os.system("ls "+nameInputDirectory+" | grep pull | grep root | grep "+options.fgen+"_"+options.fres+"  | grep jetmass | grep ttbar > list_temp.txt");
 elif options.fitjetmass and not options.ttbarcontrolregion:
  os.system("ls "+nameInputDirectory+" | grep pull | grep root | grep "+options.fgen+"_"+options.fres+"  | grep jetmass | grep -v ttbar > list_temp.txt");
 elif not options.fitjetmass and options.ttbarcontrolregion:
  os.system("ls "+nameInputDirectory+" | grep pull | grep root | grep "+options.fgen+"_"+options.fres+"  | grep -v jetmass | grep ttbar > list_temp.txt");
 elif not options.fitjetmass and not options.ttbarcontrolregion:
  os.system("ls "+nameInputDirectory+" | grep pull | grep root | grep "+options.fgen+"_"+options.fres+"  | grep -v jetmass | grep -v ttbar > list_temp.txt");

 vector_root_file = [];

 os.system("mkdir -p "+options.outputDir);

 with open("list_temp.txt") as input_list:
  for line in input_list:
      for name in line.split():
        name.replace(" ", "");
        vector_root_file.append(TFile(TString(nameInputDirectory+"/"+name).Data(),"READ"));


  ## make the histograms

  histogram_pull_vs_mass_nback = ROOT.TH1F("histogram_pull_vs_mass_nback","",len(mass),0,len(mass));
  gaussian_pull_vs_mass_nback  = ROOT.TH1F("gaussian_pull_vs_mass_nback","",len(mass),0,len(mass));

  if not options.onlybackgroundfit:
  
   histogram_pull_vs_mass_nsig  = ROOT.TH1F("histogram_pull_vs_mass_nsig","",len(mass),0,len(mass));
   gaussian_pull_vs_mass_nsig   = ROOT.TH1F("gaussian_pull_vs_mass_nsig","",len(mass),0,len(mass));

  graph_1 = ROOT.TGraphErrors();
  graph_2 = ROOT.TGraphErrors();

  canvas_pull_bkg_data = [];
  canvas_pull_sig_data = [];
  
  ## loop and fill the histograms
  for imass in range(len(mass)):

   masspoint = -1 ;
   ifilePos  = -1 ;
   
   for ifile in range(len(vector_root_file)):
    if TString(vector_root_file[ifile].GetName()).Contains("%s"%(mass[imass])):
       masspoint = mass[imass];
       ifilePos = ifile ;
       break; 
 
   vector_root_file[ifile].cd();
   
   if masspoint == -1 or ifilePos == -1 : continue;
  
   if options.isMC == 0:
    histo_pull_data_wjet    = vector_root_file[ifilePos].Get("rrv_number_data%s_fit_%s%s_data_pull"%(options.mlvjregion,options.channel,spectrum));
    gaussian_pull_data_wjet = vector_root_file[ifilePos].Get("Gaussian_pull_rrv_number_data%s_fit_%s%s_data_pull"%(options.mlvjregion,options.channel,spectrum)); 
    if not options.onlybackgroundfit:
     histo_pull_signal       = vector_root_file[ifilePos].Get("rrv_number_signal_region_fit_ggH_vbfH_data_pull");
     gaussian_pull_signal    = vector_root_file[ifilePos].Get("Gaussian_pull_rrv_number_signal_region_fit_ggH_vbfH_data_pull"); 
   else:  
    histo_pull_data_wjet    = vector_root_file[ifilePos].Get("rrv_number_WJets0%s_fit_%s%s_data_pull"%(options.mlvjregion,options.channel,spectrum));
    gaussian_pull_signal    = vector_root_file[ifilePos].Get("Gaussian_pull_rrv_number_data%s_fit_%s%s_data_pull"%(options.mlvjregion,options.channel,spectrum)); 
    if not options.onlybackgroundfit:
     histo_pull_signal       = vector_root_file[ifilePos].Get("rrv_number_signal_region_fit_ggH_vbfH_data_pull");
     gaussian_pull_signal    = vector_root_file[ifilePos].Get("Gaussian_pull_rrv_number_signal_region_fit_ggH_vbfH_data_pull"); 

   canvas_pull_bkg_data.append(TCanvas("canvas_"+histo_pull_data_wjet.GetName()+"_mH%d"%(mass[imass]),""));
   canvas_pull_bkg_data[len(canvas_pull_bkg_data)-1].cd();
   histo_pull_data_wjet.Draw();
   canvas_pull_bkg_data[len(canvas_pull_bkg_data)-1].SaveAs(options.outputDir+"/"+canvas_pull_bkg_data[len(canvas_pull_bkg_data)-1].GetName()+".png","png");
   canvas_pull_bkg_data[len(canvas_pull_bkg_data)-1].SaveAs(options.outputDir+"/"+canvas_pull_bkg_data[len(canvas_pull_bkg_data)-1].GetName()+".pdf","pdf");

   if not options.onlybackgroundfit:
    canvas_pull_sig_data.append(TCanvas("canvas_"+histo_pull_sig_wjet.GetName()+"_mH%d"%(mass[imass]),""));
    canvas_pull_sig_data[len(canvas_pull_sig_data)-1].cd();
    histo_pull_signal.Draw();
    canvas_pull_sig_data[len(canvas_pull_sig_data)-1].SaveAs(options.outputDir+"/"+canvas_pull_sig_data[len(canvas_pull_sig_data)-1].GetName()+".png","png");
    canvas_pull_sig_data[len(canvas_pull_sig_data)-1].SaveAs(options.outputDir+"/"+canvas_pull_sig_data[len(canvas_pull_sig_data)-1].GetName()+".pdf","pdf");

   histogram_pull_vs_mass_nback.SetBinContent(imass+1,histo_pull_data_wjet.GetMean());
   histogram_pull_vs_mass_nback.SetBinError(imass+1,histo_pull_data_wjet.GetRMS());
   histogram_pull_vs_mass_nback.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));
   
   gaussian_pull_vs_mass_nback.SetBinContent(imass+1,gaussian_pull_data_wjet.GetParameter(1));
   gaussian_pull_vs_mass_nback.SetBinError(imass+1,gaussian_pull_data_wjet.GetParameter(2));
   gaussian_pull_vs_mass_nback.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));
   
   if not options.onlybackgroundfit:

    histogram_pull_vs_mass_nsig.SetBinContent(imass+1,histo_pull_signal.GetMean());
    histogram_pull_vs_mass_nsig.SetBinError(imass+1,histo_pull_signal.GetRMS());
    histogram_pull_vs_mass_nsig.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));

    gaussian_pull_vs_mass_nsig.SetBinContent(imass+1,gaussian_pull_signal.GetParameter(1));
    gaussian_pull_vs_mass_nsig.SetBinError(imass+1,gaussian_pull_signal.GetParameter(2));
    gaussian_pull_vs_mass_nsig.GetXaxis().SetBinLabel(imass+1,"%d"%(masspoint));


   graph_1.SetPoint(imass+1,imass,0);
   graph_1.SetPointError(imass+1,0,0.2);

   graph_2.SetPoint(imass+1,imass,0);
   graph_2.SetPointError(imass+1,0,1);

  graph_1.SetPoint(len(mass)+1,len(mass),0);
  graph_1.SetPointError(len(mass)+1,0,0.2);

  graph_2.SetPoint(len(mass)+1,len(mass),0);
  graph_2.SetPointError(len(mass)+1,0,1);

 ## make the output directory

 setPlotStyle();
           
 canvas_1 = ROOT.TCanvas("histo_bkg", "histo_bkg")
 histogram_pull_vs_mass_nback.GetXaxis().SetLabelSize(0.065);
 histogram_pull_vs_mass_nback.GetXaxis().SetTitleSize(0.055);
 histogram_pull_vs_mass_nback.GetXaxis().SetTitle("m_{H}");
 histogram_pull_vs_mass_nback.GetXaxis().SetTitleOffset(0.85);
 histogram_pull_vs_mass_nback.GetYaxis().SetTitle("Pull");
 histogram_pull_vs_mass_nback.GetYaxis().SetTitleOffset(1.05);
 histogram_pull_vs_mass_nback.GetYaxis().SetTitleSize(0.045);
 histogram_pull_vs_mass_nback.GetYaxis().SetLabelSize(0.045);
 histogram_pull_vs_mass_nback.SetMarkerStyle(20);
 histogram_pull_vs_mass_nback.SetMarkerSize(2.0);
 histogram_pull_vs_mass_nback.SetLineColor(1);
 histogram_pull_vs_mass_nback.SetLineWidth(2);
 graph_1.SetFillColor(kYellow+1);
 graph_1.SetFillStyle(3001);
 graph_2.SetFillColor(kGreen+1);
 graph_2.SetFillStyle(3001);

 oneLine = ROOT.TF1("oneLine","0",0,100);
 oneLine.SetLineColor(ROOT.kRed);
 oneLine.SetLineWidth(2);
                     
 histogram_pull_vs_mass_nback.SetMaximum(1.5);
 histogram_pull_vs_mass_nback.SetMinimum(-1.5);
 histogram_pull_vs_mass_nback.Draw("pe");
 graph_2.Draw("e3same");
 graph_1.Draw("e3same");
 oneLine.Draw("Lsame");
 histogram_pull_vs_mass_nback.Draw("pesame");

 banner = TLatex(0.21,1.52,("CMS Preliminary,                       19.3 fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow l #nu"));
 banner.SetTextSize(0.046);
 title  = TLatex(0.21,1.13,("N_{bkg} pull from histogram vs mass"));
 title.SetTextSize(0.042);
 title.SetTextFont(42); 

 title2  = TLatex(3.56,1.30,("Generated: %s")%(options.fgen));
 title3  = TLatex(3.56,1.12,("Fitted: %s")%(options.fres));
 title2.SetTextSize(0.042);
 title2.SetTextFont(42); 
 title2.SetTextColor(1); 
 title3.SetTextSize(0.042);
 title3.SetTextFont(42); 
 title3.SetTextColor(1);
 
 banner.Draw();
 title.Draw();
 title2.Draw();
 title3.Draw();

 canvas_1.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".png","png");
 canvas_1.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".pdf","pdf");
 canvas_1.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".root","root");

 if not options.onlybackgroundfit:
  canvas_2 = ROOT.TCanvas("histo_sig", "histo_sig")
  histogram_pull_vs_mass_nsig.GetXaxis().SetLabelSize(0.065);
  histogram_pull_vs_mass_nsig.GetXaxis().SetTitleSize(0.055);
  histogram_pull_vs_mass_nsig.GetXaxis().SetTitle("m_{H}");
  histogram_pull_vs_mass_nsig.GetXaxis().SetTitleOffset(0.85);
  histogram_pull_vs_mass_nsig.GetYaxis().SetTitle("Pull");
  histogram_pull_vs_mass_nsig.GetYaxis().SetTitleOffset(1.05);
  histogram_pull_vs_mass_nsig.GetYaxis().SetLabelSize(0.045);
  histogram_pull_vs_mass_nsig.GetYaxis().SetTitleSize(0.045);
  histogram_pull_vs_mass_nsig.SetMarkerStyle(20);
  histogram_pull_vs_mass_nsig.SetMarkerSize(2.0);
  histogram_pull_vs_mass_nsig.SetLineColor(1);
  histogram_pull_vs_mass_nsig.SetLineWidth(2);

  oneLine = ROOT.TF1("oneLine","0",0,100);
  oneLine.SetLineColor(ROOT.kRed);
  oneLine.SetLineWidth(2);
                     
  histogram_pull_vs_mass_nsig.SetMaximum(1.5);
  histogram_pull_vs_mass_nsig.SetMinimum(-1.5);
  histogram_pull_vs_mass_nsig.Draw("pe");
  graph_2.Draw("e3same");
  graph_1.Draw("e3same");
  oneLine.Draw("Lsame");
  histogram_pull_vs_mass_nsig.Draw("pesame");

  banner.Draw();
  title  = TLatex(0.21,1.13,("N_{sig} pull from histogram vs mass"));
  title.SetTextSize(0.042);
  title.SetTextFont(42); 
  title.Draw();
  title2.Draw();
  title3.Draw();


  canvas_2.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".png","png");
  canvas_2.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".pdf","pdf");
  canvas_2.SaveAs(options.outputDir+"/"+histogram_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".root","root");


 canvas_3 = ROOT.TCanvas("gaus_bkg", "gaus_bkg")
 gaussian_pull_vs_mass_nback.GetXaxis().SetLabelSize(0.065);
 gaussian_pull_vs_mass_nback.GetXaxis().SetTitleSize(0.055);
 gaussian_pull_vs_mass_nback.GetXaxis().SetTitle("m_{H}");
 gaussian_pull_vs_mass_nback.GetXaxis().SetTitleOffset(0.85);
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitle("Pull");
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitleOffset(1.05);
 gaussian_pull_vs_mass_nback.GetYaxis().SetLabelSize(0.045);
 gaussian_pull_vs_mass_nback.GetYaxis().SetTitleSize(0.045);
 gaussian_pull_vs_mass_nback.SetTitle("Pull of Background Number of Events Gaussian")
 gaussian_pull_vs_mass_nback.SetMarkerStyle(20);
 gaussian_pull_vs_mass_nback.SetMarkerSize(2.0);
 gaussian_pull_vs_mass_nback.SetLineColor(1);
 gaussian_pull_vs_mass_nback.SetLineWidth(2);


 oneLine = ROOT.TF1("oneLine","0",0,100);
 oneLine.SetLineColor(ROOT.kRed);
 oneLine.SetLineWidth(2);
                     
 gaussian_pull_vs_mass_nback.SetMaximum(1.5);
 gaussian_pull_vs_mass_nback.SetMinimum(-1.5);
 gaussian_pull_vs_mass_nback.Draw("pe");
 graph_2.Draw("e3same");
 graph_1.Draw("e3same");
 oneLine.Draw("Lsame");

 gaussian_pull_vs_mass_nback.Draw("pesame");

 banner.Draw();
 title  = TLatex(0.21,1.13,("N_{back} pull from gaussian fit vs mass"));
 title.SetTextSize(0.042);
 title.SetTextFont(42); 
 title.Draw();
 title2.Draw();
 title3.Draw();

 canvas_3.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".png","png");
 canvas_3.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".pdf","pdf");
 canvas_3.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nback.GetName()+"_"+options.fgen+"_"+options.fres+".root","root");

 if not options.onlybackgroundfit:

  canvas_4 = ROOT.TCanvas("gaus_sig", "gaus_sig")
  gaussian_pull_vs_mass_nsig.GetXaxis().SetLabelSize(0.065);
  gaussian_pull_vs_mass_nsig.GetXaxis().SetTitleSize(0.055);
  gaussian_pull_vs_mass_nsig.GetXaxis().SetTitle("m_{H}");
  gaussian_pull_vs_mass_nsig.GetXaxis().SetTitleOffset(0.85);
  gaussian_pull_vs_mass_nsig.GetYaxis().SetTitle("Pull");
  gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleOffset(1.05);
  gaussian_pull_vs_mass_nsig.GetYaxis().SetLabelSize(0.045);
  gaussian_pull_vs_mass_nsig.GetYaxis().SetTitleSize(0.045);
  gaussian_pull_vs_mass_nsig.SetTitle("Pull of Signal Number of Events Gaussian")
  gaussian_pull_vs_mass_nsig.SetMarkerStyle(20);
  gaussian_pull_vs_mass_nsig.SetMarkerSize(2.0);
  gaussian_pull_vs_mass_nsig.SetLineColor(1);
  gaussian_pull_vs_mass_nsig.SetLineWidth(2);

  gaussian_pull_vs_mass_nsig.SetMaximum(1.5);
  gaussian_pull_vs_mass_nsig.SetMinimum(-1.5);
  gaussian_pull_vs_mass_nsig.Draw("pe");
  graph_2.Draw("e3same");
  graph_1.Draw("e3same");
  oneLine.Draw("Lsame");
  gaussian_pull_vs_mass_nsig.Draw("pesame");

  banner.Draw();

  title  = TLatex(0.21,1.13,("N_{sig} pull from gaussian fit vs mass"));
  title.SetTextSize(0.042);
  title.SetTextFont(42); 
  title.Draw();
  title2.Draw();
  title3.Draw();

  canvas_4.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".png","png");
  canvas_4.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".pdf","pdf");
  canvas_4.SaveAs(options.outputDir+"/"+gaussian_pull_vs_mass_nsig.GetName()+"_"+options.fgen+"_"+options.fres+".root","root");

 os.system("rm list_temp.txt");
 
