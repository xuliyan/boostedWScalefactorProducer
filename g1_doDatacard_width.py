#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import string
import subprocess
from subprocess import Popen
from optparse import OptionParser

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, TList

############################################
#              Job steering                #
############################################

parser = OptionParser()

parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="em")

parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--datacardPath', action="store",type="string",dest="datacardPath",default="cards_em_EXO_allCat_v2_ExpTail_g1_rereco_c0p5_modelIndependent/")

parser.add_option('--category', action="store",type="string",dest="category",default="HP")

(options, args) = parser.parse_args()

#######################################################################

ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf,RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

#######################################################################

##### Point of Mass to be Analyzed 
mass          = [1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000,2100,2200,2300,2400,2500,800,900]
GammaOverMass = [0.001,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40]


### how to run it: python g1_doDatacard_width.py wwlvj_BulkG_WW_inclusive_c0p2 -b  --datacardPath cards_em_EXO_allCat_v2_ExpTail_g1_rereco_c0p5/ --channel em 

    ## Set basic TDR style for canvas, pad ..etc ..
def setTDRStyle():

 gStyle.SetCanvasBorderMode(0);
 gStyle.SetCanvasColor(kWhite);
 gStyle.SetCanvasDefH(600); #Height of canvas
 gStyle.SetCanvasDefW(600); #Width of canvas
 gStyle.SetCanvasDefX(0); #POsition on screen
 gStyle.SetCanvasDefY(0);
      
 #For the Pad:
 gStyle.SetPadBorderMode(0);
 gStyle.SetPadColor(kWhite);
 gStyle.SetPadGridX(False);
 gStyle.SetPadGridY(False);
 gStyle.SetGridColor(0);
 gStyle.SetGridStyle(3);
 gStyle.SetGridWidth(1);
      
 #For the frame:
 gStyle.SetFrameBorderMode(0);
 gStyle.SetFrameBorderSize(1);
 gStyle.SetFrameFillColor(0);
 gStyle.SetFrameFillStyle(0);
 gStyle.SetFrameLineColor(1);
 gStyle.SetFrameLineStyle(1);
 gStyle.SetFrameLineWidth(1);
      
 #For the histo:
 gStyle.SetHistLineColor(1);
 gStyle.SetHistLineStyle(0);
 gStyle.SetHistLineWidth(1);
 gStyle.SetEndErrorSize(2);
 gStyle.SetErrorX(0.);
 gStyle.SetMarkerStyle(20);
      
 #For the fit/function:
 gStyle.SetOptFit(1);
 gStyle.SetFitFormat("5.4g");
 gStyle.SetFuncColor(2);
 gStyle.SetFuncStyle(1);
 gStyle.SetFuncWidth(1);
      
 #For the date:
 gStyle.SetOptDate(0);
      
 #For the statistics box:
 gStyle.SetOptFile(0);
 gStyle.SetOptStat(0); #To display the mean and RMS:
 gStyle.SetStatColor(kWhite);
 gStyle.SetStatFont(42);
 gStyle.SetStatFontSize(0.025);
 gStyle.SetStatTextColor(1);
 gStyle.SetStatFormat("6.4g");
 gStyle.SetStatBorderSize(1);
 gStyle.SetStatH(0.1);
 gStyle.SetStatW(0.15);
      
 #Margins:
 gStyle.SetPadTopMargin(0.05);
 gStyle.SetPadBottomMargin(0.13);
 gStyle.SetPadLeftMargin(0.18);
 gStyle.SetPadRightMargin(0.06);
      
 #For the Global title:
 gStyle.SetOptTitle(0);
 gStyle.SetTitleFont(42);
 gStyle.SetTitleColor(1);
 gStyle.SetTitleTextColor(1);
 gStyle.SetTitleFillColor(10);
 gStyle.SetTitleFontSize(0.05);
      
 #For the axis titles:
 gStyle.SetTitleColor(1, "XYZ");
 gStyle.SetTitleFont(42, "XYZ");
 gStyle.SetTitleSize(0.03, "XYZ");
 gStyle.SetTitleXOffset(0.9);
 gStyle.SetTitleYOffset(1.5);
      
 #For the axis labels:
 gStyle.SetLabelColor(1, "XYZ");
 gStyle.SetLabelFont(42, "XYZ");
 gStyle.SetLabelOffset(0.007, "XYZ");
 gStyle.SetLabelSize(0.03, "XYZ");
      
 #For the axis:
 gStyle.SetAxisColor(1, "XYZ");
 gStyle.SetStripDecimals(kTRUE);
 gStyle.SetTickLength(0.03, "XYZ");
 gStyle.SetNdivisions(510, "XYZ");
 gStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
 gStyle.SetPadTickY(1);
      
 #Change for log plots:
 gStyle.SetOptLogx(0);
 gStyle.SetOptLogy(0);
 gStyle.SetOptLogz(0);
      
 #Postscript options:
 gStyle.SetPaperSize(20.,20.);
 gStyle.cd();

 
##### Get Lumi for banner title
def GetLumi():

   if options.channel=="el": return 19.7;
   elif options.channel=="mu": return 19.7;
   elif options.channel=="em": return 19.7;


#### in order to make the banner on the plots
def banner4Plot(iswithpull=0):
     
 print "############### draw the banner ########################"

 if iswithpull:
    if options.channel=="el":
       banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(GetLumi())));
    elif options.channel=="mu":
       banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(GetLumi())));
    elif options.channel=="em":
       banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(GetLumi())));
    banner.SetNDC(); banner.SetTextSize(0.04);
 else:
    if options.channel=="el":
       banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(GetLumi())));
    if options.channel=="mu":
       banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(GetLumi())));
    if options.channel=="em":
       banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu,e #nu "%(GetLumi())));
    banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
 return banner;


#### just drawing canvas with no pull
def draw_canvas(in_obj,in_directory,in_file_name,is_range=0,logy=0):

  cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

  if is_range:
    h2=TH2D("h2","",100,400,1400,4,0.00001,4);
    h2.Draw();
    in_obj.Draw("same")
  else :
    in_obj.Draw()

  in_obj.GetXaxis().SetTitleSize(0.045);
  in_obj.GetXaxis().SetTitleOffset(1.15);
  in_obj.GetXaxis().SetLabelSize(0.04);

  in_obj.GetYaxis().SetTitleSize(0.045);
  in_obj.GetYaxis().SetTitleOffset(1.40);
  in_obj.GetYaxis().SetLabelSize(0.04);

  banner = banner4Plot();
  banner.Draw();
        
  Directory=TString(in_directory);
  if not Directory.EndsWith("/"):Directory=Directory.Append("/");
  if not os.path.isdir(Directory.Data()):
    os.system("mkdir -p "+Directory.Data());

  rlt_file = TString();
  rlt_file.Form("%s_%s"%(Directory,in_file_name));
  if rlt_file.EndsWith(".root"):
     rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
  else:
     rlt_file.ReplaceAll(".root","");
     rlt_file = rlt_file.Append(".png");

  cMassFit.SaveAs(rlt_file.Data());

  rlt_file.ReplaceAll(".png",".pdf");
  cMassFit.SaveAs(rlt_file.Data());

  rlt_file.ReplaceAll(".pdf",".root");
  cMassFit.SaveAs(rlt_file.Data());

  if logy:
    cMassFit.SetLogy() ;
    cMassFit.Update();
    rlt_file.ReplaceAll(".root","_log.root");
    cMassFit.SaveAs(rlt_file.Data());
    rlt_file.ReplaceAll(".root",".pdf");
    cMassFit.SaveAs(rlt_file.Data());
    rlt_file.ReplaceAll(".pdf",".png");
    cMassFit.SaveAs(rlt_file.Data());


#### draw canvas with plots with pull
def draw_canvas_with_pull(mplot, mplot_pull,in_directory, in_file_name, in_model_name="",logy=0):# mplot + pull

 mplot.GetXaxis().SetTitleOffset(1.1);
 mplot.GetYaxis().SetTitleOffset(1.3);
 mplot.GetXaxis().SetTitleSize(0.05);
 mplot.GetYaxis().SetTitleSize(0.05);
 mplot.GetXaxis().SetLabelSize(0.045);
 mplot.GetYaxis().SetLabelSize(0.045);
 mplot_pull.GetXaxis().SetLabelSize(0.15);
 mplot_pull.GetYaxis().SetLabelSize(0.15);
 mplot_pull.GetYaxis().SetTitleSize(0.15);
 mplot_pull.GetYaxis().SetNdivisions(205);

                                                                          
 cMassFit = TCanvas("cMassFit","cMassFit", 600,600);
 pad1=TPad("pad1","pad1",0.,0. ,0.99,0.24);
 pad2=TPad("pad2","pad2",0.,0.24,0.99,1. );
 pad1.Draw();
 pad2.Draw();
                                                                                                                                                                              
 pad2.cd();
 mplot.Draw();
 banner = banner4Plot(1);
 banner.Draw();

 pad1.cd();
 mplot_pull.Draw();


 ## create the directory where store the plots
 Directory = TString(in_directory);
 if not Directory.EndsWith("/"):Directory = Directory.Append("/");
 if not os.path.isdir(Directory.Data()):
   os.system("mkdir -p "+Directory.Data());

 rlt_file = TString();
 rlt_file.Form("%s_%s"%(Directory,in_file_name));
 if rlt_file.EndsWith(".root"):
    TString(in_model_name).ReplaceAll(".root","");
    rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
 else:
   TString(in_model_name).ReplaceAll(".root","");
   rlt_file.ReplaceAll(".root","");
   rlt_file=rlt_file.Append("_"+in_model_name+"_with_pull.png");

 cMassFit.SaveAs(rlt_file.Data());

 rlt_file.ReplaceAll(".png",".pdf");
 cMassFit.SaveAs(rlt_file.Data());
        
 rlt_file.ReplaceAll(".pdf",".root");
 cMassFit.SaveAs(rlt_file.Data());

 string_file_name = TString(in_file_name);

 if string_file_name.EndsWith(".root"):
   string_file_name.ReplaceAll(".root","_"+in_model_name);
 else:
   string_file_name.ReplaceAll(".root","");
   string_file_name.Append("_"+in_model_name);

 if logy:
   pad2.SetLogy() ;
   pad2.Update();
   cMassFit.Update();
   rlt_file.ReplaceAll(".root","_log.root");
   cMassFit.SaveAs(rlt_file.Data());
   rlt_file.ReplaceAll(".root",".pdf");
   cMassFit.SaveAs(rlt_file.Data());
   rlt_file.ReplaceAll(".pdf",".png");
   cMassFit.SaveAs(rlt_file.Data());

 draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy);


### in order to get the pull
def get_pull(rrv_x, mplot_orig):

  hpull = mplot_orig.pullHist();
  x = ROOT.Double(0.); y = ROOT.Double(0) ;
  for ipoint in range(0,hpull.GetN()):
    hpull.GetPoint(ipoint,x,y);
    if(y == 0): hpull.SetPoint(ipoint,x,10)
       
  mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins())));
  medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
  mplot_pull.addObject(medianLine);
  mplot_pull.addPlotable(hpull,"P");
  mplot_pull.SetTitle("PULL");
  mplot_pull.GetXaxis().SetTitle("");
  mplot_pull.GetYaxis().SetRangeUser(-5,5);
  mplot_pull.GetYaxis().SetTitleSize(0.10);
  mplot_pull.GetYaxis().SetLabelSize(0.10);
  mplot_pull.GetXaxis().SetTitleSize(0.10);
  mplot_pull.GetXaxis().SetLabelSize(0.10);
  mplot_pull.GetYaxis().SetTitleOffset(0.40);
  mplot_pull.GetYaxis().SetTitle("#frac{data-fit}{#sigma_{data}}");
  mplot_pull.GetYaxis().CenterTitle();

  return mplot_pull;



#### Main Code
if __name__ == '__main__':

  print "################# Programmed to generate datacard for large widths ";

  ## make the list of the exsisting datacard 
  command = "ls %s/%s | grep %s | grep -v other | grep -v combo | grep -v counting | grep -v _W_ | grep %s | grep txt > %s/%s/temp_datacard_list.txt"%(options.inPath,options.datacardPath,sys.argv[1],options.channel,options.inPath,options.datacardPath);
  print command ;
  os.system(command);
  

  ## make the list of the exsisting workspaces
  command = "ls %s/%s | grep %s | grep -v other | grep -v combo | grep -v counting | grep -v _W_ | grep %s | grep root > %s/%s/temp_workspace_list.txt"%(options.inPath,options.datacardPath,sys.argv[1],options.channel,options.inPath,options.datacardPath);
  print command ;
  os.system(command);

  setTDRStyle();

  ## list of datacard and workspaces
  list_of_datacard = [] ;
  list_of_workspace = [] ;

  datacardFile = open("%s/%s/temp_datacard_list.txt"%(options.inPath,options.datacardPath));

  for fileLine in datacardFile :

    if not fileLine.find("#")  : continue ;
    if fileLine.split()[0] != "" and fileLine.find("txt") !=-1 and fileLine.find("unbin") !=-1 :  
     list_of_datacard.append(fileLine.split()[0]);


  workspaceFile = open("%s/%s/temp_workspace_list.txt"%(options.inPath,options.datacardPath));

  for fileLine in workspaceFile :

    if not fileLine.find("#")  : continue ;
    if fileLine.split()[0] != "" and fileLine.find("root") !=-1 and fileLine.find("workspace") !=-1 :  
     list_of_workspace.append(fileLine.split()[0]);

  ### print the new datacard with the new workspace names

  idatacard = 0 ;
  datasetname = "";
  
  for datacard in list_of_datacard :
    for gammaVal in GammaOverMass:

        newdatacard = TString(datacard).ReplaceAll(".txt","_W%.3f.txt"%gammaVal);
        newdatacard.ReplaceAll("0.","_");
        newdatacard.ReplaceAll("_txt",".txt");

        datacard_old = open("%s/%s/%s"%(options.inPath,options.datacardPath,datacard),"r");        
        datacardfile = open("%s/%s/%s"%(options.inPath,options.datacardPath,newdatacard),"a");
        
        for line_old in datacard_old:
          if line_old.find("workspace")!=-1:
            newworkspace = TString(list_of_workspace[idatacard]).ReplaceAll(".root","_W%.3f.root"%gammaVal);
            newworkspace.ReplaceAll("0.","_");
            newworkspace.ReplaceAll("_root",".root");
            if line_old.split()[1].find("Bulk")!=-1:
             datacardfile.write("%s %s %s %s %s"%(line_old.split()[0],line_old.split()[1],line_old.split()[2],newworkspace,line_old.split()[4]));
            elif line_old.split()[1].find("data")!=-1:
             datasetname = line_old.split()[1]; 
             datacardfile.write("\n%s %s %s %s %s\n"%(line_old.split()[0],line_old.split()[1],line_old.split()[2],newworkspace,line_old.split()[4]));
            else:
             datacardfile.write("\n%s %s %s %s %s"%(line_old.split()[0],line_old.split()[1],line_old.split()[2],newworkspace,line_old.split()[4]));
          elif line_old.find("rate")!=-1:
            datacardfile.write("rate 1 %s %s %s %s\n"%(line_old.split()[2],line_old.split()[3],line_old.split()[4],line_old.split()[5]));
            continue ;
          elif line_old.find("lumi_8TeV")!=-1:
            datacardfile.write("CMS_eff_width lnN 1.15 - - - -\n");
            datacardfile.write(line_old);            
          else:  
            datacardfile.write(line_old);
    idatacard = idatacard+1;

  ### make the new workspaces --> loop on the old, create empty new ones  
  iMass = 0 ;
  for workspace in list_of_workspace :
    for gammaVal in GammaOverMass:

        newfilename = TString(workspace).ReplaceAll(".root","_W%.3f.root"%gammaVal);
        newfilename.ReplaceAll("0.","_");
        newfilename.ReplaceAll("_root",".root");

        old_file  = TFile("%s/%s/%s"%(options.inPath,options.datacardPath,workspace),"READ");
        new_file  = TFile("%s/%s/%s"%(options.inPath,options.datacardPath,newfilename),"RECREATE");

        old_workspace = old_file.Get("workspace4limit_");
        new_workspace = RooWorkspace("workspace4limit_","workspace4limit_");

        new_file.cd();

        ### copy all the pdfs + parameter
        old_workspace.Print();
        pdf_workspace = old_workspace.allPdfs();
        par = pdf_workspace.createIterator();
        par.Reset();
        param = par.Next()

        while (param):                   
          if not TString(par.GetName()).Contains("BulkWW_xww") and not TString(par.GetName()).Contains("BulkG_WW"):
            getattr(new_workspace,"import")(param);
          else:
            RooDoubleCrystalBall = param;
          param=par.Next()

        ## cycle on DoubleCB parameters -> copy them in a new set
        sig_parameters = RooDoubleCrystalBall.getParameters(old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category));
        sig_par = sig_parameters.createIterator();
        sig_par.Reset();
        sig_param = sig_par.Next();

        while (sig_param):
          if TString(sig_param.GetName()).Contains("rrv_mean_CB_BulkG_WW"): ## copy the mean value for the DoubleCB and BW
           rrv_mean_BW = RooRealVar("rrv_mean_BW_BulkG_WW_"+options.channel+"_"+options.category,"rrv_mean_BW_BulkG_WW_"+options.channel+"_"+options.category,sig_param.getVal());
           rrv_mean_CB = sig_param;           
           rrv_mean_CB.setRange(-100,100);
           rrv_mean_CB.setVal(0.);
          elif TString(sig_param.GetName()).Contains("rrv_sigma_CB_BulkG_WW"):
           rrv_sigma_CB = sig_param;
          elif TString(sig_param.GetName()).Contains("rrv_alpha1_CB_BulkG_WW"):
           rrv_alpha1_CB = sig_param;
          elif TString(sig_param.GetName()).Contains("rrv_alpha2_CB_BulkG_WW"):
           rrv_alpha2_CB = sig_param;
          elif TString(sig_param.GetName()).Contains("rrv_n1_CB_BulkG_WW"):
           rrv_n1_CB = sig_param;
          elif TString(sig_param.GetName()).Contains("rrv_n2_CB_BulkG_WW"):
           rrv_n2_CB = sig_param;
          sig_param = sig_par.Next();

        ### copy all the datasets
        getattr(new_workspace,"import")(old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category));

        ## Breit Wigner Core  --> mean and the resolution fixed
        rrv_width_BW = RooRealVar("rrv_width_BW_BulkG_WW_"+options.channel+"_"+options.category,"rrv_width_BW_BulkG_WW_"+options.channel+"_"+options.category,mass[iMass]*gammaVal);

        rrv_mean_BW.setConstant(kTRUE);
        rrv_width_BW.setConstant(kTRUE);

        ## Breit Wigner Core  --> mean sys are to be applied on the mean of BW
        rrv_mean_scale_p1 = RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);
        rrv_mean_scale_p1.setConstant(kTRUE);
        rrv_mean_scale_p1.setError(1);
         
        rrv_mean_scale_p2 = RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
        rrv_mean_scale_p2.setConstant(kTRUE);
        rrv_mean_scale_p2.setError(1);
                                    
        rrv_mean_scale_X1 = RooRealVar("rrv_mean_shift_scale_lep_BW_"+options.channel+"_"+options.category,"rrv_mean_shift_scale_lep_BW_"+options.channel+"_"+options.category,float(0.001));
        rrv_mean_scale_X1.setConstant(kTRUE);
        rrv_mean_scale_X2 = RooRealVar("rrv_mean_shift_scale_jes_BW_"+options.channel+"_"+options.category,"rrv_mean_shift_scale_jes_BW_"+options.channel+"_"+options.category,float(0.033));
        rrv_mean_scale_X2.setConstant(kTRUE);

        rrv_total_mean_CB = RooFormulaVar("rrv_total_mean_BW_"+options.channel,"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(rrv_mean_BW,rrv_mean_scale_p1,rrv_mean_scale_X1,rrv_mean_scale_p2,rrv_mean_scale_X2));

        bw = RooBreitWigner("bw_BulkWW_xww_"+options.channel,"bw_BulkG_WW_"+options.channel,old_workspace.var("rrv_mass_lvj"),rrv_total_mean_CB,rrv_width_BW);


        ### make resolution DoubleCB
        rrv_sigma_scale_p1 = RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
        rrv_sigma_scale_p1.setConstant(kTRUE);
        rrv_sigma_scale_p1.setError(1);
        
        rrv_sigma_scale_p2 = RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
        rrv_sigma_scale_p3 = RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);        
        rrv_sigma_scale_p2.setConstant(kTRUE);
        rrv_sigma_scale_p3.setConstant(kTRUE);
        rrv_sigma_scale_p2.setError(1);
        rrv_sigma_scale_p3.setError(1);

        rrv_mean_sigma_X1 = RooRealVar("rrv_sigma_shift_lep_scale_CB"+options.channel+"_"+options.category,"rrv_sigma_shift_scale_CB"+options.channel+"_"+options.category,float(0.005));
        rrv_mean_sigma_X2 = RooRealVar("rrv_sigma_shift_jes_CB"+options.channel+"_"+options.category,"rrv_sigma_shift_scale_CB"+options.channel+"_"+options.category,float(0.033));
        rrv_mean_sigma_X3 = RooRealVar("rrv_sigma_shift_res_CB"+options.channel+"_"+options.category,"rrv_sigma_shift_res_CB"+options.channel+"_"+options.category,float(0.030));
        rrv_mean_sigma_X1.setConstant(kTRUE);
        rrv_mean_sigma_X2.setConstant(kTRUE);
        rrv_mean_sigma_X3.setConstant(kTRUE);

        rrv_total_sigma_CB = RooFormulaVar("rrv_total_sigma_CB"+options.channel+"_"+options.category,"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(rrv_sigma_CB,rrv_sigma_scale_p1,rrv_mean_sigma_X1,rrv_sigma_scale_p2,rrv_mean_sigma_X2,rrv_sigma_scale_p3,rrv_mean_sigma_X3));

        new_signal = ROOT.RooDoubleCrystalBall("DoubleCB_BulkG_WW_"+options.channel+"_"+options.category+"_mlvj","DoubleCB_BulkG_WW_"+options.channel+"_"+options.category+"_mlvj",old_workspace.var("rrv_mass_lvj"),rrv_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB); 


        ### FFT ConvPdf
        model_pdf = RooFFTConvPdf("BulkWW_xww_%s_%s"%(options.channel,options.category),"BulkWW_xww_%s_%s"%(options.channel,options.category),old_workspace.var("rrv_mass_lvj"),bw,new_signal);
        model_pdf.setBufferFraction(1.0);
        
        model_pdf.SetName("BulkWW_xww_%s_%s"%(options.channel,options.category));
        getattr(new_workspace,"import")(model_pdf);
                          
        ### iterate on the workspace element parameters
        print "----------- Parameter Workspace -------------";
        parameters_workspace = new_workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
         param.Print();
         param=par.Next()
        print "---------------------------------------------";

                                                                                  
        ## create a frame with data
        mplot = old_workspace.var("rrv_mass_lvj").frame(RooFit.Title(""), RooFit.Bins(int(old_workspace.var("rrv_mass_lvj").getBins()/10)),RooFit.Range(700,3000));

        old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category).plotOn(mplot,RooFit.Name("data_invisible"),RooFit.MarkerSize(1.5),RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0));
   
        rrv_number_WJets = old_workspace.var("rate_WJets_xww_for_unbin");
        rrv_number_VV = old_workspace.var("rate_VV_xww_for_unbin");
        rrv_number_TTbar = old_workspace.var("rate_TTbar_xww_for_unbin");
        rrv_number_STop = old_workspace.var("rate_STop_xww_for_unbin");
        rrv_number_signal = old_workspace.var("rate_BulkWW_xww_for_unbin");
        
        rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC_xww","rrv_number_Total_background_MC_xww",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal());

        rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()* rrv_number_WJets.getError()+
                rrv_number_VV.getError()* rrv_number_VV.getError()+
                rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
                rrv_number_STop.getError() *rrv_number_STop.getError()
                ));

        model_Total_background_MC = RooAddPdf("model_Total_background_MC_xww","model_Total_background_MC_xww",RooArgList(old_workspace.pdf("WJets_xww_%s_%s"%(options.channel,options.category)), old_workspace.pdf("VV_xww_%s_%s"%(options.channel,options.category)),old_workspace.pdf("TTbar_xww_%s_%s"%(options.channel,options.category)),old_workspace.pdf("STop_xww_%s_%s"%(options.channel,options.category))),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_TTbar,rrv_number_STop));

        rrv_number_signal.setVal(rrv_number_signal.getVal()*6.25);

        #### scale factor in order to scale MC to data in the final plot -> in order to avoid the normalization to data which is done by default in rooFit
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category).sumEntries();
        scale_number_signal = rrv_number_signal.getVal()/old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category).sumEntries();

        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("total_MC"),RooFit.Components("WJets_xww_%s_%s,VV_xww_%s_%s,TTbar_xww_%s_%s,STop_xww_%s_%s"%(options.channel,options.category,options.channel,options.category,options.channel,options.category,options.channel,options.category)),RooFit.DrawOption("L"), RooFit.LineColor(kRed), RooFit.VLines(),RooFit.LineWidth(2));

            
        model_signal_background_MC = RooAddPdf("model_signal_background_MC_xww","model_signal_background_MC_xww",RooArgList(model_pdf,model_Total_background_MC),RooArgList(rrv_number_signal,rrv_number_Total_background_MC));

        model_signal_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC+scale_number_signal),RooFit.Name("total_SpB_MC"),RooFit.Components("BulkWW_xww_%s_%s,model_Total_background_MC_xww"%(options.channel,options.category)),RooFit.DrawOption("L"), RooFit.LineColor(kBlue), RooFit.VLines(),RooFit.LineWidth(2),RooFit.LineStyle(7));

        model_pdf.plotOn(mplot,RooFit.Name("total_S_MC"),RooFit.Normalization(scale_number_signal),RooFit.DrawOption("L"), RooFit.LineColor(kGreen+2), RooFit.VLines(),RooFit.LineWidth(2),RooFit.LineStyle(kDashed));

        os.system("mkdir -p ./plots_signal_width");
        name = TString("check_%.3f_%.3f"%(mass[iMass],gammaVal));
        name.ReplaceAll("0.","0_");
        mplot.GetYaxis().SetRangeUser(1e-3,mplot.GetMaximum()*1.2);
#        draw_canvas(mplot,"plots_signal_width/",name,0,1);
        
        new_workspace.writeToFile(new_file.GetName());

        if ( (mass[iMass] == 1000 and gammaVal == 0.05) or (mass[iMass] == 1000 and gammaVal == 0.15) or (mass[iMass] == 1000 and gammaVal == 0.30) or
             (mass[iMass] == 1500 and gammaVal == 0.05) or (mass[iMass] == 1500 and gammaVal == 0.15) or (mass[iMass] == 1500 and gammaVal == 0.30) or
             (mass[iMass] == 2100 and gammaVal == 0.05) or (mass[iMass] == 2100 and gammaVal == 0.15) or (mass[iMass] == 2100 and gammaVal == 0.30) ):

           signal_file_name = ("treeEDBR_BulkG_WW_inclusive_M%3d_W%d_xww.root"%(mass[iMass],int(mass[iMass]*gammaVal)));
           print "#### making plot for ",signal_file_name;               
           file_Directory = "AnaSigTree_new/";
           fileIn_name = TString(options.inPath+"/"+file_Directory+signal_file_name);

           if gammaVal == 0.05   : sf = 5 ;
           elif gammaVal == 0.15 : sf = 3 ;
           elif gammaVal == 0.30 : sf = 1.5 ;

           if gammaVal == 0.05 or mass[iMass] == 1000 : binWidth = 20;
           elif gammaVal == 0.15 or mass[iMass] == 1500 : binWidth = 30;
           elif gammaVal == 0.30 or mass[iMass] == 1500 : binWidth = 50;
           elif gammaVal == 0.15 or mass[iMass] == 2000 : binWidth = 50;
           elif gammaVal == 0.30 or mass[iMass] == 2000 : binWidth = 75;
           else:
             binWidth = 20;
             
           rrv_mass_lvj = RooRealVar("rrv_mass_lvj_plot","M_{WW}",float(mass[iMass]),float(mass[iMass]-mass[iMass]*gammaVal*sf),float(mass[iMass]+mass[iMass]*gammaVal*sf),"GeV/c^{2}");
           rrv_mass_lvj.setBins(int((rrv_mass_lvj.getMax()-rrv_mass_lvj.getMin())/binWidth));
           rrv_weight   = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.);

           ## open the signal file 
           fileIn = TFile(fileIn_name.Data());
           treeIn = fileIn.Get("SelectedCandidatesPlain");


           rdataset4fit_mlvj = RooDataSet("rdataset4fit_BulkG_WW_inclusive_M%3d_W%d_%s_mlvj"%(int(mass[iMass]),int(mass[iMass]*gammaVal),options.channel),"rdataset4fit_BulkG_WW_inclusive_M%3d_W%d_%s_mlvj"%(int(mass[iMass]),int(mass[iMass]*gammaVal),options.channel),RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight));
                   
           luminosity = 19700;

           for iEvent in range(treeIn.GetEntries()):

            if iEvent % 1000 == 0: print "iEvent: ",iEvent;
            treeIn.GetEntry(iEvent);

            if iEvent==0:
                tmp_scale_to_lumi = getattr(treeIn,"LumiWeight")*luminosity;
    
            tmp_jet_mass=getattr(treeIn,"mJJNoKinFit");
            isGoodEvent = 0 ;
            
            ## event in the whole range
            if options.channel == "em" :
             if (getattr(treeIn,"categories")==1 or getattr(treeIn,"categories")==3) and getattr(treeIn,"mZZ")> rrv_mass_lvj.getMin() and getattr(treeIn,"mZZ")<rrv_mass_lvj.getMax() and tmp_jet_mass>=65 and tmp_jet_mass<=105 and options.category == "HP" :
              isGoodEvent = 1 ;
             if (getattr(treeIn,"categories")==0 or getattr(treeIn,"categories")==2) and getattr(treeIn,"mZZ")> rrv_mass_lvj.getMin() and getattr(treeIn,"mZZ")<rrv_mass_lvj.getMax() and tmp_jet_mass>=65 and tmp_jet_mass<=130 and options.category == "LP" :
              isGoodEvent = 1 ;

            if isGoodEvent == 1:
                ### weigh MC events
                tmp_event_weight = getattr(treeIn,"weight")*luminosity;
                rrv_mass_lvj.setVal(getattr(treeIn,"mZZ"));
                rdataset4fit_mlvj.add(RooArgSet(rrv_mass_lvj),tmp_event_weight);

           ### signal over convolution plot
           mplot_sig = rrv_mass_lvj.frame(RooFit.Title(""),RooFit.Bins(rrv_mass_lvj.getBins()));

           rdataset4fit_mlvj.plotOn(mplot_sig,RooFit.Name("data_invisible"),RooFit.MarkerSize(1.5),RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0));


           doubleCB_sig = ROOT.RooDoubleCrystalBall("DoubleCB_BulkG_WW_mlvj","DoubleCB_BulkG_WW_mlvj",rrv_mass_lvj,rrv_mean_CB,rrv_total_sigma_CB,rrv_alpha1_CB,rrv_n1_CB,rrv_alpha2_CB,rrv_n2_CB); 
           bw_sig = RooBreitWigner("bw_"+options.channel,"bw_"+options.channel,rrv_mass_lvj,rrv_total_mean_CB,rrv_width_BW);
           model_pdf_sig = RooFFTConvPdf("sig_xww_%s_%s"%(options.channel,options.category),"sigxww_%s_%s"%(options.channel,options.category),rrv_mass_lvj,bw_sig,doubleCB_sig);
           model_pdf_sig.setBufferFraction(1.0);


           model_pdf_sig.plotOn(mplot_sig,RooFit.Name("total_MC"),RooFit.Normalization(rrv_number_signal.getVal()/(6.25*rdataset4fit_mlvj.sumEntries())),RooFit.DrawOption("L"), RooFit.LineColor(kBlue), RooFit.VLines(),RooFit.LineWidth(2),RooFit.LineStyle(1));

           mplot_pull = get_pull(rrv_mass_lvj,mplot_sig);
           mplot_sig.GetYaxis().SetRangeUser(1e-3,mplot_sig.GetMaximum()*1.2);

           draw_canvas_with_pull(mplot_sig,mplot_pull,"plots_signal_width/","signal_width_plot_%d_W%d"%(mass[iMass],mass[iMass]*gammaVal),"",0);
                   
        new_file.Close();
                          
    iMass = iMass +1 ;                              

  os.system("rm %s/%s/temp_datacard_list.txt"%(options.inPath,options.datacardPath));
  os.system("rm %s/%s/temp_workspace_list.txt"%(options.inPath,options.datacardPath));
