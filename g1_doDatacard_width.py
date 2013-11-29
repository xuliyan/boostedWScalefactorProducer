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
          elif line_old.find("lumi_8TeV")!=-1:
            datacardfile.write("CMS_eff_width lnN 1.15 - - - -\n");
            datacardfile.write(line_old);            
          else:  
            datacardfile.write(line_old);
    idatacard = idatacard+1;

  ### make the new workspaces --> loop on the old, create empty new ones  
  for workspace in list_of_workspace :
    iMass = 0 ;
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

        ### copy all the datasets
        getattr(new_workspace,"import")(old_workspace.data(datasetname+"_xww_"+options.channel+"_"+options.category));

        ### make the Breit Wigner core
        RooDoubleCrystalBall.SetName("model_pdf_DoubleCB_BulkG_WW_"+options.channel+"_"+options.category+"_mlvj"); 

        rrv_mean_BW = RooRealVar("rrv_mean_BW_BulkG_WW_"+options.channel+"_"+options.category,"rrv_mean_BW_BulkG_WW_"+options.channel+"_"+options.category,0.);
        rrv_width_BW = RooRealVar("rrv_width_BW_BulkG_WW_"+options.channel+"_"+options.category,"rrv_width_BW_BulkG_WW_"+options.channel+"_"+options.category,mass[iMass]*gammaVal);

        rrv_mean_BW.setConstant(kTRUE);
        rrv_width_BW.setConstant(kTRUE);

        bw = RooBreitWigner("bw_BulkWW_xww_"+options.channel,"bw_BulkG_WW_"+options.channel,old_workspace.var("rrv_mass_lvj"),rrv_mean_BW,rrv_width_BW);

        ### FFT ConvPdf
        model_pdf = RooFFTConvPdf("BulkWW_xww_%s_%s"%(options.channel,options.category),"BulkWW_xww_%s_%s"%(options.channel,options.category),old_workspace.var("rrv_mass_lvj"),RooDoubleCrystalBall,bw);
        model_pdf.setBufferFraction(1.0)
        getattr(new_workspace,"import")(model_pdf);
                          
        iMass = iMass +1 ;                              

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

                                                                                  
        new_workspace.writeToFile(new_file.GetName());
        new_file.Close();


  os.system("rm %s/%s/temp_datacard_list.txt"%(options.inPath,options.datacardPath));
  os.system("rm %s/%s/temp_workspace_list.txt"%(options.inPath,options.datacardPath));
