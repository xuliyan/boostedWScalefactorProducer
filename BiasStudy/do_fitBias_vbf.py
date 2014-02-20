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

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit,RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet,RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, RooMCStudy, RooGlobalFunc,RooChi2MCSModule, RooCurve, RooBernstein, RooHist

############################################
#              Job steering                #
############################################

parser = OptionParser()
parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-j','--njets',     help='number of jets: 1 or 2 ,default:single' , type=int, default = 1) 
parser.add_option('-p','--category',  help='purity category: LP (low purity) or HP (high purity) or NP (no purity selection), default:HP',type="string", default = "HP")
parser.add_option('-l','--channel',   help='lepton flavor: el or mu  or both , default:both' ,type="string", default ="mu" ) ## lepton flavour
parser.add_option('-f','--inPath',    help='directory with workspace' , default = "./" )
parser.add_option('-m','--mass',      help='test signal yield for this mass', type=int, default=-1)
parser.add_option('-n','--nexp',      help='number of toys', type=int, default=1000)
parser.add_option('-g','--fgen',      help='function to generate toys Exp,ExpTail,Pow2,ExpN)', type="string", default="ExpN")
parser.add_option('-r','--fres',      help='function to fit toys (Exp,ExpTail,Pow2,ExpN)',     type="string", default="ExpN")
parser.add_option('-s','--storeplot', help='in case of more than 10 toys just 1/3 stored, more than 100 1/10',     type=int, default=0)
parser.add_option('-z','--isMC',      help='options to run pure mc w+jets toys', type=int, default=0)
parser.add_option('-d','--pseudodata', help='use pseudodata instead of real data', type=int, default=0)
parser.add_option('-t','--shapetest',  help='make W+jets and data fit with different parametrization', type=int, default=0)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")


from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooPow3Pdf, RooErfPow3Pdf


class doBiasStudy_mlvj:

    def __init__(self,in_channel,in_ggH_sample,in_mlvj_min=400., in_mlvj_max=1400., in_mj_min=40, in_mj_max=130, generation_model="ErfExp", fit_model="ErfExp", input_workspace=None):

        self.setTDRStyle();

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        print "###################### construnctor ############################# ";
        ### set the channel type --> electron or muon
        self.channel = in_channel;

        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;

        self.BinWidth_mj = 5.;
        nbins_mj = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max = in_mj_min+nbins_mj*self.BinWidth_mj;

        self.BinWidth_mlvj = 50.;
        self.in_mlvj_min = in_mlvj_min;
        self.nbins_mlvj  = int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        self.in_mlvj_max = in_mlvj_min+self.nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j","pruned m_{J}",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV/c^{2}");
        rrv_mass_j.setBins(nbins_mj);

        self.mj_sideband_lo_min = in_mj_min;
        self.mj_sideband_lo_max = 65;
        self.mj_signal_min      = 65;
        self.mj_signal_max      = 105;
        self.mj_sideband_hi_min = 105;
        self.mj_sideband_hi_max = in_mj_max;
                                                                                            
        ### define invariant mass WW variable
        rrv_mass_lvj = RooRealVar("rrv_mass_lvj","m_{WW}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV/c^{2}");
        rrv_mass_lvj.setBins(self.nbins_mlvj);

        ### define shapes to be used to generate and fit
        self.generation_model = generation_model ;
        self.fit_model = fit_model ;

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4bias_ = RooWorkspace("workspace4bias_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model),"workspace4bias_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model));
        else:
            self.workspace4bias_ = input_workspace;
            
        getattr(self.workspace4bias_,"import")(rrv_mass_lvj);
        getattr(self.workspace4bias_,"import")(rrv_mass_j);

        ## zone definition in the jet mass
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

                                                                           
        ## set the signal sample
        self.file_Directory = options.inPath+"/trainingtrees_%s/"%(self.channel);
        
        self.ggH_sample = in_ggH_sample;
        if in_ggH_sample=="ggH600":  self.vbfhiggs_sample="vbfH600";
        if in_ggH_sample=="ggH700":  self.vbfhiggs_sample="vbfH700";
        if in_ggH_sample=="ggH800":  self.vbfhiggs_sample="vbfH800";
        if in_ggH_sample=="ggH900":  self.vbfhiggs_sample="vbfH900";
        if in_ggH_sample=="ggH1000": self.vbfhiggs_sample="vbfH1000";
        if in_ggH_sample=="ggH1500": self.vbfhiggs_sample="vbfH1500";
        if in_ggH_sample=="ggH2000": self.vbfhiggs_sample="vbfH2000";

        if options.pseudodata :
           self.file_data  = "ofile_pseudodata4higgs.root";                                                                                
        else:  
           self.file_data  = "ofile_data.root";

        self.file_ggH   = ("ofile_%s.root"%(self.ggH_sample));
        self.file_vbfH  = ("ofile_%s.root"%(self.vbfhiggs_sample));

        #WJets0 is the default PS model, WJets1 is the alternative PS model                                                                                                                
        self.file_WJets0_mc = ("ofile_WJets_exclusive_Pythia.root");
        self.file_WJets1_mc = ("ofile_WJets_Herwig.root");

        self.file_VV_mc     = ("ofile_VV.root");# WW+WZ                                                                                                                                        
        self.file_WW_EWK_mc = ("ofile_WW2jet_phantom.root");# WW+WZ                                                                                                                        
        self.file_TTbar_mc         = ("ofile_TTbar_Powheg.root");
        #self.file_TTbar_mc         = ("ofile_TTbar_mcanlo.root");
        self.file_TTbar_matchDn_mc = ("ofile_TTbar_matchDn.root");
        self.file_TTbar_matchUp_mc = ("ofile_TTbar_matchUp.root");
        self.file_TTbar_scaleDn_mc = ("ofile_TTbar_scaleDn.root");
        self.file_TTbar_scaleUp_mc = ("ofile_TTbar_scaleUp.root");
        self.file_TTbar_mcanlo_mc  = ("ofile_TTbar_mcanlo.root");
        self.file_STop_mc          = ("ofile_STop.root");#single Top                                                                                                                            
        ## event categorization as a function of the label        
        if self.wtagger_label=="HP" :

            if self.channel=="el":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="mu":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="em":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;

        if self.wtagger_label=="LP":
                self.wtagger_cut=0.75 ;    self.wtagger_cut_min=0.5 ;

        if self.wtagger_label=="NP":
                self.wtagger_cut=10000;

        ## color palette
        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : 2,
            'VV' : 4,
            'WW_EWK' : 6,
            'STop' : 7,
            'TTbar' : 210,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        ## for basic selection
        self.vpt_cut   = 200;
        self.pfMET_cut = 50;
        self.lpt_cut   = 30;

        if self.channel=="el":
            self.pfMET_cut = 70; self.lpt_cut = 35;#very tight
        self.deltaPhi_METj_cut =2.0;

        self.top_veto_had_max = 200 ;
        self.top_veto_lep_max = 200 ;
        self.top_veto_had_min = 0 ;
        self.top_veto_lep_min = 0 ;

        self.dEta_cut = 2.5 ;
        self.Mjj_cut  = 250 ;

        if self.channel=="mu" and self.wtagger_label=="HP":
           if options.pseudodata == 1:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
             self.rrv_wtagger_eff_reweight_forT.setError(0.338);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           else:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.128);
             self.rrv_wtagger_eff_reweight_forT.setError(0.338);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());


        if self.channel=="el" and self.wtagger_label=="HP":
           if options.pseudodata == 1:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           else:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.836);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());

        if self.channel=="em" and self.wtagger_label=="HP":
           if options.pseudodata == 1:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           else:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.98);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());

        
        ### create output txt file for the shape test 
        if options.isMC == 1 and options.shapetest == 1: 
         self.file_FTestFile_txt = "FTest_%s_%s_mc.txt"%(self.channel,self.wtagger_label);
         self.file_out_FTest = open(self.file_FTestFile_txt,"w");
        elif options.shapetest == 1 and options.isMC == 0 : 
         self.file_FTestFile_txt = "FTest_%s_%s_data.txt"%(self.channel,self.wtagger_label);
         self.file_out_FTest = open(self.file_FTestFile_txt,"w");

        ### create output root file for the pull plot vs mass
        if options.shapetest  == 0:
           self.outputFilePull = ROOT.TFile("pull_output_%s_%s_%s.root"%(self.ggH_sample,options.fgen,options.fres),"RECREATE");
           self.outputFilePull.cd(); 

    ## Set basic TDR style for canvas, pad ..etc ..
    def setTDRStyle(self):
        self.tdrStyle =TStyle("tdrStyle","Style for P-TDR");
        #For the canvas:
        self.tdrStyle.SetCanvasBorderMode(0);
        self.tdrStyle.SetCanvasColor(kWhite);
        self.tdrStyle.SetCanvasDefH(600); #Height of canvas
        self.tdrStyle.SetCanvasDefW(600); #Width of canvas
        self.tdrStyle.SetCanvasDefX(0); #POsition on screen
        self.tdrStyle.SetCanvasDefY(0);
      
        #For the Pad:
        self.tdrStyle.SetPadBorderMode(0);
        self.tdrStyle.SetPadColor(kWhite);
        self.tdrStyle.SetPadGridX(False);
        self.tdrStyle.SetPadGridY(False);
        self.tdrStyle.SetGridColor(0);
        self.tdrStyle.SetGridStyle(3);
        self.tdrStyle.SetGridWidth(1);
      
        #For the frame:
        self.tdrStyle.SetFrameBorderMode(0);
        self.tdrStyle.SetFrameBorderSize(1);
        self.tdrStyle.SetFrameFillColor(0);
        self.tdrStyle.SetFrameFillStyle(0);
        self.tdrStyle.SetFrameLineColor(1);
        self.tdrStyle.SetFrameLineStyle(1);
        self.tdrStyle.SetFrameLineWidth(1);
      
        #For the histo:
        self.tdrStyle.SetHistLineColor(1);
        self.tdrStyle.SetHistLineStyle(0);
        self.tdrStyle.SetHistLineWidth(1);
        self.tdrStyle.SetEndErrorSize(2);
        self.tdrStyle.SetErrorX(0.);
        self.tdrStyle.SetMarkerStyle(20);
      
        #For the fit/function:
        self.tdrStyle.SetOptFit(1);
        self.tdrStyle.SetFitFormat("5.4g");
        self.tdrStyle.SetFuncColor(2);
        self.tdrStyle.SetFuncStyle(1);
        self.tdrStyle.SetFuncWidth(1);
      
        #For the date:
        self.tdrStyle.SetOptDate(0);
      
        #For the statistics box:
        self.tdrStyle.SetOptFile(0);
        self.tdrStyle.SetOptStat(111111); #To display the mean and RMS:
        self.tdrStyle.SetStatColor(kWhite);
        self.tdrStyle.SetStatFont(42);
        self.tdrStyle.SetStatFontSize(0.025);
        self.tdrStyle.SetStatTextColor(1);
        self.tdrStyle.SetStatFormat("6.4g");
        self.tdrStyle.SetStatBorderSize(1);
        self.tdrStyle.SetStatH(0.1);
        self.tdrStyle.SetStatW(0.15);
      
        #Margins:
        self.tdrStyle.SetPadTopMargin(0.05);
        self.tdrStyle.SetPadBottomMargin(0.13);
        self.tdrStyle.SetPadLeftMargin(0.18);
        self.tdrStyle.SetPadRightMargin(0.06);
      
        #For the Global title:
        self.tdrStyle.SetOptTitle(0);
        self.tdrStyle.SetTitleFont(42);
        self.tdrStyle.SetTitleColor(1);
        self.tdrStyle.SetTitleTextColor(1);
        self.tdrStyle.SetTitleFillColor(10);
        self.tdrStyle.SetTitleFontSize(0.05);
      
        #For the axis titles:
        self.tdrStyle.SetTitleColor(1, "XYZ");
        self.tdrStyle.SetTitleFont(42, "XYZ");
        self.tdrStyle.SetTitleSize(0.03, "XYZ");
        self.tdrStyle.SetTitleXOffset(0.9);
        self.tdrStyle.SetTitleYOffset(1.5);
      
        #For the axis labels:
        self.tdrStyle.SetLabelColor(1, "XYZ");
        self.tdrStyle.SetLabelFont(42, "XYZ");
        self.tdrStyle.SetLabelOffset(0.007, "XYZ");
        self.tdrStyle.SetLabelSize(0.03, "XYZ");
      
        #For the axis:
        self.tdrStyle.SetAxisColor(1, "XYZ");
        self.tdrStyle.SetStripDecimals(kTRUE);
        self.tdrStyle.SetTickLength(0.03, "XYZ");
        self.tdrStyle.SetNdivisions(510, "XYZ");
        self.tdrStyle.SetPadTickX(1); #To get tick marks on the opposite side of the frame
        self.tdrStyle.SetPadTickY(1);
      
        #Change for log plots:
        self.tdrStyle.SetOptLogx(0);
        self.tdrStyle.SetOptLogy(0);
        self.tdrStyle.SetOptLogz(0);
      
        #Postscript options:
        self.tdrStyle.SetPaperSize(20.,20.);
        self.tdrStyle.cd();

    #### Method to make a RooAbsPdf giving label, model name, spectrum, if it is mc or not and a constraint list for the parameters
    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],ismc = 0):

        if TString(mass_spectrum).Contains("_mj"):   rrv_x = self.workspace4bias_.var("rrv_mass_j");
        if TString(mass_spectrum).Contains("_mlvj"): rrv_x = self.workspace4bias_.var("rrv_mass_lvj");
        if in_model_name == "CB_v1":
            label_tstring = TString(label);

            if label_tstring.Contains("ggH600"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,610,580,650);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,65,60,90);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1.2,-2,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,2.8,0.1,4);
            elif label_tstring.Contains("vbfH600"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,600,550,650);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,60,50,100);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1,-2,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,10,0.1,15);
            elif label_tstring.Contains("ggH700"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,700,660,740);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,100,80,120);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1.5,-2,-0.7);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,3.,1.5,5);
            elif label_tstring.Contains("vbfH700"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,700,650,750);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,70,40,120);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1,-2,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,8.,0.1,20);
            elif label_tstring.Contains("ggH800"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,800,750,850);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,130,110,150);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1.5,-3,-1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,9.,2,15);
            elif label_tstring.Contains("vbfH800"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,800,750,850);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,90,60,150);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-3,-4,-1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,8.,0.1,30);
            elif label_tstring.Contains("ggH900"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,900,800,920);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,170,130,190);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1.3,-2.5,-0.5);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,13,10,20);
            elif label_tstring.Contains("vbfH900"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,900,850,950);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,160,140,190);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1.3,-2,-0.5);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,25.,15,30);
            elif label_tstring.Contains("ggH1000"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,1000,900,1050);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,180,150,270);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1,-3,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,15.,10,25);
            elif label_tstring.Contains("vbfH1000"):
                rrv_mean_CB  = RooRealVar("rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,"rrv_mean_CB"+label+"_"+self.channel+mass_spectrum,1000,950,1050);
                rrv_sigma_CB = RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,"rrv_sigma_CB"+label+"_"+self.channel+mass_spectrum,220,200,240);
                rrv_alpha_CB = RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,"rrv_alpha_CB"+label+"_"+self.channel+mass_spectrum,-1,-3,-0.1);
                rrv_n_CB     = RooRealVar("rrv_n_CB"+label+"_"+self.channel+mass_spectrum,"rrv_n_CB"+label+"_"+self.channel+mass_spectrum,25.,15,65);   
            else:
                if label_tstring.Contains("600") and (not label_tstring.Contains("1600") ):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 600, 550, 650);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);

                elif label_tstring.Contains("700") and (not label_tstring.Contains("1700") ):
                     rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 700, 600, 800);
                     rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 30,10 ,80);
                     
                elif label_tstring.Contains("800") and (not label_tstring.Contains("1800") ):
                     rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 800, 600, 800);
                     rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 40,10 ,90);
                     
                elif label_tstring.Contains("900") and (not label_tstring.Contains("1900") ):
                    rrv_mean_CB=RooRealVaDr("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 900, 600, 800);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 40,10 ,90);

                elif label_tstring.Contains("1000"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1000, 900,1100);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1100"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1100,1000,1200);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1200"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1200,1100,1300);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1300"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1300,1200,1400);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1400"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1400,1300,1500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1500"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1500,1400,1600);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                    
                elif label_tstring.Contains("1600"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1600,1500,1700);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1700"):
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1700,1500,1800);

                elif label_tstring.Contains("1800"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1800,1500,1900);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("1900"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,1900,1500,2000);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2000"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2000,1800,2200);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2100"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2100,1800,2300);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2200"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2200,1800,2400);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2300"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2300,1800,2500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2400"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2400,1800,2600);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                elif label_tstring.Contains("2500"):
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,2500,2000,2700);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);
                else :
                    rrv_mean_CB=RooRealVar("rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_mean_CB"+label+"_"+self.channel+"_"+self.wtagger_label,700,550,2500);
                    rrv_sigma_CB=RooRealVar("rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_sigma_CB"+label+"_"+self.channel+"_"+self.wtagger_label, 50,20 ,120);

                rrv_alpha_CB=RooRealVar("rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_alpha_CB"+label+"_"+self.channel+"_"+self.wtagger_label,4,1,5);
                rrv_n_CB=RooRealVar("rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,"rrv_n_CB"+label+"_"+self.channel+"_"+self.wtagger_label,20.,10,40);

            model_pdf = RooCBShape("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean_CB,rrv_sigma_CB,rrv_alpha_CB,rrv_n_CB);

        #######################################################
        ############### Simple Falling models #################
        #######################################################

        ## single exponential e^{-[0]*x}                                                                                                                             
        if in_model_name == "Exp" :

            rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel,-0.05,-0.1,0.);

            model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp);

        ## levelled exp for W+jets bkg fit --> e^{-[0]-[1]*x}                                                                                                                
        if in_model_name == "ExpTail":

            rrv_s_ExpTail = RooRealVar("rrv_s_ExpTail"+label+"_"+self.channel,"rrv_s_ExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
            rrv_a_ExpTail = RooRealVar("rrv_a_ExpTail"+label+"_"+self.channel,"rrv_a_ExpTail"+label+"_"+self.channel, 1e-1,-1.e-1,1e6);

            model_pdf     = ROOT.RooExpTailPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail, rrv_a_ExpTail);

        ## three parameter exonential --> e^{-[0]-[1]*x-[2]*x*x}                                                                                                                
        if in_model_name == "Exp_v3":

             if TString(label).Contains("sb_lo"):

              rrv_s_Exp = RooRealVar("rrv_s_Exp_v3"+label+"_"+self.channel,"rrv_s_Exp_v3"+label+"_"+self.channel, -0.004,-0.2,0.);
              rrv_a_Exp = RooRealVar("rrv_a_Exp_v3"+label+"_"+self.channel,"rrv_a_Exp_v3"+label+"_"+self.channel, -10,-50.,-0.001);
              rrv_c_Exp = RooRealVar("rrv_c_Exp_v3"+label+"_"+self.channel,"rrv_c_Exp_v3"+label+"_"+self.channel, -2.e-6,-0.001,0.001);
 
             else:
 
              rrv_s_Exp = RooRealVar("rrv_s_Exp_v3"+label+"_"+self.channel,"rrv_s_Exp_v3"+label+"_"+self.channel, -0.001,-0.2,0.);
              rrv_a_Exp = RooRealVar("rrv_a_Exp_v3"+label+"_"+self.channel,"rrv_a_Exp_v3"+label+"_"+self.channel, -35,-50.,-0.001);
              rrv_c_Exp = RooRealVar("rrv_c_Exp_v3"+label+"_"+self.channel,"rrv_c_Exp_v3"+label+"_"+self.channel, -2.e-6,-0.001,0.001);

             model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "TMath::Exp(%s+%s*%s+%s*%s*%s)"%(rrv_a_Exp.GetName(),rrv_x.GetName(),rrv_s_Exp.GetName(), rrv_c_Exp.GetName(),rrv_x.GetName(),rrv_x.GetName()), RooArgList(rrv_x,rrv_a_Exp,rrv_s_Exp,rrv_c_Exp) );

 
        ## two exponential --> exp^{-[0]*x}+[1]*exp^{-[2]*x}
        if in_model_name == "2Exp":

            rrv_c1_Exp = RooRealVar("rrv_c1_2Exp"+label+"_"+self.channel,"rrv_c1_2Exp"+label+"_"+self.channel,-5e-4,-0.1,0.);
            rrv_c2_Exp = RooRealVar("rrv_c2_2Exp"+label+"_"+self.channel,"rrv_c2_2Exp"+label+"_"+self.channel,-0.005,-0.1,0.);

            exp1 = ROOT.RooExponential("exp1"+label+"_"+self.channel+mass_spectrum,"exp1"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c1_Exp);
            exp2 = ROOT.RooExponential("exp2"+label+"_"+self.channel+mass_spectrum,"exp2"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c2_Exp);

            rrv_frac = RooRealVar("rrv_frac_2Exp"+label+"_"+self.channel,"rrv_frac_2Exp"+label+"_"+self.channel,0.2,0.,8.);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp1,exp2),RooArgList(rrv_frac),1);

       ## For mlvj fit -> Pow funtion: [x/sqrt(s)]^{-[0]}                                                                 
        if in_model_name == "Pow" :

            rrv_c = RooRealVar("rrv_c_Pow"+label+"_"+self.channel,"rrv_c_Pow"+label+"_"+self.channel, -5, -20, 0);

            model_pdf = RooPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x, rrv_c );

        ## For mlvj fit -> Pow function: [x/sqrt{s}]^{-[0]-[1]*log(x/sqrt(s)}                                                                                     
        if in_model_name == "Pow2":

            rrv_c0 = RooRealVar("rrv_c0_Pow2"+label+"_"+self.channel,"rrv_c0_Pow2"+label+"_"+self.channel, 5, 0, 20);
            rrv_c1 = RooRealVar("rrv_c1_Pow2"+label+"_"+self.channel,"rrv_c1_Pow2"+label+"_"+self.channel, 0, -5 , 5);

            model_pdf = RooPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1 );

        ## For mlvj fit -> Pow function: [x/sqrt{s}]^{-[0]-[1]*log(x/sqrt(s)-[2]*log(x/sqrt(s))*log(x/sqrt(s))}                                             
        if in_model_name == "Pow3":

            rrv_c0 = RooRealVar("rrv_c0_Pow3"+label+"_"+self.channel,"rrv_c0_Pow3"+label+"_"+self.channel, 5, 0, 20);
            rrv_c1 = RooRealVar("rrv_c1_Pow3"+label+"_"+self.channel,"rrv_c1_Pow3"+label+"_"+self.channel, 0, -5 , 5);
            rrv_c2 = RooRealVar("rrv_c2_Pow3"+label+"_"+self.channel,"rrv_c2_Pow3"+label+"_"+self.channel, 0, -5 , 5);

            model_pdf = RooPow3Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x, rrv_c0, rrv_c1, rrv_c2);
    
        ## For mlvj fit -> 2Pow function: x^{-[0]}+[1]*x^{-[1]}                                             
        if in_model_name == "2Pow":

            rrv_c0 = RooRealVar("rrv_c0_2Pow"+label+"_"+self.channel,"rrv_c0_2Pow"+label+"_"+self.channel, 5, 0, 20);
            rrv_c1 = RooRealVar("rrv_c1_2Pow"+label+"_"+self.channel,"rrv_c1_2Pow"+label+"_"+self.channel, 0, -5 , 5);
            rrv_frac = RooRealVar("rrv_frac_2Pow"+label+"_"+self.channel,"rrv_frac_2Po2"+label+"_"+self.channel,0.5,0.,1.);
              
            model_pdf = RooGenericPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,"TMath::Power(%s,-%s)+%s*TMath::Power(%s,-%s)"%(rrv_x.GetName(),rrv_c0.GetName(),rrv_frac.GetName(),rrv_x.GetName(),rrv_c1.GetName()),RooArgList(rrv_x,rrv_c0,rrv_c1,rrv_frac));

        ##### polynomial functions   
        if in_model_name == "Chebychev_v2":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);

            model_pdf = RooChebychev("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1));

        if in_model_name == "Chebychev_v3":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
 
            model_pdf = RooChebychev("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));

        if in_model_name == "Chebychev_v4":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p3      = RooRealVar("rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
 
            model_pdf = RooChebychev("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));


        if in_model_name == "Bernstein_v3":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, 0.3,  0, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, 0.01, -1,1);
            rrv_p2      = RooRealVar("rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum, 0.01, 0,1);
 
            model_pdf = RooBernstein("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));


        if in_model_name == "Bernstein_v4":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,0,1);
            rrv_p3      = RooRealVar("rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
 
            model_pdf = RooBernstein("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));


        if in_model_name == "Bernstein_v5":#can replace erf*exp                                                                                                                           

            rrv_p0      = RooRealVar("rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p0_Pol"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p1_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p2_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,0,1);
            rrv_p3      = RooRealVar("rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p3_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p4      = RooRealVar("rrv_p4_Pol"+label+"_"+self.channel+mass_spectrum,"rrv_p4_Pol"+label+"_"+self.channel+mass_spectrum, 0.01,0,1);
 
            model_pdf = RooBernstein("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3,rrv_p4));


        ########################################################
        ############### Turn On + Falling part #################
        ########################################################

        if in_model_name == "ErfExp_v1" : #different init-value and range                                                                                                            

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,-0.005,-0.1,0.);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum,450.,420.,600.);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,60.,15.,100.);

            model_pdf  = RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        if in_model_name == "ErfExpTail" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_ErfExpTail"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfExpTail"+label+"_"+self.channel+mass_spectrum,450.,400.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_ErfExpTail"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExpTail"+label+"_"+self.channel+mass_spectrum,70.,35.,100.);

            if TString(label).Contains("signal_region"):
             rrv_s_ExpTail = RooRealVar("rrv_s_ErfExpTail"+label+"_"+self.channel,"rrv_s_ErfExpTail"+label+"_"+self.channel, 500,0.,1e6);
             rrv_a_ExpTail = RooRealVar("rrv_a_ErfExpTail"+label+"_"+self.channel,"rrv_a_ErfExpTail"+label+"_"+self.channel, -1e-1,-1.e2,1.);
            else:
             rrv_s_ExpTail = RooRealVar("rrv_s_ErfExpTail"+label+"_"+self.channel,"rrv_s_ErfExpTail"+label+"_"+self.channel, 500, 0.,1e6);
             rrv_a_ExpTail = RooRealVar("rrv_a_ErfExpTail"+label+"_"+self.channel,"rrv_a_ErfExpTail"+label+"_"+self.channel, -1e-1,-1.e2,1.);

            erf       = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf));

            expTail   = ROOT.RooExpTailPdf("expTail"+label+"_"+self.channel+mass_spectrum,"expTail"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail,rrv_a_ExpTail);

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,expTail);
          
        if in_model_name == "ErfExp_v3" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_ErfExp_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfExp_v3"+label+"_"+self.channel+mass_spectrum,450.,400.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_ErfExp_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExp_v3"+label+"_"+self.channel+mass_spectrum,60.,35.,100.);

            if TString(label).Contains("sb_lo"):

             rrv_s_ExpTail = RooRealVar("rrv_s_ErfExp_v3"+label+"_"+self.channel,"rrv_s_ErfExp_v3"+label+"_"+self.channel, -0.004,-0.2,0.);
             rrv_a_ExpTail = RooRealVar("rrv_a_ErfExp_v3"+label+"_"+self.channel,"rrv_a_ErfExp_v3"+label+"_"+self.channel, -10,-50.,-0.001);
             rrv_c_ExpTail = RooRealVar("rrv_c_ErfExp_v3"+label+"_"+self.channel,"rrv_c_ErfExp_v3"+label+"_"+self.channel, -2.e-6,-0.001,0.001);
            
            else:

             rrv_s_ExpTail = RooRealVar("rrv_s_ErfExp_v3"+label+"_"+self.channel,"rrv_s_ErfExp_v3"+label+"_"+self.channel, -0.001,-0.2,0.);
             rrv_a_ExpTail = RooRealVar("rrv_a_ErfExp_v3"+label+"_"+self.channel,"rrv_a_ErfExp_v3"+label+"_"+self.channel, -35,-50.,-0.001);
             rrv_c_ExpTail = RooRealVar("rrv_c_ErfExp_v3"+label+"_"+self.channel,"rrv_c_ErfExp_v3"+label+"_"+self.channel, -2.e-6,-0.001,0.001);


            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) );

            exp4 = RooGenericPdf("exp4"+label+"_"+self.channel+mass_spectrum,"exp4"+label+"_"+self.channel+mass_spectrum, "TMath::Exp(%s+%s*%s+%s*%s*%s)"%( rrv_a_ExpTail.GetName(),rrv_s_ExpTail.GetName(), rrv_x.GetName(), rrv_c_ExpTail.GetName(),rrv_x.GetName(),rrv_x.GetName()), RooArgList(rrv_x,rrv_s_ExpTail,rrv_a_ExpTail,rrv_c_ExpTail));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,exp4);


        if in_model_name == "Erf2Exp" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_Erf2Exp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_Erf2Exp"+label+"_"+self.channel+mass_spectrum,450.,400.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_Erf2Exp"+label+"_"+self.channel+mass_spectrum,"rrv_width_Erf2Exp"+label+"_"+self.channel+mass_spectrum,60.,40.,100.);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) );

            rrv_c1_Exp = RooRealVar("rrv_c1_Erf2Exp"+label+"_"+self.channel,"rrv_c1_Erf2Exp"+label+"_"+self.channel,-0.006,-0.1,0.);
            rrv_c2_Exp = RooRealVar("rrv_c2_Erf2Exp"+label+"_"+self.channel,"rrv_c2_Erf2Exp"+label+"_"+self.channel,-5e-4,-0.1,0.);

            rrv_frac = RooRealVar("rrv_frac_Erf2Exp"+label+"_"+self.channel,"rrv_frac_Erf2Exp"+label+"_"+self.channel,0.005,0.,0.05);

            Exp = RooGenericPdf("2exp"+label+"_"+self.channel+mass_spectrum,"2exp"+label+"_"+self.channel+mass_spectrum,"TMath::Exp(%s*%s)+%s*TMath::Exp(%s*%s)"%(rrv_x.GetName(),rrv_c1_Exp.GetName(),rrv_frac.GetName(),rrv_x.GetName(),rrv_c2_Exp.GetName()),RooArgList(rrv_x,rrv_c1_Exp,rrv_c2_Exp,rrv_frac));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,Exp);

        ##################################################################

        if in_model_name == "ErfPow_v1":#can replace erf*exp                                                                                                                           
            rrv_c      = RooRealVar("rrv_c_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfPow"+label+"_"+self.channel+mass_spectrum, -5,-10,0);
            rrv_offset = RooRealVar("rrv_offset_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPow"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width  = RooRealVar("rrv_width_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPow"+label+"_"+self.channel+mass_spectrum,50,20,100);

            model_pdf  = RooErfPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c,rrv_offset,rrv_width);

        if in_model_name == "ErfPow2_v1":#can replace erf*exp                                                                                                                           

            rrv_c0 = RooRealVar("rrv_c0_ErfPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c0_ErfPow2"+label+"_"+self.channel+mass_spectrum,14,1,30);
            rrv_c1 = RooRealVar("rrv_c1_ErfPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c1_ErfPow2"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_offset = RooRealVar("rrv_offset_ErfPow2"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPow2"+label+"_"+self.channel+mass_spectrum, 450,400,600);
            rrv_width  = RooRealVar("rrv_width_ErfPow2"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPow2"+label+"_"+self.channel+mass_spectrum,60,35,100);
            model_pdf  = RooErfPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        if in_model_name == "ErfPow3_v1":#can replace erf*exp                                                                                                                           

            rrv_c0 = RooRealVar("rrv_c0_ErfPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c0_ErfPow3"+label+"_"+self.channel+mass_spectrum,14,-1,30);
            rrv_c1 = RooRealVar("rrv_c1_ErfPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c1_ErfPow3"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_c2 = RooRealVar("rrv_c3_ErfPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c2_ErfPow3"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_offset = RooRealVar("rrv_offset_ErfPow3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPow3"+label+"_"+self.channel+mass_spectrum, 450,400,600);
            rrv_width  = RooRealVar("rrv_width_ErfPow3"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPow3"+label+"_"+self.channel+mass_spectrum,60,35,100);
            model_pdf  = RooErfPow3Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_c2,rrv_offset,rrv_width);


        if in_model_name == "Erf2Pow":#can replace erf*exp                                                                                                                           

            rrv_offset = RooRealVar("rrv_offset_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPow2"+label+"_"+self.channel+mass_spectrum, 450,400,600);
            rrv_width  = RooRealVar("rrv_width_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPow2"+label+"_"+self.channel+mass_spectrum,60,35,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );

            if TString(label).Contains("signal_region"):
             rrv_c0 = RooRealVar("rrv_c0_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_c0_Erf2Pow"+label+"_"+self.channel+mass_spectrum,11,8.,13);
             rrv_c1 = RooRealVar("rrv_c1_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_c1_Erf2Pow"+label+"_"+self.channel+mass_spectrum,3,2,7);
             rrv_frac = RooRealVar("rrv_frac_Erf2Pow"+label+"_"+self.channel,"rrv_frac_Erf2Pow"+label+"_"+self.channel,0.3,0.1,0.4);
            else :
             rrv_c0 = RooRealVar("rrv_c0_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_c0_Erf2Pow"+label+"_"+self.channel+mass_spectrum,3, 2.,7);
             rrv_c1 = RooRealVar("rrv_c1_Erf2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_c1_Erf2Pow"+label+"_"+self.channel+mass_spectrum,3, 2.,7);
             rrv_frac = RooRealVar("rrv_frac_Erf2Pow"+label+"_"+self.channel,"rrv_frac_Erf2Pow"+label+"_"+self.channel,0.2,0.1,0.3);
              
            Pow = RooGenericPdf("2pow"+label+"_"+self.channel+mass_spectrum,"2pow"+label+"_"+self.channel+mass_spectrum,"TMath::Power(%s,-%s)+%s*TMath::Power(%s,-%s)"%(rrv_x.GetName(),rrv_c0.GetName(),rrv_frac.GetName(),rrv_x.GetName(),rrv_c1.GetName()),RooArgList(rrv_x,rrv_c0,rrv_c1,rrv_frac));

            model_pdf  = RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,Pow);

        ###################################################################

        if in_model_name == "ErfPowExp_v1":#can replace erf*exp                                                                                                                         
            rrv_c0 = RooRealVar("rrv_c0_ErfPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_c0_ErfPowExp"+label+"_"+self.channel+mass_spectrum,13,5,40);
            rrv_c1 = RooRealVar("rrv_c1_ErfPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_c1_ErfPowExp"+label+"_"+self.channel+mass_spectrum, 2,0,4);
            rrv_offset = RooRealVar("rrv_offset_ErfPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPowExp"+label+"_"+self.channel+mass_spectrum, 450,400,600);
            rrv_width  = RooRealVar("rrv_width_ErfPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPowExp"+label+"_"+self.channel+mass_spectrum,50,15,150);
            model_pdf  = RooErfPowExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);


        #################################################################

        if in_model_name == "ErfChebychev_v2":#can replace erf*exp                                                                                                           

            rrv_offset  = RooRealVar("rrv_offset_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );

            rrv_p0      = RooRealVar("rrv_p0_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfChebychev_v2"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
 
            pol = RooChebychev("Chebychev_v2"+label+"_"+self.channel+mass_spectrum,"Chebychev_v2"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);
 
        if in_model_name == "ErfChebychev_v3":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p2_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfChebychev_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooChebychev("Chebychev_v3"+label+"_"+self.channel+mass_spectrum,"Chebychev_v3"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);

        if in_model_name == "ErfChebychev_v4":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p2_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);
            rrv_p3      = RooRealVar("rrv_p3_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p3_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfChebychev_v4"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooChebychev("Chebychev_v4"+label+"_"+self.channel+mass_spectrum,"Chebychev_v4"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);


        #################################################################

        if in_model_name == "ErfBernstein_v3":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum, 1, -1., 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum, -0.01,-2.,2.);
            rrv_p2      = RooRealVar("rrv_p2_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p2_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum, 0.001,-1.,1.);

            rrv_offset  = RooRealVar("rrv_offset_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v3"+label+"_"+self.channel+mass_spectrum,"Bernstein_v3"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);

        if in_model_name == "ErfBernstein_v4":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p2_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);
            rrv_p3      = RooRealVar("rrv_p3_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p3_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfBernstein_v4"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v4"+label+"_"+self.channel+mass_spectrum,"Bernstein_v4"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);

        if in_model_name == "ErfBernstein_v5":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p0_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p1_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p2_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);
            rrv_p3      = RooRealVar("rrv_p3_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p3_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,-1,1);
            rrv_p4      = RooRealVar("rrv_p4_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p4_ErfBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);

            rrv_offset  = RooRealVar("rrv_offset_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfBernstein_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            erf = RooGenericPdf("erf"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::Erf((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v5"+label+"_"+self.channel+mass_spectrum,"Bernstein_v5"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3,rrv_p4));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);


        ########################################################
        ############### Turn On + Falling part #################
        ########################################################

        if in_model_name == "AtanExp_v1" : #different init-value and range                                                                                                            

            rrv_c_AtanExp      = RooRealVar("rrv_c_AtanExp"+label+"_"+self.channel+mass_spectrum,"rrv_c_AtanExp"+label+"_"+self.channel+mass_spectrum,-0.006,-0.1,0.);
            rrv_offset_AtanExp = RooRealVar("rrv_offset_AtanExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanExp"+label+"_"+self.channel+mass_spectrum,450.,350.,600.);
            rrv_width_AtanExp  = RooRealVar("rrv_width_AtanExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanExp"+label+"_"+self.channel+mass_spectrum,60.,15.,100.);

            model_pdf  = RooAtanExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_AtanExp,rrv_offset_AtanExp,rrv_width_AtanExp);

        if in_model_name == "AtanExpTail" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_AtanExpTail"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanExpTail"+label+"_"+self.channel+mass_spectrum,450.,350.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_AtanExpTail"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanExpTail"+label+"_"+self.channel+mass_spectrum,70.,15.,100.);

            rrv_s_ExpTail = RooRealVar("rrv_s_AtanExpTail"+label+"_"+self.channel,"rrv_s_AtanExpTail"+label+"_"+self.channel, 250,-1.e6,1e6);
            rrv_a_ExpTail = RooRealVar("rrv_a_AtanExpTail"+label+"_"+self.channel,"rrv_a_AtanExpTail"+label+"_"+self.channel, 5e-1,-1.e2,1e6);

            atan       = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf));

            expTail   = ROOT.RooExpTailPdf("expTail"+label+"_"+self.channel+mass_spectrum,"expTail"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_s_ExpTail,rrv_a_ExpTail);

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,expTail);
          
        if in_model_name == "AtanExp_v3" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_AtanExp_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanExp_v3"+label+"_"+self.channel+mass_spectrum,450.,350.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_AtanExp_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanExp_v3"+label+"_"+self.channel+mass_spectrum,70.,15.,100.);

            rrv_s_ExpTail = RooRealVar("rrv_s_AtanExp_v3"+label+"_"+self.channel,"rrv_s_AtanExp_v3"+label+"_"+self.channel, -0.006,-0.2,0.);

            rrv_a_ExpTail = RooRealVar("rrv_a_AtanExp_v3"+label+"_"+self.channel,"rrv_a_AtanExp_v3"+label+"_"+self.channel, -0.1,-50.,-0.001);
            rrv_c_ExpTail = RooRealVar("rrv_c_AtanExp_v3"+label+"_"+self.channel,"rrv_c_AtanExp_v3"+label+"_"+self.channel, -0.00001,-0.001,0.);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) );

            exp4 = RooGenericPdf("exp4"+label+"_"+self.channel+mass_spectrum,"exp4"+label+"_"+self.channel+mass_spectrum, "TMath::Exp(%s+%s*%s+%s*%s*%s)"%( rrv_a_ExpTail.GetName(),rrv_s_ExpTail.GetName(), rrv_x.GetName(), rrv_c_ExpTail.GetName(),rrv_x.GetName(),rrv_x.GetName()), RooArgList(rrv_x,rrv_s_ExpTail,rrv_a_ExpTail,rrv_c_ExpTail));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,exp4);


        if in_model_name == "Atan2Exp" : #different init-value and range                                                                                                             
   
            rrv_offset_erf = RooRealVar("rrv_offset_Atan2Exp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_Atan2Exp"+label+"_"+self.channel+mass_spectrum,450.,350.,600.);
            rrv_width_erf  = RooRealVar("rrv_width_Atan2Exp"+label+"_"+self.channel+mass_spectrum,"rrv_width_Atan2Exp"+label+"_"+self.channel+mass_spectrum,70.,15.,100.);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf) );

            rrv_c1_Exp = RooRealVar("rrv_c1_Atan2Exp"+label+"_"+self.channel,"rrv_c1_Atan2Exp"+label+"_"+self.channel,-0.05,-0.1,0.);
            rrv_c2_Exp = RooRealVar("rrv_c2_Atan2Exp"+label+"_"+self.channel,"rrv_c2_Atan2Exp"+label+"_"+self.channel,-0.05,-0.1,0.);

            exp1 = ROOT.RooExponential("exp1"+label+"_"+self.channel+mass_spectrum,"exp1"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c1_Exp);
            exp2 = ROOT.RooExponential("exp2"+label+"_"+self.channel+mass_spectrum,"exp3"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c2_Exp);

            rrv_frac = RooRealVar("rrv_frac_Atan2Exp"+label+"_"+self.channel,"rrv_frac_Atan2Exp"+label+"_"+self.channel,0.5,0.,1.);

            Exp = RooAddPdf("2Exp"+label+"_"+self.channel+mass_spectrum,"2Exp"+label+"_"+self.channel+mass_spectrum,RooArgList(exp1,exp2),RooArgList(rrv_frac));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,Exp);

        ##################################################################

        if in_model_name == "AtanPow_v1":#can replace erf*exp                                                                                 
                                          
            rrv_c      = RooRealVar("rrv_c_AtanPow"+label+"_"+self.channel+mass_spectrum,"rrv_c_AtanPow"+label+"_"+self.channel+mass_spectrum, -5,-10,0);
            rrv_offset = RooRealVar("rrv_offset_AtanPow"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanPow"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width  = RooRealVar("rrv_width_AtanPow"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanPow"+label+"_"+self.channel+mass_spectrum,50,20,100);

            model_pdf  = RooAtanPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c,rrv_offset,rrv_width);

        if in_model_name == "AtanPow2_v1":#can replace erf*exp                                                                                                                           

            rrv_c0 = RooRealVar("rrv_c0_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c0_AtanPow2"+label+"_"+self.channel+mass_spectrum,14,1,30);
            rrv_c1 = RooRealVar("rrv_c1_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c1_AtanPow2"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_offset = RooRealVar("rrv_offset_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanPow2"+label+"_"+self.channel+mass_spectrum, 600,400,600);
            rrv_width  = RooRealVar("rrv_width_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanPow2"+label+"_"+self.channel+mass_spectrum,60,10,100);
            model_pdf  = RooAtanPow2Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        if in_model_name == "AtanPow3_v1":#can replace erf*exp                                                                                                                           

            rrv_c0 = RooRealVar("rrv_c0_AtanPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c0_AtanPow3"+label+"_"+self.channel+mass_spectrum,14,1,30);
            rrv_c1 = RooRealVar("rrv_c1_AtanPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c1_AtanPow3"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_c2 = RooRealVar("rrv_c3_AtanPow3"+label+"_"+self.channel+mass_spectrum,"rrv_c2_AtanPow3"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_offset = RooRealVar("rrv_offset_AtanPow3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanPow3"+label+"_"+self.channel+mass_spectrum, 600,400,600);
            rrv_width  = RooRealVar("rrv_width_AtanPow3"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanPow3"+label+"_"+self.channel+mass_spectrum,60,10,100);
            model_pdf  = RooAtanPow3Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_c2,rrv_offset,rrv_width);

        if in_model_name == "Atan2Pow":#can replace erf*exp                                                                                                                           

            rrv_offset = RooRealVar("rrv_offset_Atan2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanPow2"+label+"_"+self.channel+mass_spectrum, 600,400,600);
            rrv_width  = RooRealVar("rrv_width_Atan2Pow"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanPow2"+label+"_"+self.channel+mass_spectrum,60,10,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%(rrv_x.GetName(),rrv_offset_erf.GetName(), rrv_width_erf.GetName()), RooArgList(rrv_x,rrv_offset_erf,rrv_width_erf));

            rrv_c0 = RooRealVar("rrv_c0_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c0_AtanPow2"+label+"_"+self.channel+mass_spectrum,14,1,30);
            rrv_c1 = RooRealVar("rrv_c1_AtanPow2"+label+"_"+self.channel+mass_spectrum,"rrv_c1_AtanPow2"+label+"_"+self.channel+mass_spectrum, 5,-5,10);
            rrv_frac = RooRealVar("rrv_frac_Atan2Pow"+label+"_"+self.channel,"rrv_frac_Atan2Pow"+label+"_"+self.channel,0.5,0.,1.);
              
            Pow = RooGenericPdf("2pow"+label+"_"+self.channel+mass_spectrum,"2pow"+label+"_"+self.channel+mass_spectrum,"TMath::Pow(%s,-%s)+%s*TMath::Pow(%s,-%s)"(rrv_x.GetName(),rrv_c0.GetName(),rrv_frac.GetName(),rrv_x.GetName(),rrv_c1.GetName()),RooArgList(rrv_x,rrv_c0,rrv_c1,rrv_frac));

            model_pdf  = RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,Pow);

        ###################################################################

        if in_model_name == "AtanPowExp_v1":#can replace erf*exp                                                                                                                         

            rrv_c0 = RooRealVar("rrv_c0_AtanPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_c0_AtanPowExp"+label+"_"+self.channel+mass_spectrum,13,5,40);
            rrv_c1 = RooRealVar("rrv_c1_AtanPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_c1_AtanPowExp"+label+"_"+self.channel+mass_spectrum, 2,0,4);
            rrv_offset = RooRealVar("rrv_offset_AtanPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanPowExp"+label+"_"+self.channel+mass_spectrum, 450,400,600);
            rrv_width  = RooRealVar("rrv_width_AtanPowExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanPowExp"+label+"_"+self.channel+mass_spectrum,50,15,150);
            model_pdf  = RooAtanPowExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c0,rrv_c1,rrv_offset,rrv_width);

        #################################################################

        if in_model_name == "AtanChebychev_v2":#can replace erf*exp                                                                                                           

            rrv_offset  = RooRealVar("rrv_offset_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );

            rrv_p0      = RooRealVar("rrv_p0_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanChebychev_v2"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
 
            pol = RooChebychev("Chebychev_v2"+label+"_"+self.channel+mass_spectrum,"Chebychev_v2"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,pol);
 
        if in_model_name == "AtanChebychev_v3":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p2_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanChebychev_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width));
 
            pol = RooChebychev("Chebychev_v3"+label+"_"+self.channel+mass_spectrum,"Chebychev_v3"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,pol);

        if in_model_name == "AtanChebychev_v4":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.3, -3, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p2_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);
            rrv_p3      = RooRealVar("rrv_p3_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p3_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum, -0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanChebychev_v4"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooChebychev("Chebychev_v4"+label+"_"+self.channel+mass_spectrum,"Chebychev_v4"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,erf,pol);

        #################################################################

        if in_model_name == "AtanBernstein_v3":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_p2_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);

            rrv_offset  = RooRealVar("rrv_offset_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v3"+label+"_"+self.channel+mass_spectrum,"Bernstein_v3"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,pol);

        if in_model_name == "AtanBernstein_v4":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p2_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);
            rrv_p3      = RooRealVar("rrv_p3_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_p3_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum, 0.001,-1,1);

            rrv_offset  = RooRealVar("rrv_offset_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanBernstein_v4"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"atan"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v4"+label+"_"+self.channel+mass_spectrum,"Bernstein_v4"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,pol);

        if in_model_name == "AtanBernstein_v5":#can replace erf*exp                                                                                               

            rrv_p0      = RooRealVar("rrv_p0_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p0_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.3, 0, 3);
            rrv_p1      = RooRealVar("rrv_p1_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p1_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.01,-1,1);
            rrv_p2      = RooRealVar("rrv_p2_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p2_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);
            rrv_p3      = RooRealVar("rrv_p3_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p3_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,-1,1);
            rrv_p4      = RooRealVar("rrv_p4_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum,"rrv_p4_AtanBernstein_v5"+label+"_"+self.channel+mass_spectrum, 0.001,0,1);

            rrv_offset  = RooRealVar("rrv_offset_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_offset_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum, 450,350,600);
            rrv_width   = RooRealVar("rrv_width_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,"rrv_width_AtanBernstein_v3"+label+"_"+self.channel+mass_spectrum,50,20,100);

            atan = RooGenericPdf("atan"+label+"_"+self.channel+mass_spectrum,"erf"+label+"_"+self.channel+mass_spectrum, "(1.+TMath::ATan((%s-%s)/%s))/2."%( rrv_x.GetName(),rrv_offset.GetName(), rrv_width.GetName()), RooArgList(rrv_x,rrv_offset,rrv_width) );
 
            pol = RooBernstein("Bernstein_v5"+label+"_"+self.channel+mass_spectrum,"Bernstein_v5"+label+"_"+self.channel+mass_spectrum,rrv_x,RooArgList(rrv_p0,rrv_p1,rrv_p2,rrv_p3,rrv_p4));

            model_pdf = ROOT.RooProdPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,atan,pol);

        ###################################################################

        if in_model_name == "Keys":

            if TString(label).Contains("sb_lo"): 
               rdataset = self.workspace4bias_.data("rdataset_WJets0_sb_lo_%s_mlvj"%(self.channel))
            elif TString(label).Contains("signal_region"): 
               rdataset = self.workspace4bias_.data("rdataset_WJets0_signal_region_%s_mlvj"%(self.channel))
        
            model_pdf = RooKeysPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum, rrv_x,rdataset,RooKeysPdf.MirrorRight);
                                           
        ## return the pdf
        getattr(self.workspace4bias_,"import")(model_pdf)
        return self.workspace4bias_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)

    ### Define the Extended Pdf for and mJ fit giving: label, fit model name, list constraint and ismc
    def make_Model(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[], ismc_wjet=0, area_init_value=500):

      ##### define an extended pdf from a standard Roofit One
      print " "
      print "###############################################"
      print "## Make model : ",label," ",in_model_name,"##";
      print "###############################################"
      print " "

      rrv_number = RooRealVar("rrv_number"+label+"_"+self.channel+mass_spectrum,"rrv_number"+label+"_"+self.channel+mass_spectrum,500.,0.,1e5);
      ## call the make RooAbsPdf method
      model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList,ismc_wjet)
      print "######## Model Pdf ########"
      model_pdf.Print();
      
      ## create the extended pdf
      model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
      print "######## Model Extended Pdf ########"

      #### put all the parameters ant the shape in the workspace
      getattr(self.workspace4bias_,"import")(rrv_number)
      getattr(self.workspace4bias_,"import")(model)
      self.workspace4bias_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();

      ## return the total extended pdf
      return self.workspace4bias_.pdf("model"+label+"_"+self.channel+mass_spectrum);
  
    #### get a generic mlvj model from the workspace
    def get_mlvj_Model(self,label, mlvj_region, model):
        print "model"+label+mlvj_region+model+"_"+self.channel+"_mlvj"
        return self.workspace4bias_.pdf("model"+label+mlvj_region+model+"_"+self.channel+"_mlvj");

    #### get a general mlvj model and fiz the paramters --> for extended pdf
    def get_General_mlvj_Model(self, label, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing amd return a general mlvj model ############"
        rdataset_General_mlvj = self.workspace4bias_.data("rdataset4bias%s%s_%s_mlvj"%(label, mlvj_region,self.channel))
        model_General = self.get_mlvj_Model(label,mlvj_region,model);
        rdataset_General_mlvj.Print();
        model_General.Print();
        parameters_General = model_General.getParameters(rdataset_General_mlvj);
        par=parameters_General.createIterator(); par.Reset();
        param=par.Next()
        while (param):
          param.setConstant(kTRUE);
          param.Print();
          param=par.Next()
        return model_General;


    ###### get TTbar model mlvj in a region
    def get_TTbar_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing TTbar mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_TTbar",mlvj_region,model);

    ###### get Single Top model mlvj in a region
    def get_STop_mlvj_Model(self,  mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_STop",mlvj_region,model);

    ###### get ggH model mlvj in a region
    def get_ggH_signal_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_%s"%(self.ggH_sample),mlvj_region,model);

    ###### get ggH model mlvj in a region
    def get_vbfH_signal_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_%s"%(self.vbfhiggs_sample),mlvj_region,model);

    ###### get VV mlvj in a region
    def get_VV_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_VV",mlvj_region,model);

    ###### get VV mlvj in a region
    def get_WW_EWK_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_WW_EWK",mlvj_region,model);
                                                                                                                    

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj",model="",additional_info="",notExtended=0):
        print "########### Fixing an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
        if notExtended:
           label = "_pdf"+label;    
        model = self.get_mlvj_Model(label,mlvj_region,model+additional_info);
        rdataset.Print();
        model.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); 
        par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### clone Model Extend in a not extend skipping the normalization
    def clone_Model(self, inputPdf, label, mlvj_region="_signal_region",mass_spectrum="_mlvj", model =""):
        print "########### Cloning an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum));
        model = self.get_mlvj_Model(label,mlvj_region,model);                
        rdataset.Print();
        inputPdf.Print();
        model.Print()
        parameters  = inputPdf.getParameters(rdataset);
        parameters2 = model.getParameters(rdataset);
        par         = parameters.createIterator();
        par2        = parameters2.createIterator();
        par.Reset();
        par2.Reset();
        param  = par.Next();
        param2 = par2.Next();

        if parameters.getSize() < parameters2.getSize() :
         while (param):
             if(not TString(param2.GetName()).Contains("number")):
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param.Print(); param2.Print();
                 param  = par.Next();
                 param2 = par2.Next();
             else:
                 param2 = par2.Next();
        elif  parameters.getSize() > parameters2.getSize() :
         while (param2):
             if(not TString(param.GetName()).Contains("number")):
                 param.Print(); param2.Print();
                 param.Print(); param2.Print();
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param  = par.Next();
                 param2 = par2.Next();
             else:
                param = par.Next();

        else:
         while (param):
                 param.Print(); param2.Print(); 
                 param.setVal(param2.getVal());
                 param.setError(param2.getError());
                 param  = par.Next();
                 param2 = par2.Next();
            
    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Deco_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj"):
        print "########### Fixing an Extended Pdf for mlvj ############"
        rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum));
        model = self.workspace4bias_.pdf("model_pdf"+label+mlvj_region+"_"+self.channel+"_mlvj"+"_Deco"+label+mlvj_region+"_"+self.channel+"_"+self.wtagger_label+"_mlvj");
        rdataset.Print();
        model.Print();
        parameters = model.getParameters(rdataset);
        par=parameters.createIterator(); 
        par.Reset();
        param=par.Next()
        while (param):
            param.setConstant(kTRUE);
            param.Print();
            param=par.Next()

    ### Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
    def fit_mlvj_model_single_MC(self,in_file_name, label, in_range, mlvj_model, deco=0, show_constant_parameter=0, logy=0, ismc=0):

        print "############### Fit mlvj single MC sample ",in_file_name," ",label," ",mlvj_model," ",in_range," ##################"
        ## imporparam_generatedt variable and dataset
        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj")
        rdataset     = self.workspace4bias_.data("rdataset4bias"+label+in_range+"_"+self.channel+"_mlvj");
        constrainslist =[];

        ## make the extended pdf model
        model = self.make_Model(label+in_range+mlvj_model,mlvj_model,"_mlvj",constrainslist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4bias_,"import")(rfresult)


        ## plot the result
        mplot = rrv_mass_lvj.frame(RooFit.Title("M_{lvj"+in_range+"} fitted by "+mlvj_model), RooFit.Bins(int(rrv_mass_lvj.getBins())));
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0));
        ## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset, model, rfresult,mplot,2,"L")
        model.plotOn( mplot, RooFit.Name("model_mc") )#, RooFit.VLines()); in order to have the right pull
        rdataset.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Name("data") );

        ## get the pull
        mplot_pull = self.get_pull(rrv_mass_lvj,mplot,"data","model_mc");
        parameters_list = model.getParameters(rdataset);
        
        ##CALCULATE CHI2                                                                                                                                                    
        datahist   = rdataset.binnedClone(rdataset.GetName()+"_binnedClone",rdataset.GetName()+"_binnedClone");
        histo_data = datahist.createHistogram("histo_data",rrv_mass_lvj) ;
        histo_data.SetName("histo_data");
        histo_func = model.createHistogram("histo_func",rrv_mass_lvj) ;
        histo_func.SetName("histo_func");
        
        Nbin     = int(rrv_mass_lvj.getBins());
        rresult_param = rfresult.floatParsFinal();
        nparameters   = rresult_param.getSize();
        ChiSquare = model.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
        chi_over_ndf  = ChiSquare.getVal()/(Nbin-nparameters);

        
        residHist = mplot.residHist("data","model_mc");
        residual = 0. ;
        for iPoint in range(residHist.GetN()):
         x = ROOT.Double(0.); y = ROOT.Double(0) ;
         residHist.GetPoint(iPoint,x,y); 
         residual = residual + y**2 ;
        
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

        if options.shapetest == 1:
         self.file_out_FTest.write(" ###################### \n");
         self.file_out_FTest.write(" Model pdf %s"%(model.GetName()));
         self.file_out_FTest.write(" Chi2 Chi2var %0.2f "%(chi_over_ndf));
         self.file_out_FTest.write(" Residual %0.2f   Nbin %0.2f nparameters %0.2f \n"%(residual,Nbin,nparameters));
        
        ##Add Chisquare to mplot_pull                                                                                                                                                    
        cs2 = TLatex(0.75,0.8,"#chi^{2}/ndf = %0.2f "%(float(chi_over_ndf)));
        cs2.SetNDC();
        cs2.SetTextSize(0.12);
        cs2.AppendPad("same");
        mplot_pull.addObject(cs2);
        
        self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_g1/m_lvj_fitting_%s/"%(options.additioninformation, self.channel,self.wtagger_label,mlvj_model),in_file_name,"m_lvj"+in_range+"_"+mlvj_model, show_constant_parameter, logy);

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region                                                 
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").getVal()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel+in_range+"_mlvj").getVal());
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").setError(self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").getError()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel+in_range+"_mlvj").getVal() );
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").Print();
        
     ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region," ",mlvj_model," ##################"

        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj");
        rdataset_data_mlvj = self.workspace4bias_.data("rdataset_data%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get and fix the minor component shapes in the sb low

        model_VV_backgrounds     = self.get_VV_mlvj_Model("_sb_lo",options.fgen);
        model_STop_backgrounds   = self.get_STop_mlvj_Model("_sb_lo",options.fgen);
        model_TTbar_backgrounds  = self.get_TTbar_mlvj_Model("_sb_lo",options.fgen);
        model_WW_EWK_backgrounds = self.get_WW_EWK_mlvj_Model("_sb_lo",options.fgen);

        number_VV_sb_lo_mlvj      = self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)); ## get the normalization
        number_TTbar_sb_lo_mlvj   = self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)); ## get the normalization
        number_STop_sb_lo_mlvj    = self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)); ## get the normalization
        number_WW_EWK_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)); ## get the normalization       

        self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).Print();

        ### Make the Pdf for the WJets
        model_pdf_WJets = self.make_Pdf("%s_sb_lo%s_from_fitting"%(label,mlvj_model), mlvj_model,"_mlvj");
        model_pdf_WJets.Print();
        ### inititalize the value to what was fitted with the mc in the sideband
        number_WJets_sb_lo = self.workspace4bias_.var("rrv_number%s_sb_lo%s_%s_mlvj"%(label,options.fgen,self.channel)).clone("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+self.channel+"_mlvj");
    
        model_WJets        = RooExtendPdf("model%s_sb_lo%s_from_fitting_%s_mlvj"%(label,mlvj_model,self.channel),"model%s_sb_lo%s_from_fitting_%s_mlvj"%(label,mlvj_model,self.channel),model_pdf_WJets,number_WJets_sb_lo);
        number_WJets_sb_lo.Print();

        ## Add the other bkg component fixed to the total model --> in the extended way
        model_data = RooAddPdf("model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds, model_WW_EWK_backgrounds));
        
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.SumW2Error(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4bias_,"import")(model_data)

        model_WJets.getParameters(rdataset_data_mlvj).Print("v");

        ### data in the sideband plus error from fit
        rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_sb_lo%s_%s_mlvj"%(mlvj_model,self.channel),"rrv_number_data_sb_lo%s_%s_mlvj"%(mlvj_model,self.channel),
                                                 self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+
                                                 self.workspace4bias_.var("rrv_number%s_sb_lo%s_from_fitting_%s_mlvj"%(label,mlvj_model,self.channel)).getVal() );

        rrv_number_data_sb_lo_mlvj.setError(TMath.Sqrt(self.workspace4bias_.var("rrv_number%s_sb_lo%s_from_fitting_%s_mlvj"%(label,mlvj_model,self.channel)).getError()*
                                                       self.workspace4bias_.var("rrv_number%s_sb_lo%s_from_fitting_%s_mlvj"%(label,mlvj_model,self.channel)).getError()+
                                                       self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()*
                                                       self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()+
                                                       self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()*
                                                       self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()+
                                                       self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()*
                                                       self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()+
                                                       self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()*
                                                       self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getError()));

        rrv_number_data_sb_lo_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rrv_number_data_sb_lo_mlvj)

        ### plot for WJets default + default shape
        if label=="_WJets0":

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in M_j sideband "), RooFit.Bins(int(rrv_mass_lvj.getBins())));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Invisible(), RooFit.Name("data_invisible") );

            model_data.plotOn(mplot, RooFit.Components(model_WJets.GetName()+","+model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WJets"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WW_EWK"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()), RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_STop_backgrounds.GetName()), RooFit.Name("STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());


            #solid line
            model_data.plotOn(mplot,  RooFit.Components(model_WJets.GetName()+","+model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WW_EWK_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()), RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_TTbar_backgrounds.GetName()+","+model_STop_backgrounds.GetName()), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_STop_backgrounds.GetName()), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            ### draw the error band
            draw_error_band(rdataset_data_mlvj, model_data,self.workspace4bias_.var("rrv_number_data_sb_lo%s_%s_mlvj"%(mlvj_model,self.channel)) ,rfresult,mplot,self.color_palet["Uncertainty"],"F");
            model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
            model_data.plotOn( mplot , RooFit.Invisible(), RooFit.Name("model_mc"));

            if options.pseudodata == 1:
             rdataset_data_mlvj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Name("data"));               
            else: 
             self.getData_PoissonInterval(rdataset_data_mlvj,mplot);

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
            
            ### Add the legend to the plot
            leg = self.legend4Plot(mplot,0,1,0., 0.06, 0.16, 0.);
            mplot.addObject(leg)

            ### get the pull plot and store the canvas
            mplot_pull = self.get_pull(rrv_mass_lvj,mplot,"data","model_mc");
            parameters_list = model_data.getParameters(rdataset_data_mlvj);

            ##CALCULATE CHI2                                                                                                                                               
            datahist   = rdataset_data_mlvj.binnedClone(rdataset_data_mlvj.GetName()+"_binnedClone",rdataset_data_mlvj.GetName()+"_binnedClone");
            histo_data = datahist.createHistogram("histo_data",rrv_mass_lvj) ;
            histo_data.SetName("histo_data");
            histo_func = model_data.createHistogram("histo_func",rrv_mass_lvj) ;
            histo_func.SetName("histo_func");

            Nbin     = int(rrv_mass_lvj.getBins());
            rresult_param = rfresult.floatParsFinal();
            nparameters   = rresult_param.getSize();
            ChiSquare = model_data.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
            chi_over_ndf  = ChiSquare.getVal()/(Nbin-nparameters);

            residHist = mplot.residHist("data","model_mc");
            residual = 0. ;
            for iPoint in range(residHist.GetN()):
              x = ROOT.Double(0.); y = ROOT.Double(0) ;
              residHist.GetPoint(iPoint,x,y); 
              residual = residual + y**2 ;

            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);

            if options.shapetest == 1:
             self.file_out_FTest.write(" ###################### \n");
             self.file_out_FTest.write(" Model pdf %s"%(model_data.GetName()));
             self.file_out_FTest.write(" Chi2 Chi2var %0.2f "%(chi_over_ndf));
             self.file_out_FTest.write(" Residual %0.2f   Nbin %0.2f nparameters %0.2f \n"%(residual,Nbin,nparameters));
        
            ##Add Chisquare to mplot_pull                                                                                                                                             
            cs2 = TLatex(0.75,0.8,"#chi^{2}/ndf = %0.2f "%(float(chi_over_ndf)));
            cs2.SetNDC();
            cs2.SetTextSize(0.12);
            cs2.AppendPad("same");
            mplot_pull.addObject(cs2);
                
            self.draw_canvas_with_pull( mplot, mplot_pull,parameters_list,"plots_%s_%s_%s_g1/m_lvj_fitting_%s/"%(options.additioninformation, self.channel,self.wtagger_label,mlvj_model), "m_lvj_sb_lo%s_%s"%(label,mlvj_model),"",1,1)
            

    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass="jet_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");

        rrv_mass_j         = self.workspace4bias_.var("rrv_mass_j") 
        rrv_mass_lvj       = self.workspace4bias_.var("rrv_mass_lvj")
        rrv_weight         = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 
        
        #dataset of m_lvj -> before and after vbf cuts -> central object value
        rdataset_sb_lo_mlvj     = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_sb_hi_mlvj         = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        rdataset4bias_sb_lo_mlvj     = RooDataSet("rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4bias_signal_region_mlvj = RooDataSet("rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4bias_sb_hi_mlvj     = RooDataSet("rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        ###### Define the event categorization
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");

        combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4bias = RooDataSet("combData4bias"+label+"_"+self.channel,"combData4bias"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        print "N entries: ", treeIn.GetEntries();

        for i in range(treeIn.GetEntries()):


          if i % 100000 == 0: print "iEvent: ",i
          treeIn.GetEntry(i);

          discriminantCut = False;

          wtagger=-1;
          wtagger=getattr(treeIn,"jet_tau2tau1");

          if wtagger < self.wtagger_cut:
            discriminantCut = True;
          else:
            discriminantCut = False;
                  
          tmp_scale_to_lumi = treeIn.wSampleWeight;
                                
          jet_1 = ROOT.TLorentzVector();
          jet_2 = ROOT.TLorentzVector();

          njet = 0. ; tmp_vbf_dEta =0.; tmp_vbf_Mjj = 0.; ungroomed_jet_pt = 0.; pfMET = 0.; mass_lvj = 0. ;

          # jet mass , central value
          tmp_jet_mass = getattr(treeIn, jet_mass);
          tmp_vbf_dEta = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta")-getattr(treeIn,"vbf_maxpt_j2_eta"));
          tmp_vbf_Mjj  = getattr(treeIn, "vbf_maxpt_jj_m");
          njet         = getattr(treeIn,"numberJetBin");
          ungroomed_jet_pt = getattr(treeIn,"ungroomed_jet_pt");
          pfMET    = getattr(treeIn,"pfMET");
          mass_lvj = getattr(treeIn,"mass_lvj_type0_met");

          isFullVBF = 0 ;

          if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2:
            isFullVBF = 1 ;
          
          if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
            isFullVBF = 2 ;
   
          tmp_event_weight = 0 ;       
          tmp_event_weight4bias = 0 ;

          if isFullVBF !=0 :
              
             tmp_event_weight      = getattr(treeIn,"totalEventWeight");
             tmp_event_weight4bias = getattr(treeIn,"eff_and_pu_Weight");
             tmp_interference_weight_H600  = getattr(treeIn,"interference_Weight_H600");
             tmp_interference_weight_H700  = getattr(treeIn,"interference_Weight_H700");
             tmp_interference_weight_H800  = getattr(treeIn,"interference_Weight_H800");
             tmp_interference_weight_H900  = getattr(treeIn,"interference_Weight_H900");
             tmp_interference_weight_H1000 = getattr(treeIn,"interference_Weight_H1000");
             
             if TString(label).Contains("ggH600"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H600
                tmp_event_weight4bias = tmp_event_weight4bias*tmp_interference_weight_H600
             if TString(label).Contains("ggH700"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H700
                tmp_event_weight4bias = tmp_event_weight4bias*tmp_interference_weight_H700
             if TString(label).Contains("ggH800"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H800
                tmp_event_weight4bias = tmp_event_weight4bias*tmp_interference_weight_H800
             if TString(label).Contains("ggH900"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H900
                tmp_event_weight4bias = tmp_event_weight4bias*tmp_interference_weight_H900
             if TString(label).Contains("ggH1000"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H1000
                tmp_event_weight4bias = tmp_event_weight4bias*tmp_interference_weight_H1000

             if TString(label).Contains("vbfH600"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H600;
                tmp_event_weight4bias = tmp_event_weight4bias*treeIn.cps_Weight_H600;
             if TString(label).Contains("vbfH700"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H700;
                tmp_event_weight4bias = tmp_event_weight4bias*treeIn.cps_Weight_H700;
             if TString(label).Contains("vbfH800"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H800;
                tmp_event_weight4bias = tmp_event_weight4bias*treeIn.cps_Weight_H800;
             if TString(label).Contains("vbfH900"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H900;
                tmp_event_weight4bias = tmp_event_weight4bias*treeIn.cps_Weight_H900;
             if TString(label).Contains("vbfH1000"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H1000;
                tmp_event_weight4bias = tmp_event_weight4bias*treeIn.cps_Weight_H1000;

             if not label=="_data":
                     if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
                         tmp_event_weight      = tmp_event_weight*self.rrv_wtagger_eff_reweight_forT.getVal();
                         tmp_event_weight4bias = tmp_event_weight4bias*self.rrv_wtagger_eff_reweight_forT.getVal();
                     elif TString(label).Contains("_ggH") or TString(label).Contains("_vbfH") or TString(label).Contains("_VV") or TString(label).Contains("_WW_EWK") :
                         tmp_event_weight      = tmp_event_weight*self.rrv_wtagger_eff_reweight_forV.getVal();
                         tmp_event_weight4bias = tmp_event_weight4bias*self.rrv_wtagger_eff_reweight_forV.getVal();

             tmp_event_weight       = tmp_event_weight* getattr(treeIn,"btag_weight"); ## add the btag weight 
             tmp_event_weight4bias  = tmp_event_weight4bias* getattr(treeIn,"btag_weight"); ## add the btag weight 

             ## central values
             rrv_mass_lvj.setVal(mass_lvj);
             rrv_mass_j.setVal(tmp_jet_mass);

             if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max and isFullVBF >= 2:
                 rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                 data_category.setLabel("sideband");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);

             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max and isFullVBF >= 2:
                 rdataset_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                 data_category.setLabel("signal_region");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);
                   
             if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max and isFullVBF >= 2:
                 rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                                                     
        print "########### Nominal Value ###########";
        rrv_scale_to_lumi            = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,(rdataset_sb_lo_mlvj.sumEntries()+rdataset_signal_region_mlvj.sumEntries()+rdataset_sb_hi_mlvj.sumEntries())/(rdataset4bias_sb_lo_mlvj.sumEntries()+rdataset4bias_sb_hi_mlvj.sumEntries()+rdataset4bias_signal_region_mlvj.sumEntries()));
        rrv_scale_to_lumi_sb_lo_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj",(rdataset_sb_lo_mlvj.sumEntries())/(rdataset4bias_sb_lo_mlvj.sumEntries()));
        rrv_scale_to_lumi_sb_hi_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj",(rdataset_sb_hi_mlvj.sumEntries())/(rdataset4bias_sb_hi_mlvj.sumEntries()));
        rrv_scale_to_lumi_signal_region_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj",(rdataset_signal_region_mlvj.sumEntries())/(rdataset4bias_signal_region_mlvj.sumEntries()));
        
        rrv_scale_to_lumi.Print();
        rrv_scale_to_lumi_sb_lo_mlvj.Print();
        rrv_scale_to_lumi_sb_hi_mlvj.Print();
        rrv_scale_to_lumi_signal_region_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi);  
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_sb_lo_mlvj);  
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_sb_hi_mlvj);               
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_signal_region_mlvj);               
                              
        getattr(self.workspace4bias_,"import")(rdataset_sb_lo_mlvj); rdataset_sb_lo_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset_signal_region_mlvj); rdataset_signal_region_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset_sb_hi_mlvj); rdataset_sb_hi_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_lo_mlvj); rdataset4bias_sb_lo_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_signal_region_mlvj); rdataset4bias_signal_region_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_hi_mlvj); rdataset4bias_sb_hi_mlvj.Print();
        getattr(self.workspace4bias_,"import")(combData); combData.Print();
        getattr(self.workspace4bias_,"import")(combData4bias); combData4bias.Print();
                
    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig, dataname = "data", pdfname = "model_mc"):

        print "############### draw the pull plot ########################"
        hpull = mplot_orig.pullHist(dataname,pdfname);        
        x = ROOT.Double(0.); y = ROOT.Double(0) ;
        for ipoint in range(0,hpull.GetN()):
          hpull.GetPoint(ipoint,x,y);
          if(y == 0):
           hpull.SetPoint(ipoint,x,10)
       
        mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins())));
        medianLine = TLine(rrv_x.getMin(),0.,rrv_x.getMax(),0); medianLine.SetLineWidth(2); medianLine.SetLineColor(kRed);
        mplot_pull.addObject(medianLine);
        mplot_pull.addPlotable(hpull,"P");
        mplot_pull.SetTitle("");
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
        

    def getData_PoissonInterval(self,data_obs,mplot):
        rrv_x = self.workspace4bias_.var("rrv_mass_lvj");
        datahist   = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone");
        data_histo = datahist.createHistogram("histo_data",rrv_x) ;
        data_histo.SetName("data");
        data_plot  = RooHist(data_histo);
        data_plot.SetMarkerStyle(20);
        data_plot.SetMarkerSize(1.5);
        
        alpha = 1 - 0.6827;
        for iPoint  in range(data_plot.GetN()):
          N = data_plot.GetY()[iPoint];
          if N==0 : L = 0;
          else : L = (ROOT.Math.gamma_quantile(alpha/2,N,1.));
          U =  ROOT.Math.gamma_quantile_c(alpha/2,N+1,1);
          data_plot.SetPointEYlow(iPoint, N-L);
          data_plot.SetPointEYhigh(iPoint,U-N);
          data_plot.SetPointEXlow(iPoint,0);
          data_plot.SetPointEXhigh(iPoint,0);

        mplot.addPlotable(data_plot,"PE");

    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):
      print "############### draw the banner ########################"


    def banner4Plot(self, iswithpull=0):
      print "############### draw the banner ########################"

      if iswithpull:
       if self.channel=="el":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
       elif self.channel=="mu":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
       elif self.channel=="em":
        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow l #nu "%(self.GetLumi())));
       banner.SetNDC(); banner.SetTextSize(0.04);
      else:
       if self.channel=="el":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
       if self.channel=="mu":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
       if self.channel=="em":
        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow l #nu "%(self.GetLumi())));
       banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
      return banner;

    ### in order to make the legend                                                                                                                                                        
    def legend4Plot(self, plot, left=1, isFill=1, x_offset_low=0., y_offset_low=0., x_offset_high =0., y_offset_high =0., TwoCoulum =1.):
        print "############### draw the legend ########################"
        if left==-1:
            theLeg = TLegend(0.65+x_offset_low, 0.58+y_offset_low, 0.93+x_offset_low, 0.87+y_offset_low, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetLineColor(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.41+x_offset_low, 0.61+y_offset_low, 0.76+x_offset_high, 0.93+y_offset_high, "", "NDC");
            theLeg.SetName("theLegend");
            if TwoCoulum :
                theLeg.SetNColumns(2);

        theLeg.SetFillColor(0);
        theLeg.SetFillStyle(0);
        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextSize(0.040);
        theLeg.SetTextFont(42);

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";

        if   self.channel == "mu": legHeader="(#mu#nu, 1JHP)";
        if   self.channel == "el": legHeader="(e#nu, 1JHP)";
        if   self.channel == "em": legHeader="(l#nu, 1JHP)";

        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="data" : theLeg.AddEntry(theObj, "CMS Data "+legHeader,"PE");  objName_before=objName;
            else: objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";

        for obj in range(int(plot.numItems()) ):
          objName = plot.nameOf(obj);
          if objName == "errorband" : objName = "Uncertainty";
          print objName;
          if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or
objName ==objName_before ):
            theObj = plot.getObject(obj);
            objTitle = objName;
            drawoption= plot.getDrawOptions(objName).Data()
            if drawoption=="P":drawoption="PE"
            if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):  objName_before=objName; continue ;
            elif TString(objName).Contains("Graph") :  objName_before=objName; continue ;
            elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");  objName_before=objName;
            else:  objName_before=objName; continue ;

        entryCnt = 0;
        objName_before = "";
        objName_signal_graviton = "";
        objNameLeg_signal_graviton = "";

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            if objName == "errorband" : objName = "Uncertainty";
            print objName;
            if not ( ( (plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName ==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data()=="STop" : theLeg.AddEntry(theObj, "Single Top","F");
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t}","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ","F");
                    elif TString(objName).Data()=="data" :  objName_before=objName; entryCnt = entryCnt+1; continue ;
                    elif TString(objName).Data()=="WJets" : objName_before=objName; entryCnt = entryCnt+1; continue;
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    elif TString(objName).Contains("Uncertainty"): theLeg.AddEntry(theObj, objTitle,drawoption);
                    elif TString(objName).Contains("Bulk"):
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.6 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.7 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.8 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.9 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.3 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.4 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.5 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1600"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.6 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1700"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.7 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1800"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.8 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1900"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.9 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2000"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2100"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.1 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2200"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.2 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2300"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.3 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2400"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.4 TeV #tilde{k}=0.5 (#times100)";
                       if TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2500"):
                           objName_signal_graviton = theObj ;
                           objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.5 TeV #tilde{k}=0.5 (#times100)";
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;
        if objName_signal_graviton !="" :
           theLeg.AddEntry(objName_signal_graviton, TString(objNameLeg_signal_graviton).Data() ,"L");
        return theLeg;


    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull

        print "############### draw the canvas with pull ########################"
        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.055);
        mplot.GetYaxis().SetTitleSize(0.055);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot_pull.GetXaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetLabelSize(0.14);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);

                                                                          
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);
        # if parameters_list is empty, don't draw pad3
        par_first=parameters_list.createIterator();
        par_first.Reset();
        param_first=par_first.Next()
        doParameterPlot = 0 ;
        if param_first and doParameterPlot != 0:
         pad1=TPad("pad1","pad1",0.,0. ,0.8,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.8,1. );
         pad3=TPad("pad3","pad3",0.8,0.,1,1);
         pad1.Draw();
         pad2.Draw();
         pad3.Draw();
        else:
         pad1=TPad("pad1","pad1",0.,0. ,0.99,0.24);
         pad2=TPad("pad2","pad2",0.,0.24,0.99,1. );
         pad1.Draw();
         pad2.Draw();
                                                                                                                                                                              
        pad2.cd();
        mplot.Draw();
        banner = self.banner4Plot(1);
        banner.Draw();

        pad1.cd();
        mplot_pull.Draw();

        if param_first and doParameterPlot != 0:

            pad3.cd();
            latex=TLatex();
            latex.SetTextSize(0.1);
            par=parameters_list.createIterator();
            par.Reset();
            param=par.Next()
            i=0;
            while param:
                if (not param.isConstant() ) or show_constant_parameter:
                    param.Print();
                    icolor=1;#if a paramenter is constant, color is 2
                    if param.isConstant(): icolor=2
                    latex.DrawLatex(0,0.9-i*0.04,"#color[%s]{%s}"%(icolor,param.GetName()) );
                    latex.DrawLatex(0,0.9-i*0.04-0.02," #color[%s]{%4.3e +/- %2.1e}"%(icolor,param.getVal(),param.getError()) );
                    i=i+1;
                param=par.Next();

        ## create the directory where store the plots
        Directory = TString(in_directory+self.ggH_sample);
        if not Directory.EndsWith("/"):Directory = Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file = TString(Directory.Data()+in_file_name);
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
            mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy,1);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0, frompull=0):

        print "############### draw the canvas without pull ########################"
        cMassFit = TCanvas("cMassFit","cMassFit", 600,600);

        if frompull and logy :
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()/100)
        elif not frompull and logy :
            in_obj.GetYaxis().SetRangeUser(0.00001,in_obj.GetMaximum())
            

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

        banner = self.banner4Plot();
        banner.Draw();
        
        Directory=TString(in_directory+self.ggH_sample);
        if not Directory.EndsWith("/"):Directory=Directory.Append("/");
        if not os.path.isdir(Directory.Data()):
              os.system("mkdir -p "+Directory.Data());

        rlt_file=TString(Directory.Data()+in_file_name);
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
            in_obj.GetYaxis().SetRangeUser(1e-2,in_obj.GetMaximum()*100);
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());
       

    ##### Get Lumi for banner title
    def GetLumi(self):
        if self.channel=="el": return 19.3;
        if self.channel=="mu": return 19.3;
        if self.channel=="em": return 19.3;


    def shapeParametrizationAnalysis(self):

     ## get WJets and fit it in the sb
     if options.isMC == 1:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj                                                                                             
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfExp_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfExp_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfExp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfExp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Erf2Exp",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Erf2Exp",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfPow_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfPow_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfPow2_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfPow2_v1",1,0,1); 
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfPow3_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfPow3_v1",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Erf2Pow",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Erf2Pow",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfChebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfChebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfChebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfChebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfChebychev_v4",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ErfChebychev_v4",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Keys",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Keys",1,0,1);

      else:

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj                                                                                             

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Exp",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Exp",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","ExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Exp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Exp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","2Exp",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","2Exp",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Pow",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Pow",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Pow2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Pow2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Pow3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Pow3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","2Pow",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","2Pow",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Chebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Chebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Chebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Chebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Chebychev_v4",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Chebychev_v4",1,0,1);

       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Keys",1,0,1);
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_signal_region","Keys",1,0,1);
         
     else:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV","_sb_lo","ErfExp_v1",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                               
       self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sb_lo","ErfExp_v1",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                           
       self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sb_lo","ErfExp_v1",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                   
       self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK","_sb_lo","ErfExp_v1",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj                                                                          
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","ErfExp_v1",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_data,"_data"); ## global fit of data in the sidand fixing non dominant bkg
 
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfExp_v1");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfExpTail");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfExp_v3");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Erf2Exp");

       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfPow_v1");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfPow2_v1");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfPow3_v1");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Erf2Pow");

       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfChebychev_v2");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfChebychev_v3");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfChebychev_v4");

       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ErfPowExp_v1");

      else: 

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV","_sb_lo","Exp",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                                                                 
       self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sb_lo","Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
       self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sb_lo","Exp",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                                               
       self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK","_sb_lo","Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj                                                                                             
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo","Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_data,"_data"); ## global fit of data in the sidand fixing non dominant bkg
 
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Exp");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","ExpTail");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Exp_v3");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","2Exp");

       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Pow");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Pow2");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Pow3");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","2Pow");

       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Chebychev_v2");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Chebychev_v3");
       self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo","Chebychev_v4");

    def biasAnalysis(self):
   
     print"######################## begin the bias analysis ###########################";  
    
     ## get the signal and fit it     
     self.get_mj_and_mlvj_dataset(self.file_ggH ,"_%s"%(self.ggH_sample), "jet_mass_pr")# to get the shape of m_lvj
     self.fit_mlvj_model_single_MC(self.file_ggH,"_%s"%(self.ggH_sample),"_signal_region","CB_v1", 0, 0, 1);
     
     self.get_mj_and_mlvj_dataset(self.file_vbfH,"_%s"%(self.vbfhiggs_sample), "jet_mass_pr")# to get the shape of m_lvj
     self.fit_mlvj_model_single_MC(self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"_signal_region","CB_v1", 0, 0, 1);

     if options.isMC == 0:
      ## get diboson and fit it
      self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
      self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV","_sb_lo",options.fgen,0,0,1);
 
      ## get SingleTop and fit it
      self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                                                                 
      self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop","_sb_lo",options.fgen,0,0,1);

      ## get TTbar and fit it
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar","_sb_lo",options.fgen,0,0,1);
            
      ## get WW EWK and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                                               
      self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK","_sb_lo",options.fgen,0,0,1);

     ## get WJets and fit it in the sb
     self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
     self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0","_sb_lo",options.fgen,0,0,1);

     if options.isMC == 0:
 
      ## get data in sb and fit it                                                                                                                                              
      self.get_mj_and_mlvj_dataset(self.file_data,"_data", "jet_mass_pr"); ## global fit of data in the sidand fixing non dominant bkg             
      self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo",options.fgen,1);


     ## fix signal and bkg models that are going to be used in the generation
     self.fix_Model("_%s"%self.ggH_sample,"_signal_region","_mlvj","CB_v1");
     self.fix_Model("_%s"%self.vbfhiggs_sample,"_signal_region","_mlvj","CB_v1");
     self.fix_Model("_WJets0","_sb_lo","_mlvj",options.fgen);

     if options.isMC == 0:

      self.fix_Model("_TTbar" ,"_sb_lo","_mlvj",options.fgen);
      self.fix_Model("_STop"  ,"_sb_lo","_mlvj",options.fgen);
      self.fix_Model("_VV"    ,"_sb_lo","_mlvj",options.fgen) ;
      self.fix_Model("_WW_EWK","_sb_lo","_mlvj",options.fgen) ;
          
     ### clone the signal shape --> parameter already fixed
     fitted_signal_ggH   = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.ggH_sample,"signal_region",self.channel,"mlvj"));
     fitted_signal_vbfH  = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.vbfhiggs_sample,"signal_region",self.channel,"mlvj"));

     ### make the pdf for ggH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_ggH = [];
     model_signal_ggH = self.make_Pdf("_%s_signal_region_fit"%(self.ggH_sample),"CB_v1","_mlvj",constrainslist_signal_ggH,1);
     model_signal_ggH.Print();     
     
     self.clone_Model(model_signal_ggH,"_%s"%self.ggH_sample,"_signal_region","_mlvj","CB_v1");
     self.fix_Model("_%s"%self.ggH_sample,"_signal_region","_mlvj","","_fit",1);
      
     rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
     rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal());
     rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getError());
     rrv_number_signal_signal_fit_ggH.Print();
 
     ### make the pdf for vbfH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_vbfH = [];
     model_signal_vbfH = self.make_Pdf("_%s_signal_region_fit"%(self.vbfhiggs_sample),"CB_v1","_mlvj",constrainslist_signal_vbfH,1);
     model_signal_vbfH.Print();     
     
     self.clone_Model(model_signal_vbfH,"_%s"%self.vbfhiggs_sample,"_signal_region","_mlvj","CB_v1");
     self.fix_Model("_%s"%self.vbfhiggs_sample,"_signal_region","_mlvj","","_fit",1);
      
     rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
     rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal());
     rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getError());
     rrv_number_signal_signal_fit_vbfH.Print();

     ## total pdf: ggH + vbfH 
     rrv_fraction_ggH_vbf    = RooRealVar("rrv_fraction_ggH_vbf","rrv_fraction_ggH_vbf",rrv_number_signal_signal_fit_vbfH.getVal()/rrv_number_signal_signal_fit_ggH.getVal());
     rrv_fraction_ggH_vbf.setConstant(kTRUE);
     rrv_fraction_ggH_vbf.Print();

     model_modified_signal   = RooAddPdf("model_modified_signal","model_modified_signal",RooArgList(model_signal_vbfH,model_signal_ggH),RooArgList(rrv_fraction_ggH_vbf),1);
     model_modified_signal.Print();
     getattr(self.workspace4bias_,"import")(model_modified_signal);
  

     rrv_number_signal_region_fit_ggH = RooRealVar("rrv_number_signal_region_fit_H","rrv_number_signal_region_fit_H",0,-1e8,1e8);
     rrv_number_signal_region_fit_ggH.Print();

     model_total_signal = RooExtendPdf("model_higgs_signal_region_fit_%s_mlvj"%(self.channel),
                                       "model_higgs_signal_region_fit_%s_mlvj"%(self.channel),model_modified_signal,rrv_number_signal_region_fit_ggH);

     model_total_signal.Print();
     getattr(self.workspace4bias_,"import")(model_total_signal);
          
     ############### Make the MC analysis --> make the Entended pdf for the bkg
     if options.isMC == 1 :

      constrainslist_bkg_wjet = [];
      model_bkg_wjet    = self.make_Model("_WJets0_sb_lo_fit",options.fres,"_mlvj",constrainslist_bkg_wjet,1);
      model_bkg_wjet.Print();

      self.clone_Model(model_bkg_wjet,"_WJets0","_sb_lo","_mlvj",options.fgen);
      self.workspace4bias_.var("rrv_number_WJets0_sb_lo_fit_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number_WJets0_sb_lo"+options.fgen+"_"+self.channel+"_mlvj").getVal());

      ##### Total model for MC
      model_Total_mc    = RooAddPdf("model_Total_background_mc","model_Total_background_mc",RooArgList(model_total_signal,model_bkg_wjet));

      model_Total_mc.Print();
      getattr(self.workspace4bias_,"import")(model_Total_mc);

      ##### generate models  --> the one fixed and already fitted
      generation_model_wjet = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%("WJets0","sb_lo"+options.fgen,self.channel,"mlvj"));
      generation_model_wjet.Print();

      self.workspace4bias_.Print();             

      ## variable and RooMC study
      numevents_mc   = self.workspace4bias_.data("rdataset"+"_WJets0"+"_sb_lo"+"_"+self.channel+"_mlvj").sumEntries();
      print"########  numevents mc ",numevents_mc;
       
      
      mc_wjet = RooMCStudy(generation_model_wjet,
                           RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                           RooFit.FitModel(model_Total_mc),
                           RooFit.FitOptions(RooFit.Save(kTRUE),RooFit.SumW2Error(kTRUE),RooFit.Minimizer("Minuit2"),RooFit.Extended(kTRUE)),
                           RooFit.Extended(kTRUE),
#                           RooFit.Silence()
      );

      mc_wjet.Print();     
      ## create module for chi2 evaluation
      chi2_mc = RooChi2MCSModule();
      mc_wjet.addModule(chi2_mc);
     
      ## generate and fit storing the generated distribution for each toy
      mc_wjet.generateAndFit(options.nexp,int(numevents_mc),1);

      generatedData_wjet           = []; ## distribution of generated toy according to bkg only hypothesis
      fittedPdf_wjet               = []; ## fitted pdf signal + bkg
      fitResults_wjet              = [];
      parameterHisto_wjet          = []; ## histo of the parameters of the fitted pdf
      parameterHisto_fraction_wjet = []; ## histo of the parameters of the fitted pdf
      parameterHistoError_wjet     = []; ## errror on the parameters
      parameterHistoPull_wjet      = []; ## pull wrt the generated one
      chi2distribution_wjet        = [];
      nLLdistribution_wjet         = [];
    
      for iToy in range(options.nexp): ## loop on the toy
 
         if not mc_wjet.genData(iToy) : continue ;
         generatedData_wjet.append(mc_wjet.genData(iToy)); ## take the generated dataset ad store it
         fitResults_wjet.append(mc_wjet.fitResult(iToy));  ## take the fit result and store it        
         fittedPdf_wjet.append(model_Total_mc.Clone("model_Total_mc_toy_%d"%iToy)); ## clone the pdf and take the list to re-build the pdf shape for plotting reason

         if not mc_wjet.fitResult(iToy) : continue ;
         if fitResults_wjet[len(fitResults_wjet)-1].status() != 0 : continue ; ## skip not status zero fit
          
         parset = mc_wjet.fitParams(iToy); ## get the parameters of the fit only those non constant
         if not parset : continue ;
         parlist = RooArgList(parset);
         
         param       = fittedPdf_wjet[len(fittedPdf_wjet)-1].getParameters(generatedData_wjet[len(generatedData_wjet)-1]); ## parameter of the new pdf .. all the parameters
         if not param : continue ;
         parameters = RooArgList(param); 

         param_generated  = generation_model_wjet.getParameters(self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%("_WJets0","_sb_lo",self.channel,"_mlvj"))); ## parameters used in the generation
         ### parameters used in the generation
         
         if not param_generated : continue ;
         parameters_generated = RooArgList(param_generated);

         iparNotConstant = 0;
         iGenerated = 0 ;
         iPull      = 0 ;    

         for ipar in range(parameters.getSize()): ## to clone the pdf value for plotting reasons
              if (parameters.at(ipar).GetName() == parlist.at(ipar).GetName()):
                  parameters.at(ipar).setVal(parlist.at(ipar).getVal());
                  parameters.at(ipar).setError(parlist.at(ipar).getError()) ## copy parameters in order to have the fitted shape
                  
              if TString(parameters.at(ipar).GetName()).Contains("signal_region") and not TString(parameters.at(ipar).GetName()).Contains("number"): continue ; ## skip signal par 
              if TString(parameters.at(ipar).GetName()).Contains("rrv_fraction_ggH_vbf"): continue ;
             
              if iToy == 0 or len(parameterHisto_wjet)==0 : ## create parameter histo and pull 

               if not TString(parlist.at(ipar).GetName()).Contains("signal_region"):
                 if not TString(parlist.at(ipar).GetName()).Contains("number") :  
                  parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",75,-math.fabs(parlist.at(ipar).getVal()*5),math.fabs(parlist.at(ipar).getVal()*5)));
                  parameterHisto_fraction_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_fraction","",50,-50,50));
                  parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",40,0.,math.fabs(parlist.at(ipar).getError()*3)));
                  if options.fgen == options.fres :
                   parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",30,-5,5));
                 else : 
                  parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",75,math.fabs(parlist.at(ipar).getVal())/3,math.fabs(parlist.at(ipar).getVal())*3));
                  parameterHisto_fraction_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_fraction","",50,-50,50));
                  parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",40,0.,math.fabs(parlist.at(ipar).getError()*3)));
                  parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",30,-5,5)); 
               else:
                  parameterHisto_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName(),"",100,-50,50));
                  parameterHisto_fraction_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_fraction","",50,-50,50));
                  parameterHistoError_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_error","",50,0.,math.fabs(parlist.at(ipar).getError()*3)));
                  parameterHistoPull_wjet.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_pull","",30,-5,5)); 
                 
              if not TString(parlist.at(ipar).GetName()).Contains("signal_region"): ## fill pulls and parameters histo
               parameterHisto_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal());
               parameterHistoError_wjet[iparNotConstant].Fill(parlist.at(ipar).getError());

               if not TString(parlist.at(ipar).GetName()).Contains("number") and options.fgen == options.fres:
                print "parameter name ",parlist.at(ipar).GetName()," value ",parlist.at(ipar).getVal()," error ",parlist.at(ipar).getError()
                parameterHisto_fraction_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal()/parameters_generated.at(iGenerated).getVal());
                parameterHistoPull_wjet[iPull].Fill((parlist.at(ipar).getVal()-parameters_generated.at(iGenerated).getVal())/parlist.at(ipar).getError());
                iPull = iPull +1;
               elif TString(parlist.at(ipar).GetName()).Contains("number"):
                print "parameter name ",parlist.at(ipar).GetName()," value ",parlist.at(ipar).getVal()," error ",parlist.at(ipar).getError()
                parameterHistoPull_wjet[iPull].Fill((parlist.at(ipar).getVal()-parlist.find("ngen").getVal())/parlist.at(ipar).getError());
                parameterHisto_fraction_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal()/parlist.find("ngen").getVal());
                iPull = iPull +1;                
               iGenerated = iGenerated +1 ;
              else:
                 ### fill with the fraction of the fitted signal events (positive or negative) divided by what is given by mc -> signal strentght in a only bkg generation
                 print "parameter name ",parlist.at(ipar).GetName()," value ",parlist.at(ipar).getVal()," error ",parlist.at(ipar).getError()
                 parameterHisto_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal());
                 parameterHistoError_wjet[iparNotConstant].Fill(parlist.at(ipar).getError());  
                 parameterHistoPull_wjet[iPull].Fill(parlist.at(ipar).getVal()/parlist.at(ipar).getError());
                 parameterHisto_fraction_wjet[iparNotConstant].Fill(parlist.at(ipar).getVal()/(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal()+self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal()));                  

                 iPull = iPull +1;
                    
              iparNotConstant = iparNotConstant+1;

         ## fill chi2, NNLL
         if len(chi2distribution_wjet)==0 and parlist.find("chi2red"):
          chi2distribution_wjet.append(ROOT.TH1F("chi2distribution_wjet","",30,0.,3));
         if len(nLLdistribution_wjet)==0 and parlist.find("NLL"):
          nLLdistribution_wjet.append(ROOT.TH1F("nLLdistribution_wjet","",40,math.fabs(parlist.find("NLL").getVal())*0.5,math.fabs(parlist.find("NLL").getVal())*2));

         if parlist.find("chi2red") :
          chi2distribution_wjet[0].Fill(parlist.find("chi2red").getVal());
         if parlist.find("NLL") : 
          nLLdistribution_wjet[0].Fill(parlist.find("NLL").getVal());

      ### Plot in Canvas + Gaussian fit of each histogram 
      canvas_generatedToys_wjet        = []; ## canvas for show generated distribution + fits
      canvas_parameters_wjet           = []; ## canvas parameters
      canvas_parameters_fraction_wjet  = []; ## canvas parameters
      canvas_parameters_err_wjet       = [];
      canvas_parameters_pull_wjet      = [];
      canvas_chi2_wjet                 = [];
      canvas_chi2_wjet_frame           = [];
      canvas_nLL_wjet                  = [];

      chi2distribution_wjet_frame   = [];

      ### print the canvas of the single jobs
      for iObj in range(len(generatedData_wjet)):
       storethisPlot = 0;
       if options.nexp <= 10 :
           storethisPlot = 1;
       elif options.nexp <= 100 and iObj%3 == 0:
           storethisPlot = 1;
       elif options.nexp > 100 and iObj%10 == 0:
           storethisPlot = 1;
                   
       parameters = fittedPdf_wjet[iObj].getParameters(generatedData_wjet[iObj]).selectByAttrib("Constant",kFALSE) ; ## take only nn constant parameters
        
       if parameters :
        wjet_binned   = generatedData_wjet[iObj].binnedClone(); ## binned cloning of the unbinned dataset
        ChiSquare = fittedPdf_wjet[iObj].createChi2(wjet_binned,RooFit.Extended(kTRUE),RooFit.SumW2Error(kTRUE)); ## create the chi2
        if len(chi2distribution_wjet_frame) ==0:
         chi2distribution_wjet_frame.append(ROOT.TH1F("chi2distribution_wjet_frame",""30,0,30));
         
        chi2distribution_wjet_frame[0].Fill(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize()));

       if options.storeplot and  storethisPlot == 1:

         mplot = self.workspace4bias_.var("rrv_mass_lvj").frame(RooFit.Title("frame_generatedToys_%d"%iObj), RooFit.Bins(self.workspace4bias_.var("rrv_mass_lvj").getBins()));      

         generatedData_wjet[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
         draw_error_band_extendPdf(generatedData_wjet[iObj], fittedPdf_wjet[iObj], fitResults_wjet[iObj],mplot,2,"L");
         generatedData_wjet[iObj].plotOn(mplot,RooFit.MarkerSize(1.5),RooFit.DataError(RooAbsData.SumW2),RooFit.XErrorSize(0),RooFit.Name(generatedData_wjet[iObj].GetName()+"_curve"));
         fittedPdf_wjet[iObj].plotOn(mplot,RooFit.Name(fittedPdf_wjet[iObj].GetName()+"_curve"));
                     
         canvas_generatedToys_wjet.append(ROOT.TCanvas("canvas_generatedToys_wjet_%d"%iObj,""));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].cd();
         ROOT.SetOwnership(canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1], False);
         pad1 = ROOT.TPad("pad1_%d"%iObj,"pad1_%d"%iObj,0.,0.24,0.99,1. );
         ROOT.SetOwnership(pad1, False);
         pad1.Draw();
         pad2 = ROOT.TPad("pad2_%d"%iObj,"pad2_%d"%iObj,0.,0.,0.99,0.24 );
         ROOT.SetOwnership(pad2, False);
         pad2.Draw();
         pad1.cd();
         mplot.GetXaxis().SetTitleOffset(1.1);
         mplot.GetYaxis().SetTitleOffset(1.3);
         mplot.GetXaxis().SetTitleSize(0.05);
         mplot.GetYaxis().SetTitleSize(0.05);
         mplot.GetXaxis().SetLabelSize(0.045);
         mplot.GetYaxis().SetLabelSize(0.045);
         mplot.Draw();
         banner = self.banner4Plot() ;
         banner.Draw();
                                
         pad2.cd();
         mplot_pull = self.get_pull(self.workspace4bias_.var("rrv_mass_lvj"), mplot,generatedData_wjet[iObj].GetName()+"_curve",fittedPdf_wjet[iObj].GetName()+"_curve");
         mplot_pull.Draw();
         mplot_pull.GetXaxis().SetLabelSize(0.15);
         mplot_pull.GetYaxis().SetLabelSize(0.15);
         mplot_pull.GetYaxis().SetTitleSize(0.15);
         mplot_pull.GetYaxis().SetNdivisions(205);

         if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
          os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
          os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));

         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));

         pad1.SetLogy();
         pad1.Update();
         mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].Update()

         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s_log.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
         canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s_log.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_wjet[len(canvas_generatedToys_wjet)-1].GetName()));
        

      ### print plots of the parameters
      for ipar in range(len(parameterHisto_wjet)):
       canvas_parameters_wjet.append(ROOT.TCanvas("canvas_parameters_wjet_%s"%parameterHisto_wjet[ipar].GetName(),""));
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_wjet[len(canvas_parameters_wjet)-1], False);
       parameterHisto_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHisto_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHisto_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHisto_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHisto_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHisto_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHisto_wjet[ipar].GetXaxis().SetTitle(parameterHisto_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_par_%s"%(parameterHisto_wjet[ipar].GetName()),"gaus",parameterHisto_wjet[ipar].GetXaxis().GetXmin(),parameterHisto_wjet[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHisto_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHisto_wjet[ipar].Draw("E");

       self.outputFilePull.cd();
       Gaussian.Write();
       parameterHisto_wjet[ipar].Write();
       
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
         os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
         
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_wjet[len(canvas_parameters_wjet)-1].GetName()));
       canvas_parameters_wjet[len(canvas_parameters_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_wjet[len(canvas_parameters_wjet)-1].GetName()));

      ### print plots of the parameters ratio wrt to the fitted one or the expected yield
      for ipar in range(len(parameterHisto_fraction_wjet)):
       canvas_parameters_fraction_wjet.append(ROOT.TCanvas("canvas_parameters_wjet_%s"%parameterHisto_fraction_wjet[ipar].GetName(),""));
       canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1], False);
       parameterHisto_fraction_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHisto_fraction_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHisto_fraction_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHisto_fraction_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHisto_fraction_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHisto_fraction_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHisto_fraction_wjet[ipar].GetXaxis().SetTitle(parameterHisto_fraction_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_par_fraction%s"%(parameterHisto_fraction_wjet[ipar].GetName()),"gaus",parameterHisto_fraction_wjet[ipar].GetXaxis().GetXmin(),parameterHisto_fraction_wjet[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHisto_fraction_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHisto_fraction_wjet[ipar].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
         os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
         
       canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1].GetName()));
       canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_fraction_wjet[len(canvas_parameters_fraction_wjet)-1].GetName()));


      ### print the parameters error distribution 
      for ipar in range(len(parameterHistoError_wjet)):
       canvas_parameters_err_wjet.append(ROOT.TCanvas("canvas_parameterHistoError_wjet_%s"%parameterHistoError_wjet[ipar].GetName(),""));
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1], False);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoError_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoError_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoError_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoError_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoError_wjet[ipar].GetXaxis().SetTitle(parameterHistoError_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_error_%s"%(parameterHistoError_wjet[ipar].GetName()),"gaus",parameterHistoError_wjet[ipar].GetXaxis().GetXmin(),parameterHistoError_wjet[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoError_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHistoError_wjet[ipar].Draw("E");
        
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
         
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].GetName()));
       canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_wjet[len(canvas_parameters_err_wjet)-1].GetName()));
 
      ## print pulls of each parameter
      for ipar in range(len(parameterHistoPull_wjet)):
       canvas_parameters_pull_wjet.append(ROOT.TCanvas("canvas_parameters_pull_wjet_%s"%parameterHistoPull_wjet[ipar].GetName(),""));
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].cd();
       ROOT.SetOwnership(canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1], False);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoPull_wjet[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoPull_wjet[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoPull_wjet[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoPull_wjet[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoPull_wjet[ipar].GetXaxis().SetTitle(parameterHistoPull_wjet[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_pull_%s"%(parameterHistoPull_wjet[ipar].GetName()),"gaus",-5,5);
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoPull_wjet[ipar].Fit(Gaussian,"MSQ");
       parameterHistoPull_wjet[ipar].Draw("E");

       self.outputFilePull.cd();
       Gaussian.Write();
       parameterHistoPull_wjet[ipar].Write();


       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
         
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].GetName()));
       canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_wjet[len(canvas_parameters_pull_wjet)-1].GetName()));
  
      ## print chi2
      if len(chi2distribution_wjet) !=0:
       canvas_chi2_wjet.append(ROOT.TCanvas("canvas_chi2_wjet_%s"%chi2distribution_wjet[0].GetName(),""));
       canvas_chi2_wjet[0].cd();
       ROOT.SetOwnership(canvas_chi2_wjet[0], False);
       chi2distribution_wjet[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_wjet[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_wjet[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_wjet[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_wjet[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_wjet[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_wjet[0].GetXaxis().SetTitle(chi2distribution_wjet[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_wjet[0].GetXaxis().GetXmin(),chi2distribution_wjet[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_wjet[0].Fit(Gaussian,"MSQ");
       chi2distribution_wjet[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
        
       canvas_chi2_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet[0].GetName()));
       canvas_chi2_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet[0].GetName()));

      ## print chi2 from frame
      if len(chi2distribution_wjet_frame) !=0:
       canvas_chi2_wjet_frame.append(ROOT.TCanvas("canvas_chi2_wjet_frame_%s"%chi2distribution_wjet_frame[0].GetName(),""));
       canvas_chi2_wjet_frame[0].cd();
       ROOT.SetOwnership(canvas_chi2_wjet_frame[0], False);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_wjet_frame[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_wjet_frame[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_wjet_frame[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_wjet_frame[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_wjet_frame[0].GetXaxis().SetTitle(chi2distribution_wjet_frame[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_wjet_frame[0].GetXaxis().GetXmin(),chi2distribution_wjet_frame[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_wjet_frame[0].Fit(Gaussian,"MSQ");
       chi2distribution_wjet_frame[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
        
       canvas_chi2_wjet_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet_frame[0].GetName()));
       canvas_chi2_wjet_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_wjet_frame[0].GetName()));

      ## print -log(l)
      if len(nLLdistribution_wjet)!=0:
       canvas_nLL_wjet.append(ROOT.TCanvas("canvas_nLL_wjet_%s"%nLLdistribution_wjet[0].GetName(),""));
       canvas_nLL_wjet[0].cd();
       ROOT.SetOwnership(canvas_nLL_wjet[0], False);
       nLLdistribution_wjet[0].GetXaxis().SetTitleOffset(1.1);
       nLLdistribution_wjet[0].GetYaxis().SetTitleOffset(1.3);
       nLLdistribution_wjet[0].GetXaxis().SetTitleSize(0.04);
       nLLdistribution_wjet[0].GetYaxis().SetTitleSize(0.04);
       nLLdistribution_wjet[0].GetXaxis().SetLabelSize(0.035);
       nLLdistribution_wjet[0].GetYaxis().SetLabelSize(0.035); 
       nLLdistribution_wjet[0].GetXaxis().SetTitle(nLLdistribution_wjet[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",nLLdistribution_wjet[0].GetXaxis().GetXmin(),nLLdistribution_wjet[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       nLLdistribution_wjet[0].Fit(Gaussian,"MSQ");
       nLLdistribution_wjet[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toy_mc")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys_mc"));
        
       canvas_nLL_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_wjet[0].GetName()));
       canvas_nLL_wjet[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys_mc/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_wjet[0].GetName()));

     else: 
         
      ############### Make the Data analysis --> make the Entended pdf for the bkg
      constrainslist_bkg_data = [];

      model_VV_backgrounds     = self.get_VV_mlvj_Model("_sb_lo",options.fgen);
      model_STop_backgrounds   = self.get_STop_mlvj_Model("_sb_lo",options.fgen);
      model_TTbar_backgrounds  = self.get_TTbar_mlvj_Model("_sb_lo",options.fgen);
      model_WW_EWK_backgrounds = self.get_WW_EWK_mlvj_Model("_sb_lo",options.fgen);

      model_bkg_data    = self.make_Model("_data_sb_lo_fit",options.fres,"_mlvj",constrainslist_bkg_data,1); ## basic model used for fit in the toys
      model_bkg_data.Print();

      self.clone_Model(model_bkg_data,"_WJets0","_sb_lo","_mlvj",options.fgen+"_from_fitting");
      self.workspace4bias_.var("rrv_number_data_sb_lo_fit_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number_WJets0_sb_lo"+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal());

      ## Add the other bkg component fixed to the total model --> in the extended way
      model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
                                                                                                       
      model_Total_data.Print();
      getattr(self.workspace4bias_,"import")(model_Total_data);

      ##### generate models  --> the one fixed and already fitted
      generation_model_data = self.workspace4bias_.pdf("model_%s_%s_%s_%s_%s"%("data","WJets0","sb_lo"+options.fgen,self.channel,"mlvj"));
      generation_model_data.Print();

      self.workspace4bias_.Print();             

      ## variable and RooMC study
      numevents_data   = self.workspace4bias_.data("rdataset"+"_data"+"_sb_lo"+"_"+self.channel+"_mlvj").sumEntries();
      print"########  numevents data ",numevents_data;

      data_wjet = RooMCStudy(generation_model_data,
                             RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                             RooFit.FitModel(model_Total_data),
                             RooFit.FitOptions(RooFit.Save(kTRUE),RooFit.Minimizer("Minuit2"),RooFit.Extended(kTRUE)),
                             RooFit.Extended(kTRUE),
#                             RooFit.Silence()
                            );

      data_wjet.Print();
     
      ## create module for chi2 evaluation
      chi2_data = RooChi2MCSModule();
      data_wjet.addModule(chi2_data)

      ## generate and fit storing the generated distribution for each toy
      data_wjet.generateAndFit(options.nexp,int(numevents_data),1);

      generatedData_data         = []; ## distribution of generated toy according to bkg only hypothesis -> RooAbsData
      fittedPdf_data             = []; ## fitted pdf signal + bkg
      fitResults_data            = [];
      parameterHisto_data           = []; ## histo of the parameters of the fitted pdf
      parameterHisto_fraction_data  = []; ## histo of the parameters of the fitted pdf
      parameterHistoError_data      = []; ## errror on the parameters
      parameterHistoPull_data       = []; ## pull wrt the generated one
      chi2distribution_data         = [];
      nLLdistribution_data          = [];
            
      for iToy in range(options.nexp): ## loop on the toy

         generatedData_data.append(data_wjet.genData(iToy)); ## take the generated dataset ad store them
         fitResults_data.append(data_wjet.fitResult(iToy));       
         fittedPdf_data.append(model_Total_data.Clone("model_Total_data_toy_%d"%iToy)); ## clone the pdf and take the list to re-buil the pdf shape for plotting reason
         if not generatedData_data or not fitResults_data : continue ;
         if fitResults_data[len(fitResults_data)-1].status() != 0 : continue ;
         
         parset  = data_wjet.fitParams(iToy); ## get the parameters of the fit
         if not parset : continue ;
         parlist = RooArgList(parset);

         param   = fittedPdf_data[len(fittedPdf_data)-1].getParameters(generatedData_data[len(generatedData_data)-1]); ## parameter of the new pdf
         if not param : continue ;
         parameters = RooArgList(param); 

         param_generated  = generation_model_data.getParameters(self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%("_Wdata","_sb_lo",self.channel,"_mlvj"))); ### parameters used in the generation
         if not param_generated : continue ;
         parameters_generated = RooArgList(param_generated);

         iparNotConstant = 0;
         iGenerated = 0 ;
         iPull = 0;
             
         for ipar in range(parameters.getSize()): ## to clone the pdf value for plotting reasons

              if (parameters.at(ipar).GetName() == parlist.at(ipar).GetName()):
                  parameters.at(ipar).setVal(parlist.at(ipar).getVal());
                  parameters.at(ipar).setError(parlist.at(ipar).getError()) ## copy parameters in order to have the fitted shape
 
             
              if TString(parameters.at(ipar).GetName()).Contains("signal_region") and not TString(parameters.at(ipar).GetName()).Contains("number"): continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_VV_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_STop_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_TTbar_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("_WW_EWK_") : continue ;
              if TString(parameters.at(ipar).GetName()).Contains("rrv_fraction_ggH_vbf"): continue ;

              if iToy == 0 or len(parameterHisto_data)==0 :

               if not TString(parlist.at(ipar).GetName()).Contains("signal_region"):
                 if not TString(parlist.at(ipar).GetName()).Contains("number") :    
                  parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",50,-math.fabs(parlist.at(ipar).getVal()*2),math.fabs(parlist.at(ipar).getVal()*2)));
                  parameterHisto_fraction_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_fraction","",50,-50,50));
                  parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",45,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  if options.fgen == options.fres :
                   parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",40,-5,5));
                 else : 
                  parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",50,0.,math.fabs(parlist.at(ipar).getVal())*2));
                  parameterHisto_fraction_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_fraction","",50,-50,50));
                  parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",100,0.,math.fabs(parlist.at(ipar).getError()*2)));
                  parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",40,-5,5)); 
               else:
                parameterHisto_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data","",50,-50,50));
                parameterHisto_fraction_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_fraction","",500,-50,50));
                parameterHistoError_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_error","",45,0.,math.fabs(parlist.at(ipar).getError()*2)));
                parameterHistoPull_data.append(ROOT.TH1F(parlist.at(ipar).GetName()+"_data_pull","",40,-5,5)); 
                 
              if not TString(parlist.at(ipar).GetName()).Contains("signal_region"): ## fill pulls and parameters histo
               parameterHistoError_data[iparNotConstant].Fill(parlist.at(ipar).getError());
               
               if not TString(parlist.at(ipar).GetName()).Contains("number") :
                parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal());
                parameterHisto_fraction_data[iparNotConstant].Fill(parlist.at(ipar).getVal()/parameters_generated.at(iGenerated).getVal());
                if  options.fgen == options.fres:                                            
                 parameterHistoPull_data[iPull].Fill((parlist.at(ipar).getVal()-parameters_generated.at(iGenerated).getVal())/parlist.at(ipar).getError());
                 iPull = iPull +1;
               elif TString(parlist.at(ipar).GetName()).Contains("number"):
                 parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal()+self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal());
                 parameterHisto_fraction_data[iparNotConstant].Fill((parlist.at(ipar).getVal()+self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal())/parlist.find("ngen").getVal());
                
                 parameterHistoPull_data[iPull].Fill((parlist.at(ipar).getVal()+self.workspace4bias_.var("rrv_number_VV_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_TTbar_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_STop_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()+self.workspace4bias_.var("rrv_number_WW_EWK_sb_lo%s_%s_mlvj"%(options.fgen,self.channel)).getVal()-parlist.find("ngen").getVal())/parlist.at(ipar).getError());

                 iPull = iPull +1;                
               iGenerated = iGenerated +1 ;
              else:
                 ### fill with the fraction of the fitted signal events (positive or negative) divided by what is given by mc -> signal strentght in a only bkg generation  
                 parameterHisto_data[iparNotConstant].Fill(parlist.at(ipar).getVal());
                 parameterHistoError_data[iparNotConstant].Fill(parlist.at(ipar).getError());  
                 parameterHistoPull_data[iPull].Fill(parlist.at(ipar).getVal()/parlist.at(ipar).getError());
                 parameterHisto_fraction_data[iparNotConstant].Fill(parlist.at(ipar).getVal()/(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal()+self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+"_signal_regionCB_v1_"+self.channel+"_mlvj").getVal()));             
 
                 iPull = iPull +1;
                    
              iparNotConstant = iparNotConstant+1;

         ## fill chi2, NNLL
         if len(nLLdistribution_data)==0 and parlist.find("NLL"):
            nLLdistribution_data.append(ROOT.TH1F("nLLdistribution_data,"",40,math.fabs(parlist.find("NLL").getVal())*0.5,math.fabs(parlist.find("NLL").getVal())*2));
             
         if len(chi2distribution_data)==0 and parlist.find("chi2red"):
            chi2distribution_data.append(ROOT.TH1F("chi2distribution_data","",30,0.,3));

         if parlist.find("chi2red") :
          chi2distribution_data[0].Fill(parlist.find("chi2red").getVal());
         if parlist.find("NLL"): 
          nLLdistribution_data[0].Fill(parlist.find("NLL").getVal());

      ### Plot in Canvas + Gaussian fit of each histogram 
      canvas_generatedToys_data        = []; ## canvas for show generated distribution + fits
      canvas_parameters_data           = []; ## canvas parameters
      canvas_parameters_fraction_data  = []; ## canvas parameters
      canvas_parameters_err_data       = [];
      canvas_parameters_pull_data      = [];
      canvas_chi2_data                 = [];
      canvas_nLL_data                  = [];
      canvas_chi2_data_frame           = [];     

      chi2distribution_data_frame      = [];
      ### print the canvas of the single jobs         
      for iObj in range(len(generatedData_data)):
       storethisPlot = 0;

       if options.nexp <= 10 :
           storethisPlot = 1;
       elif options.nexp <= 100 and iObj%3 == 0:
           storethisPlot = 1;
       elif options.nexp > 100 and iObj%10 == 0:
          storethisPlot = 1;

       parameters = fittedPdf_data[iObj].getParameters(generatedData_data[iObj]).selectByAttrib("Constant",kFALSE) ;
       if parameters :
        data_binned   = generatedData_data[iObj].binnedClone();
        ChiSquare = fittedPdf_data[iObj].createChi2(data_binned,RooFit.Extended(kTRUE),RooFit.SumW2Error(kFALSE));
        if len (chi2distribution_data_frame) ==0:
         chi2distribution_data_frame.append(ROOT.TH1F("chi2distribution_data_frame","",30,0,3));
         
        chi2distribution_data_frame[0].Fill(ChiSquare.getVal()/(self.workspace4bias_.var("rrv_mass_lvj").getBins()-parameters.getSize()));
        
       if options.storeplot and  storethisPlot == 1:

        mplot = self.workspace4bias_.var("rrv_mass_lvj").frame(RooFit.Title("frame_generatedToys_data_%d"%iObj), RooFit.Bins(self.workspace4bias_.var("rrv_mass_lvj").getBins()));      
        generatedData_data[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.Invisible(), RooFit.XErrorSize(0));
        draw_error_band_extendPdf(generatedData_data[iObj], fittedPdf_data[iObj], fitResults_data[iObj],mplot,2,"L");
        fittedPdf_data[iObj].plotOn(mplot, RooFit.Name(fittedPdf_data[iObj].GetName()+"_curve"));
        if options.pseudodata == 1:
         generatedData_data[iObj].plotOn(mplot,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0),RooFit.Name(generatedData_data[iObj].GetName()+"_curve"));
        else:
         self.getData_PoissonInterval(generatedData_data[iObj],mplot);

            
        canvas_generatedToys_data.append(ROOT.TCanvas("canvas_generatedToys_data_%d"%iObj,""));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].cd();
        ROOT.SetOwnership(canvas_generatedToys_data[len(canvas_generatedToys_data)-1], False);
        pad1 = ROOT.TPad("pad1_%d"%iObj,"pad1_%d"%iObj,0.,0.24,0.99,1. );
        ROOT.SetOwnership(pad1, False);
        pad1.Draw();
        pad2 = ROOT.TPad("pad2_%d"%iObj,"pad2_%d"%iObj,0.,0.,0.99,0.24 );
        ROOT.SetOwnership(pad2, False);
        pad2.Draw();
        pad1.cd();
        mplot.GetXaxis().SetTitleOffset(1.1);
        mplot.GetYaxis().SetTitleOffset(1.3);
        mplot.GetXaxis().SetTitleSize(0.05);
        mplot.GetYaxis().SetTitleSize(0.05);
        mplot.GetXaxis().SetLabelSize(0.045);
        mplot.GetYaxis().SetLabelSize(0.045);
        mplot.Draw();
        banner = self.banner4Plot() ;
        banner.Draw();                                
        pad2.cd();

        if options.pseudodata == 1 :
          mplot_pull = self.get_pull(self.workspace4bias_.var("rrv_mass_lvj"),mplot,generatedData_data[iObj].GetName()+"_curve",fittedPdf_data[iObj].GetName()+"_curve");
        else:
         mplot_pull = self.get_pull(self.workspace4bias_.var("rrv_mass_lvj"),mplot,"data",fittedPdf_data[iObj].GetName()+"_curve");

        mplot_pull.Draw();
        mplot_pull.GetXaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetLabelSize(0.15);
        mplot_pull.GetYaxis().SetTitleSize(0.15);
        mplot_pull.GetYaxis().SetNdivisions(205);

        if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
         os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.root"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));

        pad1.SetLogy();
        pad1.Update();
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*100);
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].Update()

        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        canvas_generatedToys_data[len(canvas_generatedToys_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s_log.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_generatedToys_data[len(canvas_generatedToys_data)-1].GetName()));
        

      ### print plots of the parameters
      for ipar in range(len(parameterHisto_data)):
       canvas_parameters_data.append(ROOT.TCanvas("canvas_parameters_data_%s"%parameterHisto_data[ipar].GetName(),""));
       canvas_parameters_data[len(canvas_parameters_data)-1].cd();
       ROOT.SetOwnership(canvas_parameters_data[len(canvas_parameters_data)-1], False);
       parameterHisto_data[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHisto_data[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHisto_data[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHisto_data[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHisto_data[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHisto_data[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHisto_data[ipar].GetXaxis().SetTitle(parameterHisto_data[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_par_%s"%(parameterHisto_data[ipar].GetName()),"gaus",parameterHisto_data[ipar].GetXaxis().GetXmin(),parameterHisto_data[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHisto_data[ipar].Fit(Gaussian,"MSQ");
       parameterHisto_data[ipar].Draw("E");

       self.outputFilePull.cd();
       Gaussian.Write();
       parameterHisto_data[ipar].Write();

       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

       canvas_parameters_data[len(canvas_parameters_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_data[len(canvas_parameters_data)-1].GetName()));
       canvas_parameters_data[len(canvas_parameters_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_data[len(canvas_parameters_data)-1].GetName()));

      ### print plots of the parameters
      for ipar in range(len(parameterHisto_fraction_data)):
       canvas_parameters_fraction_data.append(ROOT.TCanvas("canvas_parameters_data_%s"%parameterHisto_fraction_data[ipar].GetName(),""));
       canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1].cd();
       ROOT.SetOwnership(canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1], False);
       parameterHisto_fraction_data[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHisto_fraction_data[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHisto_fraction_data[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHisto_fraction_data[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHisto_fraction_data[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHisto_fraction_data[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHisto_fraction_data[ipar].GetXaxis().SetTitle(parameterHisto_fraction_data[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_%s"%(parameterHisto_fraction_data[ipar].GetName()),"gaus",parameterHisto_fraction_data[ipar].GetXaxis().GetXmin(),parameterHisto_fraction_data[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHisto_fraction_data[ipar].Fit(Gaussian,"MSQ");
       parameterHisto_fraction_data[ipar].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

       canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1].GetName()));
       canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_fraction_data[len(canvas_parameters_fraction_data)-1].GetName()));


      ### print the parameters error distribution 
      for ipar in range(len(parameterHistoError_data)):
       canvas_parameters_err_data.append(ROOT.TCanvas("canvas_parameterHistoError_data_%s"%parameterHistoError_data[ipar].GetName(),""));
       canvas_parameters_err_data[len(canvas_parameters_err_data)-1].cd();
       ROOT.SetOwnership(canvas_parameters_err_data[len(canvas_parameters_err_data)-1], False);
       parameterHistoError_data[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoError_data[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoError_data[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoError_data[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoError_data[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoError_data[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoError_data[ipar].GetXaxis().SetTitle(parameterHistoError_data[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_%s"%(parameterHistoError_data[ipar].GetName()),"gaus",parameterHistoError_data[ipar].GetXaxis().GetXmin(),parameterHistoError_data[ipar].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoError_data[ipar].Fit(Gaussian,"MSQ");
       parameterHistoError_data[ipar].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

       canvas_parameters_err_data[len(canvas_parameters_err_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_data[len(canvas_parameters_err_data)-1].GetName()));
       canvas_parameters_err_data[len(canvas_parameters_err_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_err_data[len(canvas_parameters_err_data)-1].GetName()));

      ## print pulls of each parameter
      for ipar in range(len(parameterHistoPull_data)):
       canvas_parameters_pull_data.append(ROOT.TCanvas("canvas_parameters_pull_data_%s"%parameterHistoPull_data[ipar].GetName(),""));
       canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].cd();
       ROOT.SetOwnership(canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1], False);
       parameterHistoPull_data[ipar].GetXaxis().SetTitleOffset(1.1);
       parameterHistoPull_data[ipar].GetYaxis().SetTitleOffset(1.3);
       parameterHistoPull_data[ipar].GetXaxis().SetTitleSize(0.04);
       parameterHistoPull_data[ipar].GetYaxis().SetTitleSize(0.04);
       parameterHistoPull_data[ipar].GetXaxis().SetLabelSize(0.035);
       parameterHistoPull_data[ipar].GetYaxis().SetLabelSize(0.035);
       parameterHistoPull_data[ipar].GetXaxis().SetTitle(parameterHistoPull_data[ipar].GetName());
       Gaussian = ROOT.TF1("Gaussian_pull_%s"%(parameterHistoPull_data[ipar].GetName()),"gaus",-5,5);
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       parameterHistoPull_data[ipar].Fit(Gaussian,"MSQ");
       parameterHistoPull_data[ipar].Draw("E");

       self.outputFilePull.cd();
       Gaussian.Write();
       parameterHistoPull_data[ipar].Write();

       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));

       canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].GetName()));
       canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_parameters_pull_data[len(canvas_parameters_pull_data)-1].GetName()));

      ## print chi2
      if len (chi2distribution_data)!=0: 
       canvas_chi2_data.append(ROOT.TCanvas("canvas_chi2_data_%s"%chi2distribution_data[0].GetName(),""));
       canvas_chi2_data[0].cd();
       ROOT.SetOwnership(canvas_chi2_data[0], False);
       chi2distribution_data[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_data[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_data[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_data[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_data[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_data[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_data[0].GetXaxis().SetTitle(chi2distribution_data[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_data[0].GetXaxis().GetXmin(),chi2distribution_data[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_data[0].Fit(Gaussian,"MSQ");
       chi2distribution_data[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
         
       canvas_chi2_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data[0].GetName()));
       canvas_chi2_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data[0].GetName()));

      ## print -log(l)
      if len(nLLdistribution_data)!=0: 
       canvas_nLL_data.append(ROOT.TCanvas("canvas_nLL_data_%s"%nLLdistribution_data[0].GetName(),""));
       canvas_nLL_data[0].cd();
       ROOT.SetOwnership(canvas_nLL_data[0], False);
       nLLdistribution_data[0].GetXaxis().SetTitleOffset(1.1);
       nLLdistribution_data[0].GetYaxis().SetTitleOffset(1.3);
       nLLdistribution_data[0].GetXaxis().SetTitleSize(0.04);
       nLLdistribution_data[0].GetYaxis().SetTitleSize(0.04);
       nLLdistribution_data[0].GetXaxis().SetLabelSize(0.035);
       nLLdistribution_data[0].GetYaxis().SetLabelSize(0.035);
       nLLdistribution_data[0].GetXaxis().SetTitle(nLLdistribution_data[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",nLLdistribution_data[0].GetXaxis().GetXmin(),nLLdistribution_data[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       nLLdistribution_data[0].Fit(Gaussian,"MSQ");
       nLLdistribution_data[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
       canvas_nLL_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_data[0].GetName()));
       canvas_nLL_data[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_nLL_data[0].GetName()));

      ## print chi2 from frame
      if len(chi2distribution_data_frame) !=0:
       canvas_chi2_data_frame.append(ROOT.TCanvas("canvas_chi2_data_frame_%s"%chi2distribution_data_frame[0].GetName(),""));
       canvas_chi2_data_frame[0].cd();
       ROOT.SetOwnership(canvas_chi2_data_frame[0], False);
       chi2distribution_data_frame[0].GetXaxis().SetTitleOffset(1.1);
       chi2distribution_data_frame[0].GetYaxis().SetTitleOffset(1.3);
       chi2distribution_data_frame[0].GetXaxis().SetTitleSize(0.04);
       chi2distribution_data_frame[0].GetYaxis().SetTitleSize(0.04);
       chi2distribution_data_frame[0].GetXaxis().SetLabelSize(0.035);
       chi2distribution_data_frame[0].GetYaxis().SetLabelSize(0.035);
       chi2distribution_data_frame[0].GetXaxis().SetTitle(chi2distribution_data_frame[0].GetName());
       Gaussian = ROOT.TF1("GaussianChi","gaus",chi2distribution_data_frame[0].GetXaxis().GetXmin(),chi2distribution_data_frame[0].GetXaxis().GetXmax());
       Gaussian.SetLineColor(kBlue);
       Gaussian.SetLineWidth(2);
       chi2distribution_data_frame[0].Fit(Gaussian,"MSQ");
       chi2distribution_data_frame[0].Draw("E");
      
       if not os.path.isdir("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys")):
        os.system("mkdir -p plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
        os.system("mkdir plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"toys"));
        
       canvas_chi2_data_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.pdf"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data_frame[0].GetName()));
       canvas_chi2_data_frame[0].SaveAs("plots_%s_%s_%s_g1/m_lvj_fitting_%s_%s/toys/%s.png"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,canvas_chi2_data_frame[0].GetName()));
 
#### Main code     
if __name__ == "__main__":

  print "###################### begin the analysis: channel %s, signal name %s, mlvj min %s, mlvj max %s, mj min %s, mj max %s, genfunction %s, fitfunction %s"%(options.channel,sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],options.fgen,options.fres);
  
  fitBiasAnalysis = doBiasStudy_mlvj (options.channel,sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),options.fgen,options.fres)

  if options.shapetest == 0 :
     fitBiasAnalysis.biasAnalysis();
     fitBiasAnalysis.outputFilePull.Close();                                  
  else: 
     fitBiasAnalysis.shapeParametrizationAnalysis();


           
  

