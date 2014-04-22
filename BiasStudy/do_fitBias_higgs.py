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
from array import array

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit,RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet,RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, RooMCStudy, RooGlobalFunc,RooChi2MCSModule, RooCurve, RooBernstein, RooHist, RooRandom

############################################
#              Job steering                #
############################################

## general options 
parser = OptionParser()
parser.add_option('-a', '--additioninformation',action = "store",      type="string",dest="additioninformation",default="higgs")
parser.add_option('-b',                         action = 'store_true', dest='noX',   default=False,             help='no X11 windows')
parser.add_option('-f','--inPath',    help='directory with workspace' , default = "./" )

## event setup: pseudodata, only w+jets mc, cateogry and channel
parser.add_option('-m','--mass',       help='test signal yield for this mass', type=int, default=-1)
parser.add_option('-l','--channel',    help='lepton flavor: el or mu  or both , default:both' ,type="string", default ="mu" ) ## lepton flavour
parser.add_option('-p','--category',   help='purity category: LP (low purity) or HP (high purity) or NP (no purity selection), default:HP',type="string", default = "HP")
parser.add_option('-d','--pseudodata', help='use pseudodata instead of real data', type=int, default=0)
parser.add_option('-z','--isMC',       help='options to run pure mc w+jets toys', type=int, default=0)
parser.add_option('--jetBin',          action="store", type="string", dest="jetBin",    default="")

## F-Test or bias study
parser.add_option('-t','--shapetest',  help='make W+jets and data fit with different parametrization', type=int, default=0)

## Bias study info and options
parser.add_option('-n','--nexp',      help='number of toys', type=int, default=1000)
parser.add_option('-g','--fgen',      help='function to generate toys Exp,ExpTail,Pow2,ExpN)', type="string", default="ExpN")
parser.add_option('-r','--fres',      help='function to fit toys (Exp,ExpTail,Pow2,ExpN)',     type="string", default="ExpN")
parser.add_option('-s','--storeplot', help='in case of more than 10 toys just 1/3 stored, more than 100 1/10',     type=int, default=0)
parser.add_option('-k','--ttbarcontrolregion',  help='run the toy or f-test analysis in the ttbar control region', type=int, default=0)
parser.add_option('-y','--mlvjregion',  help='run the toy taking the events in the : low sb, high sb or signal region', type="string", default="_sb_lo")
parser.add_option('-q','--fitjetmass',  help='run the toy on the jet mass fit', type=int, default=0)
parser.add_option('-w','--onlybackgroundfit',  help='run only background fit',  type=int, default=0)
parser.add_option('-i','--inflatejobstatistic',  help='enlarge the generated statistics in the fit',  type=int, default=1)
parser.add_option('--scalesignalwidth', help='reduce the signal width by a factor x', type=float, default=1.)
parser.add_option('--injectSingalStrenght', help='inject a singal in the toy generation', type=float, default=1.)

(options, args) = parser.parse_args()

### Lybrary import

ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/MakePdf_cxx.so")
ROOT.gSystem.Load(options.inPath+"/BiasStudy/BiasUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/FitUtils/FitUtils_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error

from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooPow3Pdf, RooErfPow3Pdf, RooUser1Pdf

from ROOT import biasModelAnalysis, MakeGeneralPdf, MakeExtendedModel, get_TTbar_mj_Model, get_STop_mj_Model, get_VV_mj_Model, get_WW_EWK_mj_Model, get_WJets_mj_Model, get_ggH_mj_Model, get_vbfH_mj_Model, get_TTbar_mlvj_Model, get_STop_mlvj_Model, get_VV_mlvj_Model, get_WW_EWK_mlvj_Model, get_WJets_mlvj_Model, get_ggH_mlvj_Model, get_vbfH_mlvj_Model, fix_Model,  clone_Model

from ROOT import setTDRStyle, get_pull, draw_canvas, draw_canvas_with_pull, legend4Plot, GetDataPoissonInterval, GetLumi

from ROOT import fit_mj_single_MC, fit_mlvj_model_single_MC, fit_WJetsNormalization_in_Mj_signal_region, fit_mlvj_in_Mj_sideband, get_WJets_mlvj_correction_sb_lo_to_signal_region

from ROOT import *

gInterpreter.GenerateDictionary("std::map<std::string,std::string>", "map;string;string")

###############################
## doFit Class Implemetation ##
###############################

class doBiasStudy_mlvj:

    def __init__(self,in_channel,in_ggH_sample,in_mlvj_min=400., in_mlvj_max=1400., in_mj_min=40, in_mj_max=130, generation_model="ErfExp", fit_model="ErfExp", input_workspace=None):

        setTDRStyle(); #import the tdr style

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;
        ROOT.RooRandom.randomGenerator().SetSeed(0);

        ### shapes to be used in mj
        self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string) () ;
        self.mj_shape["TTbar"]   = "2Gaus_ErfExp";
        self.mj_shape["VV"]      = "2_2Gaus";
        self.mj_shape["WW_EWK"]  = "2_2Gaus";
        if options.jetBin == "_2jet":  self.mj_shape["STop"]    = "ErfExp";
        else: self.mj_shape["STop"]    = "2Gaus_ErfExp";
        self.mj_shape["WJets0"]  = "ErfExp";
        self.mj_shape["WJets1"]  = "ErfExp";
        self.mj_shape["WJets01"] = "User1";        
        self.mj_shape["ggH"]     = "2Gaus";
        self.mj_shape["vbfH"]    = "2Gaus";

        self.mlvj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string) () ;
        self.mlvj_shape["TTbar"]   = generation_model;
        self.mlvj_shape["VV"]      = generation_model;
        self.mlvj_shape["WW_EWK"]  = generation_model;
        self.mlvj_shape["STop"]    = generation_model;
        self.mlvj_shape["WJets0"]  = generation_model;
        self.mlvj_shape["WJets1"]  = generation_model;
        self.mlvj_shape["WJets01"] = generation_model;
        self.mlvj_shape["ggH"]     = "CB_v1";
        self.mlvj_shape["vbfH"]    = "CB_v1";

        self.tmpFile = TFile("tmp.root","RECREATE");
        self.tmpFile.cd();
                                                                                                                                                                
        print "###################### construnctor ############################# ";
        ### set the channel type --> electron or muon
        self.channel = in_channel;
        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;

        self.BinWidth_mj = 5.;
        nbins_mj         = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max        = in_mj_min+nbins_mj*self.BinWidth_mj;

        self.BinWidth_mlvj = 50.;
        self.in_mlvj_min   = in_mlvj_min;
        self.nbins_mlvj    = int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        self.in_mlvj_max   = in_mlvj_min+self.nbins_mlvj*self.BinWidth_mlvj;

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

        suffix = "";
        suffix = "_"+options.channel ;
        if options.jetBin == "_2jet": suffix = suffix + "_2jet";
        if options.ttbarcontrolregion : suffix = suffix+"_ttbar";
        if options.fitjetmass         : suffix = suffix+"_jetmass";
        if in_mlvj_min < 550          : suffix = suffix+"_turnOn";
        if options.scalesignalwidth !=1 : suffix = suffix+ ("_width_%0.1f")%(options.scalesignalwidth);
        if options.injectSingalStrenght !=0 : suffix = suffix+"_SB";
        else :suffix = suffix+"_B";         
        if options.onlybackgroundfit  : suffix = suffix+"_B";        
        else: suffix = suffix+"_SB";
        
        ## create the workspace and import them
        if input_workspace is None:
             self.workspace4bias_ = RooWorkspace("workspace4bias_%s_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model,suffix),"workspace4bias_%s_%s_%s_%s_%s"%(self.channel,self.wtagger_label,self.generation_model,self.fit_model,suffix));
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
        if in_ggH_sample == "ggH600":  self.vbfhiggs_sample = "vbfH600";
        if in_ggH_sample == "ggH700":  self.vbfhiggs_sample = "vbfH700";
        if in_ggH_sample == "ggH800":  self.vbfhiggs_sample = "vbfH800";
        if in_ggH_sample == "ggH900":  self.vbfhiggs_sample = "vbfH900";
        if in_ggH_sample == "ggH1000": self.vbfhiggs_sample = "vbfH1000";
        if in_ggH_sample == "ggH1500": self.vbfhiggs_sample = "vbfH1500";
        if in_ggH_sample == "ggH2000": self.vbfhiggs_sample = "vbfH2000";

        if options.pseudodata :
           self.file_data  = "ofile_pseudodata.root";                                                                                
        else:  
           self.file_data  = "ofile_data.root";

        self.file_ggH   = ("ofile_%s.root"%(self.ggH_sample));
        self.file_vbfH  = ("ofile_%s.root"%(self.vbfhiggs_sample));

        #WJets0 is the default PS model, WJets1 is the alternative PS model                                                                                       
        if options.jetBin == "_2jet" :
         self.file_WJets0_mc = ("ofile_WJets_exclusive_Pythia.root");
         self.file_WJets1_mc = ("ofile_WJets_Herwig.root");
        else:
         self.file_WJets0_mc = ("ofile_WJets_Pythia100.root");
         self.file_WJets1_mc = ("ofile_WJets_Herwig.root");
            

        self.file_VV_mc     = ("ofile_VV.root");# WW+WZ
        self.file_WW_EWK_mc = ("ofile_WW2jet_phantom.root");# WW+WZ
        
        self.file_TTbar_mc         = ("ofile_TTbar_Powheg.root");
        self.file_TTbar_mcanlo_mc  = ("ofile_TTbar_mcanlo.root");
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
            self.wtagger_cut=0.75 ; self.wtagger_cut_min=0.5 ;
        if self.wtagger_label=="NP":
            elf.wtagger_cut=10000;

        ## color palette
        self.color_palet = ROOT.std.map(ROOT.std.string, int) () ;
        self.color_palet["data"]   = 1;
        self.color_palet["WJets"]  = 2;
        self.color_palet["VV"]     = 4;
        self.color_palet["WW_EWK"] = 6;
        self.color_palet["STop"]   = 7;
        self.color_palet["TTbar"]  = 210;
        self.color_palet["ggH"]    = 1;
        self.color_palet["vbfH"]   = 12;
        self.color_palet["Signal"] = 1;
        self.color_palet["Uncertainty"] = 1;
        self.color_palet["Other_Backgrounds"] = 1;

        ## for basic selection
        self.vpt_cut   = 200;
        self.pfMET_cut = 50;
        self.lpt_cut   = 30;

        if self.channel=="el":
            self.pfMET_cut = 70; self.lpt_cut = 35;#very tight
        self.deltaPhi_METj_cut = 2.0;

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
           elif options.pseudodata == 0 and options.jetBin == "_2jet":
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.128);
             self.rrv_wtagger_eff_reweight_forT.setError(0.338);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.91);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           elif options.pseudodata == 0 and not  options.jetBin == "_2jet": 
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.93);
             self.rrv_wtagger_eff_reweight_forT.setError(0.06*self.rrv_wtagger_eff_reweight_forT.getVal());
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.91);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());

        if self.channel=="el" and self.wtagger_label=="HP":
           if options.pseudodata == 1:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           elif options.pseudodata == 0 and options.jetBin == "_2jet":
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.96);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.94);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           elif options.pseudodata == 0 and not  options.jetBin == "_2jet": 
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.93);
             self.rrv_wtagger_eff_reweight_forT.setError(0.06*self.rrv_wtagger_eff_reweight_forT.getVal());
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.94);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());

        if self.channel=="em" and self.wtagger_label=="HP":
           if options.pseudodata == 1:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
           else:
             self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.128);
             self.rrv_wtagger_eff_reweight_forT.setError(0.369);
             self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
             self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());

        if options.pseudodata == 1:
            self.mean_shift = 1.; self.sigma_scale = 1.;
        else:    
            self.mean_shift = 1.4; self.sigma_scale = 1.11;
        
        ### create output txt file for the shape test 
        if options.isMC == 1 and options.shapetest == 1:
         if options.ttbarcontrolregion == 0:    
          self.file_FTestFile_txt = "FTest%s_%s_%s_mc.txt"%(options.jetBin,self.channel,self.wtagger_label);
         else:
          self.file_FTestFile_txt = "FTest%s_%s_%s_mc_ttbar.txt"%(options.jetBin,self.channel,self.wtagger_label);
         os.system("rm "+self.file_FTestFile_txt); 
        elif options.shapetest == 1 and options.isMC == 0 :
         if options.ttbarcontrolregion == 0 :   
          self.file_FTestFile_txt = "FTest%s_%s_%s_data.txt"%(options.jetBin,self.channel,self.wtagger_label);
         else:
          self.file_FTestFile_txt = "FTest%s_%s_%s_data_ttbar.txt"%(options.jetBin,self.channel,self.wtagger_label);         
         os.system("rm "+self.file_FTestFile_txt); 
 
          
        ### create output root file for the pull plot vs mass
        if options.shapetest  == 0:
            self.outputFile  = ROOT.TFile("output_%s_%s_%s%s.root"%(self.ggH_sample,options.fgen,options.fres,suffix),"RECREATE");
            self.outputTree  = ROOT.TTree("otree","otree");
            self.outputFile.cd();
   	    
    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass ="jet_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");

        rrv_mass_j    = self.workspace4bias_.var("rrv_mass_j") 
        rrv_mass_lvj  = self.workspace4bias_.var("rrv_mass_lvj")
        rrv_weight    = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 
         
	rdataset_mj      = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4bias_mj = RooDataSet("rdataset4bias"+label+"_"+self.channel+"_mj","rdataset4bias"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        #dataset of m_lvj -> before and after vbf cuts -> central object value
        rdataset_sb_lo_mlvj  = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_sb_hi_mlvj  = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        rdataset4bias_sb_lo_mlvj     = RooDataSet("rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4bias_sb_hi_mlvj     = RooDataSet("rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4bias_signal_region_mlvj = RooDataSet("rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4bias"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
 
        ###### Define the event categorization
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");

        combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4bias = RooDataSet("combData4bias"+label+"_"+self.channel,"combData4bias"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        print "N entries: ", treeIn.GetEntries();
        
	hnum_4region = TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j   -1: 
	
        for i in range(treeIn.GetEntries()):

          if i % 100000 == 0: print "iEvent: ",i
          treeIn.GetEntry(i);

          discriminantCut = False;

          wtagger=-1;
          njet = 0. ; tmp_vbf_dEta =0.; tmp_vbf_Mjj = 0.; ungroomed_jet_pt = 0.; pfMET = 0.; mass_lvj = 0. ;

          jet_1 = ROOT.TLorentzVector();
          jet_2 = ROOT.TLorentzVector();

          tmp_scale_to_lumi = treeIn.wSampleWeight;
                                
          if options.ttbarcontrolregion == 0 :
           wtagger = getattr(treeIn,"jet_tau2tau1");
          
           if wtagger < self.wtagger_cut:
            discriminantCut = True;
           else:
            discriminantCut = False;
                  
           # jet mass , central value
           tmp_jet_mass = getattr(treeIn, jet_mass);
           tmp_vbf_dEta = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta")-getattr(treeIn,"vbf_maxpt_j2_eta"));
           tmp_vbf_Mjj  = getattr(treeIn, "vbf_maxpt_jj_m");
           njet         = getattr(treeIn,"numberJetBin");
           ungroomed_jet_pt = getattr(treeIn,"ungroomed_jet_pt");
           pfMET    = getattr(treeIn,"pfMET");
           mass_lvj = getattr(treeIn,"mass_lvj_type0_met");

          else:

           wtagger = getattr(treeIn,"ttb_ca8_tau2tau1");
          
           if wtagger < self.wtagger_cut:
            discriminantCut = True;
           else:
            discriminantCut = False;
                  
           # jet mass , central value
           tmp_jet_mass = getattr(treeIn, jet_mass);
           tmp_vbf_dEta = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta")-getattr(treeIn,"vbf_maxpt_j2_eta"));
           tmp_vbf_Mjj  = getattr(treeIn, "vbf_maxpt_jj_m");
           njet         = getattr(treeIn,"numberJetBin");
           ungroomed_jet_pt = getattr(treeIn,"ttb_ca8_ungroomed_pt");
           pfMET    = getattr(treeIn,"pfMET");
           mass_lvj = getattr(treeIn,"mass_lvj_type0_met");
          
          isPassingCut = 0 ;

          if options.jetBin == "_2jet": 

           if options.ttbarcontrolregion == 0 or TString(in_file_name).Contains("ggH") or TString(in_file_name).Contains("vbfH"):
          
            if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
             isPassingCut = 1 ;

           else:

            if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and ( getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") > 0.679 or getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV") > 0.679 ) and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
             isPassingCut = 1 ;
             

          else:

           if options.ttbarcontrolregion == 0 or TString(in_file_name).Contains("ggH") or TString(in_file_name).Contains("vbfH"):

            if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet < 2:
             isPassingCut = 1 ;
          
           else:

            if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet < 2:
             isPassingCut = 1 ;
          
              
          tmp_event_weight = 0 ;       
          tmp_event_weight4bias = 0 ;

          if isPassingCut !=0 :
              
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

             if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max and isPassingCut == 1:
                 rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                 data_category.setLabel("sideband");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);

             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max and isPassingCut == 1:
                 rdataset_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                 data_category.setLabel("signal_region");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4bias.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4bias);
                   
             if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max and isPassingCut == 1:
                 rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4bias_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4bias );
                                     
	     if isPassingCut == 1: 
              rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
              rdataset4bias_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4bias );
             
	      if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max: 
                  hnum_4region.Fill(-1,tmp_event_weight );
              if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max : 
                  hnum_4region.Fill(0,tmp_event_weight);
              if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max: 
                  hnum_4region.Fill(1,tmp_event_weight);
              hnum_4region.Fill(2,tmp_event_weight);
              
        print "########### Nominal Value ###########";        
        rrv_scale_to_lumi            = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,(rdataset_sb_lo_mlvj.sumEntries()+rdataset_signal_region_mlvj.sumEntries()+rdataset_sb_hi_mlvj.sumEntries())/(rdataset4bias_sb_lo_mlvj.sumEntries()+rdataset4bias_sb_hi_mlvj.sumEntries()+rdataset4bias_signal_region_mlvj.sumEntries()));
        rrv_scale_to_lumi.Print();
        getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi);  

        if rdataset4bias_sb_lo_mlvj.sumEntries() !=0 :
         rrv_scale_to_lumi_sb_lo_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj",(rdataset_sb_lo_mlvj.sumEntries())/(rdataset4bias_sb_lo_mlvj.sumEntries()));
         rrv_scale_to_lumi_sb_lo_mlvj.Print();
         getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_sb_lo_mlvj);  

        if rdataset4bias_sb_hi_mlvj.sumEntries()!=0 :
         rrv_scale_to_lumi_sb_hi_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj",(rdataset_sb_hi_mlvj.sumEntries())/(rdataset4bias_sb_hi_mlvj.sumEntries()));
         rrv_scale_to_lumi_sb_hi_mlvj.Print();
         getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_sb_hi_mlvj);               
 
        if rdataset4bias_signal_region_mlvj.sumEntries() !=0:
         rrv_scale_to_lumi_signal_region_mlvj = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj",(rdataset_signal_region_mlvj.sumEntries())/(rdataset4bias_signal_region_mlvj.sumEntries()));
         rrv_scale_to_lumi_signal_region_mlvj.Print();        
         getattr(self.workspace4bias_,"import")(rrv_scale_to_lumi_signal_region_mlvj);               
                              
        getattr(self.workspace4bias_,"import")(rdataset_sb_lo_mlvj); rdataset_sb_lo_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset_signal_region_mlvj); rdataset_signal_region_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset_sb_hi_mlvj); rdataset_sb_hi_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_lo_mlvj); rdataset4bias_sb_lo_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_signal_region_mlvj); rdataset4bias_signal_region_mlvj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_sb_hi_mlvj); rdataset4bias_sb_hi_mlvj.Print();
        getattr(self.workspace4bias_,"import")(combData); combData.Print();
        getattr(self.workspace4bias_,"import")(combData4bias); combData4bias.Print();
       
        print "########### nominal value ###########";
        rrv_number_dataset_sb_lo_mj = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));        
        rrv_number_dataset_sb_hi_mj = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));        


        getattr(self.workspace4bias_,"import")(rdataset_mj); rdataset_mj.Print();
        getattr(self.workspace4bias_,"import")(rdataset4bias_mj); rdataset4bias_mj.Print();
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_sb_lo_mj); rrv_number_dataset_sb_lo_mj.Print();
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_signal_region_mj); rrv_number_dataset_signal_region_mj.Print();
        getattr(self.workspace4bias_,"import")(rrv_number_dataset_sb_hi_mj); rrv_number_dataset_sb_hi_mj.Print();

    #########################
    #### Code for F-test ####
    #########################    

    def shapeParametrizationAnalysis(self,label = "_WJets0"):

     if label == "_WJets0": fileName = self.file_WJets0_mc ;
     else:                  fileName = self.file_TTbar_mc ;
     
     mlvj_region = options.mlvjregion;

     if options.isMC == 1:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(fileName,label)# to get the shape of m_lvj                                                                                             

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfExp_v1",self.channel,self.wtagger_label,0,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfExp_v1",self.channel,self.wtagger_label,0,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfExpTail",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfExpTail",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfExp_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfExp_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Erf2Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Erf2Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfPow_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfPow_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfPow2_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfPow2_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfPow3_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfPow3_v1",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Erf2Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Erf2Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfChebychev_v2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfChebychev_v2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfChebychev_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfChebychev_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ErfChebychev_v4",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ErfChebychev_v4",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

      else:

       self.get_mj_and_mlvj_dataset(fileName,label)# to get the shape of m_lvj                                                                                             

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","ExpTail",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","ExpTail",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Exp_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
#       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Exp_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
#       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","2Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","2Exp",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Pow2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Pow2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Pow3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Pow3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","2Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","2Pow",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Chebychev_v2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Chebychev_v2",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Chebychev_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Chebychev_v3",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_sb_lo","Chebychev_v4",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_model_single_MC(self.workspace4bias_,fileName,label,"_signal_region","Chebychev_v4",self.channel,self.wtagger_label,1,0,0,0,label,self.file_FTestFile_txt);
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

     else:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_VV_mc,"_VV",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,1,0,0,0,"_VV");
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                               
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,1,0,0,0,"_STop");
           
       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,1,0,0,0,"_TTbar");

       if options.jetBin == "_2jet":    
        self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj
        fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WW_EWK_mc,"_WW_EWK",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,1,0,0,0,"_WW_EWK");
           
       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,1,0,0,0,"_WJets0");

       self.get_mj_and_mlvj_dataset(self.file_data,"_data"); ## global fit of data in the sidand fixing non dominant bkg

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfExp_v1",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfExpTail",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
#       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfExp_v3",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Erf2Exp",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfPow_v1",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfPow2_v1",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfPow3_v1",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Erf2Pow",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfChebychev_v2",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfChebychev_v3",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfChebychev_v4",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ErfPowExp_v1",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

      else: 

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_VV_mc,"_VV",mlvj_region,"Exp",self.channel,self.wtagger_label,1,0,0,0,"_VV");
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop");
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",mlvj_region,"Exp",self.channel,self.wtagger_label,1,0,0,0,"_STop");

       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar");
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",mlvj_region,"Exp",self.channel,self.wtagger_label,1,0,0,0,"_TTbar");

       if options.jetBin == "_2jet":    
        self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr");
        fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WW_EWK_mc,"_WW_EWK",mlvj_region,"Exp",self.channel,self.wtagger_label,1,0,0,0,"_WW_EWK");

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0");
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",mlvj_region,"Exp",self.channel,self.wtagger_label,1,0,0,0,"_WJets0");

       self.get_mj_and_mlvj_dataset(self.file_data,"_data");
       
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Exp",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"ExpTail",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       #fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Exp_v3",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"2Exp",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Pow",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Pow2",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Pow3",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"2Pow",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);

       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Chebychev_v2",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       #fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Chebychev_v3",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);
       fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",mlvj_region,"Chebychev_v4",self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,self.file_FTestFile_txt);


    ############################
    #### Code for Bias test ####
    ############################    

    def biasAnalysis(self,label="_WJets0",fitjetmass = 0):
   
     print"######################## begin the bias analysis ###########################";  
     if fitjetmass == 1:
         options.mlvjregion = "";
     if fitjetmass == 1 and options.isMC == 1:
         print " cannot use MC to perform bias test for the mJ fit -> options not provided" ;
         return ;

     ##### fix signal and bkg models that are going to be used in the generation
     if fitjetmass :
         spectrum = "_mj"   ;
         signal_region = "";
         signal_model  = "2Gaus";
     else:
         spectrum = "_mlvj" ;
         signal_region = "_signal_region";
         signal_model  = "CB_v1";

     
     ###### get the signal and fit it     
     self.get_mj_and_mlvj_dataset(self.file_ggH ,"_%s"%(self.ggH_sample),"jet_mass_pr")# to get the shape of m_lvj
     self.get_mj_and_mlvj_dataset(self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"jet_mass_pr")# to get the shape of m_lvj
     self.workspace4bias_.writeToFile(self.tmpFile.GetName());

     if fitjetmass:
        fit_mj_single_MC(self.workspace4bias_,self.file_ggH,"_%s"%(self.ggH_sample),self.mj_shape["ggH"],self.channel,self.wtagger_label,1);                
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());
        fit_mj_single_MC(self.workspace4bias_,self.file_vbfH,"_%s"%(self.vbfhiggs_sample),self.mj_shape["vbfH"],self.channel,self.wtagger_label,1);                
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());
     else:
        fit_mlvj_model_single_MC(self.workspace4bias_,self.file_ggH,"_%s"%(self.ggH_sample),"_signal_region",self.mlvj_shape["ggH"],self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.ggH_sample)); 
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());
        fit_mlvj_model_single_MC(self.workspace4bias_,self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"_signal_region",self.mlvj_shape["vbfH"],self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));   
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());
     
     ####### Monte Carlo Analysis
     if options.isMC and options.ttbarcontrolregion:
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_TTbar");   
      self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      if options.fgen != options.fres:
          fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_TTbar");   
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());

     elif options.isMC and not options.ttbarcontrolregion:    
      self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
      fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");   
      self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      if options.fgen != options.fres:
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");   
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());

      if options.mlvjregion == "_signal_region":
       fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0","_sb_lo",options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");   
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       ## calculate the alpha factor
       #get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4bias_,label,options.fgen,spectrum,"4bias",self.channel,self.wtagger_label,1); 
       self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       if options.fgen != options.fres:
        fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0","_sb_lo",options.fres,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");   
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());
        ## calculate the alpha factor
        #get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4bias_,label,options.fres,spectrum,"4bias",self.channel,self.wtagger_label,1); 
        self.workspace4bias_.writeToFile(self.tmpFile.GetName());

     elif not options.isMC:
      
      ###### get diboson and fit it 
      self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV","jet_mass_pr");      
      if fitjetmass:
          fit_mj_single_MC(self.workspace4bias_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label,1);                
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());   
      else:
          fit_mlvj_model_single_MC(self.workspace4bias_,self.file_VV_mc,"_VV",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_VV"); 
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());   
          if options.fgen != options.fres:
           fit_mlvj_model_single_MC(self.workspace4bias_,self.file_VV_mc,"_VV",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_VV"); 
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());
   

      ####### get SingleTop and fit it
      self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                                                                 
      if fitjetmass: 
	  if options.ttbarcontrolregion:
             self.mj_shape["STop"] = "ErfExpGaus_sp" ; 
             fit_mj_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",self.mj_shape["STop"],self.channel,self.wtagger_label,1);                
             self.workspace4bias_.writeToFile(self.tmpFile.GetName());
	  else:
             fit_mj_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",self.mj_shape["STop"],self.channel,self.wtagger_label,1);                
             self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      else:
          fit_mlvj_model_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_STop"); 
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
          if options.fgen != options.fres:  
           fit_mlvj_model_single_MC(self.workspace4bias_,self.file_STop_mc,"_STop",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_STop"); 
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());

      ######## get WW EWK and fit it in the sb
      if options.jetBin == "_2jet":    
     
       self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                                             
       if fitjetmass:
          fit_mj_single_MC(self.workspace4bias_,self.file_WW_EWK_mc,"_WW_EWK",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1);                
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
       else:
          fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WW_EWK_mc,"_WW_EWK",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK"); 
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
          if options.fgen != options.fres: 
           fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WW_EWK_mc,"_WW_EWK",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK"); 
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());

      ###### get WJets and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
      if fitjetmass: 
           fit_mj_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.fgen,self.channel,self.wtagger_label,1);                
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());
           if options.fgen != options.fres:
            fit_mj_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.fres,self.channel,self.wtagger_label,1);                
            self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      else: 
           fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_WJets0"); 
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());
           if options.fgen != options.fres:
            fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_WJets0"); 
            self.workspace4bias_.writeToFile(self.tmpFile.GetName());

           if options.mlvjregion == "_signal_region":
             fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0","_sb_lo",options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_WJets0"); 
             self.workspace4bias_.writeToFile(self.tmpFile.GetName());
             ## calculate the alpha factor
             #get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4bias_,label,options.fgen,spectrum,"4bias",self.channel,self.wtagger_label,1); 
             self.workspace4bias_.writeToFile(self.tmpFile.GetName());
             if options.fgen != options.fres:
               fit_mlvj_model_single_MC(self.workspace4bias_,self.file_WJets0_mc,"_WJets0","_sb_lo",options.fres,self.channel,self.wtagger_label,0,0,0,0,"_WJets0"); 
               self.workspace4bias_.writeToFile(self.tmpFile.GetName());
               ## calculate the alpha factor
               #get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4bias_,label,options.fres,spectrum,"4bias",self.channel,self.wtagger_label,1); 
               self.workspace4bias_.writeToFile(self.tmpFile.GetName());

	
      ######## get TTbar and fit it
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      if fitjetmass: 
        if options.ttbarcontrolregion :    
          fit_mj_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.fgen,self.channel,self.wtagger_label,1);                
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
          if options.fgen != options.fres:
           fit_mj_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.fres,self.channel,self.wtagger_label,1);                
           self.workspace4bias_.writeToFile(self.tmpFile.GetName());
        else:
          fit_mj_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1);
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      else:
          fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,0,0,0,0,"_TTbar"); 
          self.workspace4bias_.writeToFile(self.tmpFile.GetName());
          if options.fgen != options.fres:
              fit_mlvj_model_single_MC(self.workspace4bias_,self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fres,self.channel,self.wtagger_label,0,0,0,0,"_TTbar"); 
              self.workspace4bias_.writeToFile(self.tmpFile.GetName());

      ##### get data in sb and fit it                                                                                                                                              
      self.get_mj_and_mlvj_dataset(self.file_data,"_data", "jet_mass_pr"); ## global fit of data in the sidand fixing non dominant bkg
      if fitjetmass :
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4bias_,self.color_palet,self.mj_shape,label,"",options.fgen,self.channel,self.wtagger_label,options.ttbarcontrolregion,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution                 
         self.workspace4bias_.writeToFile(self.tmpFile.GetName());
      else:
         self.workspace4bias_.Print(); 
         fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",options.mlvjregion,options.fgen,self.channel,self.wtagger_label,options.ttbarcontrolregion,options.pseudodata,1,options.jetBin,"");
         self.workspace4bias_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_in_Mj_sideband(self.workspace4bias_,self.color_palet,self.mlvj_shape,label,"",options.mlvjregion,options.fres,self.channel,self.wtagger_label,options.ttbarcontrolregion,1,options.pseudodata,options.jetBin,"");
         self.workspace4bias_.writeToFile(self.tmpFile.GetName());

     ####### fix the backgrund models for the generation
     print "#############################################################################################";
     print "################ Begin of the toy analysis -> fix Pdf to what is pre-fitted #################";
     print "#############################################################################################";

     if options.fitjetmass : 
      fix_Model(self.workspace4bias_,"_%s"%(self.ggH_sample),signal_region,spectrum,self.mj_shape["ggH"],self.channel,"",0);
      fix_Model(self.workspace4bias_,"_%s"%(self.vbfhiggs_sample),signal_region,spectrum,self.mj_shape["vbfH"],self.channel,"",0);
     else: 
      fix_Model(self.workspace4bias_,"_%s"%(self.ggH_sample),signal_region,spectrum,self.mlvj_shape["ggH"],self.channel,"",0);
      fix_Model(self.workspace4bias_,"_%s"%(self.vbfhiggs_sample),signal_region,spectrum,self.mlvj_shape["vbfH"],self.channel,"",0);
            
     if options.isMC == 0  and options.ttbarcontrolregion == 0 and not options.fitjetmass :

      fix_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,spectrum,options.fres,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fres,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fres,self.channel,"",0);

      if options.jetBin == "_2jet":    
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.fres,self.channel,"",0);

      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      if options.fgen != options.fres : 
       fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel,"",0);

     elif options.isMC == 0  and options.ttbarcontrolregion == 1 and not options.fitjetmass:

      fix_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,spectrum,options.fres,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fres,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fres,self.channel,"",0);
 
      if options.jetBin == "_2jet":    
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.gen,self.channel,"",0);
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.fres,self.channel,"",0);

      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      if options.fgen != options.fres : 
       fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel,"",0);
 
     elif options.fitjetmass and not options.ttbarcontrolregion:

      fix_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,spectrum,self.mj_shape["TTbar"],self.channel,"",0);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,self.mj_shape["STop"],self.channel,"",0);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,self.mj_shape["VV"],self.channel,"",0);

      if options.jetBin == "_2jet":    
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,self.mj_shape["WW_EWK"],self.channel,"",0);
 
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel,"",0);
      if options.fgen != options.fres : 
       fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel,"",0);
 
     elif options.fitjetmass and options.ttbarcontrolregion:

      fix_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,"",self.channel);

      if options.jetBin == "_2jet":
       fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,"",self.channel);

      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel);
      if options.fgen != options.fres : 
       fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel);
      
     print "#########################################################";
     print "################ Build the signal model #################";
     print "#########################################################";
     
     ### clone the signal shape --> parameter already fixed
     fitted_signal_ggH   = self.workspace4bias_.pdf("model_%s_%s_%s%s"%(self.ggH_sample,signal_region,self.channel,spectrum));
     fitted_signal_vbfH  = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.vbfhiggs_sample,signal_region,self.channel,spectrum));

     ### make the pdf for ggH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_ggH = RooArgList();
      
     if fitjetmass:

      model_signal_ggH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.ggH_sample,signal_region+self.mj_shape["ggH"]),self.mj_shape["ggH"],spectrum,self.wtagger_label,self.channel);  
      model_signal_ggH.Print();     
      clone_Model(self.workspace4bias_,model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,self.mj_shape["ggH"],self.channel);      
      getattr(self.workspace4bias_,"import")(model_signal_ggH); 
    
      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      print "rrv_number_"+self.ggH_sample+signal_region+self.mj_shape["ggH"]+self.channel+spectrum ; 
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+self.mj_shape["ggH"]+"_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+self.mj_shape["ggH"]+"_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);

     else:
      model_signal_ggH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.ggH_sample,signal_region+self.mlvj_shape["ggH"]),self.mlvj_shape["ggH"],spectrum,self.wtagger_label,self.channel);  
      model_signal_ggH.Print();     
 
      clone_Model(self.workspace4bias_,model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,self.mlvj_shape["ggH"],self.channel);      
      getattr(self.workspace4bias_,"import")(model_signal_ggH); 

      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+self.mlvj_shape["ggH"]+"_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+self.mlvj_shape["ggH"]+"_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);
     
      
     ### make the pdf for vbfH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_vbfH = RooArgList();

     if fitjetmass:
      model_signal_vbfH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.vbfhiggs_sample,signal_region+self.mj_shape["vbfH"]),self.mj_shape["vbfH"],spectrum,self.wtagger_label,self.channel);  
      model_signal_vbfH.Print();        

      clone_Model(self.workspace4bias_,model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,self.mj_shape["vbfH"],self.channel);      
      getattr(self.workspace4bias_,"import")(model_signal_vbfH); 

      rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+self.mj_shape["vbfH"]+"_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+self.mj_shape["vbfH"]+"_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_vbfH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_vbfH);

     else: 
      model_signal_vbfH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.vbfhiggs_sample,signal_region+self.mlvj_shape["vbfH"]),self.mlvj_shape["vbfH"],spectrum,self.wtagger_label,self.channel);  
      model_signal_vbfH.Print();        

      clone_Model(self.workspace4bias_,model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,self.mlvj_shape["vbfH"],self.channel);      
      getattr(self.workspace4bias_,"import")(model_signal_vbfH); 

      rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+self.mlvj_shape["vbfH"]+"_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+self.mlvj_shape["vbfH"]+"_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_vbfH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_vbfH);

     ## total pdf: ggH + vbfH --> total higgs sum --> relative fraction fixed to the expectation after all the cuts
     rrv_fraction_ggH_vbf = RooRealVar("rrv_fraction_ggH_vbf","rrv_fraction_ggH_vbf",rrv_number_signal_signal_fit_vbfH.getVal()/rrv_number_signal_signal_fit_ggH.getVal());
     rrv_fraction_ggH_vbf.setConstant(kTRUE);
     rrv_fraction_ggH_vbf.Print();
     getattr(self.workspace4bias_,"import")(rrv_fraction_ggH_vbf);

     model_modified_signal   = RooAddPdf("model_modified_signal","model_modified_signal",RooArgList(model_signal_vbfH,model_signal_ggH),RooArgList(rrv_fraction_ggH_vbf),1);
     model_modified_signal.Print();
     getattr(self.workspace4bias_,"import")(model_modified_signal);
  
     rrv_number_signal_region_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH_vbfH","rrv_number_signal_region_fit_ggH_vbfH",0,-1e8,1e8);
     rrv_number_signal_region_fit_ggH.Print();
     getattr(self.workspace4bias_,"import")(rrv_number_signal_region_fit_ggH);
     
     model_total_signal = RooExtendPdf("model_higgs_signal_region_fit_%s%s"%(self.channel,spectrum),
                                       "model_higgs_signal_region_fit_%s%s"%(self.channel,spectrum),model_modified_signal,rrv_number_signal_region_fit_ggH);

     model_total_signal.Print();
     getattr(self.workspace4bias_,"import")(model_total_signal);
        
     ############### Make the MC analysis --> make the Entended pdf for the bkg
     if options.isMC == 1 :

      print "#################################################################################";
      print "################ Start the MC analysis -> bkg model in the toy ##################";
      print "#################################################################################";

      constraintlist_bkg_wjet = RooArgList();
      
      model_bkg_wjet = MakeExtendedModel(self.workspace4bias_,label+options.mlvjregion+"_fit",options.fres,spectrum,self.channel,self.wtagger_label,constraintlist_bkg_wjet,1);
      model_bkg_wjet.Print();
      getattr(self.workspace4bias_,"import")(model_bkg_wjet); 

      if options.fres == options.fgen :
       clone_Model(self.workspace4bias_,model_bkg_wjet,label,options.mlvjregion,spectrum,options.fgen,self.channel,0);
          
      self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+"_fit_"+self.channel+spectrum).setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_"+self.channel+spectrum).getVal());

      ##### Total model for MC
      if options.onlybackgroundfit: ## fit with only mc distribution descrived by fres --> function test
       model_Total_mc    = model_bkg_wjet.Clone("model_Total_background_mc");
      else: 
       model_Total_mc    = RooAddPdf("model_Total_background_mc","model_Total_background_mc",RooArgList(model_total_signal,model_bkg_wjet));

      model_Total_mc.Print();
      getattr(self.workspace4bias_,"import")(model_Total_mc);

      ##### generate models  --> the one fixed and already fitted
      generation_model_wjet = self.workspace4bias_.pdf("model%s%s_%s%s"%(label,options.mlvjregion+options.fgen,self.channel,spectrum));
      generation_model_wjet.Print();

      self.workspace4bias_.Print();             

      ## variable and RooMC study
      numevents_mc   =( self.workspace4bias_.data("rdataset"+label+options.mlvjregion+"_"+self.channel+spectrum).sumEntries()+(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght)*options.inflatejobstatistic;
      print"########  numevents mc ",numevents_mc;

      print "###########################################################";
      print "################ Call the toy class tool ##################";
      print "###########################################################";
      os.system("rm tmp2.root");      
      self.outputFile.cd();
      mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                                           generation_model_wjet,
                                           self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                           int(options.nexp),
                                           int(options.isMC));

      mcWjetTreeResult.setTree(self.outputTree);
      mcWjetTreeResult.setFittingModel(model_Total_mc);
      mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel,label);
      mcWjetTreeResult.setBackgroundPdfCore(model_bkg_wjet);
      mcWjetTreeResult.setSignalInjection(model_total_signal,(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght,options.scalesignalwidth);
      mcWjetTreeResult.generateAndFitToys(int(numevents_mc));
      mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
      mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_,self.mlvj_shape,options.jetBin);
      self.outputTree.Write();

      ratePlotsToStore = 0 ;
      if options.nexp <= 10 :
          ratePlotsToStore = 1 ;
      elif options.nexp > 10 and options.nexp < 50:
          ratePlotsToStore = 2 ;
      elif options.nexp >= 50 and options.nexp < 100:
          ratePlotsToStore = 3 ;
      elif options.nexp >= 100:
          ratePlotsToStore = 10 ;
          
      if(options.storeplot):
          mcWjetTreeResult.saveToysPlots(int(ratePlotsToStore),options.fitjetmass); 

      self.outputFile.Close();

     else: 
         
      ############### Make the Data analysis --> make the Entended pdf for the bkg
      constraintlist_bkg_data = RooArgList();
      print "#################################################################################";
      print "################ Start the MC analysis -> bkg model in the toy ##################";
      print "#################################################################################";

      ### take the models for the background component
      if fitjetmass: 
       model_VV_backgrounds     = get_VV_mj_Model(self.workspace4bias_,"_VV",self.mj_shape["VV"],self.channel)
       model_STop_backgrounds   = get_STop_mj_Model(self.workspace4bias_,"_STop",self.mj_shape["STop"],self.channel);

       if options.jetBin == "_2jet":    
        model_WW_EWK_backgrounds = get_WW_EWK_mj_Model(self.workspace4bias_,"_WW_EWK",self.mj_shape["WW_EWK"],self.channel);
       
       ## inflate the number of events and print them
       self.workspace4bias_.var("rrv_number_VV%s_%s_mj"%(self.mj_shape["VV"],self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV%s_%s_mj"%(self.mj_shape["VV"],self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop%s_%s_mj"%(self.mj_shape["STop"],self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop%s_%s_mj"%(self.mj_shape["STop"],self.channel)).getVal()*options.inflatejobstatistic);
       print "VV number ",self.workspace4bias_.var("rrv_number_VV%s_%s_mj"%(self.mj_shape["VV"],self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "STop number ",self.workspace4bias_.var("rrv_number_STop%s_%s_mj"%(self.mj_shape["STop"],self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.jetBin == "_2jet":    
        self.workspace4bias_.var("rrv_number_WW_EWK%s_%s_mj"%(self.mj_shape["WW_EWK"],self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK%s_%s_mj"%(self.mj_shape["WW_EWK"],self.channel)).getVal()*options.inflatejobstatistic);
        print "WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK%s_%s_mj"%(self.mj_shape["WW_EWK"],self.channel)).getVal()," inflate ",options.inflatejobstatistic;


       if options.ttbarcontrolregion:

        model_TTbar_backgrounds  = get_TTbar_mj_Model(self.workspace4bias_,label,options.fgen,self.channel);
        model_WJets_backgrounds  = get_WJets_mj_Model(self.workspace4bias_,"_WJets0",self.mj_shape["WJets0"],self.channel,1);

        self.workspace4bias_.var("rrv_number_WJets0%s_%s_mj"%(self.mj_shape["WJets0"],self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0%s_%s_mj"%(self.mj_shape["WJets0"],self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);

        print "WJets number ",self.workspace4bias_.var("rrv_number_WJets0%s_%s_mj"%(self.mj_shape["WJets0"],self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       else:

        model_TTbar_backgrounds  = get_TTbar_mj_Model(self.workspace4bias_,"_TTbar",self.mj_shape["TTbar"],self.channel);           
        model_WJets_backgrounds  = get_WJets_mj_Model(self.workspace4bias_,label,options.fgen,self.channel,1,options.jetBin);

        self.workspace4bias_.var("rrv_number_TTbar%s_%s_mj"%(self.mj_shape["TTbar"],self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar%s_%s_mj"%(self.mj_shape["TTbar"],self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic)

        print "WJets number ",self.workspace4bias_.var("rrv_number_TTbar%s_%s_mj"%(self.mj_shape["TTbar"],self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

      else:	       
       ### in case of mWW analysis
       model_VV_backgrounds     = get_VV_mlvj_Model(self.workspace4bias_,"_VV",options.mlvjregion,options.fres,self.channel);
       model_STop_backgrounds   = get_STop_mlvj_Model(self.workspace4bias_,"_STop",options.mlvjregion,options.fres,self.channel);
       model_TTbar_backgrounds  = get_TTbar_mlvj_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,options.fres,self.channel);
       model_WJets_backgrounds  = get_WJets_mlvj_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,options.fres,self.channel);
       if options.jetBin == "_2jet":    
        model_WW_EWK_backgrounds = get_WW_EWK_mlvj_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,options.fres,self.channel);

       ## inflate yields
       self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()*options.inflatejobstatistic);
       print "VV number ",self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "STop number ",self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.jetBin == "_2jet":    
        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()*options.inflatejobstatistic) ## get the normalization
        print "WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()," inflate ",options.inflatejobstatistic;


       if options.ttbarcontrolregion == 0:

        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()*options.inflatejobstatistic)  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print "TTbar number ",self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "WJets number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
       else:

        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()*options.inflatejobstatistic);  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print " WJets number ",self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fres,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " TTbar number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fres+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
           
      #### make the global model for the background  
      if options.fitjetmass:

       model_bkg_data = MakeExtendedModel(self.workspace4bias_,label+signal_region+"_fit",options.fres,spectrum,self.channel,self.wtagger_label,constraintlist_bkg_data,1);
       model_bkg_data.Print();

       if options.fgen == options.fres:
        clone_Model(self.workspace4bias_,model_bkg_data,label,options.mlvjregion,spectrum,options.fgen,self.channel,0);
              
       getattr(self.workspace4bias_,"import")(model_bkg_data);    
       self.workspace4bias_.var("rrv_number"+label+signal_region+"_fit_"+self.channel+spectrum).setVal(self.workspace4bias_.var("rrv_number"+label+signal_region+options.fgen+"_"+self.channel+spectrum).getVal());
       self.workspace4bias_.var("rrv_number"+label+signal_region+"_fit_"+self.channel+spectrum).Print();

      else:
       model_bkg_data = MakeExtendedModel(self.workspace4bias_,label+options.mlvjregion+"_fit",options.fres,spectrum,self.channel,self.wtagger_label,constraintlist_bkg_data,1);
       model_bkg_data.Print();

       if options.fgen == options.fres :
        clone_Model(self.workspace4bias_,model_bkg_data,label,options.mlvjregion,spectrum,options.fgen+"_from_fitting",self.channel,0);
       else:
        clone_Model(self.workspace4bias_,model_bkg_data,label,options.mlvjregion,spectrum,options.fres,self.channel,0);
 
       getattr(self.workspace4bias_,"import")(model_bkg_data);
       self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+"_fit_"+self.channel+spectrum).setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+spectrum).getVal());
       self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+"_fit_"+self.channel+spectrum).Print();

      ## Add the other bkg component fixed to the total model --> in the extended way
      if options.onlybackgroundfit == 1 and options.ttbarcontrolregion == 0:
        if options.jetBin == "_2jet":
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
        else:
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds));
            
      elif options.onlybackgroundfit == 1 and options.ttbarcontrolregion == 1: 
        if options.jetBin == "_2jet":
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
        else:
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds));
            
      elif options.onlybackgroundfit == 0 and options.ttbarcontrolregion == 0:
        if options.jetBin == "_2jet":
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
        else:
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds));
            
      elif options.onlybackgroundfit == 0 and options.ttbarcontrolregion == 1:  
        if options.jetBin == "_2jet":
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
        else:
         model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds));
                                                                                                       
      model_Total_data.Print();
      getattr(self.workspace4bias_,"import")(model_Total_data);

      ##### generate models  --> the one fixed and already fitted
      if options.fitjetmass:
       generation_model_data = self.workspace4bias_.pdf("model_data_%s_mj"%(self.channel));
       generation_model_data.Print();
      else:    
       generation_model_data = self.workspace4bias_.pdf("model_%s%s%s_%s%s"%("data",label,options.mlvjregion+options.fgen,self.channel,spectrum));
       generation_model_data.Print();
       
      self.workspace4bias_.Print();             
      
      numevents_data   = (self.workspace4bias_.data("rdataset_data"+options.mlvjregion+"_"+self.channel+spectrum).sumEntries()+(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght)*options.inflatejobstatistic;
      print "##### number of events generated ",numevents_data ;
      
      print "###########################################################";
      print "################ Call the toy class tool ##################";
      print "###########################################################";

      if options.fitjetmass :
       self.outputFile.cd();
       mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_j")),
                                            generation_model_data,
                                            self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                            int(options.nexp),
                                            int(options.isMC));

       mcWjetTreeResult.setTree(self.outputTree);
       mcWjetTreeResult.setFittingModel(model_Total_data);
       mcWjetTreeResult.setBackgroundPdfCore(model_bkg_data);
       mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel,label);
       mcWjetTreeResult.setSignalInjection(model_total_signal,(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght,options.scalesignalwidth);
       mcWjetTreeResult.generateAndFitToys(int(numevents_data),"sb_lo,sb_hi");
       self.outputFile.cd();
       mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
       mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_,self.mj_shape,options.jetBin);
      else:
       self.outputFile.cd();
       mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                                            generation_model_data,
                                            self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                            int(options.nexp),
                                            int(options.isMC));
       mcWjetTreeResult.setTree(self.outputTree);
       mcWjetTreeResult.setFittingModel(model_Total_data);
       mcWjetTreeResult.setBackgroundPdfCore(model_bkg_data);
       mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel,label);
       mcWjetTreeResult.setSignalInjection(model_total_signal,(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght,options.scalesignalwidth);
       mcWjetTreeResult.generateAndFitToys(int(numevents_data));
       self.outputFile.cd();
       mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
       mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_,self.mlvj_shape,options.jetBin);


      ratePlotsToStore = 0 ;
      if options.nexp <= 10 :
          ratePlotsToStore = 1 ;
      elif options.nexp > 10 and options.nexp < 50:
          ratePlotsToStore = 2 ;
      elif options.nexp >= 50 and options.nexp < 100:
          ratePlotsToStore = 3 ;
      elif options.nexp >= 100:
          ratePlotsToStore = 10 ;
          
      if(options.storeplot):
       mcWjetTreeResult.saveToysPlots(int(ratePlotsToStore),options.fitjetmass); 

      self.outputTree.Write();
      self.outputFile.Close();

#### Main code     
if __name__ == "__main__":

  print "###################### begin the analysis: channel %s, signal name %s, mlvj min %s, mlvj max %s, mj min %s, mj max %s, genfunction %s, fitfunction %s"%(options.channel,sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],options.fgen,options.fres);

  fitBiasAnalysis = doBiasStudy_mlvj (options.channel,sys.argv[1],int(sys.argv[2]),int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),options.fgen,options.fres)

  if options.shapetest == 0 and options.ttbarcontrolregion == 0:
     fitBiasAnalysis.biasAnalysis("_WJets0",options.fitjetmass);
     fitBiasAnalysis.outputFile.Close();                                  
  elif options.shapetest == 1 and options.ttbarcontrolregion == 0: 
     fitBiasAnalysis.shapeParametrizationAnalysis("_WJets0");
  elif options.shapetest == 0 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.biasAnalysis("_TTbar",options.fitjetmass);
     fitBiasAnalysis.outputFile.Close();                                  
  elif options.shapetest == 1 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.shapeParametrizationAnalysis("_TTbar");                     
