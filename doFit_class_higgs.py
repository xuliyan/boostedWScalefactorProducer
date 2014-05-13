#! /Usr/bin/env python
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

from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG,RooDataSet,RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, RooChi2Var, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite, TGraph, RooMCStudy

############################################
#              Job steering                #
############################################

parser = OptionParser()

##### basic options ###########

parser.add_option('-b', '--noPlots',action='store_true', dest='noX',    default=False, help='no X11 windows')
parser.add_option('--check',  action='store_true', dest='check',  default=False, help='check the workspace for limit setting')
parser.add_option('-s', '--simple', action='store_true', dest='simple', default=False, help='pre-limit in simple mode')
parser.add_option('-m', '--multi',  action='store_true', dest='multi',  default=True,  help='pre-limit in multi mode')


#### additional information: channel. jet bin, signal properties

parser.add_option('-p', '--psmodel',             action="store", type="string", dest="psmodel",             default="pythia")
parser.add_option('-a', '--additioninformation', action="store", type="string", dest="additioninformation", default="HIGGS")
parser.add_option('-c', '--channel',             action="store", type="string", dest="channel",             default="mu")
parser.add_option('-i', '--inPath',              action="store", type="string", dest="inPath",              default="./")

parser.add_option('--cprime',      action="store", type="int",    dest="cprime",      default=10)
parser.add_option('--BRnew',       action="store", type="int",    dest="BRnew",       default=0)
parser.add_option('--closuretest', action='store', type="int",    dest='closuretest', default=0, help='closure test; 0: no test; 1: A1->A2; 2: A->B')
parser.add_option('--pseudodata',  action='store', type="int",    dest='pseudodata',  default=1, help='pseudodata 0 -> use real data, else use stack of MC backgrounds')
parser.add_option('--fitSignal',   action='store', type="int",    dest='fitsignal',   default=0, help='fit only signal lineshape with a chosen model')
parser.add_option('--category',    action="store", type="string", dest="category",    default="HP")
parser.add_option('--jetBin',      action="store", type="string", dest="jetBin",      default="")
parser.add_option('--skipJetSystematics',      action="store", type="int", dest="skipJetSystematics",  default=0)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/MakePdf_cxx.so")
ROOT.gSystem.Load(options.inPath+"/BiasStudy/BiasUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/FitUtils/FitUtils_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooAlpha4ExpTailPdf, Roo2ExpPdf, RooAlpha42ExpPdf

from ROOT import MakeGeneralPdf, MakeExtendedModel, get_TTbar_mj_Model, get_STop_mj_Model, get_VV_mj_Model, get_WW_EWK_mj_Model, get_WJets_mj_Model, get_ggH_mj_Model, get_vbfH_mj_Model, get_TTbar_mlvj_Model, get_STop_mlvj_Model, get_VV_mlvj_Model, get_WW_EWK_mlvj_Model, get_WJets_mlvj_Model, get_ggH_mlvj_Model, get_vbfH_mlvj_Model, fix_Model,  clone_Model

from ROOT import setTDRStyle, get_pull, draw_canvas, draw_canvas_with_pull, legend4Plot, GetDataPoissonInterval, GetLumi, draw_error_band_ws

from ROOT import fit_mj_single_MC, fit_mlvj_model_single_MC, fit_WJetsNormalization_in_Mj_signal_region, fit_mlvj_in_Mj_sideband, get_WJets_mlvj_correction_sb_lo_to_signal_region, get_mlvj_normalization_insignalregion, fit_genHMass, SystematicUncertaintyHiggs_2jetBin, SystematicUncertaintyHiggs_01jetBin

from ROOT import *

gInterpreter.GenerateDictionary("std::map<std::string,std::string>", "map;string;string")

###############################
## doFit Class Implemetation ##
###############################

class doFit_wj_and_wlvj:

    def __init__(self, in_channel,in_higgs_sample, in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700, in_mj_min=30, in_mj_max=140, in_mlvj_min=400., in_mlvj_max=1400., fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1", input_workspace=None):

                   
        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ### shapes to be used in mj
        self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string) () ;
        self.mj_shape["TTbar"]   = "2Gaus_ErfExp";
        self.mj_shape["VV"]      = "2_2Gaus";
        self.mj_shape["WW_EWK"]  = "2_2Gaus";
        if options.jetBin == "_2jet":  self.mj_shape["STop"]    = "ErfExp";
        else: self.mj_shape["STop"]    = "ErfExpGaus_sp"; 
        if options.jetBin == "_2jet":
            self.mj_shape["WJets0"]  = "ErfExp";
            self.mj_shape["WJets1"]  = "ErfExp";
            self.mj_shape["WJets01"] = "User1";
        else: 
            self.mj_shape["WJets0"]  = "User1";
            self.mj_shape["WJets1"]  = "User1";
            self.mj_shape["WJets01"] = "ErfExp";

        self.mlvj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string) () ;
        self.mlvj_shape["TTbar"]   = fit_model;
        self.mlvj_shape["VV"]      = fit_model;
        self.mlvj_shape["WW_EWK"]  = fit_model;
        self.mlvj_shape["STop"]    = fit_model;
        self.mlvj_shape["WJets0"]  = fit_model;
        self.mlvj_shape["WJets1"]  = fit_model;
        self.mlvj_shape["WJets01"] = fit_model_alter;
        self.mlvj_shape["ggH"]     = "CB_v1";
        self.mlvj_shape["vbfH"]    = "CB_v1";

        self.tmpFile = TFile("tmp2.root","RECREATE");
        self.tmpFile.cd();
        ### set the channel type --> electron or muon
        self.channel = in_channel;
        self.higgs_sample = in_higgs_sample;

        if in_higgs_sample == "ggH600":  self.vbfhiggs_sample = "vbfH600";
        if in_higgs_sample == "ggH700":  self.vbfhiggs_sample = "vbfH700";
        if in_higgs_sample == "ggH800":  self.vbfhiggs_sample = "vbfH800";
        if in_higgs_sample == "ggH900":  self.vbfhiggs_sample = "vbfH900";
        if in_higgs_sample == "ggH1000": self.vbfhiggs_sample = "vbfH1000";
        if in_higgs_sample == "ggH1500": self.vbfhiggs_sample = "vbfH1500";
        if in_higgs_sample == "ggH2000": self.vbfhiggs_sample = "vbfH2000";

        print "########################################################################################"
        print "######## define class: binning, variables, cuts, files and nuissance parameters ########"
        print "########################################################################################"

        ### Set the mj binning for plots
        self.BinWidth_mj = 5.;

        ## set the model used for the background parametrization
        self.MODEL_4_mlvj = fit_model;
        self.MODEL_4_mlvj_alter = fit_model_alter;

        ### Set the binning for mlvj plots as a function of the model
        if not options.fitsignal:
         if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1":
            self.BinWidth_mlvj = 35.;
         else:
            self.BinWidth_mlvj = 50.;
        else:
         if self.MODEL_4_mlvj == "ErfPowExp_v1" or self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfExp_v1":
            self.BinWidth_mlvj = 10.;
         else:
            self.BinWidth_mlvj = 10.;

        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy solution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW

        self.leg = TLegend();        
        self.narrow_factor = 10;

        ## correct the binning of mj
        self.BinWidth_mj = self.BinWidth_mj;
        self.nbins_mj    = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max        = in_mj_min+self.nbins_mj*self.BinWidth_mj;

        ## correct the binning of mlvj
        self.BinWidth_mlvj = self.BinWidth_mlvj;
        self.nbins_mlvj    = int((in_mlvj_max-in_mlvj_min)/self.BinWidth_mlvj);
        in_mlvj_max        = in_mlvj_min+self.nbins_mlvj*self.BinWidth_mlvj;

        ## define jet mass variable
        rrv_mass_j = RooRealVar("rrv_mass_j","pruned m_{J}",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV/c^{2}");
        rrv_mass_j.setBins(self.nbins_mj);

        ## define invariant mass WW variable
        rrv_mass_lvj = RooRealVar("rrv_mass_lvj","m_{WW}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV/c^{2}");
        rrv_mass_lvj.setBins(self.nbins_mlvj);

        ## generator higgs mass
        rrv_mass_gen_WW = RooRealVar("rrv_mass_gen_WW","gen_m_{WW}",(in_mlvj_min+in_mlvj_max)/2.,in_mlvj_min,in_mlvj_max,"GeV/c^{2}");
        rrv_mass_gen_WW.setBins(self.nbins_mlvj);

        ## create the workspace and import them
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_");
        else:
            self.workspace4fit_ = input_workspace;

        getattr(self.workspace4fit_,"import")(rrv_mass_j);
        getattr(self.workspace4fit_,"import")(rrv_mass_lvj);
        getattr(self.workspace4fit_,"import")(rrv_mass_gen_WW);

        #prepare workspace for unbin-Limit -> just fo the stuff on which running the limit
        self.workspace4limit_ = RooWorkspace("workspace4limit_","workspace4limit_");

        ## different code operation mode -> just normal analysis
	if options.closuretest == 0:
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 65;
            self.mj_signal_min = 65;
            self.mj_signal_max = 105;
            self.mj_sideband_hi_min = 105;
            self.mj_sideband_hi_max = in_mj_max;
        if options.closuretest == 1: ##closure test A1->A2
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 55;
            self.mj_signal_min = 55;
            self.mj_signal_max = 65;
            self.mj_sideband_hi_min = 105;
            self.mj_sideband_hi_max = in_mj_max;
        if options.closuretest == 2: #closure test A->B
            self.mj_sideband_lo_min = in_mj_min;
            self.mj_sideband_lo_max = 65;
            self.mj_signal_min = 100;
            self.mj_signal_max = 115;
            self.mj_sideband_hi_min = 115;
            self.mj_sideband_hi_max = in_mj_max;

        ## zone definition in the jet mass
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("sblo_to_sbhi",self.mj_sideband_lo_min,self.mj_sideband_hi_max);

        ## signal region definition in the mlvj variable in case of counting limit
        self.mlvj_signal_min = in_mlvj_signal_region_min
        self.mlvj_signal_max = in_mlvj_signal_region_max
        rrv_mass_lvj.setRange("signal_region",self.mlvj_signal_min,self.mlvj_signal_max);

        #prepare the data and mc files --> set the working directory and the files name
        self.file_Directory = options.inPath+"./trainingtrees_%s/"%(self.channel);

        self.PS_model = options.psmodel;

        if options.pseudodata == 1:
            self.file_data = ("ofile_pseudodata.root");
        else:
            self.file_data  = ("ofile_data.root");

        self.file_ggH   = ("ofile_%s.root"%(self.higgs_sample));
        self.file_vbfH  = ("ofile_%s.root"%(self.vbfhiggs_sample));

        #WJets0 is the default PS model, WJets1 is the alternative PS model
        if self.PS_model == "pythia":
            if options.jetBin == "_2jet" :
             self.file_WJets0_mc = ("ofile_WJets_exclusive_Pythia.root");
             self.file_WJets1_mc = ("ofile_WJets_Herwig.root");
            else: 
             self.file_WJets1_mc = ("ofile_WJets_Pythia100.root");
             self.file_WJets0_mc = ("ofile_WJets_Herwig.root");
        else:
            if options.jetBin == "_2jet" :
             self.file_WJets0_mc = ("ofile_WJets_Herwig.root");
             self.file_WJets1_mc = ("ofile_WJets_exclusive_Pythia.root");
            else:
             self.file_WJets0_mc = ("ofile_WJets_Herwig.root");
             self.file_WJets1_mc = ("ofile_WJets_Pythia100.root");
                

        self.file_VV_mc     = ("ofile_VV.root");# WW+WZ
        self.file_WW_EWK_mc = ("ofile_WW2jet_phantom.root");# WW_EWK       
        self.file_STop_mc   = ("ofile_STop.root");#STop
        self.file_TTbar_mc  = ("ofile_TTbar_Powheg.root");
        self.file_TTbar_matchDn_mc = ("ofile_TTbar_matchDn.root");
        self.file_TTbar_matchUp_mc = ("ofile_TTbar_matchUp.root");
        self.file_TTbar_scaleDn_mc = ("ofile_TTbar_scaleDn.root");
        self.file_TTbar_scaleUp_mc = ("ofile_TTbar_scaleUp.root");
        self.file_TTbar_mcanlo_mc = ("ofile_TTbar_mcanlo.root");
                                                                       
        self.PS_model= options.psmodel
 
        ## event categorization as a function of the purity and the applied selection
        self.wtagger_label = options.category;

        if self.wtagger_label=="HP" :
            if self.channel=="el":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="mu":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
            if self.channel=="em":
                self.wtagger_cut=0.5 ; self.wtagger_cut_min=0. ;
        if self.wtagger_label=="LP":
            self.wtagger_cut=0.75 ;
            self.wtagger_cut_min=0.5 ;

        if self.wtagger_label=="nocut":
            self.wtagger_cut=10000;

        #medium wtagger_eff reweight between data and mc #Wtagger_forV SF have be add to ntuple weight;
	if self.channel=="mu" and self.wtagger_label=="HP":
          if options.pseudodata == 1:
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
           self.rrv_wtagger_eff_reweight_forT.setError(0.06);
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
          elif options.pseudodata == 0 and  options.jetBin == "_2jet":
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.128);
           self.rrv_wtagger_eff_reweight_forT.setError(0.338);
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
          elif options.pseudodata == 0 and not options.jetBin == "_2jet":
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.96);
           self.rrv_wtagger_eff_reweight_forT.setError(0.06*self.rrv_wtagger_eff_reweight_forT.getVal());
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
                                                                                                                          
        if self.channel=="el" and self.wtagger_label=="HP":
          if options.pseudodata == 1:
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",1.);
           self.rrv_wtagger_eff_reweight_forT.setError(0.08);
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
          elif options.pseudodata == 0 and options.jetBin == "_2jet":
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.96);
           self.rrv_wtagger_eff_reweight_forT.setError(0.369);
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
          elif options.pseudodata == 0 and not  options.jetBin == "_2jet":
           self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT",0.89);
           self.rrv_wtagger_eff_reweight_forT.setError(0.06*self.rrv_wtagger_eff_reweight_forT.getVal());
           self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.87);
           self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
                                                  
        if self.channel=="em" and self.wtagger_label=="HP":
          if options.pseudodata == 1:
            self.rrv_wtagger_eff_reweight_forT=RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.0);
            self.rrv_wtagger_eff_reweight_forT.setError(0.265);
            self.rrv_wtagger_eff_reweight_forV=RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",1.0);
            self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());
          elif options.pseudodata ==0 and options.jetBin == "_2jet":
            self.rrv_wtagger_eff_reweight_forT = RooRealVar("rrv_wtagger_eff_reweight_forT","rrv_wtagger_eff_reweight_forT", 1.);
            self.rrv_wtagger_eff_reweight_forT.setError(0.265);
            self.rrv_wtagger_eff_reweight_forV = RooRealVar("rrv_wtagger_eff_reweight_forV","rrv_wtagger_eff_reweight_forV",0.93);
            self.rrv_wtagger_eff_reweight_forV.setError(0.097*self.rrv_wtagger_eff_reweight_forV.getVal());              


        print "wtagger efficiency correction for Top sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forT.getVal(), self.rrv_wtagger_eff_reweight_forT.getError());
        print "wtagger efficiency correction for V sample: %s +/- %s"%(self.rrv_wtagger_eff_reweight_forV.getVal(), self.rrv_wtagger_eff_reweight_forV.getError());


        #result files: The event number, parameters and error write into a txt file. The dataset and pdfs write into a root file
        if not os.path.isdir("cards_%s_%s"%(self.channel,self.mlvj_shape["WJets0"])): os.system("mkdir cards_%s_%s"%(self.channel,self.mlvj_shape["WJets0"]));
        self.rlt_DIR = "cards_%s_%s/"%(self.channel,self.mlvj_shape["WJets0"]);

        if options.jetBin == "_2jet" : 
         self.file_rlt_txt                   = self.rlt_DIR+"other_hwwlvj_%s_%s%s_%02d_%02d.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_rlt_root                  = self.rlt_DIR+"hwwlvj_%s_%s%s_%02d_%02d_workspace.root"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_unbin_ggHvbfH    = self.rlt_DIR+"hwwlvj_%s_%s%s_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_unbin_ggH        = self.rlt_DIR+"hwwlvj_%s_%s%s_ggH_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_unbin_vbfH       = self.rlt_DIR+"hwwlvj_%s_%s%s_vbfH_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_counting_ggHvbfH = self.rlt_DIR+"hwwlvj_%s_%s%s_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_counting_ggH     = self.rlt_DIR+"hwwlvj_%s_%s%s_ggH_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
         self.file_datacard_counting_vbfH    = self.rlt_DIR+"hwwlvj_%s_%s%s_vbfH_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.jetBin,options.cprime,options.BRnew)
        else:
         self.file_rlt_txt                   = self.rlt_DIR+"other_hwwlvj_%s_%s_%02d_%02d.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_rlt_root                  = self.rlt_DIR+"hwwlvj_%s_%s_%02d_%02d_workspace.root"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_unbin_ggHvbfH    = self.rlt_DIR+"hwwlvj_%s_%s_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_unbin_ggH        = self.rlt_DIR+"hwwlvj_%s_%s_ggH_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_unbin_vbfH       = self.rlt_DIR+"hwwlvj_%s_%s_vbfH_%02d_%02d_unbin.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_counting_ggHvbfH = self.rlt_DIR+"hwwlvj_%s_%s_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_counting_ggH     = self.rlt_DIR+"hwwlvj_%s_%s_ggH_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
         self.file_datacard_counting_vbfH    = self.rlt_DIR+"hwwlvj_%s_%s_vbfH_%02d_%02d_counting.txt"%(self.higgs_sample,self.channel,options.cprime,options.BRnew)
            
        self.file_out = open(self.file_rlt_txt,"w");
        self.file_out.write("Welcome:\n");
        self.file_out.close()
        self.file_out = open(self.file_rlt_txt,"a+");

        self.higgs_xs_scale=1.0; #higgs XS scale

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
                                                                        
        
        self.Lumi = 19297;
        if self.channel=="el":
            self.Lumi = 19166;

	#met cut:el 70; mu: 50
        self.pfMET_cut = 50;
        self.lpt_cut   = 30;
        self.vpt_cut   = 200;
        self.bcut      = 0.679;
        if self.channel=="el":
            self.pfMET_cut = 50; 
            self.lpt_cut   = 35;        
        #deltaPhi_METj cut
        self.deltaPhi_METj_cut = 2.0;
        self.top_veto_had = 200 ;
        self.top_veto_lep = 200 ;
        self.top_veto_had_min = 0 ;
        self.top_veto_lep_min = 0 ;
        self.dEta_cut = 2.5 ;
        self.Mjj_cut  = 250 ;

        # parameters of data-driven method to get the WJets background event number.
        self.number_WJets_insideband = -1;
        self.datadriven_alpha_WJets_unbin = -1;
        self.datadriven_alpha_WJets_counting = -1;

        #uncertainty for datacard
        self.lumi_uncertainty        = 0.026;
        self.XS_STop_uncertainty     = 0.30 ;
        self.XS_VV_uncertainty       = 0.20 ;
        self.XS_WW_EWK_uncertainty   = 0.20 ;        
        self.XS_TTbar_uncertainty    = 0.07 ;
        self.XS_TTbar_NLO_uncertainty = 0.063 ;# from AN-12-368 table8
        self.XS_STop_NLO_uncertainty  = 0.05 ; # from AN-12-368 table8
        self.XS_VV_NLO_uncertainty    = 0.10 ; # from AN-12-368 table8

        ### jet binning uncertainty 
        if options.jetBin == "_2jet":                                                
         self.QCDscale_ggH0in   = 0.000;
         self.QCDscale_ggH2in   = 0.190;
        else:
         self.QCDscale_ggH0in   = 0.260;
         self.QCDscale_ggH2in   = -0.060;
            
        self.QCDscale_qqH      = 0.0;
        self.pdf_gg            = 0.0;
        self.pdf_qqbar         = 0.0;
        self.QCDscale_ggH_ACCEPT = 0.0;
        self.QCDscale_qqH_ACCEPT = 0.0;

        # from twiki https:#twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt8TeV,  
        if self.higgs_sample == "ggH600": 
            self.QCDscale_qqH         = 0.007 ;
            self.pdf_gg               = 0.091 ;
            self.pdf_qqbar            = 0.050 ;
            self.QCDscale_ggH_ACCEPT  = 0.036 ;
            self.QCDscale_qqH_ACCEPT  = 0.007 ;

        elif self.higgs_sample == "ggH700": 
            self.QCDscale_qqH         = 0.008;
            self.pdf_gg               = 0.101;
            self.pdf_qqbar            = 0.042
            self.QCDscale_ggH_ACCEPT  = 0.038;
            self.QCDscale_qqH_ACCEPT  = 0.008

        elif self.higgs_sample == "ggH800": 
            self.QCDscale_qqH        = 0.010;
            self.pdf_gg              = 0.106;
            self.pdf_qqbar           = 0.047;
            self.QCDscale_ggH_ACCEPT = 0.040;
            self.QCDscale_qqH_ACCEPT = 0.009

        elif self.higgs_sample == "ggH900": 
            self.QCDscale_qqH        = 0.012
            self.pdf_gg              = 0.111;
            self.pdf_qqbar           = 0.053
            self.QCDscale_ggH_ACCEPT = 0.042;
            self.QCDscale_qqH_ACCEPT = 0.010

        elif self.higgs_sample == "ggH1000": 
            self.QCDscale_qqH        = 0.013;
            self.pdf_gg              = 0.121;
            self.pdf_qqbar           = 0.059; 
            self.QCDscale_ggH_ACCEPT = 0.046;
            self.QCDscale_qqH_ACCEPT = 0.011

        ### interference effect 
        self.interference_ggH_uncertainty  = 0.10;
        self.interference_vbfH_uncertainty = 0.25;

        #normalization uncertainty from jet scale
        self.WJets_normalization_uncertainty_from_jet_scale  = 0.;        
        self.VV_normalization_uncertainty_from_jet_scale     = 0.;
        self.WW_EWK_normalization_uncertainty_from_jet_scale = 0.;        
        self.STop_normalization_uncertainty_from_jet_scale   = 0.;
        self.TTbar_normalization_uncertainty_from_jet_scale  = 0.;
        self.ggH_normalization_uncertainty_from_jet_scale    = 0.;
        self.vbf_normalization_uncertainty_from_jet_scale    = 0.;

        #normalization uncertainty from jet_res
        self.WJets_normalization_uncertainty_from_jet_res  = 0.;        
        self.VV_normalization_uncertainty_from_jet_res     = 0.;
        self.WW_EWK_normalization_uncertainty_from_jet_res = 0.;        
        self.STop_normalization_uncertainty_from_jet_res   = 0.;
        self.TTbar_normalization_uncertainty_from_jet_res  = 0.;
        self.ggH_normalization_uncertainty_from_jet_res    = 0.;
        self.vbf_normalization_uncertainty_from_jet_res    = 0.;        

        #normalization uncertainty from lep scale
        if self.channel == "mu":
         self.WJets_normalization_uncertainty_from_lep_scale  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_scale     = 1.083;
         self.WW_EWK_normalization_uncertainty_from_lep_scale = 1.008;        
         self.STop_normalization_uncertainty_from_lep_scale   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_scale  = 1.008;
         self.ggH_normalization_uncertainty_from_lep_scale    = 1.028;
         self.vbf_normalization_uncertainty_from_lep_scale    = 1.015;
        elif self.channel == "el":
         self.WJets_normalization_uncertainty_from_lep_scale  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_scale     = 1.068;
         self.WW_EWK_normalization_uncertainty_from_lep_scale = 1.006;        
         self.STop_normalization_uncertainty_from_lep_scale   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_scale  = 1.000;
         self.ggH_normalization_uncertainty_from_lep_scale    = 1.014;
         self.vbf_normalization_uncertainty_from_lep_scale    = 1.004;
        elif self.channel == "em":
         self.WJets_normalization_uncertainty_from_lep_scale  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_scale     = 1.075;
         self.WW_EWK_normalization_uncertainty_from_lep_scale = 1.007;        
         self.STop_normalization_uncertainty_from_lep_scale   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_scale  = 1.004;
         self.ggH_normalization_uncertainty_from_lep_scale    = 1.021;
         self.vbf_normalization_uncertainty_from_lep_scale    = 1.010;
            

        #normalization uncertainty from lep_res
        if self.channel == "mu":
         self.WJets_normalization_uncertainty_from_lep_res  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_res     = 1.016;
         self.WW_EWK_normalization_uncertainty_from_lep_res = 1.000;        
         self.STop_normalization_uncertainty_from_lep_res   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_res  = 1.000;
         self.ggH_normalization_uncertainty_from_lep_res    = 1.001;
         self.vbf_normalization_uncertainty_from_lep_res    = 1.000;        
        elif self.channel == "el":
         self.WJets_normalization_uncertainty_from_lep_res  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_res     = 1.000;
         self.WW_EWK_normalization_uncertainty_from_lep_res = 1.000;        
         self.STop_normalization_uncertainty_from_lep_res   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_res  = 1.000;
         self.ggH_normalization_uncertainty_from_lep_res    = 1.015;
         self.vbf_normalization_uncertainty_from_lep_res    = 1.001;        
        elif self.channel == "em":
         self.WJets_normalization_uncertainty_from_lep_res  = 1.000;        
         self.VV_normalization_uncertainty_from_lep_res     = 1.008;
         self.WW_EWK_normalization_uncertainty_from_lep_res = 1.000;        
         self.STop_normalization_uncertainty_from_lep_res   = 1.000;
         self.TTbar_normalization_uncertainty_from_lep_res  = 1.000;
         self.ggH_normalization_uncertainty_from_lep_res    = 1.008;
         self.vbf_normalization_uncertainty_from_lep_res    = 1.001;        
           
        #normalization uncertainty from btag
        self.WJets_normalization_uncertainty_from_btag  = 1.000;        
        self.VV_normalization_uncertainty_from_btag     = 1.006;
        self.WW_EWK_normalization_uncertainty_from_btag = 1.007;        
        self.STop_normalization_uncertainty_from_btag   = 1.033;
        self.TTbar_normalization_uncertainty_from_btag  = 1.017;
        self.ggH_normalization_uncertainty_from_btag    = 1.005;
        self.vbf_normalization_uncertainty_from_btag    = 1.002;        

        #el and mu trigger and eff uncertainty, AN2012_368_v5 12.3
        self.lep_trigger_uncertainty = 0.01;
        self.lep_eff_uncertainty     = 0.02;

        #### increase shape uncertainty
        self.shape_para_error_WJets0 = 2.0;
        self.shape_para_error_TTbar  = 2.0;

        if self.higgs_sample == "ggH600" or self.higgs_sample == "ggH700":
            self.shape_para_error_alpha = 2.0;
        else:
            self.shape_para_error_alpha = 2.0;
        
        # shape parameter uncertainty
        self.FloatingParams = RooArgList("floatpara_list");
        self.FloatingParams_wjet = RooArgList("floatpara_list_wjet");

        ### set the TDR Style
        setTDRStyle();


    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, scaleJetMass = 0.): # to get the normalization of WJets in signal_region

        print "############### Fit mj Normalization ##################"

        ## fit the two version of pdf for Wjets shape if available
        fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0","",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets01","",self.mj_shape["WJets01"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.jetBin == "_2jet": 
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets1","",self.mj_shape["WJets1"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
          
        rrv_WJets0  = self.workspace4fit_.var("rrv_number_WJets0_in_mj_signal_region_from_fitting_%s"%(self.channel)); ## nominal parametrization for Wjets
        rrv_WJets01 = self.workspace4fit_.var("rrv_number_WJets01_in_mj_signal_region_from_fitting_%s"%(self.channel)); ## alternate descrption
        if not options.jetBin == "_2jet": rrv_WJets1 = self.workspace4fit_.var("rrv_number_WJets1_in_mj_signal_region_from_fitting_%s"%(self.channel));
                
        rrv_WJets0.Print();
        rrv_WJets01.Print();
        if not options.jetBin == "_2jet":  rrv_WJets1.Print();
        
        if options.jetBin == "_2jet":
            total_uncertainty = TMath.Sqrt(TMath.Power(rrv_WJets0.getError(),2)+TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2)); ## add in quadrature the difference
        else:
            total_uncertainty = TMath.Sqrt(TMath.Power(rrv_WJets0.getError(),2)+TMath.Power(rrv_WJets01.getVal()-rrv_WJets0.getVal(),2)+TMath.Power(rrv_WJets1.getVal()-rrv_WJets0.getVal(),2)); ## add in quadrature the difference
            
        rrv_WJets0.setError(total_uncertainty);
        rrv_WJets0.Print();
        print "Total Uncertainty in WJtes0 due to fit and shape: uncertainty ",total_uncertainty/rrv_WJets0.getVal();

        if scaleJetMass :


         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0massvbf_jes_up","massvbf_jes_up",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0massvbf_jes_dn","massvbf_jes_dn",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0massvbf_jer","massvbf_jer",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0massvbf_jer_up","massvbf_jer_up",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_WJetsNormalization_in_Mj_signal_region(self.workspace4fit_,self.color_palet,self.mj_shape,"_WJets0massvbf_jer_dn","massvbf_jer_dn",self.mj_shape["WJets0"],self.channel,self.wtagger_label,0,options.pseudodata,self.mj_signal_min,self.mj_signal_max,options.jetBin); ## fit jet mass distribution
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         
         rrv_WJetsmassvbf_jes_up = self.workspace4fit_.var("rrv_number_WJets0massvbf_jes_up_in_mj_signal_region_from_fitting_%s"%(self.channel));
         rrv_WJetsmassvbf_jes_dn = self.workspace4fit_.var("rrv_number_WJets0massvbf_jes_dn_in_mj_signal_region_from_fitting_%s"%(self.channel));        
         rrv_WJetsmassvbf_jer    = self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_in_mj_signal_region_from_fitting_%s"%(self.channel));
         rrv_WJetsmassvbf_jer_up = self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_up_in_mj_signal_region_from_fitting_%s"%(self.channel));        
         rrv_WJetsmassvbf_jer_dn = self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_dn_in_mj_signal_region_from_fitting_%s"%(self.channel));        

         print "######################### wjets scale and resolution effect " ;
         rrv_WJetsmassvbf_jes_up.Print();
         rrv_WJetsmassvbf_jes_dn.Print();        
         rrv_WJetsmassvbf_jer.Print();        
         rrv_WJetsmassvbf_jer_up.Print();        
         rrv_WJetsmassvbf_jer_dn.Print();     

         #jet mass uncertainty on WJets normalization
         if(self.workspace4fit_.var("rrv_number_WJets0massvbf_jes_up_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0massvbf_jes_dn_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_up_in_mj_signal_region_from_fitting_%s"%(self.channel)) and self.workspace4fit_.var("rrv_number_WJets0massvbf_jer_dn_in_mj_signal_region_from_fitting_%s"%(self.channel))):    

            self.WJets_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_WJetsmassvbf_jes_up.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJetsmassvbf_jes_dn.getVal()-rrv_WJets0.getVal() ) )/2.)/rrv_WJets0.getVal();         
            print "Total Uncertainty on WJtes0 due to jes: relaxed uncertainty ",self.WJets_normalization_uncertainty_from_jet_scale;
            self.WJets_normalization_uncertainty_from_jet_res   = ((TMath.Abs(rrv_WJetsmassvbf_jer.getVal()-rrv_WJets0.getVal())+TMath.Abs(rrv_WJetsmassvbf_jer_up.getVal()-rrv_WJets0.getVal() )+TMath.Abs(rrv_WJetsmassvbf_jer_dn.getVal()-rrv_WJets0.getVal() ) )/3.)/rrv_WJets0.getVal();         
            print "Total Uncertainty on WJtes0 due to jes: relaxed uncertainty ",self.WJets_normalization_uncertainty_from_jet_res;

         #jet mass uncertainty on sTop normalization
         rrv_STop                = self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_%s_mj"%(self.channel))
         rrv_STopmassvbf_jes_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jes_up_%s_mj"%(self.channel))
         rrv_STopmassvbf_jes_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jes_dn_%s_mj"%(self.channel))
         rrv_STopmassvbf_jer     = self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_%s_mj"%(self.channel))
         rrv_STopmassvbf_jer_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_up_%s_mj"%(self.channel))
         rrv_STopmassvbf_jer_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_dn_%s_mj"%(self.channel))
         
         rrv_STop.Print();
         rrv_STopmassvbf_jes_up.Print();
         rrv_STopmassvbf_jes_dn.Print();        
         rrv_STopmassvbf_jer.Print();        
         rrv_STopmassvbf_jer_up.Print();        
         rrv_STopmassvbf_jer_dn.Print();        

         #jet mass uncertainty on STop normalization
         if(self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jes_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jes_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_STopmassvbf_jer_%s_mj"%(self.channel))):    

            self.STop_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_STopmassvbf_jes_up.getVal()-rrv_STop.getVal())+TMath.Abs(rrv_STopmassvbf_jes_dn.getVal()-rrv_STop.getVal() ) )/2.)/rrv_STop.getVal();         
            print "Total Uncertainty on STop due to jes: uncertainty ",self.STop_normalization_uncertainty_from_jet_scale;
            self.STop_normalization_uncertainty_from_jet_res   = ((TMath.Abs(rrv_STopmassvbf_jer.getVal()-rrv_STop.getVal())+TMath.Abs(rrv_STopmassvbf_jer_up.getVal()-rrv_STop.getVal() )+TMath.Abs(rrv_STopmassvbf_jer_dn.getVal()-rrv_STop.getVal() ) )/3.)/rrv_STop.getVal();         
            print "Total Uncertainty on STop due to jer: uncertainty ",self.STop_normalization_uncertainty_from_jet_res;

         #jet mass uncertainty on TTbar normalization
         rrv_TTbar                 = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_%s_mj"%(self.channel))
         rrv_TTbarmassvbf_jes_up   = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jes_up_%s_mj"%(self.channel))
         rrv_TTbarmassvbf_jes_dn   = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jes_dn_%s_mj"%(self.channel))
         rrv_TTbarmassvbf_jer      = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_%s_mj"%(self.channel))
         rrv_TTbarmassvbf_jer_dn   = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_up_%s_mj"%(self.channel))
         rrv_TTbarmassvbf_jer_up   = self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_dn_%s_mj"%(self.channel))

         rrv_TTbar.Print();
         rrv_TTbarmassvbf_jes_up.Print();
         rrv_TTbarmassvbf_jes_dn.Print();        
         rrv_TTbarmassvbf_jer.Print();        
         rrv_TTbarmassvbf_jer_up.Print();        
         rrv_TTbarmassvbf_jer_dn.Print();        


         #jet mass uncertainty on TTbar normalization
         if(self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jes_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jes_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbarmassvbf_jer_%s_mj"%(self.channel))):    

            self.TTbar_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_TTbarmassvbf_jes_up.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassvbf_jes_dn.getVal()-rrv_TTbar.getVal() ) )/2.)/rrv_TTbar.getVal();         
            print "Total Uncertainty on TTbar due to jes: uncertainty ",self.TTbar_normalization_uncertainty_from_jet_scale;
            self.TTbar_normalization_uncertainty_from_jet_res   = ((TMath.Abs(rrv_TTbarmassvbf_jer.getVal()-rrv_TTbar.getVal())+TMath.Abs(rrv_TTbarmassvbf_jer_up.getVal()-rrv_TTbar.getVal() )+TMath.Abs(rrv_TTbarmassvbf_jer_dn.getVal()-rrv_TTbar.getVal() ) )/3.)/rrv_TTbar.getVal();         
            print "Total Uncertainty on TTbar due to jer: uncertainty ",self.TTbar_normalization_uncertainty_from_jet_res;


         #jet mass uncertainty on VV normalization
         rrv_VV                 = self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_%s_mj"%(self.channel))
         rrv_VVmassvbf_jes_up   = self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jes_up_%s_mj"%(self.channel))
         rrv_VVmassvbf_jes_dn   = self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jes_dn_%s_mj"%(self.channel))
         rrv_VVmassvbf_jer      = self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_%s_mj"%(self.channel))
         rrv_VVmassvbf_jer_up   = self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_up_%s_mj"%(self.channel))
         rrv_VVmassvbf_jer_dn   = self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_dn_%s_mj"%(self.channel))

         rrv_VV.Print();
         rrv_VVmassvbf_jes_up.Print();
         rrv_VVmassvbf_jes_dn.Print();        
         rrv_VVmassvbf_jer_up.Print();
         rrv_VVmassvbf_jer_dn.Print();        
         rrv_VVmassvbf_jer.Print();        


         #jet mass uncertainty on VV normalization
         if(self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jes_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jes_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_VVmassvbf_jer_%s_mj"%(self.channel))):    

            self.VV_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_VVmassvbf_jes_up.getVal()-rrv_VV.getVal())+TMath.Abs(rrv_VVmassvbf_jes_dn.getVal()-rrv_VV.getVal() ) )/2.)/rrv_VV.getVal();         
            print "Total Uncertainty on VV due to jes: uncertainty ",self.VV_normalization_uncertainty_from_jet_scale;
            self.VV_normalization_uncertainty_from_jet_res = ((TMath.Abs(rrv_VVmassvbf_jer_up.getVal()-rrv_VV.getVal())+TMath.Abs(rrv_VVmassvbf_jer_dn.getVal()-rrv_VV.getVal() )+TMath.Abs(rrv_VVmassvbf_jer.getVal()-rrv_VV.getVal() ) )/3.)/rrv_VV.getVal();         
            print "Total Uncertainty on VV due to jer: uncertainty ",self.VV_normalization_uncertainty_from_jet_res;
 
         #jet mass uncertainty on WW_EWK normalization
         if options.jetBin == "_2jet" :
          rrv_WW_EWK                = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWK_%s_mj"%(self.channel))
          rrv_WW_EWKmassvbf_jes_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jes_up_%s_mj"%(self.channel))
          rrv_WW_EWKmassvbf_jes_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jes_dn_%s_mj"%(self.channel))
          rrv_WW_EWKmassvbf_jer     = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_%s_mj"%(self.channel))
          rrv_WW_EWKmassvbf_jer_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_up_%s_mj"%(self.channel))
          rrv_WW_EWKmassvbf_jer_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_dn_%s_mj"%(self.channel))

          rrv_WW_EWK.Print();
          rrv_WW_EWKmassvbf_jes_up.Print();
          rrv_WW_EWKmassvbf_jes_dn.Print();
          rrv_WW_EWKmassvbf_jer_up.Print();
          rrv_WW_EWKmassvbf_jer_dn.Print();
          rrv_WW_EWKmassvbf_jer.Print();

          #jet mass uncertainty on WW_EWK normalization
          if(self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jes_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jes_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_up_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_dn_%s_mj"%(self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_WW_EWKmassvbf_jer_%s_mj"%(self.channel))):    

            self.WW_EWK_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_WW_EWKmassvbf_jes_up.getVal()-rrv_WW_EWK.getVal())+TMath.Abs(rrv_WW_EWKmassvbf_jes_dn.getVal()-rrv_WW_EWK.getVal() ) )/2.)/rrv_WW_EWK.getVal();         
            print "Total Uncertainty on WW_EWK due to jes: uncertainty ",self.WW_EWK_normalization_uncertainty_from_jet_scale;
            self.WW_EWK_normalization_uncertainty_from_jet_res = ((TMath.Abs(rrv_WW_EWKmassvbf_jer_up.getVal()-rrv_WW_EWK.getVal())+TMath.Abs(rrv_WW_EWKmassvbf_jer_dn.getVal()-rrv_WW_EWK.getVal() )+TMath.Abs(rrv_WW_EWKmassvbf_jer.getVal()-rrv_WW_EWK.getVal() ) )/3.)/rrv_WW_EWK.getVal();         
            print "Total Uncertainty on WW_EWK due to jer: uncertainty ",self.WW_EWK_normalization_uncertainty_from_jet_res;


         #jet mass uncertainty on ggH normalization
         rrv_ggH                  = self.workspace4fit_.var("rrv_number_dataset_signal_region_%s_%s_mj"%(self.higgs_sample,self.channel))
         rrv_ggHmassvbf_jes_up    = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_up_%s_mj"%(self.higgs_sample,self.channel))
         rrv_ggHmassvbf_jes_dn    = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_dn_%s_mj"%(self.higgs_sample,self.channel))
         rrv_ggHmassvbf_jer       = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_%s_mj"%(self.higgs_sample,self.channel))
         rrv_ggHmassvbf_jer_up    = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_up_%s_mj"%(self.higgs_sample,self.channel))
         rrv_ggHmassvbf_jer_dn    = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_dn_%s_mj"%(self.higgs_sample,self.channel))

         rrv_ggH.Print();
         rrv_ggHmassvbf_jes_up.Print();
         rrv_ggHmassvbf_jes_dn.Print();        
         rrv_ggHmassvbf_jer.Print();        
         rrv_ggHmassvbf_jer_up.Print();        
         rrv_ggHmassvbf_jer_dn.Print();        

         #jet mass uncertainty on ggH normalization
         if(self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_up_%s_mj"%(self.higgs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_up_%s_mj"%(self.higgs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_up_%s_mj"%(self.higgs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_dn_%s_mj"%(self.higgs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_%s_mj"%(self.higgs_sample,self.channel))):    

            self.ggH_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_ggHmassvbf_jes_up.getVal()-rrv_ggH.getVal())+TMath.Abs(rrv_ggHmassvbf_jes_dn.getVal()-rrv_ggH.getVal() ) )/2.)/rrv_ggH.getVal();         
            print "Total Uncertainty on ggH due to jes: uncertainty ",self.ggH_normalization_uncertainty_from_jet_scale;
            self.ggH_normalization_uncertainty_from_jet_res = ((TMath.Abs(rrv_ggHmassvbf_jer_up.getVal()-rrv_ggH.getVal())+TMath.Abs(rrv_ggHmassvbf_jer_dn.getVal()-rrv_ggH.getVal() )+TMath.Abs(rrv_ggHmassvbf_jer.getVal()-rrv_ggH.getVal() ) )/3.)/rrv_ggH.getVal();         
            print "Total Uncertainty on ggH due to jer: uncertainty ",self.ggH_normalization_uncertainty_from_jet_res;


         #jet mass uncertainty on vbf normalizatio
         rrv_vbf                = self.workspace4fit_.var("rrv_number_dataset_signal_region_%s_%s_mj"%(self.vbfhiggs_sample,self.channel))
         rrv_vbfmassvbf_jes_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_up_%s_mj"%(self.vbfhiggs_sample,self.channel))
         rrv_vbfmassvbf_jes_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_dn_%s_mj"%(self.vbfhiggs_sample,self.channel))
         rrv_vbfmassvbf_jer_up  = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_up_%s_mj"%(self.vbfhiggs_sample,self.channel))
         rrv_vbfmassvbf_jer_dn  = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_dn_%s_mj"%(self.vbfhiggs_sample,self.channel))
         rrv_vbfmassvbf_jer     = self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_%s_mj"%(self.vbfhiggs_sample,self.channel))

         rrv_vbf.Print();
         rrv_vbfmassvbf_jes_up.Print();
         rrv_vbfmassvbf_jes_dn.Print();        
         rrv_vbfmassvbf_jer_up.Print();        
         rrv_vbfmassvbf_jer_dn.Print();        
         rrv_vbfmassvbf_jer.Print();        

         #jet mass uncertainty on vbf normalization
         if(self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_up_%s_mj"%(self.vbfhiggs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jes_dn_%s_mj"%(self.vbfhiggs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_up_%s_mj"%(self.vbfhiggs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_dn_%s_mj"%(self.vbfhiggs_sample,self.channel)) and self.workspace4fit_.var("rrv_number_dataset_signal_region_%smassvbf_jer_%s_mj"%(self.vbfhiggs_sample,self.channel))):
             
            self.vbf_normalization_uncertainty_from_jet_scale = ((TMath.Abs(rrv_vbfmassvbf_jes_up.getVal()-rrv_vbf.getVal())+TMath.Abs(rrv_vbfmassvbf_jes_dn.getVal()-rrv_vbf.getVal() ) )/2.)/rrv_vbf.getVal();         
            print "Total Uncertainty on vbfH due to jes: uncertainty ",self.vbf_normalization_uncertainty_from_jet_scale;
            self.vbf_normalization_uncertainty_from_jet_res   = ((TMath.Abs(rrv_vbfmassvbf_jer_up.getVal()-rrv_vbf.getVal())+TMath.Abs(rrv_vbfmassvbf_jer_dn.getVal()-rrv_vbf.getVal() )+TMath.Abs(rrv_vbfmassvbf_jer.getVal()-rrv_vbf.getVal() ) )/3.)/rrv_vbf.getVal();         
            print "Total Uncertainty on vbfH due to jer: uncertainty ",self.vbf_normalization_uncertainty_from_jet_res;


    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass="jet_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");

        rrv_mass_j      = self.workspace4fit_.var("rrv_mass_j") 
        rrv_mass_lvj    = self.workspace4fit_.var("rrv_mass_lvj")
        rrv_mass_gen_WW = self.workspace4fit_.var("rrv_mass_gen_WW")

        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 
        
        # dataset of m_j -> before and after vbf cuts -> central object value
        rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        if TString(label).Contains("ggH") or TString(label).Contains("vbfH"):
          rdataset4fit_m_WW_gen = RooDataSet("rdataset4fit"+label+"_genHMass_"+self.channel,"rdataset4fit"+label+"_genHMass_"+self.channel,RooArgSet(rrv_mass_gen_WW,rrv_weight),RooFit.WeightVar(rrv_weight));
                            
        #dataset of m_lvj -> before and after vbf cuts -> central object value
        rdataset_sb_lo_mlvj     = RooDataSet("rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_signal_region_mlvj = RooDataSet("rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset_sb_hi_mlvj         = RooDataSet("rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        rdataset4fit_sb_lo_mlvj     = RooDataSet("rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4fit_signal_region_mlvj = RooDataSet("rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
        rdataset4fit_sb_hi_mlvj     = RooDataSet("rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        
        if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
        
         #dataset of jes_up
         rdataset_mj_jes_up  = RooDataSet("rdataset"+label+"massvbf_jes_up"+"_"+self.channel+"_mj","rdataset"+label+"massvbf_jes_up"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         rdataset4fit_mj_jes_up = RooDataSet("rdataset4fit"+label+"massvbf_jes_up"+"_"+self.channel+"_mj","rdataset4fit"+label+"massvbf_jes_up"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         
         rdataset_sb_lo_mlvj_jes_up         = RooDataSet("rdataset"+label+"massvbf_jes_up"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_up"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_signal_region_mlvj_jes_up = RooDataSet("rdataset"+label+"massvbf_jes_up"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_up"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_sb_hi_mlvj_jes_up         = RooDataSet("rdataset"+label+"massvbf_jes_up"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_up"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         rdataset4fit_sb_lo_mlvj_jes_up  = RooDataSet("rdataset4fit"+label+"massvbf_jes_up"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_up"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_signal_region_mlvj_jes_up = RooDataSet("rdataset4fit"+label+"massvbf_jes_up"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_up"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_sb_hi_mlvj_jes_up  = RooDataSet("rdataset4fit"+label+"massvbf_jes_up"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_up"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         #dataset of applying jes dn
         rdataset_mj_jes_dn  = RooDataSet("rdataset"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj","rdataset"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         rdataset4fit_mj_jes_dn = RooDataSet("rdataset4fit"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj","rdataset4fit"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

         rdataset_sb_lo_mlvj_jes_dn         = RooDataSet("rdataset"+label+"massvbf_jes_dn"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_dn"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_signal_region_mlvj_jes_dn = RooDataSet("rdataset"+label+"massvbf_jes_dn"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_dn"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_sb_hi_mlvj_jes_dn         = RooDataSet("rdataset"+label+"massvbf_jes_dn"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jes_dn"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         rdataset4fit_sb_lo_mlvj_jes_dn  = RooDataSet("rdataset4fit"+label+"massvbf_jes_dn"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_dn"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_signal_region_mlvj_jes_dn = RooDataSet("rdataset4fit"+label+"massvbf_jes_dn"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_dn"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_sb_hi_mlvj_jes_dn     = RooDataSet("rdataset4fit"+label+"massvbf_jes_dn"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jes_dn"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 


         #dataset of applying jer up
         rdataset_mj_jer_up  = RooDataSet("rdataset"+label+"massvbf_jer_up"+"_"+self.channel+"_mj","rdataset"+label+"massvbf_jer_up"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         rdataset4fit_mj_jer_up = RooDataSet("rdataset4fit"+label+"massvbf_jer_up"+"_"+self.channel+"_mj","rdataset4fit"+label+"massvbf_jer_up"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

         rdataset_sb_lo_mlvj_jer_up         = RooDataSet("rdataset"+label+"massvbf_jer_up"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_up"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_signal_region_mlvj_jer_up = RooDataSet("rdataset"+label+"massvbf_jer_up"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_up"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_sb_hi_mlvj_jer_up         = RooDataSet("rdataset"+label+"massvbf_jer_up"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_up"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         rdataset4fit_sb_lo_mlvj_jer_up  = RooDataSet("rdataset4fit"+label+"massvbf_jer_up"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_up"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_signal_region_mlvj_jer_up = RooDataSet("rdataset4fit"+label+"massvbf_jer_up"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_up"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_sb_hi_mlvj_jer_up     = RooDataSet("rdataset4fit"+label+"massvbf_jer_up"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_up"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         #dataset of applying jer dn
         rdataset_mj_jer_dn  = RooDataSet("rdataset"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj","rdataset"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         rdataset4fit_mj_jer_dn = RooDataSet("rdataset4fit"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj","rdataset4fit"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

         rdataset_sb_lo_mlvj_jer_dn         = RooDataSet("rdataset"+label+"massvbf_jer_dn"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_dn"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_signal_region_mlvj_jer_dn = RooDataSet("rdataset"+label+"massvbf_jer_dn"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_dn"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_sb_hi_mlvj_jer_dn         = RooDataSet("rdataset"+label+"massvbf_jer_dn"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer_dn"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         rdataset4fit_sb_lo_mlvj_jer_dn  = RooDataSet("rdataset4fit"+label+"massvbf_jer_dn"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_dn"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_signal_region_mlvj_jer_dn = RooDataSet("rdataset4fit"+label+"massvbf_jer_dn"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_dn"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_sb_hi_mlvj_jer_dn     = RooDataSet("rdataset4fit"+label+"massvbf_jer_dn"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer_dn"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        
         #dataset of applying jer 
         rdataset_mj_jer  = RooDataSet("rdataset"+label+"massvbf_jer"+"_"+self.channel+"_mj","rdataset"+label+"massvbf_jer"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
         rdataset4fit_mj_jer = RooDataSet("rdataset4fit"+label+"massvbf_jer"+"_"+self.channel+"_mj","rdataset4fit"+label+"massvbf_jer"+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

         rdataset_sb_lo_mlvj_jer         = RooDataSet("rdataset"+label+"massvbf_jer"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_signal_region_mlvj_jer = RooDataSet("rdataset"+label+"massvbf_jer"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset_sb_hi_mlvj_jer         = RooDataSet("rdataset"+label+"massvbf_jer"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset"+label+"massvbf_jer"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

         rdataset4fit_sb_lo_mlvj_jer  = RooDataSet("rdataset4fit"+label+"massvbf_jer"+"_sb_lo"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer"+"_sb_lo"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_signal_region_mlvj_jer = RooDataSet("rdataset4fit"+label+"massvbf_jer"+"_signal_region"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer"+"_signal_region"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 
         rdataset4fit_sb_hi_mlvj_jer     = RooDataSet("rdataset4fit"+label+"massvbf_jer"+"_sb_hi"+"_"+self.channel+"_mlvj","rdataset4fit"+label+"massvbf_jer"+"_sb_hi"+"_"+self.channel+"_mlvj",RooArgSet(rrv_mass_lvj,rrv_weight),RooFit.WeightVar(rrv_weight) ); 

        ###### Define the event categorization
        data_category = RooCategory("data_category","data_category");
        data_category.defineType("sideband");
        data_category.defineType("signal_region");

        combData = RooDataSet("combData"+label+"_"+self.channel,"combData"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
        combData4fit = RooDataSet("combData4fit"+label+"_"+self.channel,"combData4fit"+label+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );


        if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
          
         ## jes_up
         combData_jes_up = RooDataSet("combData"+label+"massvbf_jes_up"+"_"+self.channel,"combData"+label+"massvbf_jes_up"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
         combData4fit_jes_up = RooDataSet("combData4fit"+label+"massvbf_jes_up"+"_"+self.channel,"combData4fit"+label+"massvbf_jes_up"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

         ## jes_dn
         combData_jes_dn = RooDataSet("combData"+label+"massvbf_jes_dn"+"_"+self.channel,"combData"+label+"massvbf_jes_dn"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
         combData4fit_jes_dn = RooDataSet("combData4fit"+label+"massvbf_jes_dn"+"_"+self.channel,"combData4fit"+label+"massvbf_jes_dn"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

         ## jer_up
         combData_jer_up = RooDataSet("combData"+label+"massvbf_jer_up"+"_"+self.channel,"combData"+label+"massvbf_jer_up"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
         combData4fit_jer_up = RooDataSet("combData4fit"+label+"massvbf_jer_up"+"_"+self.channel,"combData4fit"+label+"massvbf_jer_up"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

         ## jer_dn
         combData_jer_dn = RooDataSet("combData"+label+"massvbf_jer_dn"+"_"+self.channel,"combData"+label+"massvbf_jer_dn"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
         combData4fit_jer_dn = RooDataSet("combData4fit"+label+"massvbf_jer_dn"+"_"+self.channel,"combData4fit"+label+"massvbf_jer_dn"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

         ## jer
         combData_jer = RooDataSet("combData"+label+"massvbf_jer"+"_"+self.channel,"combData"+label+"massvbf_jer"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );
         combData4fit_jer = RooDataSet("combData4fit"+label+"massvbf_jer"+"_"+self.channel,"combData4fit"+label+"massvbf_jer"+"_"+self.channel,RooArgSet(rrv_mass_lvj, data_category, rrv_weight),RooFit.WeightVar(rrv_weight) );

        print "N entries: ", treeIn.GetEntries();

        hnum_4region = TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j   -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total        
        hnum_2region = TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj  0: signal_region; 1: total

        if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:

         hnum_4region_jes_up = TH1D("hnum_4region"+label+"massvbf_jes_up"+"_"+self.channel,"hnum_4region"+label+"massvbf_jes_up"+"_"+self.channel,4,-1.5,2.5);  
         hnum_2region_jes_up = TH1D("hnum_2region"+label+"massvbf_jes_up"+"_"+self.channel,"hnum_2region"+label+"massvbf_jes_up"+"_"+self.channel,2,-0.5,1.5);

         hnum_4region_jes_dn = TH1D("hnum_4region"+label+"massvbf_jes_dn"+"_"+self.channel,"hnum_4region"+label+"massvbf_jes_dn"+"_"+self.channel,4,-1.5,2.5);  
         hnum_2region_jes_dn = TH1D("hnum_2region"+label+"massvbf_jes_dn"+"_"+self.channel,"hnum_2region"+label+"massvbf_jes_dn"+"_"+self.channel,2,-0.5,1.5);

         hnum_4region_jer_up = TH1D("hnum_4region"+label+"massvbf_jer_up"+"_"+self.channel,"hnum_4region"+label+"massvbf_jer_up"+"_"+self.channel,4,-1.5,2.5);  
         hnum_2region_jer_up = TH1D("hnum_2region"+label+"massvbf_jer_up"+"_"+self.channel,"hnum_2region"+label+"massvbf_jer_up"+"_"+self.channel,2,-0.5,1.5);

         hnum_4region_jer_dn = TH1D("hnum_4region"+label+"massvbf_jer_dn"+"_"+self.channel,"hnum_4region"+label+"massvbf_jer_dn"+"_"+self.channel,4,-1.5,2.5);  
         hnum_2region_jer_dn = TH1D("hnum_2region"+label+"massvbf_jer_dn"+"_"+self.channel,"hnum_2region"+label+"massvbf_jer_dn"+"_"+self.channel,2,-0.5,1.5);

         hnum_4region_jer = TH1D("hnum_4region"+label+"massvbf_jer"+"_"+self.channel,"hnum_4region"+label+"massvbf_jer"+"_"+self.channel,4,-1.5,2.5);  
         hnum_2region_jer = TH1D("hnum_2region"+label+"massvbf_jer"+"_"+self.channel,"hnum_2region"+label+"massvbf_jer"+"_"+self.channel,2,-0.5,1.5);

        tmp_scale_to_lumi = 0 ;
        
        for i in range(treeIn.GetEntries()):


          if i % 10000 == 0: print "iEvent: ",i
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

          mass_WW_gen = 0 ;
          if TString(label).Contains("ggH") or TString(label).Contains("vbfH"):
           mass_WW_gen = getattr(treeIn,"genHMass");
                              

          njet = 0. ; tmp_vbf_dEta =0.; tmp_vbf_Mjj = 0.; ungroomed_jet_pt = 0.; pfMET = 0.; mass_lvj = 0. ;

          # jet mass , central value
          tmp_jet_mass = getattr(treeIn, jet_mass);
          tmp_vbf_dEta = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta")-getattr(treeIn,"vbf_maxpt_j2_eta"));
          tmp_vbf_Mjj  = getattr(treeIn, "vbf_maxpt_jj_m");
          njet         = getattr(treeIn,"numberJetBin");
          ungroomed_jet_pt = getattr(treeIn,"ungroomed_jet_pt");
          pfMET    = getattr(treeIn,"pfMET");
          mass_lvj = getattr(treeIn,"mass_lvj_type0_met");

          if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:

           # jet mass jes_up
           tmp_jet_mass_jes_up = getattr(treeIn, "jet_mass_pr_jes_up");
           tmp_vbf_dEta_jes_up = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta_jes_up")-getattr(treeIn,"vbf_maxpt_j2_eta_jes_up"));
           jet_1.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j1_pt_jes_up"),getattr(treeIn,"vbf_maxpt_j1_eta_jes_up"),getattr(treeIn,"vbf_maxpt_j1_phi_jes_up"),getattr(treeIn,"vbf_maxpt_j1_m_jes_up"));
           jet_2.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j2_pt_jes_up"),getattr(treeIn,"vbf_maxpt_j2_eta_jes_up"),getattr(treeIn,"vbf_maxpt_j2_phi_jes_up"),getattr(treeIn,"vbf_maxpt_j2_m_jes_up"));
           tmp_vbf_Mjj_jes_up      = (jet_1+jet_2).M();
           ungroomed_jet_pt_jes_up = getattr(treeIn,"ungroomed_jet_pt_jes_up");
           njet_jes_up = 0 ;
           if(getattr(treeIn,"vbf_maxpt_j1_pt_jes_up") > 30.):
               njet_jes_up = njet_jes_up +1 ;
           if(getattr(treeIn,"vbf_maxpt_j2_pt_jes_up") > 30.):
               njet_jes_up = njet_jes_up +1 ;
           pfMET_jes_up = getattr(treeIn,"pfMET_jes_up");
           mass_lvj_jes_up = getattr(treeIn,"mass_lvj_type0_met_jes_up");

           # jet mass jes_dn
           tmp_jet_mass_jes_dn = getattr(treeIn, "jet_mass_pr_jes_dn");
           tmp_vbf_dEta_jes_dn = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta_jes_dn")-getattr(treeIn,"vbf_maxpt_j2_eta_jes_dn"));
           jet_1.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j1_pt_jes_dn"),getattr(treeIn,"vbf_maxpt_j1_eta_jes_dn"),getattr(treeIn,"vbf_maxpt_j1_phi_jes_dn"),getattr(treeIn,"vbf_maxpt_j1_m_jes_dn"));
           jet_2.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j2_pt_jes_dn"),getattr(treeIn,"vbf_maxpt_j2_eta_jes_dn"),getattr(treeIn,"vbf_maxpt_j2_phi_jes_dn"),getattr(treeIn,"vbf_maxpt_j2_m_jes_dn"));
           tmp_vbf_Mjj_jes_dn = (jet_1 + jet_2).M();
           ungroomed_jet_pt_jes_dn = getattr(treeIn,"ungroomed_jet_pt_jes_dn");
           njet_jes_dn = 0 ;
           if(getattr(treeIn,"vbf_maxpt_j1_pt_jes_dn") > 30.):
                njet_jes_dn = njet_jes_dn +1 ;
           if(getattr(treeIn,"vbf_maxpt_j2_pt_jes_dn") > 30.):
               njet_jes_dn = njet_jes_dn +1 ;
           pfMET_jes_dn = getattr(treeIn,"pfMET_jes_dn");
           mass_lvj_jes_dn = getattr(treeIn,"mass_lvj_type0_met_jes_dn");

           #jet mass jer
           tmp_jet_mass_jer = getattr(treeIn, "jet_mass_pr_jer");
           tmp_vbf_dEta_jer = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta_jer")-getattr(treeIn,"vbf_maxpt_j2_eta_jer"));
           jet_1.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j1_pt_jer"),getattr(treeIn,"vbf_maxpt_j1_eta_jer"),getattr(treeIn,"vbf_maxpt_j1_phi_jer"),getattr(treeIn,"vbf_maxpt_j1_m_jer"));
           jet_2.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j2_pt_jer"),getattr(treeIn,"vbf_maxpt_j2_eta_jer"),getattr(treeIn,"vbf_maxpt_j2_phi_jer"),getattr(treeIn,"vbf_maxpt_j2_m_jer"));
           tmp_vbf_Mjj_jer = (jet_1+jet_2).M();
           ungroomed_jet_pt_jer = getattr(treeIn,"ungroomed_jet_pt_jer");
           njet_jer = 0 ;
           if(getattr(treeIn,"vbf_maxpt_j1_pt_jer") > 30.):
               njet_jer = njet_jer +1 ;
           if(getattr(treeIn,"vbf_maxpt_j2_pt_jer") > 30.):
               njet_jer = njet_jer +1 ;
           pfMET_jer = getattr(treeIn,"pfMET_jer");
           mass_lvj_jer = getattr(treeIn,"mass_lvj_type0_met_jer");

           #jet mass jer up
           tmp_jet_mass_jer_up = getattr(treeIn, "jet_mass_pr_jer_up");
           tmp_vbf_dEta_jer_up = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta_jer_up")-getattr(treeIn,"vbf_maxpt_j2_eta_jer_up"));
           jet_1.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j1_pt_jer_up"),getattr(treeIn,"vbf_maxpt_j1_eta_jer_up"),getattr(treeIn,"vbf_maxpt_j1_phi_jer_up"),getattr(treeIn,"vbf_maxpt_j1_m_jer_up"));
           jet_2.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j2_pt_jer_up"),getattr(treeIn,"vbf_maxpt_j2_eta_jer_up"),getattr(treeIn,"vbf_maxpt_j2_phi_jer_up"),getattr(treeIn,"vbf_maxpt_j2_m_jer_up"));
           tmp_vbf_Mjj_jer_up = (jet_1+jet_2).M();
           ungroomed_jet_pt_jer_up = getattr(treeIn,"ungroomed_jet_pt_jer_up");
           njet_jer_up = 0;
           if(getattr(treeIn,"vbf_maxpt_j1_pt_jer_up") > 30.):
               njet_jer_up = njet_jer_up +1 ;
           if(getattr(treeIn,"vbf_maxpt_j2_pt_jer_up") > 30.):
               njet_jer_up = njet_jer_up +1 ;
           pfMET_jer_up = getattr(treeIn,"pfMET_jer_up");
           mass_lvj_jer_up = getattr(treeIn,"mass_lvj_type0_met_jer_up");

           #jet mass jer dn
           tmp_jet_mass_jer_dn = getattr(treeIn, "jet_mass_pr_jer_dn");
           tmp_vbf_dEta_jer_dn = math.fabs(getattr(treeIn, "vbf_maxpt_j1_eta_jer_dn")-getattr(treeIn,"vbf_maxpt_j2_eta_jer_dn"));
           jet_1.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j1_pt_jer_dn"),getattr(treeIn,"vbf_maxpt_j1_eta_jer_dn"),getattr(treeIn,"vbf_maxpt_j1_phi_jer_dn"),getattr(treeIn,"vbf_maxpt_j1_m_jer_dn"));
           jet_2.SetPtEtaPhiM(getattr(treeIn,"vbf_maxpt_j2_pt_jer_dn"),getattr(treeIn,"vbf_maxpt_j2_eta_jer_dn"),getattr(treeIn,"vbf_maxpt_j2_phi_jer_dn"),getattr(treeIn,"vbf_maxpt_j2_m_jer_dn"));
           tmp_vbf_Mjj_jer_dn = (jet_1 + jet_2).M();
           ungroomed_jet_pt_jer_dn = getattr(treeIn,"ungroomed_jet_pt_jer_dn");
           njet_jer_dn = 0 ;
           if(getattr(treeIn,"vbf_maxpt_j1_pt_jer_dn") > 30.):
               njet_jer_dn = njet_jer_dn +1 ;
           if(getattr(treeIn,"vbf_maxpt_j2_pt_jer_dn") > 30.):
               njet_jer_dn = njet_jer_dn +1 ;
           pfMET_jer_dn = getattr(treeIn,"pfMET_jer_dn");
           mass_lvj_jer_dn = getattr(treeIn,"mass_lvj_type0_met_jer_dn");

          isPassingCut = 0 ;
          isPassingCut_jes_up = 0;
          isPassingCut_jes_dn = 0;
          isPassingCut_jer = 0;
          isPassingCut_jer_up = 0;
          isPassingCut_jer_dn = 0;

          if options.jetBin == "_2jet" :
           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
            isPassingCut = 1 ;

           if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
           
            if ungroomed_jet_pt_jes_up > 200. and discriminantCut and tmp_jet_mass_jes_up >= rrv_mass_j.getMin() and tmp_jet_mass_jes_up<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj_jes_up >= rrv_mass_lvj.getMin() and mass_lvj_jes_up <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jes_up > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet_jes_up >=2 and tmp_vbf_dEta_jes_up > self.dEta_cut and tmp_vbf_Mjj_jes_up > self.Mjj_cut :
              isPassingCut_jes_up = 1 ;

          
            if ungroomed_jet_pt_jes_dn > 200. and discriminantCut and tmp_jet_mass_jes_dn >= rrv_mass_j.getMin() and tmp_jet_mass_jes_dn<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj_jes_dn >= rrv_mass_lvj.getMin() and mass_lvj_jes_dn <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jes_dn > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet_jes_dn >=2 and tmp_vbf_dEta_jes_dn > self.dEta_cut and tmp_vbf_Mjj_jes_dn > self.Mjj_cut:
             isPassingCut_jes_dn = 1 ;

            if ungroomed_jet_pt_jer > 200. and discriminantCut and tmp_jet_mass_jer >= rrv_mass_j.getMin() and tmp_jet_mass_jer<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj_jer >= rrv_mass_lvj.getMin() and mass_lvj_jer <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet_jer >=2 and tmp_vbf_dEta_jer > self.dEta_cut and tmp_vbf_Mjj_jer > self.Mjj_cut:          
             isPassingCut_jer = 1;

            if ungroomed_jet_pt_jer_up > 200. and discriminantCut and tmp_jet_mass_jer_up >= rrv_mass_j.getMin() and tmp_jet_mass_jer_up<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj_jer_up >= rrv_mass_lvj.getMin() and mass_lvj_jer_up <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer_up > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet_jer_up >=2 and tmp_vbf_dEta_jer_up > self.dEta_cut and tmp_vbf_Mjj_jer_up > self.Mjj_cut:          
             isPassingCut_jer_up = 1 ;

            if ungroomed_jet_pt_jer_dn > 200. and discriminantCut and tmp_jet_mass_jer_dn >= rrv_mass_j.getMin() and tmp_jet_mass_jer_dn<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj_jer_dn >= rrv_mass_lvj.getMin() and mass_lvj_jer_dn <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer_dn > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had and getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep and njet_jer_dn >=2 and tmp_vbf_dEta_jer_dn > self.dEta_cut and tmp_vbf_Mjj_jer_dn > self.Mjj_cut:          
             isPassingCut_jer_dn = 1 ;

          else:
              
           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet < 2:
               
            isPassingCut = 1 ;

           if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:

            if ungroomed_jet_pt_jes_up > 200. and discriminantCut and tmp_jet_mass_jes_up >= rrv_mass_j.getMin() and tmp_jet_mass_jes_up<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj_jes_up >= rrv_mass_lvj.getMin() and mass_lvj_jes_up <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jes_up > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet_jes_up < 2:

             isPassingCut_jes_up = 1 ;
          
            if ungroomed_jet_pt_jes_dn > 200. and discriminantCut and tmp_jet_mass_jes_dn >= rrv_mass_j.getMin() and tmp_jet_mass_jes_dn<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj_jes_dn >= rrv_mass_lvj.getMin() and mass_lvj_jes_dn <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jes_dn > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet_jes_dn < 2:
             isPassingCut_jes_dn = 1 ;

            if ungroomed_jet_pt_jer > 200. and discriminantCut and tmp_jet_mass_jer >= rrv_mass_j.getMin() and tmp_jet_mass_jer<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj_jer >= rrv_mass_lvj.getMin() and mass_lvj_jer <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet_jer < 2:
             isPassingCut_jer = 1;

            if ungroomed_jet_pt_jer_up > 200. and discriminantCut and tmp_jet_mass_jer_up >= rrv_mass_j.getMin() and tmp_jet_mass_jer_up<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj_jer_up >= rrv_mass_lvj.getMin() and mass_lvj_jer_up <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer_up > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet_jer_up < 2:
             isPassingCut_jer_up = 1 ;

            if ungroomed_jet_pt_jer_dn > 200. and discriminantCut and tmp_jet_mass_jer_dn >= rrv_mass_j.getMin() and tmp_jet_mass_jer_dn<=rrv_mass_j.getMax() and getattr(treeIn,"nbjets_csvm_veto") < 1 and mass_lvj_jer_dn >= rrv_mass_lvj.getMin() and mass_lvj_jer_dn <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET_jer_dn > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and njet_jer_dn < 2:
             isPassingCut_jer_dn = 1 ;

          if isPassingCut !=0 or isPassingCut_jes_up!=0 or isPassingCut_jes_dn!=0 or isPassingCut_jer!=0 or isPassingCut_jer_up!=0 or isPassingCut_jer_dn!=0 :
             
             tmp_event_weight = getattr(treeIn,"totalEventWeight");
             tmp_event_weight4fit = getattr(treeIn,"eff_and_pu_Weight");
             tmp_interference_weight_H600 = getattr(treeIn,"interference_Weight_H600");
             tmp_interference_weight_H700 = getattr(treeIn,"interference_Weight_H700");
             tmp_interference_weight_H800 = getattr(treeIn,"interference_Weight_H800");
             tmp_interference_weight_H900 = getattr(treeIn,"interference_Weight_H900");
             tmp_interference_weight_H1000 = getattr(treeIn,"interference_Weight_H1000");
             ## added by Nhan, getting additional BSM weight
             bsmWeightName = "bsmReweight_cPrime%02d_brNew%02d"%(options.cprime,options.BRnew);
             tmp_bsmWeight = getattr(treeIn, bsmWeightName);
             if tmp_bsmWeight < 0: tmp_bsmWeight = 1;
             
             if TString(label).Contains("ggH600"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H600
                tmp_event_weight4fit = tmp_event_weight4fit*tmp_interference_weight_H600
             if TString(label).Contains("ggH700"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H700
                tmp_event_weight4fit = tmp_event_weight4fit*tmp_interference_weight_H700
             if TString(label).Contains("ggH800"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H800
                tmp_event_weight4fit = tmp_event_weight4fit*tmp_interference_weight_H800
             if TString(label).Contains("ggH900"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H900
                tmp_event_weight4fit = tmp_event_weight4fit*tmp_interference_weight_H900
             if TString(label).Contains("ggH1000"):
                tmp_event_weight = tmp_event_weight*tmp_interference_weight_H1000
                tmp_event_weight4fit = tmp_event_weight4fit*tmp_interference_weight_H1000

             if TString(label).Contains("vbfH600"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H600;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.cps_Weight_H600;
             if TString(label).Contains("vbfH700"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H700;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.cps_Weight_H700;
             if TString(label).Contains("vbfH800"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H800;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.cps_Weight_H800;
             if TString(label).Contains("vbfH900"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H900;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.cps_Weight_H900;
             if TString(label).Contains("vbfH1000"):
                tmp_event_weight = tmp_event_weight*treeIn.cps_Weight_H1000;
                tmp_event_weight4fit = tmp_event_weight4fit*treeIn.cps_Weight_H1000;
 
             # for multi-sample, like STop and VV. There are two sample, and two wSampleWeight_value.Use the least wSampleWeight as scale.
             tmp_event_weight4fit = tmp_event_weight4fit*treeIn.wSampleWeight/tmp_scale_to_lumi;
             if TString(label).Contains("ggH") or TString(label).Contains("vbfH"):
                 tmp_event_weight = tmp_event_weight/self.higgs_xs_scale;
                 tmp_event_weight4fit = tmp_event_weight4fit/self.higgs_xs_scale;
                 ## added by Nhan
                 tmp_event_weight = tmp_event_weight*tmp_bsmWeight;
                 tmp_event_weight4fit=tmp_event_weight4fit*tmp_bsmWeight;
        
             if not label=="_data":
                     if TString(label).Contains("_TTbar") or TString(label).Contains("_STop") :
                         tmp_event_weight = tmp_event_weight*self.rrv_wtagger_eff_reweight_forT.getVal();
                         tmp_event_weight4fit = tmp_event_weight4fit*self.rrv_wtagger_eff_reweight_forT.getVal();
                     elif TString(label).Contains("_ggH") or TString(label).Contains("_vbfH") or TString(label).Contains("_VV") or TString(label).Contains("_WW_EWK") :
                         tmp_event_weight = tmp_event_weight*self.rrv_wtagger_eff_reweight_forV.getVal();
                         tmp_event_weight4fit = tmp_event_weight4fit*self.rrv_wtagger_eff_reweight_forV.getVal();

             tmp_event_weight        = tmp_event_weight* getattr(treeIn,"btag_weight"); ## add the btag weight 
             tmp_event_weight4fit    = tmp_event_weight4fit* getattr(treeIn,"btag_weight"); ## add the btag weight 

             if TString(label).Contains("ggH") or TString(label).Contains("vbfH"):
              rrv_mass_gen_WW.setVal(mass_WW_gen);
              rdataset4fit_m_WW_gen.add( RooArgSet( rrv_mass_gen_WW ), tmp_event_weight4fit );
                                                
             ## central values
             rrv_mass_lvj.setVal(mass_lvj);
             rrv_mass_j.setVal(tmp_jet_mass);

             if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max and isPassingCut == 1:
                 rdataset_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max and isPassingCut == 1:
                 rdataset_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region.Fill(1,tmp_event_weight);
                 if mass_lvj >=self.mlvj_signal_min and mass_lvj <self.mlvj_signal_max: 
                   hnum_2region.Fill(0,tmp_event_weight);

             if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max and isPassingCut == 1:
                 rdataset_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

             if isPassingCut == 1: 
              rdataset_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
              rdataset4fit_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
              if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max: 
                 hnum_4region.Fill(-1,tmp_event_weight );
              if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max : 
                 hnum_4region.Fill(0,tmp_event_weight);
              if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max: 
                 hnum_4region.Fill(1,tmp_event_weight);
              hnum_4region.Fill(2,tmp_event_weight);

             ## JES UP
             if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
             
              rrv_mass_lvj.setVal(mass_lvj_jes_up);
              rrv_mass_j.setVal(tmp_jet_mass_jes_up);

              if tmp_jet_mass_jes_up >= self.mj_sideband_lo_min and tmp_jet_mass_jes_up < self.mj_sideband_lo_max and isPassingCut_jes_up == 1:
                 rdataset_sb_lo_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData_jes_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jes_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

              if tmp_jet_mass_jes_up >= self.mj_signal_min and tmp_jet_mass_jes_up < self.mj_signal_max and isPassingCut_jes_up == 1:
                 rdataset_signal_region_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData_jes_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jes_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region_jes_up.Fill(1,tmp_event_weight);
                 if mass_lvj_jes_up >=self.mlvj_signal_min and mass_lvj_jes_up <self.mlvj_signal_max: 
                   hnum_2region_jes_up.Fill(0,tmp_event_weight);
                   
              if tmp_jet_mass_jes_up >= self.mj_sideband_hi_min and tmp_jet_mass_jes_up < self.mj_sideband_hi_max and isPassingCut_jes_up == 1:
                 rdataset_sb_hi_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj_jes_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

              if isPassingCut_jes_up == 1: 
               rdataset_mj_jes_up.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
               rdataset4fit_mj_jes_up.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
               
              if tmp_jet_mass_jes_up >=self.mj_sideband_lo_min and tmp_jet_mass_jes_up <self.mj_sideband_lo_max: 
                 hnum_4region_jes_up.Fill(-1,tmp_event_weight );
              if tmp_jet_mass_jes_up >=self.mj_signal_min and tmp_jet_mass_jes_up <self.mj_signal_max : 
                 hnum_4region_jes_up.Fill(0,tmp_event_weight);
              if tmp_jet_mass_jes_up >=self.mj_sideband_hi_min and tmp_jet_mass_jes_up <self.mj_sideband_hi_max: 
                 hnum_4region_jes_up.Fill(1,tmp_event_weight);
              hnum_4region_jes_up.Fill(2,tmp_event_weight);
              
 
              ### jes dn
              rrv_mass_lvj.setVal(mass_lvj_jes_dn);
              rrv_mass_j.setVal(tmp_jet_mass_jes_dn);

              if tmp_jet_mass_jes_dn >= self.mj_sideband_lo_min and tmp_jet_mass_jes_dn < self.mj_sideband_lo_max and isPassingCut_jes_dn == 1:
                 rdataset_sb_lo_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData_jes_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jes_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

              if tmp_jet_mass_jes_dn >= self.mj_signal_min and tmp_jet_mass_jes_dn < self.mj_signal_max and isPassingCut_jes_dn == 1:
                 rdataset_signal_region_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData_jes_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jes_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region_jes_dn.Fill(1,tmp_event_weight);
                 if mass_lvj_jes_dn >=self.mlvj_signal_min and mass_lvj_jes_dn <self.mlvj_signal_max: 
                   hnum_2region_jes_dn.Fill(0,tmp_event_weight);
                   

              if tmp_jet_mass_jes_dn >= self.mj_sideband_hi_min and tmp_jet_mass_jes_dn < self.mj_sideband_hi_max and isPassingCut_jes_dn == 1:
                 rdataset_sb_hi_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj_jes_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

              if isPassingCut_jes_dn == 1: 
               rdataset_mj_jes_dn.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
               rdataset4fit_mj_jes_dn.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
               if tmp_jet_mass_jes_dn >=self.mj_sideband_lo_min and tmp_jet_mass_jes_dn <self.mj_sideband_lo_max: 
                 hnum_4region_jes_dn.Fill(-1,tmp_event_weight );
               if tmp_jet_mass_jes_dn >=self.mj_signal_min and tmp_jet_mass_jes_dn <self.mj_signal_max : 
                  hnum_4region_jes_dn.Fill(0,tmp_event_weight);
               if tmp_jet_mass_jes_dn >=self.mj_sideband_hi_min and tmp_jet_mass_jes_dn <self.mj_sideband_hi_max: 
                 hnum_4region_jes_dn.Fill(1,tmp_event_weight);
               hnum_4region_jes_dn.Fill(2,tmp_event_weight);


              ########################JER
              rrv_mass_lvj.setVal(mass_lvj_jer);
              rrv_mass_j.setVal(tmp_jet_mass_jer);

              if tmp_jet_mass_jer >= self.mj_sideband_lo_min and tmp_jet_mass_jer < self.mj_sideband_lo_max and isPassingCut_jer == 1:
                 rdataset_sb_lo_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData_jer.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

              if tmp_jet_mass_jer >= self.mj_signal_min and tmp_jet_mass_jer < self.mj_signal_max and isPassingCut_jer == 1:
                 rdataset_signal_region_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData_jer.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region_jer.Fill(1,tmp_event_weight);
                 if mass_lvj_jer >=self.mlvj_signal_min and mass_lvj_jer <self.mlvj_signal_max: 
                   hnum_2region_jer.Fill(0,tmp_event_weight);
                   
              if tmp_jet_mass_jer >= self.mj_sideband_hi_min and tmp_jet_mass_jer < self.mj_sideband_hi_max and isPassingCut_jer == 1:
                 rdataset_sb_hi_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj_jer.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

              if isPassingCut_jer == 1: 
               rdataset_mj_jer.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
               rdataset4fit_mj_jer.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
               if tmp_jet_mass_jer >=self.mj_sideband_lo_min and tmp_jet_mass_jer <self.mj_sideband_lo_max: 
                 hnum_4region_jer.Fill(-1,tmp_event_weight );
               if tmp_jet_mass_jer >=self.mj_signal_min and tmp_jet_mass_jer <self.mj_signal_max : 
                 hnum_4region_jer.Fill(0,tmp_event_weight);
               if tmp_jet_mass_jer >=self.mj_sideband_hi_min and tmp_jet_mass_jer <self.mj_sideband_hi_max: 
                 hnum_4region_jer.Fill(1,tmp_event_weight);
               hnum_4region_jer.Fill(2,tmp_event_weight);

              ########################JER_UP
              rrv_mass_lvj.setVal(mass_lvj_jer_up);
              rrv_mass_j.setVal(tmp_jet_mass_jer_up);

              if tmp_jet_mass_jer_up >= self.mj_sideband_lo_min and tmp_jet_mass_jer_up < self.mj_sideband_lo_max and isPassingCut_jer_up == 1:
                 rdataset_sb_lo_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData_jer_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

              if tmp_jet_mass_jer_up >= self.mj_signal_min and tmp_jet_mass_jer_up < self.mj_signal_max and isPassingCut_jer_up == 1:
                 rdataset_signal_region_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData_jer_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer_up.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region_jer_up.Fill(1,tmp_event_weight);
                 if mass_lvj_jer_up >=self.mlvj_signal_min and mass_lvj_jer_up <self.mlvj_signal_max: 
                   hnum_2region_jer_up.Fill(0,tmp_event_weight);
                   
              if tmp_jet_mass_jer_up >= self.mj_sideband_hi_min and tmp_jet_mass_jer_up < self.mj_sideband_hi_max and isPassingCut_jer_up == 1:
                 rdataset_sb_hi_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj_jer_up.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

              if isPassingCut_jer_up == 1: 
               rdataset_mj_jer_up.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
               rdataset4fit_mj_jer_up.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
               if tmp_jet_mass_jer_up >=self.mj_sideband_lo_min and tmp_jet_mass_jer_up <self.mj_sideband_lo_max: 
                 hnum_4region_jer_up.Fill(-1,tmp_event_weight );
               if tmp_jet_mass_jer_up >=self.mj_signal_min and tmp_jet_mass_jer_up <self.mj_signal_max : 
                 hnum_4region_jer_up.Fill(0,tmp_event_weight);
               if tmp_jet_mass_jer_up >=self.mj_sideband_hi_min and tmp_jet_mass_jer_up <self.mj_sideband_hi_max: 
                 hnum_4region_jer_up.Fill(1,tmp_event_weight);
               hnum_4region_jer_up.Fill(2,tmp_event_weight);

              ########################JER_DN
              rrv_mass_lvj.setVal(mass_lvj_jer_dn);
              rrv_mass_j.setVal(tmp_jet_mass_jer_dn);

              if tmp_jet_mass_jer_dn >= self.mj_sideband_lo_min and tmp_jet_mass_jer_dn < self.mj_sideband_lo_max and isPassingCut_jer_dn == 1:
                 rdataset_sb_lo_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_lo_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("sideband");
                 combData_jer_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);

              if tmp_jet_mass_jer_dn >= self.mj_signal_min and tmp_jet_mass_jer_dn < self.mj_signal_max and isPassingCut_jer_dn == 1:
                 rdataset_signal_region_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_signal_region_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );
                 data_category.setLabel("signal_region");
                 combData_jer_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight);
                 combData4fit_jer_dn.add(RooArgSet(rrv_mass_lvj,data_category),tmp_event_weight4fit);
                 hnum_2region_jer_dn.Fill(1,tmp_event_weight);
                 if mass_lvj_jer_dn >=self.mlvj_signal_min and mass_lvj_jer_dn <self.mlvj_signal_max: 
                   hnum_2region_jer_dn.Fill(0,tmp_event_weight);
                   
              if tmp_jet_mass_jer_dn >= self.mj_sideband_hi_min and tmp_jet_mass_jer_dn < self.mj_sideband_hi_max and isPassingCut_jer_dn == 1:
                 rdataset_sb_hi_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight );
                 rdataset4fit_sb_hi_mlvj_jer_dn.add( RooArgSet( rrv_mass_lvj ), tmp_event_weight4fit );

              if isPassingCut_jer_dn == 1: 
               rdataset_mj_jer_dn.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
               rdataset4fit_mj_jer_dn.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
               if tmp_jet_mass_jer_dn >=self.mj_sideband_lo_min and tmp_jet_mass_jer_dn <self.mj_sideband_lo_max: 
                 hnum_4region_jer_dn.Fill(-1,tmp_event_weight );
               if tmp_jet_mass_jer_dn >=self.mj_signal_min and tmp_jet_mass_jer_dn <self.mj_signal_max : 
                 hnum_4region_jer_dn.Fill(0,tmp_event_weight);
               if tmp_jet_mass_jer_dn >=self.mj_sideband_hi_min and tmp_jet_mass_jer_dn <self.mj_sideband_hi_max: 
                 hnum_4region_jer_dn.Fill(1,tmp_event_weight);
               hnum_4region_jer_dn.Fill(2,tmp_event_weight);

        ## scale 4fit dataset in order to have the right luminosity normalization
        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,rdataset_mj.sumEntries()/rdataset4fit_mj.sumEntries());
        rrv_scale_to_lumi.Print();
        rrv_scale_to_lumi_sb_lo = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_lo_mlvj",0);
        if rdataset4fit_sb_lo_mlvj.sumEntries() != 0: 
         rrv_scale_to_lumi_sb_lo.setVal(rdataset_sb_lo_mlvj.sumEntries()/rdataset4fit_sb_lo_mlvj.sumEntries());
         rrv_scale_to_lumi_sb_lo.Print();

        rrv_scale_to_lumi_sb_hi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_sb_hi_mlvj",0);
        if rdataset4fit_sb_hi_mlvj.sumEntries() != 0: 
         rrv_scale_to_lumi_sb_hi.setVal(rdataset_sb_hi_mlvj.sumEntries()/rdataset4fit_sb_hi_mlvj.sumEntries());
         rrv_scale_to_lumi_sb_hi.Print();

        rrv_scale_to_lumi_signal_region = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj","rrv_scale_to_lumi"+label+"_"+self.channel+"_signal_region_mlvj",rdataset_signal_region_mlvj.sumEntries()/rdataset4fit_signal_region_mlvj.sumEntries());
        rrv_scale_to_lumi_signal_region.Print();
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi);
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_sb_lo);
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_sb_hi);
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_signal_region);

        if TString(label).Contains("ggH") or TString(label).Contains("vbfH"):
         print "######### genHMass for BW fit #########";
         getattr(self.workspace4fit_,"import")(rdataset4fit_m_WW_gen);
         rdataset4fit_m_WW_gen.Print();
                                           
        print "########### Nominal Value ###########";

        rrv_number_dataset_signal_region_mlvj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(1));       
        rrv_number_dataset_AllRange_mlvj = RooRealVar("rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"_"+self.channel+"_mlvj",hnum_2region.GetBinContent(2));
               
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj); rrv_number_dataset_signal_region_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj); rrv_number_dataset_AllRange_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj); rdataset_sb_lo_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj); rdataset_signal_region_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj); rdataset_sb_hi_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj); rdataset4fit_sb_lo_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj); rdataset4fit_signal_region_mlvj.Print();
        getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj); rdataset4fit_sb_hi_mlvj.Print();
        getattr(self.workspace4fit_,"import")(combData); combData.Print();
        getattr(self.workspace4fit_,"import")(combData4fit); combData4fit.Print();

        ###################jes_up
        if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
   
         print "########### jes up ###########";
         rrv_scale_to_lumi_jes_up = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jes_up_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jes_up_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jes_up);

         rrv_number_dataset_signal_region_mlvj_jes_up = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jes_up"+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"massvbf_jes_up"+"_"+self.channel+"_mlvj",hnum_2region_jes_up.GetBinContent(1));       
         rrv_number_dataset_AllRange_mlvj_jes_up = RooRealVar("rrv_number_dataset_AllRange"+label+"massvbf_jes_up"+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"massvbf_jes_up"+"_"+self.channel+"_mlvj",hnum_2region_jes_up.GetBinContent(2));
               
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj_jes_up); rrv_number_dataset_signal_region_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj_jes_up); rrv_number_dataset_AllRange_mlvj_jes_up.Print();
                       
         getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj_jes_up); rdataset_sb_lo_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj_jes_up); rdataset_signal_region_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj_jes_up); rdataset_sb_hi_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj_jes_up); rdataset4fit_sb_lo_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj_jes_up); rdataset4fit_signal_region_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj_jes_up); rdataset4fit_sb_hi_mlvj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(combData_jes_up); combData_jes_up.Print();
         getattr(self.workspace4fit_,"import")(combData4fit_jes_up); combData4fit_jes_up.Print();

         ####################jes_dn
         print "########### jes dn ###########";
         rrv_scale_to_lumi_jes_dn = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jes_dn_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jes_dn_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jes_dn);


         rrv_number_dataset_signal_region_mlvj_jes_dn = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jes_dn"+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"massvbf_jes_dn"+"_"+self.channel+"_mlvj",hnum_2region_jes_dn.GetBinContent(1));       
         rrv_number_dataset_AllRange_mlvj_jes_dn = RooRealVar("rrv_number_dataset_AllRange"+label+"massvbf_jes_dn"+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"massvbf_jes_dn"+"_"+self.channel+"_mlvj",hnum_2region_jes_dn.GetBinContent(2));
               
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj_jes_dn); rrv_number_dataset_signal_region_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj_jes_dn); rrv_number_dataset_AllRange_mlvj_jes_dn.Print();
                       
         getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj_jes_dn); rdataset_sb_lo_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj_jes_dn); rdataset_signal_region_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj_jes_dn); rdataset_sb_hi_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj_jes_dn); rdataset4fit_sb_lo_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj_jes_dn); rdataset4fit_signal_region_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj_jes_dn); rdataset4fit_sb_hi_mlvj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(combData_jes_dn); combData_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(combData4fit_jes_dn); combData4fit_jes_dn.Print();

         ####################jer
         print "########### jer ###########";
         rrv_scale_to_lumi_jer = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jer_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jer_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jer.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jer);

         rrv_number_dataset_signal_region_mlvj_jer = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer"+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"massvbf_jer"+"_"+self.channel+"_mlvj",hnum_2region_jer.GetBinContent(1));       
         rrv_number_dataset_AllRange_mlvj_jer = RooRealVar("rrv_number_dataset_AllRange"+label+"massvbf_jer"+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"massvbf_jer"+"_"+self.channel+"_mlvj",hnum_2region_jer.GetBinContent(2));
               
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj_jer); rrv_number_dataset_signal_region_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj_jer); rrv_number_dataset_AllRange_mlvj_jer.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj_jer); rdataset_sb_lo_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj_jer); rdataset_signal_region_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj_jer); rdataset_sb_hi_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj_jer); rdataset4fit_sb_lo_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj_jer); rdataset4fit_signal_region_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj_jer); rdataset4fit_sb_hi_mlvj_jer.Print();
         getattr(self.workspace4fit_,"import")(combData_jer); combData_jer.Print();
         getattr(self.workspace4fit_,"import")(combData4fit_jer); combData4fit_jer.Print();

         ####################jer_up
         print "########### jer up ###########";
         rrv_scale_to_lumi_jer_up = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jer_up_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jer_up_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jer_up);

         rrv_scale_to_lumi_jer_up = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jer_up_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jer_up_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jer_up);

         rrv_number_dataset_signal_region_mlvj_jer_up = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer_up"+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"massvbf_jer_up"+"_"+self.channel+"_mlvj",hnum_2region_jer_up.GetBinContent(1));       
         rrv_number_dataset_AllRange_mlvj_jer_up = RooRealVar("rrv_number_dataset_AllRange"+label+"massvbf_jer_up"+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"massvbf_jer_up"+"_"+self.channel+"_mlvj",hnum_2region_jer_up.GetBinContent(2));
               
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj_jer_up); rrv_number_dataset_signal_region_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj_jer_up); rrv_number_dataset_AllRange_mlvj_jer_up.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj_jer_up); rdataset_sb_lo_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj_jer_up); rdataset_signal_region_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj_jer_up); rdataset_sb_hi_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj_jer_up); rdataset4fit_sb_lo_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj_jer_up); rdataset4fit_signal_region_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj_jer_up); rdataset4fit_sb_hi_mlvj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(combData_jer_up); combData_jer_up.Print();
         getattr(self.workspace4fit_,"import")(combData4fit_jer_up); combData4fit_jer_up.Print();

         ####################jer_dn
         print "########### jer dn ###########";
         rrv_scale_to_lumi_jer_dn = RooRealVar("rrv_scale_to_lumi"+label+"massvbf_jer_dn_"+self.channel,"rrv_scale_to_lumi"+label+"massvbf_jer_dn_"+self.channel,tmp_scale_to_lumi);
         rrv_scale_to_lumi_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_jer_dn);

         rrv_number_dataset_signal_region_mlvj_jer_dn = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer_dn"+"_"+self.channel+"_mlvj","rrv_number_dataset_signal_region"+label+"massvbf_jer_dn"+"_"+self.channel+"_mlvj",hnum_2region_jer_dn.GetBinContent(1));       
         rrv_number_dataset_AllRange_mlvj_jer_dn = RooRealVar("rrv_number_dataset_AllRange"+label+"massvbf_jer_dn"+"_"+self.channel+"_mlvj","rrv_number_dataset_AllRange"+label+"massvbf_jer_dn"+"_"+self.channel+"_mlvj",hnum_2region_jer_dn.GetBinContent(2));
               
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mlvj_jer_dn); rrv_number_dataset_signal_region_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_AllRange_mlvj_jer_dn); rrv_number_dataset_AllRange_mlvj_jer_dn.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_sb_lo_mlvj_jer_dn); rdataset_sb_lo_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset_signal_region_mlvj_jer_dn); rdataset_signal_region_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset_sb_hi_mlvj_jer_dn); rdataset_sb_hi_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_lo_mlvj_jer_dn); rdataset4fit_sb_lo_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_signal_region_mlvj_jer_dn); rdataset4fit_signal_region_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_sb_hi_mlvj_jer_dn); rdataset4fit_sb_hi_mlvj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(combData_jer_dn); combData_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(combData4fit_jer_dn); combData4fit_jer_dn.Print();

        self.file_out.write("\n%s events number in m_lvj from dataset: %s"%(label,rdataset_signal_region_mlvj.sumEntries()))
                             
        #prepare m_j dataset       
        print "########### nominal value ###########";
        rrv_number_dataset_sb_lo_mj = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));        
        rrv_number_dataset_sb_hi_mj = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));        


        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj); rrv_number_dataset_sb_lo_mj.Print();
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj); rrv_number_dataset_signal_region_mj.Print();
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj); rrv_number_dataset_sb_hi_mj.Print();
        
        getattr(self.workspace4fit_,"import")(rdataset_mj); rdataset_mj.Print();
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj); rdataset4fit_mj.Print();

        if label != "_WJets01" and label != "_WJets1" and label !="_data" and not options.skipJetSystematics:
         #####jes_up
         #prepare m_j dataset       
         print "########### jes up ###########";
         rrv_number_dataset_sb_lo_mj_jes_up = RooRealVar("rrv_number_dataset_sb_lo"+label+"massvbf_jes_up"+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"massvbf_jes_up"+"_"+self.channel+"_mj",hnum_4region_jes_up.GetBinContent(1));
         rrv_number_dataset_signal_region_mj_jes_up = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jes_up"+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"massvbf_jes_up"+"_"+self.channel+"_mj",hnum_4region_jes_up.GetBinContent(2));        
         rrv_number_dataset_sb_hi_mj_jes_up = RooRealVar("rrv_number_dataset_sb_hi"+label+"massvbf_jes_up"+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"massvbf_jes_up"+"_"+self.channel+"_mj",hnum_4region_jes_up.GetBinContent(3));        

         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj_jes_up); rrv_number_dataset_sb_lo_mj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj_jes_up); rrv_number_dataset_signal_region_mj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj_jes_up); rrv_number_dataset_sb_hi_mj_jes_up.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_mj_jes_up); rdataset_mj_jes_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_mj_jes_up); rdataset4fit_mj_jes_up.Print();

         #####jes_dn
    
         #prepare m_j dataset       
         print "########### jes dn ###########";
         rrv_number_dataset_sb_lo_mj_jes_dn = RooRealVar("rrv_number_dataset_sb_lo"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj",hnum_4region_jes_dn.GetBinContent(1));
         rrv_number_dataset_signal_region_mj_jes_dn = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj",hnum_4region_jes_dn.GetBinContent(2));        
         rrv_number_dataset_sb_hi_mj_jes_dn = RooRealVar("rrv_number_dataset_sb_hi"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"massvbf_jes_dn"+"_"+self.channel+"_mj",hnum_4region_jes_dn.GetBinContent(3));        

         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj_jes_dn); rrv_number_dataset_sb_lo_mj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj_jes_dn); rrv_number_dataset_signal_region_mj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj_jes_dn); rrv_number_dataset_sb_hi_mj_jes_dn.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_mj_jes_dn); rdataset_mj_jes_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_mj_jes_dn); rdataset4fit_mj_jes_dn.Print();

        #####jer        
         #prepare m_j dataset       
         print "########### jer ###########";
         rrv_number_dataset_sb_lo_mj_jer = RooRealVar("rrv_number_dataset_sb_lo"+label+"massvbf_jer"+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"massvbf_jer"+"_"+self.channel+"_mj",hnum_4region_jer.GetBinContent(1));
         rrv_number_dataset_signal_region_mj_jer = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer"+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"massvbf_jer"+"_"+self.channel+"_mj",hnum_4region_jer.GetBinContent(2));        
         rrv_number_dataset_sb_hi_mj_jer = RooRealVar("rrv_number_dataset_sb_hi"+label+"massvbf_jer"+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"massvbf_jer"+"_"+self.channel+"_mj",hnum_4region_jer.GetBinContent(3));        

         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj_jer); rrv_number_dataset_sb_lo_mj_jer.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj_jer); rrv_number_dataset_signal_region_mj_jer.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj_jer); rrv_number_dataset_sb_hi_mj_jer.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_mj_jer); rdataset_mj_jer.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_mj_jer); rdataset4fit_mj_jer.Print();

        #####jer_up
         #prepare m_j dataset       
         print "########### jer up ###########";
         rrv_number_dataset_sb_lo_mj_jer_up = RooRealVar("rrv_number_dataset_sb_lo"+label+"massvbf_jer_up"+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"massvbf_jer_up"+"_"+self.channel+"_mj",hnum_4region_jer_up.GetBinContent(1));
         rrv_number_dataset_signal_region_mj_jer_up = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer_up"+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"massvbf_jer_up"+"_"+self.channel+"_mj",hnum_4region_jer_up.GetBinContent(2));        
         rrv_number_dataset_sb_hi_mj_jer_up = RooRealVar("rrv_number_dataset_sb_hi"+label+"massvbf_jer_up"+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"massvbf_jer_up"+"_"+self.channel+"_mj",hnum_4region_jer_up.GetBinContent(3));        

         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj_jer_up); rrv_number_dataset_sb_lo_mj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj_jer_up); rrv_number_dataset_signal_region_mj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj_jer_up); rrv_number_dataset_sb_hi_mj_jer_up.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_mj_jer_up); rdataset_mj_jer_up.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_mj_jer_up); rdataset4fit_mj_jer_up.Print();

         #####jer_dn
         #prepare m_j dataset       
         print "########### jer dn ###########";
         rrv_number_dataset_sb_lo_mj_jer_dn = RooRealVar("rrv_number_dataset_sb_lo"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj",hnum_4region_jer_dn.GetBinContent(1));
         rrv_number_dataset_signal_region_mj_jer_dn = RooRealVar("rrv_number_dataset_signal_region"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj",hnum_4region_jer_dn.GetBinContent(2));        
         rrv_number_dataset_sb_hi_mj_jer_dn = RooRealVar("rrv_number_dataset_sb_hi"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"massvbf_jer_dn"+"_"+self.channel+"_mj",hnum_4region_jer_dn.GetBinContent(3));        

         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj_jer_dn); rrv_number_dataset_sb_lo_mj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj_jer_dn); rrv_number_dataset_signal_region_mj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj_jer_dn); rrv_number_dataset_sb_hi_mj_jer_dn.Print();
        
         getattr(self.workspace4fit_,"import")(rdataset_mj_jer_dn); rdataset_mj_jer_dn.Print();
         getattr(self.workspace4fit_,"import")(rdataset4fit_mj_jer_dn); rdataset4fit_mj_jer_dn.Print();

    ##### Prepare the workspace for the limit and to store info to be printed in the datacard
    def prepare_limit(self,mode,isTTbarFloating = 1, isVVFloating = 0, isSTopFloating = 0, isWW_EWKFloating = 0):

        print "####################### prepare_limit for %s method ####################"%(mode);
        self.workspace4fit_.var("rrv_mass_lvj").setBins(int(self.nbins_mlvj*self.narrow_factor));
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_mass_lvj"));

        #################################################################################################################  
        ### whole number of events from the considered signal sample, WJets, VV, TTbar, STop -> counting experiment  ####
        #################################################################################################################
        
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_%s_%s_mlvj"%(self.higgs_sample,self.channel)).clone("rate_%s_for_counting"%(self.higgs_sample) ) )
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_%s_%s_mlvj"%(self.vbfhiggs_sample,self.channel)).clone("rate_%s_for_counting"%(self.vbfhiggs_sample)))

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_WJets0_%s_mlvj"%(self.channel)).clone("rate_WJets_for_counting"))
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_VV_%s_mlvj"%(self.channel)).clone("rate_VV_for_counting"))
        if options.jetBin == "_2jet" : getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_WW_EWK_%s_mlvj"%(self.channel)).clone("rate_WW_EWK_for_counting"))        
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_TTbar_%s_mlvj"%(self.channel)).clone("rate_TTbar_for_counting"))
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_fitting_signal_region_STop_%s_mlvj"%(self.channel)).clone("rate_STop_for_counting"))

        ################################################################################
        ### number of signal, Wjets, VV, TTbar and STop --> unbinned shape analysis  ###
        ################################################################################
        
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).clone("rate_%s_for_unbin"%(self.higgs_sample)));
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).clone("rate_%s_for_unbin"%(self.vbfhiggs_sample)));

        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_WJets0_signal_region%s_%s_mlvj"%(self.mlvj_shape["WJets0"],self.channel)).clone("rate_WJets_for_unbin"));
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_VV_signal_region%s_%s_mlvj"%(self.mlvj_shape["VV"],self.channel)).clone("rate_VV_for_unbin"));
        if options.jetBin == "_2jet" :
            getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_WW_EWK_signal_region%s_%s_mlvj"%(self.mlvj_shape["WW_EWK"],self.channel)).clone("rate_WW_EWK_for_unbin"));        
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_TTbar_signal_region%s_%s_mlvj"%(self.mlvj_shape["TTbar"],self.channel)).clone("rate_TTbar_for_unbin"));
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_STop_signal_region%s_%s_mlvj"%(self.mlvj_shape["STop"],self.channel)).clone("rate_STop_for_unbin"));

        ##############################################################################################################################
        ### Set the error properly -> taking into account lumi, Vtagger and theoretical uncertainty on XS -> for VV, TTbar and STop ##
        ##############################################################################################################################
        
        self.workspace4limit_.var("rate_VV_for_unbin").setError(self.workspace4limit_.var("rate_VV_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty+self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()+self.XS_VV_uncertainty*self.XS_VV_uncertainty));
        if options.jetBin == "_2jet" :
         self.workspace4limit_.var("rate_WW_EWK_for_unbin").setError(self.workspace4limit_.var("rate_WW_EWK_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty+self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()*self.rrv_wtagger_eff_reweight_forV.getError()/self.rrv_wtagger_eff_reweight_forV.getVal()+self.XS_VV_uncertainty*self.XS_VV_uncertainty));        
        self.workspace4limit_.var("rate_STop_for_unbin").setError(self.workspace4limit_.var("rate_STop_for_unbin").getVal()*TMath.Sqrt(self.lumi_uncertainty*self.lumi_uncertainty+self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()+self.XS_STop_uncertainty*self.XS_STop_uncertainty));
        self.workspace4limit_.var("rate_TTbar_for_unbin").setError(self.workspace4limit_.var("rate_TTbar_for_unbin").getVal()*TMath.Sqrt( self.lumi_uncertainty*self.lumi_uncertainty  + self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()*self.rrv_wtagger_eff_reweight_forT.getError()/self.rrv_wtagger_eff_reweight_forT.getVal()));

        ######################################        
        ######### Take the pdf ###############
        ######################################
        
        if mode == "sideband_correction_method1":
            
         getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WJets0_signal_region%s_%s_after_correct_mlvj"%(self.mlvj_shape["WJets0"],self.channel)).clone("WJets_%s"%(self.channel)));
         self.workspace4limit_.allVars().Print();

        if isTTbarFloating:
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signal_region%s_%s_mlvj_Deco_TTbar_signal_region%s_%s_%s_mlvj"%(self.mlvj_shape["TTbar"],self.channel,self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).clone("TTbar_%s"%(self.channel)));
        else :
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_TTbar_signal_region%s_%s_mlvj"%(self.mlvj_shape["TTbar"],self.channel)).clone("TTbar_%s"%(self.channel)));

        if isSTopFloating :
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_STop_signal_region%s_%s_mlvj_Deco_STop_signal_region%s_%s_%s_mlvj"%(self.mlvj_shape["STop"],self.channel,self.mlvj_shape["STop"],self.channel,self.wtagger_label)).clone("STop_%s"%(self.channel)));   
        else :
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_STop_signal_region%s_%s_mlvj"%(self.mlvj_shape["STop"],self.channel)).clone("STop_%s"%(self.channel)));

        if isVVFloating :
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signal_region%s_%s_mlvj_Deco_VV_signal_region%s_%s_%s_mlvj"%(self.mlvj_shape["VV"],self.channel,self.mlvj_shape["VV"],self.channel,self.wtagger_label)).clone("VV_%s"%(self.channel)));
        else:
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_VV_signal_region%s_%s_mlvj"%(self.mlvj_shape["VV"],self.channel)).clone("VV_%s"%(self.channel)));
             
        if isWW_EWKFloating and options.jetBin == "_2jet" :
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WW_EWK_signal_region%s_%s_mlvj_Deco_WW_EWK_signal_region%s_%s_%s_mlvj"%(self.mlvj_shape["WW_EWK"],self.channel,self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label)).clone("WW_EWK_%s"%(self.channel)));             
        elif options.jetBin == "_2jet":
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_WW_EWK_signal_region%s_%s_mlvj"%(self.mlvj_shape["WW_EWK"],self.channel)).clone("WW_EWK_%s"%(self.channel)));         

        ### signal shape
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).clone("ggH_%s"%(self.channel)))
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).clone("qqH_%s"%(self.channel)))

        ##create "fake data" for the limit
        rrv_x = self.workspace4limit_.var("rrv_mass_lvj");
        data_obs         = self.workspace4fit_.data("rdataset_data_signal_region_%s_mlvj"%(self.channel));

        model_pdf_WJets  = self.workspace4limit_.pdf("WJets_%s"%(self.channel));
        model_pdf_VV     = self.workspace4limit_.pdf("VV_%s"%(self.channel));
        model_pdf_TTbar  = self.workspace4limit_.pdf("TTbar_%s"%(self.channel));
        model_pdf_STop   = self.workspace4limit_.pdf("STop_%s"%(self.channel));
        if options.jetBin == "_2jet" :
         model_pdf_WW_EWK = self.workspace4limit_.pdf("WW_EWK_%s"%(self.channel));

        rrv_number_WJets = self.workspace4limit_.var("rate_WJets_for_unbin");
        rrv_number_VV    = self.workspace4limit_.var("rate_VV_for_unbin");
        rrv_number_TTbar = self.workspace4limit_.var("rate_TTbar_for_unbin");
        rrv_number_STop = self.workspace4limit_.var("rate_STop_for_unbin");
        if options.jetBin == "_2jet" :
            rrv_number_WW_EWK = self.workspace4limit_.var("rate_WW_EWK_for_unbin");


        if options.jetBin == "_2jet" :
         rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC","rrv_number_Total_background_MC",rrv_number_WJets.getVal()+rrv_number_VV.getVal()+rrv_number_WW_EWK.getVal()+rrv_number_TTbar.getVal()+rrv_number_STop.getVal());

         rrv_number_Total_background_MC.setError(TMath.Sqrt(rrv_number_WJets.getError()* rrv_number_WJets.getError()+rrv_number_VV.getError()* rrv_number_VV.getError()+rrv_number_WW_EWK.getError()* rrv_number_WW_EWK.getError()+rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+rrv_number_STop.getError() *rrv_number_STop.getError()));

         model_Total_background_MC = RooAddPdf("model_Total_background_MC","model_Total_background_MC",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_WW_EWK,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_WW_EWK,rrv_number_TTbar,rrv_number_STop));
        else:
         rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC","rrv_number_Total_background_MC",rrv_number_WJets.getVal()+rrv_number_VV.getVal()+rrv_number_TTbar.getVal()+rrv_number_STop.getVal());

         rrv_number_Total_background_MC.setError(TMath.Sqrt(rrv_number_WJets.getError()* rrv_number_WJets.getError()+rrv_number_VV.getError()* rrv_number_VV.getError()+rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+rrv_number_STop.getError() *rrv_number_STop.getError()));

         model_Total_background_MC = RooAddPdf("model_Total_background_MC","model_Total_background_MC",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_TTbar,rrv_number_STop));
            

        new_data_obs = model_Total_background_MC.generate(RooArgSet(rrv_x),int(data_obs.sumEntries()) );
        new_data_obs.Print();

        if options.pseudodata == 0:
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.data("rdataset_data_signal_region_%s_mlvj"%(self.channel)).Clone("data_obs_%s"%(self.channel)))
        else:
          getattr(self.workspace4limit_,"import")(new_data_obs.Clone("data_obs_%s"%(self.channel)))

        ###############
        ### FIX PDF ###
        ###############
        
        fix_Pdf(self.workspace4limit_.pdf("TTbar_%s"%(self.channel)), RooArgSet(rrv_x)); 
        fix_Pdf(self.workspace4limit_.pdf("STop_%s"%(self.channel)),  RooArgSet(rrv_x)); 
        fix_Pdf(self.workspace4limit_.pdf("VV_%s"%(self.channel)),    RooArgSet(rrv_x));
        if options.jetBin == "_2jet": fix_Pdf(self.workspace4limit_.pdf("WW_EWK_%s"%(self.channel)),RooArgSet(rrv_x));         
        fix_Pdf(self.workspace4limit_.pdf("WJets_%s"%(self.channel)), RooArgSet(rrv_x)); 
        fix_Pdf(self.workspace4limit_.pdf("ggH_%s"%(self.channel)),   RooArgSet(rrv_x));
        fix_Pdf(self.workspace4limit_.pdf("qqH_%s"%(self.channel)),   RooArgSet(rrv_x));

        ### print the workspace 4 limit 
        print " ############## Workspace for limit ";
        parameters_workspace = self.workspace4limit_.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next()
        while (param):
         param.Print();
         param = par.Next()
        self.workspace4limit_.Print()

        params_list=[];

        ###################################################  
        ### main modality for the alpha function method ###
        ###################################################
        
        if mode == "sideband_correction_method1":

            if self.MODEL_4_mlvj == "ErfExp_v1" or self.MODEL_4_mlvj == "ErfPow_v1" or self.MODEL_4_mlvj == "2Exp" :
                ### uncertainty inflation on the Wjets shape from fitting data in sb_lo
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                ### Add to the parameter list
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                ### Do the same for alpha paramter
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                ### Add to the parameter list
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))

                ### Do the same for the TTbar
                if isTTbarFloating !=0:
                  self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                  self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                  self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);

                  params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                  params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                  params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));

            if self.MODEL_4_mlvj == "ErfPow2_v1" or self.MODEL_4_mlvj == "ErfPowExp_v1" :
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig6"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig7"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig6"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig7"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))


                if isTTbarFloating !=0 :

                 self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                 self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                 self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                 self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);

                 params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                 params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                 params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                 params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));

            if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                if isTTbarFloating !=0 :
                  self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                  params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));

            if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);
                self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_WJets0);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);
                self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)).setError(self.shape_para_error_alpha);

                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))
                params_list.append(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)))

                ### TTbar use exp
                if isTTbarFloating !=0:
                     self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)).setError(self.shape_para_error_TTbar);
                     params_list.append(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));


        ### signal parameters
        if options.jetBin == "_2jet":            
          systematic = SystematicUncertaintyHiggs_2jetBin() ;
        else:
          systematic = SystematicUncertaintyHiggs_01jetBin() ;

  
        if TString(self.higgs_sample).Contains("600") :

           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_600);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_ggH_600);
                                                                                                                   
           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_600);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_600);
                                                                                                                                                                  
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_600);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_600);
                                     
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_600);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_600);

        elif TString(self.higgs_sample).Contains("700") :

           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_700);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_ggH_700);
                                                                                                                   
           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_700);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_700);
                                                                                                                                                                  

           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_700);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_700);
                                     
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_700);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_700);

        elif TString(self.higgs_sample).Contains("800") :

           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_800);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_ggH_800);
                                                                                                                   
           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_800);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_800);
                                                                                                                                                                  

           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_800);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_800);
                                     
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_800);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_800);

        elif TString(self.higgs_sample).Contains("900") :

           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_900);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_ggH_900);
                                                                                                                   
           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_900);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_900);
                                                                                                                                                                  

           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_900);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_900);
                                     
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_900);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_900);

        elif TString(self.higgs_sample).Contains("1000") :

           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_1000);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_ggH_1000);
                                                                                                                   
           self.workspace4limit_.var("rrv_mean_shift_scale_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_1000);
           self.workspace4limit_.var("rrv_mean_shift_scale_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_1000);
                                                                                                                                                                  

           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_1000);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_1000);
                                     
           self.workspace4limit_.var("rrv_sigma_shift_jes_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_1000);
           self.workspace4limit_.var("rrv_sigma_shift_jer_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_1000);

        self.workspace4limit_.var("CMS_scale_j").setError(1);
        self.workspace4limit_.var("CMS_res_j").setError(1);

        params_list.append(self.workspace4limit_.var("CMS_scale_j"));
        params_list.append(self.workspace4limit_.var("CMS_res_j"));

        ### print the datacard                                       
        self.print_limit_datacard("unbin", "ggHvbfH",params_list);

        if mode == "sideband_correction_method1":

          if self.MODEL_4_mlvj == "ErfExp_v1" or self.MODEL_4_mlvj == "ErfPow_v1" or self.MODEL_4_mlvj == "2Exp" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)));

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                if isTTbarFloating!=0:
                   self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)) );
                   self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)) );
                   self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)) );

          if self.MODEL_4_mlvj=="ErfPow2_v1" or self.MODEL_4_mlvj=="ErfPowExp_v1" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig6"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig7"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig4"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig5"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig6"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig7"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                if isTTbarFloating!=0:
                     self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                     self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                     self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));
                     self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));

          if self.MODEL_4_mlvj=="Exp" or self.MODEL_4_mlvj=="Pow" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                if isTTbarFloating!=0:
                  self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));


          if self.MODEL_4_mlvj=="ExpN" or self.MODEL_4_mlvj=="ExpTail" or self.MODEL_4_mlvj=="Pow2" :

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sb_lo%s_from_fitting_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig1"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig2"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );
                self.FloatingParams_wjet.add(self.workspace4limit_.var("Deco_WJets0_sim_%s_%s_%s_mlvj_eig3"%(self.mlvj_shape["WJets0"],self.channel,self.wtagger_label)) );

                if isTTbarFloating!=0:
                     self.FloatingParams.add(self.workspace4limit_.var("Deco_TTbar_signal_region%s_%s_%s_mlvj_eig0"%(self.mlvj_shape["TTbar"],self.channel,self.wtagger_label)));

 
        self.FloatingParams.add(self.workspace4limit_.var("CMS_scale_j"));
        self.FloatingParams.add(self.workspace4limit_.var("CMS_res_j"));

        ### Add the floating list to the combiner --> the pdf which are not fixed are floating by default
        getattr(self.workspace4limit_,"import")(self.FloatingParams);

        ### Save the workspace
        self.save_workspace_to_file();


    #### Method used in order to save the workspace in a output root file
    def save_workspace_to_file(self):
        self.workspace4limit_.writeToFile(self.file_rlt_root);
        self.file_out.close()


    #### Method used to print the general format of the datacard for both counting and unbinned analysis         
    def print_limit_datacard(self, mode, signalchannel, params_list=[] ): #mode:unbin or counting

      print "############## print_limit_datacard for %s %s ################"%(mode,signalchannel)

      if (mode == "unbin" or mode == "counting"): 
         print "print_limit_datacard use wrong mode: %s"%(mode);#raw_input("ENTER");

         ### open the datacard
         datacard_out = open(getattr(self,"file_datacard_%s_%s"%(mode, signalchannel)),"w");
         datacard_out.write( "imax 1" )
         if signalchannel == "ggH" or signalchannel == "vbfH":
            datacard_out.write( "\njmax *" )
         elif signalchannel=="ggHvbfH":
            datacard_out.write( "\njmax *" )
         else:
            raw_input("Wrong signal channel, please check!!");
            
         datacard_out.write( "\nkmax *" )
         datacard_out.write( "\n--------------- ")


         datacard_out.write( "\nbin CMS_%s"%(self.channel));
         if mode == "unbin":
            datacard_out.write("\nobservation %0.2f "%(self.workspace4limit_.data("data_obs_%s"%(self.channel)).sumEntries()));
         if mode == "counting":
            datacard_out.write("\nobservation %0.2f "%(self.workspace4limit_.var("observation_for_counting").getVal()));

         datacard_out.write( "\n------------------------------" );

         if mode == "unbin":
          fnOnly = ntpath.basename(self.file_rlt_root)
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signal_region%s_%s_mlvj"%(self.higgs_sample,self.mlvj_shape["ggH"],self.channel)).clone(self.higgs_sample+"_%s"%(self.channel)))
          getattr(self.workspace4limit_,"import")(self.workspace4fit_.pdf("model_pdf_%s_signal_region%s_%s_mlvj"%(self.vbfhiggs_sample,self.mlvj_shape["vbfH"],self.channel)).clone(self.vbfhiggs_sample+"_%s"%(self.channel)))

          if signalchannel == "ggH":
              datacard_out.write("\nshapes ggH CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          elif signalchannel == "vbfH":
              datacard_out.write("\nshapes qqH CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          elif signalchannel == "ggHvbfH":
              datacard_out.write("\nshapes ggH CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
              datacard_out.write("\nshapes qqH CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          
          datacard_out.write("\nshapes WJets CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          datacard_out.write("\nshapes TTbar CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          datacard_out.write("\nshapes STop CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          datacard_out.write("\nshapes VV CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          if options.jetBin == "_2jet" : datacard_out.write("\nshapes WW_EWK CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(),self.channel));
          datacard_out.write("\nshapes data_obs CMS_%s %s %s:$PROCESS_%s"%(self.channel,fnOnly,self.workspace4limit_.GetName(), 
self.channel));
          datacard_out.write( "\n--------------- ")
          

         if signalchannel == "ggH":
             
            if options.jetBin == "_2jet" : 
             datacard_out.write( "\nbin                CMS_%s    CMS_%s   CMS_%s   CMS_%s  CMS_%s   CMS_%s"%(self.channel,self.channel,self.channel,self.channel,self.channel,self.channel));
             datacard_out.write( "\nprocess            ggH        WJets    TTbar    STop    VV     WW_EWK");
             datacard_out.write( "\nprocess            -1          1        2        3       4       5  " );

             if mode == "unbin":
                datacard_out.write( "\nrate          %0.2f          %0.2f   %0.2f    %0.2f    %0.2f     %0.2f "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_unbin").getVal()));

             elif mode == "counting":
                datacard_out.write( "\nrate          %0.2f          %0.2f   %0.2f    %0.2f    %0.2f    %0.2f "%(self.workspace4limit_.var("rate_%s_for_counting"%(self.higgs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_counting").getVal()));
                
             datacard_out.write( "\n-------------------------------- " );

             datacard_out.write( "\nQCDscale_ggH0in     lnN   %0.3f     -             -        -       -       - "%(1.+self.QCDscale_ggH0in));

             datacard_out.write( "\nQCDscale_ggH2in     lnN   %0.3f     -             -        -       -       - "%(1.+self.QCDscale_ggH2in));   

             datacard_out.write( "\npdf_gg              lnN   %0.3f     -             -        -       -       - "%(1.+self.pdf_gg));

             datacard_out.write( "\nQCDscale_ggH_ACCEPT lnN   %0.3f     -             -        -       -       - "%(1.+self.QCDscale_ggH_ACCEPT));

             datacard_out.write( "\nintf_ggH            lnN   %0.3f     -             -        -       -       - "%(1.+self.interference_ggH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop     lnN   -         -             -      %0.3f     -        - "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV       lnN   -         -             -        -     %0.3f      - "%(1+self.XS_VV_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_WW_EWK   lnN   -         -             -        -       -     %0.3f "%(1+self.XS_WW_EWK_uncertainty));

            else:

             datacard_out.write( "\nbin                CMS_%s    CMS_%s   CMS_%s    CMS_%s   CMS_%s "%(self.channel,self.channel,self.channel,self.channel,self.channel));
             datacard_out.write( "\nprocess            ggH        WJets    TTbar     STop      VV   ");
             datacard_out.write( "\nprocess            -1          1        2        3         4    ");

             if mode == "unbin":
                datacard_out.write( "\nrate          %0.2f          %0.2f   %0.2f    %0.2f    %0.2f  "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal()));

             elif mode == "counting":
                datacard_out.write( "\nrate          %0.2f          %0.2f   %0.2f    %0.2f    %0.2f"  %(self.workspace4limit_.var("rate_%s_for_counting"%(self.higgs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal()));
                
             datacard_out.write( "\n-------------------------------- " ); 

             datacard_out.write( "\nQCDscale_ggH0in     lnN   %0.3f       -      -        -       -     "%(1.+self.QCDscale_ggH0in));

             datacard_out.write( "\nQCDscale_ggH2in     lnN   %0.3f       -      -        -       -     "%(1.+self.QCDscale_ggH2in));   

             datacard_out.write( "\npdf_gg              lnN   %0.3f       -      -        -       -     "%(1.+self.pdf_gg));

             datacard_out.write( "\nQCDscale_ggH_ACCEPT lnN   %0.3f       -      -        -       -     "%(1.+self.QCDscale_ggH_ACCEPT));

             datacard_out.write( "\nintf_ggH            lnN   %0.3f       -      -        -       -     "%(1.+self.interference_ggH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop     lnN     -         -      -      %0.3f     -     "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV       lnN     -         -      -        -     %0.3f   "%(1+self.XS_VV_uncertainty));

 
         elif signalchannel=="vbfH":
            if options.jetBin == "_2jet":  
             datacard_out.write( "\nbin                CMS_%s    CMS_%s   CMS_%s   CMS_%s   CMS_%s   CMS_%s"%(self.channel,self.channel,self.channel,self.channel,self.channel,self.channel));
             datacard_out.write( "\nprocess            qqH        WJets    TTbar    STop     VV      WW_EWK");
             datacard_out.write( "\nprocess            -1            1      2        3      4        5     ");

             if mode == "unbin":
                datacard_out.write( "\nrate            %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    %0.2f"%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_unbin").getVal()))

             elif mode == "counting":
                datacard_out.write( "\nrate            %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    %0.2f"%(self.workspace4limit_.var("rate_%s_for_counting"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_counting").getVal()))

             datacard_out.write( "\n-------------------------------- " );
            
             datacard_out.write( "\nQCDscale_qqH         lnN       %0.3f         -        -         -       -       -  "%(1.+self.QCDscale_qqH));

             datacard_out.write( "\npdf_qqbar            lnN       %0.3f         -        -         -       -       -  "%(1.+self.pdf_qqbar));

             datacard_out.write( "\nQCDscale_qqH_ACCEPT  lnN       %0.3f         -        -         -       -       -  "%(1.+self.QCDscale_qqH_ACCEPT));

             datacard_out.write( "\nintf_vbfH            lnN       %0.3f         -        -         -       -       -  "%(1.+self.interference_vbfH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop      lnN       -             -        -        %0.3f    -       -  "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV        lnN       -             -        -         -     %0.3f     -  "%(1+self.XS_VV_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_WW_EWK    lnN       -             -        -         -       -     %0.3f  "%(1+self.XS_WW_EWK_uncertainty));


            else:  
             datacard_out.write( "\nbin                CMS_%s    CMS_%s   CMS_%s   CMS_%s  CMS_%s "%(self.channel,self.channel,self.channel,self.channel,self.channel));
             datacard_out.write( "\nprocess            qqH       WJets    TTbar    STop    VV     ");
             datacard_out.write( "\nprocess            -1            1      2       3       4     ");

             if mode == "unbin":
                datacard_out.write( "\nrate            %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal()))

             elif mode == "counting":
                datacard_out.write( "\nrate            %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    "%(self.workspace4limit_.var("rate_%s_for_counting"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal()))

             datacard_out.write( "\n-------------------------------- " );
            
             datacard_out.write( "\nQCDscale_qqH        lnN    %0.3f         -        -       -       -     "%(1.+self.QCDscale_qqH) );

             datacard_out.write( "\npdf_qqbar           lnN    %0.3f         -        -       -       -     "%(1.+self.pdf_qqbar));

             datacard_out.write( "\nQCDscale_qqH_ACCEPT lnN    %0.3f         -        -       -       -     "%(1.+self.QCDscale_qqH_ACCEPT));

             datacard_out.write( "\nintf_vbfH           lnN    %0.3f         -        -       -       -     "%(1.+self.interference_vbfH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop     lnN      -           -        -     %0.3f     -     "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV       lnN      -           -        -       -      %0.3f  "%(1+self.XS_VV_uncertainty));
 

         elif signalchannel == "ggHvbfH":
             
            if options.jetBin == "_2jet":
                
             datacard_out.write( "\nbin                CMS_%s    CMS_%s     CMS_%s   CMS_%s   CMS_%s   CMS_%s    CMS_%s"%(self.channel,self.channel,self.channel,self.channel,self.channel,self.channel,self.channel));            
             datacard_out.write( "\nprocess            ggH        qqH        WJets     TTbar    STop       VV     WW_EWK ");
             datacard_out.write( "\nprocess            -1         0           1          2       3         4     5       ");

             if mode == "unbin":
                datacard_out.write( "\nrate               %0.2f    %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    %0.2f "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_%s_for_unbin"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_unbin").getVal()))

             elif mode == "counting":
                datacard_out.write( "\nrate               %0.2f    %0.2f         %0.2f   %0.2f    %0.2f    %0.2f    %0.2f"%(self.workspace4limit_.var("rate_%s_for_counting"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_%s_for_counting"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal(), self.workspace4limit_.var("rate_WW_EWK_for_counting").getVal()  ) )


             datacard_out.write( "\n-------------------------------- " );

             datacard_out.write( "\nQCDscale_ggH0in     lnN   %0.3f       -        -        -       -       -       -   "%(1.+self.QCDscale_ggH0in));

             datacard_out.write( "\nQCDscale_ggH2in      lnN   %0.3f       -        -        -       -       -       -   "%(1.+self.QCDscale_ggH2in));   

             datacard_out.write( "\npdf_gg               lnN   %0.3f       -        -        -       -       -       -   "%(1.+self.pdf_gg));

             datacard_out.write( "\nQCDscale_ggH_ACCEPT  lnN   %0.3f       -        -        -       -       -       -   "%(1.+self.QCDscale_ggH_ACCEPT));

             datacard_out.write( "\nintf_ggH             lnN   %0.3f       -        -        -       -       -       -   "%(1.+self.interference_ggH_uncertainty));
            
             datacard_out.write( "\nQCDscale_qqH         lnN    -         %0.3f     -        -       -       -       -   "%(1.+self.QCDscale_qqH));

             datacard_out.write( "\npdf_qqbar            lnN    -         %0.3f     -        -       -       -       -   "%(1.+self.pdf_qqbar));

             datacard_out.write( "\nQCDscale_qqH_ACCEPT  lnN    -         %0.3f     -        -       -       -       -   "%(1.+self.QCDscale_qqH_ACCEPT));

             datacard_out.write( "\nintf_vbfH            lnN    -         %0.3f     -        -       -       -       -   "%(1.+self.interference_vbfH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop      lnN    -          -        -        -      %0.3f    -       -   "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV        lnN    -          -        -        -       -      %0.3f    -   "%(1+self.XS_VV_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_WW_EWK    lnN    -          -        -        -       -       -     %0.3f "%(1+self.XS_WW_EWK_uncertainty));

            else:

             datacard_out.write( "\nbin                CMS_%s    CMS_%s    CMS_%s   CMS_%s   CMS_%s  CMS_%s  "%(self.channel,self.channel,self.channel,self.channel,self.channel,self.channel));            
             datacard_out.write( "\nprocess            ggH       qqH       WJets    TTbar    STop    VV   ");
             datacard_out.write( "\nprocess            -1        0          1        2        3      4    ");

             if mode == "unbin":
                datacard_out.write( "\nrate               %0.2f    %0.2f         %0.2f   %0.2f    %0.2f    %0.2f  "%(self.workspace4limit_.var("rate_%s_for_unbin"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_%s_for_unbin"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_unbin").getVal(), self.workspace4limit_.var("rate_TTbar_for_unbin").getVal(), self.workspace4limit_.var("rate_STop_for_unbin").getVal(), self.workspace4limit_.var("rate_VV_for_unbin").getVal()));

             elif mode == "counting":
                datacard_out.write( "\nrate               %0.2f    %0.2f         %0.2f   %0.2f    %0.2f    %0.2f  "%(self.workspace4limit_.var("rate_%s_for_counting"%(self.higgs_sample)).getVal(),self.workspace4limit_.var("rate_%s_for_counting"%(self.vbfhiggs_sample)).getVal(), self.workspace4limit_.var("rate_WJets_for_counting").getVal(), self.workspace4limit_.var("rate_TTbar_for_counting").getVal(), self.workspace4limit_.var("rate_STop_for_counting").getVal(), self.workspace4limit_.var("rate_VV_for_counting").getVal()));


             datacard_out.write( "\n-------------------------------- " );

             datacard_out.write( "\nQCDscale_ggH0in    lnN   %0.3f     -       -       -       -       -       "%(1.+self.QCDscale_ggH0in));

             datacard_out.write( "\nQCDscale_ggH2in     lnN   %0.3f     -       -       -       -       -       "%(1.+self.QCDscale_ggH2in));   

             datacard_out.write( "\npdf_gg              lnN   %0.3f     -       -       -       -       -       "%(1.+self.pdf_gg));

             datacard_out.write( "\nQCDscale_ggH_ACCEPT lnN   %0.3f     -       -       -       -       -       "%(1.+self.QCDscale_ggH_ACCEPT));

             datacard_out.write( "\nintf_ggH            lnN   %0.3f     -       -       -       -       -       "%(1.+self.interference_ggH_uncertainty));
            
             datacard_out.write( "\nQCDscale_qqH        lnN   -         %0.3f   -       -       -       -       "%(1.+self.QCDscale_qqH));

             datacard_out.write( "\npdf_qqbar           lnN   -         %0.3f   -       -       -       -       "%(1.+self.pdf_qqbar));

             datacard_out.write( "\nQCDscale_qqH_ACCEPT lnN   -         %0.3f   -       -       -       -       "%(1.+self.QCDscale_qqH_ACCEPT));

             datacard_out.write( "\nintf_vbfH           lnN   -         %0.3f   -       -       -       -       "%(1.+self.interference_vbfH_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_STop     lnN   -         -       -       -       %0.3f   -       "%(1+self.XS_STop_uncertainty));

             datacard_out.write( "\nCMS_hwwlvj_VV       lnN   -         -       -       -       -       %0.3f   "%(1+self.XS_VV_uncertainty));

            
         if options.jetBin == "_2jet":
             
          datacard_out.write( "\nlumi_8TeV       lnN       %0.3f     %0.3f         -        -  %0.3f   %0.3f    %0.3f"%(1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty) );

          datacard_out.write( "\nCMS_trigger_m  lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f    %0.3f"%(1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty, 1+self.lep_trigger_uncertainty ) );

          datacard_out.write( "\nCMS_trigger_e  lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f    %0.3f"%(1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty, 1+self.lep_trigger_uncertainty ) );

          datacard_out.write( "\nCMS_eff_m      lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f    %0.3f"%(1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty ) );

          datacard_out.write( "\nCMS_eff_e      lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f    %0.3f"%(1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty, 1+self.lep_eff_uncertainty ) );

          datacard_out.write( "\nCMS_TTbar_norm_2jet lnN       -          -            -        %0.3f   -   -            -"%(1+self.rrv_wtagger_eff_reweight_forT.getError()));

          datacard_out.write( "\nCMS_wtagger     lnN       %0.3f     %0.3f         -        -       -       %0.3f    %0.3f"%(1+self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError() ) );
            
          ### nousiance for the bkg
          ### WJets normaliztion due to data fit and alternate modellization
          if self.number_WJets_insideband >0:
            datacard_out.write( "\nWjet_Norm%s lnN     %0.3f       -         %0.3f      -       -        -        -  "%(options.jetBin,self.number_WJets_insideband, getattr(self,"datadriven_alpha_WJets_%s"%(mode))));
          else:              
            datacard_out.write( "\nWjet_Norm%s lnN     -           -         %0.3f      -       -        -        -  "%(options.jetBin, 1+ self.workspace4limit_.var("rate_WJets_for_unbin").getError()/self.workspace4limit_.var("rate_WJets_for_unbin").getVal() ) );


          ## jet mass systematic scaling up and down vbf jets detajj, mjj, and pt selection effect
          if self.ggH_normalization_uncertainty_from_jet_scale!=0 and self.vbf_normalization_uncertainty_from_jet_scale!=0 and self.WJets_normalization_uncertainty_from_jet_scale!=0 and self.TTbar_normalization_uncertainty_from_jet_scale!=0 and self.STop_normalization_uncertainty_from_jet_scale!=0 and self.VV_normalization_uncertainty_from_jet_scale!=0 : 

           datacard_out.write( "\nCMS_scale_j lnN   %0.3f     %0.3f     %0.3f/%0.3f    %0.3f/%0.3f   %0.3f/%0.3f   %0.3f/%0.3f    %0.3f/%0.3f"%(1+self.ggH_normalization_uncertainty_from_jet_scale, 1+self.vbf_normalization_uncertainty_from_jet_scale, 1-self.WJets_normalization_uncertainty_from_jet_scale, 1+self.WJets_normalization_uncertainty_from_jet_scale, 1+self.TTbar_normalization_uncertainty_from_jet_scale, 1-self.TTbar_normalization_uncertainty_from_jet_scale, 1+self.STop_normalization_uncertainty_from_jet_scale, 1-self.STop_normalization_uncertainty_from_jet_scale, 1+self.VV_normalization_uncertainty_from_jet_scale, 1-self.VV_normalization_uncertainty_from_jet_scale, 1+self.WW_EWK_normalization_uncertainty_from_jet_scale, 1-self.WW_EWK_normalization_uncertainty_from_jet_scale ) )        

          if self.ggH_normalization_uncertainty_from_jet_res!=0 and self.vbf_normalization_uncertainty_from_jet_res!=0 and self.WJets_normalization_uncertainty_from_jet_res!=0 and self.TTbar_normalization_uncertainty_from_jet_res!=0 and self.STop_normalization_uncertainty_from_jet_res!=0 and self.VV_normalization_uncertainty_from_jet_res!=0 :
             
           datacard_out.write( "\nCMS_res_j lnN   %0.3f     %0.3f     %0.3f/%0.3f    %0.3f/%0.3f   %0.3f/%0.3f   %0.3f/%0.3f    %0.3f/%0.3f"%(1+self.ggH_normalization_uncertainty_from_jet_res, 1+self.vbf_normalization_uncertainty_from_jet_res, 1-self.WJets_normalization_uncertainty_from_jet_res, 1+self.WJets_normalization_uncertainty_from_jet_res, 1+self.TTbar_normalization_uncertainty_from_jet_res, 1-self.TTbar_normalization_uncertainty_from_jet_res, 1+self.STop_normalization_uncertainty_from_jet_res, 1-self.STop_normalization_uncertainty_from_jet_res, 1+self.VV_normalization_uncertainty_from_jet_res, 1-self.VV_normalization_uncertainty_from_jet_res, 1+self.WW_EWK_normalization_uncertainty_from_jet_res, 1-self.WW_EWK_normalization_uncertainty_from_jet_res ) )        

          if self.ggH_normalization_uncertainty_from_lep_scale!=0 and self.vbf_normalization_uncertainty_from_lep_scale!=0 and self.WJets_normalization_uncertainty_from_lep_scale!=0 and self.TTbar_normalization_uncertainty_from_lep_scale!=0 and self.STop_normalization_uncertainty_from_lep_scale!=0 and self.VV_normalization_uncertainty_from_lep_scale!=0 : 

           datacard_out.write( "\nCMS_scale_e lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    %0.3f"%(self.ggH_normalization_uncertainty_from_lep_scale, self.vbf_normalization_uncertainty_from_lep_scale, self.WJets_normalization_uncertainty_from_lep_scale, self.TTbar_normalization_uncertainty_from_lep_scale, self.STop_normalization_uncertainty_from_lep_scale, self.VV_normalization_uncertainty_from_lep_scale, self.WW_EWK_normalization_uncertainty_from_lep_scale ) )        

           datacard_out.write( "\nCMS_scale_m lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    %0.3f"%(self.ggH_normalization_uncertainty_from_lep_scale, self.vbf_normalization_uncertainty_from_lep_scale, self.WJets_normalization_uncertainty_from_lep_scale, self.TTbar_normalization_uncertainty_from_lep_scale, self.STop_normalization_uncertainty_from_lep_scale, self.VV_normalization_uncertainty_from_lep_scale, self.WW_EWK_normalization_uncertainty_from_lep_scale ) )        


          if self.ggH_normalization_uncertainty_from_lep_res!=0 and self.vbf_normalization_uncertainty_from_lep_res!=0 and self.WJets_normalization_uncertainty_from_lep_res!=0 and self.TTbar_normalization_uncertainty_from_lep_res!=0 and self.STop_normalization_uncertainty_from_lep_res!=0 and self.VV_normalization_uncertainty_from_lep_res!=0 : 

           datacard_out.write( "\nCMS_res_m lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    %0.3f"%(self.ggH_normalization_uncertainty_from_lep_res, self.vbf_normalization_uncertainty_from_lep_res, self.WJets_normalization_uncertainty_from_lep_res, self.TTbar_normalization_uncertainty_from_lep_res, self.STop_normalization_uncertainty_from_lep_res, self.VV_normalization_uncertainty_from_lep_res, self.WW_EWK_normalization_uncertainty_from_lep_res ) )

           datacard_out.write( "\nCMS_res_e lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    %0.3f"%(self.ggH_normalization_uncertainty_from_lep_res, self.vbf_normalization_uncertainty_from_lep_res, self.WJets_normalization_uncertainty_from_lep_res, self.TTbar_normalization_uncertainty_from_lep_res, self.STop_normalization_uncertainty_from_lep_res, self.VV_normalization_uncertainty_from_lep_res, self.WW_EWK_normalization_uncertainty_from_lep_res ) )


          if self.ggH_normalization_uncertainty_from_btag!=0 and self.vbf_normalization_uncertainty_from_btag!=0 and self.WJets_normalization_uncertainty_from_btag!=0 and self.TTbar_normalization_uncertainty_from_btag!=0 and self.STop_normalization_uncertainty_from_btag!=0 and self.VV_normalization_uncertainty_from_btag!=0 : 

           datacard_out.write( "\nCMS_eff_b lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    %0.3f"%(self.ggH_normalization_uncertainty_from_btag, self.vbf_normalization_uncertainty_from_btag, self.WJets_normalization_uncertainty_from_btag, self.TTbar_normalization_uncertainty_from_btag, self.STop_normalization_uncertainty_from_btag, self.VV_normalization_uncertainty_from_btag, self.WW_EWK_normalization_uncertainty_from_btag ) )                  

         else:
 
          datacard_out.write( "\nlumi_8TeV       lnN       %0.3f     %0.3f         -        -   %0.3f   %0.3f  "%(1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty,1.+self.lumi_uncertainty) )

          if self.channel == "el" :
           datacard_out.write( "\nCMS_trigger_e  lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f   "%(1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty) );
          elif self.channel == "mu" :
           datacard_out.write( "\nCMS_trigger_m  lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f   "%(1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty,1+self.lep_trigger_uncertainty) );

          if self.channel == "el" :
           datacard_out.write( "\nCMS_eff_e      lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f   "%(1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty) );
          elif self.channel == "mu" :
           datacard_out.write( "\nCMS_eff_m      lnN       %0.3f     %0.3f         -        %0.3f   %0.3f   %0.3f   "%(1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty,1+self.lep_eff_uncertainty) );
   
          if self.channel == "el" :
           datacard_out.write( "\nCMS_TTbar_norm_e lnN         -         -           -        %0.3f   -   -      "%(1+self.rrv_wtagger_eff_reweight_forT.getError()));
          elif self.channel == "mu":
           datacard_out.write( "\nCMS_TTbar_norm_m lnN         -         -           -        %0.3f   -   -      "%(1+self.rrv_wtagger_eff_reweight_forT.getError()));
              
          datacard_out.write( "\nCMS_wtagger     lnN       %0.3f     %0.3f         -        -       -       %0.3f    "%(1+self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError(), 1+self.rrv_wtagger_eff_reweight_forV.getError()) );
            
          ### nousiance for the bkg
          ### WJets normaliztion due to data fit and alternate modellization
          if self.channel == "mu" : channel = "m" ;
          elif self.channel == "el" : channel = "e" ;
          
          if self.number_WJets_insideband >0:
            datacard_out.write( "\nWjet_Norm_%s   lnN    %0.3f      -    %0.3f     -      -        -      "%(channel,self.number_WJets_insideband, getattr(self,"datadriven_alpha_WJets_%s"%(mode))));
          else:
            datacard_out.write( "\nWjet_Norm_%s   lnN     -         -    %0.3f     -      -        -      "%(channel,1+ self.workspace4limit_.var("rate_WJets_for_unbin").getError()/self.workspace4limit_.var("rate_WJets_for_unbin").getVal() ) );


          ## jet mass systematic scaling up and down vbf jets detajj, mjj, and pt selection effect
          if self.ggH_normalization_uncertainty_from_jet_scale!=0 and self.vbf_normalization_uncertainty_from_jet_scale!=0 and self.WJets_normalization_uncertainty_from_jet_scale!=0 and self.TTbar_normalization_uncertainty_from_jet_scale!=0 and self.STop_normalization_uncertainty_from_jet_scale!=0 and self.VV_normalization_uncertainty_from_jet_scale!=0 : 

           datacard_out.write( "\nCMS_scale_j lnN   %0.3f     %0.3f     %0.3f/%0.3f    %0.3f/%0.3f   %0.3f/%0.3f   %0.3f/%0.3f    "%(1+self.ggH_normalization_uncertainty_from_jet_scale, 1+self.vbf_normalization_uncertainty_from_jet_scale, 1-self.WJets_normalization_uncertainty_from_jet_scale, 1+self.WJets_normalization_uncertainty_from_jet_scale, 1+self.TTbar_normalization_uncertainty_from_jet_scale, 1-self.TTbar_normalization_uncertainty_from_jet_scale, 1+self.STop_normalization_uncertainty_from_jet_scale, 1-self.STop_normalization_uncertainty_from_jet_scale, 1+self.VV_normalization_uncertainty_from_jet_scale, 1-self.VV_normalization_uncertainty_from_jet_scale) )        

          if self.ggH_normalization_uncertainty_from_jet_res!=0 and self.vbf_normalization_uncertainty_from_jet_res!=0 and self.WJets_normalization_uncertainty_from_jet_res!=0 and self.TTbar_normalization_uncertainty_from_jet_res!=0 and self.STop_normalization_uncertainty_from_jet_res!=0 and self.VV_normalization_uncertainty_from_jet_res!=0 :
             
           datacard_out.write( "\nCMS_res_j lnN   %0.3f     %0.3f     %0.3f/%0.3f    %0.3f/%0.3f   %0.3f/%0.3f   %0.3f/%0.3f  "%(1+self.ggH_normalization_uncertainty_from_jet_res, 1+self.vbf_normalization_uncertainty_from_jet_res, 1-self.WJets_normalization_uncertainty_from_jet_res, 1+self.WJets_normalization_uncertainty_from_jet_res, 1+self.TTbar_normalization_uncertainty_from_jet_res, 1-self.TTbar_normalization_uncertainty_from_jet_res, 1+self.STop_normalization_uncertainty_from_jet_res, 1-self.STop_normalization_uncertainty_from_jet_res, 1+self.VV_normalization_uncertainty_from_jet_res, 1-self.VV_normalization_uncertainty_from_jet_res) )        

          if self.ggH_normalization_uncertainty_from_lep_scale!=0 and self.vbf_normalization_uncertainty_from_lep_scale!=0 and self.WJets_normalization_uncertainty_from_lep_scale!=0 and self.TTbar_normalization_uncertainty_from_lep_scale!=0 and self.STop_normalization_uncertainty_from_lep_scale!=0 and self.VV_normalization_uncertainty_from_lep_scale!=0 : 

           datacard_out.write( "\nCMS_scale_%s lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    "%(channel,self.ggH_normalization_uncertainty_from_lep_scale, self.vbf_normalization_uncertainty_from_lep_scale, self.WJets_normalization_uncertainty_from_lep_scale, self.TTbar_normalization_uncertainty_from_lep_scale, self.STop_normalization_uncertainty_from_lep_scale, self.VV_normalization_uncertainty_from_lep_scale) )        


          if self.ggH_normalization_uncertainty_from_lep_res!=0 and self.vbf_normalization_uncertainty_from_lep_res!=0 and self.WJets_normalization_uncertainty_from_lep_res!=0 and self.TTbar_normalization_uncertainty_from_lep_res!=0 and self.STop_normalization_uncertainty_from_lep_res!=0 and self.VV_normalization_uncertainty_from_lep_res!=0 : 

           datacard_out.write( "\nCMS_res_%s lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f    "%(channel,self.ggH_normalization_uncertainty_from_lep_res, self.vbf_normalization_uncertainty_from_lep_res, self.WJets_normalization_uncertainty_from_lep_res, self.TTbar_normalization_uncertainty_from_lep_res, self.STop_normalization_uncertainty_from_lep_res, self.VV_normalization_uncertainty_from_lep_res) )


          if self.ggH_normalization_uncertainty_from_btag!=0 and self.vbf_normalization_uncertainty_from_btag!=0 and self.WJets_normalization_uncertainty_from_btag!=0 and self.TTbar_normalization_uncertainty_from_btag!=0 and self.STop_normalization_uncertainty_from_btag!=0 and self.VV_normalization_uncertainty_from_btag!=0 : 

           datacard_out.write( "\nCMS_eff_b lnN   %0.3f     %0.3f     %0.3f    %0.3f   %0.3f   %0.3f "%(self.ggH_normalization_uncertainty_from_btag, self.vbf_normalization_uncertainty_from_btag, self.WJets_normalization_uncertainty_from_btag, self.TTbar_normalization_uncertainty_from_btag, self.STop_normalization_uncertainty_from_btag, self.VV_normalization_uncertainty_from_btag) )                  

         if mode == "unbin":
            for i in range(len(params_list)):
                if TString(params_list[i].GetName()).Contains("Deco_TTbar_signal_region"):
                    datacard_out.write( "\n%s param  %0.1f  %0.1f "%( params_list[i].GetName(), params_list[i].getVal(), params_list[i].getError() ) ) 
                else:
                    datacard_out.write( "\n%s param  %0.1f  %0.1f "%( params_list[i].GetName(), params_list[i].getVal(), params_list[i].getError() ) ) 
         if mode == "counting":
            datacard_out.write( "\nShape    lnN       -         -             %0.3f    -       -       -      -"%(1+self.rrv_counting_uncertainty_from_shape_uncertainty.getError()))


    ### in order to get the pull
    def read_workspace(self):

        ### Taket the workspace for limits
        file      = TFile(self.file_rlt_root) ;
        workspace = file.Get("workspace4limit_") ;
        workspace.Print()

        parameters_workspace = workspace.allVars();
        par = parameters_workspace.createIterator();
        par.Reset();
        param = par.Next();
        while (param):
            param.Print();
            param = par.Next();
        print "___________________________________________________"

        workspace.data("data_obs_%s"%(self.channel)).Print()
        print "_________________ Pdf in the Workspace  __________________________________"
        pdfs_workspace = workspace.allPdfs();
        par = pdfs_workspace.createIterator();
        par.Reset();
        param = par.Next();
        while (param):
            param.Print();
            param = par.Next();
        print "___________________________________________________"


        rrv_x = workspace.var("rrv_mass_lvj");
        data_obs = workspace.data("data_obs_%s"%(self.channel));
        model_pdf_ggH   = workspace.pdf("ggH_%s"%(self.channel));
        model_pdf_vbfH  = workspace.pdf("qqH_%s"%(self.channel));
        model_pdf_WJets = workspace.pdf("WJets_%s"%(self.channel));
        model_pdf_VV = workspace.pdf("VV_%s"%(self.channel));        
        if options.jetBin == "_2jet": model_pdf_WW_EWK = workspace.pdf("WW_EWK_%s"%(self.channel));
        model_pdf_TTbar = workspace.pdf("TTbar_%s"%(self.channel));
        model_pdf_STop  = workspace.pdf("STop_%s"%(self.channel));

        model_pdf_ggH.Print();
        model_pdf_vbfH.Print();
        model_pdf_WJets.Print();
        model_pdf_VV.Print();
        model_pdf_TTbar.Print();
        model_pdf_STop.Print();
        
        if options.jetBin == "_2jet": model_pdf_WW_EWK.Print();

        rrv_number_ggH = workspace.var("rate_%s_for_unbin"%(self.higgs_sample));
        rrv_number_vbfH = workspace.var("rate_%s_for_unbin"%(self.vbfhiggs_sample));
        rrv_number_WJets = workspace.var("rate_WJets_for_unbin");
        rrv_number_VV = workspace.var("rate_VV_for_unbin");
        if options.jetBin == "_2jet": rrv_number_WW_EWK = workspace.var("rate_WW_EWK_for_unbin");        
        rrv_number_TTbar = workspace.var("rate_TTbar_for_unbin");
        rrv_number_STop = workspace.var("rate_STop_for_unbin");

        rrv_number_ggH.Print();
        rrv_number_vbfH.Print();
        rrv_number_WJets.Print();
        rrv_number_VV.Print();
        if options.jetBin == "_2jet": rrv_number_WW_EWK.Print();       
        rrv_number_TTbar.Print();
        rrv_number_STop.Print();

        if options.jetBin == "_2jet": 
         rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC","rrv_number_Total_background_MC",
                rrv_number_WJets.getVal()+
                rrv_number_VV.getVal()+
                rrv_number_WW_EWK.getVal()+                                                    
                rrv_number_TTbar.getVal()+
                rrv_number_STop.getVal());
         rrv_number_Total_background_MC.setError(TMath.Sqrt(
                rrv_number_WJets.getError()* rrv_number_WJets.getError()+
                rrv_number_VV.getError()* rrv_number_VV.getError()+
                rrv_number_WW_EWK.getError()* rrv_number_WW_EWK.getError()+                
                rrv_number_TTbar.getError()* rrv_number_TTbar.getError()+
                rrv_number_STop.getError() *rrv_number_STop.getError() 
                ));

         model_Total_background_MC = RooAddPdf("model_Total_background_MC","model_Total_background_MC",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_WW_EWK,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_WW_EWK,rrv_number_TTbar,rrv_number_STop));

        else:  

         rrv_number_Total_background_MC = RooRealVar("rrv_number_Total_background_MC","rrv_number_Total_background_MC",
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

         model_Total_background_MC = RooAddPdf("model_Total_background_MC","model_Total_background_MC",RooArgList(model_pdf_WJets,model_pdf_VV,model_pdf_TTbar,model_pdf_STop),RooArgList(rrv_number_WJets,rrv_number_VV,rrv_number_TTbar,rrv_number_STop));


        scale_number_ggH  = rrv_number_ggH.getVal()/data_obs.sumEntries();
        scale_number_vbfH = rrv_number_vbfH.getVal()/data_obs.sumEntries();
        scale_number_Total_background_MC = rrv_number_Total_background_MC.getVal()/data_obs.sumEntries();

        ######################################################### 
        ###### Final Plot of mWW invarant mass distribution #####
        #########################################################
        
        mplot = rrv_x.frame(RooFit.Title("check_workspace"),RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
        data_obs.plotOn(mplot ,RooFit.DataError(RooAbsData.SumW2), RooFit.Name("data_invisible"),RooFit.Invisible());

        #### create the frame
        if options.jetBin == "_2jet" :
         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_%s,WW_EWK_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WW_EWK"), RooFit.Components("WW_EWK_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(1), RooFit.VLines());        

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("TTbar_%s,STop_%s"%(self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop_%s"%(self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(1), RooFit.VLines());                

         #solid line
         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_%s,WW_EWK_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WW_EWK_line_invisible"), RooFit.Components("WW_EWK_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_%s,STop_%s"%(self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop_%s"%(self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        #### create the frame
        else :
         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.Components("WJets_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV"), RooFit.Components("VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar"), RooFit.Components("TTbar_%s,STop_%s"%(self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(1), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop"), RooFit.Components("STop_%s"%(self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(1), RooFit.VLines());                


         #solid line
         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets_line_invisible"), RooFit.Components("WJets_%s,VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("VV_line_invisible"), RooFit.Components("VV_%s,TTbar_%s,STop_%s"%(self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("TTbar_line_invisible"), RooFit.Components("TTbar_%s,STop_%s"%(self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

         model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("STop_line_invisible"), RooFit.Components("STop_%s"%(self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());


        if self.higgs_sample == "ggH600" or self.higgs_sample == "ggH700":
           signal_scale = 2;
        else: 
           signal_scale = 2;
        
        if self.higgs_sample == "ggH600":
            model_pdf_ggH.plotOn(mplot,RooFit.Normalization(scale_number_ggH*signal_scale),RooFit.Name("ggH#times%s, m_{H}=0.6TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["ggH"]), RooFit.LineStyle(2), RooFit.VLines());
            model_pdf_vbfH.plotOn(mplot,RooFit.Normalization(scale_number_vbfH*signal_scale),RooFit.Name("qqH#times%s, m_{H}=0.6TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["vbfH"]), RooFit.LineStyle( 9), RooFit.VLines());

        if self.higgs_sample=="ggH700":
            model_pdf_ggH.plotOn(mplot,RooFit.Normalization(scale_number_ggH*signal_scale),RooFit.Name("ggH#times%s, m_{H}=0.7TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["ggH"]), RooFit.LineStyle(2), RooFit.VLines());
            model_pdf_vbfH.plotOn(mplot,RooFit.Normalization(scale_number_vbfH*signal_scale),RooFit.Name("qqH#times%s, m_{H}=0.7TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["vbfH"]), RooFit.LineStyle( 9), RooFit.VLines());

        if self.higgs_sample=="ggH800":
            model_pdf_ggH.plotOn(mplot,RooFit.Normalization(scale_number_ggH*signal_scale),RooFit.Name("ggH#times%s, m_{H}=0.8TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["ggH"]), RooFit.LineStyle(2), RooFit.VLines());
            model_pdf_vbfH.plotOn(mplot,RooFit.Normalization(scale_number_vbfH*signal_scale),RooFit.Name("qqH#times%s, m_{H}=0.8TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["vbfH"]), RooFit.LineStyle( 9), RooFit.VLines());

        if self.higgs_sample=="ggH900":
            model_pdf_ggH.plotOn(mplot,RooFit.Normalization(scale_number_ggH*signal_scale),RooFit.Name("ggH#times%s, m_{H}=0.9TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["ggH"]), RooFit.LineStyle(2), RooFit.VLines());
            model_pdf_vbfH.plotOn(mplot,RooFit.Normalization(scale_number_vbfH*signal_scale),RooFit.Name("qqH#times%s, m_{H}=0.9TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["vbfH"]), RooFit.LineStyle( 9), RooFit.VLines());

        if self.higgs_sample=="ggH1000":
            model_pdf_ggH.plotOn(mplot,RooFit.Normalization(scale_number_ggH*signal_scale),RooFit.Name("ggH#times%s, m_{H}=1TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["ggH"]), RooFit.LineStyle(2), RooFit.VLines());
            model_pdf_vbfH.plotOn(mplot,RooFit.Normalization(scale_number_vbfH*signal_scale),RooFit.Name("qqH#times%s, m_{H}=1TeV"%(signal_scale)),RooFit.DrawOption("L"), RooFit.LineColor(self.color_palet["vbfH"]), RooFit.LineStyle( 9), RooFit.VLines());

        data_obs.plotOn(mplot, RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        model_Total_background_MC.plotOn(mplot,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Invisible(),RooFit.Name("model_mc"));

        rfresult = RooFitResult() ;
        mplot_pull = get_pull(rrv_x,mplot,data_obs,model_Total_background_MC,rfresult,"data","model_mc",0,1);

        self.FloatingParams.Print("v");
        if options.closuretest == 0:
            draw_error_band_ws(data_obs,model_Total_background_MC,rrv_x.GetName(),rrv_number_Total_background_MC,self.FloatingParams,workspace,mplot,self.color_palet["Uncertainty"],"F");
        else:
            draw_error_band_ws(data_obs,model_Total_background_MC,rrv_x.GetName(),rrv_number_Total_background_MC,self.FloatingParams,workspace,mplot,self.color_palet["Uncertainty"],"F");

        self.leg = legend4Plot(mplot,0,0.,0.06,0.16,0.,1,self.channel);
            
        mplot.addObject(self.leg);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.3);

        datahist = data_obs.binnedClone(data_obs.GetName()+"_binnedClone",data_obs.GetName()+"_binnedClone")
        Nbin = int(rrv_x.getBins()); 
        nparameters = self.FloatingParams.getSize();        
        ChiSquare = model_Total_background_MC.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
        chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters);
        ## Add Chisquare to mplot_pull
        cs = TLatex(0.75,0.8,"#chi^{2}/ndf = %0.2f "%(float(chi_over_ndf)));
        cs.SetNDC();
        cs.SetTextSize(0.12);
        cs.AppendPad("same");
        mplot_pull.addObject(cs)

        print "nPar=%s, chiSquare=%s/%s"%(nparameters, ChiSquare.getVal()*(Nbin - nparameters), (Nbin - nparameters) );

        parameters_list = RooArgList();
        draw_canvas_with_pull(mplot,mplot_pull,parameters_list,"plots_%s_%s_g1/mlvj_fitting/"%(self.channel,self.wtagger_label),"check_worksapce","",channel,0,1,GetLumi());


        ################################        
        ## Final Plot for the w+jet ####
        ################################
        mplot_wjet = rrv_x.frame(RooFit.Title("check_workspace_wjet"));
        data_obs.plotOn(mplot ,RooFit.DataError(RooAbsData.SumW2), RooFit.Name("data_invisible"),RooFit.Invisible());
        scale_number_Total_background_MC = rrv_number_WJets.getVal()/data_obs.sumEntries();
        model_pdf_WJets.plotOn(mplot_wjet,RooFit.Normalization(scale_number_Total_background_MC),RooFit.Name("WJets"), RooFit.LineColor(self.color_palet["WJets"]), RooFit.VLines());
        draw_error_band_ws(data_obs,model_pdf_WJets,rrv_x.GetName(),rrv_number_WJets,self.FloatingParams_wjet,workspace,mplot_wjet,self.color_palet["Uncertainty"],"F");
        draw_canvas(mplot_wjet,"plots_%s_%s_g1/mlvj_fitting/"%(self.channel,self.wtagger_label),"wjet_shape",self.channel,GetLumi(),0,1,0);
        

    ######## +++++++++++
    def get_data(self, signal_model="CB_v1"):
        print "############### get_data ########################"
        self.get_mj_and_mlvj_dataset(self.file_data,"_data")     
        getattr(self.workspace4limit_,"import")(self.workspace4fit_.var("rrv_number_dataset_signal_region_data_%s_mlvj"%(self.channel)).clone("observation_for_counting"))
 

    #### Define the steps to fit signal distribution in the mj and mlvj spectra
    def fit_Signal(self):
        print "############# fit_Signal #################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_ggH,"_%s"%(self.higgs_sample));

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%smassvbf_jes_up"%(self.higgs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));

         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         self.workspace4fit_.Print();
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%smassvbf_jes_dn"%(self.higgs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%smassvbf_jer"%(self.higgs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%smassvbf_jer_up"%(self.higgs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%smassvbf_jer_dn"%(self.higgs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_ggH,"_%s"%(self.higgs_sample),"_signal_region",self.mlvj_shape["ggH"],self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.higgs_sample));
        self.get_mj_and_mlvj_dataset(self.file_vbfH,"_%s"%(self.vbfhiggs_sample));

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%smassvbf_jes_up"%(self.vbfhiggs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%smassvbf_jes_dn"%(self.vbfhiggs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%smassvbf_jer"%(self.vbfhiggs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%smassvbf_jer_up"%(self.vbfhiggs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%smassvbf_jer_dn"%(self.vbfhiggs_sample),"_signal_region","CB_v1",self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"_signal_region",self.mlvj_shape["vbfH"],self.channel,self.wtagger_label,0,0,0,0,"_%s"%(self.vbfhiggs_sample));
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
    
        print "________________________________________________________________________"
        
 
    ##### Define the steps to fit WJets MC in the mj and mlvj spectra
    def fit_WJets(self):
        print "######################### fit_WJets ########################"        
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj
        if not options.jetBin == "_2jet": self.get_mj_and_mlvj_dataset(self.file_WJets1_mc,"_WJets1")# to get the shape of m_lvj
        self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets01")# to get the shape of m_lvj

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_up",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_dn",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_up",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_dn",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets01",self.mj_shape["WJets01"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        if not options.jetBin == "_2jet": fit_mj_single_MC(self.workspace4fit_,self.file_WJets1_mc,"_WJets1",self.mj_shape["WJets1"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_up","_sb_lo",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_dn","_sb_lo",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer","_sb_lo",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_up","_sb_lo",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_dn","_sb_lo",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets01","_sb_lo",self.mlvj_shape["WJets01"],self.channel,self.wtagger_label,0,0,0,0,"_WJets01");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        if not options.jetBin == "_2jet": fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets1_mc,"_WJets1","_sb_lo",self.mlvj_shape["WJets1"],self.channel,self.wtagger_label,1,0,0,0,"_WJets01");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_up","_signal_region",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jes_dn","_signal_region",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer","_signal_region",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_up","_signal_region",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0massvbf_jer_dn","_signal_region",self.MODEL_4_mlvj,self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets01","_signal_region",self.mlvj_shape["WJets01"],self.channel,self.wtagger_label,0,0,0,0,"_WJets01");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        if not options.jetBin == "_2jet": fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets1_mc,"_WJets1","_signal_region",self.mlvj_shape["WJets1"],self.channel,self.wtagger_label,0,0,0,0,"_WJets01");
        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0","_signal_region",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,0,0,0,"_WJets0");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        
        print "________________________________________________________________________"

    ##### Define the steps to fit VV MC in the mj and mlvj spectra
    def fit_VV(self):
        print "############################# fit_VV ################################"

        ### Build the dataset        
        self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV","jet_mass_pr")# to get the shape of m_lvj

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_up",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_dn",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_up",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_dn",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_up","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_dn","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_up","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_dn","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV","_sb_lo",self.mlvj_shape["VV"],self.channel,self.wtagger_label,1,0,1,0,"_VV");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_up","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jes_dn","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_up","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VVmassvbf_jer_dn","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV","_signal_region",self.mlvj_shape["VV"],self.channel,self.wtagger_label,0,0,0,0,"_VV");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
                
        print "________________________________________________________________________"


    ##### Define the steps to fit WW_EWK MC in the mj and mlvj spectra
    def fit_WW_EWK(self):
        print "############################# fit_WW_EWK ################################"

        ### Build the dataset        
        self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_up",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_dn",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_up",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_dn",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWK",self.mj_shape["WW_EWK"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_up","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_dn","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_up","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_dn","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWK","_sb_lo",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_up","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jes_dn","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_up","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWKmassvbf_jer_dn","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,0,0,0,0,"_WW_EWK");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_WW_EWK_mc,"_WW_EWK","_signal_region",self.mlvj_shape["WW_EWK"],self.channel,self.wtagger_label,1,0,0,0,"_WW_EWK");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        
        print "________________________________________________________________________"


    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_TTbar(self):

        print "################################ fit_TTbar #########################################"
        ### Build the dataset

        self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj        
        self.get_mj_and_mlvj_dataset(self.file_TTbar_mcanlo_mc,"_TTbar_mcanlo")# to get the shape of m_lvj

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_up",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata); 
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_dn",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_up",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_dn",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlo",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_up",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_dn",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_up",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_dn",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar",self.mj_shape["TTbar"],self.channel,self.wtagger_label,1);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_up","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_dn","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_up","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_dn","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        
        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlo","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_up","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_dn","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_up","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_dn","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar","_sb_lo",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_up","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jes_dn","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_up","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlomassvbf_jer_dn","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mcanlo_mc,"_TTbar_mcanlo","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar_mcanlo");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_up","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jes_dn","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_up","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbarmassvbf_jer_dn","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar","_signal_region",self.mlvj_shape["TTbar"],self.channel,self.wtagger_label,1,0,0,0,"_TTbar");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
    
        print "________________________________________________________________________"


    ##### Define the steps to fit TTbar MC in the mj and mlvj spectra
    def fit_STop(self):
        print "############################### fit_STop #########################################"
        ### Build the dataset
        self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj

        if not options.skipJetSystematics:
         fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_up",self.mj_shape["STop"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_dn",self.mj_shape["STop"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer",self.mj_shape["STop"],self.channel,self.wtagger_label,1,options.pseudodata);        
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_up",self.mj_shape["STop"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_dn",self.mj_shape["STop"],self.channel,self.wtagger_label,1,options.pseudodata);
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop",self.mj_shape["STop"],self.channel,self.wtagger_label,1);
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_up","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_dn","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_up","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_dn","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop","_sb_lo",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        if not options.skipJetSystematics:
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_up","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jes_dn","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_up","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());
         fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STopmassvbf_jer_dn","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
         self.workspace4fit_.writeToFile(self.tmpFile.GetName());

        fit_mlvj_model_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop","_signal_region",self.mlvj_shape["STop"],self.channel,self.wtagger_label,1,0,0,0,"_STop");
        self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        
        print "________________________________________________________________________"  


    ##### Fit of all the MC in both mj and mlvj : Signal, TTbar, STop, VV and Wjets
    def fit_AllSamples_Mj_and_Mlvj(self):
        print "################### fit_AllSamples_Mj_and_Mlvj #####################"
        self.fit_Signal();
        self.fit_STop();
        self.fit_VV();
        if options.jetBin == "_2jet": self.fit_WW_EWK();        
        self.fit_WJets();
        self.fit_TTbar();

        print "________________________________________________________________________"
                                                                    

    ##### Analysis with sideband alpha correction
    def analysis_sideband_correction_method1(self):
      print "##################### Start sideband correction full analysis ##############";
      ### Fit all MC components in both mj and mlvj      
      self.fit_AllSamples_Mj_and_Mlvj();
      ### take the real data
      self.get_data();
      ### fit the WJets Normalization into the signal region -> no jet mass fluctuation has been done
      if not options.skipJetSystematics:
       self.fit_WJetsNorm(1);
      else:
       self.fit_WJetsNorm(0);
          
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       
      ### fit data in the mlvj low sideband with two different models      
      if not options.skipJetSystematics:
       fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0massvbf_jes_up","massvbf_jes_up","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0massvbf_jes_dn","massvbf_jes_dn","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0massvbf_jer_up","massvbf_jer_up","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0massvbf_jer_dn","massvbf_jer_dn","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0massvbf_jer","massvbf_jer","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
        

      fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets01","","_sb_lo",self.mlvj_shape["WJets01"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());
      if not options.jetBin == "_2jet": fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets1","","_sb_lo",self.mlvj_shape["WJets1"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());
      fit_mlvj_in_Mj_sideband(self.workspace4fit_,self.color_palet,self.mlvj_shape,"_WJets0","","_sb_lo",self.mlvj_shape["WJets0"],self.channel,self.wtagger_label,0,1,options.pseudodata,options.jetBin);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());

     
      ### Prepare the workspace and datacards
      #### Call the alpha evaluation in automatic
      if not options.skipJetSystematics:
      
       get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0massvbf_jes_up",self.mlvj_shape["WJets0"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0massvbf_jes_dn",self.mlvj_shape["WJets0"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0massvbf_jer",self.mlvj_shape["WJets0"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0massvbf_jer_up",self.mlvj_shape["WJets0"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0massvbf_jer_dn",self.mlvj_shape["WJets0"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
       self.workspace4fit_.writeToFile(self.tmpFile.GetName());
       
      get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets01",self.mlvj_shape["WJets01"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());
      if not options.jetBin == "_2jet": get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets1",self.mlvj_shape["WJets1"],"","_mlvj","4fit",self.channel,self.wtagger_label,1);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());
      get_WJets_mlvj_correction_sb_lo_to_signal_region(self.workspace4fit_,"_WJets0",self.mlvj_shape["WJets0"],self.mlvj_shape["WJets01"],"_mlvj","4fit",self.channel,self.wtagger_label,1);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());

      ### Fix the pdf of signal, TTbar, STop and VV in the signal region
      
      fix_Model(self.workspace4fit_,"_%s"%(self.higgs_sample),"_signal_region","_mlvj",self.mlvj_shape["ggH"],self.channel,"",0);
      fix_Model(self.workspace4fit_,"_%s"%(self.vbfhiggs_sample),"_signal_region","_mlvj",self.mlvj_shape["vbfH"],self.channel,"",0);
      fix_Model(self.workspace4fit_,"_TTbar","_signal_region","_mlvj",self.mlvj_shape["TTbar"],self.channel,"",0);
      fix_Model(self.workspace4fit_,"_STop","_signal_region","_mlvj",self.mlvj_shape["STop"],self.channel,"",0);
      fix_Model(self.workspace4fit_,"_VV","_signal_region","_mlvj",self.mlvj_shape["VV"],self.channel,"",0);
      if options.jetBin == "_2jet" : fix_Model(self.workspace4fit_,"_WW_EWK","_signal_region","_mlvj",self.mlvj_shape["WW_EWK"],self.channel,"",0);
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());

      ### Call the evaluation of the normalization in the signal region for signal, TTbar, VV, STop, and WJets after the extrapolation via alpha
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_%s"%(self.higgs_sample),self.mlvj_shape["ggH"],"_signal_region",self.channel);
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_%s"%(self.vbfhiggs_sample),self.mlvj_shape["vbfH"],"_signal_region",self.channel);
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_TTbar",self.mlvj_shape["TTbar"],"_signal_region",self.channel);
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_STop",self.mlvj_shape["STop"],"_signal_region",self.channel);
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_VV",self.mlvj_shape["VV"],"_signal_region",self.channel);
      if options.jetBin == "_2jet":
       get_mlvj_normalization_insignalregion(self.workspace4fit_,"_WW_EWK",self.mlvj_shape["WW_EWK"],"_signal_region",self.channel);
    
      get_mlvj_normalization_insignalregion(self.workspace4fit_,"_WJets0",self.mlvj_shape["WJets0"],"_signal_region",self.channel,1);  
      self.workspace4fit_.writeToFile(self.tmpFile.GetName());

      self.prepare_limit("sideband_correction_method1",1);
      self.read_workspace();
       
    ####### +++++++++++++++
    def analysis_sideband_correction_method1_without_shape_and_psmodel_systermatic(self):
        self.fit_AllSamples_Mj_and_Mlvj()
        self.get_data();
        self.fit_WJetsNormalization_in_Mj_signal_region("_WJets0");
        self.fit_mlvj_in_Mj_sideband("_WJets0","_sb_lo", self.MODEL_4_mlvj,1);
        self.prepare_limit("sideband_correction_method1");
        self.read_workspace(1);


### funtion to run the analysis without systematics
def pre_limit_sb_correction_without_systermatic( channel, higgs_sample="HWWMH600", in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700, in_mj_min=30, in_mj_max=140, in_mlvj_min=400, in_mlvj_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1"): # the WJets M_lvj shape and normalization are from sb_correction

    print "#################### pre_limit_sb_correction_without_systermatic: channel %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min mlvj %f-f, fit model %s and alternate %s ######################"%(channel,higgs_sample,in_mlvj_signal_region_min,in_mlvj_signal_region_max,in_mj_min,in_mj_max,in_mlvj_mi,in_mlvj_max,fit_model,fit_model_alter);

    boostedW_fitter=doFit_wj_and_wlvj(channel, higgs_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max,fit_model, fit_model_alter);
    boostedW_fitter.analysis_sideband_correction_method1_without_shape_and_psmodel_systermatic();


### funtion to run the complete alpha analysis
def pre_limit_sb_correction(method, channel, higgs_sample="HWWMH600", in_mlvj_signal_region_min=500, in_mlvj_signal_region_max=700, in_mj_min=30, in_mj_max=140, in_mlvj_min=400, in_mlvj_max=1400, fit_model="ErfExp_v1", fit_model_alter="ErfPow_v1"): # the WJets M_lvj shape and normalization are from sb_correction

    print "#################### pre_limit_sb_correction: channel %s, signal %s, max and min signal region %f-%f, max and min mJ %f-%f, max and min mlvj %f-%f, fit model %s and alternate %s ######################"%(channel,higgs_sample,in_mlvj_signal_region_min,in_mlvj_signal_region_max,in_mj_min,in_mj_max,in_mlvj_min,in_mlvj_max,fit_model,fit_model_alter);

    boostedW_fitter=doFit_wj_and_wlvj(channel, higgs_sample, in_mlvj_signal_region_min, in_mlvj_signal_region_max, in_mj_min, in_mj_max, in_mlvj_min, in_mlvj_max,fit_model, fit_model_alter);
    getattr(boostedW_fitter,"analysis_sideband_correction_%s"%(method) )();


### funtion to run the analysis without systematic
def pre_limit_simple(channel):
    print "######################### pre_limit_simple for %s sampel"%(channel)
    pre_limit_sb_correction_without_systermatic(channel, "ggH600",550, 700,40,130,400,1000,"ErfPowExp_v1","ErfPow2_v1");
    
### function to check the workspace once it has already created
def check_workspace(channel, higgs):
    boostedW_fitter = doFit_wj_and_wlvj(channel,higgs);
    boostedW_fitter.read_workspace()

####################################
######### Main Programme ###########
####################################    
                            
if __name__ == '__main__':

    channel=options.channel;#mu or el; default is mu;

    if options.check:
        print '################# check workspace for %s sample'%(channel);
        check_workspace(channel,"ggH600");

    if options.simple and ( not options.multi) and ( not options.check) :
        print '################# simple mode for %s sample'%(channel)
        pre_limit_simple(channel);

    if options.multi:
        print '################# multi mode for %s sample'%(channel)
        pre_limit_sb_correction("method1",sys.argv[1],sys.argv[2],int(sys.argv[3]),int(sys.argv[4]),int(sys.argv[5]),int(sys.argv[6]),int(sys.argv[7]),int(sys.argv[8]), sys.argv[9], sys.argv[10] );
