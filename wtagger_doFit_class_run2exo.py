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


from ROOT import gROOT, TPaveLabel, gStyle, gSystem, TGaxis, TStyle, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooArgusBG, RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooChebychev, RooSimultaneous, RooGenericPdf,RooConstVar, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan, kMagenta, kWhite


############################################
# Job steering #
############################################

def foo_callback(option, opt, value, parser):
          setattr(parser.values, option.dest, value.split(','))


parser = OptionParser()

parser.add_option('-a', '--additioninformation',action="store",type="string",dest="additioninformation",default="EXO")
parser.add_option('--cprime', action="store",type="int",dest="cprime",default=10)
parser.add_option('--BRnew', action="store",type="int",dest="BRnew",default=0)

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")


parser.add_option('--fitwtagger', action='store_true', dest='fitwtagger', default=True, help='fit wtagger jet in ttbar control sample')
parser.add_option('--fitwtaggersim', action='store_true', dest='fitwtaggersim', default=False, help='fit wtagger jet in ttbar control sample with mu and el samples simultaneously')

parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")

parser.add_option('--herwig', action="store",type="int",dest="herwig",default=0)
parser.add_option('--pTbin', action="callback",callback=foo_callback,type="string",dest="pTbin",default="")
parser.add_option('--shift', action="store", type="int",dest="shift",default=0)
parser.add_option('--smear', action="store", type="int",dest="smear",default=0)

parser.add_option('--tau2tau1cutHP', action="store", type="float",dest="tau2tau1cutHP",default=0.6)
parser.add_option('--tau2tau1cutLP', action="store", type="float",dest="tau2tau1cutLP",default=0.75)

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

from ROOT import fit_mj_single_MC, fit_mlvj_model_single_MC, fit_WJetsNormalization_in_Mj_signal_region, fit_mlvj_in_Mj_sideband, get_WJets_mlvj_correction_sb_lo_to_signal_region, get_mlvj_normalization_insignalregion, fit_genHMass, SystematicUncertaintyHiggs_2jetBin, SystematicUncertaintyHiggs_01jetBin, ScaleFactorTTbarControlSampleFit,DrawScaleFactorTTbarControlSample

from ROOT import *

gInterpreter.GenerateDictionary("std::map<std::string,std::string>", "map;string;string")
gInterpreter.GenerateDictionary("std::vector<std::string>", "vector;string")


###############################
## doFit Class Implemetation ##
###############################

class doFit_wj_and_wlvj:

    ## contructor taking channel, signal name, range in mj, label and a workspace
    def __init__(self, in_channel,in_signal_sample, in_mj_min=40, in_mj_max=150, label="", input_workspace=None):

        print " ";
        print "#####################################################";
        print "## Constructor of the fit object doFit_wj_and_wlvj ##";
        print "#####################################################";
        print " ";

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9);
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9);

        ### set the channel type --> electron or muon
        self.channel = in_channel;

        ### shapes to be used in mj                                                                                                                                          
        self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string)();
        self.mj_shape["TTbar"]   = "2Gaus_ErfExp";

        if self.channel == "mu":
            self.mj_shape["STop"]             = "ErfExpGaus_sp";
            self.mj_shape["STop_fail"]        = "Exp";
            self.mj_shape["STop_extremefail"] = "Exp";
            self.mj_shape["VV"]               = "ErfExpGaus_sp";
            self.mj_shape["VV_fail"]          = "ExpGaus";
            self.mj_shape["VV_extremefail"]   = "Exp";
        else:
            self.mj_shape["STop"]              = "ExpGaus";
            self.mj_shape["STop_fail"]         = "ExpGaus";
            self.mj_shape["STop_extremefail"]  = "Exp";
            self.mj_shape["VV"]                = "Gaus";
            self.mj_shape["VV_fail"]           = "Exp";
            self.mj_shape["VV_extremefail"]    = "Exp";
            
        self.mj_shape["WJets0"]              = "ErfExp";
        self.mj_shape["WJets0_fail"]         = "Exp";
        self.mj_shape["WJets0_extremefail"]  = "Exp";

        self.mj_shape["bkg_data"]         = "ErfExp_ttbar";
        self.mj_shape["bkg_data_fail"]    = "ErfExp_ttbar_failtau2tau1cut";
        self.mj_shape["signal_data"]      = "2Gaus_ttbar" ;
        self.mj_shape["signal_data_fail"] = "GausChebychev_ttbar_failtau2tau1cut";
        self.mj_shape["bkg_mc"]           = "ErfExp_ttbar";
        self.mj_shape["bkg_mc_fail"]      = "ErfExp_ttbar_failtau2tau1cut";
        self.mj_shape["signal_mc"]        = "2Gaus_ttbar" ;
        self.mj_shape["signal_mc_fail"]   = "GausChebychev_ttbar_failtau2tau1cut";
        self.mj_shape["data_extremefail"]     = "Exp_ttbar_extremefailtau2tau1cut";
        self.mj_shape["data_bkg_extremefail"] = "Exp_bkg_extremefailtau2tau1cut";
        self.mj_shape["mc_extremefail"]       = "Exp_ttbar_extremefailtau2tau1cut";
        self.mj_shape["mc_bkg_extremefail"]   = "Exp_bkg_extremefailtau2tau1cut";

        self.Lumi=1263;

        ### Set the mj binning for plots
        self.BinWidth_mj = 5.;
        
        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor = 1.;

        ## set the range max properly in order to have a integer number of bins
        self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor;
        nbins_mj         = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max        = in_mj_min+nbins_mj*self.BinWidth_mj;

        ## declare the RooRealVar + binning
        rrv_mass_j = RooRealVar("rrv_mass_j","pruned jet mass",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV");
        rrv_mass_j.setBins(nbins_mj);

        ## create the workspace and import the variable
        if input_workspace is None:
            self.workspace4fit_ = RooWorkspace("workspace4fit"+label+"_","Workspace4fit"+label+"_");
        else:
            self.workspace4fit_ = input_workspace;
        getattr(self.workspace4fit_,"import")(rrv_mass_j);

        ## Region definition --> signal region between 65 and 105 GeV
        self.mj_sideband_lo_min = in_mj_min;
        self.mj_sideband_lo_max = 65;
        self.mj_signal_min      = 65;
        self.mj_signal_max      = 105;
        self.mj_sideband_hi_min = 135;
        self.mj_sideband_hi_max = in_mj_max;

        ## define ranges on mj
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("controlsample_fitting_range",40,150);

        ## directory where are the trees to be run
        self.file_Directory = "/afs/cern.ch/user/l/lbrianza/work/public/WWTree_3nov/WWTree_%s/"%(self.channel);

        ## taking the root file for data and mc
        self.file_data = ("WWTree_data.root");

        self.signal_sample = in_signal_sample;
        self.file_signal   = ("WWTree_%s.root"%(self.signal_sample));

        self.file_pseudodata        = ("WWTree_pseudodata.root");
        self.file_pseudodata_herwig = ("WWTree_pseudodata.root");

        if self.channel != "em":
            self.file_WJets0_mc = ("WWTree_WJets.root");
        else:
            self.file_WJets0_mc = ("WWTree_WJets.root");

        self.file_WJets1_mc = ("WWTree_WJets.root");

        self.file_VV_mc            = ("WWTree_VV.root");# WW+WZ
        self.file_TTbar_mc         = ("WWTree_TTbar.root"); ## powheg TTbar
        self.file_TTbar_herwig     = ("WWTree_TTbar_powheg.root"); ## mc@nlo TTbar
        self.file_TTbar_matchDn_mc = ("WWTree_TTbar.root"); ## madgraph matching down
        self.file_TTbar_matchUp_mc = ("WWTree_TTbar.root"); ## madgraph matching up
        self.file_TTbar_scaleDn_mc = ("WWTree_TTbar.root"); ## madgraph qcd scale down
        self.file_TTbar_scaleUp_mc = ("WWTree_TTbar.root"); ## madgraph qcd scale up
        self.file_TTbar_MG_mc      = ("WWTree_TTbar_madgraph.root"); ## madgraph ttbar
        self.file_STop_mc          = ("WWTree_STop.root"); ##single Top
 
        ## Define the workin poit on the N-subjettines cut

        self.wtagger_label = options.category;

        if self.wtagger_label == "HP" :
            if self.channel == "el"    : self.wtagger_cut = options.tau2tau1cutHP ; self.wtagger_cut_min = 0. ;
            if self.channel == "mu"    : self.wtagger_cut = options.tau2tau1cutHP ; self.wtagger_cut_min = 0. ;
            if self.channel == "em":     self.wtagger_cut = options.tau2tau1cutHP ; self.wtagger_cut_min = 0. ;

        if self.wtagger_label == "LP":
            self.wtagger_cut = options.tau2tau1cutLP ; self.wtagger_cut_min = options.tau2tau1cutHP ;

        if self.wtagger_label == "nocut":
            self.wtagger_cut = 10000;

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

                                                       
        ## Some basic cut vaule
        self.vpt_cut      = 200; ## hadronic and leptonic W cut
        self.mass_lvj_max = 5000.; ## invariant mass of 3 body max
        self.mass_lvj_min = 0.; ## invariant mass of 3 body min
        self.pfMET_cut    = 40; ## missing transverse energy
        self.lpt_cut      = 53; ## lepton pT
        self.deltaPhi_METj_cut = 2.0;

        ## binning in the W jet pT for differential SF study in bins of pT

        if options.pTbin != "" :
         self.ca8_ungroomed_pt_min = int(options.pTbin[0]);
         self.ca8_ungroomed_pt_max = int(options.pTbin[1]);
        else:
         self.ca8_ungroomed_pt_min = 200;
         self.ca8_ungroomed_pt_max = 5000;

        
        ## tighter cut for the electron channel
        if self.channel == "el" or self.channel == "em":
#        if self.channel == "em":
            self.pfMET_cut = 80; self.lpt_cut = 90;

        ## out txt file with info about fit and event couting
        self.file_ttbar_control_txt = "ttbar_control_%s_%s_wtaggercut%s.txt"%(self.signal_sample,self.channel,self.wtagger_label);
        self.file_out_ttbar_control = open(self.file_ttbar_control_txt,"w");

        ### set the TDR Style                                                                                                                                                             
        setTDRStyle();


    ############# -----------
    def fit_mj_TTbar_controlsample(self,in_file_name,label=""):

        ##### Print the final result for the number of events passing the cut and before the cut + the efficiency for the W-tagging -> dataset yields in the signal region
        self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_VV"  +label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_WJets0"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_STop"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar"+label+"_"+self.channel+"_mj").Print()

        number_dataset_signal_region_data_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region: %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj);

        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel));
        self.file_out_ttbar_control.write("event number of data in signal_region: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj));

        number_dataset_signal_region_TotalMC_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        print "event number of TotalMC %s in signal_region: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj);

        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj));


        number_dataset_signal_region_before_cut_data_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region before_cut: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj);
        self.file_out_ttbar_control.write("event number of data in signal_region before_cut: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj));

        number_dataset_signal_region_before_cut_TotalMC_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        print "event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj);

        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj));
                                                             
        # wtagger_eff reweight: only reweight the efficiency difference between MC and data
        wtagger_eff_MC   = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj;
        wtagger_eff_data = number_dataset_signal_region_data_mj/number_dataset_signal_region_before_cut_data_mj;

        wtagger_eff_reweight     = wtagger_eff_data/wtagger_eff_MC;
        wtagger_eff_reweight_err = wtagger_eff_reweight*TMath.Sqrt(number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +number_dataset_signal_region_before_cut_error2_data_mj/number_dataset_signal_region_before_cut_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_before_cut_error2_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj);
        
        print "wtagger efficiency of %s channel"%(self.channel);
        print "wtagger_eff_MC %s = %s "%(label,wtagger_eff_MC);
        print "wtagger_eff_data = %s "%(wtagger_eff_data);
        print "wtagger_eff_reweight %s = %s +/- %s"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err);

        self.file_out_ttbar_control.write("wtagger_eff_MC %s = %s \n"%(label,wtagger_eff_MC ));
        self.file_out_ttbar_control.write("wtagger_eff_data = %s \n"%(wtagger_eff_data ));
        self.file_out_ttbar_control.write("wtagger_eff_reweight %s= %s +/- %s\n"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err));


  
    ##########################################
    ## To build the dataset to be fitted  ####
    ##########################################
        
    def get_mj_and_mlvj_dataset_TTbar_controlsample(self,in_file_name, label, jet_mass="ttb_jet_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn      = TFile(fileIn_name.Data());
        treeIn      = fileIn.Get("otree");
        
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ### dataset of m_j before tau2tau1 cut : Passed
        rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight));
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
    
        ### dataset of m_j before tau2tau1 cut : Total
        rdataset_beforetau2tau1cut_mj     = RooDataSet("rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_beforetau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

   
        ### dataset of m_j failed tau2tau1 cut :
        rdataset_failtau2tau1cut_mj     = RooDataSet("rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_failtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

  
        ### dataset of m_j extreme failed tau2tau1 cut: >0.75
        rdataset_extremefailtau2tau1cut_mj     = RooDataSet("rdataset"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_extremefailtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

       
        ### create the category for the fit
        if TString(label).Contains("herwig"):
         if self.workspace4fit_.cat("category_p_f"+"_herwig_"+self.channel):
             category_p_f = self.workspace4fit_.cat("category_p_f"+"_herwig_"+self.channel);
         else:
             category_p_f = RooCategory("category_p_f"+"_herwig_"+self.channel,"category_p_f"+"_herwig_"+self.channel);
             category_p_f.defineType("pass");
             category_p_f.defineType("fail");
             getattr(self.workspace4fit_,"import")(category_p_f);

        else:
            
         if self.workspace4fit_.cat("category_p_f"+"_"+self.channel):
             category_p_f = self.workspace4fit_.cat("category_p_f"+"_"+self.channel);
         else:
             category_p_f = RooCategory("category_p_f"+"_"+self.channel,"category_p_f"+"_"+self.channel);
             category_p_f.defineType("pass");
             category_p_f.defineType("fail");
             getattr(self.workspace4fit_,"import")(category_p_f);
        
        combData_p_f = RooDataSet("combData_p_f"+label+"_"+self.channel,"combData_p_f"+label+"_"+self.channel,RooArgSet(rrv_mass_j, category_p_f, rrv_weight),RooFit.WeightVar(rrv_weight));
            
        print "N entries: ", treeIn.GetEntries()

        hnum_4region = TH1D("hnum_4region"+label+"_"+self.channel,"hnum_4region"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_4region_error2 = TH1D("hnum_4region_error2"+label+"_"+self.channel,"hnum_4region_error2"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total

        hnum_4region_before_cut = TH1D("hnum_4region_before_cut"+label+"_"+self.channel,"hnum_4region_before_cut"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
        hnum_4region_before_cut_error2 = TH1D("hnum_4region_before_cut_error2"+label+"_"+self.channel,"hnum_4region_before_cut_error2"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total

        hnum_2region = TH1D("hnum_2region"+label+"_"+self.channel,"hnum_2region"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total
        hnum_2region_error2 = TH1D("hnum_2region_error2"+label+"_"+self.channel,"hnum_2region_error2"+label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total

        rand = ROOT.TRandom3() ;
        rand.SetSeed();
        
        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0: print "iEntry: ",i
            treeIn.GetEntry(i);

            if i==0: tmp_scale_to_lumi = treeIn.wSampleWeight ## weigth for xs for the mc
                
            discriminantCut = 0;
            wtagger = getattr(treeIn,"ttb_jet_tau2tau1");
            
            if wtagger < options.tau2tau1cutHP:
                discriminantCut = 2;
            elif wtagger > options.tau2tau1cutHP and wtagger < options.tau2tau1cutLP:
                discriminantCut = 1;
            elif wtagger > options.tau2tau1cutLP :
                discriminantCut = 0;

            tmp_jet_mass = 0. ;    
            if options.shift == 1 and not TString(label).Contains("data"):
                 tmp_jet_mass = getattr(treeIn, jet_mass) + self.mean_shift;
            elif options.smear == 1 and not TString(label).Contains("data"):
                 tmp_jet_mass = getattr(treeIn, jet_mass)*(1+rand.Gaus(0,self.sigma_scale-1));
            else:
                 tmp_jet_mass = getattr(treeIn, jet_mass);

#            tmp_event_weight     = getattr(treeIn,"totalEventWeight");
            tmp_event_weight = getattr(treeIn,"wSampleWeight")*self.Lumi#*getattr(treeIn,"genWeight")*getattr(treeIn,"eff_and_pu_Weight")
            tmp_event_weight4fit = getattr(treeIn,"eff_and_pu_Weight")#*getattr(treeIn,"genWeight");
            if TString(label).Contains("TTbar"):
                tmp_event_weight4fit = getattr(treeIn,"eff_and_pu_Weight")#*getattr(treeIn,"genWeight");

#            tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"wSampleWeight")/tmp_scale_to_lumi

                 
            if not TString(label).Contains("data"):
#                  tmp_event_weight = tmp_event_weight*getattr(treeIn,"btag_weight");
#                  tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"btag_weight");                  
                  tmp_event_weight = tmp_event_weight;
                  tmp_event_weight4fit = tmp_event_weight4fit;                  
#                  tmp_event_weight = 1;
#                  tmp_event_weight4fit = 1;                  
            else:
                  tmp_event_weight = 1.;
                  tmp_event_weight4fit = 1.;
  
       

            ### Cut for the HP category
            if discriminantCut ==2 and getattr(treeIn,"mass_lvj_type0") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0") > self.mass_lvj_min and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ungroomed_jet_pt") > 100 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ungroomed_jet_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ungroomed_jet_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"nBTagJet_medium") > 0:


                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;

                rrv_mass_j.setVal(tmp_jet_mass);

                rdataset_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight);
                rdataset4fit_mj.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);
     
                if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                    hnum_2region.Fill(1,tmp_event_weight);
                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                    hnum_4region_error2.Fill(0,tmp_event_weight*tmp_event_weight);
                if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                    hnum_4region.Fill(1,tmp_event_weight);

                hnum_4region.Fill(2,tmp_event_weight);

                category_p_f.setLabel("pass");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

            ### Cut for the Total category
            if getattr(treeIn,"mass_lvj_type0") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0") > self.mass_lvj_min and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ungroomed_jet_pt") > 100 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ungroomed_jet_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ungroomed_jet_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"nBTagJet_medium") > 0:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;

                rrv_mass_j.setVal(tmp_jet_mass);

                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                   hnum_4region_before_cut.Fill(0,tmp_event_weight);
                   hnum_4region_before_cut_error2.Fill(0,tmp_event_weight*tmp_event_weight);

                rdataset_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight);
                rdataset4fit_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight4fit);
     
            ### 1-HP category
            if (discriminantCut==1 or discriminantCut==0) and getattr(treeIn,"mass_lvj_type0") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0") > self.mass_lvj_min and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ungroomed_jet_pt") > 100 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ungroomed_jet_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ungroomed_jet_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"nBTagJet_medium") > 0:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
            
                rrv_mass_j.setVal(tmp_jet_mass);

                rdataset_failtau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight);
                rdataset4fit_failtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
      
                category_p_f.setLabel("fail");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

            ### extreme fail category
            if discriminantCut==1 and getattr(treeIn,"mass_lvj_type0") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0") > self.mass_lvj_min and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ungroomed_jet_pt") > 100 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ungroomed_jet_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ungroomed_jet_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"nBTagJet_medium") > 0: #to be changed with the one below!!!!
#            if discriminantCut==0 and getattr(treeIn,"mass_lvj_type0") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0") > self.mass_lvj_min and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ungroomed_jet_pt") > 100 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ungroomed_jet_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ungroomed_jet_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"nBTagJet_medium") > 0: 
           
                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;

                rdataset_extremefailtau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight);
                rdataset4fit_extremefailtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,rdataset_mj.sumEntries()/rdataset4fit_mj.sumEntries());
        rrv_scale_to_lumi_failtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,rdataset_failtau2tau1cut_mj.sumEntries()/rdataset4fit_failtau2tau1cut_mj.sumEntries());
        rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,rdataset_extremefailtau2tau1cut_mj.sumEntries()/rdataset4fit_extremefailtau2tau1cut_mj.sumEntries());

        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi);
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut);
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut);
                                                            
        #prepare m_j dataset
        rrv_number_dataset_sb_lo_mj         = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_signal_region_error2_mj = RooRealVar("rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj",hnum_4region_error2.GetBinContent(2));

        rrv_number_dataset_signal_region_before_cut_mj        = RooRealVar("rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj",hnum_4region_before_cut.GetBinContent(2));
        rrv_number_dataset_signal_region_before_cut_error2_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj",hnum_4region_before_cut_error2.GetBinContent(2));
        rrv_number_dataset_sb_hi_mj                           = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));

        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj);
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj);
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_error2_mj);
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_mj);
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_error2_mj);
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj);
        getattr(self.workspace4fit_,"import")(combData_p_f);
        
        print "N_rdataset_mj: "
        getattr(self.workspace4fit_,"import")(rdataset_mj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj)
        getattr(self.workspace4fit_,"import")(rdataset_beforetau2tau1cut_mj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_beforetau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset_failtau2tau1cut_mj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_failtau2tau1cut_mj);
        getattr(self.workspace4fit_,"import")(rdataset_extremefailtau2tau1cut_mj);
        getattr(self.workspace4fit_,"import")(rdataset4fit_extremefailtau2tau1cut_mj)
                  
        rdataset_mj.Print();
        rdataset4fit_mj.Print();
        rdataset_failtau2tau1cut_mj.Print();
        rdataset4fit_failtau2tau1cut_mj.Print();
        rdataset_extremefailtau2tau1cut_mj.Print();
        rdataset4fit_extremefailtau2tau1cut_mj.Print();
        rrv_number_dataset_sb_lo_mj.Print();
        rrv_number_dataset_signal_region_mj.Print();
        rrv_number_dataset_signal_region_error2_mj.Print();
        rrv_number_dataset_signal_region_before_cut_mj.Print();
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print();
        rrv_number_dataset_sb_hi_mj.Print();

        rdataset_mj.Print();
        rdataset_beforetau2tau1cut_mj.Print();
        rdataset_failtau2tau1cut_mj.Print();
        rdataset_extremefailtau2tau1cut_mj.Print();
        rrv_number_dataset_signal_region_mj.Print();
        rrv_number_dataset_signal_region_error2_mj.Print();
        rrv_number_dataset_signal_region_before_cut_mj.Print();
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print();
        combData_p_f.Print("v");


    #### defines two different way to fit depending on pythia or herwig analysis
    def fit_TTbar_controlsample(self, isherwig=0, ttbarMC=0):

      if isherwig ==0 :

        print "fit_TTbar_controlsample --> Pythia samples"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");

        print "##################################################"
        print "############### Single Top DataSet ###############"
        print "##################################################"
          
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop");

        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop",self.mj_shape["STop"],self.channel,self.wtagger_label,1);
        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_failtau2tau1cut",self.mj_shape["STop_fail"],self.channel,self.wtagger_label,1);                      
        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_extremefailtau2tau1cut",self.mj_shape["STop_extremefail"],self.channel,self.wtagger_label,1);          

        ### Build WJet fit pass and fail distributions
        print "###########################################"
        print "############### WJets Pythia ##############"
        print "###########################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets0_mc,"_WJets0");

        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1);
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_failtau2tau1cut",self.mj_shape["WJets0_fail"],self.channel,self.wtagger_label,1);
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_extremefailtau2tau1cut",self.mj_shape["WJets0_extremefail"],self.channel,self.wtagger_label,1);


        ### Build VV fit pass and fail distributions
        print "#########################################"
        print "############### VV Pythia ###############"
        print "#########################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV");
        
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label,1);
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_failtau2tau1cut",self.mj_shape["VV_fail"],self.channel,self.wtagger_label,1);
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_extremefailtau2tau1cut",self.mj_shape["VV_extremefail"],self.channel,self.wtagger_label,1);

        print "#########################################"
        print "############# TTbar Powheg ##############"
        print "#########################################"
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_TTbar_mc,"_TTbar");

        print "################################################"
        print "############## Pseudo Data Powheg ##############"
        print "################################################"
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_pseudodata,"_TotalMC");


        print "#################################"
        print "############# Data ##############"
        print "#################################"
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_data,"_data");
        self.fit_mj_TTbar_controlsample(self.file_data);

        self.constrainslist_data = ROOT.std.vector(ROOT.std.string)();
        self.constrainslist_mc   = ROOT.std.vector(ROOT.std.string)();
        
        ScaleFactorTTbarControlSampleFit(self.workspace4fit_,self.mj_shape,self.color_palet,self.constrainslist_data,self.constrainslist_mc,"",self.channel,self.wtagger_label,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max);

      else:


          print "###############################################"
          print "############# Single Top DataSet ##############"
          print "###############################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop_herwig");
          fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_herwig",self.mj_shape["STop"],self.channel,self.wtagger_label,1);
          fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_herwig_failtau2tau1cut",self.mj_shape["STop_fail"],self.channel,self.wtagger_label,1);          
          fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_extremefailtau2tau1cut",self.mj_shape["STop_extremefail"],self.channel,self.wtagger_label,1);          
          
          print "##########################################"
          print "############## WJets Herwig ##############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets1_mc,"_WJets0_herwig");

          fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_herwig",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1);
          fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_herwig_failtau2tau1cut",self.mj_shape["WJets0_fail"],self.channel,self.wtagger_label,1);
          fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_herwig_extremefailtau2tau1cut",self.mj_shape["WJets0_extremefail"],self.channel,self.wtagger_label,1);


          print "##########################################"
          print "################ VV Pythia ###############"
          print "##########################################"
          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV_herwig");

          fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_herwig",self.mj_shape["VV"],self.channel,self.wtagger_label,1);
          fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_herwig_failtau2tau1cut",self.mj_shape["VV_fail"],self.channel,self.wtagger_label,1);
          fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_herwig_extremefailtau2tau1cut",self.mj_shape["VV_extremefail"],self.channel,self.wtagger_label,1);

          self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp

          print "##########################################"
          print "################ TTbar Herwig ############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_TTbar_herwig,"_TTbar_herwig");
 
          print "##############################################"
          print "############# Pseudo Data Herwig #############"
          print "##############################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_pseudodata_herwig,"_TotalMC_herwig");

          print "##################################"
          print "############## Data ##############"
          print "##################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_data,"_data_herwig");
          self.fit_mj_TTbar_controlsample(self.file_data,"_herwig");
          self.constrainslist_data = ROOT.std.vector(ROOT.std.string)();
          self.constrainslist_mc   = ROOT.std.vector(ROOT.std.string)();
          ScaleFactorTTbarControlSampleFit(self.workspace4fit_,self.mj_shape,self.color_palet,self.constrainslist_data,self.constrainslist_mc,"",self.channel,self.wtagger_label,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max);


class doFit_wj_and_wlvj_simultaneous:

    def __init__(self,isherwig=0):

        label = "";
        if isherwig==1 : label = "_herwig" ;

        self.workspace4fit_ = RooWorkspace("workspace4fit"+label+"_","workspace4fit"+label+"_"); ## create the workspace

        self.boostedW_fitter_el = doFit_wj_and_wlvj("el","ggH600",40,150,label, self.workspace4fit_); ## single object analysis for electrons
        self.boostedW_fitter_mu = doFit_wj_and_wlvj("mu","ggH600",40,150,label, self.workspace4fit_); ## single object analysis for muons
        self.boostedW_fitter_em = doFit_wj_and_wlvj("em","ggH600",40,150,label, self.workspace4fit_); ## single object analysis for muons

        self.boostedW_fitter_el.fit_TTbar_controlsample(isherwig); ## run the electron analysis
        self.boostedW_fitter_mu.fit_TTbar_controlsample(isherwig); ## run the muon analysis
        self.boostedW_fitter_em.fit_TTbar_controlsample(isherwig); ## run the muon analysis

        self.workspace4fit_.data("rdataset_data"+label+"_mu_mj").Print();
        self.workspace4fit_.data("rdataset_data"+label+"_el_mj").Print();
        self.workspace4fit_.data("rdataset_data"+label+"_em_mj").Print();

        #### Define simultaneusly 4 category
        sample_type = RooCategory("sample_type"+label,"sample_type"+label);
        sample_type.defineType("mu_pass");
        sample_type.defineType("mu_fail");
        sample_type.defineType("el_pass");
        sample_type.defineType("el_fail");
        sample_type.defineType("em_pass");
        sample_type.defineType("em_fail");
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 

        #### take the datasets 
        rdataset_data_mu_mj      = self.workspace4fit_.data("rdataset_data"+label+"_mu_mj");
        rdataset_data_el_mj      = self.workspace4fit_.data("rdataset_data"+label+"_el_mj");
        rdataset_data_em_mj      = self.workspace4fit_.data("rdataset_data"+label+"_em_mj");
        rdataset_data_mu_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_data_el_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_el_mj"); 
        rdataset_data_em_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_em_mj"); 
 
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j");

        ## combined dataset fill
#        combData_data = RooDataSet("combData_data"+label,"combData_data"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_data_mu_mj),RooFit.Import("el_pass",rdataset_data_el_mj),RooFit.Import("mu_fail",rdataset_data_mu_mj_fail),RooFit.Import("el_fail",rdataset_data_el_mj_fail) );
        combData_data = RooDataSet("combData_data"+label,"combData_data"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_data_mu_mj),RooFit.Import("el_pass",rdataset_data_el_mj),RooFit.Import("em_pass",rdataset_data_em_mj),RooFit.Import("mu_fail",rdataset_data_mu_mj_fail),RooFit.Import("el_fail",rdataset_data_el_mj_fail),RooFit.Import("em_fail",rdataset_data_em_mj_fail) );
        combData_data.Print();

        rdataset_TotalMC_mu_mj      = self.workspace4fit_.data("rdataset_TotalMC"+label+"_mu_mj");
        rdataset_TotalMC_el_mj      = self.workspace4fit_.data("rdataset_TotalMC"+label+"_el_mj");
        rdataset_TotalMC_em_mj      = self.workspace4fit_.data("rdataset_TotalMC"+label+"_em_mj");
        rdataset_TotalMC_mu_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_TotalMC_el_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_el_mj"); 
        rdataset_TotalMC_em_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_em_mj"); 

#        combData_TotalMC = RooDataSet("combData_TotalMC"+label,"combData_TotalMC"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_TotalMC_mu_mj),RooFit.Import("el_pass",rdataset_TotalMC_el_mj),RooFit.Import("mu_fail",rdataset_TotalMC_mu_mj_fail),RooFit.Import("el_fail",rdataset_TotalMC_el_mj_fail) );
        combData_TotalMC = RooDataSet("combData_TotalMC"+label,"combData_TotalMC"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_TotalMC_mu_mj),RooFit.Import("el_pass",rdataset_TotalMC_el_mj),RooFit.Import("em_pass",rdataset_TotalMC_em_mj),RooFit.Import("mu_fail",rdataset_TotalMC_mu_mj_fail),RooFit.Import("el_fail",rdataset_TotalMC_el_mj_fail),RooFit.Import("em_fail",rdataset_TotalMC_em_mj_fail) );
        combData_TotalMC.Print();

        # fit data --> import the pdf from the single fits and define the simultaneous total pdf
        model_data_mu      = self.workspace4fit_.pdf("model_data"+label+"_mu");
        model_data_fail_mu = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_mu");
        model_data_el      = self.workspace4fit_.pdf("model_data"+label+"_el");
        model_data_fail_el = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_el");
        model_data_em      = self.workspace4fit_.pdf("model_data"+label+"_em");
        model_data_fail_em = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_em");

        simPdf_data = RooSimultaneous("simPdf_data_em"+label,"simPdf_data_em"+label,sample_type);
        simPdf_data.addPdf(model_data_mu,"mu_pass");
        simPdf_data.addPdf(model_data_el,"el_pass");
        simPdf_data.addPdf(model_data_fail_mu,"mu_fail");
        simPdf_data.addPdf(model_data_fail_el,"el_fail");

        constrainslist_data_em = ROOT.std.vector(ROOT.std.string)();
        for i in range(self.boostedW_fitter_el.constrainslist_data.size()):
            constrainslist_data_em.push_back(self.boostedW_fitter_el.constrainslist_data.at(i));
        for i in range(self.boostedW_fitter_mu.constrainslist_data.size()):
            constrainslist_data_em.push_back(self.boostedW_fitter_mu.constrainslist_data.at(i));
            
        pdfconstrainslist_data_em = RooArgSet("pdfconstrainslist_data_em"+label);
        for i in range(constrainslist_data_em.size()):
            self.workspace4fit_.pdf(constrainslist_data_em.at(i)).Print();
            pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em.at(i)) );
        pdfconstrainslist_data_em.Print();

        rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em));
        rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em));

        # fit TotalMC --> define the simultaneous total pdf
        model_TotalMC_mu      = self.workspace4fit_.pdf("model_TotalMC"+label+"_mu");
        model_TotalMC_fail_mu = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_mu");
        model_TotalMC_el      = self.workspace4fit_.pdf("model_TotalMC"+label+"_el");
        model_TotalMC_fail_el = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_el");

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em"+label,"simPdf_TotalMC_em"+label,sample_type);
        simPdf_TotalMC.addPdf(model_TotalMC_mu,"mu_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_el,"el_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_mu,"mu_fail");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_el,"el_fail");

        constrainslist_TotalMC_em = ROOT.std.vector(ROOT.std.string)();
        for i in range(self.boostedW_fitter_el.constrainslist_mc.size()):
            constrainslist_TotalMC_em.push_back(self.boostedW_fitter_el.constrainslist_mc.at(i));
        for i in range(self.boostedW_fitter_mu.constrainslist_mc.size()):
            constrainslist_TotalMC_em.push_back(self.boostedW_fitter_mu.constrainslist_mc.at(i));

        pdfconstrainslist_TotalMC_em = RooArgSet("pdfconstrainslist_TotalMC_em"+label);
        for i in range(constrainslist_TotalMC_em.size()):
            self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]).Print();
            pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) );
        pdfconstrainslist_TotalMC_em.Print();

        rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em));
        rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em));

        ## draw the plots
        DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_mu.color_palet,label,"mu",self.boostedW_fitter_mu.wtagger_label,self.boostedW_fitter_mu.ca8_ungroomed_pt_min,self.boostedW_fitter_mu.ca8_ungroomed_pt_max);
        DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_el.color_palet,label,"el",self.boostedW_fitter_el.wtagger_label,self.boostedW_fitter_el.ca8_ungroomed_pt_min,self.boostedW_fitter_el.ca8_ungroomed_pt_max);
        DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_em.color_palet,label,"em",self.boostedW_fitter_em.wtagger_label,self.boostedW_fitter_em.ca8_ungroomed_pt_min,self.boostedW_fitter_em.ca8_ungroomed_pt_max);

        rfresult_TotalMC.Print();
        rfresult_data.Print();

        ### take the efficienty value from the fits in both data nad mc         
        rrv_eff_MC_el   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_el_mj");
        rrv_eff_MC_mu   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_mu_mj");
        rrv_eff_MC_em   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_em_mj");
        rrv_mean_MC_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_el_mj");
        rrv_sigma_MC_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_el_mj");
        rrv_mean_MC_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_em_mj");
        rrv_sigma_MC_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_em_mj");
        rrv_mean_MC_mu  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_mu_mj");
        rrv_sigma_MC_mu = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_mu_mj");

        rrv_eff_data_el   = self.workspace4fit_.var("eff_ttbar_data"+label+"_el_mj");
        rrv_eff_data_em   = self.workspace4fit_.var("eff_ttbar_data"+label+"_em_mj");
        rrv_eff_data_mu   = self.workspace4fit_.var("eff_ttbar_data"+label+"_mu_mj");
        
        rrv_mean_data_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_el_mj");
        rrv_sigma_data_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_el_mj");
        rrv_mean_data_mu  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_mu_mj");
        rrv_sigma_data_mu = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_mu_mj");
        rrv_mean_data_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_em_mj");
        rrv_sigma_data_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_em_mj");

        rrv_eff_MC_el.Print()   ; rrv_eff_MC_mu.Print();
        rrv_eff_data_el.Print() ; rrv_eff_data_mu.Print();
        rrv_mean_MC_el.Print()  ; rrv_mean_data_el.Print();  
        rrv_sigma_MC_el.Print() ; rrv_sigma_data_el.Print();  
        rrv_mean_MC_mu.Print()  ; rrv_mean_data_mu.Print();  
        rrv_sigma_MC_mu.Print() ; rrv_sigma_data_mu.Print();  
        rrv_mean_MC_em.Print()  ; rrv_mean_data_em.Print();  
        rrv_sigma_MC_em.Print() ; rrv_sigma_data_em.Print();  

        ## compute the scale factors and uncertainty for HP
        pure_wtagger_sf_el = rrv_eff_data_el.getVal()/rrv_eff_MC_el.getVal(); 
        pure_wtagger_sf_mu = rrv_eff_data_mu.getVal()/rrv_eff_MC_mu.getVal(); 
        pure_wtagger_sf_em = rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal(); 
        pure_wtagger_mean_shift_el    = rrv_mean_data_el.getVal()-rrv_mean_MC_el.getVal();
        pure_wtagger_sigma_enlarge_el = rrv_sigma_data_el.getVal()/rrv_sigma_MC_el.getVal();
        pure_wtagger_mean_shift_mu    = rrv_mean_data_mu.getVal()-rrv_mean_MC_mu.getVal();
        pure_wtagger_sigma_enlarge_mu = rrv_sigma_data_mu.getVal()/rrv_sigma_MC_mu.getVal();
        pure_wtagger_mean_shift_em    = rrv_mean_data_em.getVal()-rrv_mean_MC_em.getVal();
        pure_wtagger_sigma_enlarge_em = rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal();

        pure_wtagger_sf_el_err = ((rrv_eff_data_el.getError()/rrv_eff_data_el.getVal())**2 + (rrv_eff_MC_el.getError()/rrv_eff_MC_el.getVal())**2 )**0.5* pure_wtagger_sf_el;

        pure_wtagger_sf_mu_err = ((rrv_eff_data_mu.getError()/rrv_eff_data_mu.getVal())**2 + (rrv_eff_MC_mu.getError()/rrv_eff_MC_mu.getVal())**2 )**0.5* pure_wtagger_sf_mu;

        pure_wtagger_sf_em_err = ((rrv_eff_data_em.getError()/rrv_eff_data_em.getVal())**2 + (rrv_eff_MC_em.getError()/rrv_eff_MC_em.getVal())**2 )**0.5* pure_wtagger_sf_em;
        
        pure_wtagger_mean_shift_err_el = (rrv_mean_data_el.getError()**2 + rrv_mean_MC_el.getError()**2)**0.5;

        pure_wtagger_sigma_enlarge_err_el = ((rrv_sigma_data_el.getError()/rrv_sigma_data_el.getVal())**2 + (rrv_sigma_MC_el.getError()/rrv_sigma_MC_el.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_el;

        pure_wtagger_mean_shift_err_mu = (rrv_mean_data_mu.getError()**2 + rrv_mean_MC_mu.getError()**2)**0.5;

        pure_wtagger_sigma_enlarge_err_mu = ((rrv_sigma_data_mu.getError()/rrv_sigma_data_mu.getVal())**2 + (rrv_sigma_MC_mu.getError()/rrv_sigma_MC_mu.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_mu;

        pure_wtagger_mean_shift_err_em = (rrv_mean_data_em.getError()**2 + rrv_mean_MC_em.getError()**2)**0.5;

        pure_wtagger_sigma_enlarge_err_em = ((rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 + (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_em;

        print "Pure W-tagger SF of el %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err);
        print "Pure W-tagger SF of mu %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err);
        print "Pure W-tagger SF of em %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em, pure_wtagger_sf_em_err);
        print "Pure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el);
        print "Pure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el);
        print "Pure W-tagger mean shift mu %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_mu, pure_wtagger_mean_shift_err_mu);
        print "Pure W-tagger sigma enlarge mu %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_mu, pure_wtagger_sigma_enlarge_err_mu);
        print "Pure W-tagger mean shift em %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em);
        print "Pure W-tagger sigma enlarge em %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em);

        self.boostedW_fitter_el.file_out_ttbar_control.write( "\n***************************************************" )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of el %s     : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of mu %s        : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger mean shift mu %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_mu, pure_wtagger_mean_shift_err_mu));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge mu %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_mu, pure_wtagger_sigma_enlarge_err_mu));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of em %s     : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em, pure_wtagger_sf_em_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger mean shift em %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em));
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge em %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em));

        ## take the extreme fail informations
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.Print();
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_em_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_em_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.Print();

        rrv_number_total_ttbar_TotalMC_el = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_el_mj");
        rrv_number_total_ttbar_TotalMC_mu = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_mu_mj");
        rrv_number_total_ttbar_data_el = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_el_mj");
        rrv_number_total_ttbar_data_mu = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_mu_mj");
        rrv_number_total_ttbar_TotalMC_em = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_em_mj");
        rrv_number_total_ttbar_data_em = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_em_mj");
###
        rrv_number_total_ttbar_TotalMC_el.Print();
        rrv_number_total_ttbar_TotalMC_mu.Print();
        rrv_number_total_ttbar_TotalMC_em.Print();
        rrv_number_total_ttbar_data_el.Print();
        rrv_number_total_ttbar_data_mu.Print();
        rrv_number_total_ttbar_data_em.Print();

        print "el TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal());
        print "mu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal());
        print "el data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal());
        print "mu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal());
        print "em TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal());
        print "em data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_data_em.getVal());
        
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal()));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal()));

        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal()));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal()));

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal()));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_data_em.getVal()));

        tmp_eff_MC_el_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal();
        tmp_eff_MC_mu_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal();
        tmp_eff_data_el_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal();
        tmp_eff_data_mu_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal();
        tmp_eff_MC_em_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal();
        tmp_eff_data_em_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_data_em.getVal();

        tmp_eff_MC_el_extremefail_error =tmp_eff_MC_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_el.getError()/rrv_number_total_ttbar_TotalMC_el.getVal() )**2 );
        tmp_eff_MC_mu_extremefail_error =tmp_eff_MC_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_mu.getError()/rrv_number_total_ttbar_TotalMC_mu.getVal() )**2 );
        tmp_eff_data_el_extremefail_error =tmp_eff_data_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_data_el.getError()/rrv_number_total_ttbar_data_el.getVal() )**2 );
        tmp_eff_data_mu_extremefail_error =tmp_eff_data_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_data_mu.getError()/rrv_number_total_ttbar_data_mu.getVal() )**2 );
        tmp_eff_MC_em_extremefail_error =tmp_eff_MC_em_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_em.getError()/rrv_number_total_ttbar_TotalMC_em.getVal() )**2 );
        tmp_eff_data_em_extremefail_error =tmp_eff_data_em_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal() )**2+ (rrv_number_total_ttbar_data_em.getError()/rrv_number_total_ttbar_data_em.getVal() )**2 );

        print "eff_MC_el_extremefail_error %s: %f"%(label,tmp_eff_MC_el_extremefail_error);
        print "eff_MC_mu_extremefail_error %s: %f"%(label,tmp_eff_MC_mu_extremefail_error);
        print "eff_data_el_extremefail_error %s: %f"%(label,tmp_eff_data_el_extremefail_error);
        print "eff_data_mu_extremefail_error %s: %f"%(label,tmp_eff_data_mu_extremefail_error);
        print "eff_MC_em_extremefail_error %s: %f"%(label,tmp_eff_MC_em_extremefail_error);
        print "eff_data_em_extremefail_error %s: %f"%(label,tmp_eff_data_em_extremefail_error);

        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_MC_el_extremefail_error %s: %f"%(label,tmp_eff_MC_el_extremefail_error));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_data_el_extremefail_error %s: %f"%(label,tmp_eff_data_el_extremefail_error));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_MC_mu_extremefail_error %s: %f"%(label,tmp_eff_MC_mu_extremefail_error));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_data_mu_extremefail_error %s: %f"%(label,tmp_eff_data_mu_extremefail_error));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_MC_em_extremefail_error %s: %f"%(label,tmp_eff_MC_em_extremefail_error));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_data_em_extremefail_error %s: %f"%(label,tmp_eff_data_em_extremefail_error));
        
        ## Low purity scale factors LP
        tmp_eff_MC_el_LP = 1.-rrv_eff_MC_el.getVal()-tmp_eff_MC_el_extremefail;
        tmp_eff_MC_mu_LP = 1.-rrv_eff_MC_mu.getVal()-tmp_eff_MC_mu_extremefail;
        tmp_eff_data_el_LP =1.-rrv_eff_data_el.getVal()-tmp_eff_data_el_extremefail;
        tmp_eff_data_mu_LP =1.-rrv_eff_data_mu.getVal()-tmp_eff_data_mu_extremefail;
        tmp_eff_MC_em_LP = 1.-rrv_eff_MC_em.getVal()-tmp_eff_MC_em_extremefail;
        tmp_eff_data_em_LP =1.-rrv_eff_data_em.getVal()-tmp_eff_data_em_extremefail;

        tmp_eff_MC_el_LP_err = TMath.Sqrt( rrv_eff_MC_el.getError()**2 + tmp_eff_MC_el_extremefail_error**2 );
        tmp_eff_MC_mu_LP_err = TMath.Sqrt( rrv_eff_MC_mu.getError()**2 + tmp_eff_MC_mu_extremefail_error**2 );
        tmp_eff_data_el_LP_err = TMath.Sqrt( rrv_eff_data_el.getError()**2 + tmp_eff_data_el_extremefail_error**2 );
        tmp_eff_data_mu_LP_err = TMath.Sqrt( rrv_eff_data_mu.getError()**2 + tmp_eff_data_mu_extremefail_error**2 );
        tmp_eff_MC_em_LP_err = TMath.Sqrt( rrv_eff_MC_em.getError()**2 + tmp_eff_MC_em_extremefail_error**2 );
        tmp_eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2 + tmp_eff_data_em_extremefail_error**2 );

        print "LP Eff of el data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err);
        print "LP Eff of el MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err);
        print "LP Eff of mu data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err);
        print "LP Eff of mu MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err);
        print "LP Eff of em data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_em_LP, tmp_eff_data_em_LP_err);
        print "LP Eff of em MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_em_LP, tmp_eff_MC_em_LP_err);

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el data %s: %f +/- %f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el MC %s: %f +/- %f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu data %s: %f +/- %f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu MC %s: %f +/- %f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of em data %s: %f +/- %f"%(label,tmp_eff_data_em_LP, tmp_eff_data_em_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of em MC %s: %f +/- %f"%(label,tmp_eff_MC_em_LP, tmp_eff_MC_em_LP_err));
        
        pure_wtagger_sf_el_LP = tmp_eff_data_el_LP / tmp_eff_MC_el_LP;
        pure_wtagger_sf_mu_LP = tmp_eff_data_mu_LP / tmp_eff_MC_mu_LP;
        pure_wtagger_sf_el_LP_err = pure_wtagger_sf_el_LP*TMath.Sqrt( (tmp_eff_data_el_LP_err/tmp_eff_data_el_LP)**2 + (tmp_eff_MC_el_LP_err/tmp_eff_MC_el_LP)**2 );
        pure_wtagger_sf_mu_LP_err = pure_wtagger_sf_mu_LP*TMath.Sqrt( (tmp_eff_data_mu_LP_err/tmp_eff_data_mu_LP)**2 + (tmp_eff_MC_mu_LP_err/tmp_eff_MC_mu_LP)**2 );
        pure_wtagger_sf_em_LP = tmp_eff_data_em_LP / tmp_eff_MC_em_LP;
        pure_wtagger_sf_em_LP_err = pure_wtagger_sf_em_LP*TMath.Sqrt( (tmp_eff_data_em_LP_err/tmp_eff_data_em_LP)**2 + (tmp_eff_MC_em_LP_err/tmp_eff_MC_em_LP)**2 );

        print "Pure W-tagger LP SF of el %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err);
        print "Pure W-tagger LP SF of mu %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err);
        print "Pure W-tagger LP SF of em %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em_LP, pure_wtagger_sf_em_LP_err);

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nPure W-tagger LP SF of el %s: %f +/- %f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nPure W-tagger LP SF of mu %s: %f +/- %f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nPure W-tagger LP SF of em %s: %f +/- %f"%(label,pure_wtagger_sf_em_LP, pure_wtagger_sf_em_LP_err));

### function to call single channel fits
def control_sample(channel="mu",isherwig=0, ttbarMC=0):

    print "control sample "+channel;
    if isherwig == 0:
        ## create the object and do the single for pythia
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,150);
        boostedW_fitter.fit_TTbar_controlsample(isherwig);

    elif isherwig == 1:
        ## create the object and do the single for herwig
        boostedW_fitter_herwig = doFit_wj_and_wlvj(channel,"ggH600",40,150,"_herwig");
        boostedW_fitter_herwig.fit_TTbar_controlsample(isherwig);

    elif isherwig == 2:
        ## do Pythia and herwig analysis at the same time
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,150);
        boostedW_fitter.fit_TTbar_controlsample(0);
        boostedW_fitter_herwig = doFit_wj_and_wlvj(channel,"ggH600",40,150,"_herwig");
        boostedW_fitter_herwig.fit_TTbar_controlsample(1);
                                                             
### function to call simultaneous channel fits
def control_sample_simultaneous():

    print "control_sample_simultaneous";
    if options.herwig == 0 :
        boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous(options.herwig);
    elif options.herwig == 1 :
        boostedW_fitter_sim_herwig = doFit_wj_and_wlvj_simultaneous(options.herwig);
    elif options.herwig == 2:
        boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous();
        boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous(options.herwig);


### main code
if __name__ == '__main__':

    channel = options.channel; 
    if options.fitwtaggersim:
        print 'fitwtagger for el+mu sample'
        control_sample_simultaneous();

    elif options.fitwtagger:
        print 'fitwtagger for %s sample'%(channel)
        control_sample(channel,options.herwig);
         
