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
parser.add_option('-k','--ttbarcontrolregion',  help='run the toy or f-test analysis in the ttbar control region', type=int, default=0)
parser.add_option('-y','--mlvjregion',  help='run the toy taking the events in the : low sb, high sb or signal region', type="string", default="_sb_lo")
parser.add_option('-q','--fitjetmass',  help='run the toy on the jet mass fit', type=int, default=0)
parser.add_option('-w','--onlybackgroundfit',  help='run only background fit',  type=int, default=0)
parser.add_option('-i','--inflatejobstatistic',  help='enlarge the generated statistics in the fit',  type=int, default=1)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/BiasStudy/BiasUtils_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error, RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooPow3Pdf, RooErfPow3Pdf, RooUser1Pdf, biasModelAnalysis

from ROOT import setTDRStyle, get_pull, draw_canvas, draw_canvas_with_pull, legend4Plot, GetDataPoissonInterval

class doBiasStudy_mlvj:

    def __init__(self,in_channel,in_ggH_sample,in_mlvj_min=400., in_mlvj_max=1400., in_mj_min=40, in_mj_max=130, generation_model="ErfExp", fit_model="ErfExp", input_workspace=None):

        setTDRStyle(); #import the tdr style

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;
        ROOT.RooRandom.randomGenerator().SetSeed(0);

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

        suffix = "";
        if options.ttbarcontrolregion : suffix = "_ttbar";
        if options.fitjetmass         : suffix = "_jetmass";
        if in_mlvj_min < 550          : suffix = suffix+"_turnOn";                
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
        self.color_palet={ 'data' : 1, 'WJets' : 2, 'VV' : 4, 'WW_EWK' : 6, 'STop' : 7, 'TTbar' : 210, 'ggH' : 1,
                           'vbfH' : 12, 'Signal': 1, 'Uncertainty' : kBlack, 'Other_Backgrounds' : kBlue};

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
          self.file_FTestFile_txt = "FTest_%s_%s_mc.txt"%(self.channel,self.wtagger_label);
         else:
          self.file_FTestFile_txt = "FTest_%s_%s_mc_ttbar.txt"%(self.channel,self.wtagger_label);
         self.file_out_FTest = open(self.file_FTestFile_txt,"w");
        elif options.shapetest == 1 and options.isMC == 0 :
         if options.ttbarcontrolregion == 0 :   
          self.file_FTestFile_txt = "FTest_%s_%s_data.txt"%(self.channel,self.wtagger_label);
         else:
          self.file_FTestFile_txt = "FTest_%s_%s_data_ttbar.txt"%(self.channel,self.wtagger_label);         
         self.file_out_FTest = open(self.file_FTestFile_txt,"w");

        ### create output root file for the pull plot vs mass
        if options.shapetest  == 0:
            self.outputFile  = ROOT.TFile("output_%s_%s_%s%s.root"%(self.ggH_sample,options.fgen,options.fres,suffix),"RECREATE");
            self.outputTree  = ROOT.TTree("otree","otree");
            self.outputFile.cd(); 

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
	    
        ######################### models for mJ fit 
        
	if in_model_name == "2Gaus":
       
            mean1_tmp     =8.3141e+01;  mean1_tmp_err     =1.63e-01;
            deltamean_tmp =6.9129e+00;  deltamean_tmp_err =1.24e+00;
            sigma1_tmp    =7.5145e+00;  sigma1_tmp_err    =1.99e-01;
            scalesigma_tmp=3.6819e+00;  scalesigma_tmp_err=2.11e-01;
            frac_tmp      =6.7125e-01;  frac_tmp_err      =2.09e-02; 


            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel+mass_spectrum,"gaus1"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,deltamean_tmp, -4, deltamean_tmp+deltamean_tmp_err*4);
            rrv_mean2_gaus     = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel+mass_spectrum,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*8, scalesigma_tmp+scalesigma_tmp_err*8); 
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel+mass_spectrum,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel+mass_spectrum,"gaus2"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel+mass_spectrum,"rrv_frac"+label+"_"+self.channel+mass_spectrum,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)
 
        ###############	
	if in_model_name == "2_2Gaus":#for VV m_j

            mean1_tmp     =8.3141e+01;  mean1_tmp_err     =1.63e-01;
            deltamean_tmp =6.9129e+00;  deltamean_tmp_err =1.24e+00;
            sigma1_tmp    =7.5145e+00;  sigma1_tmp_err    =1.99e-01;
            scalesigma_tmp=3.6819e+00;  scalesigma_tmp_err=2.11e-01;
            frac_tmp      =6.7125e-01;  frac_tmp_err      =2.09e-02;
            
            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            
	    gaus1 = RooGaussian("gaus1"+label+"_"+self.channel+mass_spectrum,"gaus1"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,deltamean_tmp, deltamean_tmp-deltamean_tmp_err*4, deltamean_tmp+deltamean_tmp_err*4); 
            
	    rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel+mass_spectrum,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            	                
	    rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*10, scalesigma_tmp+scalesigma_tmp_err*10); 

	    rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel+mass_spectrum,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            
	    gaus2 = RooGaussian("gaus2"+label+"_"+self.channel+mass_spectrum,"gaus2"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_shift = RooRealVar("rrv_shift"+label+"_"+self.channel+mass_spectrum,"rrv_shift"+label+"_"+self.channel+mass_spectrum,10.8026)   # Z mass: 91.1876;  shift=91.1876-80.385=10.8026
            
            rrv_frac1  = RooRealVar("rrv_frac1"+label+"_"+self.channel+mass_spectrum,"rrv_frac1"+label+"_"+self.channel+mass_spectrum,frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            gausguas_1 = RooAddPdf("gausguas_1"+label+"_"+self.channel+mass_spectrum+mass_spectrum,"gausguas_1"+label+"_"+self.channel+mass_spectrum+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac1),1)

            rrv_mean3_gaus = RooFormulaVar("rrv_mean3_gaus"+label+"_"+self.channel+mass_spectrum,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_shift));
            rrv_mean4_gaus = RooFormulaVar("rrv_mean4_gaus"+label+"_"+self.channel+mass_spectrum,"@0+@1",RooArgList(rrv_mean2_gaus, rrv_shift));

            gaus3 = RooGaussian("gaus3"+label+"_"+self.channel+mass_spectrum,"gaus3"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean3_gaus,rrv_sigma1_gaus);
            gaus4 = RooGaussian("gaus4"+label+"_"+self.channel+mass_spectrum,"gaus4"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean4_gaus,rrv_sigma2_gaus);

            gausguas_2 = RooAddPdf("gausguas_2"+label+"_"+self.channel+mass_spectrum+mass_spectrum,"gausguas_2"+label+"_"+self.channel+mass_spectrum+mass_spectrum,RooArgList(gaus3,gaus4),RooArgList(rrv_frac1),1)

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel+mass_spectrum,"rrv_frac"+label+"_"+self.channel+mass_spectrum,0.74)#,0.5,1.0);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gausguas_1,gausguas_2),RooArgList(rrv_frac),1)
        ###############
        if in_model_name == "ErfExpGaus_sp":#offset == mean

            rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,-0.05,-0.2,0.);
            rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,30.,10.,200.);
            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,84,78,88);

            erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel+mass_spectrum,"erfExp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp);

            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,7,4,10);
            rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel+mass_spectrum,"rrv_high"+label+"_"+self.channel+mass_spectrum,0.5,0.,1.);

            gaus = RooGaussian("gaus"+label+"_"+self.channel+mass_spectrum,"gaus"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))

        ###############
        if in_model_name == "ErfExp" :
            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,-0.05,-0.5,-1e-4);
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum,60.,20.,200);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,30.)#,10, 80.);

            model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        if in_model_name == "ErfPow":#can replace erf*exp                                                                                                                           
            rrv_c      = RooRealVar("rrv_c_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfPow"+label+"_"+self.channel+mass_spectrum,-0.05,-15.,-1e-4);
            rrv_offset = RooRealVar("rrv_offset_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfPow"+label+"_"+self.channel+mass_spectrum,60.,20.,200);
            rrv_width  = RooRealVar("rrv_width_ErfPow"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfPow"+label+"_"+self.channel+mass_spectrum,30.)#,10, 80.);;

            model_pdf  = RooErfPowPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c,rrv_offset,rrv_width);

        ###############
        if in_model_name == "User1":
            rrv_p0 = RooRealVar("rrv_p0_User1"+label+"_"+self.channel+mass_spectrum,"rrv_p0_User1"+label+"_"+self.channel+mass_spectrum, 6, 0,100);
            rrv_p1 = RooRealVar("rrv_p1_User1"+label+"_"+self.channel+mass_spectrum,"rrv_p1_User1"+label+"_"+self.channel+mass_spectrum, -3,-30,0);

            model_pdf = RooUser1Pdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_p0,rrv_p1);

        ###############
        if in_model_name == "2Gaus_ErfExp":
                
	    mean1_tmp     =8.3141e+01;  mean1_tmp_err     =1.63e-01;
            deltamean_tmp =6.9129e+00;  deltamean_tmp_err =1.24e+00;
            sigma1_tmp    =7.5145e+00;  sigma1_tmp_err    =1.99e-01;
            scalesigma_tmp=3.6819e+00;  scalesigma_tmp_err=2.11e-01;
            frac_tmp      =6.7125e-01;  frac_tmp_err      =2.09e-02;
 
            rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_mean1_gaus"+label+"_"+self.channel+mass_spectrum,mean1_tmp, mean1_tmp-8, mean1_tmp+8);
            rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_sigma1_gaus"+label+"_"+self.channel+mass_spectrum,sigma1_tmp, sigma1_tmp-10,sigma1_tmp+10 );
            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel+mass_spectrum,"gaus1"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus  = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_deltamean_gaus"+label+"_"+self.channel+mass_spectrum,deltamean_tmp)#, deltamean_tmp, deltamean_tmp); 
            rrv_mean2_gaus      = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel+mass_spectrum,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,"rrv_scalesigma_gaus"+label+"_"+self.channel+mass_spectrum,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp); 
            rrv_sigma2_gaus     = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel+mass_spectrum,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel+mass_spectrum,"gaus2"+label+"_"+self.channel+mass_spectrum, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel+mass_spectrum,"rrv_frac_2gaus"+label+"_"+self.channel+mass_spectrum,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            c0_tmp    =   -2.9893e-02 ; c0_tmp_err     = 6.83e-03;
            offset_tmp=    7.9350e+01 ; offset_tmp_err = 9.35e+00;
            width_tmp =    3.3083e+01 ; width_tmp_err  = 2.97e+00; 

            rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_c_ErfExp"+label+"_"+self.channel+mass_spectrum,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2  );
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_offset_ErfExp"+label+"_"+self.channel+mass_spectrum, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum,"rrv_width_ErfExp"+label+"_"+self.channel+mass_spectrum, width_tmp, width_tmp-10, width_tmp+10);
            erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

            rrv_frac  = RooRealVar("rrv_frac"+label+"_"+self.channel+mass_spectrum,"rrv_frac"+label+"_"+self.channel+mass_spectrum, 0.5,0.,1.);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)

        ## 2Gaus per ttbar fit
        if in_model_name == "2Gaus_ttbar":

            mean1_tmp = 8.3141e+01; mean1_tmp_err = 1.63e-01;
            deltamean_tmp = 6.9129e+00; deltamean_tmp_err = 1.24e+00;
            sigma1_tmp = 7.5145e+00; sigma1_tmp_err = 1.99e-01;
            scalesigma_tmp = 3.6819e+00; scalesigma_tmp_err = 2.11e-01;
            frac_tmp = 6.7125e-01; frac_tmp_err = 2.09e-02;
          
            ## constrain the peak and the width among electron and muon channel to be the same in case of simultaneous fit
            if self.channel=="el":
                if self.workspace4bias_.var("rrv_mean1_gaus%s_mu"%(label)) and self.workspace4bias_.var("rrv_sigma1_gaus%s_mu"%(label)):
                    rrv_mean1_gaus = self.workspace4bias_.var("rrv_mean1_gaus%s_mu"%(label));
                    rrv_sigma1_gaus = self.workspace4bias_.var("rrv_sigma1_gaus%s_mu"%(label));
                else:
                    rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4);
            elif self.channel=="mu" :
                if self.workspace4bias_.var("rrv_mean1_gaus%s_el"%(label)) and self.workspace4bias_.var("rrv_sigma1_gaus%s_el"%(label)):
                    rrv_mean1_gaus = self.workspace4bias_.var("rrv_mean1_gaus%s_el"%(label));
                    rrv_sigma1_gaus = self.workspace4bias_.var("rrv_sigma1_gaus%s_el"%(label));
                else:
                    rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );

            gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

            rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#,deltamean_tmp-deltamean_tmp_err*5, deltamean_tmp+deltamean_tmp_err*5);
            rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
            rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#,scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
            gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

            rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp)#,frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

            model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)


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

 #### get a generic mj model from the workspace
    def get_mj_Model(self,label):
        print "model"+label+"_"+self.channel+"_mj"
        return self.workspace4bias_.pdf("model"+label+"_"+self.channel+"_mj");

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

    #### get a general mj model and fiz the paramters --> for extended pdf
    def get_General_mj_Model(self, label, model, fix = 1):
        print "########### Fixing amd return a general mj model ############"
        rdataset_General_mj = self.workspace4bias_.data("rdataset4bias%s_%s_mj"%(label,self.channel))
        model_General = self.get_mj_Model(label+model);
        rdataset_General_mj.Print();
        model_General.Print();
	if fix == 1 :
         parameters_General = model_General.getParameters(rdataset_General_mj);
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

    ###### get TTbar model mlvj in a region
    def get_TTbar_mj_Model(self, label ="_TTbar", model = "", fix = 1):
        print "########### Fixing TTbar mj model ############"
        return self.get_General_mj_Model(label,model,fix);

    ###### get Single Top model mlvj in a region
    def get_STop_mlvj_Model(self,  mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing Stop mlvj model for the region",mlvj_region," ############"
        return self.get_General_mlvj_Model("_STop",mlvj_region,model);

    ###### get Single Top model mj in a region
    def get_STop_mj_Model(self, label = "_STop", model ="", fix = 1):
        print "########### Fixing Stop mj model ############"
        return self.get_General_mj_Model(label,model,fix);

    ###### get ggH model mlvj in a region
    def get_ggH_signal_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_%s"%(self.ggH_sample),mlvj_region,model);

    ###### get ggH model mj in a region
    def get_ggH_signal_mj_Model(self, model = "2Gaus", fix = 1):
        print "########### Fixing VV mj model ############"
        return self.get_General_mj_Model("_%s"%(self.ggH_sample),model,fix);

    ###### get ggH model mlvj in a region
    def get_vbfH_signal_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_%s"%(self.vbfhiggs_sample),mlvj_region,model);

    ###### get ggH model mj in a region
    def get_vbfH_signal_mj_Model(self, model = "2Gaus", fix = 1):
        print "########### Fixing VV mj for the region ############"
        return self.get_General_mj_Model("_%s"%(self.vbfhiggs_sample),model,fix);

    ###### get VV mlvj in a region
    def get_VV_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_VV",mlvj_region,model);

    ###### get VV mlvj in a region
    def get_VV_mj_Model(self, label = "_VV", model = "", fix = 1):
        print "########### Fixing VV mj for the region model ############"
        return self.get_General_mj_Model(label,model,fix);

    ###### get VV mlvj in a region
    def get_WW_EWK_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing VV mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_WW_EWK",mlvj_region,model);

    ###### get VV mlvj in a region
    def get_WW_EWK_mj_Model(self, label = "_WW_EWK", model = "", fix = 1):
        print "########### Fixing WW_EWK mj for the region ############"
        return self.get_General_mj_Model(label,model,fix);

    ###### get WJets mlvj in a region
    def get_WJets_mlvj_Model(self, mlvj_region="_signal_region", model = "ErfExp_v1"):
        print "########### Fixing WJets mlvj for the region",mlvj_region," model",model,"############"
        return self.get_General_mlvj_Model("_WJets0",mlvj_region,model);

    ###### get WJets mj in a region
    def get_WJets_mj_Model(self, label = "_WJets0", model = "", fix = 1):
        print "########### Fixing WJets mj for the region ############"
        return self.get_General_mj_Model(label,model,fix);

    ### fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
    def fix_Model(self, label, mlvj_region="_signal_region",mass_spectrum="_mlvj",model="",additional_info="",notExtended=0):
              
        print "########### Fixing an Extended Pdf for mlvj ############"
	if mass_spectrum == "_mj" :  
         rdataset = self.workspace4bias_.data("rdataset4bias%s_%s%s"%(label,self.channel,mass_spectrum))
         if notExtended == 1 :
           label = "_pdf"+label;
         model = self.get_mj_Model(label+model+additional_info);             
        else:  
         rdataset = self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,mlvj_region,self.channel,mass_spectrum))
         if notExtended == 1:
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
	if mass_spectrum == "_mj" :  
         rdataset = self.workspace4bias_.data("rdataset4bias%s_%s%s"%(label,self.channel,mass_spectrum));
         model = self.get_mj_Model(label+model);
        else:
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
            
	    
        ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4bias_.var("rrv_mass_j");
        rdataset_mj = self.workspace4bias_.data("rdataset4bias"+label+"_"+self.channel+"_mj");
        rdataset_mj.Print();
        
        ## make the extended model
        if additioninformation == 1:
         model = self.make_Model(label+in_model_name,in_model_name);
        else:
         model = self.make_Model(label,in_model_name);
            
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE),   RooFit.SumW2Error(kTRUE));
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## Plot the result
        mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins())));
        rdataset_mj.plotOn( mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

	## draw the error band for an extend pdf
        draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,2,"L");

        ## re-draw the dataset
        rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0),RooFit.Name("data"));
        ## draw the function
        model.plotOn( mplot,RooFit.Name("model_mc")  );# remove RooFit.VLines() in order to get right pull in the 1st bin

        ## Get the pull
        mplot_pull = get_pull(rrv_mass_j,mplot,rdataset_mj,model,rfresult,"data","model_mc",0,1);
        mplot.GetYaxis().SetRangeUser(1e-5,mplot.GetMaximum()*1.2);

        ## CALCULATE CHI2
        datahist = rdataset_mj.binnedClone(rdataset_mj.GetName()+"_binnedClone",rdataset_mj.GetName()+"_binnedClone")
        Nbin = int(rrv_mass_j.getBins()); 
        rresult_param = rfresult.floatParsFinal();        
        nparameters =  rresult_param.getSize()                                         
        ChiSquare = model.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
        chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters);

        ## Add Chisquare to mplot_pull
        cs = TLatex(0.75,0.8,"#chi^{2}/ndf = %0.2f "%(float(chi_over_ndf)));
        cs.SetNDC();
        cs.SetTextSize(0.12);
        cs.AppendPad("same");
        mplot_pull.addObject(cs)

        parameters_list = model.getParameters(rdataset_mj);

        if not os.path.isdir("plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/")):
         os.system("mkdir -p plots_%s_%s_%s_g1/mj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/"));

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),label+in_file_name,in_model_name,"em",0,1,self.GetLumi());


        print "####################################################";
        print "######## Normalization Factor in mJ ################"
        print "####################################################";

        #normalize the number of total events to lumi --> correct the number to scale to the lumi
        if additioninformation == 1:
         self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setVal(self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getVal()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
         self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getError()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )

         self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").Print();
         
         if TString(label).Contains("ggH"):
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setVal( self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getVal() )
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getError() )
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").Print();

         if TString(label).Contains("vbfH"):
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setVal( self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getVal() )
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").getError() )
            self.workspace4bias_.var("rrv_number"+label+in_model_name+"_"+self.channel+"_mj").Print();
            
        else:
         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )
         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal() )

         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();
         
         if TString(label).Contains("ggH"):
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal() )
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError() )
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

         if TString(label).Contains("vbfH"):
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal() )
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError() )
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

            
        ##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV
        par=parameters_list.createIterator();
        par.Reset();
        param=par.Next()
        while (param):
              if (TString(label).Contains("VV") or TString(label).Contains("WW_EWK") or TString(label).Contains("STop") or TString(label).Contains("TTbar")):
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
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
        mplot_pull = get_pull(rrv_mass_lvj,mplot,rdataset,model,rfresult,"data","model_mc",0,1);
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

        if not os.path.isdir("plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/")):
         os.system("mkdir -p plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/"));

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),label+in_file_name,mlvj_model,"em",0,1,self.GetLumi());

        ### Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
        print "rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj";
        print "rrv_scale_to_lumi"+label+"_"+self.channel+in_range+"_mlvj";
        
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").getVal()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel+in_range+"_mlvj").getVal());
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").setError(self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").getError()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel+in_range+"_mlvj").getVal() );
        self.workspace4bias_.var("rrv_number"+label+in_range+mlvj_model+"_"+self.channel+"_mlvj").Print();
	
    #### method to fit the WJets normalization inside the mj signal region -> and write the jets mass sys if available
    def fit_WJetsNorm(self, label): # to get the normalization of WJets in signal_region

        print "############### Fit mj Normalization ##################"
        ## fit the two version of pdf for Wjets shape if available
        self.fit_WJetsNormalization_in_Mj_signal_region(label);
	 
        rrv_WJets0  = self.workspace4bias_.var("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel)); ## nominal parametrization for Wjets        
        rrv_WJets0.Print();
         
	rrv_STop  = self.workspace4bias_.var("rrv_number_dataset_signal_region_STop_%s_mj"%(self.channel));
        rrv_STop.Print();

        rrv_TTbar = self.workspace4bias_.var("rrv_number_dataset_signal_region_TTbar_%s_mj"%(self.channel));
        rrv_TTbar.Print();

        rrv_WJets  = self.workspace4bias_.var("rrv_number_dataset_signal_region_WJets0_%s_mj"%(self.channel));
        rrv_WJets.Print();
        rrv_WJets.Print();

        rrv_VV = self.workspace4bias_.var("rrv_number_dataset_signal_region_VV_%s_mj"%(self.channel));
        rrv_VV.Print();
        
        rrv_WW_EWK  = self.workspace4bias_.var("rrv_number_dataset_signal_region_WW_EWK_%s_mj"%(self.channel));
        rrv_WW_EWK.Print();
         
	rrv_ggH  = self.workspace4bias_.var("rrv_number_dataset_signal_region_%s_%s_mj"%(self.ggH_sample,self.channel))        
        rrv_ggH.Print();
         
	rrv_vbf  = self.workspace4bias_.var("rrv_number_dataset_signal_region_%s_%s_mj"%(self.vbfhiggs_sample,self.channel))
        rrv_vbf.Print();
    
    #### make the mj sideband fit on data ti get the Wjets normaliztion
    def fit_WJetsNormalization_in_Mj_signal_region(self,label):

        print "############### Fit mj Normalization: ",label," ##################"
	rrv_mass_j = self.workspace4bias_.var("rrv_mass_j")
	rdataset_data_mj = self.workspace4bias_.data("rdataset_data_%s_mj"%self.channel)

	### Fix TTbar, VV and STop
        if options.ttbarcontrolregion:
         model_WJets  = self.get_WJets_mj_Model("_WJets0");
         model_STop   = self.get_STop_mj_Model("_STop");
         model_VV     = self.get_VV_mj_Model("_VV");
         model_WW_EWK = self.get_WW_EWK_mj_Model("_WW_EWK");
	 model_TTbar  = self.get_TTbar_mj_Model(label,options.fgen,0);
        else :
         model_TTbar  = self.get_TTbar_mj_Model("_TTbar");
         model_STop   = self.get_STop_mj_Model("_STop");
         model_VV     = self.get_VV_mj_Model("_VV");
         model_WW_EWK = self.get_WW_EWK_mj_Model("_WW_EWK");
         model_WJets  = self.get_WJets_mj_Model(label,options.fgen,0);


        ## Total Pdf and fit only in sideband
        model_data = RooAddPdf("model_data_%s_mj"%(self.channel),"model_data_%s_mj"%(self.channel),RooArgList(model_WJets,model_VV,model_WW_EWK,model_TTbar,model_STop));
                                    
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(4) );
        rfresult = model_data.fitTo( rdataset_data_mj, RooFit.Save(1) , RooFit.Range("sb_lo,sb_hi") ,RooFit.Extended(kTRUE), RooFit.NumCPU(4), RooFit.Minimizer("Minuit2") );
        rfresult.Print();
	rfresult.covarianceMatrix().Print();
        getattr(self.workspace4bias_,"import")(model_data);

	## Total numver of event --> full propagation of error due to all the background sources coming from the fit
        if options.ttbarcontrolregion:
         rrv_number_data_mj = RooRealVar("rrv_number_data_%s_mj"%(self.channel),"rrv_number_data_%s_mj"%(self.channel),
                                          self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getVal()+ ## TTbar
                                          self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()+  ## STop
                                          self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()+    ## VV
                                          self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()+ ## WW_EWK
                                          self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal());  ## WJets

         rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getError()+       
                                                self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getError()));         
         getattr(self.workspace4bias_,"import")(rrv_number_data_mj);

         print "TTbar  events: ",self.workspace4bias_.var("rrv_number_%s_%s_mj"%(label+options.fgen,self.channel)).getVal();
         print "STop   events: ",self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal();
         print "VV     events: ",self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal();
         print "WW_EWK events: ",self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal();
         print "WJets  events: ",self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal();
         print "Data   events: ",self.workspace4bias_.var("rrv_number_data_%s_mj"%(self.channel)).getVal();

        else:
         rrv_number_data_mj = RooRealVar("rrv_number_data_%s_mj"%(self.channel),"rrv_number_data_%s_mj"%(self.channel),
                                          self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal()+ ## TTbar
                                          self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()+  ## STop
                                          self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()+    ## VV
                                          self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()+ ## WW_EWK
                                          self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getVal());  ## WJets

         rrv_number_data_mj.setError(TMath.Sqrt(self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getError()+
                                                self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getError()+       
                                                self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getError()*
                                                self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getError()));         
         getattr(self.workspace4bias_,"import")(rrv_number_data_mj);

         print "TTbar  events: ",self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal();
         print "STop   events: ",self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal();
         print "VV     events: ",self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal();
         print "WW_EWK events: ",self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal();
         print "WJets  events: ",self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getVal();
         print "Data   events: ",self.workspace4bias_.var("rrv_number_data_%s_mj"%(self.channel)).getVal();

        ## draw the plot for the default WJets Shape
        mplot = rrv_mass_j.frame(RooFit.Title("rrv_mass_j"), RooFit.Bins(int(rrv_mass_j.getBins())));
        rdataset_data_mj.plotOn(mplot, RooFit.Name("data_invisible"), RooFit.MarkerSize(1.5),RooFit.DataError(RooAbsData.Poisson), RooFit.XErrorSize(0),RooFit.Invisible());
        ## plot solid style
	if options.ttbarcontrolregion:
         model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
            
	 model_data.plotOn(mplot,RooFit.Name("WW_EWK"), RooFit.Components("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

	 model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj"%(self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model_STop_%s_mj%s,model_WJets0_%s_mj"%(self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model_WJets0_%s_mj"%(self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("TTbar_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("WW_EWK_invisible"), RooFit.Components("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
           
	 ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("VV_EWK_invisible"), RooFit.Components("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj"%(self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("STop_invisible"), RooFit.Components("model_STop_%s_mj,model_WJets0_%s_mj"%(self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
	 
	 ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("WJets_invisible"), RooFit.Components("model_WJets0_%s_mj"%(self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
	 
         ### solid line
         model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj"%(label+options.fgen,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj"%(self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj"%(self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("_invisible"), RooFit.Components("model%s_%s_mj,model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj,"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

        else:
         model_data.plotOn(mplot,RooFit.Name("WW_EWK"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
            
	 model_data.plotOn(mplot,RooFit.Name("VV"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label+options.fgen+options.fgen,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

	 model_data.plotOn(mplot,RooFit.Name("TTbar"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

	 model_data.plotOn(mplot,RooFit.Name("STop"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj"%(label+options.fgen,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

	 model_data.plotOn(mplot,RooFit.Name("WJets"), RooFit.Components("model%s_%s_mj"%(label+options.fgen,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("WW_EWK_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("VV_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
           
	 ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("TTbar_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("STop_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj"%(label+options.fgen,self.channel,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
	 
	 ## plot "dashed" style area
         model_data.plotOn(mplot,RooFit.Name("WJets_invisible"), RooFit.Components("model%s_%s_mj"%(label+options.fgen,self.channel)),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack),RooFit.FillStyle(3003),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());
	 
         ### solid line
         model_data.plotOn( mplot,RooFit.Name("WJets_line_invisible"), RooFit.Components("model%s_%s_mj"%(label+options.fgen,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("STop_line_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj"%(label+options.fgen,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("TTbar_line_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("VV_line_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_TTbar_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("WW_EWK_line_invisible"), RooFit.Components("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj,"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2),RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         ### dash line
         model_data.plotOn( mplot,RooFit.Name("WJets_dashed_invisible"), RooFit.Components("model%s_%s_mj"%(label+options.fgen,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("STop_dashed_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj"%(label+options.fgen,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("TTbar_dashed_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel)), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("VV_dashed_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"), RooFit.VLines());

         model_data.plotOn( mplot,RooFit.Name("WW_EWK_dashed_invisible"), RooFit.Components("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj"%(label+options.fgen,self.channel,self.channel,self.channel,self.channel,self.channel)),RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.LineStyle(kDashed) ,RooFit.NormRange("sb_lo,sb_hi"));
                                                                     
        ### draw the error band using the sum of all the entries component MC + fit and the total error == Normalization for the fixed MC, shape + normalization for W+jets
        draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,self.color_palet["Uncertainty"],"F");
        model_data.plotOn(mplot,RooFit.Name("model_mc"),RooFit.Range(rrv_mass_j.getMin(),rrv_mass_j.getMax()),RooFit.NormRange("sb_lo,sb_hi"),RooFit.Invisible());
                        
	if options.pseudodata == 1:
         rdataset_data_mj.plotOn( mplot ,RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Name("data"));               
        else: 
         GetDataPoissonInterval(rdataset_data_mj,rrv_mass_j,mplot);
		
        ### Get the pull and plot it
        mplot_pull = get_pull(rrv_mass_j,mplot,rdataset_data_mj,model_data,rfresult,"data","model_mc",1,1);
  
        ### signal window zone with vertical lines
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,mplot.GetMaximum()*0.9); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kGray+2); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,mplot.GetMaximum()*0.9); upperLine.SetLineWidth(2); upperLine.SetLineColor(kGray+2); upperLine.SetLineStyle(9);
        mplot.addObject(lowerLine);
        mplot.addObject(upperLine);

        ### legend of the plot
        leg = legend4Plot(mplot,0,-0.2,0.07,0.04,0.,1,"em");
        mplot.addObject(leg);
        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.5);

        ## CALCULATE CHI2
        datahist = rdataset_data_mj.binnedClone(rdataset_data_mj.GetName()+"_binnedClone",rdataset_data_mj.GetName()+"_binnedClone")
        Nbin = int(rrv_mass_j.getBins()); 
        rresult_param = rfresult.floatParsFinal();        
        nparameters =  rresult_param.getSize()                                         
        ChiSquare = model_data.createChi2(datahist,RooFit.Extended(kTRUE),RooFit.DataError(RooAbsData.Poisson));
        chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters);

        ## Add Chisquare to mplot_pull
        cs = TLatex(0.75,0.8,"#chi^{2}/ndf = %0.2f "%(float(chi_over_ndf)));
        cs.SetNDC();
        cs.SetTextSize(0.12);
        cs.AppendPad("same");
        mplot_pull.addObject(cs)

        parameters_list = model_data.getParameters(rdataset_data_mj);

        if not os.path.isdir("plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/")):
         os.system("mkdir -p plots_%s_%s_%s_g1/mj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/"));

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),"m_j_sideband%s"%(label),"","em",0,1,self.GetLumi());

                
        #### to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data
	if options.ttbarcontrolregion:
	 fullInt   = model_TTbar.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j) );
         signalInt = model_TTbar.createIntegral(RooArgSet(rrv_mass_j),RooArgSet(rrv_mass_j),("signal_region"));
         fullInt_val   = fullInt.getVal();
         signalInt_val = signalInt.getVal()/fullInt_val;
         ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
         rrv_number_TTbar_in_mj_signal_region_from_fitting = RooRealVar("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),"rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getVal()*signalInt_val);

         #### Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
         rrv_number_TTbar_in_mj_signal_region_from_fitting.setError( Calc_error_extendPdf(rdataset_data_mj, model_TTbar, rfresult,"signal_region") );
         print "########## error on the normaliztion due to shape + norm = %s"%(rrv_number_TTbar_in_mj_signal_region_from_fitting.getError());
         getattr(self.workspace4bias_,"import")(rrv_number_TTbar_in_mj_signal_region_from_fitting);
         rrv_number_TTbar_in_mj_signal_region_from_fitting.Print();

	else:            
         fullInt   = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooFit.NormSet(RooArgSet(rrv_mass_j)) );
         signalInt = model_WJets.createIntegral(RooArgSet(rrv_mass_j),RooFit.NormSet(RooArgSet(rrv_mass_j)),RooFit.Range("signal_region"));
         fullInt_val   = fullInt.getVal();
         signalInt_val = signalInt.getVal()/fullInt_val;
         ## take the value from the fit (normalization) and multiply it from the ratio of the integrals
         rrv_number_WJets_in_mj_signal_region_from_fitting = RooRealVar("rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),"rrv_number%s_in_mj_signal_region_from_fitting_%s"%(label,self.channel),self.workspace4bias_.var("rrv_number%s_%s_mj"%(label+options.fgen,self.channel)).getVal()*signalInt_val);

         #### Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
         rrv_number_WJets_in_mj_signal_region_from_fitting.setError( Calc_error_extendPdf(rdataset_data_mj, model_WJets, rfresult,"signal_region") );
         print "########## error on the normaliztion due to shape + norm = %s"%(rrv_number_WJets_in_mj_signal_region_from_fitting.getError());
         getattr(self.workspace4bias_,"import")(rrv_number_WJets_in_mj_signal_region_from_fitting);
         rrv_number_WJets_in_mj_signal_region_from_fitting.Print();
 
            
    ##### Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape
    def fit_mlvj_in_Mj_sideband(self, label, mlvj_region, mlvj_model,logy=0):

        print "############### Fit mlvj in mj sideband: ",label," ",mlvj_region," ",mlvj_model," ##################"

        rrv_mass_lvj = self.workspace4bias_.var("rrv_mass_lvj");
        rdataset_data_mlvj = self.workspace4bias_.data("rdataset_data%s_%s_mlvj"%(mlvj_region,self.channel))

        ## get and fix the minor component shapes in the sb low
        model_VV_backgrounds     = self.get_VV_mlvj_Model(mlvj_region,options.fgen);
        model_STop_backgrounds   = self.get_STop_mlvj_Model(mlvj_region,options.fgen);
        model_WW_EWK_backgrounds = self.get_WW_EWK_mlvj_Model(mlvj_region,options.fgen);
     
        number_VV_sb_lo_mlvj      = self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        number_STop_sb_lo_mlvj    = self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        number_WW_EWK_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization       

        if options.ttbarcontrolregion == 0:
         model_TTbar_backgrounds  = self.get_TTbar_mlvj_Model(mlvj_region,options.fgen);
         number_TTbar_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        else: 
         model_WJets_backgrounds  = self.get_WJets_mlvj_Model(mlvj_region,options.fgen);
         number_WJets_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
         

        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        print "rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel) ; 
        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();

        ### Make the Pdf for the WJets
        if options.ttbarcontrolregion == 0 :
         model_pdf_WJets = self.make_Pdf("%s%s%s_from_fitting"%(label,mlvj_region,mlvj_model), mlvj_model,"_mlvj");
         model_pdf_WJets.Print();
         ### inititalize the value to what was fitted with the mc in the sideband
         number_WJets_sb_lo = self.workspace4bias_.var("rrv_number%s%s%s_%s_mlvj"%(label,mlvj_region,options.fgen,self.channel)).clone("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+self.channel+"_mlvj");
    
         model_WJets        = RooExtendPdf("model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),model_pdf_WJets,number_WJets_sb_lo);
         number_WJets_sb_lo.Print();

         ## Add the other bkg component fixed to the total model --> in the extended way
         model_data = RooAddPdf("model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds, model_WW_EWK_backgrounds));

        else: 
         model_pdf_TTbar = self.make_Pdf("%s%s%s_from_fitting"%(label,mlvj_region,mlvj_model), mlvj_model,"_mlvj");
         model_pdf_TTbar.Print();
         ### inititalize the value to what was fitted with the mc in the sideband
         number_TTbar_signal_region = self.workspace4bias_.var("rrv_number%s%s%s_%s_mlvj"%(label,mlvj_region,options.fgen,self.channel)).clone("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+self.channel+"_mlvj");
    
         model_TTbar        = RooExtendPdf("model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),model_pdf_TTbar,number_TTbar_signal_region);
         number_TTbar_signal_region.Print();

         ## Add the other bkg component fixed to the total model --> in the extended way
         model_data = RooAddPdf("model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),RooArgList(model_TTbar,model_VV_backgrounds, model_WJets_backgrounds, model_STop_backgrounds, model_WW_EWK_backgrounds));
        

        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.SumW2Error(kTRUE));
        rfresult = model_data.fitTo( rdataset_data_mlvj, RooFit.Save(1) ,RooFit.Extended(kTRUE), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));
        rfresult.Print();
        rfresult.covarianceMatrix().Print();
        getattr(self.workspace4bias_,"import")(model_data)

        if options.ttbarcontrolregion == 0:
         model_WJets.getParameters(rdataset_data_mlvj).Print("v");

         ### data in the sideband plus error from fit        
         rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel),"rrv_number_data%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel),
                                                  self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_WJets0%s%s_from_fitting_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getVal() );

         rrv_number_data_sb_lo_mlvj.setError(TMath.Sqrt(self.workspace4bias_.var("rrv_number_WJets0%s%s_from_fitting_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_WJets0%s%s_from_fitting_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()));

         rrv_number_data_sb_lo_mlvj.Print();
         getattr(self.workspace4bias_,"import")(rrv_number_data_sb_lo_mlvj)

        else:        

         ### data in the sideband plus error from fit        
         rrv_number_data_sb_lo_mlvj = RooRealVar("rrv_number_data_%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel),"rrv_number_data%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel),
                                                  self.workspace4bias_.var("rrv_number_TTbar%s%s_from_fitting_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getVal()+
                                                  self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getVal() );

         rrv_number_data_sb_lo_mlvj.setError(TMath.Sqrt(self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_TTbar%s%s_from_fitting_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_TTbar%s%s_from_fitting_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()+
                                                        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()*
                                                        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).getError()));

         rrv_number_data_sb_lo_mlvj.Print();
         getattr(self.workspace4bias_,"import")(rrv_number_data_sb_lo_mlvj)

         model_TTbar.getParameters(rdataset_data_mlvj).Print("v");


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

        ### plot for WJets default + default shape
        if label=="_TTbar":

            mplot = rrv_mass_lvj.frame(RooFit.Title("M_lvj fitted in TTbar Control Region "), RooFit.Bins(int(rrv_mass_lvj.getBins())));

            rdataset_data_mlvj.plotOn( mplot , RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Invisible(), RooFit.Name("data_invisible") );

            model_data.plotOn(mplot, RooFit.Components(model_TTbar.GetName()+","+model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WW_EWK"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WW_EWK"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()), RooFit.Name("VV"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()), RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines());

            model_data.plotOn(mplot, RooFit.Components(model_STop_backgrounds.GetName()), RooFit.Name("STop"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines());


            #solid line
            model_data.plotOn(mplot,  RooFit.Components(model_TTbar.GetName()+","+model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()+","+model_WW_EWK_backgrounds.GetName()), RooFit.Name("WW_EWK_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()+","+model_VV_backgrounds.GetName()), RooFit.Name("VV_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_WJets_backgrounds.GetName()+","+model_STop_backgrounds.GetName()), RooFit.Name("WJets_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

            model_data.plotOn(mplot,  RooFit.Components(model_STop_backgrounds.GetName()), RooFit.Name("STop_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines());

        ### draw the error band
        draw_error_band(rdataset_data_mlvj, model_data,self.workspace4bias_.var("rrv_number_data_%s%s_%s_mlvj"%(mlvj_region,mlvj_model,self.channel)),rfresult,mplot,self.color_palet["Uncertainty"],"F");
        model_data.plotOn( mplot , RooFit.VLines(), RooFit.Invisible());
        model_data.plotOn( mplot , RooFit.Invisible(), RooFit.Name("model_mc"));

        if options.pseudodata == 1:
         rdataset_data_mlvj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0), RooFit.Name("data"));               
        else: 
         GetDataPoissonInterval(rdataset_data_mlvj,rrv_mass_lvj,mplot);

        mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
            
        ### Add the legend to the plot
        leg = legend4Plot(mplot,0,0.,0.06,0.16,0.,1,"em");
        mplot.addObject(leg);

        ### get the pull plot and store the canvas
        mplot_pull = get_pull(rrv_mass_lvj,mplot,rdataset_data_mlvj,model_data,rfresult,"data","model_mc",1,1);
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

        if not os.path.isdir("plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/")):
         os.system("mkdir -p plots_%s_%s_%s_g1/mj_fitting_%s_%s/"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres));
         os.system("mkdir plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel,self.wtagger_label,options.fgen,options.fres,"basePlot/"));

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),"m_lvj_sb_lo%s_%s"%(label,mlvj_model),"","em",0,1,self.GetLumi());

                
    ##### Method used to cycle on the events and for the dataset to be fitted
    def get_mj_and_mlvj_dataset(self,in_file_name, label, jet_mass="jet_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");

        rrv_mass_j         = self.workspace4bias_.var("rrv_mass_j") 
        rrv_mass_lvj       = self.workspace4bias_.var("rrv_mass_lvj")
        rrv_weight         = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 
         
	rdataset_mj     = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4bias_mj = RooDataSet("rdataset4bias"+label+"_"+self.channel+"_mj","rdataset4bias"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

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
          
          isFullVBF = 0 ;

          if options.ttbarcontrolregion == 0 or TString(in_file_name).Contains("ggH") or TString(in_file_name).Contains("vbfH"):

           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2:
            isFullVBF = 1 ;
          
           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") < 0.679 and getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV")<0.679 and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
            isFullVBF = 2 ;

          else:

           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and ( getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") > 0.679 or getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV") > 0.679 ) and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2:
            isFullVBF = 1 ;
          
           if ungroomed_jet_pt > 200. and discriminantCut and tmp_jet_mass >= rrv_mass_j.getMin() and tmp_jet_mass<=rrv_mass_j.getMax() and ( getattr(treeIn,"vbf_maxpt_j1_bDiscriminatorCSV") > 0.679 or getattr(treeIn,"vbf_maxpt_j2_bDiscriminatorCSV") > 0.679 ) and mass_lvj >= rrv_mass_lvj.getMin() and mass_lvj <=rrv_mass_lvj.getMax() and getattr(treeIn,"v_pt") > self.vpt_cut and pfMET > self.pfMET_cut and getattr(treeIn,"l_pt") > self.lpt_cut and getattr(treeIn,"issignal")==1 and ( getattr(treeIn,"mass_ungroomedjet_closerjet") > self.top_veto_had_max  or getattr(treeIn,"mass_ungroomedjet_closerjet") < self.top_veto_had_min ) and ( getattr(treeIn,"mass_leptonic_closerjet") > self.top_veto_lep_max or getattr(treeIn,"mass_leptonic_closerjet") < self.top_veto_lep_min) and njet >=2 and tmp_vbf_dEta > self.dEta_cut and tmp_vbf_Mjj > self.Mjj_cut:
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
                                     
	     if isFullVBF >= 2: 
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
    
    ##### Get Lumi for banner title
    def GetLumi(self):
        if self.channel=="el": return 19.3;
        if self.channel=="mu": return 19.3;
        if self.channel=="em": return 19.3;


    def shapeParametrizationAnalysis(self,label = "_WJets0"):

     if label == "_WJets0": fileName = self.file_WJets0_mc ;
     else:                  fileName = self.file_TTbar_mc ;
     mlvj_region = options.mlvjregion;

     if options.isMC == 1:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(fileName,label)# to get the shape of m_lvj                                                                                             
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfExp_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfExp_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfExp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfExp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Erf2Exp",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Erf2Exp",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfPow_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfPow_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfPow2_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfPow2_v1",1,0,1); 
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfPow3_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfPow3_v1",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Erf2Pow",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Erf2Pow",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfChebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfChebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfChebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfChebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ErfChebychev_v4",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ErfChebychev_v4",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Keys",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Keys",1,0,1);

      else:

       self.get_mj_and_mlvj_dataset(fileName,label)# to get the shape of m_lvj                                                                                             

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Exp",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Exp",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","ExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","ExpTail",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Exp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Exp_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","2Exp",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","2Exp",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Pow",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Pow",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Pow2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Pow2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Pow3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Pow3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","2Pow",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","2Pow",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Chebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Chebychev_v2",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Chebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Chebychev_v3",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Chebychev_v4",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Chebychev_v4",1,0,1);

       self.fit_mlvj_model_single_MC(fileName,label,"_sb_lo","Keys",1,0,1);
       self.fit_mlvj_model_single_MC(fileName,label,"_signal_region","Keys",1,0,1);
         
     else:

      if self.in_mlvj_min < 500:

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV",mlvj_region,"ErfExp_v1",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                               
       self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop",mlvj_region,"ErfExp_v1",0,0,1);
           
       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj
       self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",mlvj_region,"ErfExp_v1",0,0,1);
           
       self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj
       self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK",mlvj_region,"ErfExp_v1",0,0,1);
           
       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0")# to get the shape of m_lvj
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",mlvj_region,"ErfExp_v1",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_data,"_data"); ## global fit of data in the sidand fixing non dominant bkg

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfExp_v1");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfExpTail");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfExp_v3");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Erf2Exp");

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfPow_v1");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfPow2_v1");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfPow3_v1");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Erf2Pow");

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfChebychev_v2");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfChebychev_v3");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfChebychev_v4");

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ErfPowExp_v1");

      else: 

       self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV", "jet_mass_pr");
       self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV",mlvj_region,"Exp",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop");
       self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop",mlvj_region,"Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar");
       self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",mlvj_region,"Exp",0,0,1);
 
       self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr");
       self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK",mlvj_region,"Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0");
       self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",mlvj_region,"Exp",0,0,1);

       self.get_mj_and_mlvj_dataset(self.file_data,"_data");
       
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Exp");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"ExpTail");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Exp_v3");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"2Exp");

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Pow");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Pow2");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Pow3");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"2Pow");

       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Chebychev_v2");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Chebychev_v3");
       self.fit_mlvj_in_Mj_sideband(label,mlvj_region,"Chebychev_v4");

    def biasAnalysis(self,label="_WJets0",fitjetmass = 0):
   
     print"######################## begin the bias analysis ###########################";  
     if fitjetmass == 1:
         options.mlvjregion = "";
     if fitjetmass == 1 and options.isMC == 1:
         print " cannot use MC to perform bias test for the mJ fit -> options not provided" ;
         return ;

     ###### get the signal and fit it     
     self.get_mj_and_mlvj_dataset(self.file_ggH ,"_%s"%(self.ggH_sample),"jet_mass_pr")# to get the shape of m_lvj
     self.get_mj_and_mlvj_dataset(self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"jet_mass_pr")# to get the shape of m_lvj

     if fitjetmass:
	self.fit_mj_single_MC(self.file_ggH,"_%s"%(self.ggH_sample),"2Gaus",1);
	self.fit_mj_single_MC(self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"2Gaus",1);
     else:
        self.fit_mlvj_model_single_MC(self.file_ggH,"_%s"%(self.ggH_sample),"_signal_region","CB_v1", 0, 0, 1);
        self.fit_mlvj_model_single_MC(self.file_vbfH,"_%s"%(self.vbfhiggs_sample),"_signal_region","CB_v1", 0, 0, 1);

     ####### Monte Carlo Analysis
     if options.isMC and options.ttbarcontrolregion:
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fgen,0,0,1);

     elif options.isMC and not options.ttbarcontrolregion:    
      self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
      self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fgen,0,0,1);
   
     elif not options.isMC:

      ###### get diboson and fit it 
      self.get_mj_and_mlvj_dataset(self.file_VV_mc,"_VV","jet_mass_pr");      
      if fitjetmass:
          self.fit_mj_single_MC(self.file_VV_mc,"_VV","2_2Gaus");
      else:
          self.fit_mlvj_model_single_MC(self.file_VV_mc,"_VV",options.mlvjregion,options.fgen,0,0,1); 

      ####### get SingleTop and fit it
      self.get_mj_and_mlvj_dataset(self.file_STop_mc,"_STop")# to get the shape of m_lvj                                                                                                 
      if fitjetmass: 
	  if options.ttbarcontrolregion:
	     self.fit_mj_single_MC(self.file_STop_mc,"_STop","ErfExpGaus_sp");
	  else:
             self.fit_mj_single_MC(self.file_STop_mc,"_STop","2Gaus_ErfExp");		   
      else: self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop",options.mlvjregion,options.fgen,0,0,1);                  

      ######## get WW EWK and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                                             
      if fitjetmass:
          self.fit_mj_single_MC(self.file_WW_EWK_mc,"_WW_EWK","2Gaus"); 
      else:
          self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK",options.mlvjregion,options.fgen,0,0,1);
      ###### get WJets and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
      if fitjetmass: 
	   self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0","ErfExp",1);
      else: 
           self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fgen,0,0,1);
	
      ######## get TTbar and fit it
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      if fitjetmass: 
       if options.ttbarcontrolregion :    
          self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar",options.fgen,1);
       else:
          self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar","2Gaus_ErfExp");	   		   
      else: self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fgen,0,0,1);

      ##### get data in sb and fit it                                                                                                                                              
      self.get_mj_and_mlvj_dataset(self.file_data,"_data", "jet_mass_pr"); ## global fit of data in the sidand fixing non dominant bkg             
      if fitjetmass :
         self.fit_WJetsNorm(label); ## fit jet mass distribution
      else:      
         self.fit_mlvj_in_Mj_sideband(label,options.mlvjregion,options.fgen,1); ## sideband or TTbar signal region fit

     ##### fix signal and bkg models that are going to be used in the generation
     if fitjetmass :
         spectrum = "_mj"   ;
         signal_region = "";
         signal_model  = "2Gaus";
     else:
         spectrum = "_mlvj" ;
         signal_region = "_signal_region";
         signal_model  = "CB_v1";
    
     self.fix_Model("_%s"%self.ggH_sample,signal_region,spectrum,signal_model);
     self.fix_Model("_%s"%self.vbfhiggs_sample,signal_region,spectrum,signal_model);

     ###### fix the backgrund models for the generation     
     if options.isMC == 0  and options.ttbarcontrolregion == 0 and not options.fitjetmass :

      self.fix_Model("_TTbar" ,options.mlvjregion,spectrum,options.fgen);
      self.fix_Model("_STop"  ,options.mlvjregion,spectrum,options.fgen);
      self.fix_Model("_VV"    ,options.mlvjregion,spectrum,options.fgen) ;
      self.fix_Model("_WW_EWK",options.mlvjregion,spectrum,options.fgen) ;
      self.fix_Model(label,options.mlvjregion,spectrum,options.fgen);

     elif options.isMC == 0  and options.ttbarcontrolregion == 1 and not options.fitjetmass:

      self.fix_Model("_WJets0",options.mlvjregion,spectrum,options.fgen);
      self.fix_Model("_STop"  ,options.mlvjregion,spectrum,options.fgen);
      self.fix_Model("_VV"    ,options.mlvjregion,spectrum,options.fgen) ;
      self.fix_Model("_WW_EWK",options.mlvjregion,spectrum,options.fgen) ;
      self.fix_Model(label,options.mlvjregion,spectrum,options.fgen);

     elif options.fitjetmass and not options.ttbarcontrolregion:

      self.fix_Model("_STop"  ,options.mlvjregion,spectrum);
      self.fix_Model("_VV"    ,options.mlvjregion,spectrum);
      self.fix_Model("_WW_EWK",options.mlvjregion,spectrum);
      self.fix_Model("_TTbar" ,options.mlvjregion,spectrum);
      self.fix_Model(label,options.mlvjregion,spectrum,options.fgen);

     elif options.fitjetmass and options.ttbarcontrolregion:

      self.fix_Model("_STop"  ,options.mlvjregion,spectrum);
      self.fix_Model("_VV"    ,options.mlvjregion,spectrum);
      self.fix_Model("_WW_EWK",options.mlvjregion,spectrum);
      self.fix_Model("_WJets0" ,options.mlvjregion,spectrum);
      self.fix_Model(label,options.mlvjregion,spectrum,options.fgen);

     
     ### clone the signal shape --> parameter already fixed
     fitted_signal_ggH   = self.workspace4bias_.pdf("model_%s_%s_%s%s"%(self.ggH_sample,signal_region,self.channel,spectrum));
     fitted_signal_vbfH  = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.vbfhiggs_sample,signal_region,self.channel,spectrum));

     ### make the pdf for ggH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_ggH = [];

     if fitjetmass:
      model_signal_ggH = self.make_Pdf("_%s%s_fit"%(self.ggH_sample,signal_region+"2Gaus"),"2Gaus",spectrum,constrainslist_signal_ggH,1);
      model_signal_ggH.Print();     
      self.clone_Model(model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,"2Gaus");
      self.fix_Model("_%s"%self.ggH_sample,signal_region,spectrum,"2Gaus","_fit",1);
      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"2Gaus_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"2Gaus_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);

     else:
      model_signal_ggH = self.make_Pdf("_%s%s_fit"%(self.ggH_sample,signal_region+"CB_v1"),"CB_v1",spectrum,constrainslist_signal_ggH,1);
      model_signal_ggH.Print();     
      self.clone_Model(model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,"CB_v1");
      self.fix_Model("_%s"%self.ggH_sample,signal_region,spectrum,"CB_v1","_fit",1);
      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"CB_v1_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"CB_v1_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);
     
      
     ### make the pdf for vbfH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_vbfH = [];
     if fitjetmass:
      model_signal_vbfH = self.make_Pdf("_%s%s_fit"%(self.vbfhiggs_sample,signal_region+"2Gaus"),"2Gaus",spectrum,constrainslist_signal_vbfH,1);
      model_signal_vbfH.Print();        
      self.clone_Model(model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,"2Gaus");
      self.fix_Model("_%s"%self.vbfhiggs_sample,signal_region,spectrum,"2Gaus","_fit",1);      
      rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"2Gaus_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"2Gaus_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_vbfH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_vbfH);

     else: 
      model_signal_vbfH = self.make_Pdf("_%s%s_fit"%(self.vbfhiggs_sample,signal_region+"CB_v1"),"CB_v1",spectrum,constrainslist_signal_vbfH,1);
      model_signal_vbfH.Print();        
      self.clone_Model(model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,"CB_v1");
      self.fix_Model("_%s"%self.vbfhiggs_sample,signal_region,spectrum,"CB_v1","_fit",1);      
      rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"CB_v1_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"CB_v1_"+self.channel+spectrum).getError());
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

      constrainslist_bkg_wjet = [];
      model_bkg_wjet    = self.make_Model(label+options.mlvjregion+"_fit",options.fres,spectrum,constrainslist_bkg_wjet,1); ## only mWW analysis can be taken into account
      model_bkg_wjet.Print();
      if options.fres == options.fgen :
       self.clone_Model(model_bkg_wjet,label,options.mlvjregion,spectrum,options.fgen); ## clone the parameter from the old mc fit
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
      numevents_mc   = self.workspace4bias_.data("rdataset"+label+options.mlvjregion+"_"+self.channel+spectrum).sumEntries()*options.inflatejobstatistic;
      print"########  numevents mc ",numevents_mc;

      mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                                           generation_model_wjet,
                                           self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                           int(options.nexp),
                                           int(options.isMC));

      mcWjetTreeResult.setTree(self.outputTree);
      mcWjetTreeResult.setFittingModel(model_Total_mc);
      mcWjetTreeResult.setPdfInfomation(options.mlvjregion,spectrum,self.channel);
      mcWjetTreeResult.setBackgroundPdfCore(model_bkg_wjet);
      mcWjetTreeResult.generateAndFitToys(int(numevents_mc));
      self.outputFile.cd();
      mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
      mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_);

      ratePlotsToStore = 0 ;
      if options.nexp <= 10 :
          ratePlotsToStore = 1 ;
      elif options.nexp > 10 and options.nexp < 50:
          ratePlotsToStore = 2 ;
      elif options.nexp >= 50 and options.nexp < 100:
          ratePlotsToStore = 3 ;
      elif options.nexp >= 100:
          ratePlotsToStore = 10 ;
          
      mcWjetTreeResult.saveToysPlots(int(ratePlotsToStore),options.fitjetmass); 
      self.outputTree.Write();
      self.outputFile.Close();

     else: 
         
      ############### Make the Data analysis --> make the Entended pdf for the bkg
      constrainslist_bkg_data = [];

      ### take the models for the background component
      if fitjetmass:

       model_VV_backgrounds     = self.get_VV_mj_Model("_VV");
       model_STop_backgrounds   = self.get_STop_mj_Model("_STop");
       model_WW_EWK_backgrounds = self.get_WW_EWK_mj_Model("_WW_EWK");
       ## inflate the number of events and print them
       self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);

       print " VV number ",self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print " STop number ",self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print " WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.ttbarcontrolregion:
        model_TTbar_backgrounds  = self.get_TTbar_mj_Model(label,options.fgen);
        model_WJets_backgrounds  = self.get_WJets_mj_Model("_WJets0");

        self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);

        print " WJets number ",self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       else:

        model_TTbar_backgrounds  = self.get_TTbar_mj_Model("_TTbar");           
        model_WJets_backgrounds  = self.get_WJets_mj_Model(label,options.fgen);

        self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic)

        print " WJets number ",self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

      else:	       
       ### in case of mWW analysis
       model_VV_backgrounds     = self.get_VV_mlvj_Model(options.mlvjregion,options.fgen);
       model_STop_backgrounds   = self.get_STop_mlvj_Model(options.mlvjregion,options.fgen);
       model_TTbar_backgrounds  = self.get_TTbar_mlvj_Model(options.mlvjregion,options.fgen);
       model_WW_EWK_backgrounds = self.get_WW_EWK_mlvj_Model(options.mlvjregion,options.fgen);
       model_WJets_backgrounds  = self.get_WJets_mlvj_Model(options.mlvjregion,options.fgen);

       ## inflate yields
       self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic) ## get the normalization

       print " VV number ",self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print " STop number ",self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print " WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.ttbarcontrolregion == 0:

        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic)  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print " TTbar number ",self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " WJets number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
       else:

        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print " WJets number ",self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " TTbar number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
           
      #### make the global model for the background  
      if options.fitjetmass:
       model_bkg_data    = self.make_Model("_data"+signal_region+"_fit",options.fres,spectrum,constrainslist_bkg_data,1); ## basic model used for fit in the toys
       model_bkg_data.Print();

       if options.fgen == options.fres:
        self.clone_Model(model_bkg_data,label,signal_region,spectrum,options.fgen);

       self.workspace4bias_.var("rrv_number_data"+signal_region+"_fit_"+self.channel+spectrum).setVal(self.workspace4bias_.var("rrv_number"+label+signal_region+options.fgen+"_"+self.channel+spectrum).getVal());
       self.workspace4bias_.var("rrv_number_data"+signal_region+"_fit_"+self.channel+spectrum).Print();

      else:
       model_bkg_data    = self.make_Model("_data"+options.mlvjregion+"_fit",options.fres,spectrum,constrainslist_bkg_data,1); ## basic model used for fit in the toys
       model_bkg_data.Print();

       if options.fgen == options.fres :
        self.clone_Model(model_bkg_data,label,options.mlvjregion,spectrum,options.fgen+"_from_fitting");        

       self.workspace4bias_.var("rrv_number_data"+options.mlvjregion+"_fit_"+self.channel+spectrum).setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+spectrum).getVal());
       self.workspace4bias_.var("rrv_number_data"+options.mlvjregion+"_fit_"+self.channel+spectrum).Print();

      ## Add the other bkg component fixed to the total model --> in the extended way
      if options.onlybackgroundfit == 1 and options.ttbarcontrolregion == 0:
        model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
      elif options.onlybackgroundfit == 1 and options.ttbarcontrolregion == 1: 
        model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
      elif options.onlybackgroundfit == 0 and options.ttbarcontrolregion == 0:
        model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_TTbar_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
      elif options.onlybackgroundfit == 0 and options.ttbarcontrolregion == 1:  
        model_Total_data = RooAddPdf("model_Total_background_data","model_Total_background_data",RooArgList(model_total_signal,model_bkg_data,model_VV_backgrounds,model_WJets_backgrounds, model_STop_backgrounds,model_WW_EWK_backgrounds));
                                                                                                       
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
      
      numevents_data   = self.workspace4bias_.data("rdataset"+"_data"+signal_region+"_"+self.channel+spectrum).sumEntries()*options.inflatejobstatistic;

      if options.fitjetmass :
       mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_j")),
                                            generation_model_data,
                                            self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                            int(options.nexp),
                                            int(options.isMC));

       mcWjetTreeResult.setTree(self.outputTree);
       mcWjetTreeResult.setFittingModel(model_Total_data);
       mcWjetTreeResult.setBackgroundPdfCore(model_bkg_data);
       mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel);
       mcWjetTreeResult.generateAndFitToys(int(numevents_data),"sb_lo,sb_hi");
       self.outputFile.cd();
       mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
       mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_);
      else:
       mcWjetTreeResult = biasModelAnalysis(RooArgSet(self.workspace4bias_.var("rrv_mass_lvj")),
                                            generation_model_data,
                                            self.workspace4bias_.data("rdataset4bias%s%s_%s%s"%(label,options.mlvjregion,self.channel,spectrum)),
                                            int(options.nexp),
                                            int(options.isMC));
       mcWjetTreeResult.setTree(self.outputTree);
       mcWjetTreeResult.setFittingModel(model_Total_data);
       mcWjetTreeResult.setBackgroundPdfCore(model_bkg_data);
       mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel);
       mcWjetTreeResult.generateAndFitToys(int(numevents_data));
       self.outputFile.cd();
       mcWjetTreeResult.createBranches(options.fgen,options.fres,options.ttbarcontrolregion);
       mcWjetTreeResult.fillBranches(options.ttbarcontrolregion,options.fitjetmass,self.workspace4bias_);


      ratePlotsToStore = 0 ;
      if options.nexp <= 10 :
          ratePlotsToStore = 1 ;
      elif options.nexp > 10 and options.nexp < 50:
          ratePlotsToStore = 2 ;
      elif options.nexp >= 50 and options.nexp < 100:
          ratePlotsToStore = 3 ;
      elif options.nexp >= 100:
          ratePlotsToStore = 10 ;
          
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
     fitBiasAnalysis.shapeParametrizationAnalysis("_WJets0",options.fitjetmass);
  elif options.shapetest == 0 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.biasAnalysis("_TTbar",options.fitjetmass);
     fitBiasAnalysis.outputFile.Close();                                  
  elif options.shapetest == 1 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.shapeParametrizationAnalysis("_TTbar",options.fitjetmass);                     
