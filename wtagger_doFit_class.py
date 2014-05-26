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
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="el")


parser.add_option('--fitwtagger', action='store_true', dest='fitwtagger', default=True, help='fit wtagger jet in ttbar control sample')
parser.add_option('--fitwtaggersim', action='store_true', dest='fitwtaggersim', default=False, help='fit wtagger jet in ttbar control sample with mu and el samples simultaneously')

parser.add_option('--inPath', action="store",type="string",dest="inPath",default="./")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")

parser.add_option('--herwig', action="store",type="int",dest="herwig",default=0)
parser.add_option('--ttbarMC', action="store",type="int",dest="ttbarMC",default=0)
parser.add_option('--pTbin', action="callback",callback=foo_callback,type="string",dest="pTbin",default="")
parser.add_option('--shift', action="store", type="int",dest="shift",default=0)
parser.add_option('--smear', action="store", type="int",dest="smear",default=0)

(options, args) = parser.parse_args()


ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
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

    ## contructor taking channel, signal name, range in mj, label and a workspace
    def __init__(self, in_channel,in_signal_sample, in_mj_min=30, in_mj_max=140, label="", input_workspace=None):

        print " "
        print "#####################################################"
        print "## Constructor of the fit object doFit_wj_and_wlvj ##"
        print "#####################################################"
        print " "

        #set plots style
        self.setTDRStyle();

        RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9) ;
        RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9) ;

        ### set the channel type --> electron or muon
        self.channel = in_channel;

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


        ### Set the mj binning for plots
        self.BinWidth_mj = 5.;
        
        #narrow the BinWidth_mj and BinWidth_mlvj by a factor of 5. Because Higgs-Combination-Tools will generate a binned sample, so need the bin width narrow. So, as a easy selution, we will increase the bin-width by a factor of 5 when ploting m_j m_WW
        self.narrow_factor=1.;

        ## set the range max properly in order to have a integer number of bins
        self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor;
        nbins_mj = int((in_mj_max-in_mj_min)/self.BinWidth_mj);
        in_mj_max = in_mj_min+nbins_mj*self.BinWidth_mj;

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
        self.mj_signal_min = 65;
        self.mj_signal_max = 105;
        self.mj_sideband_hi_min = 105;
        self.mj_sideband_hi_max = in_mj_max;

        ## define ranges on mj
        rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max);
        rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max);
        rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max);
        rrv_mass_j.setRange("controlsample_fitting_range",40,130);

        ## directory where are the trees to be run
        self.file_Directory="trainingtrees_exo_%s/"%(self.channel);

        ## taking the root file for data and mc
        self.file_data=("ofile_data.root");

        self.signal_sample = in_signal_sample;
        self.file_signal=("ofile_%s.root"%(self.signal_sample));

        self.file_pseudodata=("ofile_pseudodata4exo.root");
        self.file_pseudodata_herwig=("ofile_pseudodata4exo_herwig.root");

        if self.channel!="merged": self.file_WJets0_mc=("ofile_WJets_Pythia180.root");
        else: self.file_WJets0_mc=("ofile_WJets_Pythia100.root");

        self.file_WJets1_mc=("ofile_WJets_Herwig.root");

        self.file_VV_mc=("ofile_VV.root");# WW+WZ
        self.file_TTbar_mc=("ofile_TTbar_Powheg.root"); ## powheg TTbar
        self.file_TTbar_herwig=("ofile_TTbar_mcanlo.root"); ## mc@nlo TTbar
        self.file_TTbar_matchDn_mc=("ofile_TTbar_matchDn.root"); ## madgraph matching down
        self.file_TTbar_matchUp_mc=("ofile_TTbar_matchUp.root"); ## madgraph matching up
        self.file_TTbar_scaleDn_mc=("ofile_TTbar_scaleDn.root"); ## madgraph qcd scale down
        self.file_TTbar_scaleUp_mc=("ofile_TTbar_scaleUp.root"); ## madgraph qcd scale up
        self.file_TTbar_MG_mc=("ofile_TTbar_MG.root"); ## madgraph ttbar
        self.file_STop_mc =("ofile_STop.root"); ##single Top
 
        ## Define the workin poit on the N-subjettines cut

        self.wtagger_label=options.category;

        if self.wtagger_label=="HP" :
            if self.channel=="el"    :self.wtagger_cut=0.75 ; self.wtagger_cut_min=0. ;
            if self.channel=="mu"    :self.wtagger_cut=0.75 ; self.wtagger_cut_min=0. ;
            if self.channel=="merged":self.wtagger_cut=0.75 ; self.wtagger_cut_min=0. ;

        if self.wtagger_label=="LP": self.wtagger_cut=0.75 ; self.wtagger_cut_min=0.5 ;

        if self.wtagger_label=="nocut": self.wtagger_cut=10000;

        self.categoryID=-1;
        if self.wtagger_label=="LP" and self.channel=="el": self.categoryID=0;
        if self.wtagger_label=="HP" and self.channel=="el": self.categoryID=1;

        if self.wtagger_label=="LP" and self.channel=="mu": self.categoryID=2;
        if self.wtagger_label=="HP" and self.channel=="mu": self.categoryID=3;

        if self.wtagger_label=="LP" and self.channel=="merged": self.categoryID=4;
        if self.wtagger_label=="HP" and self.channel=="merged": self.categoryID=5;

        self.color_palet={ #color palet
            'data' : 1,
            'WJets' : 2,
            'VV' : 4,
            'STop' : 7,
            'TTbar' : 210,
            'TTbar_herwig' :210,
            'ggH' : 1,
            'vbfH' : 12,
            'Signal': 1,
            'Uncertainty' : kBlack,
            'Other_Backgrounds' : kBlue
        }

        ## Some basic cut vaule

        self.vpt_cut = 200; ## hadronic and leptonic W cut
        self.mass_lvj_max = 2000.; ## invariant mass of 3 body max
        self.mass_lvj_min = 200.; ## invariant mass of 3 body min
        self.pfMET_cut = 50; ## missing transverse energy
        self.lpt_cut = 50; ## lepton pT
        self.deltaPhi_METj_cut = 2.0;

        ## binning in the W jet pT for differential SF study in bins of pT

        if options.pTbin != "" :
         self.ca8_ungroomed_pt_min = int(options.pTbin[0]);
         self.ca8_ungroomed_pt_max = int(options.pTbin[1]);
        else:
         self.ca8_ungroomed_pt_min = 200;
         self.ca8_ungroomed_pt_max = 10000;

        
        ## tighter cut for the electron channel
        if self.channel=="el" or self.channel=="merged":
            self.pfMET_cut= 80; self.lpt_cut = 90;

        ## out txt file with info about fit and event couting
        self.file_ttbar_control_txt = "ttbar_control_%s_%s_wtaggercut%s.txt"%(self.signal_sample,self.channel,self.wtagger_label);
        self.file_out_ttbar_control=open(self.file_ttbar_control_txt,"w");

        self.mean_shift = 1.4 ; self.sigma_scale = 1.08 ;

        ### set the TDR Style                                                                                                                                                             
        setTDRStyle();


    ### Make a generic model for the ttbar control sample fit
    def make_Model_for_ttbar_controlsample(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[],info=""):

        print " "
        print "###################################################################################"
        print " ## make_Model_for_ttbar_controlsample : ",label," ",in_model_name," ",info," "##";
        print "###################################################################################"
        print " "
        
        if TString(label).Contains("_ttbar_data") and not TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = RooRealVar("rrv_number_total_ttbar_data"+info+"_"+self.channel,"rrv_number_total_ttbar_data"+info+"_"+self.channel,500,0.,1e7);
            eff_ttbar = RooRealVar("eff_ttbar_data"+info+"_"+self.channel,"eff_ttbar_data"+info+"_"+self.channel,0.7,0.3,0.9);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "@0*@1", RooArgList(rrv_number_total,eff_ttbar));

        ## fail fit has to be done after the pass one
        elif TString(label).Contains("_ttbar_data") and TString(label).Contains("failtau2tau1cut"):
            ## correlated paramters among the two cateogry of pass and fail sample
            rrv_number_total = self.workspace4fit_.var("rrv_number_total_ttbar_data"+info+"_"+self.channel);
            eff_ttbar = self.workspace4fit_.var("eff_ttbar_data"+info+"_"+self.channel);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total));

        elif TString(label).Contains("_ttbar_TotalMC") and not TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = RooRealVar("rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel,"rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel,500,0.,1e7);
            eff_ttbar = RooRealVar("eff_ttbar_TotalMC"+info+"_"+self.channel,"eff_ttbar_TotalMC"+info+"_"+self.channel,0.7,0.3,0.9);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total));

        elif TString(label).Contains("_ttbar_TotalMC") and TString(label).Contains("failtau2tau1cut"):
            rrv_number_total = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+info+"_"+self.channel);
            eff_ttbar = self.workspace4fit_.var("eff_ttbar_TotalMC"+info+"_"+self.channel);
            rrv_number = RooFormulaVar("rrv_number"+label+"_"+self.channel+mass_spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total)) ;

        rrv_number_total.Print();
        eff_ttbar.Print();
        rrv_number.Print();
        ## make the pdf for the given model
        model_pdf = self.make_Pdf(label,in_model_name,mass_spectrum,ConstraintsList)
        print " model Pdf :"
        model_pdf.Print();
        ## make the extended pdf via rrv_number
        print " extended PDF :"
        model = RooExtendPdf("model"+label+"_"+self.channel+mass_spectrum,"model"+label+"_"+self.channel+mass_spectrum, model_pdf, rrv_number );
        ## import in the worspace
        getattr(self.workspace4fit_,"import")(model)
        print "model"+label+"_"+self.channel+mass_spectrum
        self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum).Print();
        ## return the model
        return self.workspace4fit_.pdf("model"+label+"_"+self.channel+mass_spectrum);



    ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

      ## import variable and dataset
      rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
      rdataset_mj = self.workspace4fit_.data("rdataset4fit"+label+"_"+self.channel+"_mj");
      rdataset_mj.Print();

      ## make the extended model
      model = self.make_Model(label,in_model_name);
      ## Fit Sequence setting the errors correctly
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));
      rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2"));

      ## Plot the result
      mplot = rrv_mass_j.frame(RooFit.Title(label+" fitted by "+in_model_name), RooFit.Bins(int(rrv_mass_j.getBins())));
      rdataset_mj.plotOn(mplot, RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

      ## draw the error band for an extend pdf
      draw_error_band_extendPdf(rdataset_mj, model, rfresult,mplot,6,"L");
      ## re-draw the dataset
      rdataset_mj.plotOn( mplot , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
      ## draw the function
      model.plotOn(mplot);

      ## Get the pull
      mplot_pull = self.get_pull(rrv_mass_j, mplot);
      mplot.GetYaxis().SetRangeUser(1e-2,mplot.GetMaximum()*1.2);
      
      parameters_list = model.getParameters(rdataset_mj);
      self.draw_canvas_with_pull( mplot,mplot_pull,parameters_list,"plots_%s_%s_%s/m_j_fitting%s_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation,self.channel,self.wtagger_label,additioninformation, self.wtagger_label),label+in_file_name+"_"+str(self.ca8_ungroomed_pt_min)+"_"+str(self.ca8_ungroomed_pt_max), in_model_name)
      
      #normalize the number of total events to lumi
      if options.ttbarMC == 0 :

       self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal()  )
       self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4fit_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal()  )
       if TString(label).Contains("ggH"):
            self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal( self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()  )
            self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").getError())
       self.workspace4fit_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

       #correct the mean and sigma from the ttbar contral sample
       par=parameters_list.createIterator();
       par.Reset();
       param=par.Next()
       while (param):
           if (TString(label).Contains("VV") or TString(label).Contains("STop") or TString(label).Contains("TTbar")) and (not (options.fitwtaggersim or options.fitwtagger)):
                #param.Print();
                if TString(param.GetName()).Contains("rrv_mean1_gaus"):
                    param.setRange(param.getMin()+self.mean_shift, param.getMax()+self.mean_shift);
                    param.setVal(param.getVal()+self.mean_shift);
                    #param.Print(); raw_input("mean"+label);
                if TString(param.GetName()).Contains("rrv_deltamean_gaus"):
                    param.setRange(param.getMin()-self.mean_shift, param.getMax()-self.mean_shift);
                    param.setVal(param.getVal()-self.mean_shift);
                    #param.Print(); raw_input("mean"+label);
                if TString(param.GetName()).Contains("rrv_sigma1_gaus"):
                    param.setVal(param.getVal()*self.sigma_scale);
                    param.setRange(param.getMin()*self.sigma_scale, param.getMax()*self.sigma_scale);
                    #param.Print(); raw_input("sigma"+label);
                if TString(param.GetName()).Contains("rrv_scalesigma_gaus"):
                    param.setRange(param.getMin()/self.sigma_scale, param.getMax()/self.sigma_scale);
                    param.setVal(param.getVal()/self.sigma_scale);
                    #param.Print(); raw_input("sigma"+label);
           param=par.Next()
                                       

    ############# -----------
    def fit_mj_TTbar_controlsample(self,in_file_name,label=""):

        ##### Print the final result for the number of events passing the cut and before the cut + the efficiency for the W-tagging

        self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_VV"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_WJets0"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_STop"+label+"_"+self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar"+label+"_"+self.channel+"_mj").Print()

        number_dataset_signal_region_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region: %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj);

        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel));
        self.file_out_ttbar_control.write("event number of data in signal_region: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj));

        number_dataset_signal_region_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();

        print "event number of TotalMC %s in signal_region: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj));


        number_dataset_signal_region_before_cut_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_data"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_data"+label+"_"+self.channel+"_mj").getVal();
        print "event number of data in signal_region before_cut: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj);
        self.file_out_ttbar_control.write("event number of data in signal_region before_cut: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj));

        number_dataset_signal_region_before_cut_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        number_dataset_signal_region_before_cut_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_TotalMC"+label+"_"+self.channel+"_mj").getVal();
        print "event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) "%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj);
        self.file_out_ttbar_control.write("event number of TotalMC %s in signal_region before_cut: %s +/- sqrt(%s) \n"%(label,number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj));

                                                             
        # wtagger_eff reweight: only reweight the efficiency difference between MC and data
        wtagger_eff_MC = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj;
        wtagger_eff_data= number_dataset_signal_region_data_mj/number_dataset_signal_region_before_cut_data_mj;

        wtagger_eff_reweight=wtagger_eff_data/wtagger_eff_MC;
        wtagger_eff_reweight_err=wtagger_eff_reweight*TMath.Sqrt(number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +number_dataset_signal_region_before_cut_error2_data_mj/number_dataset_signal_region_before_cut_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_before_cut_error2_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj);
        
        print "wtagger efficiency of %s channel"%(self.channel )
        print "wtagger_eff_MC %s = %s "%(label,wtagger_eff_MC )
        print "wtagger_eff_data = %s "%(wtagger_eff_data )
        print "wtagger_eff_reweight %s = %s +/- %s"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err)

        self.file_out_ttbar_control.write("wtagger_eff_MC %s = %s \n"%(label,wtagger_eff_MC ));
        self.file_out_ttbar_control.write("wtagger_eff_data = %s \n"%(wtagger_eff_data ));
        self.file_out_ttbar_control.write("wtagger_eff_reweight %s= %s +/- %s\n"%(label,wtagger_eff_reweight, wtagger_eff_reweight_err));


    ## method to do fail and pass fit for the scale factor evaluation
    def ScaleFactor_forPureWJet_TTbar_controlsample(self,in_file_name,label="",fit_or_not=1):

        print " ############################## Pass: dataset #################################### "
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_data_mj = self.workspace4fit_.data("rdataset_data"+label+"_"+self.channel+"_mj");
        rdataset_TTbar_mj = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+self.channel+"_mj");
        rdataset_STop_mj = self.workspace4fit_.data("rdataset_STop"+label+"_"+self.channel+"_mj");
        rdataset_VV_mj = self.workspace4fit_.data("rdataset_VV"+label+"_"+self.channel+"_mj");
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj)

        # inmport the pdf for each background
        print " ############################### Pass: model ################################### "
        model_histpdf_TTbar = self.workspace4fit_.pdf(rdataset_TTbar_mj.GetName()+"_histpdf")
        model_histpdf_STop = self.workspace4fit_.pdf(rdataset_STop_mj.GetName()+"_histpdf")
        model_histpdf_VV = self.workspace4fit_.pdf(rdataset_VV_mj.GetName()+"_histpdf")
        model_histpdf_WJets = self.workspace4fit_.pdf(rdataset_WJets_mj.GetName()+"_histpdf")
         
        number_TTbar = RooRealVar("rrv_number_TTbar"+label+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+self.channel,rdataset_TTbar_mj.sumEntries());
        number_STop = RooRealVar("rrv_number_STop"+label+"_"+self.channel ,"rrv_number_STop"+label+"_"+self.channel,rdataset_STop_mj.sumEntries());
        number_VV = RooRealVar("rrv_number_VV"+label+"_"+self.channel ,"rrv_number_VV"+label+"_"+self.channel,rdataset_VV_mj.sumEntries());
        number_WJets = RooRealVar("rrv_number_WJets"+label+"_"+self.channel,"rrv_number_WJets"+label+"_"+self.channel,rdataset_WJets_mj.sumEntries());

        #### Total MC Model
        model_TTbar_STop_VV_WJets = RooAddPdf("model_TTbar_STop_VV_WJets"+label+"_"+self.channel,"model_TTbar_STop_VV_WJets"+label+"_"+self.channel,RooArgList(model_histpdf_TTbar, model_histpdf_STop, model_histpdf_VV, model_histpdf_WJets), RooArgList(number_TTbar, number_STop, number_VV, number_WJets) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets);


        print " ############################## Fail: dataset #################################### "
        rdataset_data_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_fail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_fail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_fail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_fail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_fail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_fail)

        print " ############################### Fail: model ################################### "
        model_histpdf_TTbar_fail = self.workspace4fit_.pdf(rdataset_TTbar_mj_fail.GetName()+"_histpdf")
        model_histpdf_STop_fail = self.workspace4fit_.pdf(rdataset_STop_mj_fail.GetName()+"_histpdf")
        model_histpdf_VV_fail = self.workspace4fit_.pdf(rdataset_VV_mj_fail.GetName()+"_histpdf")
        model_histpdf_WJets_fail = self.workspace4fit_.pdf(rdataset_WJets_mj_fail.GetName()+"_histpdf")
 
        number_TTbar_fail = RooRealVar("rrv_number_TTbar_fail"+label+"_"+self.channel ,"rrv_number_TTbar_fail"+label+"_"+self.channel,rdataset_TTbar_mj_fail.sumEntries());
        number_STop_fail = RooRealVar("rrv_number_STop_fail"+label+"_"+self.channel ,"rrv_number_STop_fail"+label+"_"+self.channel,rdataset_STop_mj_fail.sumEntries());
        number_VV_fail = RooRealVar("rrv_number_VV_fail"+label+"_"+self.channel ,"rrv_number_VV_fail"+label+"_"+self.channel,rdataset_VV_mj_fail.sumEntries());
        number_WJets_fail = RooRealVar("rrv_number_WJets_fail"+label+"_"+self.channel,"rrv_number_WJets_fail"+label+"_"+self.channel,rdataset_WJets_mj_fail.sumEntries());

        model_TTbar_STop_VV_WJets_fail = RooAddPdf("model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel,"model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel, RooArgList(model_histpdf_TTbar_fail, model_histpdf_STop_fail, model_histpdf_VV_fail, model_histpdf_WJets_fail), RooArgList(number_TTbar_fail, number_STop_fail, number_VV_fail, number_WJets_fail) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets_fail);

        ##### inclusive scale factor from total integral

        scale_number_TTbar_STop_VV_WJets=(rdataset_TTbar_mj.sumEntries()+rdataset_STop_mj.sumEntries()+rdataset_VV_mj.sumEntries() +rdataset_WJets_mj.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )
        scale_number_TTbar_STop_VV_WJets_fail=(rdataset_TTbar_mj_fail.sumEntries()+rdataset_STop_mj_fail.sumEntries()+rdataset_VV_mj_fail.sumEntries() +rdataset_WJets_mj_fail.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )

        rrv_scale_number_TTbar_STop_VV_WJets = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets"+label,"rrv_scale_number_TTbar_STop_VV_WJets"+label,scale_number_TTbar_STop_VV_WJets);
        rrv_scale_number_TTbar_STop_VV_WJets_fail = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label,"rrv_scale_number_TTbar_STop_VV_WJets_fail"+label,scale_number_TTbar_STop_VV_WJets_fail);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets_fail);


        #### All the shape parameters and normalization are fixed

        print " ############################### Pass: Single MC model ################################### "
        
        model_STop = self.get_STop_mj_Model("_STop"+label);
        model_VV = self.get_VV_mj_Model("_VV"+label);
        model_WJets = self.get_General_mj_Model("_WJets0"+label);

        print " ############################### Fail: Single MC model ################################### "

        model_STop_fail = self.get_STop_mj_Model("_STop"+label+"_"+"failtau2tau1cut");
        model_VV_fail = self.get_VV_mj_Model("_VV"+label+"_"+"failtau2tau1cut");
        model_WJets_fail = self.get_General_mj_Model("_WJets0"+label+"_"+"failtau2tau1cut");

        self.constrainslist_data=[];
        ### Model for unmatched events passing and failing the cut --> ErfExp
        print " ############################### Pass: Data Bkg ################################### "
        model_bkg_data = self.make_Model("_bkg_data"+label,"ErfExp_ttbar","_mj",self.constrainslist_data); ## all the parameters contrained within the assigned error
        print " ############################### Fail: Data Bkg ################################### "
        model_bkg_data_fail = self.make_Model("_bkg_data"+label+"_"+"failtau2tau1cut","ErfExp_ttbar_failtau2tau1cut","_mj",self.constrainslist_data); 
## all the parameters contrained within the assigned error

        ### Model for matched events passing and failing the cut --> 2Gaus_ttbar and GausChebychev_ttbar_failtau2tau1cut
        print " ############################### Pass: Data ttbar resonant component ################################### "
        model_ttbar_data = self.make_Model_for_ttbar_controlsample("_ttbar_data"+label,"2Gaus_ttbar","_mj",self.constrainslist_data,label); ## No constraint
        print " ############################### Fail: Data ttbar resonant component ################################### "
        model_ttbar_data_fail = self.make_Model_for_ttbar_controlsample("_ttbar_data"+label+"_"+"failtau2tau1cut","GausChebychev_ttbar_failtau2tau1cut","_mj",self.constrainslist_data,label); 
## No Constraint


        print " ############################### Total Data Pdf Fail ################################### "
        model_data_fail=RooAddPdf("model_data"+label+"_"+"failtau2tau1cut"+"_"+self.channel,"model_data+"+label+"_"+"failtau2tau1cut"+"_"+self.channel, RooArgList(model_ttbar_data_fail, model_bkg_data_fail, model_STop_fail, model_VV_fail, model_WJets_fail));
        print " ############################### Total Data Pdf Pass ################################### "
        model_data=RooAddPdf("model_data"+label+"_"+self.channel,"model_data"+label+"_"+self.channel, RooArgList(model_ttbar_data, model_bkg_data, model_STop, model_VV, model_WJets));

        getattr(self.workspace4fit_,"import")(model_data);
        getattr(self.workspace4fit_,"import")(model_data_fail);

        ### take the event category from the workspace
        print " ############################### Data Event Category ################################### "
        category_p_f=self.workspace4fit_.cat("category_p_f"+label+"_"+self.channel);

        ### Define the simultaneous fit
        simPdf_data = RooSimultaneous("simPdf_data"+label+"_"+self.channel,"simPdf_data"+label+"_"+self.channel,category_p_f);
        simPdf_data.addPdf(model_data,"pass");
        simPdf_data.addPdf(model_data_fail,"fail");
        simPdf_data.Print();

        combData_p_f_data= self.workspace4fit_.data("combData_p_f_data"+label+"_"+self.channel);## combined dataset
        combData_p_f_data.Print();

        print " ############################### Data Total Fit ################################### "
        pdfconstrainslist_data=RooArgSet("pdfconstrainslist_data"+label+"_"+self.channel);
        for i in range(len(self.constrainslist_data)):
            self.workspace4fit_.pdf(self.constrainslist_data[i]).Print();
            pdfconstrainslist_data.add(self.workspace4fit_.pdf(self.constrainslist_data[i]) );
            
        if fit_or_not :
            rfresult_data=simPdf_data.fitTo(combData_p_f_data,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
            rfresult_data=simPdf_data.fitTo(combData_p_f_data,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data),RooFit.Minimizer("Minuit2"), RooFit.SumW2Error(kTRUE));
            rfresult_data.Print();


        #fit TotalMC
        self.constrainslist_TotalMC=[];

        print " ############################### Pass: MC Bkg ################################### "
        model_bkg_TotalMC = self.make_Model("_bkg_TotalMC"+label,"ErfExp_ttbar","_mj",self.constrainslist_TotalMC);
        print " ############################### Fail: MC Bkg ################################### "
        model_bkg_TotalMC_fail = self.make_Model("_bkg_TotalMC"+label+"_"+"failtau2tau1cut","ErfExp_ttbar_failtau2tau1cut","_mj",self.constrainslist_TotalMC);

        print " ############################### Pass: MC ttbar ################################### "
        model_ttbar_TotalMC = self.make_Model_for_ttbar_controlsample("_ttbar_TotalMC"+label,"2Gaus_ttbar","_mj",self.constrainslist_TotalMC,label);
        print " ############################### Fail: MC ttbar ################################### "
        model_ttbar_TotalMC_fail = self.make_Model_for_ttbar_controlsample("_ttbar_TotalMC"+label+"_"+"failtau2tau1cut","GausChebychev_ttbar_failtau2tau1cut","_mj",self.constrainslist_TotalMC,label);

        print " ############################### Fail: MC total PDF ################################### "
        model_TotalMC_fail = RooAddPdf("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel,"model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel, RooArgList(model_ttbar_TotalMC_fail, model_bkg_TotalMC_fail, model_STop_fail, model_VV_fail, model_WJets_fail));
        print " ############################### Pass: MC total PDF ################################### "
        model_TotalMC=RooAddPdf("model_TotalMC"+label+"_"+self.channel,"model_TotalMC"+label+"_"+self.channel, RooArgList(model_ttbar_TotalMC, model_bkg_TotalMC, model_STop, model_VV, model_WJets));
        getattr(self.workspace4fit_,"import")(model_TotalMC_fail);
        getattr(self.workspace4fit_,"import")(model_TotalMC);


        print " ############################### Data MC Category ################################### "
        simPdf_TotalMC=RooSimultaneous("simPdf_TotalMC"+label+"_"+self.channel,"simPdf_TotalMC"+label+"_"+self.channel,category_p_f);
        simPdf_TotalMC.addPdf(model_TotalMC,"pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail,"fail");
        
        combData_p_f_TotalMC=self.workspace4fit_.data("combData_p_f_TotalMC"+label+"_"+self.channel);


        print " ############################### MC Total Fit ################################### "
        pdfconstrainslist_TotalMC=RooArgSet("pdfconstrainslist_TotalMC"+label+"_"+self.channel);
        for i in range(len(self.constrainslist_TotalMC)):
            self.workspace4fit_.pdf(self.constrainslist_TotalMC[i]).Print();
            pdfconstrainslist_TotalMC.add(self.workspace4fit_.pdf(self.constrainslist_TotalMC[i]) );
        if fit_or_not :
            rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_p_f_TotalMC,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC), RooFit.Minimizer("Minuit2"));
            rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_p_f_TotalMC,RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC), RooFit.Minimizer("Minuit2"));
            rfresult_TotalMC.Print();


        xframe_data = rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_fail = rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_extremefail=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));

        ### plot data pass and fail
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        
        ### plot total pass mc normalizing it to data, passing + failing
        if TString(label).Contains("herwig"):
         model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar_herwig"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        else:
         model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(),model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        ### plot total fail mc normalizing it to data, passing + failing

        if TString(label).Contains("herwig"):
         model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar_herwig"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        else:
         model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())


        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        ## plot data again
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## plot mc fit function
        combData_p_f_TotalMC.plotOn(xframe_data,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible() , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("total MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
        simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit bkg"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_TotalMC.GetName(),model_STop.GetName(),model_VV.GetName(),model_WJets.GetName()) ), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
              
        ## plot data fit function
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0));

        simPdf_data.plotOn(xframe_data,RooFit.Name("total data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data,RooFit.Name("data fit bkg"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_data.GetName(),model_STop.GetName(),model_VV.GetName(),model_WJets.GetName())), RooFit.LineColor(kRed))
            
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        #fail plots -> plot MC fit
        combData_p_f_TotalMC.plotOn(xframe_data_fail,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("total MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
        simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit bkg"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_TotalMC_fail.GetName(),model_STop_fail.GetName(),model_VV_fail.GetName(),model_WJets_fail.GetName())), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))


        #fail plots -> plot data fit
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("total data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit bkg"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_data_fail.GetName(),model_STop_fail.GetName(),model_VV_fail.GetName(),model_WJets_fail.GetName())), RooFit.LineColor(kRed))

        #signal window
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data.GetMaximum()*0.7); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data.GetMaximum()*0.7); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
        xframe_data.addObject(lowerLine); xframe_data.addObject(upperLine);
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data_fail.GetMaximum()*0.7); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data_fail.GetMaximum()*0.7); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
        xframe_data_fail.addObject(lowerLine); xframe_data_fail.addObject(upperLine);

        # make the legend
        leg_data = self.legend4Plot(xframe_data,0,1,0.12)
        xframe_data.addObject(leg_data)
        leg_data_fail=self.legend4Plot(xframe_data_fail,0,1,0.12)
        xframe_data_fail.addObject(leg_data_fail)

        #add mean and width
        rrv_mean_gaus_data = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+self.channel);
        rrv_sigma_gaus_data = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_"+self.channel);
        rrv_mean_gaus_TotalMC = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+self.channel);
        rrv_sigma_gaus_TotalMC = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+self.channel);

        if rrv_mean_gaus_TotalMC:
            tl_MC_mean =TLatex(0.25 ,0.62, ("Mean_{MC } = %3.1f #pm %2.1f")%(rrv_mean_gaus_TotalMC.getVal(), rrv_mean_gaus_TotalMC.getError()) );
            tl_MC_sigma =TLatex(0.25 ,0.57, ("Sigma_{MC }= %2.1f #pm %2.1f")%(rrv_sigma_gaus_TotalMC.getVal(), rrv_sigma_gaus_TotalMC.getError()) );
            tl_MC_mean.SetNDC(); tl_MC_sigma.SetNDC();
            tl_MC_mean.SetTextSize(0.03)
            tl_MC_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_MC_mean);
            # xframe_data.addObject(tl_MC_sigma);

        if rrv_mean_gaus_data:
            tl_data_mean =TLatex(0.25 ,0.52, ("Mean_{data} = %3.1f #pm %2.1f")%(rrv_mean_gaus_data.getVal(), rrv_mean_gaus_data.getError()) );
            tl_data_sigma=TLatex(0.25 ,0.47, ("Sigma_{data}= %2.1f #pm %2.1f")%(rrv_sigma_gaus_data.getVal(), rrv_sigma_gaus_data.getError()) );
            tl_data_mean.SetNDC(); tl_data_sigma.SetNDC();
            tl_data_mean.SetTextSize(0.03)
            tl_data_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_data_mean);
            # xframe_data.addObject(tl_data_sigma);
         
        xframe_data.GetYaxis().SetRangeUser(1e-2,xframe_data.GetMaximum()*1.4);
        xframe_data_fail.GetYaxis().SetRangeUser(1e-2,xframe_data_fail.GetMaximum()*1.4);


        self.draw_canvas(xframe_data,"plots_%s_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation, self.channel, self.wtagger_label, self.wtagger_label),"control%s_%s_%s_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
        self.draw_canvas(xframe_data_fail,"plots_%s_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_fail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));


        self.ShowParam_Pdf(simPdf_data,RooArgSet(rrv_mass_j,category_p_f));
        self.ShowParam_Pdf(simPdf_TotalMC,RooArgSet(rrv_mass_j,category_p_f));
        
        if fit_or_not:
            rfresult_TotalMC.covarianceMatrix().Print();
            rfresult_data.covarianceMatrix().Print();
            rfresult_TotalMC.Print();
            rfresult_data.Print();

        self.ShowParam_Pdf(simPdf_data,RooArgSet(rrv_mass_j,category_p_f));


        # SF for LP
        print "########################## Extreme Fail Analysis ############################## "

        rdataset_data_mj_extremefail = self.workspace4fit_.data("rdataset_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_TotalMC_mj_extremefail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_extremefail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_extremefail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_extremefail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_extremefail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"extremefailtau2tau1cut_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_extremefail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_extremefail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_extremefail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_extremefail)

        model_histpdf_TTbar_extremefail = self.workspace4fit_.pdf(rdataset_TTbar_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_STop_extremefail = self.workspace4fit_.pdf(rdataset_STop_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_VV_extremefail = self.workspace4fit_.pdf(rdataset_VV_mj_extremefail.GetName()+"_histpdf")
        model_histpdf_WJets_extremefail = self.workspace4fit_.pdf(rdataset_WJets_mj_extremefail.GetName()+"_histpdf")

        number_TTbar_extremefail = RooRealVar("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+"extremefail"+"_"+self.channel,rdataset_TTbar_mj_extremefail.sumEntries());
        number_STop_extremefail = RooRealVar("rrv_number_STop"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_STop"+label+"_"+"extremefail"+"_"+self.channel,rdataset_STop_mj_extremefail.sumEntries());
        number_VV_extremefail = RooRealVar("rrv_number_VV"+label+"_"+"extremefail"+"_"+self.channel ,"rrv_number_VV"+label+"_"+"extremefail"+"_"+self.channel,rdataset_VV_mj_extremefail.sumEntries());
        number_WJets_extremefail = RooRealVar("rrv_number_WJets"+label+"_"+"extremefail"+"_"+self.channel,"rrv_number_WJets"+label+"_"+"extremefail"+"_"+self.channel,rdataset_WJets_mj_extremefail.sumEntries());
        
        model_TTbar_STop_VV_WJets_extremefail = RooAddPdf("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail"+"_"+self.channel,"model_TTbar_STop_VV_WJets"+label+"_"+"extremefail""_"+self.channel, RooArgList(model_histpdf_TTbar_extremefail, model_histpdf_STop_extremefail, model_histpdf_VV_extremefail, model_histpdf_WJets_extremefail), RooArgList(number_TTbar_extremefail, number_STop_extremefail, number_VV_extremefail, number_WJets_extremefail) );
        getattr(self.workspace4fit_,"import")(model_TTbar_STop_VV_WJets_extremefail);

        scale_number_TTbar_STop_VV_WJets_extremefail = (rdataset_TTbar_mj_extremefail.sumEntries()+rdataset_STop_mj_extremefail.sumEntries()+rdataset_VV_mj_extremefail.sumEntries() +rdataset_WJets_mj_extremefail.sumEntries() )/( rdataset_data_mj_extremefail.sumEntries() )
        rrv_scale_number_TTbar_STop_VV_WJets_extremefail = RooRealVar("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail","rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail",scale_number_TTbar_STop_VV_WJets_extremefail);
        getattr(self.workspace4fit_,"import")(rrv_scale_number_TTbar_STop_VV_WJets_extremefail);

        model_STop_extremefail = self.get_STop_mj_Model("_STop"+label+"_"+"extremefailtau2tau1cut");
        model_VV_extremefail = self.get_VV_mj_Model("_VV"+label+"_"+"extremefailtau2tau1cut");
        model_WJets_extremefail = self.get_General_mj_Model("_WJets0"+label+"_"+"extremefailtau2tau1cut");
                                                   

        tmp_constrainslist_data_LP=[];
        model_ttbar_data_extremefail = self.make_Model("_ttbar_data"+label+"_"+"extremefailtau2tau1cut","Exp_ttbar_extremefailtau2tau1cut","_mj",tmp_constrainslist_data_LP,50.);
        model_bkg_data_extremefail = self.make_Model("_bkg_data"+label+"_"+"extremefailtau2tau1cut","Exp_bkg_extremefailtau2tau1cut","_mj",tmp_constrainslist_data_LP,50);
        ## make the model for the extreme fail fit on data
        model_data_extremefail = RooAddPdf("model_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel,"model_data"+label+"_"+"extremefailtau2tau1cut_"+self.channel, RooArgList(model_ttbar_data_extremefail, model_bkg_data_extremefail, model_STop_extremefail, model_VV_extremefail, model_WJets_extremefail));
        getattr(self.workspace4fit_,"import")(model_data_extremefail);

        ## constraint the parameter in the extreme fail fit
        pdfconstrainslist_data_LP = RooArgSet("pdfconstrainslist_data_LP"+label+"_"+self.channel);
        for i in range(len(tmp_constrainslist_data_LP)):
         self.workspace4fit_.pdf(tmp_constrainslist_data_LP[i]).Print();
         pdfconstrainslist_data_LP.add(self.workspace4fit_.pdf(tmp_constrainslist_data_LP[i]) );

        if fit_or_not :
          rfresult_data_LP = model_data_extremefail.fitTo(rdataset_data_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data_LP),RooFit.Minimizer("Minuit2"));
          rfresult_data_LP=model_data_extremefail.fitTo(rdataset_data_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_data_LP),RooFit.Minimizer("Minuit2"));
          rfresult_data_LP.Print();

          rrv_number_bkg_data_extremefailtau2tau1cut = self.workspace4fit_.var("rrv_number_bkg_data%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
          rrv_number_ttbar_data_extremefailtau2tau1cut = self.workspace4fit_.var("rrv_number_ttbar_data%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel))
          rrv_number_ttbar_data_extremefailtau2tau1cut.setError( (rrv_number_ttbar_data_extremefailtau2tau1cut.getVal()+rrv_number_bkg_data_extremefailtau2tau1cut.getVal())/2. );
          rrv_number_bkg_data_extremefailtau2tau1cut.setError( (rrv_number_ttbar_data_extremefailtau2tau1cut.getVal()+rrv_number_bkg_data_extremefailtau2tau1cut.getVal())/2. );

          model_STop_extremefail.Print();
          model_VV_extremefail.Print();
          model_WJets_extremefail.Print();
          rdataset_data_mj_extremefail.Print();
          rdataset_TTbar_mj_extremefail.Print();
          rdataset_STop_mj_extremefail.Print();
          rdataset_VV_mj_extremefail.Print();
          rdataset_WJets_mj_extremefail.Print();

        ## extremefail fit of the mc
        tmp_constrainslist_TotalMC_LP=[];
        model_ttbar_TotalMC_extremefail = self.make_Model("_ttbar_TotalMC"+label+"_"+"extremefailtau2tau1cut","Exp_ttbar_extremefailtau2tau1cut","_mj",tmp_constrainslist_TotalMC_LP, 100.);
        model_bkg_TotalMC_extremefail = self.make_Model("_bkg_TotalMC"+label+"_"+"extremefailtau2tau1cut","Exp_bkg_extremefailtau2tau1cut","_mj",tmp_constrainslist_TotalMC_LP, 0.);

        model_TotalMC_extremefail=RooAddPdf("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel,"model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+self.channel, RooArgList(model_ttbar_TotalMC_extremefail, model_bkg_TotalMC_extremefail, model_STop_extremefail, model_VV_extremefail, model_WJets_extremefail));
        getattr(self.workspace4fit_,"import")(model_TotalMC_extremefail);

        pdfconstrainslist_TotalMC_LP=RooArgSet("pdfconstrainslist_TotalMC_LP_"+self.channel);
        for i in range(len(tmp_constrainslist_TotalMC_LP)):
            self.workspace4fit_.pdf(tmp_constrainslist_TotalMC_LP[i]).Print();
            pdfconstrainslist_TotalMC_LP.add(self.workspace4fit_.pdf(tmp_constrainslist_TotalMC_LP[i]) );

        if fit_or_not :
                                 
           rfresult_TotalMC_LP=model_TotalMC_extremefail.fitTo(rdataset_TotalMC_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_LP),RooFit.Minimizer("Minuit2"));
           rfresult_TotalMC_LP=model_TotalMC_extremefail.fitTo(rdataset_TotalMC_mj_extremefail, RooFit.Save(kTRUE), RooFit.Range("controlsample_fitting_range"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_LP),RooFit.Minimizer("Minuit2"));
           rfresult_TotalMC_LP.Print();
           rrv_number_bkg_TotalMC_extremefailtau2tau1cut=self.workspace4fit_.var("rrv_number_bkg_TotalMC%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
           rrv_number_ttbar_TotalMC_extremefailtau2tau1cut=self.workspace4fit_.var("rrv_number_ttbar_TotalMC%s_extremefailtau2tau1cut_%s_mj"%(label,self.channel));
           rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut.getVal())/2. );
           rrv_number_bkg_TotalMC_extremefailtau2tau1cut.setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut.getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut.getVal())/2. );
           model_STop_extremefail.Print();
           model_VV_extremefail.Print();
           model_WJets_extremefail.Print();
           rdataset_TotalMC_mj_extremefail.Print();
                                                        
        ## plot the extreme fail fit result + distribution

        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail, RooFit.Name("data invisi"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        if TString(label).Contains("herwig"):
         model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("TTbar_herwig"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines()) 
        else:
         model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines()) 

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_extremefail.GetName(), model_histpdf_VV_extremefail.GetName(), model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_extremefail.GetName(), model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_extremefail.plotOn(xframe_data_extremefail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_extremefail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())
        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail, RooFit.Name("data"), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        rdataset_TotalMC_mj_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("TotalMC invisi"), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );


        model_TotalMC_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("total MC fit"),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))

        model_TotalMC_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("MC fit bkg"),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_TotalMC_extremefail.GetName(),model_STop_extremefail.GetName(),model_VV_extremefail.GetName(),model_WJets_extremefail.GetName())), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))


        rdataset_data_mj_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data invisi"), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("total data fit"),RooFit.NormRange("controlsample_fitting_range"))
        model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit invisi"),RooFit.NormRange("controlsample_fitting_range"))
        model_data_extremefail.plotOn(xframe_data_extremefail,RooFit.Name("data fit bkg"),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%(model_bkg_data_extremefail.GetName(),model_STop_extremefail.GetName(),model_VV_extremefail.GetName(),model_WJets_extremefail.GetName())), RooFit.LineColor(kRed))


              
        xframe_data_extremefail.GetYaxis().SetRangeUser(1e-2,xframe_data_extremefail.GetMaximum()*1.5);

        leg_data_extremefail = self.legend4Plot(xframe_data_extremefail,0,1,0.12)
        xframe_data_extremefail.addObject(leg_data_extremefail)

        self.draw_canvas(xframe_data_extremefail,"plots_%s_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_extremefail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
                                                                  
  
    ### function to be used for the simultaneous fit of electron and muon channel --> just drawing the results
    def draw_ScaleFactor_forPureWJet_TTbar_controlsample(self,in_file_name, label="",fit_or_not=1):

        #dataset after tau2tau1 cut
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");
        rdataset_data_mj = self.workspace4fit_.data("rdataset_data"+label+"_"+self.channel+"_mj");
        rdataset_TTbar_mj = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+self.channel+"_mj");
        rdataset_STop_mj = self.workspace4fit_.data("rdataset_STop"+label+"_"+self.channel+"_mj");
        rdataset_VV_mj = self.workspace4fit_.data("rdataset_VV"+label+"_"+self.channel+"_mj");
        rdataset_WJets_mj = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+self.channel+"_mj");

        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj)

        model_histpdf_TTbar = self.workspace4fit_.pdf(rdataset_TTbar_mj.GetName()+"_histpdf")
        model_histpdf_STop = self.workspace4fit_.pdf(rdataset_STop_mj.GetName()+"_histpdf")
        model_histpdf_VV = self.workspace4fit_.pdf(rdataset_VV_mj.GetName()+"_histpdf")
        model_histpdf_WJets = self.workspace4fit_.pdf(rdataset_WJets_mj.GetName()+"_histpdf")

        number_TTbar = RooRealVar("rrv_number_TTbar"+label+"_"+self.channel ,"rrv_number_TTbar"+label+"_"+self.channel,rdataset_TTbar_mj.sumEntries());
        number_STop = RooRealVar("rrv_number_STop"+label+"_"+self.channel ,"rrv_number_STop"+label+"_"+self.channel,rdataset_STop_mj.sumEntries());
        number_VV = RooRealVar("rrv_number_VV"+label+"_"+self.channel ,"rrv_number_VV"+label+"_"+self.channel,rdataset_VV_mj.sumEntries());
        number_WJets = RooRealVar("rrv_number_WJets"+label+"_"+self.channel,"rrv_number_WJets"+label+"_"+self.channel,rdataset_WJets_mj.sumEntries());

        model_TTbar_STop_VV_WJets = self.workspace4fit_.pdf("model_TTbar_STop_VV_WJets"+label+"_"+self.channel);

        #dataset fail tau2tau1 cut
        rdataset_data_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_TTbar_mj_fail = self.workspace4fit_.data("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_STop_mj_fail = self.workspace4fit_.data("rdataset_STop"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_VV_mj_fail = self.workspace4fit_.data("rdataset_VV"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
        rdataset_WJets_mj_fail = self.workspace4fit_.data("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+self.channel+"_mj");
 
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_TTbar_mj_fail);
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_STop_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_VV_mj_fail)
        self.change_dataset_to_histpdf(rrv_mass_j, rdataset_WJets_mj_fail)

        model_histpdf_TTbar_fail = self.workspace4fit_.pdf(rdataset_TTbar_mj_fail.GetName()+"_histpdf")
        model_histpdf_STop_fail = self.workspace4fit_.pdf(rdataset_STop_mj_fail.GetName()+"_histpdf")
        model_histpdf_VV_fail = self.workspace4fit_.pdf(rdataset_VV_mj_fail.GetName()+"_histpdf")
        model_histpdf_WJets_fail = self.workspace4fit_.pdf(rdataset_WJets_mj_fail.GetName()+"_histpdf")

        number_TTbar_fail = RooRealVar("rrv_number_TTbar_fail"+label+"_"+self.channel ,"rrv_number_TTbar_fail"+label+"_"+self.channel,rdataset_TTbar_mj_fail.sumEntries());
        number_STop_fail = RooRealVar("rrv_number_STop_fail"+label+"_"+self.channel ,"rrv_number_STop_fail"+label+"_"+self.channel,rdataset_STop_mj_fail.sumEntries());
        number_VV_fail = RooRealVar("rrv_number_VV_fail"+label+"_"+self.channel ,"rrv_number_VV_fail"+label+"_"+self.channel,rdataset_VV_mj_fail.sumEntries());
        number_WJets_fail = RooRealVar("rrv_number_WJets_fail"+label+"_"+self.channel,"rrv_number_WJets_fail"+label+"_"+self.channel,rdataset_WJets_mj_fail.sumEntries());

        model_TTbar_STop_VV_WJets_fail = self.workspace4fit_.pdf("model_TTbar_STop_VV_WJets_fail"+label+"_"+self.channel);

        ## rescale factor for plots
        scale_number_TTbar_STop_VV_WJets = (rdataset_TTbar_mj.sumEntries()+rdataset_STop_mj.sumEntries()+rdataset_VV_mj.sumEntries() +rdataset_WJets_mj.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )
        scale_number_TTbar_STop_VV_WJets_fail = (rdataset_TTbar_mj_fail.sumEntries()+rdataset_STop_mj_fail.sumEntries()+rdataset_VV_mj_fail.sumEntries() +rdataset_WJets_mj_fail.sumEntries() )/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries() )


        model_STop  = self.get_STop_mj_Model("_STop"+label);
        model_VV    = self.get_VV_mj_Model("_VV"+label);
        model_WJets = self.get_General_mj_Model("_WJets0"+label);

        model_STop_fail  = self.get_STop_mj_Model("_STop"+label+"_"+"failtau2tau1cut");
        model_VV_fail    = self.get_VV_mj_Model("_VV"+label+"_"+"failtau2tau1cut");
        model_WJets_fail = self.get_General_mj_Model("_WJets0"+label+"_"+"failtau2tau1cut");

        model_data_fail = self.workspace4fit_.pdf("model_data"+label+"_"+"failtau2tau1cut"+"_"+self.channel)
        model_data      = self.workspace4fit_.pdf("model_data"+label+"_"+self.channel);

        category_p_f = self.workspace4fit_.cat("category_p_f"+label+"_"+self.channel);

        simPdf_data = RooSimultaneous("simPdf_data"+label+"_"+self.channel,"simPdf_data"+label+"_"+self.channel,category_p_f);
        simPdf_data.addPdf(model_data,"pass");
        simPdf_data.addPdf(model_data_fail,"fail");
        combData_p_f_data = self.workspace4fit_.data("combData_p_f_data"+label+"_"+self.channel);

        simPdf_data.Print();
        combData_p_f_data.Print();

        model_TotalMC_fail = self.workspace4fit_.pdf("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+self.channel);
        model_TotalMC = self.workspace4fit_.pdf("model_TotalMC"+label+"_"+self.channel);

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC"+label+"_"+self.channel,"simPdf_TotalMC"+label+"_"+self.channel,category_p_f);
        simPdf_TotalMC.addPdf(model_TotalMC,"pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail,"fail");
        combData_p_f_TotalMC = self.workspace4fit_.data("combData_p_f_TotalMC"+label+"_"+self.channel);


        xframe_data=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));
        xframe_data_fail=rrv_mass_j.frame( RooFit.Bins(int(rrv_mass_j.getBins())));

        ## plot data and mc on the frame
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        
        #pass plot
        if TString(label).Contains("herwig"):     
         model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar_herwig"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        else:
         model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar"),RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #pass plots
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("STop_line_invisible"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("VV_line_invisible"),RooFit.Components("%s,%s"%(model_histpdf_VV.GetName(), model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets.plotOn(xframe_data,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),RooFit.Name("WJets_line_invisible"),RooFit.Components("%s"%(model_histpdf_WJets.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        #fail plot
        if TString(label).Contains("herwig"):     
         model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar_herwig"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        else:
         model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar"), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["TTbar"+label]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["STop"]), RooFit.LineColor(kBlack), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["VV"]), RooFit.LineColor(kBlack), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.DrawOption("F"), RooFit.FillColor(self.color_palet["WJets"]), RooFit.LineColor(kBlack), RooFit.VLines())

        #solid line
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("TTbar_line_invisible"), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())
        
        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("STop_line_invisible"),RooFit.Components("%s,%s,%s"%(model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("VV_line_invisible"),RooFit.Components("%s,%s"%(model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())

        model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit.Name("WJets_line_invisible"),RooFit.Components("%s"%( model_histpdf_WJets_fail.GetName()) ), RooFit.LineColor(kBlack), RooFit.LineWidth(2), RooFit.VLines())


        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );
        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        ## plot mc fit
        combData_p_f_TotalMC.plotOn(xframe_data,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible() , RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("total MC fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
        simPdf_TotalMC.plotOn(xframe_data,RooFit.Name("MC fit bkg"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC), RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%("model_bkg_TotalMC_%s_mj"%(self.channel),"model_STop_%s_mj"%(self.channel),"model_VV_%s_mj"%(self.channel),"model_WJets0_%s_mj"%(self.channel)) ), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))


        ## plot data fit
        combData_p_f_data.plotOn(xframe_data,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::pass"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_data.plotOn(xframe_data,RooFit.Name("total data fit"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data,RooFit.Name("data fit bkg"),RooFit.Slice(category_p_f,"pass"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%("model_bkg_data_%s_mj"%(self.channel),"model_STop_%s_mj"%(self.channel),"model_VV_%s_mj"%(self.channel),"model_WJets0_%s_mj"%(self.channel))), RooFit.LineColor(kRed))
         

        combData_p_f_TotalMC.plotOn(xframe_data_fail,RooFit.Name("TotalMC invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerColor(kWhite), RooFit.LineColor(kWhite), RooFit.Invisible(), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("total MC fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.LineStyle(kDashed))
        simPdf_TotalMC.plotOn(xframe_data_fail,RooFit.Name("MC fit bkg"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_TotalMC),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%("model_bkg_TotalMC_failtau2tau1cut_%s_mj"%(self.channel),"model_STop_failtau2tau1cut_%s_mj"%(self.channel),"model_VV_failtau2tau1cut_%s_mj"%(self.channel),"model_WJets0_failtau2tau1cut_%s_mj"%(self.channel))), RooFit.LineColor(kRed), RooFit.LineStyle(kDashed))
            

        combData_p_f_data.plotOn(xframe_data_fail,RooFit.Name("data invisi"), RooFit.Cut("category_p_f%s_%s==category_p_f%s_%s::fail"%(label,self.channel,label,self.channel)), RooFit.MarkerSize(1.5), RooFit.DataError(RooAbsData.SumW2), RooFit.XErrorSize(0) );

        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("total data fit"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit invisi"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"))
        simPdf_data.plotOn(xframe_data_fail,RooFit.Name("data fit bkg"),RooFit.Slice(category_p_f,"fail"),RooFit.ProjWData(RooArgSet(category_p_f),combData_p_f_data),RooFit.NormRange("controlsample_fitting_range"), RooFit.Components("%s,%s,%s,%s"%("model_bkg_data_failtau2tau1cut_%s_mj"%(self.channel),"model_STop_failtau2tau1cut_%s_mj"%(self.channel),"model_VV_failtau2tau1cut_%s_mj"%(self.channel),"model_WJets0_failtau2tau1cut_%s_mj"%(self.channel))), RooFit.LineColor(kRed))
      
        #signal window
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data.GetMaximum()*0.7); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data.GetMaximum()*0.7); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
        xframe_data.addObject(lowerLine); xframe_data.addObject(upperLine);
        lowerLine = TLine(self.mj_signal_min,0.,self.mj_signal_min,xframe_data_fail.GetMaximum()*0.7); lowerLine.SetLineWidth(2); lowerLine.SetLineColor(kBlack); lowerLine.SetLineStyle(9);
        upperLine = TLine(self.mj_signal_max,0.,self.mj_signal_max,xframe_data_fail.GetMaximum()*0.7); upperLine.SetLineWidth(2); upperLine.SetLineColor(kBlack); upperLine.SetLineStyle(9);
        xframe_data_fail.addObject(lowerLine); xframe_data_fail.addObject(upperLine);

        #legend
        leg_data=self.legend4Plot(xframe_data,0,1, 0.12)
        xframe_data.addObject(leg_data)
        leg_data_fail=self.legend4Plot(xframe_data_fail,0,1,0.12)
        xframe_data_fail.addObject(leg_data_fail)

        #add mean and width
        tmp_channel="el";
        if self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+"el"): tmp_channel="el";
        else: tmp_channel="mu";

        rrv_mean_gaus_data =self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_"+tmp_channel);
        rrv_sigma_gaus_data =self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_"+tmp_channel);
        rrv_mean_gaus_TotalMC =self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel);
        rrv_sigma_gaus_TotalMC=self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel);

        if rrv_mean_gaus_TotalMC:
            tl_MC_mean =TLatex(0.25 ,0.62, ("Mean_{MC } = %3.1f #pm %2.1f")%(rrv_mean_gaus_TotalMC.getVal(), rrv_mean_gaus_TotalMC.getError()) );
            tl_MC_sigma =TLatex(0.25 ,0.57, ("Sigma_{MC }= %2.1f #pm %2.1f")%(rrv_sigma_gaus_TotalMC.getVal(), rrv_sigma_gaus_TotalMC.getError()) );
            tl_MC_mean.SetNDC(); tl_MC_sigma.SetNDC();
            tl_MC_mean.SetTextSize(0.03)
            tl_MC_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_MC_mean);
            # xframe_data.addObject(tl_MC_sigma);

        if rrv_mean_gaus_data:
            tl_data_mean =TLatex(0.25 ,0.52, ("Mean_{data} = %3.1f #pm %2.1f")%(rrv_mean_gaus_data.getVal(), rrv_mean_gaus_data.getError()) );
            tl_data_sigma=TLatex(0.25 ,0.47, ("Sigma_{data}= %2.1f #pm %2.1f")%(rrv_sigma_gaus_data.getVal(), rrv_sigma_gaus_data.getError()) );
            tl_data_mean.SetNDC(); tl_data_sigma.SetNDC();
            tl_data_mean.SetTextSize(0.03)
            tl_data_sigma.SetTextSize(0.03)
            # xframe_data.addObject(tl_data_mean); xframe_data.addObject(tl_data_sigma);
         
        xframe_data.GetYaxis().SetRangeUser(1e-2,xframe_data.GetMaximum()*1.2);
        xframe_data_fail.GetYaxis().SetRangeUser(1e-2,xframe_data_fail.GetMaximum()*1.2);

        self.draw_canvas(xframe_data,"plots_%s_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation, self.channel,self.wtagger_label, self.wtagger_label),"control%s_%s_%s_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));
        self.draw_canvas(xframe_data_fail,"plots_%s_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/"%(options.additioninformation, self.channel, self.wtagger_label,self.wtagger_label),"control%s_%s_%s_fail_pTbin_%d_%d"%(label,self.wtagger_label,self.channel,self.ca8_ungroomed_pt_min,self.ca8_ungroomed_pt_max));


    ## To build the dataset to be fitted
    def get_mj_and_mlvj_dataset_TTbar_controlsample(self,in_file_name, label, jet_mass="ttb_ca8_mass_pr"):# to get the shape of m_lvj

        # read in tree
        fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");
        
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ### dataset of m_j before tau2tau1 cut : Passed
        rdataset_mj = RooDataSet("rdataset"+label+"_"+self.channel+"_mj","rdataset"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj = RooDataSet("rdataset4fit"+label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j before tau2tau1 cut : Total
        rdataset_beforetau2tau1cut_mj = RooDataSet("rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_beforetau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j failed tau2tau1 cut :
        rdataset_failtau2tau1cut_mj = RooDataSet("rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_failtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j extreme failed tau2tau1 cut: >0.75
        rdataset_extremefailtau2tau1cut_mj = RooDataSet("rdataset"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_extremefailtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

                 
        #### combine of dataset before and after tau2tau1 cut
        if TString(label).Contains("herwig"):
         category_cut = RooCategory("category_cut"+"_herwig_"+self.channel,"category_cut"+"_herwig_"+self.channel);
         category_cut.defineType("cut",1);
         category_cut.defineType("beforecut",2);
         combData4cut = RooDataSet("combData4cut"+"_herwig_"+self.channel,"combData4cut"+"_"+self.channel,RooArgSet(rrv_mass_j, category_cut, rrv_weight),RooFit.WeightVar(rrv_weight) );
        else:
         category_cut = RooCategory("category_cut"+"_"+self.channel,"category_cut"+"_"+self.channel);
         category_cut.defineType("cut",1);
         category_cut.defineType("beforecut",2);
         combData4cut = RooDataSet("combData4cut"+"_"+self.channel,"combData4cut"+"_"+self.channel,RooArgSet(rrv_mass_j, category_cut, rrv_weight),RooFit.WeightVar(rrv_weight) );
            
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
        
        combData_p_f=RooDataSet("combData_p_f"+label+"_"+self.channel,"combData_p_f"+label+"_"+self.channel,RooArgSet(rrv_mass_j, category_p_f, rrv_weight),RooFit.WeightVar(rrv_weight));
            
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

            if i==0: tmp_scale_to_lumi=treeIn.wSampleWeight ## weigth for xs for the mc
                
            discriminantCut = 0;
            wtagger = getattr(treeIn,"ttb_ca8_tau2tau1");
            
            if wtagger < self.wtagger_cut:
                discriminantCut=2;
            elif wtagger > self.wtagger_cut and wtagger < 0.75:
                discriminantCut=1;
            elif wtagger > 0.75 :
                iscriminantCut=0;

            tmp_jet_mass = 0. ;    
            if options.shift == 1 and not TString(label).Contains("data"):
                 tmp_jet_mass = getattr(treeIn, jet_mass) + self.mean_shift;
            elif options.smear == 1 and not TString(label).Contains("data"):
                 tmp_jet_mass = getattr(treeIn, jet_mass)*(1+rand.Gaus(0,self.sigma_scale-1));
            else:
                 tmp_jet_mass = getattr(treeIn, jet_mass);

            tmp_event_weight = getattr(treeIn,"totalEventWeight");
            tmp_event_weight4fit = getattr(treeIn,"eff_and_pu_Weight");
            tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"wSampleWeight")/tmp_scale_to_lumi

            if not TString(label).Contains("data"):
                  tmp_event_weight = tmp_event_weight*getattr(treeIn,"btag_weight");
                  tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"btag_weight");
            else:
                  tmp_event_weight = 1.;
                  tmp_event_weight4fit = 1.;


            ### Cut for the HP category
            if discriminantCut ==2 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2:

            # and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2.0 and ( getattr(treeIn,"mass_ungroomedjet_closerjet")<230 and getattr(treeIn,"mass_ungroomedjet_closerjet")>150) and (getattr(treeIn,"mass_leptonic_closerjet") < 230 and getattr(treeIn,"mass_leptonic_closerjet") > 130) :


                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight#sum_of_weight_HP/sum_of_event_HP;

                rrv_mass_j.setVal(tmp_jet_mass);

                rdataset_mj.add( RooArgSet(rrv_mass_j ),tmp_event_weight);
                rdataset4fit_mj.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);

                if tmp_jet_mass >=self.mj_sideband_lo_min and tmp_jet_mass <self.mj_sideband_lo_max:
                    hnum_4region.Fill(-1,tmp_event_weight );
                if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                    hnum_2region.Fill(1,tmp_event_weight);
                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                    hnum_4region.Fill(0,tmp_event_weight);
                    hnum_4region_error2.Fill(0,tmp_event_weight*tmp_event_weight);
                if tmp_jet_mass >=self.mj_sideband_hi_min and tmp_jet_mass <self.mj_sideband_hi_max:
                      hnum_4region.Fill(1,tmp_event_weight);

                hnum_4region.Fill(2,tmp_event_weight);

                category_cut.setLabel("cut");
                combData4cut.add(RooArgSet(rrv_mass_j,category_cut),tmp_event_weight4fit);
                category_p_f.setLabel("pass");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

            ### Cut for the Total category
            if getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2:

#        and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2.0  and ( getattr(treeIn,"mass_ungroomedjet_closerjet")<230 and getattr(treeIn,"mass_ungroomedjet_closerjet")>150) and (getattr(treeIn,"mass_leptonic_closerjet") < 230 and getattr(treeIn,"mass_leptonic_closerjet") > 130) :


                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight

                rrv_mass_j.setVal( tmp_jet_mass );

                if tmp_jet_mass >=self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                   hnum_4region_before_cut.Fill(0,tmp_event_weight);
                   hnum_4region_before_cut_error2.Fill(0,tmp_event_weight*tmp_event_weight);

                rdataset_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight);
                rdataset4fit_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight4fit);

                category_cut.setLabel("beforecut");
                combData4cut.add(RooArgSet(rrv_mass_j,category_cut),tmp_event_weight4fit);

            ### 1-HP category
            if (discriminantCut==1 or discriminantCut==0) and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2:

# and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2.0 and ( getattr(treeIn,"mass_ungroomedjet_closerjet")<230 and getattr(treeIn,"mass_ungroomedjet_closerjet")>150) and (getattr(treeIn,"mass_leptonic_closerjet") < 230 and getattr(treeIn,"mass_leptonic_closerjet") > 130):


                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;
            
                rrv_mass_j.setVal( tmp_jet_mass );

                rdataset_failtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_failtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                category_p_f.setLabel("fail");
                combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight);

            ### extreme fail category
            if discriminantCut==0 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and getattr(treeIn,"ttb_dR_ca8_jet_closer") < 2:

# and getattr(treeIn,"ttb_dR_ca8_bjet_closer") < 2.0 and ( getattr(treeIn,"mass_ungroomedjet_closerjet")<230 and getattr(treeIn,"mass_ungroomedjet_closerjet")>150) and (getattr(treeIn,"mass_leptonic_closerjet") < 230 and getattr(treeIn,"mass_leptonic_closerjet") > 130) :

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   Tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;

                rdataset_extremefailtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_extremefailtau2tau1cut_mj.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );
    
        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi);
        rrv_scale_to_lumi_failtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
        rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
                
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut)

        #prepare m_j dataset

        rrv_number_dataset_sb_lo_mj = RooRealVar("rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1));
        rrv_number_dataset_signal_region_mj = RooRealVar("rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2));
        rrv_number_dataset_signal_region_error2_mj = RooRealVar("rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj",hnum_4region_error2.GetBinContent(2));

        rrv_number_dataset_signal_region_before_cut_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut"+label+"_"+self.channel+"_mj",hnum_4region_before_cut.GetBinContent(2));
        rrv_number_dataset_signal_region_before_cut_error2_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj",hnum_4region_before_cut_error2.GetBinContent(2));
        rrv_number_dataset_sb_hi_mj = RooRealVar("rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"+label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3));

        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_error2_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_error2_mj)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj)
        getattr(self.workspace4fit_,"import")(combData_p_f);
        
        print "N_rdataset_mj: "
        getattr(self.workspace4fit_,"import")(rdataset_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj)
        getattr(self.workspace4fit_,"import")(rdataset_beforetau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_beforetau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset_failtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_failtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset_extremefailtau2tau1cut_mj)
        getattr(self.workspace4fit_,"import")(rdataset4fit_extremefailtau2tau1cut_mj)

                
        rdataset_mj.Print();
        rdataset4fit_mj.Print();
        rrv_number_dataset_sb_lo_mj.Print()
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_signal_region_error2_mj.Print()
        rrv_number_dataset_signal_region_before_cut_mj.Print()
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
        rrv_number_dataset_sb_hi_mj.Print()

        rdataset_mj.Print();
        rdataset_beforetau2tau1cut_mj.Print();
        rdataset_failtau2tau1cut_mj.Print();
        rdataset_extremefailtau2tau1cut_mj.Print();
        rrv_number_dataset_signal_region_mj.Print()
        rrv_number_dataset_signal_region_error2_mj.Print()
        rrv_number_dataset_signal_region_before_cut_mj.Print()
        rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
        combData_p_f.Print("v");

    ## To build the dataset to be fitted in case of mc ttbar matched truth analysis
    def get_mj_and_mlvj_dataset_TTbar_controlsample_matched_with_truth(self,in_file_name, label, jet_mass="ttb_ca8_mass_pr"):

        # read in tree
        fileIn_name = TString(options.inPath+"/"+self.file_Directory+in_file_name);
        fileIn = TFile(fileIn_name.Data());
        treeIn = fileIn.Get("otree");
        
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

        ### dataset of m_j : Passed matched dataset
        rdataset_mj_matched_pass = RooDataSet("rdataset"+label+"_matched_pass_"+self.channel+"_mj","rdataset"+label+"_matched_pass_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_matched_pass = RooDataSet("rdataset4fit"+label+"_matched_pass_"+self.channel+"_mj","rdataset4fit"+label+"_matched_pass_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j : Passed unmatched dataset
        rdataset_mj_unmatched_pass = RooDataSet("rdataset"+label+"_unmatched_pass_"+self.channel+"_mj","rdataset"+label+"_unmatched_pass_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_unmatched_pass = RooDataSet("rdataset4fit"+label+"_unmatched_pass_"+self.channel+"_mj","rdataset4fit"+label+"_unmatched_pass_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j : Fail matched dataset
        rdataset_mj_matched_fail = RooDataSet("rdataset"+label+"_matched_fail_"+self.channel+"_mj","rdataset"+label+"_matched_fail_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_matched_fail = RooDataSet("rdataset4fit"+label+"_matched_fail_"+self.channel+"_mj","rdataset4fit"+label+"_matched_fail_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j : Fail unmatched dataset
        rdataset_mj_unmatched_fail = RooDataSet("rdataset"+label+"_unmatched_fail_"+self.channel+"_mj","rdataset"+label+"_unmatched_fail_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_unmatched_fail = RooDataSet("rdataset4fit"+label+"_unmatched_fail_"+self.channel+"_mj","rdataset4fit"+label+"_unmatched_fail_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j : extremefail matched dataset
        rdataset_mj_extreme_matched = RooDataSet("rdataset"+label+"_extreme_matched_"+self.channel+"_mj","rdataset"+label+"_extreme_matched_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_extreme_matched = RooDataSet("rdataset4fit"+label+"_extreme_matched_"+self.channel+"_mj","rdataset4fit"+label+"_extreme_matched_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );

        ### dataset of m_j : Fail unmatched dataset
        rdataset_mj_extreme_unmatched = RooDataSet("rdataset"+label+"_extreme_unmatched_"+self.channel+"_mj","rdataset"+label+"_extreme_unmatched_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
        rdataset4fit_mj_extreme_unmatched = RooDataSet("rdataset4fit"+label+"_extreme_unmatched_"+self.channel+"_mj","rdataset4fit"+label+"_extreme_unmatched_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) );
                            
        print "N entries: ", treeIn.GetEntries()

        hnum_6region = TH1D("hnum_6region"+label+"_"+self.channel,"hnum_6region"+label+"_"+self.channel,5,0,6);# 0 pass matched, 1 pass unmatched, 2 fail matched, 3 fail unmatched, 4 extreme matched, 5 extreme unmatched
        hnum_6region_error2 = TH1D("hnum_6region_error2"+label+"_"+self.channel,"hnum_6region_error2"+label+"_"+self.channel,4,-1.5,2.5);

        for i in range(treeIn.GetEntries()):
            if i % 100000 == 0: print "iEntry: ",i

            treeIn.GetEntry(i);
            if i==0: tmp_scale_to_lumi=treeIn.wSampleWeight
                
            discriminantCut = 0;
            wtagger = getattr(treeIn,"ttb_ca8_tau2tau1");
            
            if wtagger < self.wtagger_cut:
                discriminantCut=2;
            elif wtagger > self.wtagger_cut and wtagger < 0.75:
                discriminantCut=1;
            elif wtagger > 0.75 :
                iscriminantCut=0;
                
            tmp_jet_mass=getattr(treeIn, jet_mass);

            tmp_event_weight = getattr(treeIn,"totalEventWeight");
            tmp_event_weight4fit = getattr(treeIn,"eff_and_pu_Weight");
            tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"wSampleWeight")/tmp_scale_to_lumi

            if not TString(label).Contains("data"):
                  tmp_event_weight = tmp_event_weight*getattr(treeIn,"btag_weight");
                  tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"btag_weight");
            else:
                  tmp_event_weight = 1.;
                  tmp_event_weight4fit = 1.;

            ########## Generator level matching -- shpuld be correct once we will have W parton in the trees ##########
          
            if getattr(treeIn,"gen_W_jet_pt") > 1000. or getattr(treeIn,"gen_W_jet_pt") < 0. : continue ;

            delta_eta = math.fabs(getattr(treeIn,"ttb_ca8_ungroomed_eta")-getattr(treeIn,"gen_W_jet_eta"));

            delta_phi = 0.;
            if(math.fabs(getattr(treeIn,"ttb_ca8_ungroomed_phi")-getattr(treeIn,"gen_W_jet_phi")) <= TMath.Pi() ):
               delta_phi = math.fabs(getattr(treeIn,"ttb_ca8_ungroomed_phi")-getattr(treeIn,"gen_W_jet_phi"));
            else:
               delta_phi = 2*TMath.Pi() - math.fabs(getattr(treeIn,"ttb_ca8_ungroomed_phi")-getattr(treeIn,"gen_W_jet_phi"));

            deltaR = TMath.Sqrt(delta_phi*delta_phi+delta_eta*delta_eta);

            ### Cut for the HP matched
            if discriminantCut ==2 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR <= 0.5:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight

                rrv_mass_j.setVal(tmp_jet_mass);

                rdataset_mj_matched_pass.add( RooArgSet(rrv_mass_j ),tmp_event_weight);
                rdataset4fit_mj_matched_pass.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);

                hnum_6region.Fill(0,tmp_event_weight)
                hnum_6region_error2.Fill(0,tmp_event_weight*tmp_event_weight)

            ### Cut for the HP unmatched
            if discriminantCut ==2 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR > 0.5:


                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight

                rrv_mass_j.setVal(tmp_jet_mass);

                rdataset_mj_unmatched_pass.add( RooArgSet(rrv_mass_j ),tmp_event_weight);
                rdataset4fit_mj_unmatched_pass.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);

                hnum_6region.Fill(1,tmp_event_weight)
                hnum_6region_error2.Fill(1,tmp_event_weight*tmp_event_weight)


            ### 1-HP category matched
            if (discriminantCut==1 or discriminantCut==0) and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR <= 0.5:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;
            
                rrv_mass_j.setVal( tmp_jet_mass );

                rdataset_mj_matched_fail.add( RooArgSet(rrv_mass_j ),tmp_event_weight);
                rdataset4fit_mj_matched_fail.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);

                hnum_6region.Fill(2,tmp_event_weight)
                hnum_6region_error2.Fill(2,tmp_event_weight*tmp_event_weight)

            ### 1-HP category matched
            if (discriminantCut==1 or discriminantCut==0) and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR > 0.5:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;
            
                rrv_mass_j.setVal( tmp_jet_mass );

                rdataset_mj_unmatched_fail.add( RooArgSet(rrv_mass_j ),tmp_event_weight);
                rdataset4fit_mj_unmatched_fail.add(RooArgSet( rrv_mass_j ),tmp_event_weight4fit);

                hnum_6region.Fill(3,tmp_event_weight)
                hnum_6region_error2.Fill(3,tmp_event_weight*tmp_event_weight)


            ### extreme fail category matched
            if discriminantCut==0 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR <= 0.5:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   Tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;

                rdataset_mj_extreme_matched.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_mj_extreme_matched.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                hnum_6region.Fill(4,tmp_event_weight)
                hnum_6region_error2.Fill(4,tmp_event_weight*tmp_event_weight)

                ### extreme fail category unmatched
            if discriminantCut==0 and getattr(treeIn,"mass_lvj_type0_met") < self.mass_lvj_max and getattr(treeIn,"mass_lvj_type0_met") > self.mass_lvj_min and (getattr(treeIn,"ttb_nak5_same_csvm") > 0 or getattr(treeIn,"ttb_nak5_oppoveto_csvm") > 0) and getattr(treeIn,"isttbar") > 0 and getattr(treeIn,"v_pt") > self.vpt_cut and getattr(treeIn,"l_pt") >= self.lpt_cut and getattr(treeIn,"pfMET") > self.pfMET_cut and getattr(treeIn,"ttb_ca8_ungroomed_pt") > 200 and tmp_jet_mass>rrv_mass_j.getMin() and tmp_jet_mass<rrv_mass_j.getMax() and getattr(treeIn,"ttb_ca8_ungroomed_pt") > self.ca8_ungroomed_pt_min and getattr(treeIn,"ttb_ca8_ungroomed_pt") < self.ca8_ungroomed_pt_max and deltaR > 0.5:

                if TString(label).Contains("herwig") and not TString(label).Contains("data") :
                   tmp_event_weight = tmp_event_weight*treeIn.event_weight;
                   Tmp_event_weight4fit = tmp_event_weight4fit*treeIn.event_weight;

                rdataset_mj_extreme_unmatched.add( RooArgSet( rrv_mass_j ), tmp_event_weight );
                rdataset4fit_mj_extreme_unmatched.add( RooArgSet( rrv_mass_j ), tmp_event_weight4fit );

                hnum_6region.Fill(5,tmp_event_weight)
                hnum_6region_error2.Fill(5,tmp_event_weight*tmp_event_weight)

        rrv_scale_to_lumi = RooRealVar("rrv_scale_to_lumi"+label+"_"+self.channel,"rrv_scale_to_lumi"+label+"_"+self.channel,tmp_scale_to_lumi);
        rrv_scale_to_lumi_failtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
        rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,tmp_scale_to_lumi);
                
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut)
        getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut)

        #prepare m_j dataset
        rrv_number_dataset_matched_pass = RooRealVar("rrv_number_dataset"+label+"_matched_pass_"+self.channel+"_mj","rrv_number_dataset"+label+"_matched_pass_"+self.channel+"_mj",hnum_6region.GetBinContent(1));
        rrv_number_dataset_unmatched_pass = RooRealVar("rrv_number_dataset"+label+"_"+"_unmatched_pass_"+self.channel+"_mj","rrv_number_dataset"+label+"_unmatched_pass_"+self.channel+"_mj",hnum_6region.GetBinContent(2));
        rrv_number_dataset_matched_fail = RooRealVar("rrv_number_dataset"+label+"_matched_fail_"+self.channel+"_mj","rrv_number_dataset"+label+"_matched_fail_"+self.channel+"_mj",hnum_6region.GetBinContent(3));
        rrv_number_dataset_unmatched_fail = RooRealVar("rrv_number_dataset"+label+"_"+"_unmatched_fail_"+self.channel+"_mj","rrv_number_dataset"+label+"_unmatched_fail_"+self.channel+"_mj",hnum_6region.GetBinContent(4));
        rrv_number_dataset_extreme_matched = RooRealVar("rrv_number_dataset"+label+"_extreme_matched_"+self.channel+"_mj","rrv_number_dataset"+label+"_extreme_matched"+self.channel+"_mj",hnum_6region.GetBinContent(5));
        rrv_number_dataset_extreme_unmatched = RooRealVar("rrv_number_dataset"+label+"_"+"_extreme_unmatched_"+self.channel+"_mj","rrv_number_dataset"+label+"_extreme_unmatched_"+self.channel+"_mj",hnum_6region.GetBinContent(5));


        getattr(self.workspace4fit_,"import")(rrv_number_dataset_matched_pass)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_unmatched_pass)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_matched_fail)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_unmatched_fail)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_extreme_matched)
        getattr(self.workspace4fit_,"import")(rrv_number_dataset_extreme_unmatched)
                                                        
        
        print "N_rdataset_mj: "
        getattr(self.workspace4fit_,"import")(rdataset_mj_matched_pass)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_matched_pass)
        getattr(self.workspace4fit_,"import")(rdataset_mj_matched_fail)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_matched_fail)
        getattr(self.workspace4fit_,"import")(rdataset_mj_unmatched_pass)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_unmatched_pass)
        getattr(self.workspace4fit_,"import")(rdataset_mj_unmatched_fail)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_unmatched_fail)
        getattr(self.workspace4fit_,"import")(rdataset_mj_extreme_matched)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_extreme_matched)
        getattr(self.workspace4fit_,"import")(rdataset_mj_extreme_unmatched)
        getattr(self.workspace4fit_,"import")(rdataset4fit_mj_extreme_unmatched)

        rdataset_mj_matched_pass.Print();
        rdataset4fit_mj_matched_pass.Print();
        rdataset_mj_unmatched_pass.Print();
        rdataset4fit_mj_unmatched_pass.Print();
        rdataset_mj_matched_fail.Print();
        rdataset4fit_mj_matched_fail.Print();
        rdataset_mj_unmatched_fail.Print();
        rdataset4fit_mj_unmatched_fail.Print();
        rdataset_mj_extreme_matched.Print();
        rdataset4fit_mj_extreme_matched.Print();
        rdataset_mj_extreme_unmatched.Print();
        rdataset4fit_mj_extreme_unmatched.Print();
    

    ### in order to get the pull
    def get_pull(self, rrv_x, mplot_orig):

       hpull = mplot_orig.pullHist();
       x = ROOT.Double(0.); y = ROOT.Double(0) ;
       for ipoint in range(0,hpull.GetN()):
        hpull.GetPoint(ipoint,x,y);
        if(y == 0):
          hpull.SetPoint(ipoint,x,10)
       
       mplot_pull = rrv_x.frame(RooFit.Title("Pull Distribution"), RooFit.Bins(int(rrv_x.getBins()/self.narrow_factor)));
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

    #### in order to make the banner on the plots
    def banner4Plot(self, iswithpull=0):

      if iswithpull:
       if self.channel=="el":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e#nu "%(self.GetLumi()));                   
       elif self.channel=="mu":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu#nu "%(self.GetLumi()));                   
       elif self.channel=="merged":
#        banner = TLatex(0.3,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu/e #nu "%(self.GetLumi()));                   
       banner.SetNDC(); banner.SetTextSize(0.04);
      else:
       if self.channel=="el":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e#nu "%(self.GetLumi()));                   
       elif self.channel=="mu":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu#nu "%(self.GetLumi()));                   
       elif self.channel=="merged":
#        banner = TLatex(0.22,0.96,("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu "%(self.GetLumi())));
        banner = TLatex(0.187919,0.960069,"CMS                       L = %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu/e #nu "%(self.GetLumi()));                   
       banner.SetNDC(); banner.SetTextSize(0.033);
                                                                                                         
      return banner;

    ### in order to make the legend
    def legend4Plot(self, plot, left=1, isFill=1, xoffset=0., yoffset=0., x_right_offset=0., y_upper_offset=0.):

        if left==-1:
            theLeg = TLegend(0.65+xoffset, 0.65+yoffset, 0.95+xoffset, 0.87+yoffset, "", "NDC");
            theLeg.SetName("theLegend");
            theLeg.SetBorderSize(0);
            theLeg.SetLineColor(0);
            theLeg.SetFillColor(0);
            if isFill:
                theLeg.SetFillStyle(1001);
            else:
                theLeg.SetFillStyle(0);
            theLeg.SetLineWidth(0);
            theLeg.SetLineStyle(0);
            theLeg.SetTextFont(42);
            theLeg.SetTextSize(.04);
        else:
            theLeg = TLegend(0.1+xoffset, 0.65+yoffset, 0.8+xoffset+x_right_offset, 0.92+yoffset+y_upper_offset, "", "NDC");
            theLeg.SetFillColor(0);
            theLeg.SetFillStyle(0);
            theLeg.SetNColumns(2);
            theLeg.SetTextSize(.033);

        theLeg.SetBorderSize(0);
        theLeg.SetLineColor(0);
        theLeg.SetLineWidth(0);
        theLeg.SetLineStyle(0);
        theLeg.SetTextFont(42);
        
        entryCnt = 0;
        objName_before = "";

        for obj in range(int(plot.numItems()) ):
            objName = plot.nameOf(obj);
            print objName;
            if not (((plot.getInvisible(objName)) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisi") or TString(objName).Contains("TLine") or objName==objName_before ):
                theObj = plot.getObject(obj);
                objTitle = objName;
                drawoption= plot.getDrawOptions(objName).Data()
                if drawoption=="P":drawoption="PE"
                if TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"):
                    theLeg.AddEntry(theObj, objName,"F");
                elif TString(objName).Contains("Graph") :
                    if not (objName_before=="Graph" or objName_before=="Uncertainty"): theLeg.AddEntry(theObj, "Uncertainty","F");
                else:
                    if TString(objName).Data() =="STop" : theLeg.AddEntry(theObj, "single Top","F");
                    elif TString(objName).Data()=="TTbar" : theLeg.AddEntry(theObj, "t#bar{t} powheg,pythia6","F");
                    elif TString(objName).Data()=="TTbar_herwig" : theLeg.AddEntry(theObj, "t#bar{t} mc@nlo,herwig++","F");
                    elif TString(objName).Data()=="VV" : theLeg.AddEntry(theObj, "WW/WZ/ZZ","F");
                    elif TString(objName).Data()=="WJets" : theLeg.AddEntry(theObj, "W+jets","F");
                    elif TString(objName).Contains("vbfH"): theLeg.AddEntry(theObj, (TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
                    else : theLeg.AddEntry(theObj, objTitle,drawoption);
                entryCnt=entryCnt+1;
            objName_before=objName;

        return theLeg;

    #### draw canvas with plots with pull
    def draw_canvas_with_pull(self, mplot, mplot_pull,parameters_list,in_directory, in_file_name, in_model_name="", show_constant_parameter=0, logy=0):# mplot + pull

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
        Directory = TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));        
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
            pad2.SetLogy() ;
            pad2.Update();
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());

        self.draw_canvas(mplot,in_directory,string_file_name.Data(),0,logy);

    #### jusr drawing canvas with no pull
    def draw_canvas(self, in_obj,in_directory, in_file_name, is_range=0, logy=0):

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

        banner = self.banner4Plot();
        banner.Draw();
        
        Directory=TString(in_directory+self.signal_sample+"_%02d_%02d/"%(options.cprime,options.BRnew));
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
            cMassFit.SetLogy() ;
            cMassFit.Update();
            rlt_file.ReplaceAll(".root","_log.root");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".root",".pdf");
            cMassFit.SaveAs(rlt_file.Data());
            rlt_file.ReplaceAll(".pdf",".png");
            cMassFit.SaveAs(rlt_file.Data());



    ### Get the liminosity value
    def GetLumi(self):

        if self.channel=="el": return 19.5
        if self.channel=="mu": return 19.5
        if self.channel=="merged": return 19.5


    #### defines two different way to fit depending on pythia or herwig analysis
    def fit_TTbar_controlsample(self, isherwig=0, ttbarMC=0):

     if ttbarMC ==0 :

      if isherwig ==0:

        print "fit_TTbar_controlsample --> Pythia samples"
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j");

        print "##################################################"
        print "############### Single Top DataSet ###############"
        print "##################################################"
          
        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop");
            

        if self.channel=="mu":
         self.fit_mj_single_MC(self.file_STop_mc,"_STop","ErfExpGaus_sp","_TTbar_controlsample"); ## Erf*exp + Gaus pass sample
         self.fit_mj_single_MC(self.file_STop_mc,"_STop_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp+Gaus fail sample
        else :
         self.fit_mj_single_MC(self.file_STop_mc,"_STop","ExpGaus","_TTbar_controlsample"); ## Exp + Gaus pass sample
         self.fit_mj_single_MC(self.file_STop_mc,"_STop_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp+Gaus fail sample
            
        self.fit_mj_single_MC(self.file_STop_mc,"_STop_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp for extreme fail

        ### Build WJet fit pass and fail distributions
        print "###########################################"
        print "############### WJets Pythia ##############"
        print "###########################################"


        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets0_mc,"_WJets0");

        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0","ErfExp","_TTbar_controlsample"); ## Erf*Exp pass sample
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp fail sample
        self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp extreme fail


        ### Build VV fit pass and fail distributions
        print "#########################################"
        print "############### VV Pythia ###############"
        print "#########################################"

        self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV");

        if self.channel == "mu" :
              self.fit_mj_single_MC(self.file_VV_mc,"_VV","ErfExpGaus_sp","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp +Gaus fail
        else:
              self.fit_mj_single_MC(self.file_VV_mc,"_VV","Gaus","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp +Gaus fail

        self.fit_mj_single_MC(self.file_VV_mc,"_VV_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp


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
        self.ScaleFactor_forPureWJet_TTbar_controlsample(self.file_data);

      else:


          print "###############################################"
          print "############# Single Top DataSet ##############"
          print "###############################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop_herwig");
          if self.channel=="mu":
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig","ErfExpGaus_sp","_TTbar_controlsample"); ## Erf*exp + Gaus pass sample
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp+Gaus fail sample
          else :
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig","ExpGaus","_TTbar_controlsample"); ## Exp + Gaus pass sample
           self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp+Gaus fail sample

          self.fit_mj_single_MC(self.file_STop_mc,"_STop_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp for extreme fail
          
          print "##########################################"
          print "############## WJets Herwig ##############"
          print "##########################################"

          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets1_mc,"_WJets0_herwig");

          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig","ErfExp","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig_failtau2tau1cut","Exp","_TTbar_controlsample");
          self.fit_mj_single_MC(self.file_WJets1_mc,"_WJets0_herwig_extremefailtau2tau1cut","Exp","_TTbar_controlsample");


          print "##########################################"
          print "################ VV Pythia ###############"
          print "##########################################"
          self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV_herwig");

          if self.channel == "mu" :
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig","ErfExpGaus_sp","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig_failtau2tau1cut","ExpGaus","_TTbar_controlsample"); ## Exp +Gaus fail
          else:
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig","Gaus","_TTbar_controlsample"); ## Exp+Gaus pass
              self.fit_mj_single_MC(self.file_VV_mc,"_VV_herwig_failtau2tau1cut","Exp","_TTbar_controlsample"); ## Exp +Gaus fail

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
          self.ScaleFactor_forPureWJet_TTbar_controlsample(self.file_data,"_herwig");

     else :

      if isherwig ==0:

            print "#####################################################"
            print "############# TTbar Powheg MC Analysis ##############"
            print "#####################################################"
            self.get_mj_and_mlvj_dataset_TTbar_controlsample_matched_with_truth(self.file_TTbar_mc,"_TTbar");

            print "##################### Print Matched Pass ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_matched_pass","2Gaus","_TTbar_controlsample");
            print "##################### Print Matched Fail ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_matched_fail","GausChebychev","_TTbar_controlsample");
            print "##################### Print UnMatched Pass ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_unmatched_pass","ErfExp","_TTbar_controlsample");
            print "##################### Print UnMatched Fail ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_unmatched_fail","ErfExp","_TTbar_controlsample");
            print "##################### Print Extreme Matched ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_extreme_matched","Exp","_TTbar_controlsample");
            print "##################### Print Extrme UnMatched ################################"
            self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar_extreme_unmatched","Exp","_TTbar_controlsample");

      else:

            print "#####################################################"
            print "############# TTbar Powheg MC Analysis ##############"
            print "#####################################################"
            self.get_mj_and_mlvj_dataset_TTbar_controlsample_matched_with_truth(self.file_TTbar_herwig,"_TTbar_herwig");

            print "##################### Print Herwig Matched Pass ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_herwig_matched_pass","2Gaus","_TTbar_controlsample");
            print "##################### Print Herwig Matched Fail ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_herwig_matched_fail","GausChebychev","_TTbar_controlsample");
            print "##################### Print Herwig UnMatched Pass ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_herwig_unmatched_pass","ErfExp","_TTbar_controlsample");
            print "##################### Print Herwig UnMatched Fail ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_herwig_unmatched_fail","ErfExp","_TTbar_controlsample");
            print "##################### Print Herwig Extreme Matched ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_extreme_matched","Exp","_TTbar_controlsample");
            print "##################### Print Herwig Extrme UnMatched ################################"
            self.fit_mj_single_MC(self.file_TTbar_herwig,"_TTbar_extreme_unmatched","Exp","_TTbar_controlsample");


class doFit_wj_and_wlvj_simultaneous:
    def __init__(self,isherwig=0):

        label = "";
        if isherwig==1 : label = "_herwig" ;

        self.workspace4fit_ = RooWorkspace("workspace4fit"+label+"_","workspace4fit"+label+"_");

        self.boostedW_fitter_el = doFit_wj_and_wlvj("el","ggH600",40,130,label, self.workspace4fit_)
        self.boostedW_fitter_mu = doFit_wj_and_wlvj("mu","ggH600",40,130,label, self.workspace4fit_)

        self.boostedW_fitter_el.fit_TTbar_controlsample(isherwig);
        self.boostedW_fitter_mu.fit_TTbar_controlsample(isherwig);

        self.workspace4fit_.data("rdataset_data"+label+"_mu_mj").Print();
        self.workspace4fit_.data("rdataset_data"+label+"_el_mj").Print();

        #### Define simultaneusly 4 category
        sample_type=RooCategory("sample_type"+label,"sample_type"+label);
        sample_type.defineType("mu_pass");
        sample_type.defineType("mu_fail");
        sample_type.defineType("el_pass");
        sample_type.defineType("el_fail");
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.) 


        rdataset_data_mu_mj      = self.workspace4fit_.data("rdataset_data"+label+"_mu_mj");
        rdataset_data_el_mj      = self.workspace4fit_.data("rdataset_data"+label+"_el_mj");
        rdataset_data_mu_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_data_el_mj_fail = self.workspace4fit_.data("rdataset_data"+label+"_failtau2tau1cut_el_mj"); 
 
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j");

        combData_data=RooDataSet("combData_data"+label,"combData_data"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_data_mu_mj),RooFit.Import("el_pass",rdataset_data_el_mj),RooFit.Import("mu_fail",rdataset_data_mu_mj_fail),RooFit.Import("el_fail",rdataset_data_el_mj_fail) );
        combData_data.Print();

        rdataset_TotalMC_mu_mj = self.workspace4fit_.data("rdataset_TotalMC"+label+"_mu_mj");
        rdataset_TotalMC_el_mj = self.workspace4fit_.data("rdataset_TotalMC"+label+"_el_mj");
        rdataset_TotalMC_mu_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_mu_mj"); 
        rdataset_TotalMC_el_mj_fail = self.workspace4fit_.data("rdataset_TotalMC"+label+"_failtau2tau1cut_el_mj"); 

        combData_TotalMC = RooDataSet("combData_TotalMC"+label,"combData_TotalMC"+label,RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("mu_pass",rdataset_TotalMC_mu_mj),RooFit.Import("el_pass",rdataset_TotalMC_el_mj),RooFit.Import("mu_fail",rdataset_TotalMC_mu_mj_fail),RooFit.Import("el_fail",rdataset_TotalMC_el_mj_fail) );
        combData_TotalMC.Print();

        # fit data
        model_data_mu      = self.workspace4fit_.pdf("model_data"+label+"_mu");
        model_data_fail_mu = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_mu");
        model_data_el      = self.workspace4fit_.pdf("model_data"+label+"_el");
        model_data_fail_el = self.workspace4fit_.pdf("model_data"+label+"_failtau2tau1cut_el");

        simPdf_data = RooSimultaneous("simPdf_data_em"+label,"simPdf_data_em"+label,sample_type);
        simPdf_data.addPdf(model_data_mu,"mu_pass");
        simPdf_data.addPdf(model_data_el,"el_pass");
        simPdf_data.addPdf(model_data_fail_mu,"mu_fail");
        simPdf_data.addPdf(model_data_fail_el,"el_fail");

        constrainslist_data_em=self.boostedW_fitter_el.constrainslist_data +self.boostedW_fitter_mu.constrainslist_data
        pdfconstrainslist_data_em=RooArgSet("pdfconstrainslist_data_em"+label);
        for i in range(len(constrainslist_data_em)):
            self.workspace4fit_.pdf(constrainslist_data_em[i]).Print();
            pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em[i]) );
        pdfconstrainslist_data_em.Print();

        rfresult_data=simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em) )
        rfresult_data=simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em) )

        # fit TotalMC
        model_TotalMC_mu      = self.workspace4fit_.pdf("model_TotalMC"+label+"_mu");
        model_TotalMC_fail_mu = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_mu");
        model_TotalMC_el      = self.workspace4fit_.pdf("model_TotalMC"+label+"_el");
        model_TotalMC_fail_el = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_el");

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em"+label,"simPdf_TotalMC_em"+label,sample_type);
        simPdf_TotalMC.addPdf(model_TotalMC_mu,"mu_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_el,"el_pass");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_mu,"mu_fail");
        simPdf_TotalMC.addPdf(model_TotalMC_fail_el,"el_fail");

        constrainslist_TotalMC_em=self.boostedW_fitter_el.constrainslist_TotalMC +self.boostedW_fitter_mu.constrainslist_TotalMC
        pdfconstrainslist_TotalMC_em=RooArgSet("pdfconstrainslist_TotalMC_em"+label);

        for i in range(len(constrainslist_TotalMC_em)):
            self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]).Print();
            pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) );
        pdfconstrainslist_TotalMC_em.Print();

        rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em) )
        rfresult_TotalMC=simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em) )

        self.boostedW_fitter_el.draw_ScaleFactor_forPureWJet_TTbar_controlsample(self.boostedW_fitter_el.file_data,label,0);
        self.boostedW_fitter_mu.draw_ScaleFactor_forPureWJet_TTbar_controlsample(self.boostedW_fitter_mu.file_data,label,0);


        rfresult_TotalMC.Print();
        rfresult_data.Print();
         
        rrv_eff_MC_el   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_el");
        rrv_eff_MC_mu   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_mu");
        rrv_mean_MC_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_el");
        rrv_sigma_MC_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_el");

        rrv_eff_data_el   = self.workspace4fit_.var("eff_ttbar_data"+label+"_el");
        rrv_eff_data_mu   = self.workspace4fit_.var("eff_ttbar_data"+label+"_mu");
        rrv_mean_data_el  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_el");
        rrv_sigma_data_el = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_el");

        rrv_eff_MC_el.Print()   ; rrv_eff_MC_mu.Print();
        rrv_eff_data_el.Print() ; rrv_eff_data_mu.Print();
        rrv_mean_MC_el.Print()  ; rrv_mean_data_el.Print();  
        rrv_sigma_MC_el.Print() ; rrv_sigma_data_el.Print();  

        pure_wtagger_sf_el=rrv_eff_data_el.getVal()/ rrv_eff_MC_el.getVal(); 
        pure_wtagger_sf_mu=rrv_eff_data_mu.getVal()/ rrv_eff_MC_mu.getVal(); 
        pure_wtagger_mean_shift_el= rrv_mean_data_el.getVal()-rrv_mean_MC_el.getVal();
        pure_wtagger_sigma_enlarge_el= rrv_sigma_data_el.getVal()/rrv_sigma_MC_el.getVal();

        pure_wtagger_sf_el_err= ( (rrv_eff_data_el.getError()/rrv_eff_data_el.getVal())**2 + (rrv_eff_MC_el.getError()/rrv_eff_MC_el.getVal())**2 )**0.5* pure_wtagger_sf_el

        pure_wtagger_sf_mu_err= ( (rrv_eff_data_mu.getError()/rrv_eff_data_mu.getVal())**2 + (rrv_eff_MC_mu.getError()/rrv_eff_MC_mu.getVal())**2 )**0.5* pure_wtagger_sf_mu
        
        pure_wtagger_mean_shift_err_el= ( rrv_mean_data_el.getError()**2 + rrv_mean_MC_el.getError()**2 )**0.5

        pure_wtagger_sigma_enlarge_err_el= ( (rrv_sigma_data_el.getError()/rrv_sigma_data_el.getVal())**2 + (rrv_sigma_MC_el.getError()/rrv_sigma_MC_el.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_el

        print "Pure W-tagger SF of el %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err)
        print "Pure W-tagger SF of mu %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err)
        print "Pure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el)
        print "Pure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el)

        self.boostedW_fitter_el.file_out_ttbar_control.write( "\n***************************************************" )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of el %s     : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el, pure_wtagger_sf_el_err) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger SF of mu %s        : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu, pure_wtagger_sf_mu_err) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger mean shift el %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_el, pure_wtagger_mean_shift_err_el) )
        self.boostedW_fitter_el.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge el %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_el, pure_wtagger_sigma_enlarge_err_el) )

        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_el_mj");
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_mu_mj");
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.Print();
        rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.Print();

        rrv_number_total_ttbar_TotalMC_el = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_el");
        rrv_number_total_ttbar_TotalMC_mu = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_mu");
        rrv_number_total_ttbar_data_el = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_el");
        rrv_number_total_ttbar_data_mu = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_mu");

        rrv_number_total_ttbar_TotalMC_el.Print();
        rrv_number_total_ttbar_TotalMC_mu.Print();
        rrv_number_total_ttbar_data_el.Print();
        rrv_number_total_ttbar_data_mu.Print();

        print "el TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal());
        print "mu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal());
        print "el data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal());
        print "mu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal());
        
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal()));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nel data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal()));

        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal()));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nmu data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal()));

        tmp_eff_MC_el_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal();
        tmp_eff_MC_mu_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_TotalMC_mu.getVal();
        tmp_eff_data_el_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal();
        tmp_eff_data_mu_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() / rrv_number_total_ttbar_data_mu.getVal();

        tmp_eff_MC_el_extremefail_error =tmp_eff_MC_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_el.getError()/rrv_number_total_ttbar_TotalMC_el.getVal() )**2 );
        tmp_eff_MC_mu_extremefail_error =tmp_eff_MC_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_mu.getError()/rrv_number_total_ttbar_TotalMC_mu.getVal() )**2 );
        tmp_eff_data_el_extremefail_error =tmp_eff_data_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_data_el.getError()/rrv_number_total_ttbar_data_el.getVal() )**2 );
        tmp_eff_data_mu_extremefail_error =tmp_eff_data_mu_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_mu_mj.getVal() )**2+ (rrv_number_total_ttbar_data_mu.getError()/rrv_number_total_ttbar_data_mu.getVal() )**2 );

        print "eff_MC_el_extremefail_error %s: %f"%(label,tmp_eff_MC_el_extremefail_error)
        print "eff_MC_mu_extremefail_error %s: %f"%(label,tmp_eff_MC_mu_extremefail_error)
        print "eff_data_el_extremefail_error %s: %f"%(label,tmp_eff_data_el_extremefail_error)
        print "eff_data_mu_extremefail_error %s: %f"%(label,tmp_eff_data_mu_extremefail_error)

        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_MC_el_extremefail_error %s: %f"%(label,tmp_eff_MC_el_extremefail_error) );
        self.boostedW_fitter_el.file_out_ttbar_control.write("\neff_data_el_extremefail_error %s: %f"%(label,tmp_eff_data_el_extremefail_error) );
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_MC_mu_extremefail_error %s: %f"%(label,tmp_eff_MC_mu_extremefail_error) );
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\neff_data_mu_extremefail_error %s: %f"%(label,tmp_eff_data_mu_extremefail_error) );
        
        rrv_eff_MC_el.Print()
        rrv_eff_MC_mu.Print()
        rrv_eff_data_el.Print()
        rrv_eff_data_mu.Print()

        tmp_eff_MC_el_LP =1. - rrv_eff_MC_el.getVal() - tmp_eff_MC_el_extremefail;
        tmp_eff_MC_mu_LP =1. - rrv_eff_MC_mu.getVal() - tmp_eff_MC_mu_extremefail;
        tmp_eff_data_el_LP =1. - rrv_eff_data_el.getVal() - tmp_eff_data_el_extremefail;
        tmp_eff_data_mu_LP =1. - rrv_eff_data_mu.getVal() - tmp_eff_data_mu_extremefail;

        tmp_eff_MC_el_LP_err = TMath.Sqrt( rrv_eff_MC_el.getError()**2 + tmp_eff_MC_el_extremefail_error**2 );
        tmp_eff_MC_mu_LP_err = TMath.Sqrt( rrv_eff_MC_mu.getError()**2 + tmp_eff_MC_mu_extremefail_error**2 );
        tmp_eff_data_el_LP_err = TMath.Sqrt( rrv_eff_data_el.getError()**2 + tmp_eff_data_el_extremefail_error**2 );
        tmp_eff_data_mu_LP_err = TMath.Sqrt( rrv_eff_data_mu.getError()**2 + tmp_eff_data_mu_extremefail_error**2 );

        print "LP Eff of el data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err);
        print "LP Eff of el MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err);
        print "LP Eff of mu data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err);
        print "LP Eff of mu MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err);

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el data %s: %f +/- %f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err));
        self.boostedW_fitter_el.file_out_ttbar_control.write("\nLP Eff of el MC %s: %f +/- %f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu data %s: %f +/- %f"%(label,tmp_eff_data_mu_LP, tmp_eff_data_mu_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nLP Eff of mu MC %s: %f +/- %f"%(label,tmp_eff_MC_mu_LP, tmp_eff_MC_mu_LP_err));
        
        pure_wtagger_sf_el_LP = tmp_eff_data_el_LP / tmp_eff_MC_el_LP;
        pure_wtagger_sf_mu_LP = tmp_eff_data_mu_LP / tmp_eff_MC_mu_LP;
        pure_wtagger_sf_el_LP_err = pure_wtagger_sf_el_LP*TMath.Sqrt( (tmp_eff_data_el_LP_err/tmp_eff_data_el_LP)**2 + (tmp_eff_MC_el_LP_err/tmp_eff_MC_el_LP)**2 );
        pure_wtagger_sf_mu_LP_err = pure_wtagger_sf_mu_LP*TMath.Sqrt( (tmp_eff_data_mu_LP_err/tmp_eff_data_mu_LP)**2 + (tmp_eff_MC_mu_LP_err/tmp_eff_MC_mu_LP)**2 );

        print "Pure W-tagger LP SF of el %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err)
        print "Pure W-tagger LP SF of mu %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err)

        self.boostedW_fitter_el.file_out_ttbar_control.write("\nPure W-tagger LP SF of el %s: %f +/- %f"%(label,pure_wtagger_sf_el_LP, pure_wtagger_sf_el_LP_err));
        self.boostedW_fitter_mu.file_out_ttbar_control.write("\nPure W-tagger LP SF of mu %s: %f +/- %f"%(label,pure_wtagger_sf_mu_LP, pure_wtagger_sf_mu_LP_err));



### function to call single channel fits
def control_sample(channel="mu",isherwig=0, ttbarMC=0):

    print "control sample "+channel;

    if ttbarMC ==0 :

     if isherwig ==0:
        ## create the object and do the single for pythia
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(isherwig,ttbarMC);

     elif isherwig ==1:
        ## create the object and do the single for herwig
        boostedW_fitter_herwig=doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(isherwig,ttbarMC);

     elif isherwig==2:
        ## do Pythia and herwig analysis at the same time
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(0,ttbarMC);
        boostedW_fitter_herwig = doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(1,ttbarMC);

    else:

       if isherwig ==0:
        ## create the object and do the single for pythia
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(isherwig,ttbarMC);
  
       elif isherwig ==1:
        ## create the object and do the single for herwig
        boostedW_fitter_herwig=doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(isherwig,ttbarMC);

       elif isherwig==2:
        ## do Pythia and herwig analysis at the same time
        boostedW_fitter = doFit_wj_and_wlvj(channel,"ggH600",40,130)
        boostedW_fitter.fit_TTbar_controlsample(0,ttbarMC);
        boostedW_fitter_herwig = doFit_wj_and_wlvj(channel,"ggH600",40,130,"_herwig")
        boostedW_fitter_herwig.fit_TTbar_controlsample(1,ttbarMC);
                                                             


### function to call simultaneous channel fits
def control_sample_simultaneous():

    print "control_sample_simultaneous";

    if options.herwig == 0 :
        boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous(options.herwig)
    elif options.herwig == 1 :
        boostedW_fitter_sim_herwig = doFit_wj_and_wlvj_simultaneous(options.herwig)
    elif options.herwig == 2:
        boostedW_fitter_sim=doFit_wj_and_wlvj_simultaneous()
        boostedW_fitter_sim=doFit_wj_and_wlvj_simultaneous(options.herwig)

### main code
if __name__ == '__main__':

    channel=options.channel;
 
    if options.fitwtaggersim:
        print 'fitwtagger for el+mu sample'
        control_sample_simultaneous();

    elif options.fitwtagger:
        print 'fitwtagger for %s sample'%(channel)
        control_sample(channel,options.herwig,options.ttbarMC);
 
        


