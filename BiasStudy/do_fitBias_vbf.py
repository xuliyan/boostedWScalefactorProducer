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
parser.add_option('--scalesignalwidth', help='reduce the signal width by a factor x', type=float, default=1.)
parser.add_option('--injectSingalStrenght', help='inject a singal in the toy generation', type=float, default=1.)

(options, args) = parser.parse_args()

ROOT.gSystem.Load(options.inPath+"/PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PDFs/MakePdf_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(options.inPath+"/PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(options.inPath+"/BiasStudy/BiasUtils_cxx.so")

from ROOT import draw_error_band, draw_error_band_extendPdf, draw_error_band_Decor, draw_error_band_shape_Decor, Calc_error_extendPdf, Calc_error
from ROOT import RooErfExpPdf, RooAlpha, RooAlpha4ErfPowPdf, RooAlpha4ErfPow2Pdf, RooAlpha4ErfPowExpPdf, PdfDiagonalizer, RooPowPdf, RooPow2Pdf, RooErfPowExpPdf, RooErfPowPdf, RooErfPow2Pdf, RooQCDPdf, RooUser1Pdf, RooBWRunPdf, RooAnaExpNPdf,RooExpNPdf, RooAlpha4ExpNPdf, RooExpTailPdf, RooPow3Pdf, RooErfPow3Pdf, RooUser1Pdf
from ROOT import biasModelAnalysis, MakeGeneralPdf, MakeExtendedModel, get_TTbar_mj_Model, get_STop_mj_Model, get_VV_mj_Model, get_WW_EWK_mj_Model, get_WJets_mj_Model, get_ggH_mj_Model, get_vbfH_mj_Model, get_TTbar_mlvj_Model, get_STop_mlvj_Model, get_VV_mlvj_Model, get_WW_EWK_mlvj_Model, get_WJets_mlvj_Model, get_ggH_mlvj_Model, get_vbfH_mlvj_Model, fix_Model,  clone_Model

from ROOT import setTDRStyle, get_pull, draw_canvas, draw_canvas_with_pull, legend4Plot, GetDataPoissonInterval, GetLumi

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

   	    
        ### Method for a single MC fit of the mj spectra giving: file name, label, model name
    def fit_mj_single_MC(self,in_file_name, label, in_model_name, additioninformation=""):

        print "############### Fit mj single MC sample",in_file_name," ",label,"  ",in_model_name," ##################"
        ## import variable and dataset
        rrv_mass_j = self.workspace4bias_.var("rrv_mass_j");
        rdataset_mj = self.workspace4bias_.data("rdataset4bias"+label+"_"+self.channel+"_mj");
        rdataset_mj.Print();
        
        ## make the extended model
        if additioninformation == 1:
         model = MakeExtendedModel(self.workspace4bias_,label+in_model_name,in_model_name,"_mj",self.channel,self.wtagger_label);
        else:
         model = MakeExtendedModel(self.workspace4bias_,label,in_model_name,"_mj",self.channel,self.wtagger_label);

        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.Extended(kTRUE),   RooFit.SumW2Error(kTRUE));
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult = model.fitTo(rdataset_mj,RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        getattr(self.workspace4bias_,"import")(model);
        
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

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),label+in_file_name,in_model_name,"em",0,1,GetLumi());


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
         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())
         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError()*self.workspace4bias_.var("rrv_scale_to_lumi"+label+"_"+self.channel).getVal())

         self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();
         
         if TString(label).Contains("ggH"):
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal());
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError());
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").Print();

         if TString(label).Contains("vbfH"):
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setVal(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getVal());
            self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").setError(self.workspace4bias_.var("rrv_number"+label+"_"+self.channel+"_mj").getError());
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
        constraintlist = RooArgList();

        ## make the extended pdf model
        model = MakeExtendedModel(self.workspace4bias_,label+in_range+mlvj_model,mlvj_model,"_mlvj",self.channel,self.wtagger_label,constraintlist,ismc);

        ## make the fit
        model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE) );
        rfresult = model.fitTo( rdataset, RooFit.Save(1), RooFit.SumW2Error(kTRUE) ,RooFit.Extended(kTRUE), RooFit.Minimizer("Minuit2") );
        rfresult.Print();

        ## set the name of the result of the fit and put it in the workspace
        rfresult.SetName("rfresult"+label+in_range+"_"+self.channel+"_mlvj")
        getattr(self.workspace4bias_,"import")(model);
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

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),label+in_file_name,mlvj_model,"em",0,1,GetLumi());

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
         model_WJets  = get_WJets_mj_Model(self.workspace4bias_,"_WJets0","",self.channel);
         model_STop   = get_STop_mj_Model(self.workspace4bias_,"_STop","",self.channel);
         model_VV     = get_VV_mj_Model(self.workspace4bias_,"_VV","",self.channel);
         model_WW_EWK = get_WW_EWK_mj_Model(self.workspace4bias_,"_WW_EWK","",self.channel);
	 model_TTbar  = get_TTbar_mj_Model(self.workspace4bias_,label,options.fgen,self.channel,0);
        else :
         model_TTbar  = get_TTbar_mj_Model(self.workspace4bias_,"_TTbar","",self.channel);
         model_STop   = get_STop_mj_Model(self.workspace4bias_,"_STop","",self.channel);
         model_VV     = get_VV_mj_Model(self.workspace4bias_,"_VV","",self.channel);
         model_WW_EWK = get_WW_EWK_mj_Model(self.workspace4bias_,"_WW_EWK","",self.channel);
	 model_WJets  = get_WJets_mj_Model(self.workspace4bias_,label,options.fgen,self.channel,0);


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

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),"m_j_sideband%s"%(label),"","em",0,1,GetLumi());

                
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
        model_VV_backgrounds     = get_VV_mlvj_Model(self.workspace4bias_,"_VV",mlvj_region,options.fgen,self.channel);
        model_STop_backgrounds   = get_STop_mlvj_Model(self.workspace4bias_,"_STop",mlvj_region,options.fgen,self.channel);
        model_WW_EWK_backgrounds = get_WW_EWK_mlvj_Model(self.workspace4bias_,"_WW_EWK",mlvj_region,options.fgen,self.channel);
     
        number_VV_sb_lo_mlvj      = self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        number_STop_sb_lo_mlvj    = self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        number_WW_EWK_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization       

        if options.ttbarcontrolregion == 0:
         model_TTbar_backgrounds  = get_TTbar_mlvj_Model(self.workspace4bias_,"_TTbar",mlvj_region,options.fgen,self.channel);
         number_TTbar_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
        else: 
         model_WJets_backgrounds  = get_WJets_mlvj_Model(self.workspace4bias_,"_WJets0",mlvj_region,options.fgen,self.channel);
         number_WJets_sb_lo_mlvj  = self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)); ## get the normalization
         

        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();
        self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(mlvj_region,options.fgen,self.channel)).Print();

        ### Make the Pdf for the WJets
        constraint = RooArgList(); 
        if options.ttbarcontrolregion == 0 :
            
         model_pdf_WJets = MakeGeneralPdf(self.workspace4bias_,"%s%s%s_from_fitting"%(label,mlvj_region,mlvj_model),mlvj_model,"_mlvj",self.wtagger_label,self.channel);
         model_pdf_WJets.Print();
         ### inititalize the value to what was fitted with the mc in the sideband
         number_WJets_sb_lo = self.workspace4bias_.var("rrv_number%s%s%s_%s_mlvj"%(label,mlvj_region,options.fgen,self.channel)).clone("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+self.channel+"_mlvj");
    
         model_WJets        = RooExtendPdf("model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model%s%s%s_from_fitting_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),model_pdf_WJets,number_WJets_sb_lo);
         number_WJets_sb_lo.Print();

         ## Add the other bkg component fixed to the total model --> in the extended way
         model_data = RooAddPdf("model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),"model_data%s%s%s_%s_mlvj"%(label,mlvj_region,mlvj_model,self.channel),RooArgList(model_WJets,model_VV_backgrounds, model_TTbar_backgrounds, model_STop_backgrounds, model_WW_EWK_backgrounds));

        else: 
         model_pdf_TTbar = MakeGeneralPdf(self.workspace4bias_,"%s%s%s_from_fitting"%(label,mlvj_region,mlvj_model),mlvj_model,"_mlvj",self.wtagger_label,self.channel);
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

        draw_canvas_with_pull(mplot,mplot_pull,RooArgList(parameters_list),"plots_%s_%s_%s_g1/mlvj_fitting_%s_%s/%s"%(options.additioninformation, self.channel, self.wtagger_label,options.fgen,options.fres,"basePlot/"),"m_lvj_sb_lo%s_%s"%(label,mlvj_model),"","em",0,1,GetLumi());

                
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
      else:
          self.fit_mlvj_model_single_MC(self.file_STop_mc,"_STop",options.mlvjregion,options.fgen,0,0,1);                  

      ######## get WW EWK and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WW_EWK_mc,"_WW_EWK","jet_mass_pr")# to get the shape of m_lvj                                                                             
      if fitjetmass:
          self.fit_mj_single_MC(self.file_WW_EWK_mc,"_WW_EWK","2Gaus"); 
      else:
          self.fit_mlvj_model_single_MC(self.file_WW_EWK_mc,"_WW_EWK",options.mlvjregion,options.fgen,0,0,1);
      ###### get WJets and fit it in the sb
      self.get_mj_and_mlvj_dataset(self.file_WJets0_mc,"_WJets0","jet_mass_pr")# to get the shape of m_lvj                                                                               
      if fitjetmass: 
	   self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0",options.fgen,1);
           if options.fgen != options.fres: self.fit_mj_single_MC(self.file_WJets0_mc,"_WJets0",options.fres,1);
      else: 
           self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fgen,0,0,1);
           if options.fgen != options.fres: self.fit_mlvj_model_single_MC(self.file_WJets0_mc,"_WJets0",options.mlvjregion,options.fres,0,0,1);
	
      ######## get TTbar and fit it
      self.get_mj_and_mlvj_dataset(self.file_TTbar_mc,"_TTbar")# to get the shape of m_lvj                                                                                               
      if fitjetmass: 
       if options.ttbarcontrolregion :    
          self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar",options.fgen,1);
          if options.fgen != options.fres: self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar",options.fres,1);
       else:
          self.fit_mj_single_MC(self.file_TTbar_mc,"_TTbar","2Gaus_ErfExp");	   		   
      else:
          self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fgen,0,0,1);
          if options.fgen != options.fres: self.fit_mlvj_model_single_MC(self.file_TTbar_mc,"_TTbar",options.mlvjregion,options.fres,0,0,1);

      ##### get data in sb and fit it                                                                                                                                              
      self.get_mj_and_mlvj_dataset(self.file_data,"_data", "jet_mass_pr"); ## global fit of data in the sidand fixing non dominant bkg             
      if fitjetmass :
         self.fit_WJetsNorm(label); ## fit jet mass distribution
      else:      
         self.fit_mlvj_in_Mj_sideband(label,options.mlvjregion,options.fgen,1); ## sideband or TTbar signal region fit
         self.fit_mlvj_in_Mj_sideband(label,options.mlvjregion,options.fres,1); ## sideband or TTbar signal region fit
     
     ##### fix signal and bkg models that are going to be used in the generation
     if fitjetmass :
         spectrum = "_mj"   ;
         signal_region = "";
         signal_model  = "2Gaus";
     else:
         spectrum = "_mlvj" ;
         signal_region = "_signal_region";
         signal_model  = "CB_v1";
    
     ###### fix the backgrund models for the generation
     print "#############################################################################################";
     print "################ Begin of the toy analysis -> fix Pdf to what is pre-fitted #################";
     print "#############################################################################################";
     
     fix_Model(self.workspace4bias_,"_%s"%self.ggH_sample,signal_region,spectrum,signal_model,self.channel);
     fix_Model(self.workspace4bias_,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,signal_model,self.channel);
     
     if options.isMC == 0  and options.ttbarcontrolregion == 0 and not options.fitjetmass :

      fix_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel);

     elif options.isMC == 0  and options.ttbarcontrolregion == 1 and not options.fitjetmass:

      fix_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fgen,self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,options.fres,self.channel);

     elif options.fitjetmass and not options.ttbarcontrolregion:

      fix_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,"",self.channel);

     elif options.fitjetmass and options.ttbarcontrolregion:

      fix_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_STop",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_VV",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,spectrum,"",self.channel);
      fix_Model(self.workspace4bias_,label,options.mlvjregion,spectrum,"",self.channel);
      
     print "#########################################################";
     print "################ Build the signal model #################";
     print "#########################################################";
     
     ### clone the signal shape --> parameter already fixed
     fitted_signal_ggH   = self.workspace4bias_.pdf("model_%s_%s_%s%s"%(self.ggH_sample,signal_region,self.channel,spectrum));
     fitted_signal_vbfH  = self.workspace4bias_.pdf("model_%s_%s_%s_%s"%(self.vbfhiggs_sample,signal_region,self.channel,spectrum));

     ### make the pdf for ggH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_ggH = RooArgList();
      
     if fitjetmass:
      model_signal_ggH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.ggH_sample,signal_region+"2Gaus"),"2Gaus",spectrum,self.wtagger_label,self.channel);  
      model_signal_ggH.Print();     
      clone_Model(self.workspace4bias_,model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,"2Gaus",self.channel);
      
      getattr(self.workspace4bias_,"import")(model_signal_ggH); 
    
      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"2Gaus_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"2Gaus_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);

     else:
      model_signal_ggH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.ggH_sample,signal_region+"CB_v1"),"CB_v1",spectrum,self.wtagger_label,self.channel);  
      model_signal_ggH.Print();     
 
      clone_Model(self.workspace4bias_,model_signal_ggH,"_%s"%self.ggH_sample,signal_region,spectrum,"CB_v1",self.channel);
      
      getattr(self.workspace4bias_,"import")(model_signal_ggH); 

      rrv_number_signal_signal_fit_ggH = RooRealVar("rrv_number_signal_region_fit_ggH","rrv_number_signal_region_fit_ggH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_ggH.setVal(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"CB_v1_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_ggH.setError(self.workspace4bias_.var("rrv_number_"+self.ggH_sample+signal_region+"CB_v1_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_ggH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_ggH);
     
      
     ### make the pdf for vbfH signal not in the extend way, cloning the parameter from the fitted value and then keep the absolute number of signal yields free
     constrainslist_signal_vbfH = RooArgList();
     if fitjetmass:
      model_signal_vbfH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.vbfhiggs_sample,signal_region+"2Gaus"),"2Gaus",spectrum,self.wtagger_label,self.channel);  
      model_signal_vbfH.Print();        

      clone_Model(self.workspace4bias_,model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,"2Gaus",self.channel);
      
      getattr(self.workspace4bias_,"import")(model_signal_vbfH); 

      rrv_number_signal_signal_fit_vbfH = RooRealVar("rrv_number_signal_region_fit_vbfH","rrv_number_signal_region_fit_vbfH",0,-1e8,1e8);
      rrv_number_signal_signal_fit_vbfH.setVal(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"2Gaus_"+self.channel+spectrum).getVal());
      rrv_number_signal_signal_fit_vbfH.setError(self.workspace4bias_.var("rrv_number_"+self.vbfhiggs_sample+signal_region+"2Gaus_"+self.channel+spectrum).getError());
      rrv_number_signal_signal_fit_vbfH.Print();
      getattr(self.workspace4bias_,"import")(rrv_number_signal_signal_fit_vbfH);

     else: 
      model_signal_vbfH = MakeGeneralPdf(self.workspace4bias_,"_%s%s_fit"%(self.vbfhiggs_sample,signal_region+"CB_v1"),"CB_v1",spectrum,self.wtagger_label,self.channel);  
      model_signal_vbfH.Print();        

      clone_Model(self.workspace4bias_,model_signal_vbfH,"_%s"%self.vbfhiggs_sample,signal_region,spectrum,"CB_v1",self.channel);
      
      getattr(self.workspace4bias_,"import")(model_signal_vbfH); 

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
          
      if(options.storeplot):
          mcWjetTreeResult.saveToysPlots(int(ratePlotsToStore),options.fitjetmass); 

      self.outputTree.Write();
      self.outputFile.Close();

     else: 
         
      ############### Make the Data analysis --> make the Entended pdf for the bkg
      constraintlist_bkg_data = RooArgList();
      print "#################################################################################";
      print "################ Start the MC analysis -> bkg model in the toy ##################";
      print "#################################################################################";

      ### take the models for the background component
      if fitjetmass: 

       model_VV_backgrounds     = get_VV_mj_Model(self.workspace4bias_,"_VV","",self.channel)
       model_STop_backgrounds   = get_STop_mj_Model(self.workspace4bias_,"_STop","",self.channel);
       model_WW_EWK_backgrounds = get_WW_EWK_mj_Model(self.workspace4bias_,"_WW_EWK","",self.channel);
       
       ## inflate the number of events and print them
       self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic);

       print "VV number ",self.workspace4bias_.var("rrv_number_VV_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "STop number ",self.workspace4bias_.var("rrv_number_STop_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.ttbarcontrolregion:
        model_TTbar_backgrounds  = get_TTbar_mj_Model(self.workspace4bias_,label,options.fgen,self.channel);
        model_WJets_backgrounds  = get_WJets_mj_Model(self.workspace4bias_,"_WJets0","",self.channel);

        self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);

        print "WJets number ",self.workspace4bias_.var("rrv_number_WJets0_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       else:

        model_TTbar_backgrounds  = get_TTbar_mj_Model(self.workspace4bias_,"_TTbar","",self.channel);           
        model_WJets_backgrounds  = get_WJets_mj_Model(self.workspace4bias_,label,options.fgen,"",self.channel);

        self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal()*options.inflatejobstatistic)
        self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()*options.inflatejobstatistic)

        print "WJets number ",self.workspace4bias_.var("rrv_number_TTbar_%s_mj"%(self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "TTbar number ",self.workspace4bias_.var("rrv_number%s%s_%s_mj"%(label,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

      else:	       
       ### in case of mWW analysis
       model_VV_backgrounds     = get_VV_mlvj_Model(self.workspace4bias_,"_VV",options.mlvjregion,options.fgen,self.channel);
       model_STop_backgrounds   = get_STop_mlvj_Model(self.workspace4bias_,"_STop",options.mlvjregion,options.fgen,self.channel);
       model_TTbar_backgrounds  = get_TTbar_mlvj_Model(self.workspace4bias_,"_TTbar",options.mlvjregion,options.fgen,self.channel);
       model_WW_EWK_backgrounds = get_WW_EWK_mlvj_Model(self.workspace4bias_,"_WW_EWK",options.mlvjregion,options.fgen,self.channel);
       model_WJets_backgrounds  = get_WJets_mlvj_Model(self.workspace4bias_,"_WJets0",options.mlvjregion,options.fgen,self.channel);

       ## inflate yields
       self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);
       self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic) ## get the normalization

       print "VV number ",self.workspace4bias_.var("rrv_number_VV%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "STop number ",self.workspace4bias_.var("rrv_number_STop%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
       print "WW_EWK number ",self.workspace4bias_.var("rrv_number_WW_EWK%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;

       if options.ttbarcontrolregion == 0:

        self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic)  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print "TTbar number ",self.workspace4bias_.var("rrv_number_TTbar%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print "WJets number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
       else:

        self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).setVal(self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()*options.inflatejobstatistic);  ## get the normalization
        self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").setVal(self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()*options.inflatejobstatistic);

        print " WJets number ",self.workspace4bias_.var("rrv_number_WJets0%s%s_%s_mlvj"%(options.mlvjregion,options.fgen,self.channel)).getVal()," inflate ",options.inflatejobstatistic;
        print " TTbar number ",self.workspace4bias_.var("rrv_number"+label+options.mlvjregion+options.fgen+"_from_fitting_"+self.channel+"_mlvj").getVal()," inflate ",options.inflatejobstatistic;
           
      #### make the global model for the background  
      if options.fitjetmass:

       model_bkg_data = MakeExtendedModel(self.workspace4bias_,label+signal_region+"_fit",options.fres,spectrum,self.channel,self.wtagger_label,constraintlist_bkg_data,1);
       model_bkg_data.Print();

       if options.fgen == options.fres:
        clone_Model(self.workspace4bias_,model_bkg_data,options.mlvjregion,spectrum,options.gen,self.channel,0);
       
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
      
      numevents_data   = (self.workspace4bias_.data("rdataset_data"+options.mlvjregion+"_"+self.channel+spectrum).sumEntries()+(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght)*options.inflatejobstatistic;
      print "##### number of events generated ",numevents_data ;
      
      print "###########################################################";
      print "################ Call the toy class tool ##################";
      print "###########################################################";

      if options.fitjetmass :
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
       mcWjetTreeResult.setPdfInformation(options.mlvjregion,spectrum,self.channel,label);
       mcWjetTreeResult.setSignalInjection(model_total_signal,(rrv_number_signal_signal_fit_ggH.getVal()+rrv_number_signal_signal_fit_vbfH.getVal())*options.injectSingalStrenght,options.scalesignalwidth);
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
     fitBiasAnalysis.shapeParametrizationAnalysis("_WJets0",options.fitjetmass);
  elif options.shapetest == 0 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.biasAnalysis("_TTbar",options.fitjetmass);
     fitBiasAnalysis.outputFile.Close();                                  
  elif options.shapetest == 1 and options.ttbarcontrolregion == 1:
     fitBiasAnalysis.shapeParametrizationAnalysis("_TTbar",options.fitjetmass);                     
