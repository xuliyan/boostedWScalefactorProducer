# A class which takes histograms and plots them in a versatile way
# inputs are file names which can b "data" or "MC"
import ROOT
import os
import sys
import math

from array import array
from optparse import OptionParser

from ROOT import RooTrace
from ROOT import TTree
from ROOT import gROOT

from bsmReweighter import *

import sys
from DataFormats.FWLite import Events, Handle

ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptFit(0)

ROOT.TH1.SetDefaultSumw2()
ROOT.TH2.SetDefaultSumw2()

ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);


# FWLITE stuff
ROOT.gSystem.Load('libCondFormatsJetMETObjects')
ROOT.gSystem.Load('libFWCoreFWLite');
ROOT.gSystem.Load('libFWCoreUtilities');  

### ------------ h e l p e r s --------------------
def deltaphi(phi1,phi2):
    if math.fabs(phi1-phi2) > ROOT.TMath.Pi() :
       return (2*ROOT.TMath.Pi()-math.fabs(phi1-phi2)); 
    else: 
        return math.fabs(phi1-phi2);

def getListRMS(list):
    mean = sum(list)/float(len(list));
    return math.sqrt(sum((n-mean)*(n-mean) for n in list)/len(list));

def getListMean(list):
    return sum(list)/float(len(list));

### ----------- class implementation -----------

class sampleWrapperClass:

    ### ------------------------------------------------
    def __init__(self, label, file, channel, sampleEffLumi, lumi, treename, isData,outputfiledirectory):

        self.IsData_    = isData; ### flag when running on data
        self.FileName_  = file; ### input file name 
        self.File_      = ROOT.TFile(file); ## get the root file
        self.InputTree_ = self.File_.Get(treename); ## get the input tree
        
        self.SampleWeight_ = lumi/sampleEffLumi; ## take the lumi re-weight factor
        
        self.JetPrefix_ = "GroomedJet_CA8";
        self.Label_ = label;
        self.Channel_ = channel
        self.OFileName_ = outputfiledirectory+"trainingtrees_"+channel+"/ofile_"+label+".root";
    
        # initialization for doing BSM reweighting --> set Higgs properties when running on higgs samples
        self.SignalMass_ = -1;
      
        if TString(self.FileName_).Contains("HWW")  and TString(self.FileName_).Contains("600") :  self.SignalMass_ = 600;
        if TString(self.FileName_).Contains("HWW")  and TString(self.FileName_).Contains("700") :  self.SignalMass_ = 700;
        if TString(self.FileName_).Contains("HWW")  and TString(self.FileName_).Contains("800") :  self.SignalMass_ = 800;
        if TString(self.FileName_).Contains("HWW")  and TString(self.FileName_).Contains("900") :  self.SignalMass_ = 900;    
        if TString(self.FileName_).Contains("HWW")  and TString(self.FileName_).Contains("1000") : self.SignalMass_ = 1000;   

        self.FitSMSignal = False;
        self.FitSMSignal_mean = -1;
        self.FitSMSignal_gamma = -1;
        self.isVBF_ = False;
        if TString(self.FileName_).Contains("VBFHWW") : self.isVBF_ = True;

        ###### reweight file --> open the root file for cps re-weighting at posteriori and prepare the information for VBF interference reweight
        if self.SignalMass_ > 0:
            self.rwName  = "H_CPSrw_%03d"%(self.SignalMass_);
            self.rwF     = ROOT.TFile("CPSrw/"+self.rwName+".root");
            self.h_rwCPS = self.rwF.Get(self.rwName);
            self.x_rwCPS = self.h_rwCPS.GetXaxis();

            self.signal_parameter_1              = [] ;
            self.signal_interference_parameter_1 = [] ;
            self.signal_parameter_05             = [] ;
            self.signal_interference_parameter_05= [] ;
            self.signal_parameter_2              = [] ;
            self.signal_interference_parameter_2 = [] ;

            mass = ["500","650","800","1000"];
            
            for ipar in range(7):
                self.signal_parameter_1.append(ROOT.TGraph(4));
                self.signal_parameter_2.append(ROOT.TGraph(4));
                self.signal_parameter_05.append(ROOT.TGraph(4));

            for ipar in range(9):
                self.signal_interference_parameter_1.append(ROOT.TGraph(4));
                self.signal_interference_parameter_2.append(ROOT.TGraph(4));
                self.signal_interference_parameter_05.append(ROOT.TGraph(4));
               
            iPoint = 0 ; 

            for imass in range(len(mass)):
                    
                inputFile_1  = ROOT.TFile("InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".1.root");
                inputFile_1.cd();
                f_sig_1  =      inputFile_1.Get("func_mg_1").Clone("f_sig_"+mass[imass]+"_1");
                f_sig_int_1  =  inputFile_1.Get("func_ph_1").Clone("f_sig_int_"+mass[imass]+"_05");
                print " input file ","InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".1.root" ;                
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ;
                 if ipar == 0:
                  y = ROOT.Double(f_sig_1.GetParameter(ipar));
                  self.signal_parameter_1[ipar].SetPoint(iPoint,x,y);
                 else:
                  y = ROOT.Double(f_sig_1.GetParameter(ipar));                         
                  self.signal_parameter_1[ipar].SetPoint(iPoint,x,y);
            
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ; y = ROOT.Double(0.);
                 if ipar == 0:   
                  y = ROOT.Double(f_sig_int_1.GetParameter(ipar)) ;   
                  self.signal_interference_parameter_1[ipar].SetPoint(iPoint,x,y);  
                 else:
                  y = ROOT.Double(f_sig_int_1.GetParameter(ipar));                        
                  self.signal_interference_parameter_1[ipar].SetPoint(iPoint,x,y);  
                inputFile_1.Close();
                
                inputFile_05 = ROOT.TFile("InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".0.5.root");
                f_sig_05 = inputFile_05.Get("func_mg_1").Clone("f_sig_"+mass[imass]+"_05");
                f_sig_int_05 = inputFile_05.Get("func_ph_1").Clone("f_sig_int_"+mass[imass]+"_05");
                print " input file ","InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".05.root" ;                
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ;
                 if ipar == 0:
                  y = ROOT.Double(f_sig_05.GetParameter(ipar));
                  self.signal_parameter_05[ipar].SetPoint(iPoint,x,y);

                 else:
                  y = ROOT.Double(f_sig_05.GetParameter(ipar));                         
                  self.signal_parameter_05[ipar].SetPoint(iPoint,x,y);
            
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ; y = ROOT.Double(0.);
                 if ipar == 0:   
                  y = ROOT.Double(f_sig_int_05.GetParameter(ipar)) ;   
                  self.signal_interference_parameter_05[ipar].SetPoint(iPoint,x,y);  
                 else:
                  y = ROOT.Double(f_sig_int_05.GetParameter(ipar));                        
                  self.signal_interference_parameter_05[ipar].SetPoint(iPoint,x,y);  
                inputFile_05.Close();
                
                inputFile_2  = ROOT.TFile("InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".2.root");
                f_sig_2  = inputFile_2.Get("func_mg_1").Clone("f_sig_"+mass[imass]+"_2");            
                f_sig_int_2  = inputFile_2.Get("func_ph_1").Clone("f_sig_int_"+mass[imass]+"_2");
                print " input file ","InterferenceRW/corr_woLeftRise/results_interference."+mass[imass]+".2.root" ;                
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ;
                 if ipar == 0:
                  y = ROOT.Double(f_sig_2.GetParameter(ipar));
                  self.signal_parameter_2[ipar].SetPoint(iPoint,x,y);

                 else:
                  y = ROOT.Double(f_sig_2.GetParameter(ipar));                         
                  self.signal_parameter_2[ipar].SetPoint(iPoint,x,y);
            
                for ipar in range(7):
                 x = ROOT.Double(mass[imass]) ; y = ROOT.Double(0.);
                 if ipar == 0:   
                  y = ROOT.Double(f_sig_int_2.GetParameter(ipar)) ;   
                  self.signal_interference_parameter_2[ipar].SetPoint(iPoint,x,y);  
                 else:
                  y = ROOT.Double(f_sig_int_2.GetParameter(ipar));                        
                  self.signal_interference_parameter_2[ipar].SetPoint(iPoint,x,y);  
                inputFile_2.Close();
                
                iPoint = iPoint + 1;

        ############ ---------- Set up jet corrections on the fly of R >= 0.7 jets
        fDir = "JECs/"      

        jecUncStr       = ROOT.std.string(fDir + "GR_R_53_V10_Uncertainty_AK7PFchs.txt");
        self.jecUnc_    = ROOT.JetCorrectionUncertainty(jecUncStr);

        jecUncStrAK5    = ROOT.std.string(fDir + "GR_R_53_V10_Uncertainty_AK5PFchs.txt");
        self.jecUncAK5_ = ROOT.JetCorrectionUncertainty(jecUncStrAK5);

        self.jerUncStr_       = ROOT.std.string(fDir + "pfJetResolutionMCtoDataCorrLUT.root");
        self.jerHistName_     = ROOT.std.string("pfJetResolutionMCtoDataCorrLUT");
        self.inputJERCFile_   = ROOT.TFile(ROOT.TString(self.jerUncStr_).Data());
        self.histoJERC_       = self.inputJERCFile_.Get(ROOT.TString(self.jerHistName_).Data());

        self.jetResolutionMC_CA8_  = [0.051,0.051,0.061,0.062,0.072,0.1]
        self.jetResolutionMC_AK5_  = [0.053,0.054,0.067,0.072,0.082,0.1]
        

        self.Random_ = ROOT.TRandom3();
        self.Random_.SetSeed();
        
        if self.Channel_ == "mu":
         self.LeptonScaleUnc_ = 0.002 ;
         self.LeptonScaleUncHighPT_ = math.sqrt(0.002**2 + 0.05**2) ;
         self.LeptonResolutionUnc_ = 0.006 ;
        else:
         self.LeptonScaleUnc_ = 0.006 ;
         self.LeptonScaleUncHighPT_ = 0.006 ;
         self.LeptonResolutionUnc_ = 0.014 ;
            
            
    def getLabel(self):
        return self.Label_

    def getTrainingTreeName(self):
        return self.OFileName_
    
    def getSampleLabel(self):
        return self.Label_

    def crystalBallLowHigh( self, x, par):
        xx = x[0];
        mean = par[1];
        sigma = par[2];
        alpha = par[3];
        n = par[4];
        alpha2 = par[5];
        n2 = par[6];
        
        if (xx-mean)/sigma > math.fabs(alpha) :
            A = ROOT.TMath.Power(n/math.fabs(alpha), n) * ROOT.TMath.Exp(-0.5 * alpha*alpha)
            B = n/math.fabs(alpha) - math.fabs(alpha);
            return par[0] * A * math.pow(B + (xx-mean)/sigma, -1.*n);                                   

        elif (xx-mean)/sigma < -1.*math.fabs(alpha2):
                               
            A = ROOT.TMath.Power(n2/math.fabs(alpha2), n2) * ROOT.TMath.Exp(-0.5 * alpha2*alpha2);
            B = n2/math.fabs(alpha2) - math.fabs(alpha2);

            return par[0] * A * ROOT.TMath.Power(B - (xx-mean)/sigma, -1.*n2);

        else:
            return par[0] * ROOT.TMath.Exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
                                       

    ### main function used to create the final otree
    def createTrainingTree(self):
        
        print self.FileName_
        self.NTree_ = self.InputTree_.GetEntries();
        
        print "Turning off branches...", self.FileName_

        self.InputTree_.SetBranchStatus("vbf*Charged*",0); 
        self.InputTree_.SetBranchStatus("vbf*Neutral*",0); 
        self.InputTree_.SetBranchStatus("vbf*Neutral*",0); 
        self.InputTree_.SetBranchStatus("vbf*Photon*",0); 
        self.InputTree_.SetBranchStatus("vbf*Electron*",0); 
        self.InputTree_.SetBranchStatus("vbf*HF*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Charged*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Neutral*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Neutral*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*Electron*",0); 
        self.InputTree_.SetBranchStatus("JetPFCor*HF*",0); 
        self.InputTree_.SetBranchStatus("Hadronic*",0); 
        self.InputTree_.SetBranchStatus("boosted*ang*",0); 
        self.InputTree_.SetBranchStatus("fit_*",0); 
        self.InputTree_.SetBranchStatus("W_tb*",0); 
        self.InputTree_.SetBranchStatus("W_Parton*",0); 
        self.InputTree_.SetBranchStatus("JetGen*",0); 

        print "Initializing sample: ", self.FileName_
        print "Nentries = ", self.NTree_

        # fill histograms
        self.createBRDTree();
        self.File_.Close();


    ### produce weights for alternative models --> input are the generator mass value for a given event, c' and BRnew
    def GetInteferenceWeights(self, mass, Cprime, BRnew=0, rwCPS=1):
        
        massin = self.SignalMass_;
        if massin < 0: return;

        ## define the fit range min and max as a function of the higgs mass 
        massmin = {600:200,700:200,800:400,900:400,1000:400};
        massmax = {600:1200,700:1200,800:1500,900:1600,1000:1800};

        #### read in original file and get lineshape from fit --> just for the first event, the file is read again and the lineshape is fitted just one time, taking
        #### the generated value, max and min range --> mean and gamma for BWRunning fit is returned from the function     
        if not self.FitSMSignal: 
            fitSM = FitMassPoint(self.FileName_, massin, massmin[massin], massmax[massin],rwCPS);
            self.FitSMSignal_mean  = fitSM[0];
            self.FitSMSignal_gamma = fitSM[1];
            self.FitSMSignal = True;
    
        #### create a weight for the given event after switch is turned --> taking the fitted mean, Gamma and the factor realted to c' and BRnew
        weight_width = lineshapeWidthReweight(mass, self.FitSMSignal_mean, self.FitSMSignal_gamma, Cprime/(1.-BRnew), massmin[massin], massmax[massin]); ## new gamma
        weight_xs = Cprime*(1-BRnew);    ## new xs
        return (weight_width*weight_xs); ## the weight is the product of the two

    ### ------------------------------------------------
    def createBRDTree(self):

        #### Output file --> create
        fname = self.OFileName_;        
        self.OFile_ = ROOT.TFile(fname,"RECREATE");        
        self.otree = ROOT.TTree("otree","otree");

        ########## Initialize Variables 
        self.InitializeVariables();

        ########## Create Output branches
        self.createBranches();
        
        ###### loop preparation
        prefix = self.JetPrefix_;        
        NLoop  = min(self.NTree_,10e9);
        NLoopWeight = self.NTree_/NLoop;
        
        wSampleWeight = NLoopWeight*self.SampleWeight_;

        RooTrace.active(ROOT.kTRUE);
        RooTrace.mark();

        ################################################
        ########## Start Loop On the Evvents ###########
        ################################################
        
        for iEvent in range(NLoop/100) :
            if iEvent  %10000 ==  0: print "iEvent = ", iEvent

            self.InputTree_.GetEntry(iEvent);

            if self.Channel_ == 'mu': lepLabel = "muon";
            if self.Channel_ == 'el': lepLabel = "electron";

            ###############################################
            ########## Fill TTbar Control Region Info  ####
            ###############################################
            
            ttbarlike  = 0;
            index_ca8_in_oppoHemi = [];
            
            for iJet in range(6): ### loop on ca8 jet over pt > 200, take dR(jet,lep) > Pi/2 --> opposite hemisphere 
                if getattr( self.InputTree_, "GroomedJet_CA8_pt" )[iJet] > 200 and math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_eta" )[iJet])<2.4:
                    j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[iJet];
                    j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[iJet];
                    l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                    l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );
                    l_charge = getattr( self.InputTree_, "W_"+lepLabel+"_charge" );
                    dR_lj = math.sqrt( (l_eta - j_ca8_eta)**2 + deltaphi(l_phi,j_ca8_phi)**2 ); ## delta R between lepton and hadronic W candidate over 200 GeV
                    if dR_lj > ROOT.TMath.Pi()/2.: index_ca8_in_oppoHemi.append(iJet);             ## opposite hemishpere

            minMass = -1;
            theca8Index = -1;
            ## loop on the vector of opposite hemisphere and take the jet with the mass closer to the nominal W mass
            for iJet in range(len(index_ca8_in_oppoHemi)):
                curmass = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[index_ca8_in_oppoHemi[iJet]]; ## search for the jet with mass closer to the W
                if math.fabs(curmass-80.385) < math.fabs(minMass-80.385): 
                    minMass = curmass;
                    theca8Index = index_ca8_in_oppoHemi[iJet];

            ### count ak5 jets
            index_ak5_in_sameHemi = [];
            index_ak5_in_oppoHemi = [];
            index_ak5_in_sameHemi_vetoca8 = [];
            index_ak5_in_oppoHemi_vetoca8 = [];
            
            index_ak5_in_sameHemi_csvl = [];
            index_ak5_in_oppoHemi_csvl = [];
            index_ak5_in_sameHemi_vetoca8_csvl = [];
            index_ak5_in_oppoHemi_vetoca8_csvl = [];
            
            index_ak5_in_sameHemi_csvm = [];
            index_ak5_in_oppoHemi_csvm = [];
            index_ak5_in_sameHemi_vetoca8_csvm = [];
            index_ak5_in_oppoHemi_vetoca8_csvm = [];
            
            index_ak5_in_sameHemi_csvt = [];
            index_ak5_in_oppoHemi_csvt = [];
            index_ak5_in_sameHemi_vetoca8_csvt = [];
            index_ak5_in_oppoHemi_vetoca8_csvt = [];

            ## HT of the event --> pt scalar sum 
            ttb_ht = getattr( self.InputTree_, "W_"+lepLabel+"_pt" );
            ttb_ht += getattr( self.InputTree_, "event_met_pfmet" ); 

            dR_ca8_bjet = -999 ;
            dR_ca8_jet  = -999 ;
            
            if theca8Index >= 0: ## if a W candidate is found --> loop on ak5 jet in the JetPFCor collection -> central jets
                for iJet in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 30: ## loop on the ak5 with pt > 30 only central jets 
                        ttb_ht += getattr( self.InputTree_, "JetPFCor_Pt" )[iJet];
                        j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                        j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[iJet];
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );               
                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + deltaphi(j_ak5_phi,j_ca8_phi)**2 ); ## delta R W-jet jet 
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + deltaphi(l_phi,j_ak5_phi)**2 );         ## delta R jet lep
                        
                        if dR_lj < ROOT.TMath.Pi()/2. :  ## same hemisphere wrt to the lepton, no cleaning with W-jet --> same hemisphere -> btag counting
                            index_ak5_in_sameHemi.append( iJet );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_sameHemi_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_sameHemi_csvm.append(iJet);                        
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_sameHemi_csvt.append(iJet);
                                                                 
                        elif dR_lj > ROOT.TMath.Pi()/2. : ## opposite hemisphere wrt to the lepton, no cleaning with W-jet --> opposite hemisphere -> btag counting
                            index_ak5_in_oppoHemi.append( iJet );    
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_oppoHemi_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_oppoHemi_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_oppoHemi_csvt.append(iJet);
                            
                        if dR_lj > ROOT.TMath.Pi()/2. and dR_jj > 0.8: ### veto ca8 counter --> cleaning inside W-jet -> a W never go to b, but some btag fake rate is possible 
                            index_ak5_in_oppoHemi_vetoca8.append( iJet );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_oppoHemi_vetoca8_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_oppoHemi_vetoca8_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_oppoHemi_vetoca8_csvt.append(iJet);

                        elif dR_lj < ROOT.TMath.Pi()/2. and  dR_jj > 0.8:
                            index_ak5_in_sameHemi_vetoca8.append( iJet );
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_sameHemi_vetoca8_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_sameHemi_vetoca8_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_sameHemi_vetoca8_csvt.append(iJet);

                        ### other topological info dR between W-jet and closer jet or bjet outside its cone
                        if dR_jj > 0.8 : 
                            if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.679 and dR_jj < math.fabs(dR_ca8_bjet) :
                                dR_ca8_bjet = dR_jj ;
                                dR_ca8_jet  = dR_jj ;
                            elif getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] <= 0.679 and dR_jj < math.fabs(dR_ca8_jet) :
                                dR_ca8_jet = dR_jj ;

                for iJet in range(6): ## same loop on forward jets
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 30: ## loop on the ak5 with pt > 30 only central jets 
                        ttb_ht += getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet];

                        j_ca8_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                        j_ca8_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[iJet];
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );               

                        dR_jj = math.sqrt( (j_ak5_eta - j_ca8_eta)**2 + deltaphi(j_ak5_phi,j_ca8_phi)**2 ); ## delta R W-jet jet 
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + deltaphi(l_phi,j_ak5_phi)**2 );         ## delta R jet lep
                        
                        if dR_lj < ROOT.TMath.Pi()/2. :  ## same hemisphere wrt to lepton --> no cleaning with ak5 inside W-jet cone  
                            index_ak5_in_sameHemi.append( iJet );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_sameHemi_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_sameHemi_csvm.append(iJet);                        
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_sameHemi_csvt.append(iJet);                        
                        elif dR_lj > ROOT.TMath.Pi()/2. : ## opposite hemisphere wrt to lepton --> no cleaning with ak5 inside W-jet cone
                            index_ak5_in_oppoHemi.append( iJet );    
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_oppoHemi_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_oppoHemi_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_oppoHemi_csvt.append(iJet);

                        if dR_lj > ROOT.TMath.Pi()/2. and dR_jj > 0.8: ## same hemisphere wrt to lepton -->  cleaning with ak5 inside W-jet cone 
                            index_ak5_in_oppoHemi_vetoca8.append( iJet );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_oppoHemi_vetoca8_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_oppoHemi_vetoca8_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_oppoHemi_vetoca8_csvt.append(iJet);
                        elif dR_lj < ROOT.TMath.Pi()/2. and  dR_jj > 0.8:
                            index_ak5_in_sameHemi_vetoca8.append( iJet );
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.244: index_ak5_in_sameHemi_vetoca8_csvl.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.679: index_ak5_in_sameHemi_vetoca8_csvm.append(iJet);
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_in_sameHemi_vetoca8_csvt.append(iJet);

                        ### other topological info dR between W-jet and closer jet or bjet outside its cone
                        if dR_jj > 0.8: 
                            if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] > 0.679 and dR_jj < math.fabs(dR_ca8_bjet) :
                                dR_ca8_bjet = dR_jj ;
                                dR_ca8_jet  = dR_jj ;
                            elif getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] <= 0.679 and dR_jj < math.fabs(dR_ca8_jet) :
                                dR_ca8_jet = dR_jj ;


                ### number of jets in the same emisphere of the lepton 
                self.ttb_nak5_same_[0] = int(len(index_ak5_in_sameHemi));
                ### number of jets in the same emisphere of the lepton passing csvl                 
                self.ttb_nak5_same_csvl_[0] = int(len(index_ak5_in_sameHemi_csvl));
                ### number of jets in the same emisphere of the lepton passing csvm                                 
                self.ttb_nak5_same_csvm_[0] = int(len(index_ak5_in_sameHemi_csvm));
                ### number of jets in the same emisphere of the lepton passing csvt                                                 
                self.ttb_nak5_same_csvt_[0] = int(len(index_ak5_in_sameHemi_csvt));

                ### number of jets in the opposite emisphere of the lepton                 
                self.ttb_nak5_oppo_[0] = int(len(index_ak5_in_oppoHemi));
                ### number of jets in the opposite emisphere of the lepton passing csvl                
                self.ttb_nak5_oppo_csvl_[0] = int(len(index_ak5_in_oppoHemi_csvl));
                ### number of jets in the opposite emisphere of the lepton passing csvm                
                self.ttb_nak5_oppo_csvm_[0] = int(len(index_ak5_in_oppoHemi_csvm));
                ### number of jets in the opposite emisphere of the lepton passing csvt                
                self.ttb_nak5_oppo_csvt_[0] = int(len(index_ak5_in_oppoHemi_csvt));

                ### number of jets in the opposite emisphere of the lepton cleaned wrt to the hadronic candidate                                
                self.ttb_nak5_oppoveto_[0] = int(len(index_ak5_in_oppoHemi_vetoca8));
                self.ttb_nak5_oppoveto_csvl_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvl));
                self.ttb_nak5_oppoveto_csvm_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvm));
                self.ttb_nak5_oppoveto_csvt_[0] = int(len(index_ak5_in_oppoHemi_vetoca8_csvt));

                ### number of jets in the same emisphere of the lepton cleaned wrt to the hadronic candidate                                
                self.ttb_nak5_sameveto_[0] = int(len(index_ak5_in_sameHemi_vetoca8));
                self.ttb_nak5_sameveto_csvl_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvl));
                self.ttb_nak5_sameveto_csvm_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvm));
                self.ttb_nak5_sameveto_csvt_[0] = int(len(index_ak5_in_sameHemi_vetoca8_csvt));

                ### jet quantities 
                self.ttb_ht_[0] = ttb_ht;
                self.ttb_ca8_mass_pr_[0]       = getattr( self.InputTree_, "GroomedJet_CA8_mass_pr" )[theca8Index];
                self.ttb_ca8_ungroomed_pt_[0]  = getattr( self.InputTree_, "GroomedJet_CA8_pt" )[theca8Index];
                self.ttb_ca8_ungroomed_eta_[0] = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];
                self.ttb_ca8_ungroomed_phi_[0] = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];
                self.ttb_ca8_ungroomed_e_[0]   = getattr( self.InputTree_, "GroomedJet_CA8_e" )[theca8Index];

                ### dR between W-jet and closer jet or bjet outside its cone
                self.ttb_dR_ca8_bjet_closer_[0] = dR_ca8_bjet ;
                self.ttb_dR_ca8_jet_closer_[0]  = dR_ca8_jet ;

                ### some generator info
                if not self.IsData_ :
                 self.ttb_ca8_ungroomed_gen_pt_[0]  = getattr( self.InputTree_, "GenGroomedJet_CA8_pt" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_eta_[0] = getattr( self.InputTree_, "GenGroomedJet_CA8_eta" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_phi_[0] = getattr( self.InputTree_, "GenGroomedJet_CA8_phi" )[theca8Index];
                 self.ttb_ca8_ungroomed_gen_e_[0]   = getattr( self.InputTree_, "GenGroomedJet_CA8_e" )[theca8Index];

                self.ttb_ca8_charge_[0]     = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge" )[theca8Index];
                self.ttb_ca8_charge_k05_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k05" )[theca8Index];
                self.ttb_ca8_charge_k07_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k07" )[theca8Index];
                self.ttb_ca8_charge_k10_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetcharge_k10" )[theca8Index];
                 
                self.ttb_ca8_tau2tau1_[0]       = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1" )[theca8Index];
                self.ttb_ca8_tau2tau1_exkT_[0]  = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1_exkT" )[theca8Index];
                self.ttb_ca8_tau2tau1_pr_[0]    = getattr( self.InputTree_, "GroomedJet_CA8_tau2tau1_pr" )[theca8Index];
                self.ttb_ca8_GeneralizedECF_[0] = getattr( self.InputTree_, "GroomedJet_CA8_jetGeneralizedECF" )[theca8Index];
                self.ttb_ca8_mu_[0]             = getattr( self.InputTree_, "GroomedJet_CA8_massdrop_pr" )[theca8Index];

                ### invariant mass of the full system
                ttb_ca8J_p4  = ROOT.TLorentzVector();        
                ttb_ca8J_pt  = getattr( self.InputTree_, "GroomedJet_CA8_pt" )[theca8Index];    
                ttb_ca8J_eta = getattr( self.InputTree_, "GroomedJet_CA8_eta" )[theca8Index];    
                ttb_ca8J_phi = getattr( self.InputTree_, "GroomedJet_CA8_phi" )[theca8Index];    
                ttb_ca8J_e   = getattr( self.InputTree_, "GroomedJet_CA8_e" )[theca8Index];                        

                ttb_ca8J_p4.SetPtEtaPhiE(ttb_ca8J_pt, ttb_ca8J_eta, ttb_ca8J_phi, ttb_ca8J_e)

                ttb_V_p4 = ROOT.TLorentzVector(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type0"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type0_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                ttb_V_p4.SetPxPyPzE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type2"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type2_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                ttb_V_p4.SetPxPyPzE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type0_met"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type0_met_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                ttb_V_p4.SetPxPyPzE(getattr(self.InputTree_,"W_px"),getattr(self.InputTree_,"W_py"),getattr(self.InputTree_,"W_pz_type2_met"),getattr(self.InputTree_,"W_e")); 
                self.ttb_ca8_mlvj_type2_met_[0] = (ttb_V_p4+ttb_ca8J_p4).M();

                self.ttb_ca8_px_[0] = ttb_ca8J_p4.Px();
                self.ttb_ca8_py_[0] = ttb_ca8J_p4.Py();
                self.ttb_ca8_pz_[0] = ttb_ca8J_p4.Pz();
                self.ttb_ca8_e_[0] =  ttb_ca8J_p4.E();

                ## preselection for the ttbar control region selection: at least one btag loose same or opposite hemisphere wrt to the lepton
                oppo1same1 = (self.ttb_nak5_same_csvl_[0] > 0  or self.ttb_nak5_oppo_csvl_[0] > 0) or (self.ttb_nak5_sameveto_csvl_[0]>0 or self.ttb_nak5_oppoveto_csvl_[0]>0);
                oppo2same0 = (self.ttb_nak5_same_csvl_[0] == 0 and self.ttb_nak5_oppo_csvl_[0] > 1) or (self.ttb_nak5_sameveto_csvl_[0]==0 or self.ttb_nak5_oppoveto_csvl_[0] > 1);

                if oppo1same1 or oppo2same0:
                    self.isttbar_[0] = 1 ;

            else: ## default value for non ttbar events

                self.ttb_nak5_same_[0] = -1;
                self.ttb_nak5_same_csvl_[0] = -1;
                self.ttb_nak5_same_csvm_[0] = -1;
                self.ttb_nak5_same_csvt_[0] = -1;

                self.ttb_nak5_oppo_[0] = -1;
                self.ttb_nak5_oppo_csvl_[0] = -1;
                self.ttb_nak5_oppo_csvm_[0] = -1;
                self.ttb_nak5_oppo_csvt_[0] = -1;
                self.ttb_nak5_oppoveto_[0] = -1;


                self.ttb_dR_ca8_bjet_closer_[0] = -999 ;
                self.ttb_dR_ca8_jet_closer_[0] = -999 ;

                self.ttb_nak5_oppoveto_csvl_[0] = -1;
                self.ttb_nak5_oppoveto_csvm_[0] = -1;
                self.ttb_nak5_oppoveto_csvt_[0] = -1;

                self.ttb_nak5_sameveto_csvl_[0] = -1;
                self.ttb_nak5_sameveto_csvm_[0] = -1;
                self.ttb_nak5_sameveto_csvt_[0] = -1;

                self.ttb_ca8_mass_pr_[0]       = -1;
                self.ttb_ca8_ungroomed_pt_[0]  = -1;
                self.ttb_ca8_ungroomed_eta_[0] = -999;
                self.ttb_ca8_ungroomed_phi_[0] = -999;
                self.ttb_ca8_ungroomed_e_[0]   = -1;
                
                self.ttb_ca8_ungroomed_gen_pt_[0]  = -1; 
                self.ttb_ca8_ungroomed_gen_eta_[0] = -999; 
                self.ttb_ca8_ungroomed_gen_phi_[0] = -999; 
                self.ttb_ca8_ungroomed_gen_e_[0]   = -1;

                self.ttb_ca8_charge_[0]     = -999; 
                self.ttb_ca8_charge_k05_[0] = -999; 
                self.ttb_ca8_charge_k07_[0] = -999; 
                self.ttb_ca8_charge_k10_[0] = -999; 
                 
                self.ttb_ca8_tau2tau1_[0]       = -999; 
                self.ttb_ca8_tau2tau1_exkT_[0]  = -999; 
                self.ttb_ca8_tau2tau1_pr_[0]    = -999;
                self.ttb_ca8_GeneralizedECF_[0] = -999;
                self.ttb_ca8_mu_[0]             = -999;

                self.ttb_ca8_mlvj_type0_[0] = -999;
                self.ttb_ca8_mlvj_type2_[0] = -999;

                self.ttb_ca8_mlvj_type0_met_[0] = -999;
                self.ttb_ca8_mlvj_type2_met_[0] = -999;

                self.ttb_ca8_px_[0] = -999; 
                self.ttb_ca8_py_[0] = -999;
                self.ttb_ca8_pz_[0] = -999;
                self.ttb_ca8_e_[0] =  -999;

            #########################################################################################################################                   
            ##################### Cuts used to define flag for final signal region and ttbar control region definition ##############
            #########################################################################################################################
                
            leptonCut = 30 ; 
            leptonCutString = "W_muon_pt";
            metCut = 40;
            if self.Channel_ == "el": 
                leptonCut = 35;    
                leptonCutString = "W_electron_pt";
                metCut    = 40;

            signallike = 0;
            if ( getattr( self.InputTree_, "W_pt" ) > 200 and
                 getattr( self.InputTree_, "GroomedJet_CA8_pt" )[0] > 200 and
                 math.fabs(getattr( self.InputTree_, "GroomedJet_CA8_eta" )[0]) < 2.4 and
                 getattr( self.InputTree_, "event_met_pfmet" ) > metCut   and
                 getattr( self.InputTree_, leptonCutString ) > leptonCut  and
                 getattr( self.InputTree_, "GroomedJet_CA8_deltaR_lca8jet") > 1.57) :

                if ( self.Channel_ == "mu" and
                     math.fabs(getattr( self.InputTree_, "W_muon_dz000")) < 0.02 and
                     math.fabs(getattr( self.InputTree_, "W_muon_dzPV"))  < 0.5  and
                     math.fabs( getattr( self.InputTree_, "W_muon_eta" )) < 2.1 ) :
                    signallike = 1 ; 
                elif self.Channel_ == "el" : signallike = 1 ; 

            ttbarlike = 0;
            if self.isttbar_[0] == 1 and getattr( self.InputTree_, "event_met_pfmet" ) > metCut and getattr( self.InputTree_, leptonCutString ) > leptonCut:
                ttbarlike = 1;

            ##################### store the info just for ttbar like or signal like events
            if ttbarlike == 1 or signallike == 1 :

             #####################################################################
             ############### Fill Event property only for interesting events #####
             #####################################################################
                
             self.event_[0]            = getattr( self.InputTree_, "event_evtNo");
             self.event_runNo_[0]      = getattr( self.InputTree_, "event_runNo" );
             self.event_lumi_[0]       = getattr( self.InputTree_, "event_lumi" );
             
             effwt = getattr( self.InputTree_, "effwt" );
             puwt  = getattr( self.InputTree_, "puwt" );
             
             totSampleWeight = 1.;
             if self.IsData_: totSampleWeight = wSampleWeight;
             else: totSampleWeight = wSampleWeight*effwt*puwt; ## total weight takes into account also pileUp and efficiency (not the btag)

             self.totalEventWeight_[0]  = totSampleWeight;
             self.eff_and_pu_Weight_[0] = effwt*puwt;
             self.wSampleWeight_[0]     = wSampleWeight;

             self.btag_weight_[0]    = getattr( self.InputTree_, "eff_btag" );
             self.btag_weight_up_[0] = getattr( self.InputTree_, "eff_btag_up" );
             self.btag_weight_dn_[0] = getattr( self.InputTree_, "eff_btag_dw" );
             self.btag_weight_up_dn_[0] = getattr( self.InputTree_, "eff_btag_up_dw" );
             self.btag_weight_dn_up_[0] = getattr( self.InputTree_, "eff_btag_dw_up" );

             self.nPV_[0]      = getattr( self.InputTree_, "event_nPV" );

             self.issignal_[0] = signallike;
             self.isttbar_[0]  = ttbarlike;

             self.numberJetBin_[0] = getattr(self.InputTree_, "numberJetBin")[0];
             self.numberJetBin2_[0] = getattr(self.InputTree_, "numberJetBin")[1];
             self.numberJetBin3_[0] = getattr(self.InputTree_, "numberJetBin")[2];
             self.numberJetBin4_[0] = getattr(self.InputTree_, "numberJetBin")[3];

             if not self.IsData_:
              self.numberJetBinGen_[0] = getattr(self.InputTree_, "numberJetBinGen")[0];
              self.numberJetBinGen2_[0] = getattr(self.InputTree_, "numberJetBinGen")[1];
              self.numberJetBinGen3_[0] = getattr(self.InputTree_, "numberJetBinGen")[2];
              self.numberJetBinGen4_[0] = getattr(self.InputTree_, "numberJetBinGen")[3];

             if self.Label_ == "TTbar_mcatnlo" :
              self.event_weight_[0] = getattr( self.InputTree_, "event_weight" )/math.fabs(getattr( self.InputTree_, "event_weight" )) ;
             else: self.event_weight_[0] = 1.;

             #### CPS part -> take value from histos external file
             rwCPS = 1;
             if self.SignalMass_ > 0:
              binVal = self.x_rwCPS.FindBin(getattr(self.InputTree_,"W_H_mass_gen"));
              if binVal > self.h_rwCPS.GetNbinsX(): binVal = self.h_rwCPS.GetNbinsX();
              if binVal < 1: binVal = 1;
              rwCPS = self.h_rwCPS.GetBinContent( binVal );

              ############## interference correction for vbf
              funz_sig_1  = ROOT.TF1("funz_sig_1",  self.crystalBallLowHigh, 200, 2000, 7);    
              funz_sAi_1  = ROOT.TF1("funz_sAi_1",  self.crystalBallLowHigh, 200, 2000, 9); 
              funz_sig_05 = ROOT.TF1("funz_sig_05", self.crystalBallLowHigh, 200, 2000, 7);    
              funz_sAi_05 = ROOT.TF1("funz_sAi_05", self.crystalBallLowHigh, 200, 2000, 9); 
              funz_sig_2  = ROOT.TF1("funz_sig_2",  self.crystalBallLowHigh, 200, 2000, 7);    
              funz_sAi_2  = ROOT.TF1("funz_sAi_2",  self.crystalBallLowHigh, 200, 2000, 9); 

              funz_sig_1.SetParameter(0,self.signal_parameter_1[0].Eval(self.SignalMass_));
              funz_sig_05.SetParameter(0,self.signal_parameter_05[0].Eval(self.SignalMass_));
              funz_sig_2.SetParameter(0,self.signal_parameter_2[0].Eval(self.SignalMass_));

              funz_sAi_1.SetParameter(0,self.signal_interference_parameter_1[0].Eval(self.SignalMass_));
              funz_sAi_05.SetParameter(0,self.signal_interference_parameter_05[0].Eval(self.SignalMass_));
              funz_sAi_2.SetParameter(0,self.signal_interference_parameter_2[0].Eval(self.SignalMass_));

              for ipar in range(6):
               funz_sig_1.SetParameter(ipar+1,self.signal_parameter_1[ipar+1].Eval(self.SignalMass_));
               funz_sig_05.SetParameter(ipar+1,self.signal_parameter_05[ipar+1].Eval(self.SignalMass_));
               funz_sig_2.SetParameter(ipar+1,self.signal_parameter_2[ipar+1].Eval(self.SignalMass_));
               
              for ipar in range(6):
               funz_sAi_1.SetParameter(ipar+1,self.signal_interference_parameter_1[ipar+1].Eval(self.SignalMass_));
               funz_sAi_05.SetParameter(ipar+1,self.signal_interference_parameter_05[ipar+1].Eval(self.SignalMass_));
               funz_sAi_2.SetParameter(ipar+1,self.signal_interference_parameter_2[ipar+1].Eval(self.SignalMass_));


              if funz_sig_1.Eval(getattr(self.InputTree_,"W_H_mass_gen")) != 0:
                 self.interferencevbf_1_[0]  = funz_sAi_1.Eval(getattr(self.InputTree_,"W_H_mass_gen"))/funz_sig_1.Eval(getattr(self.InputTree_,"W_H_mass_gen"));
              else:
                 self.interferencevbf_1_[0]  = 1 ;

              if funz_sig_05.Eval(getattr(self.InputTree_,"W_H_mass_gen")) != 0:
                 self.interferencevbf_05_[0]  = funz_sAi_05.Eval(getattr(self.InputTree_,"W_H_mass_gen"))/funz_sig_05.Eval(getattr(self.InputTree_,"W_H_mass_gen"));
              else:
                 self.interferencevbf_05_[0]  = 1 ;
                 
              if funz_sig_2.Eval(getattr(self.InputTree_,"W_H_mass_gen")) != 0:
                 self.interferencevbf_2_[0]  = funz_sAi_2.Eval(getattr(self.InputTree_,"W_H_mass_gen"))/funz_sig_2.Eval(getattr(self.InputTree_,"W_H_mass_gen"));
              else:
                 self.interferencevbf_2_[0]  = 1 ;
                 
             
             ############## interference weight and cps weight
             self.complexpolewtggH600    = getattr(self.InputTree_,"complexpolewtggH600")*rwCPS;
             self.interferencewtggH600   = getattr(self.InputTree_,"interferencewtggH600");
             self.avecomplexpolewtggH600 = getattr(self.InputTree_,"avecomplexpolewtggH600"); 
             self.interference_Weight_H600_[0] = self.complexpolewtggH600*self.interferencewtggH600/self.avecomplexpolewtggH600;  ## complete weight for standard higgs

             self.complexpolewtggH700    = getattr(self.InputTree_,"complexpolewtggH700")*rwCPS; 
             self.interferencewtggH700   = getattr(self.InputTree_,"interferencewtggH700");
             self.avecomplexpolewtggH700 = getattr(self.InputTree_,"avecomplexpolewtggH700"); 
             self.interference_Weight_H700_[0] = self.complexpolewtggH700*self.interferencewtggH700/self.avecomplexpolewtggH700;  ## complete weight for standard higgs

             self.complexpolewtggH800    = getattr(self.InputTree_,"complexpolewtggH800")*rwCPS; 
             self.interferencewtggH800   = getattr(self.InputTree_,"interferencewtggH800");
             self.avecomplexpolewtggH800 = getattr(self.InputTree_,"avecomplexpolewtggH800"); 
             self.interference_Weight_H800_[0] = self.complexpolewtggH800*self.interferencewtggH800/self.avecomplexpolewtggH800;  ## complete weight for standard higgs

             self.complexpolewtggH900    = getattr(self.InputTree_,"complexpolewtggH900")*rwCPS; 
             self.interferencewtggH900   = getattr(self.InputTree_,"interferencewtggH900");
             self.avecomplexpolewtggH900 = getattr(self.InputTree_,"avecomplexpolewtggH900"); 
             self.interference_Weight_H900_[0] = self.complexpolewtggH900*self.interferencewtggH900/self.avecomplexpolewtggH900;  ## complete weight for standard higgs

             self.complexpolewtggH1000    = getattr(self.InputTree_,"complexpolewtggH1000")*rwCPS; 
             self.interferencewtggH1000   = getattr(self.InputTree_,"interferencewtggH1000");
             self.avecomplexpolewtggH1000 = getattr(self.InputTree_,"avecomplexpolewtggH1000"); 
             self.interference_Weight_H1000_[0] = self.complexpolewtggH1000*self.interferencewtggH1000/self.avecomplexpolewtggH1000;  ## complete weight for standard higgs

             self.cps_Weight_H600_[0] = self.complexpolewtggH600/self.avecomplexpolewtggH600;
             self.cps_Weight_H700_[0] = self.complexpolewtggH700/self.avecomplexpolewtggH700;
             self.cps_Weight_H800_[0] = self.complexpolewtggH800/self.avecomplexpolewtggH800;
             self.cps_Weight_H900_[0] = self.complexpolewtggH900/self.avecomplexpolewtggH900;
             self.cps_Weight_H1000_[0] = self.complexpolewtggH1000/self.avecomplexpolewtggH1000;

             #### produce weights for alternative models
        
             if self.SignalMass_ > 0:

                curIntfRw = getattr(self.InputTree_,"interferencewtggH%03d"%(self.SignalMass_));  ## take the interference value filled in the tree for the right mass
                
                self.genHMass_[0] = getattr(self.InputTree_,"W_H_mass_gen");
                self.genHphi_[0]  = getattr(self.InputTree_,"W_H_phi_gen");
                self.genHeta_[0]  = getattr(self.InputTree_,"W_H_eta_gen");
                self.genHpt_[0]   = getattr(self.InputTree_,"W_H_pt_gen"); ## generator level higgs properties will be stored in the otree
               
                if self.isVBF_: ## if is vbf signal, also tag jet quark info stored 
                 self.genTagQuark1E_[0]    = getattr(self.InputTree_,"W_TagQuark_E")[0];
                 self.genTagQuark1eta_[0]  = getattr(self.InputTree_,"W_TagQuark_eta")[0];
                 self.genTagQuark1phi_[0]  = getattr(self.InputTree_,"W_TagQuark_phi")[0];
                 self.genTagQuark1pt_[0]   = getattr(self.InputTree_,"W_TagQuark_pt")[0];

                 self.genTagQuark2E_[0]    = getattr(self.InputTree_,"W_TagQuark_E")[1];
                 self.genTagQuark2eta_[0]  = getattr(self.InputTree_,"W_TagQuark_eta")[1];
                 self.genTagQuark2phi_[0]  = getattr(self.InputTree_,"W_TagQuark_phi")[1];
                 self.genTagQuark2pt_[0]   = getattr(self.InputTree_,"W_TagQuark_pt")[1];
                                                             

                for iPar in range(len(self.cprimeVals)): ## run over the possible value of the new couplig constant and BR
                 for jPar in range(len(self.brnewVals)): 
                  curCprime = float(self.cprimeVals[iPar])/10.; ## take the c' value
                  curBRnew  = float(self.brnewVals[jPar])/10.;   ## take the BRnew value
                  if self.isVBF_:
                    ### for VBF signal the bsmReweights is given just by the value returned by GetInteferenceWeights
                    self.bsmReweights[iPar][jPar][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew, rwCPS);  
                  else:
                    ### for ggH signal the bsmReweights is given by the value returned by GetInteferenceWeights * IntfRescale  
                    self.bsmReweights[iPar][jPar][0] = self.GetInteferenceWeights( getattr(self.InputTree_,"W_H_mass_gen"), curCprime, curBRnew, rwCPS)*IntfRescale(curIntfRw,curCprime,curBRnew);
             else:
                    ### default values for generator info 
                    self.genHMass_[0] = -1; 
                    self.genHphi_[0]  = -999;
                    self.genHeta_[0]  = -999;
                    self.genHpt_[0]   = -1;
                    self.genTagQuark1E_ = -1;
                    self.genTagQuark1eta_  = -999.;
                    self.genTagQuark1phi_  = -999.;
                    self.genTagQuark1pt_   = -1.;
                    self.genTagQuark2E_ = -1;
                    self.genTagQuark2eta_  = -999.;
                    self.genTagQuark2phi_  = -999.;
                    self.genTagQuark2pt_   = -1.;
                                                             
                    for iPar in range(len(self.cprimeVals)):  ## set to -1
                        for jPar in range(len(self.brnewVals)): 
                            self.bsmReweights[iPar][jPar][0] = -1;


             ################ lepton and met side                
             if self.Channel_ == "mu" :
                    self.l_pt_[0]     = getattr( self.InputTree_, "W_muon_pt" );
                    self.l_eta_[0]    = getattr( self.InputTree_, "W_muon_eta" );
                    self.l_phi_[0]    = getattr( self.InputTree_, "W_muon_phi" );
                    self.l_charge_[0] = getattr( self.InputTree_, "W_muon_charge" );
                    
             elif self.Channel_ == "el":
                    self.l_pt_[0]     = getattr( self.InputTree_, "W_electron_pt" );
                    self.l_eta_[0]    = getattr( self.InputTree_, "W_electron_eta" );
                    self.l_phi_[0]    = getattr( self.InputTree_, "W_electron_phi" );
                    self.l_charge_[0] = getattr( self.InputTree_, "W_electron_charge" );

             self.pfMET_[0]     = getattr( self.InputTree_, "event_met_pfmet" );        
             self.pfMET_Phi_[0] = getattr( self.InputTree_, "event_met_pfmetPhi" );        

               
             ################# Other basic Observables
             self.mass_lvj_type0_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type0" );
             self.mass_lvj_type2_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type2" );

             self.mass_lvj_type0_met_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type0_met" );
             self.mass_lvj_type2_met_[0] = getattr( self.InputTree_, "boostedW_lvj_m_type2_met" );

             self.mass_lv_subj_type0_[0] = getattr( self.InputTree_, "boosted_lvj_m_type0" );
             self.mass_lv_subj_type2_[0] = getattr( self.InputTree_, "boosted_lvj_m_type2" );

             self.mass_lv_subj_type0_met_[0] = getattr( self.InputTree_, "boosted_lvj_m_type0_met" );
             self.mass_lv_subj_type2_met_[0] = getattr( self.InputTree_, "boosted_lvj_m_type2_met" );

             self.v_pt_[0]   = getattr( self.InputTree_, "W_pt" );
             self.v_mt_[0]   = getattr( self.InputTree_, "W_mt" );
             self.v_eta_[0]  = getattr( self.InputTree_, "W_eta" );
             self.v_phi_[0]  = getattr( self.InputTree_, "W_phi" );

             self.nu_pz_type0_[0] = getattr( self.InputTree_, "W_nu1_pz_type0" );
             self.nu_pz_type2_[0] = getattr( self.InputTree_, "W_nu1_pz_type2" );

             self.nu_pz_type0_met_[0] = getattr( self.InputTree_, "W_nu1_pz_type0_met" );
             self.nu_pz_type2_met_[0] = getattr( self.InputTree_, "W_nu1_pz_type2_met" );

             self.W_pz_type0_[0] = getattr( self.InputTree_, "W_pz_type0" );
             self.W_pz_type2_[0] = getattr( self.InputTree_, "W_pz_type2" );

             self.W_pz_type0_met_[0] = getattr( self.InputTree_, "W_pz_type0_met" );
             self.W_pz_type2_met_[0] = getattr( self.InputTree_, "W_pz_type2_met" );

             W_Lepton_gen = ROOT.TLorentzVector() ;
             
             if self.IsData_ == False  :
                 self.nu_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" );                 
                 if self.Channel_ == "mu" :
                     W_Lepton_gen.SetPxPyPzE(getattr(self.InputTree_, "W_muon_px_gen" )+getattr(self.InputTree_, "W_neutrino_px_gen"),
                                             getattr(self.InputTree_, "W_muon_py_gen" )+getattr(self.InputTree_, "W_neutrino_py_gen" ),
                                             getattr(self.InputTree_, "W_muon_pz_gen" )+getattr(self.InputTree_, "W_neutrino_pz_gen" ),
                                             getattr(self.InputTree_, "W_muon_e_gen" )+getattr(self.InputTree_, "W_neutrino_e_gen" ));                     

                     self.W_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_muon_pz_gen" );
                     self.W_pt_gen_[0]  =  W_Lepton_gen.Pt();
                 
                 elif self.Channel_ == "el" :
                     W_Lepton_gen.SetPxPyPzE(getattr(self.InputTree_, "W_electron_px_gen" )+getattr(self.InputTree_, "W_neutrino_px_gen"),
                                             getattr(self.InputTree_, "W_electron_py_gen" )+getattr(self.InputTree_, "W_neutrino_py_gen" ),
                                             getattr(self.InputTree_, "W_electron_pz_gen" )+getattr(self.InputTree_, "W_neutrino_pz_gen" ),
                                             getattr(self.InputTree_, "W_electron_e_gen" )+getattr(self.InputTree_, "W_neutrino_e_gen" ));

                     self.W_pz_gen_[0]  =  getattr( self.InputTree_, "W_neutrino_pz_gen" ) + getattr( self.InputTree_, "W_electron_pz_gen" );
                     self.W_pt_gen_[0]  =  W_Lepton_gen.Pt();

                 
             ######## W-Jet stuff
             self.ungroomed_jet_pt_[0]  = getattr( self.InputTree_, prefix+"_pt" )[0];
             self.ungroomed_jet_eta_[0] = getattr( self.InputTree_, prefix+"_eta" )[0];
             self.ungroomed_jet_phi_[0] = getattr( self.InputTree_, prefix+"_phi" )[0];
             self.ungroomed_jet_e_[0]   = getattr( self.InputTree_, prefix+"_e" )[0];

             self.jet_mass_pr_[0] = getattr( self.InputTree_, prefix + "_mass_pr" )[0];
             self.jet_pt_pr_[0]   = getattr( self.InputTree_, prefix + "_pt_pr" )[0];

             self.jet_charge_[0]     = getattr( self.InputTree_, prefix + "_jetcharge" )[0];
             self.jet_charge_k05_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k05" )[0];
             self.jet_charge_k07_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k07" )[0];
             self.jet_charge_k10_[0] = getattr( self.InputTree_, prefix + "_jetcharge_k10" )[0];

             self.jet_grsens_ft_[0]   = getattr( self.InputTree_, prefix + "_mass_ft" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             self.jet_grsens_tr_[0]   = getattr( self.InputTree_, prefix + "_mass_tr" )[0] / getattr( self.InputTree_, prefix + "_mass" )[0];
             self.jet_massdrop_pr_[0] = getattr( self.InputTree_, prefix + "_massdrop_pr" )[0];    

             qjetmassdistribution = getattr( self.InputTree_, prefix+"_qjetmass" );
             qjetvol              = getListRMS(qjetmassdistribution)/getListMean(qjetmassdistribution);
             self.jet_qjetvol_[0] = qjetvol;

             self.jet_tau2tau1_[0]        = getattr( self.InputTree_, prefix + "_tau2tau1" )[0];     
             self.jet_tau2tau1_exkT_[0]   = getattr( self.InputTree_, prefix + "_tau2tau1_exkT" )[0];     
             self.jet_tau2tau1_pr_[0]     = getattr( self.InputTree_, prefix + "_tau2tau1_pr" )[0];     
             self.jet_GeneralizedECF_[0]  = getattr( self.InputTree_, prefix + "_jetGeneralizedECF" )[0];     
             self.jet_jetconstituents_[0] = getattr( self.InputTree_, prefix + "_jetconstituents" )[0];     

             self.jet_rcore4_[0] = getattr( self.InputTree_, prefix + "_rcores")[3*6 + 0];
             self.jet_rcore5_[0] = getattr( self.InputTree_, prefix + "_rcores")[4*6 + 0];
             self.jet_rcore6_[0] = getattr( self.InputTree_, prefix + "_rcores")[5*6 + 0];
             self.jet_rcore7_[0] = getattr( self.InputTree_, prefix + "_rcores")[6*6 + 0];

             self.jet_planarlow04_[0] = getattr( self.InputTree_, prefix + "_planarflow04");
             self.jet_planarlow05_[0] = getattr( self.InputTree_, prefix + "_planarflow05");
             self.jet_planarlow06_[0] = getattr( self.InputTree_, prefix + "_planarflow06");
             self.jet_planarlow07_[0] = getattr( self.InputTree_, prefix + "_planarflow07");

             self.pt1FracVal = max( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );
             self.pt2FracVal = min( getattr( self.InputTree_, prefix + "_prsubjet1ptoverjetpt" ), getattr( self.InputTree_, prefix + "_prsubjet2ptoverjetpt" ) );

             self.jet_pt1frac_[0] = self.pt1FracVal;
             self.jet_pt2frac_[0] = self.pt2FracVal;

             self.jet_sjdr_[0] = getattr( self.InputTree_, prefix + "_prsubjet1subjet2_deltaR" );       
                
             self.deltaR_lca8jet_[0]         = getattr( self.InputTree_, prefix + "_deltaR_lca8jet" );       
             self.deltaphi_METca8jet_[0]     = getattr( self.InputTree_, prefix + "_deltaphi_METca8jet_type2" );       
             self.deltaphi_Vca8jet_[0]       = getattr( self.InputTree_, prefix + "_deltaphi_Vca8jet_type2" );       
             self.deltaphi_METca8jet_met_[0] = getattr( self.InputTree_, prefix + "_deltaphi_METca8jet_type2_met" );       
             self.deltaphi_Vca8jet_met_[0]   = getattr( self.InputTree_, prefix + "_deltaphi_Vca8jet_type2_met" );       

             #### some generator information
             if not self.IsData_ :
              self.ungroomed_gen_jet_pt_[0]  = getattr( self.InputTree_, "Gen"+prefix+"_pt" )[0];
              self.ungroomed_gen_jet_eta_[0] = getattr( self.InputTree_, "Gen"+prefix+"_eta" )[0];
              self.ungroomed_gen_jet_phi_[0] = getattr( self.InputTree_, "Gen"+prefix+"_phi" )[0];
              self.ungroomed_gen_jet_e_[0]   = getattr( self.InputTree_, "Gen"+prefix+"_e" )[0];
 
              self.gen_jet_mass_pr_       = getattr( self.InputTree_, "Gen"+prefix + "_mass_pr" )[0];
              self.gen_jet_pt_pr_         = getattr( self.InputTree_, "Gen"+prefix + "_pt_pr" )[0];

#              self.gen_jet_grsens_ft_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_mass_ft" )[0] / getattr( self.InputTree_, "Gen"+prefix + "_mass" )[0];
#              self.gen_jet_grsens_tr_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_mass_tr" )[0] / getattr( self.InputTree_, "Gen"+prefix + "_mass" )[0];
#              self.gen_jet_massdrop_pr_[0] = getattr( self.InputTree_, "Gen"+prefix + "_massdrop_pr" )[0];    

#              qjetmassdistribution     = getattr( self.InputTree_, "Gen"+prefix+"_qjetmass" );
#              if qjetmassdistribution != 0 :
#               qjetvol                  = getListRMS(qjetmassdistribution)/getListMean(qjetmassdistribution);
#               self.gen_jet_qjetvol_[0] = qjetvol;

              self.gen_jet_tau2tau1_[0]        = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1" )[0];     
              self.gen_jet_tau2tau1_exkT_[0]   = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1_exkT" )[0];     
              self.gen_jet_tau2tau1_pr_[0]     = getattr( self.InputTree_, "Gen"+prefix + "_tau2tau1_pr" )[0];     
              self.gen_jet_jetconstituents_[0] = getattr( self.InputTree_, "Gen"+prefix + "_jetconstituents" )[0];     

#              self.gen_jet_rcore4_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[3*6 + 0];
#              self.gen_jet_rcore5_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[4*6 + 0];
#              self.gen_jet_rcore6_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[5*6 + 0];
#              self.gen_jet_rcore7_[0] = getattr( self.InputTree_, "Gen"+prefix + "_rcores")[6*6 + 0];
 
             ######################################################################################
             ### CA8 jet collection -> scale up and down for the jet energy scale uncertainty  ####
             ######################################################################################

             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0]);
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0]);                        
             self.j_jecfactor_up_[0] = self.jecUnc_.getUncertainty( True );

             self.jecUnc_.setJetEta( getattr( self.InputTree_, prefix+"_eta" )[0]);
             self.jecUnc_.setJetPt( getattr( self.InputTree_, prefix+"_pt" )[0]);               
             self.j_jecfactor_dn_[0] = self.jecUnc_.getUncertainty( False ) ;

             self.j_jecfactor_up_[0] = math.sqrt( self.j_jecfactor_up_[0]**2 + 0.02**2 ); ## inflate the uncertainty for the difference between ca8 and ak7
             self.j_jecfactor_dn_[0] = math.sqrt( self.j_jecfactor_dn_[0]**2 + 0.02**2 );

             self.curjes_up = 1 + self.j_jecfactor_up_[0];
             self.curjes_dn = 1 - self.j_jecfactor_dn_[0];

             jorig_pt  = getattr( self.InputTree_, prefix + "_pt_pr" )[0];    
             jorig_eta = getattr( self.InputTree_, prefix + "_eta_pr" )[0];    
             jorig_phi = getattr( self.InputTree_, prefix + "_phi_pr" )[0];    
             jorig_e   = getattr( self.InputTree_, prefix + "_e_pr" )[0];                        

             jdef_ptetaphie = ROOT.TLorentzVector();
             jdef_ptetaphie.SetPtEtaPhiE(jorig_pt, jorig_eta, jorig_phi, jorig_e); ## original jet 4V after pruning

             jdef_up = ROOT.TLorentzVector(jdef_ptetaphie.Px() * self.curjes_up, jdef_ptetaphie.Py() * self.curjes_up,
                                           jdef_ptetaphie.Pz() * self.curjes_up, jdef_ptetaphie.E()  * self.curjes_up); ## scaled up jet
             jdef_dn = ROOT.TLorentzVector(jdef_ptetaphie.Px() * self.curjes_dn, jdef_ptetaphie.Py() * self.curjes_dn,
                                           jdef_ptetaphie.Pz() * self.curjes_dn, jdef_ptetaphie.E()  * self.curjes_dn); ## scaled down jet


             self.jet_mass_pr_jes_up_[0] = jdef_up.M();  ## mass pruned up
             self.jet_mass_pr_jes_dn_[0] = jdef_dn.M();  ## mass pruned dwon

             jorig_pt  = getattr( self.InputTree_, prefix + "_pt" )[0];    
             jorig_eta = getattr( self.InputTree_, prefix + "_eta" )[0];    
             jorig_phi = getattr( self.InputTree_, prefix + "_phi" )[0];    
             jorig_e   = getattr( self.InputTree_, prefix + "_e" )[0];                        

             jdef_ptetaphie.SetPtEtaPhiE(jorig_pt, jorig_eta, jorig_phi, jorig_e); ## original jet 4V after pruning

             jdef_up.SetPxPyPzE(jdef_ptetaphie.Px() * self.curjes_up, jdef_ptetaphie.Py() * self.curjes_up,
                                jdef_ptetaphie.Pz() * self.curjes_up, jdef_ptetaphie.E()  * self.curjes_up); ## scaled up jet
             jdef_dn.SetPxPyPzE(jdef_ptetaphie.Px() * self.curjes_dn, jdef_ptetaphie.Py() * self.curjes_dn,
                                jdef_ptetaphie.Pz() * self.curjes_dn, jdef_ptetaphie.E()  * self.curjes_dn); ## scaled down jet

             self.ungroomed_jet_pt_jes_up_[0] = jdef_up.Pt(); ## ungroomed pt up
             self.ungroomed_jet_pt_jes_dn_[0] = jdef_dn.Pt(); ## ungroomed pt down

             ############################################################################
             ### CA8 jet collection -> smear for the jet energy resolution uncertainty  #
             ############################################################################

             ## Match the reco selected W candidate with the gen ones:
             smearFactor     = 1.;
             smearError      = 0.;
             smearFactor_up  = 1.;
             smearFactor_dn  = 1.;
             
             x = math.fabs(jdef_ptetaphie.Eta());
             y = jdef_ptetaphie.Pt();
             if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                 y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor += self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

             smearedEnergy    = jdef_ptetaphie.E();
             smearedEnergy_up = jdef_ptetaphie.E();
             smearedEnergy_dn = jdef_ptetaphie.E();

             if(math.fabs(jdef_ptetaphie.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 0.5 and math.fabs(jdef_ptetaphie.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 1.0 and math.fabs(jdef_ptetaphie.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 1.5 and math.fabs(jdef_ptetaphie.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 2.0 and math.fabs(jdef_ptetaphie.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 2.5 and math.fabs(jdef_ptetaphie.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));

             jet_smeared_byJERC    = ROOT.TLorentzVector();
             jet_smeared_byJERC_up = ROOT.TLorentzVector();
             jet_smeared_byJERC_dn = ROOT.TLorentzVector();

             jet_smeared_byJERC    = jdef_ptetaphie*(smearedEnergy/jdef_ptetaphie.E());
             jet_smeared_byJERC_up = jdef_ptetaphie*(smearedEnergy_up/jdef_ptetaphie.E());
             jet_smeared_byJERC_dn = jdef_ptetaphie*(smearedEnergy_dn/jdef_ptetaphie.E());

             self.ungroomed_jet_pt_jer_[0]    = jet_smeared_byJERC.Pt();
             self.ungroomed_jet_pt_jer_up_[0] = jet_smeared_byJERC_up.Pt();
             self.ungroomed_jet_pt_jer_dn_[0] = jet_smeared_byJERC_dn.Pt();

             jorig_pt  = getattr( self.InputTree_, prefix + "_pt_pr" )[0];    
             jorig_eta = getattr( self.InputTree_, prefix + "_eta_pr" )[0];    
             jorig_phi = getattr( self.InputTree_, prefix + "_phi_pr" )[0];    
             jorig_e   = getattr( self.InputTree_, prefix + "_e_pr" )[0];                        

             jdef_ptetaphie.SetPtEtaPhiE(jorig_pt, jorig_eta, jorig_phi, jorig_e); ## original jet 4V after pruning

             if(math.fabs(jdef_ptetaphie.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 0.5 and math.fabs(jdef_ptetaphie.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 1.0 and math.fabs(jdef_ptetaphie.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 1.5 and math.fabs(jdef_ptetaphie.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 2.0 and math.fabs(jdef_ptetaphie.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
             elif(math.fabs(jdef_ptetaphie.Eta())>= 2.5 and math.fabs(jdef_ptetaphie.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jdef_ptetaphie.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_CA8_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));


             self.jet_mass_pr_jer_[0]    = (jdef_ptetaphie*(smearedEnergy/jdef_ptetaphie.E())).M();
             self.jet_mass_pr_jer_up_[0] = (jdef_ptetaphie*(smearedEnergy_up/jdef_ptetaphie.E())).M();
             self.jet_mass_pr_jer_dn_[0] = (jdef_ptetaphie*(smearedEnergy_dn/jdef_ptetaphie.E())).M();
                        
             #################################################################                  
             ######## VBF jet Stuff --> maxpt pair  ##########################
             #################################################################
             
             vbf_maxpt_j1    = ROOT.TLorentzVector(0,0,0,0);
             vbf_maxpt_j2    = ROOT.TLorentzVector(0,0,0,0);

             vbf_maxpt_j1_up = ROOT.TLorentzVector(0,0,0,0);
             vbf_maxpt_j2_up = ROOT.TLorentzVector(0,0,0,0);
             vbf_maxpt_j1_dn = ROOT.TLorentzVector(0,0,0,0);
             vbf_maxpt_j2_dn = ROOT.TLorentzVector(0,0,0,0);
             
             vbf_maxpt_j1_smeared_byJERC    = ROOT.TLorentzVector();
             vbf_maxpt_j1_smeared_byJERC_up = ROOT.TLorentzVector();
             vbf_maxpt_j1_smeared_byJERC_dn = ROOT.TLorentzVector();
             vbf_maxpt_j2_smeared_byJERC    = ROOT.TLorentzVector();
             vbf_maxpt_j2_smeared_byJERC_up = ROOT.TLorentzVector();
             vbf_maxpt_j2_smeared_byJERC_dn = ROOT.TLorentzVector();

             if self.numberJetBin_[0] ==0 : ### no extra hard emission --> fill j1 and j2 with default values
              self.vbf_maxpt_jj_m_[0]     = 0;                               
              self.vbf_maxpt_jj_pt_[0]    = 0;                               
              self.vbf_maxpt_jj_eta_[0]   = 0;                               
              self.vbf_maxpt_jj_phi_[0]   = 0;                               

              self.vbf_maxpt_j1_m_[0]     = 0;   self.vbf_maxpt_j2_m_[0]   = 0;                               
              self.vbf_maxpt_j1_pt_[0]    = 0;   self.vbf_maxpt_j2_pt_[0]  = 0;                               
              self.vbf_maxpt_j1_eta_[0]   = 0;   self.vbf_maxpt_j2_eta_[0] = 0;                               
              self.vbf_maxpt_j1_phi_[0]   = 0;   self.vbf_maxpt_j2_phi_[0] = 0;                               

              self.vbf_maxpt_j1_QGLikelihood_[0]    = 0;  self.vbf_maxpt_j2_QGLikelihood_[0]    = 0;                                
              self.vbf_maxpt_j1_isPileUpMedium_[0]  = 0;  self.vbf_maxpt_j2_isPileUpMedium_[0]  = 0;                             
              self.vbf_maxpt_j1_isPileUpTight_[0]   = 0;  self.vbf_maxpt_j2_isPileUpTight_[0]   = 0;                               

              self.vbf_maxpt_j1_bDiscriminatorCSV_[0] = 0; self.vbf_maxpt_j2_bDiscriminatorCSV_[0] = 0; 

              if self.numberJetBinGen_[0] ==0 and not self.IsData_:

               self.vbf_maxpt_jj_m_gen_[0]     = 0;                               
               self.vbf_maxpt_jj_pt_gen_[0]    = 0;                               
               self.vbf_maxpt_jj_eta_gen_[0]   = 0;                               
               self.vbf_maxpt_jj_phi_gen_[0]   = 0;                               

               self.vbf_maxpt_j1_m_gen_[0]     = 0; self.vbf_maxpt_j2_m_gen_[0]    = 0;                              
               self.vbf_maxpt_j1_pt_gen_[0]    = 0; self.vbf_maxpt_j2_pt_gen_[0]   = 0;                              
               self.vbf_maxpt_j1_eta_gen_[0]   = 0; self.vbf_maxpt_j2_eta_gen_[0]  = 0;                              
               self.vbf_maxpt_j1_phi_gen_[0]   = 0; self.vbf_maxpt_j2_phi_gen_[0]  = 0;                              

               self.vbf_maxpt_j1_bDiscriminatorCSV_gen_[0]  = 0; self.vbf_maxpt_j2_bDiscriminatorCSV_gen_[0]  = 0;

              self.vbf_maxpt_j1_m_jes_up_[0]   = 0;                                
              self.vbf_maxpt_j1_pt_jes_up_[0]  = 0;                              
              self.vbf_maxpt_j1_eta_jes_up_[0] = 0;                                
              self.vbf_maxpt_j1_phi_jes_up_[0] = 0;                                

              self.vbf_maxpt_j1_m_jes_dn_[0]   = 0;                                
              self.vbf_maxpt_j1_pt_jes_dn_[0]  = 0;                              
              self.vbf_maxpt_j1_eta_jes_dn_[0] = 0;                                
              self.vbf_maxpt_j1_phi_jes_dn_[0] = 0;                                

              self.vbf_maxpt_j2_m_jes_dn_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jes_dn_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jes_dn_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jes_dn_[0] = 0;                                

              self.vbf_maxpt_j2_m_jes_up_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jes_up_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jes_up_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jes_up_[0] = 0;                                

              self.vbf_maxpt_j1_m_jer_[0]   = 0;                                
              self.vbf_maxpt_j1_pt_jer_[0]  = 0;                              
              self.vbf_maxpt_j1_eta_jer_[0] = 0;                                
              self.vbf_maxpt_j1_phi_jer_[0] = 0;                                

              self.vbf_maxpt_j2_m_jer_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_[0] = 0;                                

              self.vbf_maxpt_j1_m_jer_up_[0]   = 0;                                
              self.vbf_maxpt_j1_pt_jer_up_[0]  = 0;                              
              self.vbf_maxpt_j1_eta_jer_up_[0] = 0;                                
              self.vbf_maxpt_j1_phi_jer_up_[0] = 0;                                

              self.vbf_maxpt_j1_m_jer_dn_[0]   = 0;                                
              self.vbf_maxpt_j1_pt_jer_dn_[0]  = 0;                              
              self.vbf_maxpt_j1_eta_jer_dn_[0] = 0;                                
              self.vbf_maxpt_j1_phi_jer_dn_[0] = 0;                                

              self.vbf_maxpt_j2_m_jer_dn_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_dn_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_dn_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_dn_[0] = 0;                                

              self.vbf_maxpt_j2_m_jer_up_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_up_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_up_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_up_[0] = 0;                                

             ### one jet bin category --> fill j1 with the full information
             elif self.numberJetBin_[0] ==1 :

              self.vbf_maxpt_jj_m_[0]     = 0;                               
              self.vbf_maxpt_jj_pt_[0]    = 0;                               
              self.vbf_maxpt_jj_eta_[0]   = 0;                               
              self.vbf_maxpt_jj_phi_[0]   = 0;                               

              self.vbf_maxpt_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m");                               
              self.vbf_maxpt_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt");                               
              self.vbf_maxpt_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta");                               
              self.vbf_maxpt_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi");                               
 
              self.vbf_maxpt_j1_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j1_QGLikelihood");                                
 
              self.vbf_maxpt_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpMedium");                               
              self.vbf_maxpt_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpTight");                               

              self.vbf_maxpt_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV");

              self.vbf_maxpt_j2_m_[0]     = 0;                              
              self.vbf_maxpt_j2_pt_[0]    = 0;                               
              self.vbf_maxpt_j2_eta_[0]   = 0;                              
              self.vbf_maxpt_j2_phi_[0]   = 0;                               
 
              self.vbf_maxpt_j2_QGLikelihood_[0] = 0;                                
 
              self.vbf_maxpt_j2_isPileUpMedium_[0]  = 0;                               
              self.vbf_maxpt_j2_isPileUpTight_[0]   = 0;                               

              self.vbf_maxpt_j2_bDiscriminatorCSV_[0]  = 0;


              if self.numberJetBinGen_[0] ==1 and not self.IsData_:

               self.vbf_maxpt_jj_m_gen_[0]     = 0;                               
               self.vbf_maxpt_jj_pt_gen_[0]    = 0;                               
               self.vbf_maxpt_jj_eta_gen_[0]   = 0;                               
               self.vbf_maxpt_jj_phi_gen_[0]   = 0;                               

               self.vbf_maxpt_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m_gen");                               
               self.vbf_maxpt_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt_gen");                               
               self.vbf_maxpt_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta_gen");                               
               self.vbf_maxpt_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi_gen");                               

               self.vbf_maxpt_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV_gen");

               self.vbf_maxpt_j2_m_gen_[0]     = 0;                               
               self.vbf_maxpt_j2_pt_gen_[0]    = 0;                               
               self.vbf_maxpt_j2_eta_gen_[0]   = 0;                               
               self.vbf_maxpt_j2_phi_gen_[0]   = 0;                               

               self.vbf_maxpt_j2_bDiscriminatorCSV_gen_[0]  = 0;

              ## scale the hardest jet up and down for JES systematic  

              j_jecfactorAK5_up_ = array('f',[0.]);
              j_jecfactorAK5_dn_ = array('f',[0.]);
              self.jecUncAK5_.setJetEta(getattr(self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr(self.InputTree_, "vbf_maxpt_j1_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta(getattr( self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt(getattr( self.InputTree_, "vbf_maxpt_j1_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              vbf_maxpt_j1.SetPtEtaPhiM(self.vbf_maxpt_j1_pt_[0],self.vbf_maxpt_j1_eta_[0],self.vbf_maxpt_j1_phi_[0],self.vbf_maxpt_j1_m_[0]); ## original 4V

              vbf_maxpt_j1_up.SetPxPyPzE(vbf_maxpt_j1.Px()*(1+j_jecfactorAK5_up_[0]), 
                                         vbf_maxpt_j1.Py()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j1.Pz()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j1.E()*(1+j_jecfactorAK5_up_[0])); ## scaled up
                                                   

              vbf_maxpt_j1_dn.SetPxPyPzE(vbf_maxpt_j1.Px()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.Py()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.Pz()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.E()*(1-j_jecfactorAK5_dn_[0])); ## scaled down
                                                    
              self.vbf_maxpt_j1_m_jes_up_[0]   = vbf_maxpt_j1_up.M();                                
              self.vbf_maxpt_j1_pt_jes_up_[0]  = vbf_maxpt_j1_up.Pt();                              
              self.vbf_maxpt_j1_eta_jes_up_[0] = vbf_maxpt_j1_up.Eta();                                
              self.vbf_maxpt_j1_phi_jes_up_[0] = vbf_maxpt_j1_up.Phi();                                

              self.vbf_maxpt_j1_m_jes_dn_[0]   = vbf_maxpt_j1_dn.M();                                
              self.vbf_maxpt_j1_pt_jes_dn_[0]  = vbf_maxpt_j1_dn.Pt();                              
              self.vbf_maxpt_j1_eta_jes_dn_[0] = vbf_maxpt_j1_dn.Eta();                                
              self.vbf_maxpt_j1_phi_jes_dn_[0] = vbf_maxpt_j1_dn.Phi();                                

              self.vbf_maxpt_j2_m_jes_up_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jes_up_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jes_up_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jes_up_[0] = 0;                                

              self.vbf_maxpt_j2_m_jes_dn_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jes_dn_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jes_dn_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jes_dn_[0] = 0;                                

              ## smeared jets
              smearFactor     = 1.;
              smearError      = 0.;
              smearFactor_up  = 1.;
              smearFactor_dn  = 1.;
             
              x = math.fabs(vbf_maxpt_j1.Eta());
              y = vbf_maxpt_j1.Pt();
              if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                  y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor +=  self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

              smearedEnergy    = vbf_maxpt_j1.E();
              smearedEnergy_up = vbf_maxpt_j1.E();
              smearedEnergy_dn = vbf_maxpt_j1.E();

              if(math.fabs(vbf_maxpt_j1.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 0.5 and math.fabs(vbf_maxpt_j1.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 1.0 and math.fabs(vbf_maxpt_j1.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 1.5 and math.fabs(vbf_maxpt_j1.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 2.0 and math.fabs(vbf_maxpt_j1.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 2.5 and math.fabs(vbf_maxpt_j1.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));

              vbf_j1_smeared_byJERC    = vbf_maxpt_j1*(smearedEnergy/vbf_maxpt_j1.E());
              vbf_j1_smeared_byJERC_up = vbf_maxpt_j1*(smearedEnergy_up/vbf_maxpt_j1.E());
              vbf_j1_smeared_byJERC_dn = vbf_maxpt_j1*(smearedEnergy_dn/vbf_maxpt_j1.E());

              self.vbf_maxpt_j1_m_jer_[0]   = vbf_j1_smeared_byJERC.M();                                
              self.vbf_maxpt_j1_pt_jer_[0]  = vbf_j1_smeared_byJERC.Pt();                              
              self.vbf_maxpt_j1_eta_jer_[0] = vbf_j1_smeared_byJERC.Eta();                                
              self.vbf_maxpt_j1_phi_jer_[0] = vbf_j1_smeared_byJERC.Phi();                                

              self.vbf_maxpt_j2_m_jer_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_[0] = 0;                                

              self.vbf_maxpt_j1_m_jer_up_[0]   = vbf_j1_smeared_byJERC_up.M();                                
              self.vbf_maxpt_j1_pt_jer_up_[0]  = vbf_j1_smeared_byJERC_up.Pt();                              
              self.vbf_maxpt_j1_eta_jer_up_[0] = vbf_j1_smeared_byJERC_up.Eta();                                
              self.vbf_maxpt_j1_phi_jer_up_[0] = vbf_j1_smeared_byJERC_up.Phi();                                

              self.vbf_maxpt_j1_m_jer_dn_[0]   = vbf_j1_smeared_byJERC_dn.M();                                
              self.vbf_maxpt_j1_pt_jer_dn_[0]  = vbf_j1_smeared_byJERC_dn.Pt();                              
              self.vbf_maxpt_j1_eta_jer_dn_[0] = vbf_j1_smeared_byJERC_dn.Eta();                                
              self.vbf_maxpt_j1_phi_jer_dn_[0] = vbf_j1_smeared_byJERC_dn.Phi();                                

              self.vbf_maxpt_j2_m_jer_dn_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_dn_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_dn_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_dn_[0] = 0;                                

              self.vbf_maxpt_j2_m_jer_up_[0]   = 0;                                
              self.vbf_maxpt_j2_pt_jer_up_[0]  = 0;                              
              self.vbf_maxpt_j2_eta_jer_up_[0] = 0;                                
              self.vbf_maxpt_j2_phi_jer_up_[0] = 0;                                

             elif self.numberJetBin_[0] >=2 : ## 2 jet bin or vbf category

              self.vbf_maxpt_jj_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_jj_m");                               
              self.vbf_maxpt_jj_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_jj_pt");                               
              self.vbf_maxpt_jj_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_eta");                               
              self.vbf_maxpt_jj_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_phi");                               

              self.vbf_maxpt_j1_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m");                               
              self.vbf_maxpt_j1_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt");                               
              self.vbf_maxpt_j1_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta");                               
              self.vbf_maxpt_j1_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi");                               
 
              self.vbf_maxpt_j2_m_[0]     = getattr(self.InputTree_, "vbf_maxpt_j2_m");                               
              self.vbf_maxpt_j2_pt_[0]    = getattr(self.InputTree_, "vbf_maxpt_j2_pt");                               
              self.vbf_maxpt_j2_eta_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_eta");                               
              self.vbf_maxpt_j2_phi_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_phi");                               

              self.vbf_maxpt_j1_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j1_QGLikelihood");                                
              self.vbf_maxpt_j2_QGLikelihood_[0] = getattr(self.InputTree_, "vbf_maxpt_j2_QGLikelihood");                               
 
              self.vbf_maxpt_j1_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpMedium");                               
              self.vbf_maxpt_j1_isPileUpTight_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_isPileUpTight");                               

              self.vbf_maxpt_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_isPileUpMedium");                               
              self.vbf_maxpt_j2_isPileUpMedium_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_isPileUpMedium");                               

              self.vbf_maxpt_j1_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV");
              self.vbf_maxpt_j2_bDiscriminatorCSV_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_bDiscriminatorCSV");                                

              if not self.IsData_ and self.numberJetBinGen_[0]>= 2:
                  
               self.vbf_maxpt_jj_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_jj_m_gen");                               
               self.vbf_maxpt_jj_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_jj_pt_gen");                               
               self.vbf_maxpt_jj_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_eta_gen");                               
               self.vbf_maxpt_jj_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_jj_phi_gen");                               

               self.vbf_maxpt_j1_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j1_m_gen");                               
               self.vbf_maxpt_j1_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j1_pt_gen");                               
               self.vbf_maxpt_j1_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_eta_gen");                               
               self.vbf_maxpt_j1_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j1_phi_gen");                               
 
               self.vbf_maxpt_j2_m_gen_[0]     = getattr(self.InputTree_, "vbf_maxpt_j2_m_gen");                               
               self.vbf_maxpt_j2_pt_gen_[0]    = getattr(self.InputTree_, "vbf_maxpt_j2_pt_gen");                               
               self.vbf_maxpt_j2_eta_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_eta_gen");                               
               self.vbf_maxpt_j2_phi_gen_[0]   = getattr(self.InputTree_, "vbf_maxpt_j2_phi_gen");                               
 
               self.vbf_maxpt_j1_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j1_bDiscriminatorCSV_gen");
               self.vbf_maxpt_j2_bDiscriminatorCSV_gen_[0]  = getattr(self.InputTree_, "vbf_maxpt_j2_bDiscriminatorCSV_gen");                                

              ### scaling jets
              j_jecfactorAK5_up_ = array('f',[0.]);
              j_jecfactorAK5_dn_ = array('f',[0.]);
              self.jecUncAK5_.setJetEta(getattr(self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt( getattr(self.InputTree_, "vbf_maxpt_j1_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );

              self.jecUncAK5_.setJetEta(getattr( self.InputTree_, "vbf_maxpt_j1_eta" ));
              self.jecUncAK5_.setJetPt(getattr( self.InputTree_, "vbf_maxpt_j1_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              vbf_maxpt_j1.SetPtEtaPhiM(self.vbf_maxpt_j1_pt_[0],self.vbf_maxpt_j1_eta_[0],self.vbf_maxpt_j1_phi_[0],self.vbf_maxpt_j1_m_[0]); ## original 4V
              vbf_maxpt_j1_up.SetPxPyPzE(vbf_maxpt_j1.Px()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j1.Py()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j1.Pz()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j1.E()*(1+j_jecfactorAK5_up_[0])); ## scaled up
                                                   

              self.vbf_maxpt_j1_m_jes_up_[0]   = vbf_maxpt_j1_up.M();                                
              self.vbf_maxpt_j1_pt_jes_up_[0]  = vbf_maxpt_j1_up.Pt();                              
              self.vbf_maxpt_j1_eta_jes_up_[0] = vbf_maxpt_j1_up.Eta();                                
              self.vbf_maxpt_j1_phi_jes_up_[0] = vbf_maxpt_j1_up.Phi();                                

              vbf_maxpt_j1_dn.SetPxPyPzE(vbf_maxpt_j1.Px()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.Py()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.Pz()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j1.E()*(1-j_jecfactorAK5_dn_[0])); ## scaled down 
                                                    

              self.vbf_maxpt_j1_m_jes_dn_[0]   = vbf_maxpt_j1_dn.M();                                
              self.vbf_maxpt_j1_pt_jes_dn_[0]  = vbf_maxpt_j1_dn.Pt();                              
              self.vbf_maxpt_j1_eta_jes_dn_[0] = vbf_maxpt_j1_dn.Eta();                                
              self.vbf_maxpt_j1_phi_jes_dn_[0] = vbf_maxpt_j1_dn.Phi();                                
 
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxpt_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxpt_j2_pt" ));                        
              j_jecfactorAK5_up_[0] = self.jecUncAK5_.getUncertainty( True );
  
              self.jecUncAK5_.setJetEta( getattr( self.InputTree_, "vbf_maxpt_j2_eta" ));
              self.jecUncAK5_.setJetPt( getattr( self.InputTree_, "vbf_maxpt_j2_pt" ));               
              j_jecfactorAK5_dn_[0] = self.jecUncAK5_.getUncertainty( False ) ;

              vbf_maxpt_j2.SetPtEtaPhiM(self.vbf_maxpt_j2_pt_[0],self.vbf_maxpt_j2_eta_[0],self.vbf_maxpt_j2_phi_[0],self.vbf_maxpt_j2_m_[0]); ## orginal 4V 
              vbf_maxpt_j2_up.SetPxPyPzE(vbf_maxpt_j2.Px()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j2.Py()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j2.Pz()*(1+j_jecfactorAK5_up_[0]),
                                         vbf_maxpt_j2.E()*(1+j_jecfactorAK5_up_[0])); ## scaled up
                                                   

              self.vbf_maxpt_j2_m_jes_up_[0]   = vbf_maxpt_j2_up.M();                                
              self.vbf_maxpt_j2_pt_jes_up_[0]  = vbf_maxpt_j2_up.Pt();                              
              self.vbf_maxpt_j2_eta_jes_up_[0] = vbf_maxpt_j2_up.Eta();                                
              self.vbf_maxpt_j2_phi_jes_up_[0] = vbf_maxpt_j2_up.Phi();                                

              vbf_maxpt_j2_dn.SetPxPyPzE(vbf_maxpt_j2.Px()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j2.Py()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j2.Pz()*(1-j_jecfactorAK5_dn_[0]),
                                         vbf_maxpt_j2.E()*(1-j_jecfactorAK5_dn_[0])); ## scaled down

              self.vbf_maxpt_j2_m_jes_dn_[0]   = vbf_maxpt_j2_dn.M();                                
              self.vbf_maxpt_j2_pt_jes_dn_[0]  = vbf_maxpt_j2_dn.Pt();                              
              self.vbf_maxpt_j2_eta_jes_dn_[0] = vbf_maxpt_j2_dn.Eta();                                
              self.vbf_maxpt_j2_phi_jes_dn_[0] = vbf_maxpt_j2_dn.Phi();                                


              ## smeared jets
              smearFactor     = 1.;
              smearError      = 0.;
              smearFactor_up  = 1.;
              smearFactor_dn  = 1.;
             
              x = math.fabs(vbf_maxpt_j1.Eta());
              y = vbf_maxpt_j1.Pt();
              if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                  y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor += self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

              smearedEnergy    = vbf_maxpt_j1.E();
              smearedEnergy_up = vbf_maxpt_j1.E();
              smearedEnergy_dn = vbf_maxpt_j1.E();

              if(math.fabs(vbf_maxpt_j1.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 0.5 and math.fabs(vbf_maxpt_j1.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 1.0 and math.fabs(vbf_maxpt_j1.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 1.5 and math.fabs(vbf_maxpt_j1.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 2.0 and math.fabs(vbf_maxpt_j1.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j1.Eta())>= 2.5 and math.fabs(vbf_maxpt_j1.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j1.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));


              vbf_j1_smeared_byJERC    = vbf_maxpt_j1*(smearedEnergy/vbf_maxpt_j1.E());
              vbf_j1_smeared_byJERC_up = vbf_maxpt_j1*(smearedEnergy_up/vbf_maxpt_j1.E());
              vbf_j1_smeared_byJERC_dn = vbf_maxpt_j1*(smearedEnergy_dn/vbf_maxpt_j1.E());

              self.vbf_maxpt_j1_m_jer_[0]   = vbf_j1_smeared_byJERC.M();                                
              self.vbf_maxpt_j1_pt_jer_[0]  = vbf_j1_smeared_byJERC.Pt();                              
              self.vbf_maxpt_j1_eta_jer_[0] = vbf_j1_smeared_byJERC.Eta();                                
              self.vbf_maxpt_j1_phi_jer_[0] = vbf_j1_smeared_byJERC.Phi();                                

              self.vbf_maxpt_j1_m_jer_up_[0]   = vbf_j1_smeared_byJERC_up.M();                                
              self.vbf_maxpt_j1_pt_jer_up_[0]  = vbf_j1_smeared_byJERC_up.Pt();                              
              self.vbf_maxpt_j1_eta_jer_up_[0] = vbf_j1_smeared_byJERC_up.Eta();                                
              self.vbf_maxpt_j1_phi_jer_up_[0] = vbf_j1_smeared_byJERC_up.Phi();                                

              self.vbf_maxpt_j1_m_jer_dn_[0]   = vbf_j1_smeared_byJERC_dn.M();                                
              self.vbf_maxpt_j1_pt_jer_dn_[0]  = vbf_j1_smeared_byJERC_dn.Pt();                              
              self.vbf_maxpt_j1_eta_jer_dn_[0] = vbf_j1_smeared_byJERC_dn.Eta();                                
              self.vbf_maxpt_j1_phi_jer_dn_[0] = vbf_j1_smeared_byJERC_dn.Phi();                                

              smearFactor     = 1.;
              smearError      = 0.;
              smearFactor_up  = 1.;
              smearFactor_dn  = 1.;

              x = math.fabs(vbf_maxpt_j2.Eta());
              y = vbf_maxpt_j2.Pt();
              if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                  y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor += self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

              smearedEnergy    = vbf_maxpt_j2.E();
              smearedEnergy_up = vbf_maxpt_j2.E();
              smearedEnergy_dn = vbf_maxpt_j2.E();

              if(math.fabs(vbf_maxpt_j2.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j2.Eta())>= 0.5 and math.fabs(vbf_maxpt_j2.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j2.Eta())>= 1.0 and math.fabs(vbf_maxpt_j2.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j2.Eta())>= 1.5 and math.fabs(vbf_maxpt_j2.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j2.Eta())>= 2.0 and math.fabs(vbf_maxpt_j2.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
              elif(math.fabs(vbf_maxpt_j2.Eta())>= 2.5 and math.fabs(vbf_maxpt_j2.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = vbf_maxpt_j2.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));


              vbf_j2_smeared_byJERC    = vbf_maxpt_j2*(smearedEnergy/vbf_maxpt_j2.E());
              vbf_j2_smeared_byJERC_up = vbf_maxpt_j2*(smearedEnergy_up/vbf_maxpt_j2.E());
              vbf_j2_smeared_byJERC_dn = vbf_maxpt_j2*(smearedEnergy_dn/vbf_maxpt_j2.E());

              self.vbf_maxpt_j2_m_jer_[0]   = vbf_j2_smeared_byJERC.M();                                
              self.vbf_maxpt_j2_pt_jer_[0]  = vbf_j2_smeared_byJERC.Pt();                              
              self.vbf_maxpt_j2_eta_jer_[0] = vbf_j2_smeared_byJERC.Eta();                                
              self.vbf_maxpt_j2_phi_jer_[0] = vbf_j2_smeared_byJERC.Phi();                                

              self.vbf_maxpt_j2_m_jer_up_[0]   = vbf_j2_smeared_byJERC_up.M();                                
              self.vbf_maxpt_j2_pt_jer_up_[0]  = vbf_j2_smeared_byJERC_up.Pt();                              
              self.vbf_maxpt_j2_eta_jer_up_[0] = vbf_j2_smeared_byJERC_up.Eta();                                
              self.vbf_maxpt_j2_phi_jer_up_[0] = vbf_j2_smeared_byJERC_up.Phi();                                

              self.vbf_maxpt_j2_m_jer_dn_[0]   = vbf_j2_smeared_byJERC_dn.M();                                
              self.vbf_maxpt_j2_pt_jer_dn_[0]  = vbf_j2_smeared_byJERC_dn.Pt();                              
              self.vbf_maxpt_j2_eta_jer_dn_[0] = vbf_j2_smeared_byJERC_dn.Eta();                                
              self.vbf_maxpt_j2_phi_jer_dn_[0] = vbf_j2_smeared_byJERC_dn.Phi();                                

             ##############################################################################################################################  
             ####### build met and lepton -> scale up and down the met according to the jet energy -> buld the final invariant mass  ######
             ##############################################################################################################################
              
             met_vector_type0_met_up = ROOT.TLorentzVector();
             met_vector_type0_met_dn = ROOT.TLorentzVector();
             met_vector_type0_met_up.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr( self.InputTree_, "W_nu1_pz_type0_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type0_met")*getattr( self.InputTree_,"W_nu1_pz_type0_met"))); ## original met 4V
             met_vector_type0_met_dn = met_vector_type0_met_up ;

             met_vector_type2_met_up = ROOT.TLorentzVector();
             met_vector_type2_met_dn = ROOT.TLorentzVector();
             met_vector_type2_met_up.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr( self.InputTree_, "W_nu1_pz_type2_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type2_met")*getattr( self.InputTree_,"W_nu1_pz_type2_met"))); ## original met 4V, different pZ solution
             met_vector_type2_met_dn = met_vector_type2_met_up ;

             for iJet in range(6): ## loop on all the ak5 jet central collection

               if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 0. : ## loop on the central collection
                jet_vector = ROOT.TLorentzVector();
                jet_vector.SetPtEtaPhiE(getattr( self.InputTree_,"JetPFCor_Pt")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_Eta")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_Phi")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_E")[iJet]); ## original 4V

                self.jecUncAK5_.setJetEta( getattr(self.InputTree_, "JetPFCor_Eta")[iJet]);
                self.jecUncAK5_.setJetPt(  getattr(self.InputTree_, "JetPFCor_Pt" )[iJet]);                        
                j_jecfactorAK5_up = self.jecUncAK5_.getUncertainty( True );

                self.jecUncAK5_.setJetEta( getattr(self.InputTree_, "JetPFCor_Eta")[iJet]);
                self.jecUncAK5_.setJetPt(  getattr(self.InputTree_, "JetPFCor_Pt")[iJet]);               
                j_jecfactorAK5_dn = self.jecUncAK5_.getUncertainty( False ) ;


                jet_vector_up = ROOT.TLorentzVector(jet_vector.Px() * (1+j_jecfactorAK5_up), jet_vector.Py() * (1+j_jecfactorAK5_up),
                                                    jet_vector.Pz() * (1+j_jecfactorAK5_up), jet_vector.E()  * (1+j_jecfactorAK5_up)); ## vector up

                jet_vector_dn = ROOT.TLorentzVector(jet_vector.Px() * (1-j_jecfactorAK5_dn), jet_vector.Py() * (1-j_jecfactorAK5_dn),
                                                    jet_vector.Pz() * (1-j_jecfactorAK5_dn), jet_vector.E()  * (1-j_jecfactorAK5_dn)); ## vector down

                met_vector_type0_met_up = met_vector_type0_met_up + (jet_vector - jet_vector_up)  ; ## new met 4V up
                met_vector_type0_met_dn = met_vector_type0_met_dn + (jet_vector - jet_vector_dn)  ; ## new met 4V dn

                met_vector_type2_met_up = met_vector_type2_met_up + (jet_vector - jet_vector_up)  ;
                met_vector_type2_met_dn = met_vector_type2_met_dn + (jet_vector - jet_vector_dn)  ;

               if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 0. : ## loop on the forward collection
                   
                jet_vector = ROOT.TLorentzVector();
                jet_vector.SetPtEtaPhiE(getattr( self.InputTree_,"JetPFCorVBFTag_Pt")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_Eta")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_Phi")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_E")[iJet]); ## original 4V 
  
                self.jecUncAK5_.setJetEta( getattr(self.InputTree_, "JetPFCorVBFTag_Eta")[iJet]);
                self.jecUncAK5_.setJetPt(  getattr(self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet]);                        
                j_jecfactorAK5_up = self.jecUncAK5_.getUncertainty( True );

                self.jecUncAK5_.setJetEta( getattr(self.InputTree_, "JetPFCorVBFTag_Eta")[iJet]);
                self.jecUncAK5_.setJetPt(  getattr(self.InputTree_, "JetPFCorVBFTag_Pt")[iJet]);               
                j_jecfactorAK5_dn = self.jecUncAK5_.getUncertainty( False ) ;


                jet_vector_up = ROOT.TLorentzVector(jet_vector.Px() * (1+j_jecfactorAK5_up), jet_vector.Py() * (1+j_jecfactorAK5_up),
                                                    jet_vector.Pz() * (1+j_jecfactorAK5_up), jet_vector.E()  * (1+j_jecfactorAK5_up)); ## scaled up

                jet_vector_dn = ROOT.TLorentzVector(jet_vector.Px() * (1-j_jecfactorAK5_dn), jet_vector.Py() * (1-j_jecfactorAK5_dn),
                                                    jet_vector.Pz() * (1-j_jecfactorAK5_dn), jet_vector.E()  * (1-j_jecfactorAK5_dn)); ## scaled down

                met_vector_type0_met_up = met_vector_type0_met_up + (jet_vector - jet_vector_up)  ;
                met_vector_type0_met_dn = met_vector_type0_met_dn + (jet_vector - jet_vector_dn)  ;

                met_vector_type2_met_up = met_vector_type2_met_up + (jet_vector - jet_vector_up)  ;
                met_vector_type2_met_dn = met_vector_type2_met_dn + (jet_vector - jet_vector_dn)  ;

              
             self.pfMET_jes_up_[0]     = met_vector_type0_met_up.Pt();  ## scaled up met
             self.pfMET_Phi_jes_up_[0] = met_vector_type0_met_up.Phi(); 

             self.pfMET_jes_dn_[0]     = met_vector_type0_met_dn.Pt(); ## scaled down met
             self.pfMET_Phi_jes_dn_[0] = met_vector_type0_met_dn.Phi(); 

             ### final invariant mass --> shape and normalization sys 
             lepton_vector = ROOT.TLorentzVector();
             lepton_vector.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px"),
                                      getattr(self.InputTree_, "W_"+lepLabel+"_py"),
                                      getattr(self.InputTree_, "W_"+lepLabel+"_pz"),
                                      getattr(self.InputTree_, "W_"+lepLabel+"_e"));


             jet_pt  = getattr( self.InputTree_, prefix + "_pt" )[0];    
             jet_eta = getattr( self.InputTree_, prefix + "_eta" )[0];    
             jet_phi = getattr( self.InputTree_, prefix + "_phi" )[0];    
             jet_e   = getattr( self.InputTree_, prefix + "_e" )[0];                        

             jet_ptetaphie = ROOT.TLorentzVector();
             jet_ptetaphie.SetPtEtaPhiE(jet_pt, jet_eta, jet_phi, jet_e) ## original ungroomed W-jet 4V

             jet_up = ROOT.TLorentzVector(jet_ptetaphie.Px() * self.curjes_up, jet_ptetaphie.Py() * self.curjes_up,
                                          jet_ptetaphie.Pz() * self.curjes_up, jet_ptetaphie.E()  * self.curjes_up);
             jet_dn = ROOT.TLorentzVector(jet_ptetaphie.Px() * self.curjes_dn, jet_ptetaphie.Py() * self.curjes_dn,
                                          jet_ptetaphie.Pz() * self.curjes_dn, jet_ptetaphie.E()  * self.curjes_dn);


             self.mass_lvj_type0_met_jes_up_[0] = (lepton_vector + jet_up + met_vector_type0_met_up).M(); ## final invariant mass up
             self.mass_lvj_type0_met_jes_dn_[0] = (lepton_vector + jet_dn + met_vector_type0_met_dn).M(); ## final invariant mass dn
  
             self.mass_lvj_type2_met_jes_up_[0] = (lepton_vector + jet_up + met_vector_type2_met_up).M();
             self.mass_lvj_type2_met_jes_dn_[0] = (lepton_vector + jet_dn + met_vector_type2_met_dn).M();


             ##############################################################################################################################  
             ####### build met and lepton -> smeared the met according to the jet energy resolution-> buld the final invariant mass  ######
             ##############################################################################################################################
              
             met_vector_type0_met = ROOT.TLorentzVector();
             met_vector_type0_met.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr( self.InputTree_, "W_nu1_pz_type0_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type0_met")*getattr( self.InputTree_,"W_nu1_pz_type0_met"))); ## original met 4V
             met_vector_type0_met_up = met_vector_type0_met ;
             met_vector_type0_met_dn = met_vector_type0_met ;

             met_vector_type2_met = ROOT.TLorentzVector();
             met_vector_type2_met.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr( self.InputTree_, "W_nu1_pz_type2_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type2_met")*getattr( self.InputTree_,"W_nu1_pz_type2_met"))); ## original met 4V, different pZ solution
             met_vector_type2_met_up = met_vector_type2_met ;
             met_vector_type2_met_dn = met_vector_type2_met ;


             for iJet in range(6): ## loop on all the ak5 jet central collection

               if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 0. : ## loop on the central collection
                jet_vector = ROOT.TLorentzVector();
                jet_vector.SetPtEtaPhiE(getattr( self.InputTree_,"JetPFCor_Pt")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_Eta")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_Phi")[iJet],
                                        getattr( self.InputTree_,"JetPFCor_E")[iJet]); ## original 4V

                smearFactor     = 1.;
                smearError      = 0.;
                smearFactor_up  = 1.;
                smearFactor_dn  = 1.;

                x = math.fabs(jet_vector.Eta());
                y = jet_vector.Pt();
                if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                    y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor += self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

                smearedEnergy    = jet_vector.E();
                smearedEnergy_up = jet_vector.E();
                smearedEnergy_dn = jet_vector.E();

                if(math.fabs(jet_vector.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 0.5 and math.fabs(jet_vector.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 1.0 and math.fabs(jet_vector.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 1.5 and math.fabs(jet_vector.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 2.0 and math.fabs(jet_vector.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 2.5 and math.fabs(jet_vector.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));

                jet_smeared    = ROOT.TLorentzVector();
                jet_smeared_up = ROOT.TLorentzVector();
                jet_smeared_dn = ROOT.TLorentzVector();

                jet_smeared    = jet_vector*(smearedEnergy/jet_vector.E());
                jet_smeared_up = jet_vector*(smearedEnergy_up/jet_vector.E());
                jet_smeared_dn = jet_vector*(smearedEnergy_dn/jet_vector.E());

                met_vector_type0_met    = met_vector_type0_met    + (jet_vector - jet_smeared)  ; ## new met 4V up
                met_vector_type0_met_up = met_vector_type0_met_up + (jet_vector - jet_smeared_up)  ; ## new met 4V dn
                met_vector_type0_met_dn = met_vector_type0_met_dn + (jet_vector - jet_smeared_dn)  ; ## new met 4V dn

                met_vector_type2_met    = met_vector_type2_met    + (jet_vector - jet_smeared)  ;
                met_vector_type2_met_up = met_vector_type2_met_up + (jet_vector - jet_smeared_up)  ;
                met_vector_type2_met_dn = met_vector_type2_met_dn + (jet_vector - jet_smeared_dn)  ;

               if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 0. : ## loop on the forward collection
                   
                jet_vector = ROOT.TLorentzVector();
                jet_vector.SetPtEtaPhiE(getattr( self.InputTree_,"JetPFCorVBFTag_Pt")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_Eta")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_Phi")[iJet],
                                        getattr( self.InputTree_,"JetPFCorVBFTag_E")[iJet]); ## original 4V 

                jet_smeared    = ROOT.TLorentzVector();
                jet_smeared_up = ROOT.TLorentzVector();
                jet_smeared_dn = ROOT.TLorentzVector();

                smearFactor     = 1.;
                smearError      = 0.;
                smearFactor_up  = 1.;
                smearFactor_dn  = 1.;

                x = math.fabs(jet_vector.Eta());
                y = jet_vector.Pt();
                if( x > self.histoJERC_.GetXaxis().GetXmin() and x < self.histoJERC_.GetXaxis().GetXmax() and
                    y > self.histoJERC_.GetYaxis().GetXmin() and y < self.histoJERC_.GetYaxis().GetXmax()):
                  bin = self.histoJERC_.FindBin(x,y);
                  smearFactor += self.histoJERC_.GetBinContent(bin)-1.;
                  smearError  = self.histoJERC_.GetBinError(bin);
                  smearFactor_up = smearFactor+smearFactor*smearError ;
                  smearFactor_dn = smearFactor-smearFactor*smearError ;

                smearedEnergy    = jet_vector.E();
                smearedEnergy_up = jet_vector.E();
                smearedEnergy_dn = jet_vector.E();

                if(math.fabs(jet_vector.Eta())<0.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[0]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 0.5 and math.fabs(jet_vector.Eta())<1.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[1]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 1.0 and math.fabs(jet_vector.Eta())<1.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[2]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 1.5 and math.fabs(jet_vector.Eta())<2.0):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[3]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 2.0 and math.fabs(jet_vector.Eta())<2.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[4]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                elif(math.fabs(jet_vector.Eta())>= 2.5 and math.fabs(jet_vector.Eta())<4.5):
                  if(smearFactor > 1 ):
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor*smearFactor-1)));
                  else : 
                   smearedEnergy    = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor*smearFactor)));
                  if(smearFactor_up > 1 ): 
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_up*smearFactor_up-1)));
                  else :
                   smearedEnergy_up = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_up*smearFactor_up)));
                  if(smearFactor_dn > 1 ): 
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(smearFactor_dn*smearFactor_dn-1)));
                  else :
                   smearedEnergy_dn = jet_vector.E()*(1+self.Random_.Gaus(0.,self.jetResolutionMC_AK5_[5]*math.sqrt(1-smearFactor_dn*smearFactor_dn)));
                

                jet_smeared    = jet_vector*(smearedEnergy/jet_vector.E());
                jet_smeared_up = jet_vector*(smearedEnergy_up/jet_vector.E());
                jet_smeared_dn = jet_vector*(smearedEnergy_dn/jet_vector.E());

                met_vector_type0_met    = met_vector_type0_met    + (jet_vector - jet_smeared)  ; ## new met 4V up
                met_vector_type0_met_up = met_vector_type0_met_up + (jet_vector - jet_smeared_up)  ; ## new met 4V dn
                met_vector_type0_met_dn = met_vector_type0_met_dn + (jet_vector - jet_smeared_dn)  ; ## new met 4V dn

                met_vector_type2_met    = met_vector_type2_met    + (jet_vector - jet_smeared)  ;
                met_vector_type2_met_up = met_vector_type2_met_up + (jet_vector - jet_smeared_up)  ;
                met_vector_type2_met_dn = met_vector_type2_met_dn + (jet_vector - jet_smeared_dn)  ;


             self.pfMET_jer_[0]     = met_vector_type0_met.Pt(); 
             self.pfMET_Phi_jer_[0] = met_vector_type0_met.Phi(); 
              
             self.pfMET_jer_up_[0]     = met_vector_type0_met_up.Pt();  ## scaled up met
             self.pfMET_Phi_jer_up_[0] = met_vector_type0_met_up.Phi(); 

             self.pfMET_jer_dn_[0]     = met_vector_type0_met_dn.Pt(); ## scaled down met
             self.pfMET_Phi_jer_dn_[0] = met_vector_type0_met_dn.Phi(); 

             ### final invariant mass --> shape and normalization sys 

             self.mass_lvj_type0_met_jer_[0]    = (lepton_vector + jet_smeared_byJERC + met_vector_type0_met).M(); ## final invariant mass up
             self.mass_lvj_type0_met_jer_up_[0] = (lepton_vector + jet_smeared_byJERC_up + met_vector_type0_met_up).M(); ## final invariant mass up
             self.mass_lvj_type0_met_jer_dn_[0] = (lepton_vector + jet_smeared_byJERC_dn + met_vector_type0_met_dn).M(); ## final invariant mass dn
 
             self.mass_lvj_type2_met_jer_[0]    = (lepton_vector + jet_smeared_byJERC + met_vector_type2_met).M();
             self.mass_lvj_type2_met_jer_up_[0] = (lepton_vector + jet_smeared_byJERC_up + met_vector_type2_met_up).M();
             self.mass_lvj_type2_met_jer_dn_[0] = (lepton_vector + jet_smeared_byJERC_dn + met_vector_type2_met_dn).M();


             #####################################################################################################################################  
             ######## build met and lepton -> scale up and down the lepton according to the lepton energy -> buld the final invariant mass  ######
             #####################################################################################################################################
             
             met_vector_type0_met_up.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr( self.InputTree_, "W_nu1_pz_type0_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type0_met")*getattr( self.InputTree_,"W_nu1_pz_type0_met"))); ## original met 4V
             met_vector_type0_met_dn = met_vector_type0_met_up ;

             met_vector_type2_met_up.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                                getattr( self.InputTree_, "W_nu1_pz_type2_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type2_met")*getattr( self.InputTree_,"W_nu1_pz_type2_met"))); ## original met 4V, different pZ solution
             met_vector_type2_met_dn = met_vector_type2_met_up ;


             lepton_vector_scale_up = ROOT.TLorentzVector();
             lepton_vector_scale_dn = ROOT.TLorentzVector();

             if  getattr(self.InputTree_, "W_"+lepLabel+"_pt") < 200. :
              lepton_vector_scale_up.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px")*(1+self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_py")*(1+self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_pz")*(1+self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_e")*(1+self.LeptonScaleUnc_));

              lepton_vector_scale_dn.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px")*(1-self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_py")*(1-self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_pz")*(1-self.LeptonScaleUnc_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_e")*(1-self.LeptonScaleUnc_));
             else: 
              lepton_vector_scale_up.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px")*(1+self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_py")*(1+self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_pz")*(1+self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_e")*(1+self.LeptonScaleUncHighPT_));

              lepton_vector_scale_dn.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px")*(1-self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_py")*(1-self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_pz")*(1-self.LeptonScaleUncHighPT_),
                                                getattr(self.InputTree_, "W_"+lepLabel+"_e")*(1-self.LeptonScaleUncHighPT_));

             self.l_pt_scale_up_[0] = lepton_vector_scale_up.Pt(); 
             self.l_pt_scale_dn_[0] = lepton_vector_scale_dn.Pt(); 

             self.l_eta_scale_up_[0] = lepton_vector_scale_up.Eta(); 
             self.l_eta_scale_dn_[0] = lepton_vector_scale_dn.Eta(); 

             self.l_phi_scale_up_[0] = lepton_vector_scale_up.Phi(); 
             self.l_phi_scale_dn_[0] = lepton_vector_scale_dn.Phi(); 

             self.l_e_scale_up_[0] = lepton_vector_scale_up.E(); 
             self.l_e_scale_dn_[0] = lepton_vector_scale_dn.E(); 

             met_vector_type0_met_up = met_vector_type0_met_up + (lepton_vector - lepton_vector_scale_up)  ;
             met_vector_type0_met_dn = met_vector_type0_met_dn + (lepton_vector - lepton_vector_scale_dn)  ;

             met_vector_type2_met_up = met_vector_type2_met_up + (lepton_vector - lepton_vector_scale_up)  ;
             met_vector_type2_met_dn = met_vector_type2_met_dn + (lepton_vector - lepton_vector_scale_dn)  ;

             self.pfMET_lep_scale_up_[0]     = met_vector_type0_met_up.Pt();  ## scaled up met
             self.pfMET_Phi_lep_scale_up_[0] = met_vector_type0_met_up.Phi(); 

             self.pfMET_lep_scale_dn_[0]     = met_vector_type0_met_dn.Pt(); ## scaled down met
             self.pfMET_Phi_lep_scale_dn_[0] = met_vector_type0_met_dn.Phi(); 

             self.mass_lvj_type0_met_lep_scale_up_[0] = (lepton_vector_scale_up + jet_ptetaphie + met_vector_type0_met_up).M(); ## final invariant mass up
             self.mass_lvj_type0_met_lep_scale_dn_[0] = (lepton_vector_scale_dn + jet_ptetaphie + met_vector_type0_met_dn).M(); ## final invariant mass dn

             self.mass_lvj_type2_met_lep_scale_up_[0] = (lepton_vector_scale_up + jet_ptetaphie + met_vector_type2_met_up).M();
             self.mass_lvj_type2_met_lep_scale_dn_[0] = (lepton_vector_scale_dn + jet_ptetaphie + met_vector_type2_met_dn).M();

             #####################################################################################################################################  
             ######## build met and lepton -> smear the lepton according to the lepton energy resolution -> buld the final invariant mass  ######
             ######################################################################################################################################

             met_vector_type0_met.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr( self.InputTree_, "W_nu1_pz_type0_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type0_met")*getattr( self.InputTree_,"W_nu1_pz_type0_met"))); ## original met 4V

             met_vector_type2_met.SetPxPyPzE(getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                             getattr( self.InputTree_, "W_nu1_pz_type2_met" ),ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+getattr(self.InputTree_, "W_nu1_pz_type2_met")*getattr( self.InputTree_,"W_nu1_pz_type2_met"))); ## original met 4V

             lepton_vector_res = ROOT.TLorentzVector();
             
             smearedEnergy     = lepton_vector.E()*(1+self.Random_.Gaus(0.,self.LeptonResolutionUnc_));
             lepton_vector_res = lepton_vector*(smearedEnergy/lepton_vector.E());
             
             self.l_pt_res_[0]  = lepton_vector_res.Pt(); 
             self.l_eta_res_[0] = lepton_vector_res.Eta(); 
             self.l_phi_res_[0] = lepton_vector_res.Phi(); 
             self.l_e_res_[0]   = lepton_vector_res.M(); 

             met_vector_type0_met = met_vector_type0_met + (lepton_vector - lepton_vector_res)  ;
             met_vector_type2_met = met_vector_type2_met + (lepton_vector - lepton_vector_res)  ;

             self.pfMET_lep_res_[0]     = met_vector_type0_met.Pt();  ## scaled up met
             self.pfMET_Phi_lep_res_[0] = met_vector_type0_met.Phi(); 

             self.mass_lvj_type0_met_lep_res_[0] = (lepton_vector_res + jet_ptetaphie + met_vector_type0_met).M(); ## final invariant mass up
             self.mass_lvj_type2_met_lep_res_[0] = (lepton_vector_res + jet_ptetaphie + met_vector_type2_met).M();

             ######### btag stuff  
             index_ak5_cvst = array( 'f', [0.] );
             index_ak5_cvsm = array( 'f', [0.] );
             index_ak5_cvsl = array( 'f', [0.] );

             self.nbjets_csvt_veto_cleaned_[0] = 0. ;
             self.nbjets_csvm_veto_cleaned_[0] = 0. ;
             self.nbjets_csvl_veto_cleaned_[0] = 0. ;
             self.nbjets_csvl_veto_[0] = 0. ;
             self.nbjets_csvm_veto_[0] = 0. ;
             self.nbjets_csvt_veto_[0] = 0. ;
                 
             ## take the delta R between leading AK5 jet and the lepton 
             dR_lj = 0. ;
             for iJet in range(6):
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 30:

                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[iJet];
                                                                                
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );                
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );                

                        dR_jj = math.sqrt( (j_ak5_eta - getattr(self.InputTree_, prefix+"_eta")[0])**2 + deltaphi(j_ak5_phi,getattr(self.InputTree_, prefix+"_phi")[0])**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + deltaphi(l_phi,j_ak5_phi)**2 );

                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.244: self.nbjets_csvl_veto_[0] = self.nbjets_csvl_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.679: self.nbjets_csvm_veto_[0] = self.nbjets_csvm_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.898: self.nbjets_csvt_veto_[0] = self.nbjets_csvt_veto_[0]+1;
                        
                        if dR_jj > 0.8: 
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.244: index_ak5_cvsl.append(iJet);
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.679: index_ak5_cvsm.append(iJet);
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_cvst.append(iJet);

             for iJet in range(6):
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 30:
                        
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[iJet];

                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );                
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );                

                        dR_jj = math.sqrt( (j_ak5_eta - getattr(self.InputTree_, prefix+"_eta")[0])**2 + deltaphi(j_ak5_phi,getattr(self.InputTree_, prefix+"_phi")[0])**2 );
                        dR_lj = math.sqrt( (l_eta - j_ak5_eta)**2 + deltaphi(l_phi,j_ak5_phi)**2 );

                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] >=0.244: self.nbjets_csvl_veto_[0] = self.nbjets_csvl_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] >=0.679: self.nbjets_csvm_veto_[0] = self.nbjets_csvm_veto_[0]+1;
                        if getattr( self.InputTree_, "JetPFCorVBFTag_bDiscriminatorCSV" )[iJet] >=0.898: self.nbjets_csvt_veto_[0] = self.nbjets_csvt_veto_[0]+1;

                        if dR_jj > 0.8: 
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.244: index_ak5_cvsl.append(iJet);
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] >=0.679: index_ak5_cvsm.append(iJet);
                          if getattr( self.InputTree_, "JetPFCor_bDiscriminatorCSV" )[iJet] > 0.898: index_ak5_cvst.append(iJet);

               
             self.nbjets_csvl_veto_cleaned_[0] = len(index_ak5_cvsl);
             self.nbjets_csvm_veto_cleaned_[0] = len(index_ak5_cvsm);
             self.nbjets_csvt_veto_cleaned_[0] = len(index_ak5_cvst);

             self.nbjets_ssvhem_veto_  = getattr( self.InputTree_, "numPFCorJetBTags");
             self.njets_[0] = getattr( self.InputTree_, "numPFCorJets" ) + getattr( self.InputTree_, "numPFCorVBFTagJets" ) ;
               
             ############## top mass veto selection  -- > hadronic leg
             mass_top_veto_j1  = ROOT.TLorentzVector();
             mass_top_veto_j2  = ROOT.TLorentzVector();
             mass_top_veto     = ROOT.TLorentzVector();

             
             mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt")[0],getattr(self.InputTree_, prefix+"_eta")[0],getattr(self.InputTree_, prefix+"_phi")[0],
                                           getattr(self.InputTree_, prefix+"_mass")[0]);

             mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_,  "vbf_maxpt_j1_pt"),getattr(self.InputTree_, "vbf_maxpt_j1_eta"),getattr(self.InputTree_, "vbf_maxpt_j1_phi"),
                                           getattr(self.InputTree_,   "vbf_maxpt_j1_m"));
             mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
             self.mass_ungroomedjet_vbf_j1_[0]  = mass_top_veto.M();


             mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt_pr")[0],getattr(self.InputTree_, prefix+"_eta_pr")[0],getattr(self.InputTree_, prefix+"_phi_pr")[0],
                                           getattr(self.InputTree_, prefix+"_mass_pr")[0]);
             mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
             self.mass_ungroomedjet_vbf_j1_pr_[0]  = mass_top_veto.M();


             mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt")[0],getattr(self.InputTree_, prefix+"_eta")[0],getattr(self.InputTree_, prefix+"_phi")[0],
                                           getattr(self.InputTree_, prefix+"_mass")[0]);
             mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_,  "vbf_maxpt_j2_pt"),getattr(self.InputTree_, "vbf_maxpt_j2_eta"),getattr(self.InputTree_, "vbf_maxpt_j2_phi"),
                                           getattr(self.InputTree_,   "vbf_maxpt_j2_m"));
             mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
             self.mass_ungroomedjet_vbf_j2_[0]  = mass_top_veto.M();


             mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt_pr")[0],getattr(self.InputTree_, prefix+"_eta_pr")[0],getattr(self.InputTree_, prefix+"_phi_pr")[0],
                                           getattr(self.InputTree_, prefix+"_mass_pr")[0]);
             mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
             self.mass_ungroomedjet_vbf_j2_pr_[0]  = mass_top_veto.M();


             dR_jj_min = 100 ;
             dR_jj_pr_min = 100 ;
             ijet_index = -1 ;
             ijet_index_pr = -1 ;

             for iJet in range(6): ## loop on central jets 
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 30:

                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[iJet];
                                                                                
                        ca8_eta = getattr(self.InputTree_, prefix+"_eta")[0];
                        ca8_phi = getattr(self.InputTree_, prefix+"_phi")[0];                

                        ca8_eta_pr = getattr(self.InputTree_, prefix+"_eta_pr")[0];                
                        ca8_phi_pr = getattr(self.InputTree_, prefix+"_phi_pr")[0];

                        dR_jj = math.sqrt( (j_ak5_eta-ca8_eta)**2 + deltaphi(j_ak5_phi,ca8_phi)**2 );

                        if dR_jj < 0.8 : continue ; ## only jets outside ca8 cone
                        if(dR_jj < dR_jj_min): ## take the closer to the W
                            ijet_index = iJet;
                            dR_jj_min = dR_jj ;
                            
                        dR_jj_pr = math.sqrt( (ca8_eta_pr-j_ak5_eta)**2 + deltaphi(ca8_phi_pr,j_ak5_phi)**2 );
                        if dR_jj_pr < 0.8 : continue ;
                        if(dR_jj_pr < dR_jj_pr_min):                                
                            ijet_index_pr = iJet;
                            dR_jj_pr_min = dR_jj_pr ;

             for iJet in range(6): ## same thing looping on the forward jets
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 30:
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[iJet];
                                                                                
                        ca8_eta = getattr(self.InputTree_, prefix+"_eta")[0];
                        ca8_phi = getattr(self.InputTree_, prefix+"_phi")[0];                

                        ca8_eta_pr = getattr(self.InputTree_, prefix+"_eta_pr")[0];                
                        ca8_phi_pr = getattr(self.InputTree_, prefix+"_phi_pr")[0];

                        dR_jj = math.sqrt( (j_ak5_eta-ca8_eta)**2 + deltaphi(j_ak5_phi,ca8_phi)**2 );

                        if dR_jj < 0.8 : continue ;
                        if(dR_jj < dR_jj_min):
                            ijet_index = iJet+7;
                            dR_jj_min = dR_jj ;

                        dR_jj_pr = math.sqrt( (ca8_eta_pr-j_ak5_eta)**2 + deltaphi(ca8_phi_pr,j_ak5_phi)**2 );
                        if dR_jj_pr < 0.8 : continue ;
                        if(dR_jj_pr < dR_jj_pr_min):                                
                            ijet_index_pr = iJet+7;
                            dR_jj_pr_min = dR_jj_pr ;


             if(ijet_index < 7 and ijet_index !=-1): ## found a good jet in the central region

              mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt")[0],getattr(self.InputTree_, prefix+"_eta")[0],getattr(self.InputTree_, prefix+"_phi")[0],
                                            getattr(self.InputTree_, prefix+"_mass")[0]);

              mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_, "JetPFCor_Pt")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Eta")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Phi")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Mass")[ijet_index]);
    
              mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
              self.mass_ungroomedjet_closerjet_[0]  = mass_top_veto.M();


              mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt_pr")[0],getattr(self.InputTree_, prefix+"_eta_pr")[0],getattr(self.InputTree_, prefix+"_phi_pr")[0],
                                            getattr(self.InputTree_, prefix+"_mass_pr")[0]);
              mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
              self.mass_ungroomedjet_closerjet_pr_[0]  = mass_top_veto.M();


             elif (ijet_index !=-1 and ijet_index >= 7) : ## found a good jet in the forward region

              mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt")[0],getattr(self.InputTree_, prefix+"_eta")[0],getattr(self.InputTree_, prefix+"_phi")[0],
                                           getattr(self.InputTree_, prefix+"_mass")[0]);

              mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_, "JetPFCorVBFTag_Pt")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Eta")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Phi")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Mass")[ijet_index-7]);
    
              mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
              self.mass_ungroomedjet_closerjet_[0]  = mass_top_veto.M();

              mass_top_veto_j1.SetPtEtaPhiM(getattr(self.InputTree_, prefix+"_pt_pr")[0],getattr(self.InputTree_, prefix+"_eta_pr")[0],getattr(self.InputTree_, prefix+"_phi_pr")[0],
                                            getattr(self.InputTree_, prefix+"_mass_pr")[0]);
              mass_top_veto = mass_top_veto_j1 + mass_top_veto_j2;
              self.mass_ungroomedjet_closerjet_pr_[0]  = mass_top_veto.M();


             ### top mass veto --> leptonic leg
             dR_lj_min = 100 ;
             ijet_index = -1 ;

             leptonic_W = ROOT.TLorentzVector();
             leptonic_W.SetPxPyPzE(getattr(self.InputTree_, "W_"+lepLabel+"_px")+getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Cos(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                   getattr(self.InputTree_, "W_"+lepLabel+"_py")+getattr(self.InputTree_,"event_met_pfmet")*ROOT.TMath.Sin(getattr(self.InputTree_,"event_met_pfmetPhi")),
                                   getattr(self.InputTree_, "W_"+lepLabel+"_pz")+getattr( self.InputTree_, "W_nu1_pz_type0_met" ),
                                   getattr(self.InputTree_, "W_"+lepLabel+"_e")+ROOT.TMath.Sqrt(getattr(self.InputTree_, "event_met_pfmet")*getattr(self.InputTree_,"event_met_pfmet")+
                                   getattr(self.InputTree_, "W_nu1_pz_type0_met")*getattr( self.InputTree_,"W_nu1_pz_type0_met")));

             for iJet in range(6): ## loop on central jets
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCor_Pt" )[iJet] > 30:
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCor_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCor_Phi" )[iJet];
                                                                                
                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );                

                        ca8_eta = getattr(self.InputTree_, prefix+"_eta")[0];
                        ca8_phi = getattr(self.InputTree_, prefix+"_phi")[0];                

                        dR_jj = math.sqrt( (j_ak5_eta-ca8_eta)**2 + deltaphi(j_ak5_phi,ca8_phi)**2 );
                        dR_lj = math.sqrt( (j_ak5_eta-leptonic_W.Eta())**2 + deltaphi(j_ak5_phi,leptonic_W.Phi())**2 );

                        if dR_jj < 0.8 : continue ;
                        if(dR_lj < dR_lj_min):
                            ijet_index = iJet;
                            dR_lj_min = dR_lj ;


             for iJet in range(6):
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] < 0 : break ;                        
                    if getattr( self.InputTree_, "JetPFCorVBFTag_Pt" )[iJet] > 30:
                        j_ak5_eta = getattr( self.InputTree_, "JetPFCorVBFTag_Eta" )[iJet];
                        j_ak5_phi = getattr( self.InputTree_, "JetPFCorVBFTag_Phi" )[iJet];
                                                                                
                        ca8_eta = getattr(self.InputTree_, prefix+"_eta")[0];
                        ca8_phi = getattr(self.InputTree_, prefix+"_phi")[0];                

                        l_eta = getattr( self.InputTree_, "W_"+lepLabel+"_eta" );
                        l_phi = getattr( self.InputTree_, "W_"+lepLabel+"_phi" );                

                        dR_jj = math.sqrt( (j_ak5_eta-ca8_eta)**2 + deltaphi(j_ak5_phi,ca8_phi)**2 );
                        dR_lj = math.sqrt( (j_ak5_eta-leptonic_W.Eta())**2 + deltaphi(j_ak5_phi,leptonic_W.Phi())**2 );

                        if dR_jj < 0.8 : continue ;
                        if(dR_lj < dR_lj_min):
                            ijet_index = iJet;
                            dR_lj_min = dR_lj ;

             if(ijet_index < 7 and ijet_index !=-1):

              mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_, "JetPFCor_Pt")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Eta")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Phi")[ijet_index],
                                            getattr(self.InputTree_, "JetPFCor_Mass")[ijet_index]);
              mass_top_veto = leptonic_W + mass_top_veto_j2;
              self.mass_leptonic_closerjet_[0]  = mass_top_veto.M();


             elif (ijet_index !=-1 and ijet_index >= 7) :
              mass_top_veto_j2.SetPtEtaPhiM(getattr(self.InputTree_, "JetPFCorVBFTag_Pt")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Eta")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Phi")[ijet_index-7],
                                            getattr(self.InputTree_,  "JetPFCorVBFTag_Mass")[ijet_index-7]);
              mass_top_veto = leptonic_W + mass_top_veto_j2;
              self.mass_leptonic_closerjet_[0]  = mass_top_veto.M();

              
             ### Fill the output tree
             self.otree.Fill();

        ### close the file 
        self.OFile_.cd();
        self.otree.Write();
        self.OFile_.Close();

    ### Initialize the variables used to fill the output branches
    def InitializeVariables(self):

        ### event Property and weights
        self.nPV_               = array( 'f', [ 0. ] );                        
        self.event_runNo_       = array( 'i', [ 0 ] );
        self.event_lumi_        = array( 'i', [ 0 ] );
        self.event_             = array( 'i', [0] );

        self.totalEventWeight_  = array( 'f', [ 0. ] );                                
        self.eff_and_pu_Weight_ = array( 'f', [ 0. ] );                                
        self.wSampleWeight_     = array( 'f', [ 0. ] );
        self.event_weight_      = array( 'f', [ 0. ] );
        
        self.btag_weight_    = array( 'f', [0.]); 
        self.btag_weight_up_ = array( 'f', [0.]); 
        self.btag_weight_dn_ = array( 'f', [0.]); 
        self.btag_weight_up_dn_ = array( 'f', [0.]); 
        self.btag_weight_dn_up_ = array( 'f', [0.]); 

        self.interference_Weight_H600_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H700_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H800_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H900_  = array( 'f', [ 0. ] );                                
        self.interference_Weight_H1000_ = array( 'f', [ 0. ] );                                

        self.interferencevbf_1_ = array( 'f', [ 0. ] );
        self.interferencevbf_05_ = array( 'f', [ 0. ] );
        self.interferencevbf_2_ = array( 'f', [ 0. ] );
        
        self.cps_Weight_H600_  = array( 'f', [ 0. ] );                                
        self.cps_Weight_H700_  = array( 'f', [ 0. ] );                                
        self.cps_Weight_H800_  = array( 'f', [ 0. ] );                                
        self.cps_Weight_H900_  = array( 'f', [ 0. ] );                                
        self.cps_Weight_H1000_ = array( 'f', [ 0. ] );                                

        ########## bsm branches for ewk signlet 
        self.cprimeVals   = [1,2,3,4,5,6,7,8,9,10] ;
        self.brnewVals    = [00, 01, 02, 03, 04, 05] ;
        self.bsmReweights = [];
        
        for iPar in range(len(self.cprimeVals)): 
          col_bsmReweights = [];
          for jPar in range(len(self.brnewVals)): 
            col_bsmReweights.append( array( 'f', [ 0. ] ) );
          self.bsmReweights.append( col_bsmReweights );

        ########## event topology variables
        self.issignal_      = array( 'i', [ 0 ] );             
        self.numberJetBin_  = array( 'i', [ 0 ] );
        self.numberJetBin2_ = array( 'i', [ 0 ] );
        self.numberJetBin2_ = array( 'i', [ 0 ] );
        self.numberJetBin3_ = array( 'i', [ 0 ] );
        self.numberJetBin4_ = array( 'i', [ 0 ] );

        self.numberJetBinGen_  = array( 'i', [ 0 ] );
        self.numberJetBinGen2_ = array( 'i', [ 0 ] );
        self.numberJetBinGen2_ = array( 'i', [ 0 ] );
        self.numberJetBinGen3_ = array( 'i', [ 0 ] );
        self.numberJetBinGen4_ = array( 'i', [ 0 ] );

        ### Leptonic W, Lepton and Nueutrino
        self.l_pt_     = array( 'f', [ 0. ] );
        self.l_eta_    = array( 'f', [ 0. ] );
        self.l_phi_    = array( 'f', [ 0. ] );
        self.l_charge_ = array( 'f', [ 0. ] );

        self.l_pt_res_     = array( 'f', [ 0. ] );
        self.l_eta_res_    = array( 'f', [ 0. ] );
        self.l_phi_res_    = array( 'f', [ 0. ] );
        self.l_e_res_      = array( 'f', [ 0. ] );

        self.l_pt_scale_up_     = array( 'f', [ 0. ] );
        self.l_eta_scale_up_    = array( 'f', [ 0. ] );
        self.l_phi_scale_up_    = array( 'f', [ 0. ] );
        self.l_e_scale_up_      = array( 'f', [ 0. ] );

        self.l_pt_scale_dn_     = array( 'f', [ 0. ] );
        self.l_eta_scale_dn_    = array( 'f', [ 0. ] );
        self.l_phi_scale_dn_    = array( 'f', [ 0. ] );
        self.l_e_scale_dn_      = array( 'f', [ 0. ] );

        self.pfMET_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_ = array( 'f', [ 0. ] );

        self.pfMET_jes_up_       = array( 'f', [ 0. ] );
        self.pfMET_Phi_jes_up_   = array( 'f', [ 0. ] );

        self.pfMET_jes_dn_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_jes_dn_ = array( 'f', [ 0. ] );

        self.pfMET_jer_       = array( 'f', [ 0. ] );
        self.pfMET_Phi_jer_   = array( 'f', [ 0. ] );

        self.pfMET_jer_up_       = array( 'f', [ 0. ] );
        self.pfMET_Phi_jer_up_   = array( 'f', [ 0. ] );

        self.pfMET_jer_dn_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_jer_dn_ = array( 'f', [ 0. ] );

        self.pfMET_lep_scale_up_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_lep_scale_up_ = array( 'f', [ 0. ] );

        self.pfMET_lep_scale_dn_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_lep_scale_dn_ = array( 'f', [ 0. ] );

        self.pfMET_lep_res_     = array( 'f', [ 0. ] );
        self.pfMET_Phi_lep_res_ = array( 'f', [ 0. ] );

        self.nu_pz_type0_ = array( 'f', [ 0. ] );
        self.nu_pz_type2_ = array( 'f', [ 0. ] );

        self.nu_pz_type0_met_ = array( 'f', [ 0. ] );
        self.nu_pz_type2_met_ = array( 'f', [ 0. ] );

        self.W_pz_type0_ = array( 'f', [ 0. ] );
        self.W_pz_type2_ = array( 'f', [ 0. ] );

        self.W_pz_type0_met_ = array( 'f', [ 0. ] );
        self.W_pz_type2_met_ = array( 'f', [ 0. ] );

        self.nu_pz_gen_   = array( 'f', [ 0. ] );
        self.W_pz_gen_    = array( 'f', [ 0. ] );
        self.W_pt_gen_    = array( 'f', [ 0. ] );


        ###### new variables for tree branches         
        self.mass_lvj_type0_met_    = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_    = array( 'f', [ 0. ] );

        self.mass_lvj_type0_    = array( 'f', [ 0. ] );
        self.mass_lvj_type2_    = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_jes_up_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_jes_up_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_jes_dn_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_jes_dn_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_jer_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_jer_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_jer_up_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_jer_up_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_jer_dn_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_jer_dn_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_lep_scale_up_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_lep_scale_up_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_lep_scale_dn_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_lep_scale_dn_  = array( 'f', [ 0. ] );

        self.mass_lvj_type0_met_lep_res_  = array( 'f', [ 0. ] );
        self.mass_lvj_type2_met_lep_res_  = array( 'f', [ 0. ] );

        self.mass_lv_subj_type0_met_    = array( 'f', [ 0. ] );
        self.mass_lv_subj_type2_met_    = array( 'f', [ 0. ] );

        self.mass_lv_subj_type0_    = array( 'f', [ 0. ] );
        self.mass_lv_subj_type2_    = array( 'f', [ 0. ] );

        ### leptonic W transverse mass and pt 
        self.v_pt_             = array( 'f', [ 0. ] );
        self.v_mt_             = array( 'f', [ 0. ] );
        self.v_eta_            = array( 'f', [ 0. ] );
        self.v_phi_            = array( 'f', [ 0. ] );

        ######### Hadronic W Variables ###############
        self.ungroomed_jet_eta_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_pt_  = array( 'f', [ 0. ] );
        self.ungroomed_jet_phi_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_e_   = array( 'f', [ 0. ] );

        self.ungroomed_gen_jet_pt_  = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_eta_ = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_phi_ = array( 'f', [ 0. ] );
        self.ungroomed_gen_jet_e_   = array( 'f', [ 0. ] );

        self.jet_mass_pr_       = array( 'f', [ 0. ] );
        self.jet_pt_pr_         = array( 'f', [ 0. ] );
        self.jet_charge_        = array( 'f', [ 0. ] );
        self.jet_charge_k05_    = array( 'f', [ 0. ] );
        self.jet_charge_k07_    = array( 'f', [ 0. ] );
        self.jet_charge_k10_    = array( 'f', [ 0. ] );

        self.gen_jet_mass_pr_ = array( 'f', [ 0. ] );
        self.gen_jet_pt_pr_   = array( 'f', [ 0. ] );

        self.jet_grsens_ft_   = array( 'f', [ 0. ] );
        self.jet_grsens_tr_   = array( 'f', [ 0. ] );
        self.jet_massdrop_pr_ = array( 'f', [ 0. ] );    
        self.jet_qjetvol_     = array( 'f', [ 0. ] ); 

        self.gen_jet_grsens_ft_   = array( 'f', [ 0. ] );
        self.gen_jet_grsens_tr_   = array( 'f', [ 0. ] );
        self.gen_jet_massdrop_pr_ = array( 'f', [ 0. ] );    
        self.gen_jet_qjetvol_     = array( 'f', [ 0. ] ); 

        self.jet_tau2tau1_       = array( 'f', [ 0. ] );     
        self.jet_tau2tau1_exkT_  = array( 'f', [ 0. ] );     
        self.jet_tau2tau1_pr_    = array( 'f', [ 0. ] );
        self.jet_GeneralizedECF_ = array( 'f', [ 0. ] );

        self.gen_jet_tau2tau1_       = array( 'f', [ 0. ] );     
        self.gen_jet_tau2tau1_exkT_  = array( 'f', [ 0. ] );     
        self.gen_jet_tau2tau1_pr_    = array( 'f', [ 0. ] );
        self.gen_jet_GeneralizedECF_ = array( 'f', [ 0. ] );

        self.jet_jetconstituents_      = array( 'f', [ 0. ] );     
        self.gen_jet_jetconstituents_  = array( 'f', [ 0. ] );     
        
        self.jet_rcore4_ = array( 'f', [ 0. ] );
        self.jet_rcore5_ = array( 'f', [ 0. ] );
        self.jet_rcore6_ = array( 'f', [ 0. ] );
        self.jet_rcore7_ = array( 'f', [ 0. ] );

        self.gen_jet_rcore4_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore5_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore6_ = array( 'f', [ 0. ] );
        self.gen_jet_rcore7_ = array( 'f', [ 0. ] );

        self.jet_pt1frac_  = array ('f',[ 0. ]);
        self.jet_pt2frac_  = array ('f',[ 0. ]);
        self.jet_sjdr_     = array ('f',[ 0. ]);
        
        self.j_jecfactor_up_ = array( 'f', [ 0. ] );
        self.j_jecfactor_dn_ = array( 'f', [ 0. ] );

        self.jet_mass_pr_jes_up_ = array( 'f', [ 0. ] );
        self.jet_mass_pr_jes_dn_ = array( 'f', [ 0. ] );

        self.ungroomed_jet_pt_jes_up_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_pt_jes_dn_ = array( 'f', [ 0. ] );

        self.jet_mass_pr_jer_    = array( 'f', [ 0. ] );
        self.jet_mass_pr_jer_up_ = array( 'f', [ 0. ] );
        self.jet_mass_pr_jer_dn_ = array( 'f', [ 0. ] );

        self.ungroomed_jet_pt_jer_    = array( 'f', [ 0. ] );
        self.ungroomed_jet_pt_jer_up_ = array( 'f', [ 0. ] );
        self.ungroomed_jet_pt_jer_dn_ = array( 'f', [ 0. ] );

        self.jet_planarlow04_ = array( 'f', [0.] );
        self.jet_planarlow05_ = array( 'f', [0.] );
        self.jet_planarlow06_ = array( 'f', [0.] );
        self.jet_planarlow07_ = array( 'f', [0.] );

        ###### top mass veto variables
        self.mass_ungroomedjet_vbf_j1_     = array( 'f', [0.] );
        self.mass_ungroomedjet_vbf_j2_     = array( 'f', [0.] );
        self.mass_ungroomedjet_vbf_j1_pr_  = array( 'f', [0.] );
        self.mass_ungroomedjet_vbf_j2_pr_  = array( 'f', [0.] );
        self.mass_ungroomedjet_closerjet_  = array( 'f', [0.] );
        self.mass_ungroomedjet_closerjet_pr_ = array( 'f', [0.] );

        self.mass_leptonic_closerjet_  = array( 'f', [0.] );

        ###### vbf variables
        self.vbf_maxpt_jj_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_jj_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_jj_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_jes_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_jes_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_jes_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_jes_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_jes_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_jes_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_jes_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_jes_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_jer_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_jer_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_jer_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_jer_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_jer_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_jer_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_jer_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_jer_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_m_jer_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_pt_jer_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_eta_jer_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j1_phi_jer_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_gen_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_gen_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_gen_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_jes_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_jes_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_jes_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_jes_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_jes_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_jes_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_jes_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_jes_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_jer_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_jer_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_jer_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_jer_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_jer_up_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_jer_up_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_jer_up_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_jer_up_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j2_m_jer_dn_   = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_pt_jer_dn_  = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_eta_jer_dn_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_phi_jer_dn_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_QGLikelihood_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_QGLikelihood_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j1_isPileUpTight_  = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j2_isPileUpMedium_ = array( 'i', [ 0 ] );                                
        self.vbf_maxpt_j2_isPileUpTight_  = array( 'i', [ 0 ] );                                

        self.vbf_maxpt_j1_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_bDiscriminatorCSV_ = array( 'f', [ 0. ] );                                

        self.vbf_maxpt_j1_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                
        self.vbf_maxpt_j2_bDiscriminatorCSV_gen_ = array( 'f', [ 0. ] );                                

        ###### btag counters
        self.nbjets_csvl_veto_   = array( 'f', [ 0. ] );
        self.nbjets_csvm_veto_   = array( 'f', [ 0. ] );
        self.nbjets_csvt_veto_   = array( 'f', [ 0. ] );
        self.nbjets_ssvhem_veto_ = array( 'f', [ 0. ] );

        self.nbjets_csvl_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_csvm_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_csvt_veto_cleaned_   = array( 'f', [ 0. ] );
        self.nbjets_ssvhem_veto_cleaned_ = array( 'f', [ 0. ] );

        self.njets_ = array( 'f', [ 0. ] );

        self.deltaR_lca8jet_          = array( 'f', [0.] );
        self.deltaphi_METca8jet_      = array( 'f', [0.] );
        self.deltaphi_Vca8jet_        = array( 'f', [0.] );
        self.deltaphi_METca8jet_met_  = array( 'f', [0.] );
        self.deltaphi_Vca8jet_met_    = array( 'f', [0.] );

        ### some generator level quantites
        self.genHMass_ = array( 'f', [ 0. ] );              
        self.genHphi_  = array( 'f', [ 0. ] );              
        self.genHeta_  = array( 'f', [ 0. ] );              
        self.genHpt_   = array( 'f', [ 0. ] );              

        self.genTagQuark1E_    = array( 'f', [ 0. ] );              
        self.genTagQuark1phi_  = array( 'f', [ 0. ] );              
        self.genTagQuark1eta_  = array( 'f', [ 0. ] );              
        self.genTagQuark1pt_   = array( 'f', [ 0. ] );              

        self.genTagQuark2E_    = array( 'f', [ 0. ] );              
        self.genTagQuark2phi_  = array( 'f', [ 0. ] );              
        self.genTagQuark2eta_  = array( 'f', [ 0. ] );              
        self.genTagQuark2pt_   = array( 'f', [ 0. ] );              
                                                                        
        # variables for ttbar control region
        self.ttb_nak5_same_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_same_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_oppo_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppo_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_oppoveto_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_oppoveto_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_nak5_sameveto_       = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvl_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvm_  = array( 'f', [ 0. ] );        
        self.ttb_nak5_sameveto_csvt_  = array( 'f', [ 0. ] );        

        self.ttb_ht_           = array( 'f', [ 0. ] );                
        self.ttb_ca8_mass_pr_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_   = array( 'f', [ 0. ] );
        
        self.ttb_ca8_charge_k05_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_k07_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_charge_k10_  = array( 'f', [ 0. ] );                        

        self.ttb_ca8_px_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_py_  = array( 'f', [ 0. ] );                                
        self.ttb_ca8_pz_  = array( 'f', [ 0. ] );                                
        self.ttb_ca8_e_   = array( 'f', [ 0. ] );                                
        
        self.ttb_ca8_ungroomed_pt_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_eta_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_phi_  = array( 'f', [ 0. ] );
        self.ttb_ca8_ungroomed_e_    = array( 'f', [ 0. ] );

        self.ttb_ca8_ungroomed_gen_pt_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_gen_eta_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_ungroomed_gen_phi_  = array( 'f', [ 0. ] );
        self.ttb_ca8_ungroomed_gen_e_    = array( 'f', [ 0. ] );

        self.ttb_ca8_tau2tau1_        = array( 'f', [ 0. ] );                        
        self.ttb_ca8_tau2tau1_exkT_   = array( 'f', [ 0. ] );                        
        self.ttb_ca8_tau2tau1_pr_     = array( 'f', [ 0. ] );                        
        self.ttb_ca8_GeneralizedECF_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mu_              = array( 'f', [ 0. ] );                        

        self.ttb_ca8_mlvj_type0_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mlvj_type2_  = array( 'f', [ 0. ] );                        

        self.ttb_ca8_mlvj_type0_met_  = array( 'f', [ 0. ] );                        
        self.ttb_ca8_mlvj_type2_met_  = array( 'f', [ 0. ] );                        

        self.isttbar_ = array( 'i', [ 0 ] );            

        self.ttb_dR_ca8_bjet_closer_ = array( 'f', [ 0. ] );            
        self.ttb_dR_ca8_jet_closer_  = array( 'f', [ 0. ] );            

        self.gen_parton1_px_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_py_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_pz_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton1_e_fromttbar_  = array( 'f', [ 0. ] );
        self.gen_parton1_id_fromttbar_ = array( 'i', [ 0 ] );
        self.gen_parton2_px_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_py_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_pz_fromttbar_ = array( 'f', [ 0. ] );
        self.gen_parton2_e_fromttbar_  = array( 'f', [ 0. ] );
        self.gen_parton2_id_fromttbar_ = array( 'i', [ 0 ] );

    def createBranches(self):

        ################ Branch output Tree
        self.otree.Branch("nPV", self.nPV_,"nPV/F");
        self.otree.Branch("event", self.event_ , "event/I");
        self.otree.Branch("event_lumi", self.event_lumi_ , "event_lumi/I");
        self.otree.Branch("event_runNo", self.event_runNo_ , "event_runNo/I");

        self.otree.Branch("issignal", self.issignal_ , "issignal/I");
        self.otree.Branch("numberJetBin", self.numberJetBin_ , "numberJetBin/I");
        self.otree.Branch("numberJetBin2", self.numberJetBin2_ , "numberJetBin2/I");
        self.otree.Branch("numberJetBin3", self.numberJetBin3_ , "numberJetBin3/I");
        self.otree.Branch("numberJetBin4", self.numberJetBin4_ , "numberJetBin4/I");

        self.otree.Branch("numberJetBinGen", self.numberJetBinGen_ , "numberJetBinGen/I");
        self.otree.Branch("numberJetBinGen2", self.numberJetBinGen2_ , "numberJetBinGen2/I");
        self.otree.Branch("numberJetBinGen3", self.numberJetBinGen3_ , "numberJetBinGen3/I");
        self.otree.Branch("numberJetBinGen4", self.numberJetBinGen4_ , "numberJetBinGen4/I");

        self.otree.Branch("totalEventWeight",  self.totalEventWeight_ , "totalEventWeight/F");   ## total xs * pile Up * trigger * lepton ID
        self.otree.Branch("eff_and_pu_Weight", self.eff_and_pu_Weight_ , "eff_and_pu_Weight/F"); ## product of pileUp, trigger and lepton ID
        self.otree.Branch("wSampleWeight",     self.wSampleWeight_ , "wSampleWeight/F");         ## only cross section weight
        self.otree.Branch("event_weight",      self.event_weight_ , "event_weight/F");           ## in case the MC produce weighted events at LHE level

        self.otree.Branch("btag_weight", self.btag_weight_ , "btag_weight/F"); ## btag weight according to btag pog recipe
        self.otree.Branch("btag_weight_up", self.btag_weight_up_ , "btag_weight_up/F"); 
        self.otree.Branch("btag_weight_dn", self.btag_weight_dn_ , "btag_weight_dn/F");
        self.otree.Branch("btag_weight_dn_up", self.btag_weight_dn_up_ , "btag_weight_dn_up/F");
        self.otree.Branch("btag_weight_up_dn", self.btag_weight_up_dn_ , "btag_weight_up_dn/F");

        self.otree.Branch("interference_Weight_H600", self.interference_Weight_H600_ , "interference_Weight_H600/F");
        self.otree.Branch("interference_Weight_H700", self.interference_Weight_H700_ , "interference_Weight_H700/F");
        self.otree.Branch("interference_Weight_H800", self.interference_Weight_H800_ , "interference_Weight_H800/F");
        self.otree.Branch("interference_Weight_H900", self.interference_Weight_H900_ , "interference_Weight_H900/F");
        self.otree.Branch("interference_Weight_H1000", self.interference_Weight_H1000_ , "interference_Weight_H1000/F");

        self.otree.Branch("interferencevbf_1", self.interferencevbf_1_ , "interferencevbf_1/F");
        self.otree.Branch("interferencevbf_2", self.interferencevbf_2_ , "interferencevbf_2/F");
        self.otree.Branch("interferencevbf_05", self.interferencevbf_05_ , "interferencevbf_05/F");


        self.otree.Branch("cps_Weight_H600", self.cps_Weight_H600_ , "cps_Weight_H600/F");
        self.otree.Branch("cps_Weight_H700", self.cps_Weight_H700_ , "cps_Weight_H700/F");
        self.otree.Branch("cps_Weight_H800", self.cps_Weight_H800_ , "cps_Weight_H800/F");
        self.otree.Branch("cps_Weight_H900", self.cps_Weight_H900_ , "cps_Weight_H900/F");
        self.otree.Branch("cps_Weight_H1000", self.cps_Weight_H1000_ , "cps_Weight_H1000/F");

        ### branches for bsm models 
        for iPar in range(len(self.cprimeVals)): 
            for jPar in range(len(self.brnewVals)): 
                brname = "bsmReweight_cPrime%02d_brNew%02d"%(self.cprimeVals[iPar],self.brnewVals[jPar]);
                self.otree.Branch(brname,self.bsmReweights[iPar][jPar],brname);        


        ##############
        self.otree.Branch("mass_lvj_type0_met", self.mass_lvj_type0_met_ , "mass_lvj_type0_met/F");
        self.otree.Branch("mass_lvj_type2_met", self.mass_lvj_type2_met_ , "mass_lvj_type2_met/F");

        self.otree.Branch("mass_lvj_type0", self.mass_lvj_type0_ , "mass_lvj_type0/F");
        self.otree.Branch("mass_lvj_type2", self.mass_lvj_type2_ , "mass_lvj_type2/F");

        self.otree.Branch("mass_lvj_type0_met_jes_up", self.mass_lvj_type0_met_jes_up_ , "mass_lvj_type0_met_jes_up/F");
        self.otree.Branch("mass_lvj_type2_met_jes_up", self.mass_lvj_type2_met_jes_up_ , "mass_lvj_type2_met_jes_up/F");

        self.otree.Branch("mass_lvj_type0_met_jes_dn", self.mass_lvj_type0_met_jes_dn_ , "mass_lvj_type0_met_jes_dn/F");
        self.otree.Branch("mass_lvj_type2_met_jes_dn", self.mass_lvj_type2_met_jes_dn_ , "mass_lvj_type2_met_jes_dn/F");

        self.otree.Branch("mass_lvj_type0_met_jer", self.mass_lvj_type0_met_jer_ , "mass_lvj_type0_met_jer/F");
        self.otree.Branch("mass_lvj_type2_met_jer", self.mass_lvj_type2_met_jer_ , "mass_lvj_type2_met_jer/F");

        self.otree.Branch("mass_lvj_type0_met_jer_up", self.mass_lvj_type0_met_jer_up_ , "mass_lvj_type0_met_jer_up/F");
        self.otree.Branch("mass_lvj_type2_met_jer_up", self.mass_lvj_type2_met_jer_up_ , "mass_lvj_type2_met_jer_up/F");

        self.otree.Branch("mass_lvj_type0_met_jer_dn", self.mass_lvj_type0_met_jer_dn_ , "mass_lvj_type0_met_jer_dn/F");
        self.otree.Branch("mass_lvj_type2_met_jer_dn", self.mass_lvj_type2_met_jer_dn_ , "mass_lvj_type2_met_jer_dn/F");

        self.otree.Branch("mass_lvj_type0_met_lep_scale_up", self.mass_lvj_type0_met_lep_scale_up_ , "mass_lvj_type0_met_lep_scale_up/F");
        self.otree.Branch("mass_lvj_type2_met_lep_scale_up", self.mass_lvj_type2_met_lep_scale_up_ , "mass_lvj_type2_met_lep_scale_up/F");

        self.otree.Branch("mass_lvj_type0_met_lep_scale_dn", self.mass_lvj_type0_met_lep_scale_dn_ , "mass_lvj_type0_met_lep_scale_dn/F");
        self.otree.Branch("mass_lvj_type2_met_lep_scale_dn", self.mass_lvj_type2_met_lep_scale_dn_ , "mass_lvj_type2_met_lep_scale_dn/F");

        self.otree.Branch("mass_lvj_type0_met_lep_res", self.mass_lvj_type0_met_lep_res_ , "mass_lvj_type0_met_lep_res/F");
        self.otree.Branch("mass_lvj_type2_met_lep_res", self.mass_lvj_type2_met_lep_res_ , "mass_lvj_type2_met_lep_res/F");

        self.otree.Branch("mass_lv_subj_type0_met", self.mass_lv_subj_type0_met_ , "mass_lv_subj_type0_met/F");
        self.otree.Branch("mass_lv_subj_type2_met", self.mass_lv_subj_type2_met_ , "mass_lv_subj_type2_met/F");

        self.otree.Branch("mass_lv_subj_type0", self.mass_lv_subj_type0_ , "mass_lv_subj_type0/F");
        self.otree.Branch("mass_lv_subj_type2", self.mass_lv_subj_type2_ , "mass_lv_subj_type2/F");

        ###########
        self.otree.Branch("mass_ungroomedjet_vbf_j1", self.mass_ungroomedjet_vbf_j1_ , "mass_ungroomedjet_vbf_j1/F");
        self.otree.Branch("mass_ungroomedjet_vbf_j2", self.mass_ungroomedjet_vbf_j2_ , "mass_ungroomedjet_vbf_j2/F");
        self.otree.Branch("mass_ungroomedjet_vbf_j1_pr", self.mass_ungroomedjet_vbf_j1_pr_ , "mass_ungroomedjet_vbf_j1_pr/F");
        self.otree.Branch("mass_ungroomedjet_vbf_j2_pr", self.mass_ungroomedjet_vbf_j2_pr_ , "mass_ungroomedjet_vbf_j2_pr/F");
        self.otree.Branch("mass_ungroomedjet_closerjet", self.mass_ungroomedjet_closerjet_ , "mass_ungroomedjet_closerjet/F");
        self.otree.Branch("mass_ungroomedjet_closerjet_pr", self.mass_ungroomedjet_closerjet_pr_ , "mass_ungroomedjet_closerjet_pr/F");

        self.otree.Branch("mass_leptonic_closerjet", self.mass_leptonic_closerjet_ , "mass_leptonic_closerjet/F");

        ########### Leptonic W, Lepton and Nueutrino
        self.otree.Branch("l_pt", self.l_pt_ , "l_pt/F");
        self.otree.Branch("l_eta", self.l_eta_ , "l_eta/F");
        self.otree.Branch("l_charge", self.l_charge_ , "l_charge/F");
        self.otree.Branch("l_phi", self.l_phi_ , "l_phi/F");

        self.otree.Branch("l_pt_scale_up", self.l_pt_scale_up_ , "l_pt_scale_up/F");
        self.otree.Branch("l_eta_scale_up", self.l_eta_scale_up_ , "l_eta_scale_up/F");
        self.otree.Branch("l_phi_scale_up", self.l_phi_scale_up_ , "l_phi_scale_up/F");

        self.otree.Branch("l_pt_scale_dn", self.l_pt_scale_dn_ , "l_pt_scale_dn/F");
        self.otree.Branch("l_eta_scale_dn", self.l_eta_scale_dn_ , "l_eta_scale_dn/F");
        self.otree.Branch("l_phi_scale_dn", self.l_phi_scale_dn_ , "l_phi_scale_dn/F");

        self.otree.Branch("l_pt_res", self.l_pt_res_ , "l_pt_res/F");
        self.otree.Branch("l_eta_res", self.l_eta_res_ , "l_eta_res/F");
        self.otree.Branch("l_phi_res", self.l_phi_res_ , "l_phi_res/F");

        self.otree.Branch("pfMET", self.pfMET_ , "pfMET/F");
        self.otree.Branch("pfMET_Phi", self.pfMET_Phi_ , "pfMET_Phi/F");

        self.otree.Branch("pfMET_jes_up", self.pfMET_jes_up_ , "pfMET_jes_up/F");
        self.otree.Branch("pfMET_Phi_jes_up", self.pfMET_Phi_jes_up_ , "pfMET_Phi_jes_up/F");

        self.otree.Branch("pfMET_jes_dn", self.pfMET_jes_dn_ , "pfMET_jes_dn/F");
        self.otree.Branch("pfMET_Phi_jes_dn", self.pfMET_Phi_jes_dn_ , "pfMET_Phi_jes_dn/F");

        self.otree.Branch("pfMET_jer_", self.pfMET_jer_ , "pfMET_jer/F");
        self.otree.Branch("pfMET_Phi_jer_", self.pfMET_Phi_jer_ , "pfMET_Phi_jer/F");

        self.otree.Branch("pfMET_jer_up", self.pfMET_jer_up_ , "pfMET_jer_up/F");
        self.otree.Branch("pfMET_Phi_jer_up", self.pfMET_Phi_jer_up_ , "pfMET_Phi_jer_up/F");

        self.otree.Branch("pfMET_jer_dn", self.pfMET_jer_dn_ , "pfMET_jer_dn/F");
        self.otree.Branch("pfMET_Phi_jer_dn", self.pfMET_Phi_jer_dn_ , "pfMET_Phi_jer_dn/F");

        self.otree.Branch("pfMET_lep_scale_up", self.pfMET_lep_scale_up_ , "pfMET_lep_scale_up/F");
        self.otree.Branch("pfMET_Phi_lep_scale_up", self.pfMET_Phi_lep_scale_up_ , "pfMET_Phi_lep_scale_up/F");

        self.otree.Branch("pfMET_lep_scale_dn", self.pfMET_lep_scale_dn_ , "pfMET_lep_scale_dn/F");
        self.otree.Branch("pfMET_Phi_lep_scale_dn", self.pfMET_Phi_lep_scale_dn_ , "pfMET_Phi_lep_scale_dn/F");

        self.otree.Branch("pfMET_lep_res", self.pfMET_lep_res_ , "pfMET_lep_res/F");
        self.otree.Branch("pfMET_Phi_lep_res", self.pfMET_Phi_lep_res_ , "pfMET_Phi_lep_res/F");

        self.otree.Branch("nu_pz_type0", self.nu_pz_type0_ , "nu_pz_type0/F");
        self.otree.Branch("nu_pz_type2", self.nu_pz_type2_ , "nu_pz_type2/F");

        self.otree.Branch("nu_pz_type0_met", self.nu_pz_type0_met_ , "nu_pz_type0_met/F");
        self.otree.Branch("nu_pz_type2_met", self.nu_pz_type2_met_ , "nu_pz_type2_met/F");

        self.otree.Branch("W_pz_type0", self.W_pz_type0_ , "W_pz_type0/F");
        self.otree.Branch("W_pz_type2", self.W_pz_type2_ , "W_pz_type2/F");

        self.otree.Branch("W_pz_type0_met", self.W_pz_type0_met_ , "W_pz_type0_met/F");
        self.otree.Branch("W_pz_type2_met", self.W_pz_type2_met_ , "W_pz_type2_met/F");

        self.otree.Branch("nu_pz_gen", self.nu_pz_gen_ , "nu_pz_gen/F");
        self.otree.Branch("W_pz_gen", self.W_pz_gen_ , "W_pz_gen/F");
        self.otree.Branch("W_pt_gen", self.W_pt_gen_ , "W_pt_gen/F");

        self.otree.Branch("v_pt", self.v_pt_ , "v_pt/F");
        self.otree.Branch("v_mt", self.v_mt_ , "v_mt/F");
        self.otree.Branch("v_eta", self. v_eta_ , "v_eta/F");
        self.otree.Branch("v_phi", self.v_phi_ , "v_phi/F");

        ###################
        self.otree.Branch("ungroomed_jet_eta", self.ungroomed_jet_eta_ , "ungroomed_jet_eta/F");
        self.otree.Branch("ungroomed_jet_phi", self.ungroomed_jet_phi_ , "ungroomed_jet_phi/F");
        self.otree.Branch("ungroomed_jet_pt", self.ungroomed_jet_pt_ , "ungroomed_jet_pt/F");
        self.otree.Branch("ungroomed_jet_e", self.ungroomed_jet_e_ , "ungroomed_jet_e/F");

        self.otree.Branch("ungroomed_gen_jet_eta", self.ungroomed_gen_jet_eta_ , "ungroomed_gen_jet_eta/F");
        self.otree.Branch("ungroomed_gen_jet_phi", self.ungroomed_gen_jet_phi_ , "ungroomed_gen_jet_phi/F");
        self.otree.Branch("ungroomed_gen_jet_pt", self.ungroomed_gen_jet_pt_ , "ungroomed_gen_jet_pt/F");
        self.otree.Branch("ungroomed_gen_jet_e", self.ungroomed_gen_jet_e_ , "ungroomed_gen_jet_e/F");

        self.otree.Branch("jet_mass_pr", self.jet_mass_pr_ , "jet_mass_pr/F");
        self.otree.Branch("jet_pt_pr", self.jet_pt_pr_ , "jet_pt_pr/F");
        self.otree.Branch("jet_charge", self.jet_charge_ , "jet_charge/F");
        self.otree.Branch("jet_charge_k05", self.jet_charge_k05_ , "jet_charge_k05/F");
        self.otree.Branch("jet_charge_k07", self.jet_charge_k07_ , "jet_charge_k07/F");
        self.otree.Branch("jet_charge_k10", self.jet_charge_k10_ , "jet_charge_k10/F");

        self.otree.Branch("gen_jet_mass_pr", self.gen_jet_mass_pr_ , "gen_jet_mass_pr/F");
        self.otree.Branch("gen_jet_pt_pr", self.gen_jet_pt_pr_ , "gen_jet_pt_pr/F");

        self.otree.Branch("jet_grsens_ft", self.jet_grsens_ft_ , "jet_grsens_ft/F");
        self.otree.Branch("jet_grsens_tr", self.jet_grsens_tr_ , "jet_grsens_tr/F");
        self.otree.Branch("jet_massdrop_pr", self.jet_massdrop_pr_ , "jet_massdrop_pr/F");
        self.otree.Branch("jet_qjetvol", self.jet_qjetvol_ , "jet_qjetvol/F");

        self.otree.Branch("gen_jet_grsens_ft", self.gen_jet_grsens_ft_ , "gen_jet_grsens_ft/F");
        self.otree.Branch("gen_jet_grsens_tr", self.gen_jet_grsens_tr_ , "gen_jet_grsens_tr/F");
        self.otree.Branch("gen_jet_massdrop_pr", self.gen_jet_massdrop_pr_ , "gen_jet_massdrop_pr/F");
        self.otree.Branch("gen_jet_qjetvol", self.gen_jet_qjetvol_ , "gen_jet_qjetvol/F");

        self.otree.Branch("jet_tau2tau1", self.jet_tau2tau1_ , "jet_tau2tau1/F");
        self.otree.Branch("jet_tau2tau1_exkT", self.jet_tau2tau1_exkT_ , "jet_tau2tau1_exkT/F");
        self.otree.Branch("jet_tau2tau1_pr", self.jet_tau2tau1_pr_ , "jet_tau2tau1_pr/F");
        self.otree.Branch("jet_GeneralizedECF", self.jet_GeneralizedECF_ , "jet_GeneralizedECF/F");

        self.otree.Branch("gen_jet_tau2tau1", self.gen_jet_tau2tau1_ , "gen_jet_tau2tau1/F");
        self.otree.Branch("gen_jet_tau2tau1_exkT", self.gen_jet_tau2tau1_exkT_ , "gen_jet_tau2tau1_exkT/F");
        self.otree.Branch("gen_jet_tau2tau1_pr", self.gen_jet_tau2tau1_pr_ , "gen_jet_tau2tau1_pr/F");
        self.otree.Branch("gen_jet_GeneralizedECF", self.gen_jet_GeneralizedECF_ , "gen_jet_GeneralizedECF/F");

        self.otree.Branch("jet_jetconstituents", self.jet_jetconstituents_ , "jet_jetconstituents/F");
        self.otree.Branch("gen_jet_jetconstituents", self.gen_jet_jetconstituents_ , "gen_jet_jetconstituents/F");

        self.otree.Branch("jet_rcore4", self.jet_rcore4_ , "jet_rcore4/F");
        self.otree.Branch("jet_rcore5", self.jet_rcore5_ , "jet_rcore5/F");
        self.otree.Branch("jet_rcore6", self.jet_rcore6_ , "jet_rcore6/F");
        self.otree.Branch("jet_rcore7", self.jet_rcore7_ , "jet_rcore7/F");

        self.otree.Branch("gen_jet_rcore4", self.gen_jet_rcore4_ , "gen_jet_rcore4/F");
        self.otree.Branch("gen_jet_rcore5", self.gen_jet_rcore5_ , "gen_jet_rcore5/F");
        self.otree.Branch("gen_jet_rcore6", self.gen_jet_rcore6_ , "gen_jet_rcore6/F");
        self.otree.Branch("gen_jet_rcore7", self.gen_jet_rcore7_ , "gen_jet_rcore7/F");

        self.otree.Branch("jet_pt1frac", self.jet_pt1frac_ , "jet_pt1frac/F");
        self.otree.Branch("jet_pt2frac", self.jet_pt2frac_ , "jet_pt2frac/F");
        self.otree.Branch("jet_sjdr", self.jet_sjdr_ , "jet_sjdr/F");


        self.otree.Branch("j_jecfactor_up", self.j_jecfactor_up_ , "j_jecfactor_up/F");
        self.otree.Branch("j_jecfactor_dn", self.j_jecfactor_dn_ , "j_jecfactor_dn/F");

        self.otree.Branch("jet_mass_pr_jes_up", self.jet_mass_pr_jes_up_ , "jet_mass_pr_jes_up/F");
        self.otree.Branch("jet_mass_pr_jes_dn", self.jet_mass_pr_jes_dn_ , "jet_mass_pr_jes_dn/F");

        self.otree.Branch("jet_mass_pr_jer", self.jet_mass_pr_jer_ , "jet_mass_pr_jer/F");
        self.otree.Branch("jet_mass_pr_jer_up", self.jet_mass_pr_jer_up_ , "jet_mass_pr_jer_up/F");
        self.otree.Branch("jet_mass_pr_jer_dn", self.jet_mass_pr_jer_dn_ , "jet_mass_pr_jer_dn/F");

        self.otree.Branch("ungroomed_jet_pt_jes_dn", self.ungroomed_jet_pt_jes_dn_ , "ungroomed_jet_pt_jes_dn/F");
        self.otree.Branch("ungroomed_jet_pt_jes_up", self.ungroomed_jet_pt_jes_up_ , "ungroomed_jet_pt_jes_up/F");

        self.otree.Branch("ungroomed_jet_pt_jer", self.ungroomed_jet_pt_jer_ , "ungroomed_jet_pt_jer/F");
        self.otree.Branch("ungroomed_jet_pt_jer_dn", self.ungroomed_jet_pt_jer_dn_ , "ungroomed_jet_pt_jer_dn/F");
        self.otree.Branch("ungroomed_jet_pt_jer_up", self.ungroomed_jet_pt_jer_up_ , "ungroomed_jet_pt_jer_up/F");

        self.otree.Branch("jet_planarlow04", self.jet_planarlow04_ , "jet_planarlow04/F");
        self.otree.Branch("jet_planarlow05", self.jet_planarlow05_ , "jet_planarlow05/F");
        self.otree.Branch("jet_planarlow06", self.jet_planarlow06_ , "jet_planarlow06/F");
        self.otree.Branch("jet_planarlow07", self.jet_planarlow07_ , "jet_planarlow07/F");

        ########## vbf jets 
        self.otree.Branch("vbf_maxpt_jj_m", self.vbf_maxpt_jj_m_ , "vbf_maxpt_jj_m/F");
        self.otree.Branch("vbf_maxpt_jj_pt", self.vbf_maxpt_jj_pt_ , "vbf_maxpt_jj_pt/F");
        self.otree.Branch("vbf_maxpt_jj_eta", self.vbf_maxpt_jj_eta_ , "vbf_maxpt_jj_eta/F");
        self.otree.Branch("vbf_maxpt_jj_phi", self.vbf_maxpt_jj_phi_ , "vbf_maxpt_jj_phi/F");

        self.otree.Branch("vbf_maxpt_j1_m", self.vbf_maxpt_j1_m_ , "vbf_maxpt_j1_m/F");
        self.otree.Branch("vbf_maxpt_j1_pt", self.vbf_maxpt_j1_pt_ , "vbf_maxpt_j1_pt/F");
        self.otree.Branch("vbf_maxpt_j1_eta", self.vbf_maxpt_j1_eta_ , "vbf_maxpt_j1_eta/F");
        self.otree.Branch("vbf_maxpt_j1_phi", self.vbf_maxpt_j1_phi_ , "vbf_maxpt_j1_phi/F");

        self.otree.Branch("vbf_maxpt_j2_m", self.vbf_maxpt_j2_m_ , "vbf_maxpt_j2_m/F");
        self.otree.Branch("vbf_maxpt_j2_pt", self.vbf_maxpt_j2_pt_ , "vbf_maxpt_j2_pt/F");
        self.otree.Branch("vbf_maxpt_j2_eta", self.vbf_maxpt_j2_eta_ , "vbf_maxpt_j2_eta/F");
        self.otree.Branch("vbf_maxpt_j2_phi", self.vbf_maxpt_j2_phi_ , "vbf_maxpt_j2_phi/F");

        self.otree.Branch("vbf_maxpt_jj_m_gen", self.vbf_maxpt_jj_m_gen_ , "vbf_maxpt_jj_m_gen/F");
        self.otree.Branch("vbf_maxpt_jj_pt_gen", self.vbf_maxpt_jj_pt_gen_ , "vbf_maxpt_jj_pt_gen/F");
        self.otree.Branch("vbf_maxpt_jj_eta_gen", self.vbf_maxpt_jj_eta_gen_ , "vbf_maxpt_jj_eta_gen/F");
        self.otree.Branch("vbf_maxpt_jj_phi_gen", self.vbf_maxpt_jj_phi_gen_ , "vbf_maxpt_jj_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j1_m_gen", self.vbf_maxpt_j1_m_gen_ , "vbf_maxpt_j1_m_gen/F");
        self.otree.Branch("vbf_maxpt_j1_pt_gen", self.vbf_maxpt_j1_pt_gen_ , "vbf_maxpt_j1_pt_gen/F");
        self.otree.Branch("vbf_maxpt_j1_eta_gen", self.vbf_maxpt_j1_eta_gen_ , "vbf_maxpt_j1_eta_gen/F");
        self.otree.Branch("vbf_maxpt_j1_phi_gen", self.vbf_maxpt_j1_phi_gen_ , "vbf_maxpt_j1_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j2_m_gen", self.vbf_maxpt_j2_m_gen_ , "vbf_maxpt_j2_m_gen/F");
        self.otree.Branch("vbf_maxpt_j2_pt_gen", self.vbf_maxpt_j2_pt_gen_ , "vbf_maxpt_j2_pt_gen/F");
        self.otree.Branch("vbf_maxpt_j2_eta_gen", self.vbf_maxpt_j2_eta_gen_ , "vbf_maxpt_j2_eta_gen/F");
        self.otree.Branch("vbf_maxpt_j2_phi_gen", self.vbf_maxpt_j2_phi_gen_ , "vbf_maxpt_j2_phi_gen/F");

        self.otree.Branch("vbf_maxpt_j1_m_jes_up", self.vbf_maxpt_j1_m_jes_up_ , "vbf_maxpt_j1_m_jes_up/F");
        self.otree.Branch("vbf_maxpt_j1_pt_jes_up", self.vbf_maxpt_j1_pt_jes_up_ , "vbf_maxpt_j1_pt_jes_up/F");
        self.otree.Branch("vbf_maxpt_j1_eta_jes_up", self.vbf_maxpt_j1_eta_jes_up_ , "vbf_maxpt_j1_eta_jes_up/F");
        self.otree.Branch("vbf_maxpt_j1_phi_jes_up", self.vbf_maxpt_j1_phi_jes_up_ , "vbf_maxpt_j1_phi_jes_up/F");

        self.otree.Branch("vbf_maxpt_j1_m_jes_dn", self.vbf_maxpt_j1_m_jes_dn_ , "vbf_maxpt_j1_m_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j1_pt_jes_dn", self.vbf_maxpt_j1_pt_jes_dn_ , "vbf_maxpt_j1_pt_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j1_eta_jes_dn", self.vbf_maxpt_j1_eta_jes_dn_ , "vbf_maxpt_j1_eta_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j1_phi_jes_dn", self.vbf_maxpt_j1_phi_jes_dn_ , "vbf_maxpt_j1_phi_jes_dn/F");

        self.otree.Branch("vbf_maxpt_j1_m_jer", self.vbf_maxpt_j1_m_jer_ , "vbf_maxpt_j1_m_jer/F");
        self.otree.Branch("vbf_maxpt_j1_pt_jer", self.vbf_maxpt_j1_pt_jer_ , "vbf_maxpt_j1_pt_jer/F");
        self.otree.Branch("vbf_maxpt_j1_eta_jer", self.vbf_maxpt_j1_eta_jer_ , "vbf_maxpt_j1_eta_jer/F");
        self.otree.Branch("vbf_maxpt_j1_phi_jer", self.vbf_maxpt_j1_phi_jer_ , "vbf_maxpt_j1_phi_jer/F");

        self.otree.Branch("vbf_maxpt_j1_m_jer_up", self.vbf_maxpt_j1_m_jer_up_ , "vbf_maxpt_j1_m_jer_up/F");
        self.otree.Branch("vbf_maxpt_j1_pt_jer_up", self.vbf_maxpt_j1_pt_jer_up_ , "vbf_maxpt_j1_pt_jer_up/F");
        self.otree.Branch("vbf_maxpt_j1_eta_jer_up", self.vbf_maxpt_j1_eta_jer_up_ , "vbf_maxpt_j1_eta_jer_up/F");
        self.otree.Branch("vbf_maxpt_j1_phi_jer_up", self.vbf_maxpt_j1_phi_jer_up_ , "vbf_maxpt_j1_phi_jer_up/F");

        self.otree.Branch("vbf_maxpt_j1_m_jer_dn", self.vbf_maxpt_j1_m_jer_dn_ , "vbf_maxpt_j1_m_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j1_pt_jer_dn", self.vbf_maxpt_j1_pt_jer_dn_ , "vbf_maxpt_j1_pt_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j1_eta_jer_dn", self.vbf_maxpt_j1_eta_jer_dn_ , "vbf_maxpt_j1_eta_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j1_phi_jer_dn", self.vbf_maxpt_j1_phi_jer_dn_ , "vbf_maxpt_j1_phi_jer_dn/F");

        self.otree.Branch("vbf_maxpt_j2_m_jes_up", self.vbf_maxpt_j2_m_jes_up_ , "vbf_maxpt_j2_m_jes_up/F");
        self.otree.Branch("vbf_maxpt_j2_pt_jes_up", self.vbf_maxpt_j2_pt_jes_up_ , "vbf_maxpt_j2_pt_jes_up/F");
        self.otree.Branch("vbf_maxpt_j2_eta_jes_up", self.vbf_maxpt_j2_eta_jes_up_ , "vbf_maxpt_j2_eta_jes_up/F");
        self.otree.Branch("vbf_maxpt_j2_phi_jes_up", self.vbf_maxpt_j2_phi_jes_up_ , "vbf_maxpt_j2_phi_jes_up/F");

        self.otree.Branch("vbf_maxpt_j2_m_jes_dn", self.vbf_maxpt_j2_m_jes_dn_ , "vbf_maxpt_j2_m_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j2_pt_jes_dn", self.vbf_maxpt_j2_pt_jes_dn_ , "vbf_maxpt_j2_pt_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j2_eta_jes_dn", self.vbf_maxpt_j2_eta_jes_dn_ , "vbf_maxpt_j2_eta_jes_dn/F");
        self.otree.Branch("vbf_maxpt_j2_phi_jes_dn", self.vbf_maxpt_j2_phi_jes_dn_ , "vbf_maxpt_j2_phi_jes_dn/F");

        self.otree.Branch("vbf_maxpt_j2_m_jer", self.vbf_maxpt_j2_m_jer_ , "vbf_maxpt_j2_m_jer/F");
        self.otree.Branch("vbf_maxpt_j2_pt_jer", self.vbf_maxpt_j2_pt_jer_ , "vbf_maxpt_j2_pt_jer/F");
        self.otree.Branch("vbf_maxpt_j2_eta_jer", self.vbf_maxpt_j2_eta_jer_ , "vbf_maxpt_j2_eta_jer/F");
        self.otree.Branch("vbf_maxpt_j2_phi_jer", self.vbf_maxpt_j2_phi_jer_ , "vbf_maxpt_j2_phi_jer/F");

        self.otree.Branch("vbf_maxpt_j2_m_jer_up", self.vbf_maxpt_j2_m_jer_up_ , "vbf_maxpt_j2_m_jer_up/F");
        self.otree.Branch("vbf_maxpt_j2_pt_jer_up", self.vbf_maxpt_j2_pt_jer_up_ , "vbf_maxpt_j2_pt_jer_up/F");
        self.otree.Branch("vbf_maxpt_j2_eta_jer_up", self.vbf_maxpt_j2_eta_jer_up_ , "vbf_maxpt_j2_eta_jer_up/F");
        self.otree.Branch("vbf_maxpt_j2_phi_jer_up", self.vbf_maxpt_j2_phi_jer_up_ , "vbf_maxpt_j2_phi_jer_up/F");

        self.otree.Branch("vbf_maxpt_j2_m_jer_dn", self.vbf_maxpt_j2_m_jer_dn_ , "vbf_maxpt_j2_m_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j2_pt_jer_dn", self.vbf_maxpt_j2_pt_jer_dn_ , "vbf_maxpt_j2_pt_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j2_eta_jer_dn", self.vbf_maxpt_j2_eta_jer_dn_ , "vbf_maxpt_j2_eta_jer_dn/F");
        self.otree.Branch("vbf_maxpt_j2_phi_jer_dn", self.vbf_maxpt_j2_phi_jer_dn_ , "vbf_maxpt_j2_phi_jer_dn/F");

        self.otree.Branch("vbf_maxpt_j1_QGLikelihood", self.vbf_maxpt_j1_QGLikelihood_ , "vbf_maxpt_j1_QGLikelihood/F");
        self.otree.Branch("vbf_maxpt_j2_QGLikelihood", self.vbf_maxpt_j2_QGLikelihood_ , "vbf_maxpt_j2_QGLikelihood/F");

        self.otree.Branch("vbf_maxpt_j1_isPileUpMedium", self.vbf_maxpt_j1_isPileUpMedium_ , "vbf_maxpt_j1_isPileUpMedium/I");
        self.otree.Branch("vbf_maxpt_j2_isPileUpMedium", self.vbf_maxpt_j2_isPileUpMedium_ , "vbf_maxpt_j2_isPileUpMedium/I");

        self.otree.Branch("vbf_maxpt_j1_isPileUpTight", self.vbf_maxpt_j1_isPileUpTight_ , "vbf_maxpt_j1_isPileUpTight/I");
        self.otree.Branch("vbf_maxpt_j2_isPileUpTight", self.vbf_maxpt_j2_isPileUpTight_ , "vbf_maxpt_j2_isPileUpTight/I");

        self.otree.Branch("vbf_maxpt_j1_bDiscriminatorCSV", self.vbf_maxpt_j1_bDiscriminatorCSV_ , "vbf_maxpt_j1_bDiscriminatorCSV/F");
        self.otree.Branch("vbf_maxpt_j2_bDiscriminatorCSV", self.vbf_maxpt_j2_bDiscriminatorCSV_ , "vbf_maxpt_j2_bDiscriminatorCSV/F");

        self.otree.Branch("vbf_maxpt_j1_bDiscriminatorCSV_gen", self.vbf_maxpt_j1_bDiscriminatorCSV_gen_ , "vbf_maxpt_j1_bDiscriminatorCSV_gen/F");
        self.otree.Branch("vbf_maxpt_j2_bDiscriminatorCSV_gen", self.vbf_maxpt_j2_bDiscriminatorCSV_gen_ , "vbf_maxpt_j2_bDiscriminatorCSV_gen/F");


        ###### btag counters
        self.otree.Branch("nbjets_csvl_veto", self.nbjets_csvl_veto_ , "nbjets_csvl_veto/F");
        self.otree.Branch("nbjets_csvm_veto", self.nbjets_csvm_veto_ , "nbjets_csvm_veto/F");
        self.otree.Branch("nbjets_csvt_veto", self.nbjets_csvt_veto_ , "nbjets_csvt_veto/F");
        self.otree.Branch("nbjets_ssvhem_veto", self.nbjets_ssvhem_veto_ , "nbjets_ssvhem_veto/F");

        self.otree.Branch("nbjets_csvl_veto_cleaned", self.nbjets_csvl_veto_cleaned_ , "nbjets_csvl_veto_cleaned/F");
        self.otree.Branch("nbjets_csvm_veto_cleaned", self.nbjets_csvm_veto_cleaned_ , "nbjets_csvm_veto_cleaned/F");
        self.otree.Branch("nbjets_csvt_veto_cleaned", self.nbjets_csvt_veto_cleaned_ , "nbjets_csvt_veto_cleaned/F");
        self.otree.Branch("nbjets_ssvhem_veto_cleaned", self.nbjets_ssvhem_veto_cleaned_ , "nbjets_ssvhem_veto_cleaned/F");

        self.otree.Branch("njets", self.njets_ , "njets/F");

        self.otree.Branch("deltaR_lca8jet",     self.deltaR_lca8jet_ , "deltaR_lca8jet/F")
        self.otree.Branch("deltaphi_METca8jet", self.deltaphi_METca8jet_ , "deltaphi_METca8jet/F")
        self.otree.Branch("deltaphi_Vca8jet",   self.deltaphi_Vca8jet_ , "deltaphi_Vca8jet/F")
        self.otree.Branch("deltaphi_METca8jet_met", self.deltaphi_METca8jet_met_ , "deltaphi_METca8jet_met/F")
        self.otree.Branch("deltaphi_Vca8jet_met", self.deltaphi_Vca8jet_met_ , "deltaphi_Vca8jet_met/F")

        ### some generator level quantites
        self.otree.Branch("genHMass",     self.genHMass_ , "genHMass/F")
        self.otree.Branch("genHphi",     self.genHphi_ , "genHphi/F")
        self.otree.Branch("genHeta",     self.genHeta_ , "genHeta/F")
        self.otree.Branch("genHpt",     self.genHpt_ , "genHpt/F")
                
        self.otree.Branch("genTagQuark1W",       self.genTagQuark1E_ , "genTagQuark1E/F")
        self.otree.Branch("genTagQuark1phi",     self.genTagQuark1phi_ , "genTagQuark1phi/F")
        self.otree.Branch("genTagQuark1eta",     self.genTagQuark1eta_ , "genTagQuark1eta/F")
        self.otree.Branch("genTagQuark1pt",      self.genTagQuark1pt_ , "genTagQuark1pt/F")

        self.otree.Branch("genTagQuarkE",     self.genTagQuark2E_ , "genTagQuark2E/F")
        self.otree.Branch("genTagQuark2phi",     self.genTagQuark2phi_ , "genTagQuark2phi/F")
        self.otree.Branch("genTagQuark2eta",     self.genTagQuark2eta_ , "genTagQuark2eta/F")
        self.otree.Branch("genTagQuark2pt",     self.genTagQuark2pt_ , "genTagQuark2pt/F")
                                                                        
        # variables for ttbar control region
        self.otree.Branch("ttb_nak5_same",     self.ttb_nak5_same_ , "ttb_nak5_same/F")
        self.otree.Branch("ttb_nak5_same_csvl",     self.ttb_nak5_same_csvl_ , "ttb_nak5_same_csvl/F")
        self.otree.Branch("ttb_nak5_same_csvm",     self.ttb_nak5_same_csvm_ , "ttb_nak5_same_csvm/F")
        self.otree.Branch("ttb_nak5_same_csvt",     self.ttb_nak5_same_csvt_ , "ttb_nak5_same_csvt/F")

        self.otree.Branch("ttb_nak5_oppo",     self.ttb_nak5_oppo_ , "ttb_nak5_oppo/F")
        self.otree.Branch("ttb_nak5_oppo_csvl",     self.ttb_nak5_oppo_csvl_ , "ttb_nak5_oppo_csvl/F")
        self.otree.Branch("ttb_nak5_oppo_csvm",     self.ttb_nak5_oppo_csvm_ , "ttb_nak5_oppo_csvm/F")
        self.otree.Branch("ttb_nak5_oppo_csvt",     self.ttb_nak5_oppo_csvt_ , "ttb_nak5_oppo_csvt/F")

        self.otree.Branch("ttb_nak5_oppoveto",     self.ttb_nak5_oppoveto_ , "ttb_nak5_oppoveto/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvl",     self.ttb_nak5_oppoveto_csvl_ , "ttb_nak5_oppoveto_csvl/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvm",     self.ttb_nak5_oppoveto_csvm_ , "ttb_nak5_oppoveto_csvm/F")
        self.otree.Branch("ttb_nak5_oppoveto_csvt",     self.ttb_nak5_oppoveto_csvt_ , "ttb_nak5_oppoveto_csvt/F")

        self.otree.Branch("ttb_nak5_sameveto",     self.ttb_nak5_sameveto_ , "ttb_nak5_sameveto/F")
        self.otree.Branch("ttb_nak5_sameveto_csvl",     self.ttb_nak5_sameveto_csvl_ , "ttb_nak5_sameveto_csvl/F")
        self.otree.Branch("ttb_nak5_sameveto_csvm",     self.ttb_nak5_sameveto_csvm_ , "ttb_nak5_sameveto_csvm/F")
        self.otree.Branch("ttb_nak5_sameveto_csvt",     self.ttb_nak5_sameveto_csvt_ , "ttb_nak5_sameveto_csvt/F")

        self.otree.Branch("ttb_ht",     self.ttb_ht_ , "ttb_ht/F")
        self.otree.Branch("ttb_ca8_mass_pr",     self.ttb_ca8_mass_pr_ , "ttb_ca8_mass_pr/F")
        self.otree.Branch("ttb_ca8_charge",     self.ttb_ca8_charge_ , "ttb_ca8_charge/F")
        self.otree.Branch("ttb_ca8_charge_k05",     self.ttb_ca8_charge_k05_ , "ttb_ca8_charge_k05/F")
        self.otree.Branch("ttb_ca8_charge_k07",     self.ttb_ca8_charge_k07_ , "ttb_ca8_charge_k07/F")
        self.otree.Branch("ttb_ca8_charge_k10",     self.ttb_ca8_charge_k10_ , "ttb_ca8_charge_k10/F")

        self.otree.Branch("ttb_ca8_ungroomed_pt",     self.ttb_ca8_ungroomed_pt_ , "ttb_ca8_ungroomed_pt/F")
        self.otree.Branch("ttb_ca8_ungroomed_eta",     self.ttb_ca8_ungroomed_eta_ , "ttb_ca8_ungroomed_eta/F")
        self.otree.Branch("ttb_ca8_ungroomed_phi",     self.ttb_ca8_ungroomed_phi_ , "ttb_ca8_ungroomed_phi/F")
        self.otree.Branch("ttb_ca8_ungroomed_e",     self.ttb_ca8_ungroomed_e_ , "ttb_ca8_ungroomed_e/F")
        

        self.otree.Branch("ttb_ca8_ungroomed_gen_pt",     self.ttb_ca8_ungroomed_gen_pt_ , "ttb_ca8_ungroomed_gen_pt/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_eta",     self.ttb_ca8_ungroomed_gen_eta_ , "ttb_ca8_ungroomed_gen_eta/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_phi",     self.ttb_ca8_ungroomed_gen_phi_ , "ttb_ca8_ungroomed_gen_phi/F")
        self.otree.Branch("ttb_ca8_ungroomed_gen_e",     self.ttb_ca8_ungroomed_gen_e_ , "ttb_ca8_ungroomed_gen_e/F")
        
        self.otree.Branch("ttb_ca8_tau2tau1",     self.ttb_ca8_tau2tau1_ , "ttb_ca8_tau2tau1/F")
        self.otree.Branch("ttb_ca8_tau2tau1_exkT",     self.ttb_ca8_tau2tau1_exkT_ , "ttb_ca8_tau2tau1_exkT/F")
        self.otree.Branch("ttb_ca8_tau2tau1_pr",     self.ttb_ca8_tau2tau1_pr_ , "ttb_ca8_tau2tau1_pr/F")
        self.otree.Branch("ttb_ca8_GeneralizedECF",     self.ttb_ca8_GeneralizedECF_ , "ttb_ca8_GeneralizedECF/F")
        self.otree.Branch("ttb_ca8_mu",     self.ttb_ca8_mu_ , "ttb_ca8_mu/F")

        self.otree.Branch("ttb_ca8_mlvj_type0",     self.ttb_ca8_mlvj_type0_ , "ttb_ca8_mlvj_type0/F")
        self.otree.Branch("ttb_ca8_mlvj_type2",     self.ttb_ca8_mlvj_type2_ , "ttb_ca8_mlvj_type2/F")

        self.otree.Branch("ttb_ca8_mlvj_type0_met",     self.ttb_ca8_mlvj_type0_met_ , "ttb_ca8_mlvj_type0_met/F")
        self.otree.Branch("ttb_ca8_mlvj_type2_met",     self.ttb_ca8_mlvj_type2_met_ , "ttb_ca8_mlvj_type2_met/F")

        self.otree.Branch("ttb_dR_ca8_bjet_closer",     self.ttb_dR_ca8_bjet_closer_ , "ttb_dR_ca8_bjet_closer/F")
        self.otree.Branch("ttb_dR_ca8_jet_closer",     self.ttb_dR_ca8_jet_closer_ , "ttb_dR_ca8_jet_closer/F")

        self.otree.Branch("isttbar",     self.isttbar_ , "isttbar/I")

        self.otree.Branch("gen_parton1_px_fromttbar",     self.gen_parton1_px_fromttbar_ , "gen_parton1_px_fromttbar/F")
        self.otree.Branch("gen_parton1_py_fromttbar",     self.gen_parton1_py_fromttbar_ , "gen_parton1_py_fromttbar/F")
        self.otree.Branch("gen_parton1_pz_fromttbar",     self.gen_parton1_pz_fromttbar_ , "gen_parton1_pz_fromttbar/F")
        self.otree.Branch("gen_parton1_e_fromttbar",     self.gen_parton1_e_fromttbar_ , "gen_parton1_e_fromttbar/F")
        self.otree.Branch("gen_parton1_id_fromttbar",     self.gen_parton1_id_fromttbar_ , "gen_parton1_id_fromttbar/F")

        self.otree.Branch("gen_parton2_px_fromttbar",     self.gen_parton2_px_fromttbar_ , "gen_parton2_px_fromttbar/F")
        self.otree.Branch("gen_parton2_py_fromttbar",     self.gen_parton2_py_fromttbar_ , "gen_parton2_py_fromttbar/F")
        self.otree.Branch("gen_parton2_pz_fromttbar",     self.gen_parton2_pz_fromttbar_ , "gen_parton2_pz_fromttbar/F")
        self.otree.Branch("gen_parton2_e_fromttbar",     self.gen_parton2_e_fromttbar_ , "gen_parton2_e_fromttbar/F")
        self.otree.Branch("gen_parton2_id_fromttbar",     self.gen_parton2_id_fromttbar_ , "gen_parton2_id_fromttbar/F")

