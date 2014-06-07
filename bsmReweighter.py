#! /Usr/bin/env python
import os
import glob
import math
import array
import ROOT
import subprocess
from array import array
from subprocess import Popen

from ROOT import gROOT, gStyle, gSystem, TLatex, TString, TF1,TFile,TLine, TLegend, TH1D,TH2D,THStack,TChain, TCanvas, TMatrixDSym, TMath, TText, TPad, RooFit, RooArgSet, RooArgList, RooArgSet, RooAbsData, RooAbsPdf, RooAddPdf, RooWorkspace, RooExtendPdf,RooCBShape, RooLandau, RooFFTConvPdf, RooGaussian, RooBifurGauss, RooDataSet, RooExponential,RooBreitWigner, RooVoigtian, RooNovosibirsk, RooRealVar,RooFormulaVar, RooDataHist, RooHistPdf,RooCategory, RooSimultaneous, RooGenericPdf, RooKeysPdf, RooHistPdf, RooEffProd, RooProdPdf, TIter, kTRUE, kFALSE, kGray, kRed, kDashed, kGreen,kAzure, kOrange, kBlack,kBlue,kYellow,kCyan


from ROOT import RooTrace

ROOT.gStyle.SetPadTopMargin(0.09);
ROOT.gStyle.SetPadLeftMargin(0.16);
ROOT.gStyle.SetPalette(1);

ROOT.gSystem.Load("PDFs/RooRelBWRunningWidth_cxx.so")

############################################################

def FitMassPoint(filename, massin, massmin, massmax, rwCPS=1, nbins=50):

    ### take the input tree from the file
    inputFile = ROOT.TFile(filename);
    tree      = inputFile.Get("WJet");
    
    print "Lineshape Higgs mass -->  n entries: ", tree.GetEntries();
    
    # RooFitting
    rrv_mass   = RooRealVar("rrv_mass","rrv_mass",massin,massmin,massmax);
    rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.); 

    rrv_mH1    = RooRealVar("rrv_mH1","rrv_mH1", massin, massmin, massmax );
    rrv_gamma1 = RooRealVar("rrv_gamma1","rrv_gamma1",20.,0.,massmax);

    rrv_mH2    = RooRealVar("rrv_mH2","rrv_mH2", massin, massmin, massmax )
    rrv_gamma2 = RooRealVar("rrv_gamma2","rrv_gamma2",20.,0.,massmax)

    rrv_mH3    = RooRealVar("rrv_mH3","rrv_mH3", massin, massmin, massmax );
    rrv_gamma3 = RooRealVar("rrv_gamma3","rrv_gamma3",20.,0.,massmax);

    rds_raw      = RooDataSet("rds_raw","rds_raw",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## created the raw dataset --> raw lineshape
    rds_cps      = RooDataSet("rds_cps","rds_cps",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## create the cps dataset  --> cps lineshape
    rds_cps_intf = ROOT.RooDataSet("rds_cps_intf","rds_cps_intf",RooArgSet(rrv_mass,rrv_weight),RooFit.WeightVar(rrv_weight)); ## create the cps + interference dataset

    model1_pdf = ROOT.RooRelBWRunningWidth("model1_pdf","model1_pdf",rrv_mass,rrv_mH1,rrv_gamma1); ## BWRunningWidth Pdf from the dedicated library
    model2_pdf = ROOT.RooRelBWRunningWidth("model2_pdf","model2_pdf",rrv_mass,rrv_mH2,rrv_gamma2); ## BWRunningWidth Pdf from the dedicated library
    model3_pdf = ROOT.RooRelBWRunningWidth("model3_pdf","model3_pdf",rrv_mass,rrv_mH3,rrv_gamma3); ## BWRunningWidth Pdf from the dedicated library
    
    ### loop on Higgs signal events 
    for i in range(tree.GetEntries()):
        
        if i % 10000 == 0: print "reweighter, i: ", i;
        tree.GetEntry(i);
        
        curmass = getattr(tree,"W_H_mass_gen"); ## take the higgs generated mass inside the window; point not in the window are neglected
        if curmass < massmax and curmass > massmin:
            
            rrv_mass.setVal( curmass ); ## set the value of the RooRealVar            
            tmpweight_cps = getattr(tree,"complexpolewtggH"+str(massin))*rwCPS/getattr(tree,"avecomplexpolewtggH"+str(massin)); ## take the cps weight from the tree
            tmpweight_cps_intf = getattr(tree,"complexpolewtggH"+str(massin))*rwCPS*getattr(tree,"interferencewtggH"+str(massin))/getattr(tree,"avecomplexpolewtggH"+str(massin)); ## cps*int

            rds_raw.add( RooArgSet( rrv_mass ), 1. );
            rds_cps.add( RooArgSet( rrv_mass ), tmpweight_cps );
            rds_cps_intf.add( RooArgSet( rrv_mass ), tmpweight_cps_intf );


    print ">>>>"
    RooTrace.dump(ROOT.cout,ROOT.kTRUE);
    RooTrace.mark();                
    print "<<<<"

    #### final fit  

    model1_pdf.fitTo(rds_raw,RooFit.Save(1), RooFit.SumW2Error(kTRUE));  
    model1_pdf.fitTo(rds_raw,RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));  

    model2_pdf.fitTo(rds_cps,RooFit.Save(1), RooFit.SumW2Error(kTRUE));  
    model2_pdf.fitTo(rds_cps,RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));  

    model3_pdf.fitTo(rds_cps_intf,RooFit.Save(1), RooFit.SumW2Error(kTRUE));  
    model3_pdf.fitTo(rds_cps_intf,RooFit.Save(1), RooFit.SumW2Error(kTRUE), RooFit.Minimizer("Minuit2"));  

    ## plot of the fits
    mplot = rrv_mass.frame(RooFit.Title("mass plot"));
    rds_raw.plotOn(mplot, RooFit.MarkerColor(kBlack), RooFit.LineColor(kBlack), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps.plotOn(mplot, RooFit.MarkerColor(kRed), RooFit.LineColor(kRed), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps_intf.plotOn(mplot, RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    model1_pdf.plotOn(mplot, RooFit.LineColor(kBlack));
    model2_pdf.plotOn(mplot, RooFit.LineColor(kRed), RooFit.LineStyle(2) );
    model3_pdf.plotOn(mplot, RooFit.LineColor(kBlue), RooFit.LineStyle(3) );
    rds_raw.plotOn(mplot, RooFit.MarkerColor(kBlack), RooFit.LineColor(kBlack), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps.plotOn(mplot, RooFit.MarkerColor(kRed), RooFit.LineColor(kRed), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );
    rds_cps_intf.plotOn(mplot, RooFit.MarkerColor(kBlue), RooFit.LineColor(kBlue), RooFit.Binning(nbins,massmin,massmax), RooFit.DataError(RooAbsData.SumW2) );

    print "rds_raw.sumEntries() = ", rds_raw.sumEntries()
    print "model1_pdf: mH = ", rrv_mH1.getVal(), ", gamma = ", rrv_gamma1.getVal();    
    print "rds_cps.sumEntries() = ", rds_cps.sumEntries()
    print "model2_pdf: mH = ", rrv_mH2.getVal(), ", gamma = ", rrv_gamma2.getVal();    
    print "rds_cps_intf.sumEntries() = ", rds_cps_intf.sumEntries()
    print "model3_pdf: mH = ", rrv_mH3.getVal(), ", gamma = ", rrv_gamma3.getVal();    

    dummy_h1 = ROOT.TH1F("dummy_h1","dummy_h1",1,0,1); 
    dummy_h1.SetMarkerColor( ROOT.kBlack );
    dummy_h2 = ROOT.TH1F("dummy_h2","dummy_h2",1,0,1); 
    dummy_h2.SetMarkerColor( ROOT.kRed );
    dummy_h3 = ROOT.TH1F("dummy_h3","dummy_h3",1,0,1); 
    dummy_h3.SetMarkerColor( ROOT.kBlue );


    L = TLegend(0.65,0.60,0.93,0.85);
    L.SetFillStyle(0);
    L.AddEntry(dummy_h1,"Powheg","p");
    L.AddEntry(dummy_h2,"w/CPS weight","p");
    L.AddEntry(dummy_h3,"w/CPS,Intf weight","p");

    can2 = ROOT.TCanvas("can2","can2",800,800);
    mplot.Draw();
    L.Draw();
    
    os.system("mkdir -p massFits");
    can2.SaveAs("massFits/mass_rf_"+str(massin)+".pdf");
    can2.SaveAs("massFits/mass_rf_"+str(massin)+".png");

    outputpar_1 = [];
    outputpar_1.append( rrv_mH1.getVal() );
    outputpar_1.append( rrv_gamma1.getVal() );
    outputpar_2 = [];
    outputpar_2.append( rrv_mH2.getVal() );
    outputpar_2.append( rrv_gamma2.getVal() );
    outputpar_3 = [];
    outputpar_3.append( rrv_mH3.getVal() );
    outputpar_3.append( rrv_gamma3.getVal() );


    model1_pdf.Delete();
    model2_pdf.Delete();
    model3_pdf.Delete();

    rds_raw.Delete(),
    rds_cps.Delete(),
    rds_cps_intf.Delete();

    rrv_mH1.Delete(),
    rrv_gamma1.Delete();
    rrv_mH2.Delete(),
    rrv_gamma2.Delete();
    rrv_mH3.Delete(),
    rrv_gamma3.Delete();

    rrv_mass.Delete();

    return outputpar_2

### definition of the BW relativistic --> given nominal higgs mass mH and Gamma 
def MyRunningWidthBW(mww, mH, gamma):

    numerator = (mww*mww*gamma/mH) ;
    denominator = (mww*mww-mH*mH)*(mww*mww-mH*mH) + (mww*mww*gamma/mH)*(mww*mww*gamma/mH);
    pdf = numerator/denominator;

    return pdf;

#### function for extract the reweight of each event according to the ratio of the two lineshape: the original one (denominator) and the rescaled one (numerator) 
def lineshapeWidthReweight(point, meanSM, gammaSM, Cprime, massmin, massmax):

    ## in the lineshape rescale, the mean value stay the same, the new Gamma is = Cprime/(1.-BRnew) --> in this way we can obtain the new lineshape with a simple re-weighting 
    weight = MyRunningWidthBW(point,meanSM,gammaSM*Cprime)/MyRunningWidthBW(point,meanSM,gammaSM);
    return weight;

#### used only for ggH samples ; input are the current value for the interference, the new lineshape parameters c' and BRnew
def IntfRescale(curIntfRw,cPrime,BRnew):

     curIoverS = curIntfRw - 1;
     newWeight = 1 + ((1-BRnew)*curIoverS)/cPrime; ## scale as the inverse of the xs ; ratio of the two is the new weight to go to the new lineshape for ggH
     if newWeight < 0: newWeight = 0.0;
     ratio = newWeight/curIntfRw;

     return ratio;
