#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TString.h"
#include "TGaxis.h"
#include "TGraph2D.h"
#include "TFile.h"
#include "TF1.h"

#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooGaussian.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooVoigtian.h"
#include "RooCBShape.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooArgusBG.h"
#include "RooWorkspace.h"
#include "RooChebychev.h"
#include "RooBernstein.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooChi2Var.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooProdPdf.h"

#include "../PDFs/MakePdf.h"
#include "../PlotStyle/PlotUtils.h"
#include "../PDFs/PdfDiagonalizer.h"

void fit_genHMass(RooWorkspace* , const std::string & = "", const std::string & = "", const std::string & = "", const std::string & = "",const std::string & = "", const int & = 0);

void fit_mj_single_MC(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "",const std::string & = "em", const std::string & = "HP", const int & = 0, const int & = 0, const int & = 0);

void fit_mlvj_model_single_MC(RooWorkspace*, const std::string & ="", const std::string & = " ", const std::string & = "_signal_region", const std::string & = "", const std::string & = "em", const std::string & = "HP", const int & = 0, const int & = 0, const int & = 1, const int & = 0, const std::string & = "", const std::string & = "");

void fit_WJetsNormalization_in_Mj_signal_region(RooWorkspace*,  std::map<std::string,int>, std::map<std::string,std::string>, const std::string & = "", const std::string & = "", const std::string & = "", const std::string & = "", const std::string = "", const int  & = 0 , const int & = 0, const float & = 65, const float & = 105, const std::string & = "");

void get_mj_normalization_insignalregion(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "");

void fit_mlvj_in_Mj_sideband(RooWorkspace*, std::map<std::string,int>, std::map<std::string,std::string>, const std::string & = "", const std::string & = "", const std::string & = "_sb_lo", const std::string & = "em", const std::string & = "HP", const std::string & = "", const int & = 0, const int & = 1, const int & = 0,const std::string & = "", const std::string & = "");

void get_WJets_mlvj_correction_sb_lo_to_signal_region(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "", const std::string & = "_sb_lo", const std::string & = "4fit_", const std::string & = "em", const std::string & = "HP", const int & = 1);

void get_mlvj_normalization_insignalregion(RooWorkspace* ,const std::string & = "", const std::string & = "", const std::string & = "", const std::string & = "", const int & = 0);

void ScaleFactorTTbarControlSampleFit(RooWorkspace*, std::map<std::string,std::string >, std::map<std::string,int>, std::vector<std::string>* = NULL, std::vector<std::string>* = NULL, const std::string & ="", const std::string & ="mu", const std::string & wtagger ="HP", const double & = 200, const double & = 2000);

void DrawScaleFactorTTbarControlSample(RooWorkspace*,  std::map<std::string,int>,  const std::string & ="", const std::string & ="mu", const std::string & ="HP",const double & = 200, const double & = 2000);
					
///////////////////

double crystalBallLowHigh (double* x, double* par);

double doubleGausCrystalBallLowHighPlusExp (double* x, double* par);

Double_t CrystalBallLowHighPlusExpDividedByCrystalBallLowHigh(Double_t *x,Double_t *par);

Double_t getIntWght(std::string wFile , double realMass, double Hmass = 350, double cprime = 1.0, double BRnew = 0.0);
