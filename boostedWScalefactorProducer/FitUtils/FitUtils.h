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
// #include "../PDFs/PdfDiagonalizer.h"


void fit_mj_single_MC(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "",const std::string & = "em", const std::string & = "HP");

void ScaleFactorTTbarControlSampleFit(RooWorkspace*, std::map<std::string,std::string >, std::map<std::string,int>, std::vector<std::string>* = NULL, std::vector<std::string>* = NULL, const std::string & ="mu", const std::string & wtagger ="HP");

void DrawScaleFactorTTbarControlSample(RooWorkspace*,  std::map<std::string,int>,  const std::string & ="", const std::string & ="mu", const std::string & ="HP",const double & = 200, const double & = 2000,   const std::string & ="");