#include <iostream>
#include <algorithm>
#include <vector>
#include <string>

#include "TH1F.h"
#include "TChain.h"
#include "TString.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TIterator.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "Math/QuantFuncMathCore.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"

#include "RooPlot.h"
#include "RooHist.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooCurve.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
#include "RooExtendPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"

using namespace std;
using namespace ROOT;
using namespace TMath;

//////////////////////////////////////////////////////////////

void draw_error_band(RooAbsData*, RooAbsPdf*,  RooRealVar*, RooFitResult*, RooPlot*, const int & = 6, const std::string & ="F", const int & = 100, const int & = 2000);

void draw_error_band_extendPdf(RooAbsData *, RooAbsPdf*, RooFitResult*, RooPlot*, const int & = 6, const std::string & = "F", const int & = 100,  const int & = 2000);

void draw_error_band_Decor( std::string, std::string, RooArgList *,  RooWorkspace *, RooRealVar *, RooPlot *, const int & = 6, const std::string & = "F", const int & = 100, const int & = 2000);

void draw_error_band_shape_Decor(std::string, std::string,  RooArgList *,  RooWorkspace *, double , RooPlot *, const int & = 6, const std::string & = "F", const int & = 3013, const std::string & = "", const int & = 100, const int & = 2000);

std::vector<TGraphAsymmErrors*>* draw_error_band_ws(RooAbsData*, RooAbsPdf *, const std::string & = "", RooRealVar * = NULL,  RooArgList* = NULL, RooWorkspace * = NULL, RooPlot * = NULL,const int & = 6, const std::string & = "F", const int & = 100, const int & = 1, const int & = 2000, const int & = 0);


//////////////////////////////////////////////////////////////

void draw_error_band_pull(RooAbsData *, RooAbsPdf *, RooRealVar *, RooFitResult*, RooPlot*, const int & = 6, const int & = 100, const int & = 2000);

void draw_error_band_ratio(RooAbsData *,RooAbsPdf *, RooRealVar *, RooFitResult*, RooPlot*, const int & = 6, const int & = 100, const int & = 2000);

void draw_error_band_extendPdf_pull(RooAbsData*, RooAbsPdf*, RooFitResult*, RooPlot*, const int & = 6, const int & = 100, const int & = 2000);

void draw_error_band_extendPdf_ratio(RooAbsData*, RooAbsPdf*, RooFitResult*, RooPlot*, const int & = 6, const int & = 100, const int & = 2000);

//////////////////////////////////////////////////////////////

double Calc_error_extendPdf(RooAbsData*,RooExtendPdf*,RooFitResult*,const std::string & = "",const int & = 2000);

double Calc_error(std::string,std::string,RooArgList*,RooWorkspace*,const std::string & = "",const int & = 2000);
