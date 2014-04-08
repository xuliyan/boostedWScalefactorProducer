#include <iostream>
#include <string>
#include <vector>

#include "TString.h"
#include "TGaxis.h"

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

void fit_mj_single_MC(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "",
const std::string & = "em", const std::string & = "HP", const int & = 0);

void fit_mlvj_model_single_MC(RooWorkspace*, const std::string & ="", const std::string & = " ", const std::string & = "_signal_region", const std::string & = "", const std::string & = "em", const std::string & = "HP", const int & = 0, const int & = 0, const int & = 1, const int & = 0, const std::string & = "");

void fit_WJetsNormalization_in_Mj_signal_region(RooWorkspace*,  std::map<std::string,int>, std::map<std::string,std::string>, const std::string & = "", const std::string & = "", const std::string & = "", const std::string & = "", const std::string = "", const int  & = 0 , const int & = 0, const float & = 65, const float & = 105, const std::string & = "");

void get_mj_normalization_insignalregion(RooWorkspace*, const std::string & = "", const std::string & = "");

void fit_mlvj_in_Mj_sideband(RooWorkspace*, std::map<std::string,int>, const std::string & = "", const std::string & = "_sb_lo", const std::string & = "", const std::string & = "em", const std::string & = "HP", const std::string & = "", const int & = 0, const int & = 1, const int & = 0,const std::string & = "");

void get_WJets_mlvj_correction_sb_lo_to_signal_region(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "_sb_lo", const std::string & = "4fit_", const std::string & = "em", const std::string & = "HP", const int & = 1);
