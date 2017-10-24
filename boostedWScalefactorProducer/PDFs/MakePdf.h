#include <iostream>
#include <string>
#include <vector>

#include "TString.h"

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

#include "HWWLVJRooPdfs.h"

RooAbsPdf* MakeGeneralPdf(RooWorkspace* ,const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & = "HP",  std::vector<std::string>* = NULL, const int & = 0);

RooExtendPdf* MakeExtendedModel(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & wtagger_label = "HP",  std::vector<std::string>* = NULL, const int & = 0, const int & = 500);

RooGaussian* addConstraint(RooRealVar*, double, double,  std::vector<std::string>*);

void change_dataset_to_histpdf(RooWorkspace*,RooRealVar*,RooDataSet*);

RooAbsPdf* MakeModelTTbarControlSample(RooWorkspace* workspace ,const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & = "HP", const std::string & = "", std::vector<std::string>* = NULL);

RooAbsPdf* get_mj_Model        (RooWorkspace*,const std::string & = "",const std::string & = "em");
RooAbsPdf* get_General_mj_Model(RooWorkspace*,const std::string & = "",const std::string & = "",const std::string & = "em", const int & = 1);