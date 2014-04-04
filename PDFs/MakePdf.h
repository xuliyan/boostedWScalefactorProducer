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

struct SystematicUncertaintyEXO{ const double mean_signal_uncertainty_jet_scale = 0.013;
                                 const double mean_signal_uncertainty_lep_scale = 0.001;
                                 const double sigma_signal_uncertainty_jet_scale = 0.033;
                                 const double sigma_signal_uncertainty_jet_res = 0.030;
                                 const double sigma_signal_uncertainty_lep_scale = 0.005;
                                 const double mean_mj_shift = 1.4 ;
                                 const double sigma_mj_smear = 1.11 ;
};       

struct SystematicUncertaintyHiggs{const double mean_signal_uncertainty_jet_scale = 0.013;
                                  const double mean_signal_uncertainty_lep_scale = 0.001;
                                  const double sigma_signal_uncertainty_jet_scale = 0.033;
                                  const double sigma_signal_uncertainty_jet_res = 0.030;
                                  const double sigma_signal_uncertainty_lep_scale = 0.005;
                                  const double mean_mj_shift = 1.4 ;
                                  const double sigma_mj_smear = 1.11 ;
};       

RooAbsPdf* MakeGeneralPdf(RooWorkspace* ,const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & = "HP", RooArgList* = NULL, const int & = 0);

RooExtendPdf* MakeExtendedModel(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & wtagger_label = "HP", RooArgList* = NULL, const int & = 0, const int & = 500);

RooGaussian* addConstraint(RooRealVar*, RooRealVar*, RooRealVar*, std::vector<std::string> &);

void change_dataset_to_histpdf(RooWorkspace*,RooRealVar*,RooDataSet*);

TH1F* change_dataset_to_histogram(RooRealVar*, RooDataSet*, const std::string & = "", const int & = 1);

//////////

RooAbsPdf* get_mj_Model        (RooWorkspace*,const std::string & = "",const std::string & = "em");
RooAbsPdf* get_General_mj_Model(RooWorkspace*,const std::string & = "",const std::string & = "",const std::string & = "em", const int & = 1);
RooAbsPdf* get_TTbar_mj_Model  (RooWorkspace*,const std::string & = "_TTbar",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_STop_mj_Model   (RooWorkspace*,const std::string & = "_STop",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_VV_mj_Model     (RooWorkspace*,const std::string & = "_VV",const std::string & model = "", const std::string & = "em", const int & = 1);
RooAbsPdf* get_WW_EWK_mj_Model (RooWorkspace*,const std::string & = "_WW_EWK",const std::string & model = "", const std::string & = "em", const int & = 1);
RooAbsPdf* get_WJets_mj_Model  (RooWorkspace*,const std::string & = "_WJets0",const std::string & model = "",const std::string & = "em", const int & = 1, const std::string & = "");
RooAbsPdf* get_ggH_mj_Model    (RooWorkspace*,const std::string & = "_ggH600",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_vbfH_mj_Model   (RooWorkspace*,const std::string & = "_vbfH600",const std::string & model = "",const std::string & = "em",const int & = 1);


RooAbsPdf* get_mlvj_Model        (RooWorkspace*,const std::string & = "",const std::string & = "",const std::string & = "", const std::string & = "em");
RooAbsPdf* get_General_mlvj_Model(RooWorkspace*,const std::string & = "",const std::string & = "",const std::string & = "",const std::string & = "em", const int & = 1);
RooAbsPdf* get_TTbar_mlvj_Model  (RooWorkspace*,const std::string & = "_TTbar",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_STop_mlvj_Model   (RooWorkspace*,const std::string & = "_STop",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_VV_mlvj_Model     (RooWorkspace*,const std::string & = "_VV",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_WW_EWK_mlvj_Model (RooWorkspace*,const std::string & = "_WW_EWK",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_WJets_mlvj_Model  (RooWorkspace*, const std::string & = "_WJets0",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_ggH_mlvj_Model    (RooWorkspace*, const std::string & = "_ggH600",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);
RooAbsPdf* get_vbfH_mlvj_Model   (RooWorkspace*, const std::string & = "_vbfH600",const std::string & = "",const std::string & model = "",const std::string & = "em",const int & = 1);

///////////////////////
void fix_Model(RooWorkspace*, const std::string & = "", const std::string & = "_signal_region", const std::string & = "_mlvj", const std::string & = "",const std::string & = "em", const std::string & = "", const int & = 0);

void fix_Pdf (RooAbsPdf*, RooArgSet*);

void ShowParam_Pdf(RooAbsPdf*,RooArgSet*);

////////////////////////
void clone_Model(RooWorkspace*, RooAbsPdf*, const std::string & = "", const std::string & = "_signal_region", const std::string & = "_mlvj", const std::string & = "", const std::string & = "em", const int & = 1);
