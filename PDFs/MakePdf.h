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

RooAbsPdf* MakeGeneralPdf(RooWorkspace* ,const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & = "HP", RooArgList* = NULL, const int & = 0);

RooExtendPdf* MakeExtendedModel(RooWorkspace*, const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & wtagger_label = "HP", RooArgList* = NULL, const int & = 0, const int & = 500);

RooGaussian* addConstraint(RooRealVar*, RooRealVar*, RooRealVar*, std::vector<std::string> &);

void change_dataset_to_histpdf(RooWorkspace*,RooRealVar*,RooDataSet*);

TH1F* change_dataset_to_histogram(RooRealVar*, RooDataSet*, const std::string & = "", const int & = 1);

RooAbsPdf* MakeModelTTbarControlSample(RooWorkspace* workspace ,const std::string & = "", const std::string & = "", const std::string & = "_mj", const std::string & = "em", const std::string & = "HP", const std::string & = "", RooArgList* = NULL);

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

struct SystematicUncertaintyEXO{ const double mean_signal_uncertainty_jet_scale  = 0.013;
                                 const double mean_signal_uncertainty_lep_scale  = 0.001;
                                 const double sigma_signal_uncertainty_jet_scale = 0.033;
                                 const double sigma_signal_uncertainty_jet_res   = 0.030;
                                 const double sigma_signal_uncertainty_lep_scale = 0.005;
                                 const double mean_mj_shift = 1.4 ;
                                 const double sigma_mj_smear = 1.16 ;
};       

struct SystematicUncertaintyHiggs_01jetBin{

                                  /////// effect on the mean --> from exo numbers
                                  const double mean_signal_uncertainty_jet_scale_ggH_600  = 0.013;
                                  const double mean_signal_uncertainty_jet_scale_ggH_700  = 0.013;
                                  const double mean_signal_uncertainty_jet_scale_ggH_800  = 0.013;
                                  const double mean_signal_uncertainty_jet_scale_ggH_900  = 0.013;
                                  const double mean_signal_uncertainty_jet_scale_ggH_1000 = 0.013;

                                  const double mean_signal_uncertainty_jet_scale_vbfH_600  = 0.011;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_700  = 0.009;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_800  = 0.013;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_900  = 0.012;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_1000 = 0.001;

                                  const double mean_signal_uncertainty_jet_res_ggH_600  = 0.005;
                                  const double mean_signal_uncertainty_jet_res_ggH_700  = 0.005;
                                  const double mean_signal_uncertainty_jet_res_ggH_800  = 0.005;
                                  const double mean_signal_uncertainty_jet_res_ggH_900  = 0.005;
                                  const double mean_signal_uncertainty_jet_res_ggH_1000 = 0.005;

                                  const double mean_signal_uncertainty_jet_res_vbfH_600  = 0.006;
                                  const double mean_signal_uncertainty_jet_res_vbfH_700  = 0.005;
                                  const double mean_signal_uncertainty_jet_res_vbfH_800  = 0.006;
                                  const double mean_signal_uncertainty_jet_res_vbfH_900  = 0.012;
                                  const double mean_signal_uncertainty_jet_res_vbfH_1000 = 0.014;


                                  /////// effect on the sigma 
                                  const double sigma_signal_uncertainty_jet_scale_ggH_600  = 0.030;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_700  = 0.030;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_800  = 0.030;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_900  = 0.030;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_1000 = 0.030;

                                  const double sigma_signal_uncertainty_jet_scale_vbfH_600  = 0.051;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_700  = 0.022;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_800  = 0.030;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_900  = 0.057;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_1000 = 0.015;

                                  const double sigma_signal_uncertainty_jet_res_ggH_600  = 0.050;
                                  const double sigma_signal_uncertainty_jet_res_ggH_700  = 0.050;
                                  const double sigma_signal_uncertainty_jet_res_ggH_800  = 0.050;
                                  const double sigma_signal_uncertainty_jet_res_ggH_900  = 0.050;
                                  const double sigma_signal_uncertainty_jet_res_ggH_1000 = 0.050;

                                  const double sigma_signal_uncertainty_jet_res_vbfH_600  = 0.145;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_700  = 0.101;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_800  = 0.096;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_900  = 0.099;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_1000 = 0.056;

                                  const double mean_mj_shift = 1.4 ;
                                  const double sigma_mj_smear = 1.16 ;
};       


struct SystematicUncertaintyHiggs_2jetBin{

                                  /////// effect on the mean 
                                  const double mean_signal_uncertainty_jet_scale_ggH_600  = 0.032;
                                  const double mean_signal_uncertainty_jet_scale_ggH_700  = 0.029;
                                  const double mean_signal_uncertainty_jet_scale_ggH_800  = 0.025;
                                  const double mean_signal_uncertainty_jet_scale_ggH_900  = 0.023;
                                  const double mean_signal_uncertainty_jet_scale_ggH_1000 = 0.020;

                                  const double mean_signal_uncertainty_jet_scale_vbfH_600  = 0.045;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_700  = 0.033;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_800  = 0.016;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_900  = 0.012;
                                  const double mean_signal_uncertainty_jet_scale_vbfH_1000 = 0.030;

                                  const double mean_signal_uncertainty_jet_res_ggH_600  = 0.006;
                                  const double mean_signal_uncertainty_jet_res_ggH_700  = 0.004;
                                  const double mean_signal_uncertainty_jet_res_ggH_800  = 0.002;
                                  const double mean_signal_uncertainty_jet_res_ggH_900  = 0.007;
                                  const double mean_signal_uncertainty_jet_res_ggH_1000 = 0.002;

                                  const double mean_signal_uncertainty_jet_res_vbfH_600  = 0.019;
                                  const double mean_signal_uncertainty_jet_res_vbfH_700  = 0.012;
                                  const double mean_signal_uncertainty_jet_res_vbfH_800  = 0.031;
                                  const double mean_signal_uncertainty_jet_res_vbfH_900  = 0.036;
                                  const double mean_signal_uncertainty_jet_res_vbfH_1000 = 0.046;


                                  /////// effect on the sigma 
                                  const double sigma_signal_uncertainty_jet_scale_ggH_600  = 0.080;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_700  = 0.040;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_800  = 0.021;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_900  = 0.023;
                                  const double sigma_signal_uncertainty_jet_scale_ggH_1000 = 0.017;

                                  const double sigma_signal_uncertainty_jet_scale_vbfH_600  = 0.036;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_700  = 0.049;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_800  = 0.085;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_900  = 0.123;
                                  const double sigma_signal_uncertainty_jet_scale_vbfH_1000 = 0.086;

                                  const double sigma_signal_uncertainty_jet_res_ggH_600  = 0.153;
                                  const double sigma_signal_uncertainty_jet_res_ggH_700  = 0.163;
                                  const double sigma_signal_uncertainty_jet_res_ggH_800  = 0.100;
                                  const double sigma_signal_uncertainty_jet_res_ggH_900  = 0.077;
                                  const double sigma_signal_uncertainty_jet_res_ggH_1000 = 0.057;

                                  const double sigma_signal_uncertainty_jet_res_vbfH_600  = 0.272;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_700  = 0.309;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_800  = 0.263;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_900  = 0.191;
                                  const double sigma_signal_uncertainty_jet_res_vbfH_1000 = 0.170;

                                  const double mean_mj_shift = 1.4 ;
                                  const double sigma_mj_smear = 1.16 ;
};       
