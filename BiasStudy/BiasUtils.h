#include <algorithm>
#include <vector>
#include <string>
#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooFitResult.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooAddPdf.h"
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooMCStudy.h"
#include "RooChi2Var.h"
#include "RooChi2MCSModule.h"
#include "TTree.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TLine.h"

#include "../PlotStyle/PlotUtils.h"

class biasModelAnalysis{

 public:
       
  biasModelAnalysis(){};
  ~biasModelAnalysis();
  biasModelAnalysis( RooArgSet *, RooAbsPdf* , RooDataSet*, const int &, const int &);
  biasModelAnalysis( const biasModelAnalysis &);
  biasModelAnalysis* clone() {return new biasModelAnalysis(*this); };

  void generateAndFitToys(int nevents, const std::string & = "" );
  void createBranches(const std::string &, const std::string &, const int &);
  void fillBranches(const int &, const int &, RooWorkspace&, std::map<std::string,std::string> &,const std::string & ="");

  void setFittingModel(RooAbsPdf*);
  void setTree(TTree*);
  void setNToys(const int &);
  void setIsMC(const int &);
  void setPdfInformation(const std::string &, const std::string &, const std::string &, const std::string &);
  void setBackgroundPdfCore(RooAbsPdf*);
  void setSignalInjection(RooAbsPdf*, const float & = 0, const float & = 1);
  void saveToysPlots(const int &, const int &, const int & = 0 );
 
  private: 

   TTree* tree_ ;

   RooMCStudy*   mc_study_ ;
   RooAbsPdf*    model_generation_ ;
   RooAbsPdf*    model_fit_;
   RooAbsPdf*    model_bkg_data_;
   RooArgSet*    observables_;
   RooDataSet*   generated_dataset_;
   RooChi2MCSModule chi2_module_ ;

   RooArgList* parlist_;
   RooArgList* param_;
   RooArgList* param_generated_;

   std::string fres_ ;
   std::string fgen_ ;
   std::string fitRange_ ;
 
   std::string mlvjregion_ ;
   std::string channel_ ;
   std::string spectrum_ ;
   std::string label_ ;
  
   std::vector<const RooAbsData*>   generatedData_; 
   std::vector<RooAbsPdf*>          fittedPdf_;
   std::vector<const RooFitResult*> fitResults_;

   float* parameter_ ;
   float* parameterResidual_ ;                                                                      
   float* parameterError_;
   float* parameterPull_ ;

   std::vector<TCanvas*> canvasVector_ ;
  
   float chi2_;
   float nLL_;
   float chi2_frame_;
   float numberSignalEvents_ ;
   float scalesignalwidth_;
   int   nexp_ ;
   bool  isMC_ ;
   int   nevents_ ;
};

