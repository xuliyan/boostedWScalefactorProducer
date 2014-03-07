#include "BiasUtils.h"
#include "../PDFs/Util.cxx"

biasModelAnalysis::biasModelAnalysis( RooArgSet* observables, RooAbsPdf * generation_model, RooDataSet* generated_dataset, const int & nexp, const int & isMC){
  
  if(observables !=0 && observables!=NULL)  observables_ = observables ;
  else{ std::cout<<" null observable set --> terminate"<<std::endl; std::terminate(); }

  if(generation_model !=0 && generation_model!=NULL) model_generation_ = generation_model;
  else{ std::cout<<" null fitting model --> terminate"<<std::endl; std::terminate(); }

  param_generated_ = new RooArgList(*(model_generation_->getParameters(generated_dataset))); 
  generated_dataset_ = generated_dataset ;
 
  if(nexp >0) nexp_ = nexp ;
  else{ std::cout<<" null number of toys --> terminate"<<std::endl; std::terminate(); }

  isMC_        = isMC;  
}

biasModelAnalysis::biasModelAnalysis( const biasModelAnalysis & other){

  (*this).observables_ = other.observables_ ;
  (*this).model_fit_   = other.model_fit_;
  (*this).model_generation_   = other.model_generation_;
  (*this).nexp_        = other.nexp_ ;
  (*this).isMC_        = other.isMC_;
  (*this).tree_        = other.tree_;
} 

biasModelAnalysis::~biasModelAnalysis(){
  
  std::vector<const RooAbsData*>::const_iterator itData = generatedData_.begin();
  for( ; itData != generatedData_.end() ; ++itData)
   delete (*itData);
  
  std::vector<RooAbsPdf*>::const_iterator  itPdf = fittedPdf_.begin();
  for( ; itPdf != fittedPdf_.end() ; ++itPdf)
    delete (*itPdf);
 
  std::vector<const RooFitResult*>::const_iterator itRes = fitResults_.begin();
  for( ; itRes != fitResults_.end() ; ++itRes)
    delete (*itRes);
  
  //  if(mc_study_)         delete mc_study_ ;

  /*  
  if(parlist_)          delete parlist_ ;
  if(param_)            delete param_ ;
  if(param_generated_)  delete param_generated_ ;  
  if(generated_dataset_) delete generated_dataset_ ;*/
}


void biasModelAnalysis::setTree(TTree* tree){
 tree_ = tree ;
}


void biasModelAnalysis::setNToys(const int & nexp){
 nexp_ = nexp ;
}

void biasModelAnalysis::setIsMC(const int & isMC){
 isMC_ = isMC ;
}

void biasModelAnalysis::setFittingModel( RooAbsPdf* fitting_model){

  if(fitting_model !=0 && fitting_model!=NULL) model_fit_ = fitting_model;
  else{ std::cout<<" null fitting model --> terminate"<<std::endl; std::terminate(); }

}

void biasModelAnalysis::generateAndFitToys(int nevents, const std::string & fitRange){

  if(nevents >0) nevents_ = nevents;
  else{ std::cout<<" no event in the generation --> terminate"<<std::endl; std::terminate(); }
  
  fitRange_ = fitRange ;
  

  if((*this).isMC_){
   mc_study_ = new RooMCStudy(*((*this).model_generation_), 
                              *((*this).observables_),                                                                                                            
	 	              RooFit::FitModel(*((*this).model_fit_)),
                              RooFit::FitOptions(RooFit::Save(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Minimizer("Minuit2"),RooFit::Extended(kTRUE),RooFit::Range(fitRange.c_str())),
                              RooFit::Extended(kTRUE));
  }
  else{
    mc_study_ = new RooMCStudy(*((*this).model_generation_), 
                               *((*this).observables_),                                                                                                            
	   	               RooFit::FitModel(*((*this).model_fit_)),
                               RooFit::FitOptions(RooFit::Save(kTRUE),RooFit::SumW2Error(kFALSE),RooFit::Minimizer("Minuit2"),RooFit::Extended(kTRUE),RooFit::Range(fitRange.c_str())),
                               RooFit::Extended(kTRUE));
 
  }
  
  mc_study_->addModule(chi2_module_);
  mc_study_->generateAndFit((*this).nexp_,nevents_,1);                                                                                              
  
  TString name ;

  for( int iToy = 0 ; iToy < (*this).nexp_ ; iToy++){
    if ((*this).isMC_ == 1) name.Form("model_Total_mc_toy_%d",iToy);
    else name.Form("model_Total_data_toy_%d",iToy);
    generatedData_.push_back(mc_study_->genData(iToy));   
    fitResults_.push_back(mc_study_->fitResult(iToy));
    fittedPdf_.push_back(dynamic_cast<RooAbsPdf*>((*this).model_fit_->Clone(name.Data())));
  }
  
}


void biasModelAnalysis::createBranches(const std::string & fgen, const std::string & fres){

  fgen_ = fgen ;
  fres_ = fres ;

  std::string suffix = "";
  if( isMC_ == 0) suffix = "_wjet";  
  else suffix = "_data";

  parameter_ = NULL ;
  parameterResidual_ = NULL ;
  parameterError_  = NULL ;
  parameterPull_ = NULL ;

  TString branchName;
  bool branchCreated = false ;
  for (unsigned int iToy = 0; iToy < (*this).generatedData_.size() && branchCreated != true; iToy++){

    chi2_ = 0.;
    nLL_  = 0.;
    chi2_frame_ = 0.;

    if(!(*this).generatedData_.at(iToy) || !(*this).fitResults_.at(iToy)) continue ;             
    if((*this).fitResults_.at(iToy)->status()!= 0) continue;
    parlist_ = (RooArgList*)(mc_study_->fitParams(iToy));
    if(!parlist_) continue ;

    param_ = new RooArgList(*(fittedPdf_.at(iToy)->getParameters(generatedData_.at(iToy))));
    if(!param_) continue ;
 
    if(!param_generated_){                                                                                                                   
      std::cout<<" Not Find generation Model --> Exit from the code "<<std::endl;
      break ;
    }

    if(parameter_==NULL && parameterResidual_==NULL && parameterError_==NULL && parameterPull_==NULL){
     parameter_         = new float[int(parlist_->getSize())];
     parameterResidual_ = new float[int(parlist_->getSize())];
     parameterError_    = new float[int(parlist_->getSize())];
     parameterPull_     = new float[int(parlist_->getSize())];
    }

    int iPull = 0;    
    int iparNotConstant = 0;                                                                                                                                                              

    for( int ipar = 0; ipar < param_->getSize() ; ipar++){    

     if(param_->at(ipar)->isConstant()) continue;

     if (param_->at(ipar)->GetName() == parlist_->at(ipar)->GetName()){
         dynamic_cast<RooRealVar*>(param_->at(ipar))->setVal(dynamic_cast<RooRealVar*>(parlist_->at(ipar))->getVal());
         dynamic_cast<RooRealVar*>(param_->at(ipar))->setError(dynamic_cast<RooRealVar*>(parlist_->at(ipar))->getError());
     }

     if ((TString(param_->at(ipar)->GetName()).Contains("ggH") || TString(param_->at(ipar)->GetName()).Contains("vbfH")) && 
         !TString(param_->at(ipar)->GetName()).Contains("number")) continue;
     if ( TString(param_->at(ipar)->GetName()).Contains("rrv_fraction_ggH_vbf")) continue ;

     if( !TString(parlist_->at(ipar)->GetName()).Contains("ggH") || !TString(parlist_->at(ipar)->GetName()).Contains("vbfH")){
      if (!TString(parlist_->at(ipar)->GetName()).Contains("number")){

       branchName.Form("%s%s",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameter_[iparNotConstant],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_error",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterError_[iparNotConstant],std::string(branchName+"/F").c_str());
       std::cout<<" branchName "<<branchName<<" position "<<iparNotConstant<<" iPull "<<iPull<<std::endl;      
       if ( fgen_ == fres_){
         branchName.Form("%s%s_residual",parlist_->at(ipar)->GetName(),suffix.c_str());
         tree_->Branch(branchName.Data(),&parameterResidual_[iPull],std::string(branchName+"/F").c_str());
         branchName.Form("%s%s_pull",parlist_->at(ipar)->GetName(),suffix.c_str());
         tree_->Branch(branchName.Data(),&parameterPull_[iPull],std::string(branchName+"/F").c_str());
         iPull = iPull +1 ;
       }
      }
       
      else{

       branchName.Form("%s%s",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameter_[iparNotConstant],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_error",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterError_[iparNotConstant],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_residual",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterResidual_[iPull],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_pull",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterPull_[iPull],std::string(branchName+"/F").c_str());
       std::cout<<" branchName "<<branchName<<" position "<<iparNotConstant<<" iPull "<<iPull<<std::endl;      
       iPull = iPull +1 ;
      }

    }
    else{

       branchName.Form("%s%s",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameter_[iparNotConstant],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_error",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterError_[iparNotConstant],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_residual",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterResidual_[iPull],std::string(branchName+"/F").c_str());
       branchName.Form("%s%s_pull",parlist_->at(ipar)->GetName(),suffix.c_str());
       tree_->Branch(branchName.Data(),&parameterPull_[iPull],std::string(branchName+"/F").c_str());
       std::cout<<" branchName "<<branchName<<" position "<<iparNotConstant<<" iPull "<<iPull<<std::endl;      
       iPull = iPull +1 ;
    }
     iparNotConstant = iparNotConstant +1 ;
   }

   if(parlist_->find("NLL")){
       branchName.Form("nLL");
       tree_->Branch(branchName.Data(),&nLL_,std::string(branchName+"/F").c_str());
   }
   if(parlist_->find("chi2red")){
       branchName.Form("chi2red");
       tree_->Branch(branchName.Data(),&chi2_,std::string(branchName+"/F").c_str());
   }
  
  branchName.Form("chi2red_frame");
  tree_->Branch(branchName.Data(),&chi2_frame_,std::string(branchName+"/F").c_str());
       
  branchCreated = true ;     
 
  }

  return;
}


void biasModelAnalysis::fillBranches(const int & ttbarcontrolregion, const int & fitjetmass){
                                                                                                                                                   
  RooArgList*  notconstantparameters = NULL;   
  RooDataHist* data_binned = NULL; 
  RooAbsReal*  chi2Var = NULL ; 
   
  for (unsigned int iToy = 0; iToy < generatedData_.size() ; iToy++){

    chi2_ = 0.;
    nLL_  = 0.;
    chi2_frame_ = 0.;

    int iparNotConstant = 0;                                                                                                                                                              
    int iGenerated = 0 ;                                                                                                                                                                  
    int iPull = 0;    

    if(!(*this).generatedData_.at(iToy) or !(*this).fitResults_.at(iToy)) continue ;                                                                                                      
    if((*this).fitResults_.at(iToy)->status()!= 0) continue;

    parlist_ = (RooArgList*)(mc_study_->fitParams(iToy));
    if(!parlist_) continue ;

    for( int iparameter = 0 ; iparameter < parlist_->getSize() ; iparameter++){
     parameter_[iparameter] = 0 ;
     parameterResidual_[iparameter] = 0 ;
     parameterError_[iparameter] = 0 ;
     parameterPull_[iparameter] = 0 ;
    }

    param_ = new RooArgList(*(fittedPdf_.at(iToy)->getParameters(generatedData_.at(iToy))));
    if(!param_) continue ;
 
    if(!param_generated_){                                                                                                                   
      std::cout<<" Not Find generation Model --> Exit from the code "<<std::endl;
      break ;
    }

    for( int ipar = 0; ipar < param_->getSize() ; ipar++){    

     if(param_->at(ipar)->isConstant()) continue;

     RooRealVar* rrv_param   = dynamic_cast<RooRealVar*>(param_->at(ipar));
     RooRealVar* rrv_parlist = dynamic_cast<RooRealVar*>(parlist_->at(ipar));

     if (param_->at(ipar)->GetName() == parlist_->at(ipar)->GetName()){
         rrv_param->setVal(rrv_parlist->getVal());
         rrv_param->setError(rrv_parlist->getError());
     }

     if ((TString(param_->at(ipar)->GetName()).Contains("ggH") || TString(param_->at(ipar)->GetName()).Contains("vbfH")) && 
         !TString(param_->at(ipar)->GetName()).Contains("number")) continue;
     if ( TString(param_->at(ipar)->GetName()).Contains("rrv_fraction_ggH_vbf")) continue ;

     if( !TString(parlist_->at(ipar)->GetName()).Contains("ggH") || !TString(parlist_->at(ipar)->GetName()).Contains("vbfH")){
       parameterError_[iparNotConstant] = rrv_parlist->getError(); 
 
       if(!TString(parlist_->at(ipar)->GetName()).Contains("number")){ 
     
	if (ttbarcontrolregion == 0){
         while ( ( TString(param_generated_->at(iGenerated)->GetName()).Contains("number") || TString(param_generated_->at(iGenerated)->GetName()).Contains("rrv_mass_j") || 
		   TString(param_generated_->at(iGenerated)->GetName()).Contains("rrv_mass_lvj")) && iGenerated <= param_generated_->getSize())                              
	     iGenerated = iGenerated +1 ;                                                                                                                                           
	}
        else{
         while ( ( TString(param_generated_->at(iGenerated)->GetName()).Contains("number") || TString(param_generated_->at(iGenerated)->GetName()).Contains("rrv_mass_j") || 
		   TString(param_generated_->at(iGenerated)->GetName()).Contains("rrv_mass_lvj")) && iGenerated <= param_generated_->getSize())                              
             iGenerated = iGenerated +1 ;                                                                                                                                            
       }	
          
	parameter_[iparNotConstant] = rrv_parlist->getVal();                                                                 
        if( fgen_ == fres_) {
	  parameterResidual_[iPull] = rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(param_generated_->at(iGenerated))->getVal();    
          parameterPull_[iPull] = (rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(param_generated_->at(iGenerated))->getVal())/rrv_parlist->getError();
	  iPull ++ ;
        }
        
       }
       else if(TString(parlist_->at(ipar)->GetName()).Contains("number")){
	 if(ttbarcontrolregion == 0){
          parameter_[iparNotConstant] = rrv_parlist->getVal();                                                                                                          
	  parameterResidual_[iPull] = rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal();                   
          parameterPull_[iPull] = (rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal())/rrv_parlist->getError();
         }      
	 else{
 	  parameter_[iparNotConstant] = rrv_parlist->getVal();                                                                                                          
	  parameterResidual_[iPull] = rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal();                   
          parameterPull_[iPull] = (rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal())/rrv_parlist->getError();
         }
	 iPull = iPull +1;                                                                                                                                                       
       }
       iGenerated = iGenerated +1 ;                                                                                                                                                
     }     
     else{

       parameterError_[iparNotConstant] = rrv_parlist->getError(); 
       parameter_[iparNotConstant] = rrv_parlist->getVal();                                                                                                          
       parameterResidual_[iPull] = rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal();                   
       parameterPull_[iPull] = (rrv_parlist->getVal()-dynamic_cast<RooRealVar*>(parlist_->find("ngen"))->getVal())/rrv_parlist->getError();
       iPull = iPull +1;                                                                                                                                                       
     }

     iparNotConstant = iparNotConstant+1;

    }

    
    if(parlist_->find("NLL"))
      nLL_ =  dynamic_cast<RooRealVar*>(parlist_->find("NLL"))->getVal();                                                                                                    

    if(parlist_->find("chi2red"))                                                                                                                                 
      chi2_ = dynamic_cast<RooRealVar*>(parlist_->find("chi2red"))->getVal();
                                       
    notconstantparameters = dynamic_cast<RooArgList*>(param_->selectByAttrib("Constant",kFALSE));                           
    
    if(notconstantparameters){
      TString name ; name.Form("data_%d",iToy) ;
      if(data_binned!=NULL) delete data_binned; 
      data_binned    = new RooDataHist(name.Data(),name.Data(),(*(*this).observables_),*generatedData_.at(iToy));        
      if(isMC_ && ! fitjetmass){         
       chi2Var = fittedPdf_.at(iToy)->createChi2(*data_binned,RooFit::Extended(kTRUE),RooFit::SumW2Error(kTRUE));   
       chi2_frame_    = chi2Var->getVal()/(dynamic_cast<RooRealVar*>(observables_->find("rrv_mass_lvj"))->getBins()-notconstantparameters->getSize());
      }
      else if(!isMC_ && !fitjetmass){
       chi2Var = fittedPdf_.at(iToy)->createChi2(*data_binned,RooFit::Extended(kTRUE),RooFit::SumW2Error(kFALSE));   
       chi2_frame_    = chi2Var->getVal()/(dynamic_cast<RooRealVar*>(observables_->find("rrv_mass_lvj"))->getBins()-notconstantparameters->getSize());
      }
      else if(!isMC_ && fitjetmass){
       chi2Var = fittedPdf_.at(iToy)->createChi2(*data_binned,RooFit::Extended(kTRUE),RooFit::SumW2Error(kFALSE));   
       chi2_frame_    = chi2Var->getVal()/(dynamic_cast<RooRealVar*>(observables_->find("rrv_mass_j"))->getBins()-notconstantparameters->getSize());
      }
    }                                                                                                                 
    
    tree_->Fill();                                                                                                                                                           
         
  }

  return; 
  
}


void biasModelAnalysis::saveToysPlots(const int & nPlots, const int & fitjetmass){

  RooRealVar*   rrv_x = NULL ;
  RooFitResult* fres = NULL ;

  TLatex*  banner = new TLatex(0.3,0.96,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow l #nu "));
  banner->SetNDC(); 
  banner->SetTextSize(0.04);

  TString Title ;

  for(unsigned int iToy = 0 ; iToy < (*this).generatedData_.size() ; iToy++){
    if(iToy%nPlots != 0 ) continue ; 

    if(!(*this).generatedData_.at(iToy) or !(*this).fitResults_.at(iToy)) continue ;
    if((*this).fitResults_.at(iToy)->status()!= 0) continue;

    Title.Form("frame_generatedToys_wjet_%d",iToy);    
    
    if(fitjetmass) rrv_x = dynamic_cast<RooRealVar*>((*this).observables_->find("rrv_mass_j"));
    else     if(fitjetmass) rrv_x = dynamic_cast<RooRealVar*>((*this).observables_->find("rrv_mass_mlvj"));

    RooPlot* mplot = rrv_x->frame(RooFit::Title(Title.Data()), RooFit::Bins(rrv_x->getBins()));
    (*this).generatedData_[iToy]->plotOn(mplot,RooFit::MarkerSize(1.5), RooFit::Invisible(), RooFit::XErrorSize(0));
    fres = dynamic_cast<RooFitResult*>((*this).fitResults_[iToy]->Clone("fres"));
    draw_error_band_extendPdf(*((RooAbsData*)(*this).generatedData_[iToy]), *((*this).fittedPdf_[iToy]),fres,mplot,2,"L");
    Title.Form("%s_curve",(*this).fittedPdf_[iToy]->GetName());
    (*this).fittedPdf_[iToy]->plotOn(mplot,RooFit::Name(Title.Data()));
    Title.Form("%s_curve",generatedData_[iToy]->GetName());   
    (*this).generatedData_[iToy]->plotOn(mplot,RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),RooFit::Name(Title.Data()));
  
    Title.Form("canvas_generatedToys_wjet_%d",iToy);
    canvasVector_.push_back(new TCanvas(Title.Data(),""));
    canvasVector_.back()->cd();
    Title.Form("pad1_%d",iToy);
    TPad* pad1 = new TPad(Title.Data(),Title.Data(),0.,0.24,0.99,1.);
    pad1->Draw();
    Title.Form("pad2_%d",iToy);
    TPad* pad2 = new TPad(Title.Data(),Title.Data(),0.,0.,0.99,0.24);
    pad2->Draw();
    pad1->cd();

    mplot->GetXaxis()->SetTitleOffset(1.1);
    mplot->GetYaxis()->SetTitleOffset(1.3);
    mplot->GetXaxis()->SetTitleSize(0.05);
    mplot->GetYaxis()->SetTitleSize(0.05);                                                                                                                                               
    mplot->GetXaxis()->SetLabelSize(0.045);                                                                                                                                              
    mplot->GetYaxis()->SetLabelSize(0.045);                                                                                                                                              
    mplot->Draw();                                                                                                                                                                      
    banner->Draw();

    pad2->cd();                                                                                                                                                                         
    RooHist* hpull = mplot->pullHist();
    double x = 0. ;
    double y = 0. ;

    for( int ipoint = 0 ; ipoint < hpull->GetN() ; ipoint++){
     hpull->GetPoint(ipoint,x,y);
     if(y == 0) hpull->SetPoint(ipoint,x,10);
    }
  
    RooPlot* mplot_pull = rrv_x->frame(RooFit::Title("Pull Distribution"), RooFit::Bins(int(rrv_x->getBins())));
    TLine* medianLine = new TLine(rrv_x->getMin(),0.,rrv_x->getMax(),0); 
    medianLine->SetLineWidth(2); 
    medianLine->SetLineColor(kRed);
    mplot_pull->addObject(medianLine);
    mplot_pull->addPlotable(hpull,"P");
    mplot_pull->SetTitle("");
    mplot_pull->GetXaxis()->SetTitle("");
    mplot_pull->GetYaxis()->SetRangeUser(-5,5);
    mplot_pull->GetYaxis()->SetTitleSize(0.10);
    mplot_pull->GetYaxis()->SetLabelSize(0.10);
    mplot_pull->GetXaxis()->SetTitleSize(0.10);
    mplot_pull->GetXaxis()->SetLabelSize(0.10);
    mplot_pull->GetYaxis()->SetTitleOffset(0.40);
    mplot_pull->GetYaxis()->SetTitle("#frac{data-fit}{#sigma_{data}}");
    mplot_pull->GetYaxis()->CenterTitle();

    mplot_pull->Draw();                                                                                                                                                                 
    mplot_pull->GetXaxis()->SetLabelSize(0.15);                                                                                                                                          
    mplot_pull->GetYaxis()->SetLabelSize(0.15);                                                                                                                                          
    mplot_pull->GetYaxis()->SetTitleSize(0.15);                                                                                                                                          
    mplot_pull->GetYaxis()->SetNdivisions(205);                                                                                                                           
    canvasVector_.back()->Write();
  }
  
  return ;

}
