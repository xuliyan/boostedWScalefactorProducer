#include "MakePdf.h"

////////////////////////////////////////////
RooExtendPdf* MakeExtendedModel(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & channel, const std::string & wtagger_label, std::vector<std::string>* constraint, const int & ismc_wjet, const int & area_init_value){

  std::cout<<" "<<std::endl;
  std::cout<<"#########################################"<<std::endl;
  std::cout<<"## Make model : "<<label<<" "<<model<<" "<<area_init_value<<" ##"<<std::endl;
  std::cout<<"#########################################"<<std::endl;
  std::cout<<" "<<std::endl;

  RooRealVar* rrv_number = new  RooRealVar(("rrv_number"+label+"_"+channel+spectrum).c_str(),("rrv_number"+label+"_"+channel+spectrum).c_str(),500.,0.,1e5);
  // call the make RooAbsPdf method
  RooAbsPdf* model_pdf = NULL ;
  model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label,channel,constraint,ismc_wjet);
  std::cout<< "######## Model Pdf ########"<<std::endl;
  model_pdf->Print();
      
  // create the extended pdf
  RooExtendPdf* model_extended = new RooExtendPdf(("model"+label+"_"+channel+spectrum).c_str(),("model"+label+"_"+channel+spectrum).c_str(),*model_pdf, *rrv_number);
  std::cout<< "######## Model Extended Pdf ########"<<std::endl;
  model_extended->Print();
  // return the total extended pdf
  return model_extended;
}

/////////////////////////////


RooGaussian* addConstraint(RooRealVar* rrv_x, double x_mean, double x_sigma, std::vector<std::string>* ConstraintsList){

  //########### Gaussian contraint of a parameter of a pdf
  std::cout<<"########### Add to Contraint List some parameters  ############"<<std::endl;
  RooRealVar* rrv_x_mean = new RooRealVar((std::string(rrv_x->GetName())+"_mean").c_str(),(std::string(rrv_x->GetName())+"_mean").c_str(),x_mean);
  RooRealVar* rrv_x_sigma = new RooRealVar((std::string(rrv_x->GetName())+"_sigma").c_str(),(std::string(rrv_x->GetName())+"_sigma").c_str(),x_sigma);
  RooGaussian* constrainpdf_x = new RooGaussian(("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),*rrv_x,*rrv_x_mean, *rrv_x_sigma);
  //## import in the workspace and save the name of constriant pdf
  ConstraintsList->push_back(constrainpdf_x->GetName());
  return constrainpdf_x ;
}

//////////////////

// ---------------------------------------------                                                                                                                                    
RooAbsPdf* MakeModelTTbarControlSample(RooWorkspace* workspace ,const std::string & label, const std::string & modelName, const std::string & spectrum, const std::string & channel, const std::string & wtagger, const std::string & information, std::vector<std::string>* constraint){

  std::cout<<" "<<std::endl;
  std::cout<<"###################################################################################"<<std::endl;
  std::cout<<" ## make_Model_for_ttbar_controlsample : "<<label<<" "<<modelName<<" "<<information<<" "<<std::endl;
  std::cout<<"###################################################################################"<<std::endl;
  std::cout<<" "<<std::endl;

  RooRealVar* rrv_number_total = NULL , *eff_ttbar = NULL ; 
  RooFormulaVar *rrv_number = NULL ;

  if(TString(label).Contains("_ttbar_data") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    eff_ttbar        = new RooRealVar(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.9);
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "@0*@1", RooArgList(*rrv_number_total,*eff_ttbar));
  }

  else if(TString(label).Contains("_ttbar_data") and TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = workspace->var(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str());
    eff_ttbar        = workspace->var(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str());
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "(1-@0)*@1", RooArgList(*eff_ttbar,*rrv_number_total));
  }
  else if(TString(label).Contains("_ttbar_TotalMC") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    eff_ttbar        = new RooRealVar(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.9);
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "@0*@1", RooArgList(*eff_ttbar,*rrv_number_total));
  }
  else if(TString(label).Contains("_ttbar_TotalMC") and TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = workspace->var(("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str());
    eff_ttbar        = workspace->var(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str());
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "(1-@0)*@1", RooArgList(*eff_ttbar,*rrv_number_total)) ;
  }

  rrv_number_total->Print();
  eff_ttbar->Print();
  rrv_number->Print();

  RooAbsPdf* model_pdf = MakeGeneralPdf(workspace,label,modelName,spectrum,wtagger,channel,constraint,false);

  model_pdf->Print();
  RooExtendPdf* model = new RooExtendPdf(("model"+label+"_"+channel+spectrum).c_str(),("model"+label+"_"+channel+spectrum).c_str(),*model_pdf,*rrv_number);

  workspace->import(*model);
  workspace->pdf(("model"+label+"_"+channel+spectrum).c_str())->Print();
  return workspace->pdf(("model"+label+"_"+channel+spectrum).c_str());

}


//## change a dataset to a histpdf roofit object
void change_dataset_to_histpdf(RooWorkspace* workspace,RooRealVar* x,RooDataSet* dataset){
 
  std::cout<<"######## change the dataset into a histpdf  ########"<<std::endl;
  RooDataHist* datahist = dataset->binnedClone((std::string(dataset->GetName())+"_binnedClone").c_str(),(std::string(dataset->GetName())+"_binnedClone").c_str());
  RooHistPdf* histpdf = new RooHistPdf((std::string(dataset->GetName())+"_histpdf").c_str(),(std::string(dataset->GetName())+"_histpdf").c_str(),RooArgSet(*x),*datahist);
  dataset->Print();
  histpdf->Print();
  workspace->import(*histpdf);
}
  
// change from a dataset to a histogramm of Roofit
TH1F* change_dataset_to_histogram(RooRealVar* x, RooDataSet* dataset, const std::string & label, const int & BinWidth){

  std::cout<<"######## change the dataset into a histogramm for mj distribution ########"<<std::endl;
  RooDataHist* datahist = dataset->binnedClone((std::string(dataset->GetName())+"_binnedClone").c_str(),(std::string(dataset->GetName())+"_binnedClone").c_str());
  int nbin = int( (x->getMax()-x->getMin())/BinWidth);
  TString Name ;
  if(label ==""){
    Name.Form("histo_%s",dataset->GetName());
    return (TH1F*) datahist->createHistogram(Name.Data(),*x,RooFit::Binning(nbin,x->getMin(),x->getMax()));
  }
  else{
    Name.Form("histo_%s",label.c_str());
    return (TH1F*) datahist->createHistogram(Name.Data(),*x,RooFit::Binning(nbin,x->getMin(),x->getMax()));
  }

}

////////////////////////////////////

RooAbsPdf* get_mj_Model (RooWorkspace* workspace, const std::string & label, const std::string & channel){
  std::cout<<"model"+label+"_"+channel+"_mj"<<std::endl;
  return workspace->pdf(("model"+label+"_"+channel+"_mj").c_str());

}

RooAbsPdf* get_General_mj_Model(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing amd return a general mj model ############"<<std::endl;
  TString name ;

  if(TString(workspace->GetName()).Contains("4bias")) name.Form("rdataset4bias%s_%s_mj",label.c_str(),channel.c_str());
  else if (TString(workspace->GetName()).Contains("4fit")) name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());
  else name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());

  RooAbsData* rdataset_General_mj = workspace->data(name.Data());
  RooAbsPdf* model_General = get_mj_Model(workspace,label+model,channel);
  rdataset_General_mj->Print();
  model_General->Print();

  RooArgSet* parameters_General = model_General->getParameters(*rdataset_General_mj);
  TIter par = parameters_General->createIterator(); par.Reset();
  RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());

  if(fix == 1){
    while(param){
      param->setConstant(kTRUE);
      param->Print();
      param = dynamic_cast<RooRealVar*>(par.Next());
    }
  }
  else if(TString(label).Contains("WJets")){
    while(param){
      if(TString(param->GetName()).Contains("width")){
        param->setConstant(kTRUE);
        param->Print();
      }
      param = dynamic_cast<RooRealVar*>(par.Next());
    }
  }

  return model_General;
}

RooAbsPdf* get_TTbar_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing TTbar mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}

RooAbsPdf* get_STop_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing STop mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}

RooAbsPdf* get_VV_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing VV mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}
 
RooAbsPdf* get_WW_EWK_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing WW_EWK mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}

RooAbsPdf* get_ggH_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing ggH mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}

RooAbsPdf* get_vbfH_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing vbfH mj model ############"<<std::endl;
  return get_General_mj_Model(workspace,label,model,channel,fix);
}

RooAbsPdf* get_WJets_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix, const std::string & jetBin){

  std::cout<<"########### Fixing WJets mj model ############"<<std::endl;
  if(jetBin == "_2jet") return get_General_mj_Model(workspace,label,model,channel,fix);
  else{
    
    TString name ;
    if(TString(workspace->GetName()).Contains("4bias")) name.Form("rdataset4bias%s_%s_mj",label.c_str(),channel.c_str());
    else if (TString(workspace->GetName()).Contains("4fit")) name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());
    else name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());

    RooDataSet* rdataset_WJets_mj = (RooDataSet*) workspace->data(name.Data());
    rdataset_WJets_mj->Print();
    RooAbsPdf* model_WJets = get_mj_Model(workspace,label+model,channel);
    model_WJets->Print();
    RooArgSet* parameters_WJets = model_WJets->getParameters(*rdataset_WJets_mj);
    TIter par = parameters_WJets->createIterator(); par.Reset();
    RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
    while (param){
      TString paraName; paraName = Form("%s",param->GetName());
      if (paraName.Contains("rrv_width_ErfExp_WJets") or paraName.Contains("rrv_offset_ErfExp_WJets") or paraName.Contains("rrv_p1_User1_WJets")){
        param->setConstant(kTRUE);
      }
      param =  dynamic_cast<RooRealVar*>(par.Next());
    }
    return model_WJets ;
  }       
   
}

/////////////////////////////////////////

RooAbsPdf* get_mlvj_Model (RooWorkspace* workspace,const std::string & label,const std::string & region,const std::string & model, const std::string & channel){
   
  //#### get a generic mlvj model from the workspace
  std::cout<<"model"+label+region+model+"_"+channel+"_mlvj"<<std::endl;
  return workspace->pdf(("model"+label+region+model+"_"+channel+"_mlvj").c_str());
}

//#### get a general mlvj model and fiz the paramters --> for extended pdf
RooAbsPdf* get_General_mlvj_Model        (RooWorkspace* workspace,const std::string & label,const std::string & region,const std::string & model, const std::string & channel, const int & fix){
  std::cout<<"########### Fixing amd return a general mlvj model ############"<<std::endl;
  TString Name ; 
  if(TString(workspace->GetName()).Contains("4bias")) Name.Form("rdataset4bias%s%s_%s_mlvj",label.c_str(),region.c_str(),channel.c_str());
  else if(TString(workspace->GetName()).Contains("4fit")) Name.Form("rdataset%s%s_%s_mlvj",label.c_str(),region.c_str(),channel.c_str());

  RooAbsData* rdataset_General_mlvj = workspace->data(Name.Data());
  RooAbsPdf* model_General = get_mlvj_Model(workspace,label,region,model,channel);
  rdataset_General_mlvj->Print();
  model_General->Print();
  RooArgSet* parameters_General = model_General->getParameters(*rdataset_General_mlvj);
  TIter par = parameters_General->createIterator(); par.Reset();
  RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
  if(fix == 1){
    while(param){
      param->setConstant(kTRUE);
      param->Print();
      param = dynamic_cast<RooRealVar*>(par.Next());
    }
  }
  return model_General;
}

//###### get TTbar model mlvj in a region
RooAbsPdf* get_TTbar_mlvj_Model (RooWorkspace* workspace, const std::string & label, const std::string & region, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing TTbar mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}

//###### get STop model mlvj in a region
RooAbsPdf* get_STop_mlvj_Model (RooWorkspace* workspace, const std::string & label, const std::string & region, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing STop mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}

//###### get VV model mlvj in a region
RooAbsPdf* get_VV_mlvj_Model (RooWorkspace* workspace, const std::string & label,const std::string & region, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing VV mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}
 

//###### get WW_EWK model mlvj in a region
RooAbsPdf* get_WW_EWK_mlvj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & region,const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing WW_EWK mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}

//###### get ggH model mlvj in a region
RooAbsPdf* get_ggH_mlvj_Model  (RooWorkspace* workspace, const std::string & label,const std::string & region, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing ggH mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}


//###### get vbfH model mlvj in a region
RooAbsPdf* get_vbfH_mlvj_Model  (RooWorkspace* workspace, const std::string & label,const std::string & region, const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing vbfH mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}

//###### get Wjets
RooAbsPdf* get_WJets_mlvj_Model  (RooWorkspace* workspace, const std::string & label,const std::string & region,const std::string & model, const std::string & channel, const int & fix){

  std::cout<<"########### Fixing WJets mlvj model ############"<<std::endl;
  return get_General_mlvj_Model(workspace,label,region,model,channel,fix);
}

/////////////////////////////////////////////////
// fix a given model taking the label, and the region --> for extended pdf --> all the parameter of the pdf + normalization
void fix_Model(RooWorkspace* workspace, const std::string & label, const std::string & region, const std::string & spectrum, const std::string & model,const std::string & channel, const std::string & additional_info, const int & notExtended){
  
  RooAbsData* rdataset = NULL ;
  RooAbsPdf* model_pdf = NULL ;             
  std::cout<<"########### Fixing an Extended Pdf for mlvj ############"<<std::endl;
  if(spectrum == "_mj"){  
    TString Name ; Name.Form("rdataset%s%s_%s%s",label.c_str(),region.c_str(),channel.c_str(),spectrum.c_str()); 
    std::cout<<"DataSet Name "<<Name<<std::endl;
    rdataset = workspace->data(Name.Data());
    std::string label2 ; 
    if(notExtended == 1)
      label2 = std::string("_pdf")+label;
    else label2 = label;
    model_pdf = get_mj_Model(workspace,label2+model+additional_info,channel);        
  }
  else{  
    TString Name ; Name.Form("rdataset%s%s_%s%s",label.c_str(),region.c_str(),channel.c_str(),spectrum.c_str());
    std::cout<<"DataSet Name "<<Name<<std::endl;
    rdataset = workspace->data(Name.Data());
    std::string label2 ; 
    if(notExtended == 1)
      label2 = std::string("_pdf")+label;
    else label2 = label;
    model_pdf = get_mlvj_Model(workspace,label2.c_str(),region.c_str(),(model+additional_info).c_str(),channel.c_str());
  }

  rdataset->Print();
  model_pdf->Print();
  RooArgSet* parameters = model_pdf->getParameters(*rdataset);
  TIter par = parameters->createIterator(); 
  par.Reset();
  RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
  while (param){
    param->setConstant(kTRUE);
    param->Print();       
    param = dynamic_cast<RooRealVar*>(par.Next());
  }
  return ;
}

////////////////////////////////////////////////
void fix_Pdf(RooAbsPdf* model_pdf, RooArgSet* argset_notparameter){

  std::cout<<"########### Fixing a RooAbsPdf for mlvj or mj  ############"<<std::endl;     
  RooArgSet* parameters = model_pdf->getParameters(*argset_notparameter);
  TIter par = parameters->createIterator(); par.Reset();
  RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
  while (param){
    param->setConstant(kTRUE);
    param->Print();
    param = dynamic_cast<RooRealVar*>(par.Next());
  }
}

////// print the parameters of a given pdf --> only non constant ones
void ShowParam_Pdf(RooAbsPdf* model_pdf, RooArgSet* argset_notparameter){
 
  std::cout<<"########### Show Parameters of a input model  ############"<<std::endl;        
  model_pdf->Print();
  RooArgSet* parameters = model_pdf->getParameters(*argset_notparameter);
  TIter par = parameters->createIterator(); par.Reset();
  RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
  while (param){
    if(not param->isConstant()) param->Print();        
    param = dynamic_cast<RooRealVar*>(par.Next());
  }
}

////////////////////////////////////////////////
// clone Model Extend in a not extend skipping the normalization
void clone_Model(RooWorkspace* workspace, RooAbsPdf* inputPdf, const std::string & label, const std::string & region, const std::string & spectrum, const std::string & model, const std::string & channel, const int & fixParam){

  std::cout<<"########### Cloning an Extended Pdf for mlvj ############"<<std::endl;
  RooAbsData* rdataset = NULL ;
  RooAbsPdf*  model_pdf = NULL;
  if(spectrum == "_mj"){  
    TString Name ; 
    Name.Form("rdataset%s%s_%s%s",label.c_str(),region.c_str(),channel.c_str(),spectrum.c_str());
    std::cout<<"DataSet Name "<<Name<<std::endl;
    rdataset = workspace->data(Name.Data());
    model_pdf  = get_mj_Model(workspace,label+model,channel);
  }
  else{

    TString Name ; 
    Name.Form("rdataset%s%s_%s%s",label.c_str(),region.c_str(),channel.c_str(),spectrum.c_str());
    std::cout<<"DataSet Name "<<Name<<std::endl;
    rdataset  = workspace->data(Name.Data());
    model_pdf = get_mlvj_Model(workspace,label,region,model,channel);
  }

  rdataset->Print();
  inputPdf->Print();
  model_pdf->Print();
  RooArgSet* parameters  = inputPdf->getParameters(*rdataset);
  RooArgSet* parameters2 = model_pdf->getParameters(*rdataset);

  TIter par         = parameters->createIterator();
  TIter par2        = parameters2->createIterator();
  par.Reset();
  par2.Reset();

  RooRealVar* param  = dynamic_cast<RooRealVar*>(par.Next());
  RooRealVar* param2 = dynamic_cast<RooRealVar*>(par2.Next());
  if(parameters->getSize() < parameters2->getSize()){
    while(param){
      if(not TString(param2->GetName()).Contains("number")){
        param->setVal(param2->getVal());
        param->setError(param2->getError());
        if(fixParam) param->setConstant(kTRUE);
        if(TString(label.c_str()).Contains("WJets") and TString(model.c_str()).Contains("Erf") and spectrum == "_mj" and TString(param->GetName()).Contains("width"))
          param->setConstant(kTRUE);
        param->Print(); param2->Print();
        param  = dynamic_cast<RooRealVar*>(par.Next());
        param2 = dynamic_cast<RooRealVar*>(par2.Next());
      }
      else param2 = dynamic_cast<RooRealVar*>(par2.Next());
    }
  }
  else if(parameters->getSize() > parameters2->getSize()){
    while(param2){
      if(not TString(param->GetName()).Contains("number")){
        param->setVal(param2->getVal());
        param->setError(param2->getError());
        param->Print(); param2->Print();
        if(fixParam) param->setConstant(kTRUE);
        if(TString(label.c_str()).Contains("WJets") and TString(model.c_str()).Contains("Erf") and spectrum == "_mj" and TString(param->GetName()).Contains("width"))
          param->setConstant(kTRUE);
        param  = dynamic_cast<RooRealVar*>(par.Next());
        param2 = dynamic_cast<RooRealVar*>(par2.Next());
      }
      else param = dynamic_cast<RooRealVar*>(par.Next());
    }
  }
  else{
    while (param){
      param->setVal(param2->getVal());
      param->setError(param2->getError());
      param->Print(); param2->Print(); 
      if(fixParam) param->setConstant(kTRUE);
      if(TString(label.c_str()).Contains("WJets") and TString(model.c_str()).Contains("Erf") and spectrum == "_mj" and TString(param->GetName()).Contains("width"))
        param->setConstant(kTRUE);
      param  = dynamic_cast<RooRealVar*>(par.Next());
      param2 = dynamic_cast<RooRealVar*>(par2.Next());
    }

  }  
}
///////////////////////////////////////////////
RooAbsPdf* MakeGeneralPdf(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & wtagger_label, const std::string & channel, std::vector<std::string>* constraint, const int & ismc){
  
 
  RooRealVar* rrv_x = NULL ; 
  if(TString(spectrum).Contains("_mj"))   rrv_x = workspace->var("rrv_mass_j");
  if(TString(spectrum).Contains("_mlvj")) rrv_x = workspace->var("rrv_mass_lvj");
  if(TString(spectrum).Contains("_genHMass")) rrv_x = workspace->var("rrv_mass_gen_WW");

  // Gaus for the W peak
  if(model == "Gaus"){
    std::cout<< "########### Gaus for W peak  ############"<<std::endl;
    RooRealVar* rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
    RooRealVar* rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),7,1,15);
    RooGaussian* model_pdf = new  RooGaussian(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            
    return model_pdf ;
  }

  // sum of two exponential 
  if(model == "Exp" or model == "Exp_sr"){
    std::cout<< "######### Exp = levelled exp funtion for W+jets mlvj ############" <<std::endl;
    RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,0.);
    RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
            
    return model_pdf ;
  }

  if ( model == "ErfExp" ){
    std::cout<< "########### Erf*Exp for mj fit  ############"<<std::endl;
    RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.035,-0.05,-0.01);
    RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),115.,110.,120);
    RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),50.,40, 60.);
    RooErfExpPdf* model_pdf       = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
            
    return model_pdf ;
  }     

  // Exp+Gaus or mj spectrum
  if ( model == "ExpGaus"){
    std::cout<< "########### Exp + Gaus for mj  fit  ############"<<std::endl;
    RooRealVar* rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),0.05,-0.2,0.2);
    RooExponential* exp         = new RooExponential(("exp"+label+"_"+channel+spectrum).c_str(),("exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);

    RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
    RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);
    RooRealVar* rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
    RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

    RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp,*gaus),RooArgList(*rrv_high));
            
    return model_pdf ;
  }

  // Erf*Exp + Gaus for mj spectrum with offset == mean
  if(model == "ErfExpGaus_sp"){
    std::cout<< "########### Erf*Exp + Gaus for mj  fit  ############"<<std::endl;
    RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.2,0.);
    RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,200.);
           
    //            RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
    RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),85,78,95);
    RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);

    RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);

    RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
            
    RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
    RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
    return model_pdf ;
  }

  // Erf*Exp + 2Gaus for mj spectrum
  if(model == "2Gaus_ErfExp"){

    std::cout<< "########### 2Gaus + Erf*Exp for mj fit  ############"<<std::endl;
    //          double mean1_tmp      = 8.3141e+01; 
    //          double mean1_tmp      = 80; 
    double mean1_tmp      = 85; 
    double deltamean_tmp  = 6.9129e+00; 
    double sigma1_tmp     = 7.5145e+00; 
    double scalesigma_tmp = 3.6819e+00;           
    double frac_tmp       = 6.7125e-01; 

    RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
    RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
    RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

    RooRealVar* rrv_deltamean_gaus  = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp);//, deltamean_tmp, deltamean_tmp);
    RooFormulaVar* rrv_mean2_gaus   = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
    RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp);//, scalesigma_tmp, scalesigma_tmp);
    RooFormulaVar* rrv_sigma2_gaus  = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
    RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

    RooRealVar* rrv_frac_2gaus = new RooRealVar(("rrv_frac_2gaus"+label+"_"+channel+spectrum).c_str(),("rrv_frac_2gaus"+label+"_"+channel+spectrum).c_str(),frac_tmp);//, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

    double c0_tmp     = -2.9893e-02 ; 
    double offset_tmp = 7.9350e+01  ; 
    double width_tmp  = 3.3083e+01  ; 

    RooRealVar* rrv_c_ErfExp     = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
    RooRealVar* rrv_offset_ErfExp= new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(), offset_tmp);//, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
    RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp, width_tmp-10, width_tmp+10);
    RooErfExpPdf* erfexp = new RooErfExpPdf(("erfexp"+label+"_"+channel+spectrum).c_str(),("erfexp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

    RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(), 0.5,0.,1.);
    RooAddPdf*  model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfexp,*gaus1,*gaus2),RooArgList(*rrv_frac,*rrv_frac_2gaus),1);
          
    return model_pdf ;

  }

  // 2Gaus per ttbar fit
  if( model == "2Gaus_ttbar"){

    double mean1_tmp      = 8.3141e+01;     
    double deltamean_tmp  = 6.9129e+00;  
    double sigma1_tmp     = 7.5145e+00;     
    double scalesigma_tmp = 3.6819e+00; 
    double frac_tmp       = 6.7125e-01;       

    if( TString(label).Contains("herwig") and not TString(label).Contains("data")){

      mean1_tmp      = 8.3141e+01; 
      deltamean_tmp  = 6.9129e+00; 
      sigma1_tmp     = 7.5145e+00; 
      scalesigma_tmp = 3.6819e+00; 
      frac_tmp       = 6.7125e-01; 
    }
           
    RooRealVar* rrv_mean1_gaus  = NULL;
    RooRealVar* rrv_sigma1_gaus = NULL;

    // Constrain the peak and the width among electron and muon channel+spectrum to be the same in case of simultaneous fit
    if(channel == "el" and spectrum == "_mj"){                
      if(workspace->var(("rrv_mean1_gaus"+label+"_mu_mj").c_str()) and workspace->var(("rrv_sigma1_gaus+"+label+"_mu_mj").c_str())){
        rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus"+label+"_mu_mj").c_str());
        rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus"+label+"_mu_mj").c_str());
      }
      else{
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4);
      }
    }
    if((channel == "mu" or channel == "em") and spectrum == "_mj" ){
      if( workspace->var(("rrv_mean1_gaus"+label+"_el_mj").c_str()) and workspace->var(("rrv_sigma1_gaus+"+label+"_el_mj").c_str())){
        rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus"+label+"_el_mj").c_str());
        rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus"+label+"_el_mj").c_str());
      }
      else{
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
      }
    }
	    

    RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
    RooRealVar* rrv_deltamean_gaus = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp);//,deltamean_tmp-deltamean_tmp_err*5, deltamean_tmp+deltamean_tmp_err*5);
    RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp);//#,scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
    RooFormulaVar* rrv_mean2_gaus = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
    RooFormulaVar* rrv_sigma2_gaus = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
    RooGaussian*   gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

    RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);//,frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
    RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);
            
    return model_pdf ;
  }

  // Gaus + pol 2 background for resonant part in the fail sample
  if( model == "GausChebychev_ttbar_failtau2tau1cut"){

    double p0_tmp = -3.5459e-01; 
    double p1_tmp = -1.2790e-01; 
    double frac_tmp = 2.7324e-01; 

    if( TString(label).Contains("herwig") and not TString(label).Contains("data") and channel+spectrum == "mu"){

      p0_tmp = -6.6824e-01;  
      p1_tmp = -5.3365e-02;  
      frac_tmp = 2.0421e-01; 
    }
    // take the same gaussian used in the pass sample
    RooAbsPdf* gaus = NULL ;
    if( TString(label).Contains("data") and not TString(label).Contains("herwig"))
      gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
    if( TString(label).Contains("TotalMC") and not TString(label).Contains("herwig"))
      gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());

    if( TString(label).Contains("data") and TString(label).Contains("herwig"))
      gaus = workspace->pdf(("gaus1_ttbar_data_herwig_"+channel+spectrum).c_str());
    if( TString(label).Contains("TotalMC") and TString(label).Contains("herwig"))
      gaus = workspace->pdf(("gaus1_ttbar_TotalMC_herwig_"+channel+spectrum).c_str());
                
    RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);//,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
    RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);//,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
    RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );

    RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);//,frac_tmp-frac_tmp_err*4,frac_tmp+frac_tmp_err*4);
    RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);
          
    return model_pdf ;
  }

  // double gauss for fail sample
  if( model == "2Gaus_ttbar_failtau2tau1cut"){

    double mean1_tmp = 8.3209e+01;     
    double deltamean_tmp = 1.1427e+00; double deltamean_tmp_err = 1.03e+00;
    double sigma1_tmp = 7.4932e+00;    
    double scalesigma_tmp = 4.5922e+00; double scalesigma_tmp_err = 2.87e-01;
    double frac_tmp = 5.7910e-01;       double frac_tmp_err = 1.59e-02;

    RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
    RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
    RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

    RooRealVar* rrv_deltamean_gaus = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp, deltamean_tmp-deltamean_tmp_err*4, deltamean_tmp+deltamean_tmp_err*4);
    RooFormulaVar* rrv_mean2_gaus = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
    RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
    RooFormulaVar* rrv_sigma2_gaus = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
    RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

    RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
    RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);
            
    return model_pdf ;
  }

  // Erf*Exp in the ttbar sample
  if( model == "ErfExp_ttbar"){

    double c0_tmp         = -2.9893e-02 ; double c0_tmp_err     = 6.83e-03;
    double offset_tmp     = 7.9350e+01 ;  double offset_tmp_err = 9.35e+00;
    double width_tmp      = 3.3083e+01 ;  double width_tmp_err  = 2.97e+00;

    if( TString(label).Contains("herwig") and not TString(label).Contains("data")){
      c0_tmp     = -2.9357e-02 ; 
      offset_tmp = 7.9350e+01 ; 
      width_tmp  = 3.3216e+01 ; 
    }
                                                        
    RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2,c0_tmp+4e-2);
    RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
    RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,width_tmp-10, width_tmp+10);
                
    RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

    RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
    RooGaussian* gaus2 = addConstraint(rrv_offset_ErfExp,rrv_offset_ErfExp->getVal(),offset_tmp_err,constraint);
    RooGaussian* gaus3 = addConstraint(rrv_width_ErfExp,rrv_width_ErfExp->getVal(),width_tmp_err,constraint);

    workspace->import(*gaus1);
    workspace->import(*gaus2);
    workspace->import(*gaus3);


    return model_pdf ;
  }
	      
  // Erf*Exp for failing events
  if( model == "ErfExp_ttbar_failtau2tau1cut"){

    double c0_tmp     = -1.0143e-01;  double c0_tmp_err = 1.46e-02;
    double offset_tmp = 2.7718e+02 ; 
    double width_tmp  = 7.1891e+01 ;  

    if( TString(label).Contains("herwig") and not TString(label).Contains("data")){

      c0_tmp = -1.0141e-01 ; 
      offset_tmp = 2.6730e+02 ; 
      width_tmp = 7.1505e+01 ; 

    }
                  
    RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
    RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);//#offset_tmp-offset_tmp_err*5,offset_tmp+offset_tmp_err*5);
    RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp);//#width_tmp-width_tmp_err*5, width_tmp+width_tmp_err*5);
    RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
    // constrain just the slope of the exponential
    RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
    workspace->import(*gaus1);
    return model_pdf ;
  }
  
  // extreme fail cut
  if( model == "Exp_ttbar_extremefailtau2tau1cut"){

    double c0_tmp = -3.0278e-02 ; double c0_tmp_err = 5.16e-03;

    RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
    RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp);
    // constraint within one sigma
    RooGaussian* gaus = addConstraint(rrv_c_ErfExp,c0_tmp,c0_tmp_err,constraint);
    workspace->import(*gaus);
    return model_pdf ;
  }

  if( model == "Exp_bkg_extremefailtau2tau1cut"){
    double c0_tmp = -4.2105e-02 ; double c0_tmp_err = 2.61e-03;
    RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
    RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp);
    RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
    workspace->import(*gaus1);
            
    return model_pdf ;
  }

  // Erf*Exp + 2Gaus
  if( model == "ErfExp2Gaus_ttbar"){
    RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(), -9.72533e-02, -9.72533e-02-2*6.05691e-03, -9.72533e-02+2*6.05691e-03 );
    RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),2.21236e+02, 2.21236e+02-2*9.53939e+00, 2.21236e+02+2*9.53939e+00 );
    RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),5.79251e+01,5.79251e+01-2*1.22221e+00,5.79251e+01+2*1.22221e+00);
    RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

    RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),8.36227e+01,8.36227e+01-2*1.60148e-01,8.36227e+01+2*1.60148e-01);
    RooRealVar* rrv_mean2_gaus = new RooRealVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),8.45644e+01 ,8.45644e+01-2*9.21003e-01,8.45644e+01+2*9.21003e-01 );
    RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7.45239e+00,7.45239e+00-2*2.05052e-01,7.45239e+00+2*2.05052e-01);
    RooRealVar* rrv_sigma2_gaus = new RooRealVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),3.39773e+01,3.39773e+01-2*2.44179e+00,3.39773e+01+2*2.44179e+00);
    RooRealVar* rrv_high1 = new RooRealVar(("rrv_high1"+label+"_"+channel+spectrum).c_str(),("rrv_high1"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
    RooRealVar* rrv_high2 = new RooRealVar(("rrv_high2"+label+"_"+channel+spectrum).c_str(),("rrv_high2"+label+"_"+channel+spectrum).c_str(),0.5622,0.51,0.61);
    RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
    RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

    RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus1,*gaus2),RooArgList(*rrv_high1,*rrv_high2),1);
            
    return model_pdf ;

  }
  
  // Erf*Exp + 2Gaus fail sample
  if( model == "ErfExp2Gaus_ttbar_failtau2tau1cut"){
    RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-2.81053e-02 ,-2.81053e-02-2*5.82712e-03 ,-2.81053e-02+2*5.82712e-03 );
    RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),7.63158e+01,7.63158e+01-2*8.66322e+00,7.63158e+01+2*8.66322e+00);
    RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),3.38156e+01,3.38156e+01-2*3.02392e+00,3.38156e+01+2*3.02392e+00);
    RooErfExpPdf *erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

    RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),8.34758e+01,8.34758e+01-2*1.62468e-01,8.34758e+01+2*1.62468e-01);
    RooRealVar* rrv_mean2_gaus = new RooRealVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),8.99116e+01,8.99116e+01-2*7.47952e-01,8.99116e+01+2*7.47952e-01);
    RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7.45827e+00,7.45827e+00-2*1.95598e-01,7.45827e+00+2*1.95598e-01);
    RooRealVar* rrv_sigma2_gaus = new RooRealVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),2.90745e+01,2.90745e+01-2*1.77035e+00,2.90745e+01+2*1.77035e+00);
    RooRealVar* rrv_high1 = new RooRealVar(("rrv_high1"+label+"_"+channel+spectrum).c_str(),("rrv_high1"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
    RooRealVar* rrv_high2 = new RooRealVar(("rrv_high2"+label+"_"+channel+spectrum).c_str(),("rrv_high2"+label+"_"+channel+spectrum).c_str(),0.6323,0.57,0.7);
    RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
    RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

    RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus1,*gaus2),RooArgList(*rrv_high1,*rrv_high2),1);
            
    return model_pdf ;
  }

  return NULL ;

}
      
	

