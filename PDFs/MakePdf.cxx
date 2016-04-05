#include "MakePdf.h"

////////////////////////////////////////////
RooExtendPdf* MakeExtendedModel(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & channel, const std::string & wtagger_label, std::vector<std::string>* constraint, const int & ismc_wjet, const int & area_init_value){

  std::cout<<" "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<"## Make model : "<<label<<" "<<model<<" "<<area_init_value<<" ##"<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" "<<std::endl;

  RooRealVar* rrv_number = new  RooRealVar(("rrv_number"+label+"_"+channel+spectrum).c_str(),("rrv_number"+label+"_"+channel+spectrum).c_str(),500.,0.,1e5);
  // call the make RooAbsPdf method
  RooAbsPdf* model_pdf = NULL ;
  model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label,channel,constraint,ismc_wjet);
  // std::cout<< "######## Model Pdf ########"<<std::endl;
  // model_pdf->Print();
      
  // create the extended pdf
  RooExtendPdf* model_extended = new RooExtendPdf(("model"+label+"_"+channel+spectrum).c_str(),("model"+label+"_"+channel+spectrum).c_str(),*model_pdf, *rrv_number);
  std::cout<< "######## Model Extended Pdf ########"<<std::endl;
  // model_extended->Print();
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
  std::cout<<""<<std::endl;
  std::cout<<" make_Model_for_ttbar_controlsample (label, modelName, information) : "<<label<<" "<<modelName<<" "<<information<<" "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" "<<std::endl;

  RooRealVar* rrv_number_total = NULL , *eff_ttbar = NULL ; 
  RooFormulaVar *rrv_number = NULL ;

  if(TString(label).Contains("_ttbar_data") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    eff_ttbar        = new RooRealVar(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.99);
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "@0*@1", RooArgList(*rrv_number_total,*eff_ttbar));
  }

  else if(TString(label).Contains("_ttbar_data") and TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = workspace->var(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str());
    eff_ttbar        = workspace->var(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str());
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "(1-@0)*@1", RooArgList(*eff_ttbar,*rrv_number_total));
  }
  else if(TString(label).Contains("_ttbar_TotalMC") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    eff_ttbar        = new RooRealVar(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.99);
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "@0*@1", RooArgList(*eff_ttbar,*rrv_number_total));
  }
  else if(TString(label).Contains("_ttbar_TotalMC") and TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = workspace->var(("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str());
    eff_ttbar        = workspace->var(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str());
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "(1-@0)*@1", RooArgList(*eff_ttbar,*rrv_number_total)) ;
  }


  RooAbsPdf* model_pdf = MakeGeneralPdf(workspace,label,modelName,spectrum,wtagger,channel,constraint,false);

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
  // dataset->Print();
  // histpdf->Print();
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

  std::cout<<" Fixing and return a general mj model "<<std::endl;
  TString name ;

  if(TString(workspace->GetName()).Contains("4bias")) name.Form("rdataset4bias%s_%s_mj",label.c_str(),channel.c_str());
  else if (TString(workspace->GetName()).Contains("4fit")) name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());
  else name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());

  RooAbsData* rdataset_General_mj = workspace->data(name.Data());
  RooAbsPdf* model_General = get_mj_Model(workspace,label,channel);
  // rdataset_General_mj->Print();
  // model_General->Print();

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
 

RooAbsPdf* get_WJets_mj_Model  (RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & channel, const int & fix, const std::string & jetBin){

  std::cout<<"########### Fixing WJets mj model ############"<<std::endl;
  if(jetBin == "_2jet") return get_General_mj_Model(workspace,label,model,channel,fix);
  else{
    
    TString name ;
    if(TString(workspace->GetName()).Contains("4bias")) name.Form("rdataset4bias%s_%s_mj",label.c_str(),channel.c_str());
    else if (TString(workspace->GetName()).Contains("4fit")) name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());
    else name.Form("rdataset%s_%s_mj",label.c_str(),channel.c_str());

    RooDataSet* rdataset_WJets_mj = (RooDataSet*) workspace->data(name.Data());
    // rdataset_WJets_mj->Print();
    RooAbsPdf* model_WJets = get_mj_Model(workspace,label+model,channel);
    // model_WJets->Print();
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
  return workspace->pdf(("model"+label+region+model+"_"+channel+"_mlvj").c_str());
}

//#### get a general mlvj model and fiz the paramters --> for extended pdf
RooAbsPdf* get_General_mlvj_Model        (RooWorkspace* workspace,const std::string & label,const std::string & region,const std::string & model, const std::string & channel, const int & fix){
  std::cout<<"Fixing and return a general mlvj model "<<std::endl;
  TString Name ; 
  if(TString(workspace->GetName()).Contains("4bias")) Name.Form("rdataset4bias%s%s_%s_mlvj",label.c_str(),region.c_str(),channel.c_str());
  else if(TString(workspace->GetName()).Contains("4fit")) Name.Form("rdataset%s%s_%s_mlvj",label.c_str(),region.c_str(),channel.c_str());

  RooAbsData* rdataset_General_mlvj = workspace->data(Name.Data());
  RooAbsPdf* model_General = get_mlvj_Model(workspace,label,region,model,channel);
  // rdataset_General_mlvj->Print();
  //   model_General->Print();
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

  // rdataset->Print();
  // model_pdf->Print();
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

///////////////////////////////////////////////
RooAbsPdf* MakeGeneralPdf(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & wtagger_label, const std::string & channel, std::vector<std::string>* constraint, const int & ismc){
  
  std::cout<< "Making general PDF for wtagger_label = "<<wtagger_label.c_str() << " model = "<< model.c_str() << " label = "<< label.c_str()<< " ismc = "<< ismc << " channel = "<< channel.c_str()<< std::endl;

  RooRealVar* rrv_x = NULL ; 
  if(TString(spectrum).Contains("_mj"))   rrv_x = workspace->var("rrv_mass_j");
  if(TString(spectrum).Contains("_mlvj")) rrv_x = workspace->var("rrv_mass_lvj");
  
  // FOR 0v45 ONLY!!!!!
  if(TString(wtagger_label.c_str()).Contains("0v45") or TString(wtagger_label.c_str()).Contains("0v50") or TString(wtagger_label.c_str()).Contains("0v55")){
    
     // sum of two exponential 
    if(model == "Exp"){
      std::cout<< "######### Exp = levelled exp funtion for W+jets mlvj ############" <<std::endl;
      RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.030, -0.2, 0.05);
      RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);       
      return model_pdf ;
    }

    if ( model == "ErfExp" ){
      std::cout<< "########### Erf*Exp for mj fit  ############"<<std::endl;
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.026);//,-0.05, 0.05);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),41.,0.,100);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,1.,100.);
      
      if(TString(label).Contains("_WJets0") ) {
        rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0279,-0.5,0.);
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,60.,75.);
        rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),23.8,20.,30.);
        
        if(TString(label).Contains("fail") ) {
          rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0294,-0.05,0.05);
          rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),39.,30.,50.);
          rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),22.5,10.,30.);
          
        }
      }
      
      RooErfExpPdf* model_pdf       = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);          
      return model_pdf ;
    }     

    // Exp+Gaus or mj spectrum
    if ( model == "ExpGaus"){
      
      RooRealVar* rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,70,90);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,40);
      RooRealVar* rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.,0.,1.);
      
      if( TString(label).Contains("_STop_failtau2tau1cut" ) ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.03,-0.5,0.5);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,60,120); //Too narrow limits here often lead to error!! eg max 80
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,60);
      }
      if( TString(label).Contains("_VV_") ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.2,0.2);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),90,78.,95.);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4.,50.);
      }if( TString(label).Contains("_VV_fail") ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.02,-0.2,0.2);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),82.,78.,95.);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),10,4.,20.);
      }
      
      
      RooExponential* exp         = new RooExponential(("exp"+label+"_"+channel+spectrum).c_str(),("exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp,*gaus),RooArgList(*rrv_high));
            
      return model_pdf ;
    }

    // Erf*Exp + Gaus for mj spectrum with offset == mean
    if(model == "ErfExpGaus_sp"){
      std::cout<< "########### Erf*Exp + Gaus for mj  fit  ############"<<std::endl;
      RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.04,-0.2,0.);
      RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,300.);
           
      //            RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
      RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,40,100);
      RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,0.,40);

      RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);

      RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
            
      RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
      return model_pdf ;
    }


    // 2Gaus per ttbar fit
    if( model == "2Gaus_ttbar"){

      double mean1_tmp = 8.406717e+01;   
      double sigma1_tmp = 7.7633e+00;      
      double mean2_tmp = 90.;   
      double sigma2_tmp = 30.;        
      float rangeMean = 10. ;
      float rangeWidth = 4. ;
      double frac_tmp = 6.24209e-01; double frac_tmp_err = 0.1;

      RooRealVar* rrv_mean1_gaus  = NULL;
      RooRealVar* rrv_sigma1_gaus = NULL;
      
      rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
      rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );    
      
      if( TString(label.c_str()).Contains("realW_fail")){       
        mean2_tmp = 70.;   
        sigma2_tmp = 40.;     
        frac_tmp = 0.29; frac_tmp_err = 0.1;
        sigma1_tmp = 5.3633e+00;   
        rangeMean = 10;
        rangeWidth = 20. ;
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );
       }
       
       if( TString(wtagger_label.c_str()).Contains("DDT")){
         rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp);//, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
         rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp);//, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );
         
       }
         
       if( TString(label.c_str()).Contains("fail") ){
         mean2_tmp = 70.;   
         sigma2_tmp = 40.;     
         // frac_tmp = 0.29; frac_tmp_err = 0.29;
         sigma1_tmp = 5.3633e+00;   
         rangeMean = 10;
         rangeWidth = 20. ;
         if( TString(label).Contains("data") ) {
           rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus_ttbar_data_"+channel+spectrum).c_str());
           rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus_ttbar_data_"+channel+spectrum).c_str());
         }
         else if( TString(label).Contains("TotalMC") ) {
           rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus_ttbar_TotalMC_"+channel+spectrum).c_str());
           rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus_ttbar_TotalMC_"+channel+spectrum).c_str());
         }
       }
      
      RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
      
      RooRealVar* rrv_mean2_gaus  = NULL;
      RooRealVar* rrv_sigma2_gaus = NULL;
      rrv_mean2_gaus  = new RooRealVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),96.1,20,130.);
      rrv_sigma2_gaus = new RooRealVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),26.5,0,700);
      
      RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);//,frac_tmp-frac_tmp_err, frac_tmp+frac_tmp_err);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);
            
      return model_pdf ;
    }
    
// //Double Crystal ball shape 
//      if(model == "Double_CB"){
//        std::cout<< "########### Double Cystal Ball ############"<<std::endl;
//
//       RooRealVar* rrv_mean_CB = NULL ;
//       RooRealVar* rrv_sigma_CB = NULL ;
//       RooRealVar* rrv_alpha_CB = NULL ;
//       RooRealVar* rrv_n_CB = NULL ;
//       RooRealVar* rrv_n2_CB = NULL ;
//       RooRealVar* rrv_alpha2_CB = NULL ;
//
//
//       if(TString(label).Contains("realW")){
//          rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),80.,78.,88.);
//          rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),8.,6.,15.);
//          rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),1.,0.8,3.);
//          rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),4.6,2.,8. );
//
//          rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),1.,0.,2.);
//          rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.2,0.8,3.);
//
//        }
//
//       RooDoubleCrystalBall* model_CB = new RooDoubleCrystalBall(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB,*rrv_alpha2_CB,*rrv_n2_CB);
//       return model_CB;
//     }
//
   //  if( model == "GausChebychev_ttbar_failtau2tau1cut_2"){
   //
   //    RooRealVar* rrv_mean1_gaus  = NULL;
   //    RooRealVar* rrv_sigma1_gaus = NULL;
   //
   //
   //    RooAbsPdf* model_pdf = NULL ;
   //
   //    RooAbsPdf* gaus = NULL ;
   //
   //    double p0_tmp = -4.73757e-01;
   //    double p1_tmp = 9.57687e-02;
   //    double frac_tmp = 2.90629e-01;
   //
   //      // take the same gaussian used in the pass sample
   //    if( TString(label).Contains("data") ) {
   //      rrv_mean1_gaus  = workspace->var(("rrv_mean_CB_data_"+channel+spectrum).c_str());
   //      rrv_sigma1_gaus = workspace->var(("rrv_sigma_CB_data_"+channel+spectrum).c_str());
   //    }
   //    else if( TString(label).Contains("TotalMC") ){
   //      rrv_mean1_gaus  = workspace->var(("rrv_mean_CB_TotalMC_"+channel+spectrum).c_str());
   //      rrv_sigma1_gaus = workspace->var(("rrv_sigma_CB_TotalMC_"+channel+spectrum).c_str());
   //    }
   //    else if( (TString(label).Contains("realW") or TString(label).Contains("fakeW")) ) {
   //      rrv_mean1_gaus  = workspace->var(("rrv_mean_CB_TTbar_realW_"+channel+spectrum).c_str());
   //      rrv_sigma1_gaus = workspace->var(("rrv_sigma_CB_TTbar_realW_"+channel+spectrum).c_str());
   //    }
   //
   //    gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
   //
   //    RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);
   //    RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);
   //    if( TString(label).Contains("realW") ) { // && TString(label).Contains("realW")
   //      p0_tmp = -3.03757e-01;
   //      p1_tmp = 1.307687e-01;
   //      rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);//,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
   //      rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);//,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
   //      frac_tmp = 2.90629e-02;
   //    }
   //
   //    RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );
   //
   //    RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,frac_tmp-frac_tmp*4,frac_tmp+frac_tmp*4);
   //    model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);
   //
   //    return model_pdf ;
   //  }
   //
   //  if(model == "Double_CB_Exp"){
   //    std::cout<< "########### Double Crystal Ball with exponential ############"<<std::endl;
   //
   //
   //    double c0_tmp = -4.18408e-02;     double c0_tmp_err = 3.06900e-03;
   //    RooRealVar* rrv_alpha_CB = NULL ;
   //    RooRealVar* rrv_n_CB = NULL ;
   //    RooRealVar* rrv_n2_CB = NULL ;
   //    RooRealVar* rrv_alpha2_CB = NULL ;
   //    RooRealVar* rrv_mean_CB  = NULL;
   //    RooRealVar* rrv_sigma_CB = NULL;
   //
   //    RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.01, c0_tmp-4e-1, c0_tmp+4e-1 );
   //
   //    RooRealVar* rrv_frac = NULL ;
   //
   //    rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),1.,0.0,100.);
   //    rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),4.6,0.,800. );
   //    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),1.,0.,1000.);
   //    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.2,0.,1000.);
   //
   //
   //    if( TString(label).Contains("data") and  workspace->var(("rrv_mean_CB_data_"+channel+spectrum).c_str()) ) {
   //      rrv_mean_CB  = workspace->var(("rrv_mean_CB_data_"+channel+spectrum).c_str());
   //      rrv_sigma_CB = workspace->var(("rrv_sigma_CB_data_"+channel+spectrum).c_str());
   //    }
   //    else if( TString(label).Contains("TotalMC") ){
   //      rrv_mean_CB  = workspace->var(("rrv_mean_CB_TotalMC_"+channel+spectrum).c_str());
   //      rrv_sigma_CB = workspace->var(("rrv_sigma_CB_TotalMC_"+channel+spectrum).c_str());
   //    }
   //    else if( (TString(label).Contains("realW") or TString(label).Contains("fakeW")) ) {
   //      rrv_mean_CB  = workspace->var(("rrv_mean_CB_TTbar_realW_"+channel+spectrum).c_str());
   //      rrv_sigma_CB = workspace->var(("rrv_sigma_CB_TTbar_realW_"+channel+spectrum).c_str());
   //    }
   //
   //    RooDoubleCrystalBall* model_CB = new RooDoubleCrystalBall(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB,*rrv_alpha2_CB,*rrv_n2_CB);
   //    RooExponential* model_Exp = new RooExponential(("model_Exp"+label+"_"+channel+spectrum).c_str(),("model_Exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
   //
   //    rrv_frac     = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.05,0.,1.);
   //
   //      // rrv_frac->setConstant(kTRUE);
   //
   //    RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), RooArgList(*model_Exp,*model_CB), RooArgList(*rrv_frac));
   //
   //    return model_pdf;
   // }
      

    // Gaus + pol 2 background for resonant part in the fail sample
    if( model == "GausChebychev_ttbar_failtau2tau1cut"){
    
      RooAbsPdf* model_pdf = NULL ;
  
      double p0_tmp = -3.73757e-01;     double p0_tmp_err = 1.73757e-01;     
      double p1_tmp = 7.57687e-02;      double p1_tmp_err = 4.73757e-02;       
      double frac_tmp = 2.90629e-01;   
      

      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      
        if( TString(label).Contains("data"   ) ) gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
        if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());
        if( TString(label).Contains("realW")   ) gaus = workspace->pdf(("gaus1_TTbar_realW_"+channel+spectrum).c_str());

      RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);
      RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);
      if( TString(label).Contains("realW")  ) { // && TString(label).Contains("realW")
        p0_tmp = -3.03757e-01;     
        p1_tmp = 1.307687e-01;     
        rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
        rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
        frac_tmp = 0.5;
      }
      if( TString(wtagger_label.c_str()).Contains("PuppiSD")){// and TString(label).Contains("realW") ){

        RooRealVar* rrv_mean1_gaus  = NULL;
        RooRealVar* rrv_sigma1_gaus = NULL;
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80, 80-5, 80+9);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),8.);//, 8-5,8+40 );
        gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

        rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);//,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
        rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);//,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
      }
      

       
      RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,0.,0.5);
      model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);

          
      return model_pdf ;
    }
    


    // Erf*Exp in the ttbar sample
    if( model == "ErfExp_ttbar"){

      double c0_tmp = -4.18408e-02;     double c0_tmp_err = 3.06900e-03;
      double offset_tmp = 8.27928e+01;  double offset_tmp_err = 3.20278e+00;
      double width_tmp = 2.94111e+01;   double width_tmp_err = 8.81193e-01;
      
      // if(  TString(wtagger_label.c_str()).Contains("herwig")  ){
      //
      //   c0_tmp         = -2.9893e-02 ;  c0_tmp_err     = 6.83e-03;
      //   offset_tmp     = 7.9350e+01 ;  offset_tmp_err = 9.35e+00;
      //   width_tmp      = 3.3083e+01 ;  width_tmp_err  = 2.97e+00;
      // }
      
      if(TString(wtagger_label.c_str()).Contains("76X")){
         c0_tmp = -4.18408e-02;      c0_tmp_err = 3.06900e-03;
         offset_tmp = 8.27928e+01;   offset_tmp_err = 3.20278e+00;
         width_tmp = 2.94111e+01;    width_tmp_err = 8.81193e-01;
      }
      
                                                        
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2,c0_tmp+4e-2);
      
     
      // if(TString(wtagger_label.c_str()).Contains("76X")) rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,width_tmp-10, width_tmp+10);
      
      
       // if( TString(wtagger_label.c_str()).Contains("PuppiSD") ) rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);

   
      RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

      RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
      RooGaussian* gaus2 = addConstraint(rrv_offset_ErfExp,rrv_offset_ErfExp->getVal(),offset_tmp_err,constraint);
      RooGaussian* gaus3 = addConstraint(rrv_width_ErfExp,rrv_width_ErfExp->getVal(),width_tmp_err,constraint);
      workspace->import(*gaus1);
      workspace->import(*gaus2);
      workspace->import(*gaus3);

      return model_pdf ;
    }
    

    if( model == "ErfExp_ttbar_failtau2tau1cut"){
      
      RooErfExpPdf* model_pdf = NULL ;
    
      double c0_tmp = -5.44643e-02;      double c0_tmp_err = 1.43399e-03;
      double offset_tmp = 6.594e+01;     double offset_tmp_err = 3.40396e+00;
      double width_tmp = 3.90e+01 ;      double width_tmp_err = 1.92162e+00;

      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-1, c0_tmp+4e-1);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);//, offset_tmp-offset_tmp_err*20,offset_tmp+offset_tmp_err*20);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp);//,width_tmp-width_tmp_err*20, width_tmp+width_tmp_err*20);
      
  
      
      if( TString(label).Contains("data") and !( TString(wtagger_label.c_str()).Contains("PuppiSD") ) ) {
        offset_tmp = 7.294e+01;
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);// offset_tmp-offset_tmp_err*20,offset_tmp+offset_tmp_err*20);
      }
      
      if( TString(label).Contains("TotalMC" ) ){
        offset_tmp = 8.294e+01;
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);//,offset_tmp-offset_tmp_err*20,offset_tmp+offset_tmp_err*20);
        rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp);//,width_tmp-width_tmp_err*20, width_tmp+width_tmp_err*20);
        if( TString(label.c_str()).Contains("fail") ) {
          rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp);//,-1, 1.);
          rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), 41.);//,0., 100.);
          rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),93.);//,0.,120.);
        }
      }
      
      model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
      workspace->import(*gaus1);
      
      return model_pdf ;
    }
    
    if(model == "2Voig"){
      
      RooRealVar* rrv_mean_voig = new RooRealVar(("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(),("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(),84,78,88);
      RooRealVar* rrv_shift_2Voig  = new RooRealVar(("rrv_shift_2Voig"+label+"_"+channel+spectrum).c_str(),("rrv_shift_2Voig"+label+"_"+channel+spectrum).c_str(),10.8026);
      RooFormulaVar* rrv_mean_shifted = new  RooFormulaVar(("rrv_mean_voig2"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean_voig,*rrv_shift_2Voig));

      RooRealVar* rrv_width_voig = new RooRealVar(("rrv_width_voig"+label+"_"+channel+spectrum).c_str(),("rrv_width_voig"+label+"_"+channel+spectrum).c_str(),16.,6,26);
      RooRealVar* rrv_sigma_voig = new RooRealVar(("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(),5.,0.,10.);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.8,0.5,1.);

      RooVoigtian* model_voig1 = new RooVoigtian(("model_voig1"+label+"_"+channel+spectrum).c_str(),("model_voig1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_voig,*rrv_width_voig,*rrv_sigma_voig);

      RooVoigtian* model_voig2 = new RooVoigtian(("model_voig2"+label+"_"+channel+spectrum).c_str(),("model_voig2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_shifted,*rrv_width_voig,*rrv_sigma_voig);

      RooAddPdf* model_pdf = new  RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), RooArgList(*model_voig1,*model_voig2),RooArgList(*rrv_frac));

      
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
  }
  
  // 0v60 STARTS HERE!!!!!!!!!!
  else{
    if ( model == "ExpGaus"){
      std::cout<< "########### Exp + Gaus for mj  fit  ############"<<std::endl;
      RooRealVar* rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),0.05,-0.2,0.2);
      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);
      RooRealVar* rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
      
      
      if( TString(label).Contains("VV_failtau2tau1cut") ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.2,0.2);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78.,88);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4.,10);
        rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
      }
      
      RooExponential* exp         = new RooExponential(("exp"+label+"_"+channel+spectrum).c_str(),("exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp,*gaus),RooArgList(*rrv_high));

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
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.026,-0.1, -1.e-4);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),60.,5.,120);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10.,100.);
      
      if( TString(label).Contains("VV_failtau2tau1cut")) {

        rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.006,-0.1,0.);
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),380.,360.,500.);
        rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,10.,100.);
      }
      if(TString(label).Contains("WJets0_failtau2tau1cut") ) {

        rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.01,-0.1,0.);
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),380.,360.,500.);
        rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),95.,70.,105.);
      }
      
      RooErfExpPdf* model_pdf       = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
            
      return model_pdf ;
    }     
   


    if( model == "ErfExp_ttbar_failtau2tau1cut"){

      double c0_tmp     = -1.0143e-01;  double c0_tmp_err = 1.46e-02;
      double offset_tmp = 2.7718e+02 ;
      double width_tmp  = 7.1891e+01 ;


      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp);
      RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      // constrain just the slope of the exponential
      RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
      workspace->import(*gaus1);
      return model_pdf ;
    }

    // Gaus + pol 2 background for resonant part in the fail sample
    if( model == "GausChebychev_ttbar_failtau2tau1cut"){

      double p0_tmp = -4.73757e-01;   double p0_tmp_err = 2.42754e-02;
      double p1_tmp = 9.57687e-02;    double p1_tmp_err = 2.93041e-02;
      double frac_tmp = 2.90629e-01;  //double frac_tmp_err = 1.86392e-02;


      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      if( TString(label).Contains("data") )    gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
      if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());


      RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp,p0_tmp-p0_tmp_err*7,p0_tmp+p0_tmp_err*7);
      RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
        
      RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );


      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);//,frac_tmp-frac_tmp_err*4,frac_tmp+frac_tmp_err*4);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    // Erf*Exp in the ttbar sample
    if( model == "ErfExp_ttbar"){

      double c0_tmp         = -2.9893e-02 ; double c0_tmp_err     = 6.83e-03;
      double offset_tmp     = 7.9350e+01 ;  double offset_tmp_err = 9.35e+00;
      double width_tmp      = 3.3083e+01 ;  double width_tmp_err  = 2.97e+00;

      // if( TString(label).Contains("herwig") and not TString(label).Contains("data")){
      //   c0_tmp     = -2.9357e-02 ;
      //   offset_tmp = 7.9350e+01 ;
      //   width_tmp  = 3.3216e+01 ;
      // }

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

    // 2Gaus per ttbar fit
    if( model == "2Gaus_ttbar"){

      double mean1_tmp      = 8.3141e+01;
      double deltamean_tmp  = 6.9129e+00;
      double sigma1_tmp     = 7.5145e+00;
      double scalesigma_tmp = 3.6819e+00;
      double frac_tmp       = 6.7125e-01;


      RooRealVar* rrv_mean1_gaus  = NULL;
      RooRealVar* rrv_sigma1_gaus = NULL;
      rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
      rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-2,sigma1_tmp+4 );



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
  }

  if( model == "GausExp_ttbar_failtau2tau1cut"){

    RooAbsPdf* gaus = NULL ;
    if( TString(label).Contains("data")    ) gaus = workspace->pdf(("gaus1_ttbar_data_"   +channel+spectrum).c_str());
    if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());

    RooRealVar* rrv_c_Exp = new RooRealVar( ("rrv_c_ExpGaus"+label+"_"+channel).c_str(), ("rrv_c_ExpGaus"+label+"_"+channel).c_str(),-0.05,-0.6,0.);
    
    if( TString(label).Contains("data") ){
      rrv_c_Exp = new RooRealVar( ("rrv_c_ExpGaus"+label+"_"+channel).c_str(), ("rrv_c_ExpGaus"+label+"_"+channel).c_str(),-0.001,-0.01,0.000);
    }

    RooExponential* exp   = new RooExponential(("ExpGaus"+label+"_"+channel+spectrum).c_str(),("ExpGaus"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
    
    RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel).c_str(), ("rrv_high"+label+"_"+channel).c_str(),0.5,0.,1.);

    RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*exp),RooArgList(*rrv_high),1);

    return model_pdf ;
  }
  
  if(model == "ErfExpGaus_sp"){
    RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,0.1);
    RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),26.,10,50.);
         
    RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80.,20.,105.);
    RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),8);//,2.,60);

    RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);

    RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
          
    RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
    RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
          
    return model_pdf ;
  }
  
}