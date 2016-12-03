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
  std::cout<<"########### Add to Constraint List some parameters  ############"<<std::endl;
  RooRealVar* rrv_x_mean  = new RooRealVar((std::string(rrv_x->GetName())+"_mean").c_str() ,(std::string(rrv_x->GetName())+"_mean").c_str() ,x_mean);
  RooRealVar* rrv_x_sigma = new RooRealVar((std::string(rrv_x->GetName())+"_sigma").c_str(),(std::string(rrv_x->GetName())+"_sigma").c_str(),x_sigma);
  RooGaussian* constrainpdf_x = new RooGaussian(("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),*rrv_x,*rrv_x_mean, *rrv_x_sigma);
  //import in the workspace and save the name of constraint pdf
  ConstraintsList->push_back(constrainpdf_x->GetName());
  return constrainpdf_x ;
}

//////////////////

// ---------------------------------------------                                                                                                                                    
RooAbsPdf* MakeModelTTbarControlSample(RooWorkspace* workspace ,const std::string & label, const std::string & modelName, const std::string & spectrum, const std::string & channel, const std::string & wtagger, const std::string & information, std::vector<std::string>* constraint){

  std::cout<<" "<<std::endl;
  std::cout<<"Defining number of real W-jets(rrv_number_total), W-tagging efficiency (eff_ttbar) and total yield after scaling with efficiency (rrv_number)"<<std::endl;
  std::cout<<" Making model for (label, modelName, information) : "<<label<<" "<<modelName<<" "<<information<<" "<<std::endl;
  std::cout<<""<<std::endl;
  std::cout<<" "<<std::endl;

  RooRealVar* rrv_number_total = NULL , *eff_ttbar = NULL ; 
  RooFormulaVar *rrv_number = NULL ;

  if(TString(label).Contains("_ttbar_data") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    // eff_ttbar        = new RooRealVar(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.99);
    eff_ttbar        = new RooRealVar(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_data"+information+"_"+channel+spectrum).c_str(),0.7,0.0,0.99);
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "@0*@1", RooArgList(*rrv_number_total,*eff_ttbar));
  }

  else if(TString(label).Contains("_ttbar_data") and TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = workspace->var(("rrv_number_total_ttbar_data"+information+"_"+channel+spectrum).c_str());
    eff_ttbar        = workspace->var(("eff_ttbar_data"+information+"_"+channel+spectrum).c_str());
    rrv_number       = new RooFormulaVar(("rrv_number"+label+"_"+channel+spectrum+spectrum).c_str(), "(1-@0)*@1", RooArgList(*eff_ttbar,*rrv_number_total));
  }
  else if(TString(label).Contains("_ttbar_TotalMC") and not TString(label).Contains("failtau2tau1cut")){
    rrv_number_total = new RooRealVar(("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("rrv_number_total_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),500,0.,1e7);
    // eff_ttbar        = new RooRealVar(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),0.7,0.3,0.99);
    eff_ttbar        = new RooRealVar(("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),("eff_ttbar_TotalMC"+information+"_"+channel+spectrum).c_str(),0.7,0.0,0.99);
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
  
  
  // Fits for data and MC (not matched tt)
  if( !(TString(label).Contains("realW")) and !(TString(label).Contains("fakeW")) ){

    if( model == "Exp"){
      std::cout<< "######### Exp = levelled exp funtion for W+jets mlvj ############" <<std::endl;
      RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.030, -2., 0.05);
      RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      return model_pdf ;
    }
    
    if( model == "Gaus"){

      RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,40,100);
      RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,0.,15);  
      RooGaussian* model_pdf       = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      return model_pdf ;
    }
    

    if( model == "ErfExp" ){
      std::cout<< "########### Erf*Exp for mj fit  ############"<<std::endl;
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.026);//,-0.05, 0.05);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),41.,0.,100);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,1.,100.);
      
      if(!TString(wtagger_label.c_str()).Contains("0v60")){
      if(TString(label).Contains("_WJets0") ) {
        rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0279,-0.5,0.);
        rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,60.,75.);
        rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),23.8,20.,30.);

        if(TString(label).Contains("fail") ) {
          rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0294,-0.05,0.05);
          rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),39.,30.,50.);
          rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),22.5,10.,30.);
        }
        if(TString(wtagger_label.c_str()).Contains("Puppi")){
          
          rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0279,-0.5,0.);
          rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,10.,75.);
          rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),23.8,20.,80.);
        
          if(TString(label).Contains("fail") ) {
            rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.0279,-0.5,0.);
            rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,10.,75.);
            rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),23.8,20.,80.);
          }
        }
      }
    }

      RooErfExpPdf* model_pdf       = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      return model_pdf ;
    }

    if( model == "ExpGaus"){

      RooRealVar* rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,70,90);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,40);
      RooRealVar* rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.,0.,1.);

      if( TString(label).Contains("_STop_failtau2tau1cut" ) ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.03,-0.5,0.5);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,60,120); //Too narrow limits here often lead to error!! eg max 80
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,60);
        
        if(TString(wtagger_label.c_str()).Contains("Puppi")){
          rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-1.7882e-02,-1.,0.);
          rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),8.4346e+01 ,75,89); 
          rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),1.1533e+01,4,12);
          rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),7.3428e-01,0.,1.);
        }
        
        if(TString(wtagger_label.c_str()).Contains("DDT")){
          rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),8.9813e-03);//,-5.,5.);
          rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,75,89); 
          rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,12);
        }
        
      }
      if( TString(label).Contains("_VV") ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.2,0.2);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),90,0.,150.);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4.,50.);
        
        if(TString(wtagger_label.c_str()).Contains("0v60")){
          rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-1.7882e-02,-1.,0.);
          rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),8.4346e+01 ,75,89); 
          rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),1.1533e+01,4,12);
          rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),7.3428e-01,0.,1.);
        }
        
        if(TString(wtagger_label.c_str()).Contains("DDT")){
          rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
          rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,70.,95.);
          rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4.,12.);   
        }
        if(TString(wtagger_label.c_str()).Contains("Puppi")){
          rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
          rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,70.,95.);
          rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,6.,200.);   
          rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.,0.0,1.);
          
          if( TString(label).Contains("fail") ) {
            rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
            rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,70.,95.);
            rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,6.,20.);
            rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.,0.6,1.);
          }
        }
        
      }
      RooExponential* exp         = new RooExponential(("exp"+label+"_"+channel+spectrum).c_str(),("exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp,*gaus),RooArgList(*rrv_high));
      return model_pdf ;
    }

    if( model == "ErfExpGaus_sp"){

      RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.04,-0.2,0.);
      RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,300.);
      RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,40,100);
      RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,0.,40);
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-4.0352e-02);//,-.,0.);
        rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,80.);
        rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,75,90);
        rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,5.,11);
        
        if( TString(label).Contains("_VV") ) {
          rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-4.0352e-02,-0.09,0.);
          rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,80.);
          rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,75,87);
          rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,5.,9.);
          
        }
      }
      

      RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);
      RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
      if( TString(label).Contains("_VV") )rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.8,0.5,1.);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));

      return model_pdf ;
    }

    if( model == "GausErfExp_ttbar"){

      double mean1_tmp = 8.2653e+01;
      double sigma1_tmp = 7.5932e+00;
      float rangeMean = 6. ;
      float rangeWidth = 5. ;
      double frac_tmp = 1.0;
      double c0_tmp     = -2.7180e-02 ;      //double c0_tmp_err     = 6.83e-03;
      double offset_tmp =  8.6888e+01 ;      //double offset_tmp_err = 9.35e+00;
      double width_tmp  =  2.9860e+01 ;      //double width_tmp_err  = 2.97e+00;  
      
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        mean1_tmp = 8.2402e+01;   
        sigma1_tmp = 7.5645e+00;            
        rangeMean = 6. ;
        rangeWidth = 5. ;
        c0_tmp     = -4.2672e-02 ;
        offset_tmp =  8.5656e+01 ;
        width_tmp  =  2.5308e+01  ;
      }
      if( TString(wtagger_label.c_str()).Contains("PuppiSD") ){
        mean1_tmp = 8.0857e+01;
        sigma1_tmp = 8.4035e+00;
        rangeMean = 5. ;
        rangeWidth = 5. ;
        c0_tmp     = -2.1046e-02  ;
        offset_tmp =  8.0122e+01 ;
        width_tmp  =  2.9595e+01 ;
          
        if(TString(wtagger_label.c_str()).Contains("0v56")){
          mean1_tmp = 8.2402e+01;
          sigma1_tmp = 7.5645e+00;
          rangeMean = 6. ;
          rangeWidth = 5. ;
          c0_tmp     = -4.2672e-02 ;
          offset_tmp =  8.5656e+01 ;
          width_tmp  =  2.5308e+01  ;
        }
      }
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        
        if(TString(wtagger_label.c_str()).Contains("0v38")){
          mean1_tmp = 8.6000e+01;
          sigma1_tmp = 7.2154e+00;
          c0_tmp     = -1.7193e-01 ;
          offset_tmp =  1.1819e+02 ;
          width_tmp  =   2.3273e+01  ;
        }
        
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          mean1_tmp = 86.;   
          sigma1_tmp = 1.0008e+01 ;            
          rangeMean = 6. ;
          rangeWidth = 5. ;
          c0_tmp     = -1.7193e-01 ;
          offset_tmp =  1.1819e+02 ;
          width_tmp  =  2.3273e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v44")){
          mean1_tmp = 8.6138e+01 ;   
          sigma1_tmp = 9.9019e+00 ;            
          rangeMean = 6. ;
          rangeWidth = 5. ;
          c0_tmp     = -1.2532e-01 ;
          offset_tmp =  1.1819e+02  ;
          width_tmp  =  2.3299e+00  ;
          // frac_tmp = 9.9757e-01
        }
      }

      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );
      RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      float offset_tmp_err = 20;;
      
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),width_tmp, width_tmp-10, width_tmp+10);
      
     
      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);
      

      rrv_c_ErfExp     ->setConstant(kTRUE);
      rrv_offset_ErfExp->setConstant(kTRUE);
      rrv_width_ErfExp ->setConstant(kTRUE);
      
      RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*erfExp),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "2Gaus_ttbar"){

      double mean1_tmp = 80.5;
      double sigma1_tmp = 8.1149e+00;
      float rangeMean = 8. ;
      float rangeWidth = 5. ;
      float frac_tmp = 8.3460e-01 ;
      float deltamean_tmp  = 9.499 ;
      float scalesigma_tmp = 2.5752;



      if(TString(wtagger_label.c_str()).Contains("76X")){
        frac_tmp = 7.3739e-01 ; 
        deltamean_tmp  = 7.81599999999998829e+00 ;
        scalesigma_tmp = 2.74416583258705149e+00;
      }
      if( TString(wtagger_label.c_str()).Contains("PuppiSD") ){
        mean1_tmp = 80.85;
        sigma1_tmp = 9.65;
        // frac_tmp = 6.9510e-01 ;
        deltamean_tmp  = 3.499 ;
        scalesigma_tmp = 2.5752;
      }
      
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        mean1_tmp = 8.5423e+01;   
        sigma1_tmp = 7.1694e+00 ;            
        rangeMean = 6. ;
        rangeWidth = 5. ;
        frac_tmp = 8.2613e-01 ;
        deltamean_tmp  = 4.4964e-13 ;
        scalesigma_tmp = 2.5898e+00;
      }
      

      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );

      if( TString(label.c_str()).Contains("fail") ){
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

      RooRealVar* rrv_deltamean_gaus  = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp);
      RooFormulaVar* rrv_mean2_gaus   = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
      RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp);
      RooFormulaVar* rrv_sigma2_gaus  = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
      RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);//,frac_tmp-frac_tmp_err,frac_tmp+frac_tmp_err);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "GausChebychev_ttbar_failtau2tau1cut"){

      RooAbsPdf* model_pdf = NULL ;

      double p0_tmp = 3.1099e-01 ; double p0_tmp_err = 1.86e-01;
      double p1_tmp = -2.2128e-01; double p1_tmp_err = 3.02e-01;
      double frac_tmp =  4.6400e-01  ; double frac_tmp_err = 1.20e-01;
      double mean1_tmp = 8.7682e+01;
      double sigma1_tmp = 9.2789e+00;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        p0_tmp = 2.1208e-01 ;
        p1_tmp = -3.2198e-01;
        frac_tmp = 4.3714e-01;
      }
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        if(TString(wtagger_label.c_str()).Contains("0v38")){
        p0_tmp = 3.7004e-01 ;
        p1_tmp = -4.5610e-01;
        frac_tmp = 6.8410e-01;
        mean1_tmp = 8.7682e+01;
        sigma1_tmp = 9.2789e+00;
      }
      }

      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      if( TString(label).Contains("data"   ) ) gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
      if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());
      if( TString(label).Contains("realW")   ) gaus = workspace->pdf(("gaus1_TTbar_realW_"+channel+spectrum).c_str());

      // RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, 75., 90.);
    //   RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp,5.,15. );
    //   RooGaussian* gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp);
      RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp);

      RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);
      model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "GausErfExp_ttbar_failtau2tau1cut"){

      double c0_tmp     = -3.0626e-02 ;      double c0_tmp_err = 1.46e-02;
      double offset_tmp = 5.1636e+01  ;      double offset_tmp_err = 4.92e+01;
      double width_tmp  = 3.6186e+01  ;      double width_tmp_err = 4.69e+00;
      double frac_tmp   = 1.          ;
      double mean1_tmp  = 8.7486e+01;
      double sigma1_tmp = 8.7456e+00;
      
      RooAbsPdf* model_pdf = NULL ;

      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp = -7.3294e-03    ;     
        offset_tmp = 3.7914e+01 ;    
        width_tmp = 2.6943e+01  ;
      }
      if( TString(wtagger_label.c_str()).Contains("PuppiSD") ){
        
        c0_tmp      = -5.5604e-02;
        offset_tmp  = 1.9915e+02 ;
        width_tmp   = 6.4888e+01 ;
        // frac_tmp = 4.5680e-01 ;
        mean1_tmp   = 8.0056e+01 ;
        sigma1_tmp  = 9.4160e+00 ;
          
        if(TString(wtagger_label.c_str()).Contains("0v56")){
          mean1_tmp = 8.4429e+01 ;
          sigma1_tmp = 7.8370e+00;
          c0_tmp     = -4.5251e-02 ;
          offset_tmp =  1.7294e+02 ;
          width_tmp  =  6.7105e+01  ;
          // frac_tmp = 2.5562e-01 ;
        }
      }
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        if(TString(wtagger_label.c_str()).Contains("0v38")){
          mean1_tmp = 8.7486e+01;
          sigma1_tmp = 8.7456e+00;
          c0_tmp     = -9.3552e-02 ;
          offset_tmp =  1.4301e+02 ;
          width_tmp  =   3.5879e+01  ;
          frac_tmp  = 8.3460e-01;
        }
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          mean1_tmp = 8.7107e+01;
          sigma1_tmp = 8.9305e+00;
          c0_tmp     =  -1.7229e-01 ;
          offset_tmp =  2.3192e+02 ;
          width_tmp  =  4.0557e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v44")){
          c0_tmp = -7.6941e-02    ;     
          offset_tmp = 1.2437e+02 ;    
          width_tmp = 3.3459e+01  ;
          // frac_tmp  = 5.5513e-01 ;
          mean1_tmp  = 8.7495e+01;
          sigma1_tmp = 8.8403e+00;
           
        }
      }
      
      
   
      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      if( TString(label).Contains("data"   ) ) gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
      if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());
      
      // RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, 75., 90.);
      //     // if( TString(label).Contains("data"   ) ) rrv_mean1_gaus = workspace->var( ("rrv_mean1_gaus_ttbar_data_"+channel+spectrum).c_str() );
      //     // if( TString(label).Contains("TotalMC") ) rrv_mean1_gaus = workspace->var( ("rrv_mean1_gaus_ttbar_TotalMC_"+channel+spectrum).c_str() );
      //     RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp,0.,sigma1_tmp+10 );
      //     RooGaussian* gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,width_tmp-width_tmp_err, width_tmp+width_tmp_err);

      RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp);
      
      rrv_c_ErfExp     ->setConstant(kTRUE);
      rrv_offset_ErfExp->setConstant(kTRUE);
      rrv_width_ErfExp ->setConstant(kTRUE);
      
      model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*erfExp),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "ErfExp_ttbar"){

      double c0_tmp     = -2.7180e-02 ;      double c0_tmp_err     = 6.83e-03;
      double offset_tmp =  8.6888e+01 ;      double offset_tmp_err = 9.35e+00;
      double width_tmp  =  2.9860e+01 ;      double width_tmp_err  = 2.97e+00;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        c0_tmp     = -1.9681e-02 ;
        offset_tmp =  7.4564e+01 ;
        width_tmp  =  2.6397e+01 ;
      }
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp     = -2.9373e-02 ;  
        offset_tmp =  8.2980e+01 ;
        width_tmp  =  3.6133e+01 ;
        if(TString(wtagger_label.c_str()).Contains("76X")){
          c0_tmp     = -2.2689e-02 ;
          offset_tmp =  6.9803e+01 ;
          width_tmp  =  3.2077e+01 ;
        }
      }
      
      if( TString(wtagger_label.c_str()).Contains("PuppiSD") ){
        c0_tmp     =  -2.1046e-02 ;     
        offset_tmp =  8.0122e+01 ;      
        width_tmp  =  2.9595e+01 ; 
      }
      
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        if(TString(wtagger_label.c_str()).Contains("0v38")){
          c0_tmp     = -6.0079e-02  ;
          offset_tmp =  1.4980e+02 ;
          width_tmp  =  4.8168e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          c0_tmp     =     -3.1793e-02  ;
          offset_tmp =  9.7997e+01 ;
          width_tmp  =  3.8820e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v44")){
          c0_tmp     = -4.1927e-02 ;
          offset_tmp =  1.1240e+02  ;
          width_tmp  =  3.9190e+01  ;
        }
      }
        
        
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),width_tmp, width_tmp-10, width_tmp+10);
      

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

      double c0_tmp = -3.0626e-02    ;      double c0_tmp_err = 1.46e-02;
      double offset_tmp = 5.1636e+01 ;      double offset_tmp_err = 4.92e+01;
      double width_tmp = 3.6186e+01  ;      double width_tmp_err = 4.69e+00;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        c0_tmp = -2.8654e-02    ;
        offset_tmp = 4.5140e+01 ;
        width_tmp = 3.4083e+01  ;
      }
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp = -3.4131e-02    ;  
        offset_tmp = 3.3961e+01 ;
        width_tmp = 3.2876e+01  ;
        if(TString(wtagger_label.c_str()).Contains("76X")){
          c0_tmp = -3.3221e-02    ; 
          offset_tmp = 2.9445e+01 ;
          width_tmp = 2.5680e+01  ;
        }
      }
      if( TString(wtagger_label.c_str()).Contains("PuppiSD") ){     
        c0_tmp = -3.5160e-02    ;       c0_tmp_err = 1.46e-02;
        offset_tmp = 6.0935e+01 ;       offset_tmp_err = 4.92e+01;
        width_tmp = 3.9583e+01  ;       width_tmp_err = 4.69e+00;
      }
      
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        if(TString(wtagger_label.c_str()).Contains("0v38")){
          c0_tmp     = -5.9405e-02 ;
          offset_tmp =  1.4790e+02 ;
          width_tmp  =  5.3544e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          c0_tmp     = -4.0722e-02  ;
          offset_tmp =  1.0108e+02 ;
          width_tmp  =  4.9796e+01  ;
        }
        if(TString(wtagger_label.c_str()).Contains("0v44")){
          c0_tmp     = -2.0670e-02 ;
          offset_tmp =  6.1299e+01  ;
          width_tmp  =  3.0247e+01  ;
          
        }
      }
         

      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp);
      
      model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
      workspace->import(*gaus1);
      return model_pdf ;
    }
  }


  // FOR MC FITS TO MATCHED TT MC!!!    
  else if( TString(label).Contains("realW") or TString(label).Contains("fakeW")){

    if( model == "Exp"){
      std::cout<< "######### Exp = levelled exp funtion for W+jets mlvj ############" <<std::endl;
      RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.030, -0.2, 0.05);
      RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      return model_pdf ;
    }

    if( model == "ErfExp" ){
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

    if( model == "ExpGaus"){

      RooRealVar* rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.5,0.5);
      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,70,90);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,40);
      RooRealVar* rrv_high        = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.,0.,1.);

      if( TString(label).Contains("_STop_failtau2tau1cut" ) ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.03,-0.5,0.5);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,60,120); //Too narrow limits here often lead to error!! eg max 80
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,60);
      }
      if( TString(label).Contains("_VV") ) {
        rrv_c_Exp       = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.2,0.2);
        rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),90,0.,150.);
        rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4.,50.);
      }
      RooExponential* exp         = new RooExponential(("exp"+label+"_"+channel+spectrum).c_str(),("exp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
      RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp,*gaus),RooArgList(*rrv_high));
      return model_pdf ;
    }

    if( model == "ErfExpGaus_sp"){

      RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.04,-0.2,0.);
      RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,300.);
      RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),80,40,100);
      RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,0.,40);

      RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);
      RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
      RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));

      return model_pdf ;
    }

    if( model == "GausErfExp_ttbar"){

      double mean1_tmp = 8.02653e+01;
      double sigma1_tmp = 7.5932e+00;
      float rangeMean = 6. ;
      float rangeWidth = 5. ;
      double frac_tmp = 0.6;
      double c0_tmp     = -2.7180e-02 ;      //double c0_tmp_err     = 6.83e-03;
      double offset_tmp =  8.6888e+01 ;      //double offset_tmp_err = 9.35e+00;
      double width_tmp  =  2.9860e+01 ;      //double width_tmp_err  = 2.97e+00;  
      
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        mean1_tmp = 8.2402e+01;   
        sigma1_tmp = 7.5645e+00;            
        rangeMean = 10. ;
        rangeWidth = 10. ;
        c0_tmp     = -4.2672e-02 ;
        offset_tmp =  8.5656e+01 ;
        width_tmp  =  2.5308e+01  ;
      }
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        mean1_tmp = 86.;   
        sigma1_tmp = 7.2154e+0;            
        rangeMean = 6. ;
        rangeWidth = 5. ;
        c0_tmp     = -1.7193e-01 ;
        offset_tmp =  1.1819e+02 ;
        width_tmp  =  2.3273e+01  ;
        
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          mean1_tmp = 8.7637e+01;   
          sigma1_tmp = 8.2736e+00;            
          rangeMean = 6. ;
          rangeWidth = 5. ;
          c0_tmp     = -2.5333e-01 ;
          offset_tmp =  1.6746e+02 ;
          width_tmp  =  2.6201e+01  ;
        }
      }

      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );
      RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,0.4,1);
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,-10,10.);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,0.,100.);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),width_tmp,0.,100.);


      
      RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*erfExp),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "2Gaus_ttbar"){

      double mean1_tmp = 8.5934e+01;
      double sigma1_tmp = 8.1149e+00;
      float rangeMean = 5. ;
      float rangeWidth = 5. ;
      double frac_tmp = 7.3994e-01 ;
      double deltamean_tmp  = 9.499 ;
      double scalesigma_tmp = 2.5752;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        frac_tmp = 7.3739e-01 ; 
        deltamean_tmp  = 7.81599999999998829e+00 ;
        scalesigma_tmp = 2.74416583258705149e+00;
      }

      RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean);
      RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth );

      if( TString(label.c_str()).Contains("fail") ){
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

      

      RooRealVar* rrv_deltamean_gaus  = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp,0,30);
      RooFormulaVar* rrv_mean2_gaus   = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
      RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp,0,5);
      RooFormulaVar* rrv_sigma2_gaus  = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
      RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,0,1);//,frac_tmp-frac_tmp_err,frac_tmp+frac_tmp_err);
      RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "GausChebychev_ttbar_failtau2tau1cut"){

      RooAbsPdf* model_pdf = NULL ;

      double p0_tmp = 3.1099e-01 ; double p0_tmp_err = 1.86e-01;
      double p1_tmp = -2.2128e-01; double p1_tmp_err = 3.02e-01;
      double frac_tmp =  4.6400e-01  ; double frac_tmp_err = 1.20e-01;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        p0_tmp = 2.1208e-01 ;
        p1_tmp = -3.2198e-01;
        frac_tmp = 4.3714e-01;
      }

      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      if( TString(label).Contains("data"   ) ) gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
      if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());
      if( TString(label).Contains("realW")   ) gaus = workspace->pdf(("gaus1_TTbar_realW_"+channel+spectrum).c_str());

      RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),p0_tmp,0.,1.);
      RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),p1_tmp,-1.,0.);//,p1_tmp-p1_tmp_err,p1_tmp+p1_tmp_err);

      RooChebychev* cheb = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(), *rrv_x, RooArgList(*rrv_p0_cheb,*rrv_p1_cheb) );

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,0.,1.);
      model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "GausErfExp_ttbar_failtau2tau1cut"){

      double c0_tmp     = -5.5276e-02;      double c0_tmp_err = 1.46e-02;
      double offset_tmp = 198. ;      double offset_tmp_err = 4.92e+01;
      double width_tmp  = 64.8  ;      double width_tmp_err = 4.69e+00;
      double frac_tmp   = 0.4;
      
      RooAbsPdf* model_pdf = NULL ;
      
      if(TString(wtagger_label.c_str()).Contains("0v40")){
        c0_tmp = -7.2072e-02    ;     
        offset_tmp = 1.8361e+02  ;    
        width_tmp = 5.3578e+01  ;
        frac_tmp = 4.2724e-01;
        
      }
      
      
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp = -7.3294e-03    ;     
        offset_tmp = 3.7914e+01 ;    
        width_tmp = 2.6943e+01  ;
      }
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        c0_tmp = -9.3552e-02    ;     
        offset_tmp = 1.4301e+02 ;    
        width_tmp = 3.5879e+01  ;
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          c0_tmp = -1.7229e-01    ;     
          offset_tmp = 2.3192e+02 ;    
          width_tmp = 4.0557e+01  ;
        }
      }
   
      // take the same gaussian used in the pass sample
      RooAbsPdf* gaus = NULL ;
      if( TString(label).Contains("data"   ) ) gaus = workspace->pdf(("gaus1_ttbar_data_"+channel+spectrum).c_str());
      if( TString(label).Contains("TotalMC") ) gaus = workspace->pdf(("gaus1_ttbar_TotalMC_"+channel+spectrum).c_str());
      if( TString(label).Contains("realW")   ) gaus = workspace->pdf(("gaus1_TTbar_realW_"+channel+spectrum).c_str());
      
      // double mean1_tmp  = 8.7486e+01;
      //     double sigma1_tmp = 8.7456e+00;
      //     RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-15, mean1_tmp+15);
      //     RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-3,sigma1_tmp+3 );
      //     RooGaussian* gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
      //
      

      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,-0.1,0.);//,c0_tmp-4e-2, c0_tmp+4e-2);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,0.,200.);//,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,30.,100);//,width_tmp-width_tmp_err, width_tmp+width_tmp_err);

      RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

      RooRealVar* rrv_frac = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp,0.,1.);
      
      model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*erfExp),RooArgList(*rrv_frac),1);

      return model_pdf ;
    }

    if( model == "ErfExp_ttbar"){

      double c0_tmp     = -2.7180e-02 ;      double c0_tmp_err     = 6.83e-03;
      double offset_tmp =  8.6888e+01 ;      double offset_tmp_err = 9.35e+00;
      double width_tmp  =  2.9860e+01 ;      double width_tmp_err  = 2.97e+00;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        c0_tmp     = -1.9681e-02 ;
        offset_tmp =  7.4564e+01 ;
        width_tmp  =  2.6397e+01 ;
      }
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp     = -2.9373e-02 ;  
        offset_tmp =  8.2980e+01 ;
        width_tmp  =  3.6133e+01 ;
        if(TString(wtagger_label.c_str()).Contains("76X")){
          c0_tmp     = -2.2689e-02 ;
          offset_tmp =  6.9803e+01 ;
          width_tmp  =  3.2077e+01 ;
        }
      }
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        c0_tmp     = -6.0079e-02 ;
        offset_tmp =  1.4980e+02 ;
        width_tmp  =  4.8168e+01 ;
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          c0_tmp = -3.1793e-02    ;     
          offset_tmp = 9.7997e+01 ;    
          width_tmp = 3.8820e+01  ;
        }
        
      }
      
      RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
      RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),width_tmp, width_tmp-10, width_tmp+10);

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

      double c0_tmp = -3.0626e-02    ;      double c0_tmp_err = 1.46e-02;
      double offset_tmp = 5.1636e+01 ;      double offset_tmp_err = 4.92e+01;
      double width_tmp = 3.6186e+01  ;      double width_tmp_err = 4.69e+00;

      if(TString(wtagger_label.c_str()).Contains("76X")){
        c0_tmp = -2.8654e-02    ;
        offset_tmp = 4.5140e+01 ;
        width_tmp = 3.4083e+01  ;
      }
      if(TString(wtagger_label.c_str()).Contains("0v60")){
        c0_tmp = -3.4131e-02    ;  
        offset_tmp = 3.3961e+01 ;
        width_tmp = 3.2876e+01  ;
        if(TString(wtagger_label.c_str()).Contains("76X")){
          c0_tmp = -3.3221e-02    ; 
          offset_tmp = 2.9445e+01 ;
          width_tmp = 2.5680e+01  ;
        }
      }
      if(TString(wtagger_label.c_str()).Contains("DDT")){
        c0_tmp     = -5.9405e-02 ;
        offset_tmp =  1.4790e+02 ;
        width_tmp  =  5.3544e+01 ;
        if(TString(wtagger_label.c_str()).Contains("0v52")){
          c0_tmp = -4.0722e-02    ;     
          offset_tmp = 1.0108e+02 ;    
          width_tmp = 4.9796e+01  ;
        }
      }

      RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2);
      RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),offset_tmp,0.,200.);
      RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,0,100.);
      
      model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
      RooGaussian* gaus1 = addConstraint(rrv_c_ErfExp,rrv_c_ErfExp->getVal(),c0_tmp_err,constraint);
      workspace->import(*gaus1);
      return model_pdf ;
    }
  
  }
} //End
