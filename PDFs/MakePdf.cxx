#include "MakePdf.h"

////////////////////////////////////////////
RooExtendPdf* MakeExtendedModel(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & channel, const std::string & wtagger_label, RooArgList* constraint, const int & ismc_wjet, const int & area_init_value){

      std::cout<<" "<<std::endl;
      std::cout<<"#########################################"<<std::endl;
      std::cout<<"## Make model : "<<label<<" "<<model<<"##"<<std::endl;
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


RooGaussian* addConstraint(RooRealVar* rrv_x, RooRealVar* x_mean, RooRealVar* x_sigma, std::vector<std::string> & ConstraintsList){

   //########### Gaussian contraint of a parameter of a pdf
   std::cout<<"########### Add to Contraint List some parameters  ############"<<std::endl;
   RooRealVar* rrv_x_mean = new RooRealVar((std::string(rrv_x->GetName())+"_mean").c_str(),(std::string(rrv_x->GetName())+"_mean").c_str(),x_mean->getVal());
   RooRealVar* rrv_x_sigma = new RooRealVar((std::string(rrv_x->GetName())+"_sigma").c_str(),(std::string(rrv_x->GetName())+"_sigma").c_str(),x_sigma->getVal());
   RooGaussian* constrainpdf_x = new RooGaussian(("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),("constrainpdf_"+std::string(rrv_x->GetName())).c_str(),*rrv_x,*rrv_x_mean, *rrv_x_sigma);
     //## import in the workspace and save the name of constriant pdf
   ConstraintsList.push_back(constrainpdf_x->GetName());
   return constrainpdf_x ;
}

//////////////////

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
       if(not param->isConstant()){
          param->Print();
          if ((param->getVal()-param->getMin())< param->getError()*1 or (param->getMax()- param->getVal())< param->getError()*1)
                    param->Print();
        }
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
RooAbsPdf* MakeGeneralPdf(RooWorkspace* workspace, const std::string & label, const std::string & model, const std::string & spectrum, const std::string & wtagger_label, const std::string & channel, RooArgList* constraint, const int & ismc){
  
 
  RooRealVar* rrv_x = NULL ; 
  if(TString(spectrum).Contains("_mj"))   rrv_x = workspace->var("rrv_mass_j");
  if(TString(spectrum).Contains("_mlvj")) rrv_x = workspace->var("rrv_mass_lvj");
  if(TString(spectrum).Contains("_genHMass")) rrv_x = workspace->var("rrv_mass_gen_WW");
    
   if(model == "Voig"){
            std::cout<<"########### Voigtian Pdf for mJ ############"<<std::endl;
            RooRealVar* rrv_mean_voig = new RooRealVar(("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(), ("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(), 84, 78, 88);
            RooRealVar* rrv_width_voig = new RooRealVar(("rrv_width_voig"+label+"_"+channel+spectrum).c_str(), ("rrv_width_voig"+label+"_"+channel+spectrum).c_str(), 7., 1, 40);
            RooRealVar* rrv_sigma_voig = new RooRealVar(("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(), ("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(), 5, 0.01, 20);
            RooVoigtian* model_pdf = new RooVoigtian(("model_pdf"+label+"_"+channel+spectrum).c_str(), ("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_voig,*rrv_width_voig,*rrv_sigma_voig);

            
            return model_pdf ;
   }
  
   if(model == "Voig_v1"){
            std::cout<<"########### Voigtian Pdf for Higgs mlvj ############"<<std::endl;
            RooRealVar* rrv_mean_voig = new RooRealVar(("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(), ("rrv_mean_voig"+label+"_"+channel+spectrum).c_str(),650,550,1200);
            RooRealVar* rrv_width_voig = new RooRealVar(("rrv_width_voig"+label+"_"+channel+spectrum).c_str(),("rrv_width_voig"+label+"_"+channel+spectrum).c_str(),100.,10,600);
            RooRealVar* rrv_sigma_voig = new RooRealVar(("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_voig"+label+"_"+channel+spectrum).c_str(),200,10,400);
            RooVoigtian* model_pdf = new RooVoigtian(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_voig,*rrv_width_voig,*rrv_sigma_voig);

            
            return model_pdf ;

   }

   if (model == "BWRUN"){
            RooRealVar* rrv_mean_BWRUN = NULL;
            RooRealVar* rrv_width_BWRUN = NULL;
 
            if( TString(label).Contains("H600")){                         
                rrv_mean_BWRUN  = new RooRealVar(("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),600,450,750);
                rrv_width_BWRUN = new RooRealVar(("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),60,5,300);
	    }	    
            else if( TString(label).Contains("H700")){                
                rrv_mean_BWRUN  = new RooRealVar(("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),700,550,850);
                rrv_width_BWRUN  = new RooRealVar(("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),50,5,350);
	    }
            else if(TString(label).Contains("H800")){                          
                rrv_mean_BWRUN  = new RooRealVar(("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),800,650,950);
                rrv_width_BWRUN = new RooRealVar(("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),150,5,400);
            }
            else if(TString(label).Contains("H900")){                          
                rrv_mean_BWRUN  = new RooRealVar(("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),900,700,1100);
                rrv_width_BWRUN = new RooRealVar(("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),400,70,500);
            }
            else if(TString(label).Contains("H1000")){                         
                rrv_mean_BWRUN  = new RooRealVar(("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BWRUN"+label+"_"+channel+spectrum).c_str(),1000,750,1250);
                rrv_width_BWRUN = new RooRealVar(("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),("rrv_width_BWRUN"+label+"_"+channel+spectrum).c_str(),100,2,570); 
            }

            RooBWRunPdf * model_pdf    = new RooBWRunPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_BWRUN,*rrv_width_BWRUN);
          
            //RooRealVar* rrv_mean_cb  = new RooRealVar(("rrv_mean_cb"+label+"_"+channel+spectrum).c_str(),("rrv_mean_cb"+label+"_"+channel+spectrum).c_str(),0);
            //RooRealVar* rrv_sigma_cb = new RooRealVar(("rrv_sigma_cb"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_cb"+label+"_"+channel+spectrum).c_str(),50,10,300);
            //RooGaussian* cbshape     = new RooGaussian(("cbshape"+label+"_"+channel+spectrum).c_str(),("cbshape"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_cb,*rrv_sigma_cb);
            //RooFFTConvPdf* model_pdf = new RooFFTConvPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*bwrun,*cbshape);
            
            return model_pdf ;
    }


    if(model == "CBBW_v1"){ //FFT: BreitWigner*CBShape                                                                                                                          

      // fix the breit wigner core from the generated higgs mass
     std::string mass_label = "" ;
     if(TString(label).Contains("ggH") and TString(label).Contains("600")) 
        mass_label = "ggH600";
     else if (TString(label).Contains("ggH") and TString(label).Contains("700"))
        mass_label = "ggH700";
     else if (TString(label).Contains("ggH") and TString(label).Contains("800"))
        mass_label = "ggH800";
     else if (TString(label).Contains("ggH") and TString(label).Contains("900"))
        mass_label = "ggH900";
     else if (TString(label).Contains("ggH") and TString(label).Contains("1000"))
        mass_label = "ggH1000";
     else if (TString(label).Contains("vbfH") and TString(label).Contains("600"))
        mass_label = "vbfH600";
     else if (TString(label).Contains("vbfH") and TString(label).Contains("700"))
        mass_label = "vbfH700";
     else if (TString(label).Contains("vbfH") and TString(label).Contains("800"))
        mass_label = "vbfH800";
     else if (TString(label).Contains("vbfH") and TString(label).Contains("900"))
        mass_label = "vbfH900";
     else if (TString(label).Contains("vbfH") and TString(label).Contains("1000"))
        mass_label = "vbfH1000";

     RooRealVar* rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),0);
     RooRealVar* rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),0);
     rrv_mean_BW->setConstant(kTRUE);
     rrv_width_BW->setConstant(kTRUE);

     RooBWRunPdf* bw          = new RooBWRunPdf(("bwrun_shape"+label+"_"+channel+spectrum).c_str(),("bwrun_shape"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_BW,*rrv_width_BW);

     RooRealVar* rrv_mean_CB  = NULL ; 
     RooRealVar* rrv_sigma_CB = NULL ; 
     RooRealVar* rrv_alpha_CB = NULL ; 
     RooRealVar* rrv_n_CB     = NULL ; 

     if(mass_label == "ggH600"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,0,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5,0.1,20);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());

     }
     else if(mass_label == "vbfH600"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,0,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5,0.1,20);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "ggH700"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),16,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,20,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-6,-8,-4);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),10,0.5,40);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "vbfH700"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),15,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,10,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),30,10,80);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "ggH800"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),16,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,10,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),10,0.1,30);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "vbfH800"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),16,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,20,100); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20,0.1,40);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "ggH900"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),10,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,30,120); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5,0.1,10);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "vbfH900"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),16,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,70,120); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20,0.1,40);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "ggH1000"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),10,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,30,120); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5,0.1,10);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else if(mass_label == "vbfH1000"){
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),16,-50,50);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,70,120); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20,0.1,40);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     else {
       rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),10,-20,20);
       rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,30,120); 
       rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-5,-8,-2);
       rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5,0.1,10);
       rrv_mean_BW->setVal(workspace->var(("rrv_mean_BWRUN"+label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
       rrv_width_BW->setVal(workspace->var(("rrv_width_BWRUN_"+mass_label+"BWRUN_"+channel+"_genHMass").c_str())->getVal());
     }
     //////
     std::string systematic_label ;
     if( TString(mass_label).Contains("ggH")) systematic_label = "_ggH" ;
     else if ( TString(mass_label).Contains("vbfH")) systematic_label = "_vbfH";

     RooRealVar* rrv_mean_scale_p1 = new RooRealVar(("CMS_sig_p1_jes"+systematic_label).c_str(),("CMS_sig_p1_jes"+systematic_label).c_str(),0);
     rrv_mean_scale_p1->setConstant(kTRUE);

     RooRealVar* rrv_mean_scale_p2 = new RooRealVar(("CMS_sig_p1_jer"+systematic_label).c_str(),("CMS_sig_p1_jer"+systematic_label).c_str(),0);
     rrv_mean_scale_p2->setConstant(kTRUE);

     SystematicUncertaintyHiggs_01jetBin systematic ;
     
     RooRealVar* rrv_mean_scale_X1 = new RooRealVar(("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),0);
     RooRealVar* rrv_mean_scale_X2 = new RooRealVar(("rrv_mean_shift_scale_jer"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jer"+label+"_"+channel+spectrum).c_str(),0);

     if (TString(label).Contains("ggH") and TString(label).Contains("600")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_600);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_600);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("700")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_700);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_700);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("800")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_800);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_800);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("900")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_900);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_900);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("1000")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_1000);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_1000);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("600")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_600);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_600);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("700")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_700);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_700);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("800")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_800);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_800);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("900")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_900);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_900);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("1000")){
       rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_1000);
       rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_1000);
     }

     rrv_mean_scale_X1->setConstant(kTRUE);
     rrv_mean_scale_X2->setConstant(kTRUE);

     RooFormulaVar* rrv_total_mean_CB = new RooFormulaVar(("rrv_total_mean_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_mean_CB,*rrv_mean_scale_p1,*rrv_mean_scale_X1,*rrv_mean_scale_p2,*rrv_mean_scale_X2));

     ////////////////////
            
     RooRealVar* rrv_sigma_scale_p1 = new RooRealVar(("CMS_sig_p2_jes"+systematic_label).c_str(),("CMS_sig_p2_jes"+systematic_label).c_str(),0);
     RooRealVar* rrv_sigma_scale_p2 = new RooRealVar(("CMS_sig_p2_jer"+systematic_label).c_str(),("CMS_sig_p2_jer"+systematic_label).c_str(),0);
     rrv_sigma_scale_p1->setConstant(kTRUE);
     rrv_sigma_scale_p2->setConstant(kTRUE);
 
     RooRealVar* rrv_sigma_scale_X1 = new RooRealVar(("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),0);
     RooRealVar* rrv_sigma_scale_X2 = new RooRealVar(("rrv_sigma_shift_jer"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_jer"+label+"_"+channel+spectrum).c_str(),0);

     if (TString(label).Contains("ggH") and TString(label).Contains("600")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_600);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_600);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("700")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_700);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_700);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("800")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_800);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_800);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("900")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_900);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_900);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("1000")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_1000);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_1000);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("600")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_600);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_600);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("700")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_700);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_700);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("800")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_800);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_800);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("900")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_900);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_900);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("1000")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_1000);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_1000);
     }

     rrv_sigma_scale_X1->setConstant(kTRUE);
     rrv_sigma_scale_X2->setConstant(kTRUE);

     RooFormulaVar* rrv_total_sigma_CB = new RooFormulaVar(("rrv_total_sigma_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_sigma_CB,*rrv_sigma_scale_p1,*rrv_sigma_scale_X1,*rrv_sigma_scale_p2,*rrv_sigma_scale_X2));        

     RooCBShape* cbshape  = new RooCBShape(("cbshape"+label+"_"+channel+spectrum).c_str(),("cbshape"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_total_mean_CB,*rrv_total_sigma_CB,*rrv_alpha_CB,*rrv_n_CB);

     rrv_x->setBins(5000);
     RooFFTConvPdf* model_pdf = new RooFFTConvPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*bw,*cbshape);
     model_pdf->setBufferFraction(2.0);

     return model_pdf;
   }

    //Voig for W mass peak
    if(model == "2Voig"){
 
            std::cout<< "########### Double Voigtian for mj fit ############"<<std::endl;
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

    // Gaus for the W peak
    if(model == "Gaus"){
            std::cout<< "########### Gaus for W peak  ############"<<std::endl;
            RooRealVar* rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
            RooRealVar* rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),7,1,15);
            RooGaussian* model_pdf = new  RooGaussian(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            
            return model_pdf ;
    }

    // Gaus for the higgs lineshape
    if(model == "Gaus_v1"){

            std::cout<< "########### Gaus for Higgs mlvj ############"<<std::endl;
            RooRealVar* rrv_mean_gaus = NULL ;
            RooRealVar* rrv_sigma_gaus = NULL ;
 
            if(label == "_ggH600_signal_region" or label == "_ggH600_sb_lo"){
                rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),580,550,620);
                rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),65,40,80);
            }
            if(label == "_ggH700_signal_region" or label == "_ggH700_sb_lo"){
                rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),700,650,750);
                rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),100,40,150);
            }
            if(label == "_ggH800_signal_region" or label == "_ggH800_sb_lo"){
                rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),800,750,850);
                rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),130,120,140);
            }
            if(label == "_ggH900_signal_region" or label=="_ggH900_sb_lo"){
                rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),900,850,900);
                rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),160,140,180);
            }
            if(label == "_ggH1000_signal_region" or label == "_ggH1000_sb_lo"){
                rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),920,900,1000);
                rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),200,100,300);
            }

            RooGaussian* model_pdf = new RooGaussian(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            
            return model_pdf ;
    }
 
    // Crystal Ball for the W mass peak
    if(model == "CB"){
            std::cout<< "########### Cystal Ball for mj fit ############"<<std::endl;
            RooRealVar* rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),84,78,88);
            RooRealVar* rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),7,4,10);
            RooRealVar* rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-2,-4,-0.5);
            RooRealVar* rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),2,0.,4);
            RooCBShape* model_pdf = new RooCBShape(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB);

            
            return model_pdf ;
     }

     // Sum of two CB 
     if(model == "SCB_v1"){
            std::cout<< "########### Cystal Ball + Crystall Ball ############"<<std::endl;
            RooRealVar* rrv_mean_SCB   = new RooRealVar(("rrv_mean_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_SCB"+label+"_"+channel+spectrum).c_str(),800,550,1000);
            RooRealVar* rrv_sigma_SCB  = new RooRealVar(("rrv_sigma_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_SCB"+label+"_"+channel+spectrum).c_str(),70,40,300);
            RooRealVar* rrv_alpha1_SCB = new RooRealVar(("rrv_alpha1_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_SCB"+label+"_"+channel+spectrum).c_str(),-2,-4,-0.5);
            RooRealVar* rrv_alpha2_SCB = new RooRealVar(("rrv_alpha2_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_SCB"+label+"_"+channel+spectrum).c_str(),2,0.5,4);
            RooRealVar* rrv_n1_SCB     = new RooRealVar(("rrv_n1_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_SCB"+label+"_"+channel+spectrum).c_str(),2,0.,4);
            RooRealVar* rrv_n2_SCB     = new RooRealVar(("rrv_n2_SCB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_SCB"+label+"_"+channel+spectrum).c_str(),2,0.,4);
            RooRealVar* frac           = new RooRealVar(("rrv_frac_SSCB"+label+"_"+channel+spectrum).c_str(),("rrv_frac_SSCB"+label+"_"+channel+spectrum).c_str(),0.5);
            RooCBShape* scb1 = new RooCBShape(("model_pdf_scb1"+label+"_"+channel+spectrum).c_str(),("model_pdf_scb1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_SCB,*rrv_sigma_SCB,*rrv_alpha1_SCB,*rrv_n1_SCB);
            RooCBShape* scb2 = new  RooCBShape(("model_pdf_scb2"+label+"_"+channel+spectrum).c_str(),("model_pdf_scb2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_SCB,*rrv_sigma_SCB,*rrv_alpha2_SCB,*rrv_n2_SCB);
            RooAddPdf* model_pdf = new  RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*scb1,*scb2),RooArgList(*frac));

            
            return model_pdf ;
     }
 
      // Crystal  ball shape for Bulk GR samples and higgs 
      if(model == "CB_v1"){
            std::cout<<"########### Crystal Ball for Higgs and  Bulk GR  mlvj ############"<<std::endl;

            RooRealVar* rrv_mean_CB = NULL ;
            RooRealVar* rrv_sigma_CB = NULL ;
            RooRealVar* rrv_alpha_CB = NULL ;
            RooRealVar* rrv_n_CB = NULL ;
 
            if(TString(label).Contains("H600")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),600,580,620);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),67,10,80);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,0.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20.,10,80 );
            }
            else if (TString(label).Contains("H700")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),700,650,750);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,40,150);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1,-3,-0.1);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20.,10,40);
            }
            else if (TString(label).Contains("ggH800")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),780,700,850);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),140,90,160);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,0.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5 , 2, 7);
            }
            else if (TString(label).Contains("vbfH800")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),800,750,850);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),140,10,160);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,0.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),5 , 2, 7);
            }
            else if (TString(label).Contains("ggH900")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),880,820,950);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),170,90,200);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,4.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(), 2., 0.5,5);
            }
            else if (TString(label).Contains("vbfH900")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),900,880,950);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),170,10,200);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,4.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(), 2., 0.5,5);
            }
            else if (TString(label).Contains("ggH1000")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),920,800,1150);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),200,10,300);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-7.,4.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),2.,0.5,4);
            }
            else if (TString(label).Contains("vbfH1000")){
                rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1000,800,1150);
                rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),200,10,300);
                rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-1.,-4.,4.);
                rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),2.,0.5,4);
            }
            else{
                if( TString(label).Contains("M600") and TString(label).Contains("BulkG_WW") and not          TString(label).Contains("M600_W")){
                    rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 600, 550, 650);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 30,10 ,80);
                                       
                }
                else if ( TString(label).Contains("M700") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M700_W")){
                    rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 700, 600, 800);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 30,10 ,80);
                                        
                }
                else if( TString(label).Contains("M800") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M800_W")){
                   rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 800, 600, 800);
                   rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 40,10 ,90);
                                        
                }
                else if (TString(label).Contains("M900") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M900_W") ){
                   rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 900, 600, 800);
                   rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 40,10 ,90);
                                        
                }
                else if (TString(label).Contains("M1000") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1000_W") ){
                   rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1000, 900,1100);
                   rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,20 ,120);
                                        
                }
                else if (TString(label).Contains("M1100") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1100_W") ){
                   rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1100,1000,1200);
                   rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,20 ,120);
                                        
                }
                else if (TString(label).Contains("M1200") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1200_W") ){
                   rrv_mean_CB = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1200,1100,1300);
                   rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,20 ,120);
                                       
                }
                else if (TString(label).Contains("M1300") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1300_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1300,1200,1400);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),65,45,120);
                         
                }
                else if (TString(label).Contains("M1400") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1400_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1400,1300,1500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,45,130);
                         
                }
                else if (TString(label).Contains("M1500") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1500_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1400,1600);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),75,50,145);
                                        
                }
                else if (TString(label).Contains("M1600") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1600_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1600,1500,1700);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),80,55,160);
                         
                }
                else if (TString(label).Contains("M1700") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1700_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1700,1600,1800);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),85,60,175);
                                            
                }
                else if (TString(label).Contains("M1800") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1800_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1800,1700,1900);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,65,190);
                      
                }
                else if (TString(label).Contains("M1900") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1900_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1900,1800,2000);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),95,70,200);
                       
                }
                else if (TString(label).Contains("M2000") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2000_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2000,1900,2100);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,75,210);
                       
                }
                else if (TString(label).Contains("M2100") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2100_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,2000,2200);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),105,80,225);
                      
                }
                else if (TString(label).Contains("M2200") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2200_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2200,2100,2300);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),110,85,245);
                      
                }
                else if (TString(label).Contains("M2300") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2300_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2300,2200,2400);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),115,90,250);
                        
                }
                else if (TString(label).Contains("M2400") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2400_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2400,2300,2500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),120,90,275);
                       
                }
                else if (TString(label).Contains("M2500") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2500_W") ){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2500,2400,2600);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),135,90,285);

                }
                else if(TString(label).Contains("M1000_W150") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1000,500,1500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 150,50 ,500);

                }
                else if (TString(label).Contains("M1000_W300") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1000,500,1500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 300,50 ,800);

                }
                else if (TString(label).Contains("M1000_W50") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1000,500,1500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,10 ,200);

                }
                else if (TString(label).Contains("M1500_W75") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 75,50 ,250);

                }
                else if (TString(label).Contains("M1500_W225") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new   RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 225,150 ,450);

                }
                else if (TString(label).Contains("M1500_W450") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 450,400 ,700);

                }
                else if (TString(label).Contains("M2100_W105") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,1500,2500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 105,90 ,300);

                }
                else if (TString(label).Contains("M2100_W315") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,1500,2500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 315,250 ,600);

                }
                else if (TString(label).Contains("M2100_W630") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,1500,2500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 630,550 ,900);

                }
                else{
                    rrv_mean_CB  = new  RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),700,550,2500);
                    rrv_sigma_CB = new  RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,20 ,120);
               }
            
            rrv_alpha_CB = new  RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),4,1,5);
            rrv_n_CB     = new   RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),20.,10,40);
            }                                                                            

            // experimental systematic uncertainty
           std::string systematic_label ;
          if( TString(label).Contains("ggH")) systematic_label = "_ggH" ;
          else if ( TString(label).Contains("vbfH")) systematic_label = "_vbfH";

          RooRealVar* rrv_mean_scale_p1 = new RooRealVar(("CMS_sig_p1_jes"+systematic_label).c_str(),("CMS_sig_p1_jes"+systematic_label).c_str(),0);
          rrv_mean_scale_p1->setConstant(kTRUE);

          RooRealVar* rrv_mean_scale_p2 = new RooRealVar(("CMS_sig_p1_jer"+systematic_label).c_str(),("CMS_sig_p1_jer"+systematic_label).c_str(),0);
          rrv_mean_scale_p2->setConstant(kTRUE);

          SystematicUncertaintyHiggs_01jetBin systematic ;
     
          RooRealVar* rrv_mean_scale_X1 = new RooRealVar(("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),0);
          RooRealVar* rrv_mean_scale_X2 = new RooRealVar(("rrv_mean_shift_scale_jer"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jer"+label+"_"+channel+spectrum).c_str(),0);

          if (TString(label).Contains("ggH") and TString(label).Contains("600")){
             rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_600);
             rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_600);
          }
          else if (TString(label).Contains("ggH") and TString(label).Contains("700")){
            rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_700);
            rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_700);
          }
          else if (TString(label).Contains("ggH") and TString(label).Contains("800")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_800);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_800);
          }
          else if (TString(label).Contains("ggH") and TString(label).Contains("900")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_900);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_900);
          }
          else if (TString(label).Contains("ggH") and TString(label).Contains("1000")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_ggH_1000);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_ggH_1000);
          }
          else if (TString(label).Contains("vbfH") and TString(label).Contains("600")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_600);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_600);
          }
          else if (TString(label).Contains("vbfH") and TString(label).Contains("700")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_700);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_700);
          }
          else if (TString(label).Contains("vbfH") and TString(label).Contains("800")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_800);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_800);
          }
          else if (TString(label).Contains("vbfH") and TString(label).Contains("900")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_900);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_900);
          }
          else if (TString(label).Contains("vbfH") and TString(label).Contains("1000")){
           rrv_mean_scale_X1->setVal(systematic.mean_signal_uncertainty_jet_scale_vbfH_1000);
           rrv_mean_scale_X2->setVal(systematic.mean_signal_uncertainty_jet_res_vbfH_1000);
          }

         rrv_mean_scale_X1->setConstant(kTRUE);
         rrv_mean_scale_X2->setConstant(kTRUE);

         RooFormulaVar* rrv_total_mean_CB = new RooFormulaVar(("rrv_total_mean_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_mean_CB,*rrv_mean_scale_p1,*rrv_mean_scale_X1,*rrv_mean_scale_p2,*rrv_mean_scale_X2));

     ////////////////////
            
     RooRealVar* rrv_sigma_scale_p1 = new RooRealVar(("CMS_sig_p2_jes"+systematic_label).c_str(),("CMS_sig_p2_jes"+systematic_label).c_str(),0);
     RooRealVar* rrv_sigma_scale_p2 = new RooRealVar(("CMS_sig_p2_jer"+systematic_label).c_str(),("CMS_sig_p2_jer"+systematic_label).c_str(),0);
     rrv_sigma_scale_p1->setConstant(kTRUE);
     rrv_sigma_scale_p2->setConstant(kTRUE);
 
     RooRealVar* rrv_sigma_scale_X1 = new RooRealVar(("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),0);
     RooRealVar* rrv_sigma_scale_X2 = new RooRealVar(("rrv_sigma_shift_jer"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_jer"+label+"_"+channel+spectrum).c_str(),0);

     if (TString(label).Contains("ggH") and TString(label).Contains("600")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_600);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_600);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("700")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_700);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_700);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("800")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_800);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_800);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("900")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_900);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_900);
     }
     else if (TString(label).Contains("ggH") and TString(label).Contains("1000")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_ggH_1000);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_ggH_1000);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("600")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_600);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_600);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("700")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_700);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_700);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("800")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_800);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_800);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("900")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_900);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_900);
     }
     else if (TString(label).Contains("vbfH") and TString(label).Contains("1000")){
       rrv_sigma_scale_X1->setVal(systematic.sigma_signal_uncertainty_jet_scale_vbfH_1000);
       rrv_sigma_scale_X2->setVal(systematic.sigma_signal_uncertainty_jet_res_vbfH_1000);
     }

     rrv_sigma_scale_X1->setConstant(kTRUE);
     rrv_sigma_scale_X2->setConstant(kTRUE);

     RooFormulaVar* rrv_total_sigma_CB = new RooFormulaVar(("rrv_total_sigma_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_sigma_CB,*rrv_sigma_scale_p1,*rrv_sigma_scale_X1,*rrv_sigma_scale_p2,*rrv_sigma_scale_X2));        
           
     RooCBShape* model_pdf = new RooCBShape(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_total_mean_CB,*rrv_total_sigma_CB,*rrv_alpha_CB,*rrv_n_CB);

            
      return model_pdf ;
            
    }

    // Crystal  ball shape for Bulk GR samples and higgs 
    if(model == "BWCB"){

            std::cout<<"########### Crystal Ball x Breit Wigner for Bulk Graviton width ############"<<std::endl;
            RooRealVar*  rrv_mean_BW = NULL ;
            RooRealVar*  rrv_width_BW = NULL ;
            RooRealVar*  rrv_mean_CB = NULL ;
            RooRealVar*  rrv_sigma_CB = NULL ;

            if(TString(label).Contains("M1000_W50") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(), 1000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 50);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,50.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,0,200);
	    }
            else if ( TString(label).Contains("M1000_W150") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(), 1000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 150);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,50.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,0,200);
            }
            else if ( TString(label).Contains("M1000_W300") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(), 1000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 300);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,50.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,0,200);
            }
            else if ( TString(label).Contains("M1500_W75") and TString(label).Contains("BulkG_WW")){ 
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 175,50 ,300);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,50.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,0,300);
            }
            else if ( TString(label).Contains("M1500_W225") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 225,150 ,450);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,80.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,0,250);
            }
            else if ( TString(label).Contains("M1500_W450") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 450,150 ,450);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,80.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,0,250);
            }
            else if ( TString(label).Contains("M2100_W105") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100,1000,2000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 105,150 ,450);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,80.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),150,0,250);
            }
            else if ( TString(label).Contains("M2100_W315") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100,1000,2000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 315,150 ,450);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,80.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),150,0,250);
            }
            else if ( TString(label).Contains("M2100_W630") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100,1000,2000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(), 630,150 ,450);
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,0.,80.);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),150,0,250);
	    }

            rrv_mean_BW->setConstant(kTRUE);
            rrv_width_BW->setConstant(kTRUE);
                    
            RooRealVar* rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),2,0,4);
	    RooRealVar* rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),1.,0.,4.);

            RooBreitWigner* bw      = new RooBreitWigner(("bw"+label+"_"+channel+spectrum).c_str(),("bw"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_BW,*rrv_width_BW);
            RooCBShape* cbshape     = new RooCBShape(("cbshape"+label+"_"+channel+spectrum).c_str(),("cbshape"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB);

            RooFFTConvPdf* model_pdf = new RooFFTConvPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x, *cbshape, *bw);

            
            return model_pdf ;
	    
	   }

      if(model == "ArgusBW_v1"){

   	   RooRealVar* rrv_width_BW = NULL ;
   	   RooRealVar* rrv_m0_Argus = NULL ;
   	   RooRealVar* rrv_c_Argus = NULL ;
   	   RooRealVar* rrv_frac = NULL ;

	    if(TString(label).Contains("ggH1000")){
	      rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),100,50,600);
                rrv_m0_Argus = new RooRealVar(("rrv_m0_Argus"+label+"_"+channel+spectrum).c_str(),("rrv_m0_Argus"+label+"_"+channel+spectrum).c_str(), 950 );
                rrv_c_Argus  = new RooRealVar(("rrv_c_Argus"+label+"_"+channel+spectrum).c_str(),("rrv_c_Argus"+label+"_"+channel+spectrum).c_str(),-1,-2,-1e-1);
                rrv_frac     = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.5,0.0,1.);
	    }
            else{
                rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),200,50,400);
                rrv_m0_Argus = new RooRealVar(("rrv_m0_Argus"+label+"_"+channel+spectrum).c_str(),("rrv_m0_Argus"+label+"_"+channel+spectrum).c_str(),1000);
                rrv_c_Argus  = new RooRealVar(("rrv_c_Argus"+label+"_"+channel+spectrum).c_str(),("rrv_c_Argus"+label+"_"+channel+spectrum).c_str(),-1,-2,0.1);
                rrv_frac     = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.5,0.0,1.);
	    }

            RooBreitWigner* bw    =  new RooBreitWigner(("bw"+label+"_"+channel+spectrum).c_str(),("bw"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_m0_Argus,*rrv_width_BW);
            RooArgusBG* argus     = new RooArgusBG(("argus"+label+"_"+channel+spectrum).c_str(),("argus"+label+"_"+channel+spectrum).c_str(), *rrv_x, *rrv_m0_Argus,*rrv_c_Argus);
            RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), RooArgList(*bw,*argus), RooArgList(*rrv_frac));

            
            return model_pdf ;
      }
	  
      if (model  == "CBBW"){
	    std::cout<<"########### Crystal Ball x Breit Wigner for W mass peak ############"<<std::endl;
            RooRealVar* rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),84.0,78,88);
            RooRealVar* rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),7,4,10);
            RooRealVar* rrv_alpha_CB = new RooRealVar(("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha_CB"+label+"_"+channel+spectrum).c_str(),-2,-4,-1);
            RooRealVar* rrv_n_CB     = new RooRealVar(("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n_CB"+label+"_"+channel+spectrum).c_str(),0.5,0.,2);
            RooRealVar* rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),0);
            RooRealVar* rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),10,5,20);
            RooCBShape* cbshape      = new RooCBShape(("cbshape"+label+"_"+channel+spectrum).c_str(),("cbshape"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_CB,*rrv_sigma_CB,*rrv_alpha_CB,*rrv_n_CB);
            RooBreitWigner* bw       = new RooBreitWigner(("bw"+label+"_"+channel+spectrum).c_str(),("bw"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean_BW,*rrv_width_BW);
            RooFFTConvPdf * model_pdf    = new RooFFTConvPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x, *cbshape, *bw);
            
            return model_pdf ;

      }
      
      
	    
      // Crystal  ball shape for Bulk GR samples and higgs 
      if(model == "DoubleCB_v1"){

	std::cout<<"########### Double CB for Bulk graviton mlvj ############"<<std::endl;
        RooRealVar* rrv_mean_CB = NULL ;
        RooRealVar* rrv_sigma_CB = NULL ;
        RooRealVar* rrv_n1_CB = NULL ;
        RooRealVar* rrv_n2_CB = NULL ;
        RooRealVar* rrv_alpha1_CB = NULL ;
        RooRealVar* rrv_alpha2_CB = NULL ;

	if ( TString(label).Contains("M600") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M600_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 600, 550, 650);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 30,10 ,80);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
	}             
        else if ( TString(label).Contains("M700") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M700_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(), 700, 600, 800);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 30,10 ,80);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }                                    
        else if ( TString(label).Contains("M800") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M800_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),820,790,880);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,40,70);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 15.,5.,25.);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.64,1.,1.9);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),15.,5.,25.);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,1.,1.9);
	}
        else if ( TString(label).Contains("M900") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M900_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),920,850,950);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),59,45,70);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 25.,2,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),25.,0.1,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.25,0.5,3.);
	}  
        else if ( TString(label).Contains("M1000") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1000_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1020,970,1070);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),55,40,65);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
        }                
        else if ( TString(label).Contains("M1100") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1100_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1120,1080,1150);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),65,55,75);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,25);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,25);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
	}
        else if ( TString(label).Contains("M1200") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1200_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1220,1200,1250);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),65,55,75);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,30);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,5.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,30);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,5.);
        } 
        else if ( TString(label).Contains("M1300") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1300_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1320,1300,1350);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),70,60,75);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.3,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.3,0.5,3.);
        }                 
        else if ( TString(label).Contains("M1400") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1400_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1420,1400,1440);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),77,65,85);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
        }
        else if ( TString(label).Contains("M1500") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1500_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1515,1500,1530);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),81,71,91);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 15.,0.01,25);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),15.,0.01,25);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
	}
        else if ( TString(label).Contains("M1600") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1600_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1620,1600,1640);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),81,70,90);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }                 
        else if ( TString(label).Contains("M1700") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1700_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1720,1700,1740);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,75,96);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }                                    
        else if ( TString(label).Contains("M1800") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1800_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1820,1800,1840);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,75,100);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }              
        else if ( TString(label).Contains("M1900") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M1900_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1920,1900,1940);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),95,80,115);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }               
        else if ( TString(label).Contains("M2000") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2000_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2020,2000,2040);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),100,80,115);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }               
        else if ( TString(label).Contains("M2100") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2100_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2120,2100,2140);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),105,85,115);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }              
        else if ( TString(label).Contains("M2200") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2200_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2220,2200,2250);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),115,75,140);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }              
        else if ( TString(label).Contains("M2300") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2300_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2320,2300,2340);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),115,95,120);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 15.,0.2,30);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),15.,0.2,20);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }                        
        else if ( TString(label).Contains("M2400") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2400_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2420,2400,2440);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),115,100,125);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
        }               
        else if ( TString(label).Contains("M2500") and TString(label).Contains("BulkG_WW") and not TString(label).Contains("M2500_W") ){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2520,2500,2540);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),125,90,145);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.);
	}
        else if ( TString(label).Contains("M1000_W150") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1020,970,1070);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),150,130,175);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
	}
        else if ( TString(label).Contains("M1000_W300") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1020,970,1070);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 300,50 ,800);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
        }
        else if ( TString(label).Contains("M1000_W50") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1020,970,1070);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),50,25,1000);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
        }
        else if ( TString(label).Contains("M1500_W75") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 75,50 ,250);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }
        else if ( TString(label).Contains("M1500_W225") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 225,150 ,450);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }
        else if ( TString(label).Contains("M1500_W450") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),1500,1000,2000);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 450,400 ,700);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }
        else if ( TString(label).Contains("M2100_W105") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,1500,2500);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 105,90 ,300);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }
        else if ( TString(label).Contains("M2100_W315") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,1500,2500);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 315,250 ,600);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
        }
        else if ( TString(label).Contains("M2100_W630") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),2100,2000,2200);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 630,500, 900);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);
	}
        else{
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),700,550,2500);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(), 50,20 ,120);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,35);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.,0.5,6.);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3,0.5,6.);

	}
         
        RooRealVar* rrv_mean_scale_p1 = new RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);
        RooRealVar* rrv_mean_scale_p2 = NULL ;
        rrv_mean_scale_p1->setConstant(kTRUE);
        if (channel+spectrum == "mu" ){             
	     rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_m","CMS_sig_p1_scale_m",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
	}
        else if( channel+spectrum == "el"){
	     rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_e","CMS_sig_p1_scale_e",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
	}
        else if( channel+spectrum == "em"){
   	     rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
	}

        SystematicUncertaintyEXO systematic ;

        RooRealVar* rrv_mean_scale_X1 = new RooRealVar(("rrv_mean_shift_scale_lep"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_lep"+label+"_"+channel+spectrum).c_str(),float(systematic.mean_signal_uncertainty_lep_scale));
        rrv_mean_scale_X1->setConstant(kTRUE);

        RooRealVar* rrv_mean_scale_X2 = new RooRealVar(("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),float(systematic.mean_signal_uncertainty_jet_scale));
        rrv_mean_scale_X2->setConstant(kTRUE);

        RooFormulaVar* rrv_total_mean_CB = new RooFormulaVar(("rrv_total_mean_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_mean_CB,*rrv_mean_scale_p1,*rrv_mean_scale_X1,*rrv_mean_scale_p2,*rrv_mean_scale_X2));
            
        RooRealVar* rrv_sigma_scale_p1 = NULL ;
        if(channel+spectrum == "mu"){
  	     rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_m","CMS_sig_p2_scale_m",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	}
        else if(channel+spectrum == "el"){
  	     rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_e","CMS_sig_p2_scale_e",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	}
        else if(channel+spectrum == "em"){
	     rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	}

        RooRealVar* rrv_sigma_scale_p2 = new RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
        RooRealVar* rrv_sigma_scale_p3 = new RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);
        rrv_sigma_scale_p2->setConstant(kTRUE);
        rrv_sigma_scale_p3->setConstant(kTRUE);

        RooRealVar* rrv_mean_sigma_X1 = new RooRealVar(("rrv_sigma_shift_lep_scale"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_scale"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_lep_scale));
        RooRealVar* rrv_mean_sigma_X2 = new RooRealVar(("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_scale"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_jet_scale));

        RooRealVar* rrv_mean_sigma_X3 = new RooRealVar(("rrv_sigma_shift_res"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_res"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_jet_res));
        rrv_mean_sigma_X1->setConstant(kTRUE);
        rrv_mean_sigma_X2->setConstant(kTRUE);
        rrv_mean_sigma_X3->setConstant(kTRUE);

        RooFormulaVar* rrv_total_sigma_CB = new RooFormulaVar(("rrv_total_sigma_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(*rrv_sigma_CB,*rrv_sigma_scale_p1,*rrv_mean_sigma_X1,*rrv_sigma_scale_p2,*rrv_mean_sigma_X2,*rrv_sigma_scale_p3,*rrv_mean_sigma_X3));        

        RooDoubleCrystalBall* model_pdf = new RooDoubleCrystalBall(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_total_mean_CB,*rrv_total_sigma_CB,*rrv_alpha1_CB,*rrv_n1_CB,*rrv_alpha2_CB,*rrv_n2_CB);

        
        return model_pdf ;
      }
       
      // Crystal  ball shape for Bulk GR samples and higgs 
      if (model == "BWDoubleCB"){

	std::cout<<"########### Double CB x BW for Bulk graviton width ############"<<std::endl;
        RooRealVar* rrv_mean_CB = NULL ;
        RooRealVar* rrv_sigma_CB = NULL ;
        RooRealVar* rrv_n1_CB = NULL ;
        RooRealVar* rrv_alpha2_CB = NULL ;
        RooRealVar* rrv_n2_CB = NULL ;
        RooRealVar* rrv_alpha1_CB = NULL ;
        RooRealVar* rrv_mean_BW = NULL ;
        RooRealVar* rrv_width_BW = NULL ;
     
        if(TString(label).Contains("M1000_W50") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB   = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB  = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),55,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.2,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.2,3.5);
                    rrv_mean_BW   = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1000);
                    rrv_width_BW  = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),50);
	}
       else if( TString(label).Contains("M1000_W150") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),55,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),150);
	}
        else if( TString(label).Contains("M1000_W300") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),55,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,35);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.5,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1000);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),300);
	}
        else if( TString(label).Contains("M1500_W75") and TString(label).Contains("BulkG_WW")){ 
                    rrv_mean_CB   = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB  = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),75,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.2,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.2,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),75);
        }           
        else if( TString(label).Contains("M1500_W225") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB   = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB  = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),75,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.2,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.2,3.5);
                    rrv_mean_BW   = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500);
                    rrv_width_BW  = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),225);
        }
        else if( TString(label).Contains("M1500_W450") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB   = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0,-100,100);
                    rrv_sigma_CB  = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),75,0,200);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.4,0.2,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.,0.2,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),1500);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),450);
	}
        else if( TString(label).Contains("M2100_W105") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,-100,100);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,20,250);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01, 45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),105);
	}
        else if( TString(label).Contains("M2100_W315") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,-100,100);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,20,250);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 10.,0.01,45);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,45);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),1.5,0.5,3.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),315);
	}
        else if( TString(label).Contains("M2100_W630") and TString(label).Contains("BulkG_WW")){
                    rrv_mean_CB  = new RooRealVar(("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),("rrv_mean_CB"+label+"_"+channel+spectrum).c_str(),0.,-100,100);
                    rrv_sigma_CB = new RooRealVar(("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_CB"+label+"_"+channel+spectrum).c_str(),90,20,250);
                    rrv_n1_CB     = new RooRealVar(("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n1_CB"+label+"_"+channel+spectrum).c_str(), 20.,0.01,105);
                    rrv_alpha2_CB = new RooRealVar(("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha2_CB"+label+"_"+channel+spectrum).c_str(),3.5,0.5,50.5);
                    rrv_n2_CB     = new RooRealVar(("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),("rrv_n2_CB"+label+"_"+channel+spectrum).c_str(),20.,0.01,105);
                    rrv_alpha1_CB = new RooRealVar(("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),("rrv_alpha1_CB"+label+"_"+channel+spectrum).c_str(),3.5,0.5,50.5);
                    rrv_mean_BW  = new RooRealVar(("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),("rrv_mean_BW"+label+"_"+channel+spectrum).c_str(),2100);
                    rrv_width_BW = new RooRealVar(("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),("rrv_width_BW"+label+"_"+channel+spectrum).c_str(),630);
	}
	// fix the Breit-Wigner core to the generated one  
        rrv_mean_BW->setConstant(kTRUE);
        rrv_width_BW->setConstant(kTRUE);                    
        RooBreitWigner* bw  = new  RooBreitWigner(("bw"+label+"_"+channel+spectrum).c_str(),("bw"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_BW,*rrv_width_BW);

        // Double Crystall ball term --> add parameters in order to do systematic on the signal shape inside combiner
        RooRealVar* rrv_mean_scale_p1 = new RooRealVar("CMS_sig_p1_jes","CMS_sig_p1_jes",0);  // jes effect on the mean
        rrv_mean_scale_p1->setConstant(kTRUE);

        SystematicUncertaintyEXO systematic;

        RooRealVar* rrv_mean_scale_p2 = NULL ;
        if( channel+spectrum == "mu" ){  /// lep scale effect on the mean       
             rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_m","CMS_sig_p1_scale_m",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
	}
        else if( channel+spectrum == "el" ){
             rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_e","CMS_sig_p1_scale_e",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
	}
        else if( channel+spectrum == "em"){
             rrv_mean_scale_p2 = new RooRealVar("CMS_sig_p1_scale_em","CMS_sig_p1_scale_em",0);
             rrv_mean_scale_p2->setConstant(kTRUE);
        }        
        //set the uncertainty value in other two independent variables 
        RooRealVar* rrv_mean_scale_X1 = new RooRealVar(("rrv_mean_shift_scale_lep"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_lep"+label+"_"+channel+spectrum).c_str(),float(systematic.mean_signal_uncertainty_lep_scale));
        rrv_mean_scale_X1->setConstant(kTRUE);
            
        RooRealVar* rrv_mean_scale_X2 = new RooRealVar(("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),("rrv_mean_shift_scale_jes"+label+"_"+channel+spectrum).c_str(),float(systematic.mean_signal_uncertainty_jet_scale));
        rrv_mean_scale_X2->setConstant(kTRUE);

	// total mean
        RooFormulaVar* rrv_total_mean_CB = new RooFormulaVar(("rrv_total_mean_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)", RooArgList(*rrv_mean_CB,*rrv_mean_scale_p1,*rrv_mean_scale_X1,*rrv_mean_scale_p2,*rrv_mean_scale_X2));

        ///lepton scale effect on the resolution 
        RooRealVar* rrv_sigma_scale_p1 = NULL ;
        if (channel+spectrum == "mu"){
             rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_m","CMS_sig_p2_scale_m",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	}
        else if( channel+spectrum == "el"){
             rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_e","CMS_sig_p2_scale_e",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	}
        else if( channel+spectrum == "em"){
             rrv_sigma_scale_p1 = new RooRealVar("CMS_sig_p2_scale_em","CMS_sig_p2_scale_em",0);
             rrv_sigma_scale_p1->setConstant(kTRUE);
	} 
	// Jes and jer effect on the resolution             
        RooRealVar* rrv_sigma_scale_p2 = new RooRealVar("CMS_sig_p2_jer","CMS_sig_p2_jer",0);
        RooRealVar* rrv_sigma_scale_p3 = new RooRealVar("CMS_sig_p2_jes","CMS_sig_p2_jes",0);

        rrv_sigma_scale_p2->setConstant(kTRUE);
        rrv_sigma_scale_p3->setConstant(kTRUE);

        RooRealVar* rrv_mean_sigma_X1 = new RooRealVar(("rrv_sigma_shift_lep_scale"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_scale"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_lep_scale));
        RooRealVar* rrv_mean_sigma_X2 = new RooRealVar(("rrv_sigma_shift_jes"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_scale"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_jet_scale));
        RooRealVar* rrv_mean_sigma_X3 = new RooRealVar(("rrv_sigma_shift_res"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_shift_res"+label+"_"+channel+spectrum).c_str(),float(systematic.sigma_signal_uncertainty_jet_res));

        rrv_mean_sigma_X1->setConstant(kTRUE);
        rrv_mean_sigma_X2->setConstant(kTRUE);
        rrv_mean_sigma_X3->setConstant(kTRUE);

        // total resolution 
        RooFormulaVar* rrv_total_sigma_CB = new RooFormulaVar(("rrv_total_sigma_CB"+label+"_"+channel+spectrum).c_str(),"@0*(1+@1*@2)*(1+@3*@4)*(1+@5*@6)", RooArgList(*rrv_sigma_CB,*rrv_sigma_scale_p1,*rrv_mean_sigma_X1,*rrv_sigma_scale_p2,*rrv_mean_sigma_X2,*rrv_sigma_scale_p3,*rrv_mean_sigma_X3));        

        RooDoubleCrystalBall* cbshape = new RooDoubleCrystalBall(("DoubleCB"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_total_mean_CB,*rrv_total_sigma_CB,*rrv_alpha1_CB,*rrv_n1_CB,*rrv_alpha2_CB,*rrv_n2_CB);

	/// numerical convolution via FFT
        RooFFTConvPdf* model_pdf = new RooFFTConvPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x,*bw,*cbshape);
        model_pdf->setBufferFraction(1.0);

        
        return model_pdf ;

      }
      
      // ExpN pdf for W+jets bkg fit
      if(model == "ExpN"){
            
	std::cout<< "########### ExpN funtion for W+jets mlvj ############"<<std::endl;
        RooRealVar* rrv_c_ExpN = new RooRealVar(("rrv_c_ExpN"+label+"_"+channel+spectrum).c_str(),("rrv_c_ExpN"+label+"_"+channel+spectrum).c_str(),-3e-3,-1e-1,-1e-5);
        RooRealVar* rrv_n_ExpN = NULL ;
        if(ismc==1){
	  rrv_n_ExpN = new RooRealVar(("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(),("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(), 1e3, -1e2, 1e4);
	}
	else {
	  if( channel+spectrum == "el" )
	    rrv_n_ExpN = new RooRealVar(("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(),("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(), 1e3, -1e2, 1e4);
	  else     
	    rrv_n_ExpN = new RooRealVar(("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(),("rrv_n_ExpN"+label+"_"+channel+spectrum).c_str(), 5e2, 0, 1e3);
	}
            
        RooExpNPdf* model_pdf = new RooExpNPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ExpN,*rrv_n_ExpN);
        
        return model_pdf ;

      }

      // levelled exp for W+jets bkg fit
      if( model == "ExpTail" ){
            
        std::cout<<"########### ExpTail = Levelled exp funtion for W+jets mlvj ############"<<std::endl;

        RooRealVar* rrv_s_ExpTail = new RooRealVar(("rrv_s_ExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_s_ExpTail"+label+"_"+channel+spectrum).c_str(), 250,-1.e6,1e6);
        RooRealVar* rrv_a_ExpTail = new RooRealVar(("rrv_a_ExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_a_ExpTail"+label+"_"+channel+spectrum).c_str(), 1e-1,-1.e-2,1e6);

        
	RooExpTailPdf* model_pdf = new RooExpTailPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_s_ExpTail,*rrv_a_ExpTail);
        
        return model_pdf ;

      }
 
      if(model == "Exp_v3"){

             RooRealVar* rrv_s_Exp = NULL ;
             RooRealVar* rrv_a_Exp = NULL ;
             RooRealVar* rrv_c_Exp = NULL ;

             if( TString(label).Contains("sb_lo")){
              rrv_s_Exp = new RooRealVar(("rrv_s_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_s_Exp_v3"+label+"_"+channel+spectrum).c_str(), -0.004,-0.2,0.);
              rrv_a_Exp = new RooRealVar(("rrv_a_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_a_Exp_v3"+label+"_"+channel+spectrum).c_str(), -10,-50.,-0.001);
              rrv_c_Exp = new RooRealVar(("rrv_c_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp_v3"+label+"_"+channel+spectrum).c_str(), -2.e-6,-0.001,0.001);
             }
             else{
 
              rrv_s_Exp = new RooRealVar(("rrv_s_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_s_Exp_v3"+label+"_"+channel+spectrum).c_str(), -0.001,-0.2,0.);
              rrv_a_Exp = new RooRealVar(("rrv_a_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_a_Exp_v3"+label+"_"+channel+spectrum).c_str(), -35,-50.,-0.001);
              rrv_c_Exp = new RooRealVar(("rrv_c_Exp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp_v3"+label+"_"+channel+spectrum).c_str(), -2.e-6,-0.001,0.001);
             }
             TString formula;
             formula.Form("TMath::Exp(%s+%s*%s+%s*%s*%s)",rrv_a_Exp->GetName(),rrv_x->GetName(),rrv_s_Exp->GetName(), rrv_c_Exp->GetName(),rrv_x->GetName(),rrv_x->GetName());
             RooGenericPdf* model_pdf = new RooGenericPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),"", formula.Data(), RooArgList(*rrv_x,*rrv_a_Exp,*rrv_s_Exp,*rrv_c_Exp) );

        
        return model_pdf ;

      }
      // two exponential --> exp^{-[0]*x}+[1]*exp^{-[2]*x}
      if (model == "2Exp"){

        RooRealVar* rrv_c1_Exp = new RooRealVar(("rrv_c1_2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_2Exp"+label+"_"+channel+spectrum).c_str(),-5e-4,-0.1,0.);
        RooRealVar* rrv_c2_Exp = new RooRealVar(("rrv_c2_2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c2_2Exp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.1,0.);

        RooExponential* exp1 = new RooExponential(("exp1"+label+"_"+channel+spectrum).c_str(),("exp1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c1_Exp);
        RooExponential* exp2 = new RooExponential(("exp2"+label+"_"+channel+spectrum).c_str(),("exp2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c2_Exp);

        RooRealVar* rrv_frac = new RooRealVar(("rrv_frac_2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_frac_2Exp"+label+"_"+channel+spectrum).c_str(),0.2,0.,8.);

        RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp1,*exp2),RooArgList(*rrv_frac),1);

        
        return model_pdf ;
      }

      // sum of two exponential 
      if(model == "Exp" or model == "Exp_sr"){
 	    std::cout<< "######### Exp = levelled exp funtion for W+jets mlvj ############" <<std::endl;
            RooRealVar* rrv_c_Exp = new RooRealVar(("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c_Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,0.);
            RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_Exp);
            
            return model_pdf ;
      }
      
      // Erf times for mj spectrum
      if ( model == "ErfExp" ){
   	    std::cout<< "########### Erf*Exp for mj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,-1e-4);
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),60.,30.,120);
            RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10, 60.);
            RooErfExpPdf* model_pdf       = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
            
           return model_pdf ;
      }
      // dif (ferent initial values -> for mlvj
      if ( model == "ErfExp_v1" ){
            std::cout<< "########### Erf*Exp for mlvj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp     = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.006,-0.1,0.);
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),450.,400.,550.);
            RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),70.,10,100.);
            RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);
            
            return model_pdf ;
      }
      // dif (ferent initial values -> for mlvj
        if ( model == "ErfExp_v2" ){ 
            std::cout<< "########### Erf*Exp for mlvj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.1,0.);
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),450.,400.,500.);
            RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), 50.,10,100.);
            RooRealVar* rrv_residue_ErfExp= new RooRealVar(("rrv_residue_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_residue_ErfExp"+label+"_"+channel+spectrum).c_str(),0.,0.,1.);

            TString formula ;
            formula.Form("(TMath::Exp(%s*%s) + %s)*(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_c_ErfExp->GetName(),rrv_x->GetName(),rrv_residue_ErfExp->GetName(),rrv_x->GetName(),rrv_offset_ErfExp->GetName(),rrv_width_ErfExp->GetName());

            RooGenericPdf* model_pdf = new RooGenericPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp,*rrv_residue_ErfExp));
            
            return model_pdf ;
	}

        // dif (ferent initial values -> for mlvj
        if ( model == "ErfExp_v3" ){ 
            std::cout<< "########### Erf*Exp for mlvj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.005,-0.1,0.);
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),450.,400,500.);
            RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), 50.,10,100.);
            RooRealVar* rrv_residue_ErfExp= new RooRealVar(("rrv_residue_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_residue_ErfExp"+label+"_"+channel+spectrum).c_str(),0.,0.,1.);
            RooRealVar* rrv_high_ErfExp    = new RooRealVar(("rrv_high_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_high_ErfExp"+label+"_"+channel+spectrum).c_str(),1.,0.,400);
            rrv_high_ErfExp->setConstant(kTRUE);
            TString formula ;
	    formula.Form("(TMath::Exp(%s*%s) + %s)* TMath::Power( ((1+TMath::Erf((%s-%s)/%s))/2.), %s )",rrv_c_ErfExp->GetName(),rrv_x->GetName(), rrv_residue_ErfExp->GetName(),rrv_x->GetName(),rrv_offset_ErfExp->GetName(), rrv_width_ErfExp->GetName(), rrv_high_ErfExp->GetName());
            RooGenericPdf* model_pdf = new RooGenericPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_high_ErfExp,*rrv_width_ErfExp,*rrv_residue_ErfExp) );
            
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
        // Erf*Exp + Gaus for mj spectrum 
        if(model == "ErfExpGaus"){
            std::cout<< "########### Erf*Exp + Gaus for mj  fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp     = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.4,0.);
            RooRealVar* rrv_offset_ErfExp= new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),100.,10.,300.);
            RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,100.);

            RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

            RooRealVar* rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),82,78,87);
            RooRealVar* rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);
            RooGaussian* gaus =  new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.7,0.,1.);
            RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
            return model_pdf ;
	}
        // Erf*Exp + Gaus for mj spectrum with offset == mean
        if(model == "ErfExpGaus_sp"){
 	    std::cout<< "########### Erf*Exp + Gaus for mj  fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.2,0.);
            RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,200.);
           
            RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
            RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);

            RooErfExpPdf* erfExp        = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean1_gaus,*rrv_width_ErfExp);

            RooGaussian* gaus            = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
            
            RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
            RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
            return model_pdf ;
	}
        // Erf*Exp+Gaus or mj spectrum
        if(model == "ErfExpGaus_v0"){
 	    std::cout<< "########### Erf*Exp + Gaus for mj  fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp     = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.2,0.);
            RooRealVar* rrv_offset_ErfExp= new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),100.,10.,140.);
            RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,100.);
            RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

            RooRealVar* rrv_mean_gaus     = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
            RooRealVar* rrv_sigma_gaus    = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);
            RooGaussian* gaus             = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            RooRealVar* rrv_high   = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.7,0.,1.);
            RooAddPdf* model_pdf  = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
            return model_pdf ;
	}
        // Erf*Exp+Gaus or mj spectrum
        if(model == "ErfExpGaus_v1"){
	  std::cout<< "########### Erf*Exp + Gaus for mlvj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.007,-0.1,0.);
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),800.,10.,1400.);
            RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),24.,10,150.);
            RooErfExpPdf* erfExp          = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

            RooRealVar* rrv_mean_gaus   = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),700,500,1200);
            RooRealVar* rrv_sigma_gaus  = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),150,10,300);
            RooGaussian* gaus           = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.1,0.,1.);
            RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
            return model_pdf ;

	}
        // Erf*Exp+Gaus or mj spectrum
        if(model == "ErfExpGaus_sp_v1"){
  	    std::cout<< "########### Erf*Exp + Gaus for mlvj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp    = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.007,-0.1,0.);
            RooRealVar* rrv_width_ErfExp= new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),24.,10,150.);
            RooRealVar* rrv_mean_gaus    = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),900,860,1200);
            RooErfExpPdf* erfExp         = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_mean_gaus,*rrv_width_ErfExp);

            RooRealVar* rrv_sigma_gaus   = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),150,10,300);
            RooGaussian* gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);

            RooRealVar* rrv_high  = new RooRealVar(("rrv_high"+label+"_"+channel+spectrum).c_str(),("rrv_high"+label+"_"+channel+spectrum).c_str(),0.1,0.,1.);
            RooAddPdf* model_pdf = new  RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus),RooArgList(*rrv_high));
            
            return model_pdf ;

	}

        // Erf*Exp + 2Gaus  
        if(model == "ErfExp2Gaus"){
 	    std::cout<< "########### Erf*Exp + 2Gaus for mj fit  ############"<<std::endl;
            RooRealVar* rrv_c_ErfExp     = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.2,0.);
            RooRealVar* rrv_offset_ErfExp= new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),100.,10.,140.);
            RooRealVar* rrv_width_ErfExp = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),30.,10,100.);
            RooErfExpPdf* erfExp = new RooErfExpPdf(("erfExp"+label+"_"+channel+spectrum).c_str(),("erfExp"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

            RooRealVar* rrv_mean1_gaus   = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),84,78,88);
            RooRealVar* rrv_mean2_gaus   = new RooRealVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),180,170,190);
            RooRealVar* rrv_sigma1_gaus  = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),7,4,10);
            RooRealVar* rrv_sigma2_gaus  = new RooRealVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),10,7,15);
            RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);
            RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

            RooRealVar* rrv_high1 = new RooRealVar(("rrv_high1"+label+"_"+channel+spectrum).c_str(),("rrv_high1"+label+"_"+channel+spectrum).c_str(),0.6,0.,1.);
            RooRealVar* rrv_high2 = new RooRealVar(("rrv_high2"+label+"_"+channel+spectrum).c_str(),("rrv_high2"+label+"_"+channel+spectrum).c_str(),0.4,0.,1.);
            RooAddPdf*  model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*erfExp,*gaus1,*gaus2),RooArgList(*rrv_high1,*rrv_high2));
            
            return model_pdf ;
	}
        // Gaus + Gaus for mj spectrum
        if(model == "2Gaus"){
	    std::cout<< "########### 2Gaus for mj fit  ############"<<std::endl;
            double mean1_tmp      = 8.3141e+01; 
            double deltamean_tmp  = 6.9129e+00; double deltamean_tmp_err  = 1.24e+00;
            double sigma1_tmp     = 7.5145e+00; 
            double scalesigma_tmp = 3.6819e+00; double scalesigma_tmp_err = 2.11e-01;
            double frac_tmp       = 6.7125e-01; double frac_tmp_err       = 2.09e-02;

            RooRealVar* rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

            RooRealVar* rrv_deltamean_gaus  = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp,deltamean_tmp-deltamean_tmp_err*4 ,deltamean_tmp+deltamean_tmp_err*4);
            RooFormulaVar* rrv_mean2_gaus   = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
            RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*8, scalesigma_tmp+scalesigma_tmp_err*8);
            RooFormulaVar* rrv_sigma2_gaus  = new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
            RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

            RooRealVar* rrv_frac  = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac),1);
            
            return model_pdf ;
	}
        // 2Gaus+2Gaus for VV mj spectrum -> WZ and WW
        if(model == "2_2Gaus"){

	    std::cout<< "########### 2Gaus +2Gaus for mj fit  ############"<<std::endl;
            double mean1_tmp      = 8.3141e+01; 
            double sigma1_tmp     = 7.5145e+00; 
            double scalesigma_tmp = 3.6819e+00; double scalesigma_tmp_err = 2.11e-01;
            double frac_tmp       = 6.7125e-01; double frac_tmp_err       = 2.09e-02;

            RooRealVar* rrv_shift = new RooRealVar(("rrv_shift"+label+"_"+channel+spectrum).c_str(),("rrv_shift"+label+"_"+channel+spectrum).c_str(),10.8026);

            RooRealVar* rrv_mean1_gaus = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
            RooRealVar* rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
            RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

            RooRealVar* rrv_deltamean_gaus     = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),0.,-8,10);
            RooFormulaVar* rrv_mean2_gaus      = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus, *rrv_deltamean_gaus));
            RooRealVar* rrv_scalesigma_gaus    = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp, scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
            RooFormulaVar* rrv_sigma2_gaus     =  new RooFormulaVar(("rrv_sigma2_gaus"+label+"_"+channel+spectrum).c_str(),"@0*@1", RooArgList(*rrv_sigma1_gaus,*rrv_scalesigma_gaus));
            RooGaussian* gaus2 = new RooGaussian(("gaus2"+label+"_"+channel+spectrum).c_str(),("gaus2"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean2_gaus,*rrv_sigma2_gaus);

            RooRealVar* rrv_frac1 = new RooRealVar(("rrv_frac1"+label+"_"+channel+spectrum).c_str(),("rrv_frac1"+label+"_"+channel+spectrum).c_str(),frac_tmp, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
            RooAddPdf* gausguas_1 = new RooAddPdf(("gausguas_1"+label+"_"+channel+spectrum).c_str(),("gausguas_1"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus1,*gaus2),RooArgList(*rrv_frac1),1);

            RooFormulaVar* rrv_mean3_gaus = new RooFormulaVar(("rrv_mean3_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_shift));
            RooFormulaVar* rrv_mean4_gaus = new RooFormulaVar(("rrv_mean4_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean2_gaus,*rrv_shift));
            RooGaussian* gaus3 = new RooGaussian(("gaus3"+label+"_"+channel+spectrum).c_str(),("gaus3"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean3_gaus,*rrv_sigma1_gaus);
            RooGaussian* gaus4 = new RooGaussian(("gaus4"+label+"_"+channel+spectrum).c_str(),("gaus4"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean4_gaus,*rrv_sigma2_gaus);
            RooAddPdf* gausguas_2 = new RooAddPdf(("gausguas_2"+label+"_"+channel+spectrum).c_str(),("gausguas_2"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus3,*gaus4),RooArgList(*rrv_frac1),1);

	    RooRealVar* rrv_frac  = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.74);//#,0.5,1.0)
	  	 
            RooAddPdf* model_pdf = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gausguas_1,*gausguas_2),RooArgList(*rrv_frac),1);
            
            return model_pdf ;
	}

        // Erf*Exp + 2Gaus for mj spectrum
        if(model == "2Gaus_ErfExp"){

	  std::cout<< "########### 2Gaus + Erf*Exp for mj fit  ############"<<std::endl;
          double mean1_tmp      = 8.3141e+01; 
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

        // User1 function 
        if(model == "User1"){
            std::cout<< "########### User 1 Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_p0     = new RooRealVar(("rrv_p0_User1"+label+"_"+channel+spectrum).c_str(),("rrv_p0_User1"+label+"_"+channel+spectrum).c_str(), 30, 10, 90);
            RooRealVar* rrv_p1 = NULL ;
            if (wtagger_label=="HP")
                rrv_p1 = new RooRealVar(("rrv_p1_User1"+label+"_"+channel+spectrum).c_str(),("rrv_p1_User1"+label+"_"+channel+spectrum).c_str(), -4, -9, -2);
            else
                rrv_p1 = new RooRealVar(("rrv_p1_User1"+label+"_"+channel+spectrum).c_str(),("rrv_p1_User1"+label+"_"+channel+spectrum).c_str(), -2, -4, 0.);
            RooUser1Pdf* model_pdf = new RooUser1Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_p0,*rrv_p1);
            
            return model_pdf ;
	}
        // QCD pdf  
        if (model == "QCD"){
            std::cout<< "########### QCD Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_p0 = new RooRealVar(("rrv_p0_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p0_QCD"+label+"_"+channel+spectrum).c_str(), 0,-200,200);
            RooRealVar* rrv_p1 = new RooRealVar(("rrv_p1_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p1_QCD"+label+"_"+channel+spectrum).c_str(), 0,-200,200);
            RooRealVar* rrv_p2 = new RooRealVar(("rrv_p2_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p2_QCD"+label+"_"+channel+spectrum).c_str(), 0,-200,200);
            RooQCDPdf*  model_pdf = new RooQCDPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_p0,*rrv_p1,*rrv_p2);
            
            return model_pdf ;
	}

        if(model == "QCD_v2"){//can replace exp
            std::cout<< "########### QCD Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_p0 = new RooRealVar(("rrv_p0_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p0_QCD"+label+"_"+channel+spectrum).c_str(), -15.,-50,0.);
            RooRealVar* rrv_p1 = new RooRealVar(("rrv_p1_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p1_QCD"+label+"_"+channel+spectrum).c_str(), 20,0,250);
            RooRealVar* rrv_p2 = new RooRealVar(("rrv_p2_QCD"+label+"_"+channel+spectrum).c_str(),("rrv_p2_QCD"+label+"_"+channel+spectrum).c_str(),0,-20,20);
            RooQCDPdf* model_pdf = new RooQCDPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_p0,*rrv_p1,*rrv_p2);
            
            return model_pdf ;
	}

        // For mlvj fit -> Pow function can replace exp
        if( model == "Pow" or model == "Pow_sr" ){
	    std::cout<< "########### Pow Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_c = new RooRealVar(("rrv_c_Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c_Pow"+label+"_"+channel+spectrum).c_str(), -5., -20., 0.);
            RooPowPdf* model_pdf = new RooPowPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c );
            
            return model_pdf ;
	}

        // For mlvj fit -> Pow function can replace exp
        if( model == "Pow2"){
            std::cout<< "########### Pow2 Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_Pow2"+label+"_"+channel+spectrum).c_str(),("rrv_c0_Pow2"+label+"_"+channel+spectrum).c_str(), 5, 0, 20);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_Pow2"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Pow2"+label+"_"+channel+spectrum).c_str(), 0, -5 , 5);
            RooPow2Pdf* model_pdf = new RooPow2Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(), *rrv_x, *rrv_c0, *rrv_c1);
            
            return model_pdf ;
	}

        // For mlvj fit ->Erf*Pow can replace Erf*Exp
        if( model == "ErfPow_v1"){
            std::cout<< "########### Erf*Pow Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_c      = new RooRealVar(("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(),-5.,-10.,0.);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(), 450,350,550);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),50,20,90);
            RooErfPowPdf* model_pdf  = new RooErfPowPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}

        if( model == "ErfPow"){
            std::cout<< "########### Erf*Pow Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_c      = new RooRealVar(("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(),-5.,-10.,0.);
            RooRealVar* rrv_offset_ErfPow = new RooRealVar(("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(),60.,30.,120);
            RooRealVar* rrv_width_ErfPow  = new RooRealVar(("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),30.,10, 60.);
            RooErfPowPdf* model_pdf  = new RooErfPowPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c,*rrv_offset_ErfPow,*rrv_width_ErfPow);
            
            return model_pdf ;
	}

        // For mlvj fit ->Erf*Pow can replace Erf*Exp -> in the sideband
        if( model == "ErfPow2_v1"){
            std::cout<< "########### Erf*Pow2 Pdf  for mlvj fit ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),14,1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(), 450,400,520);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),30,10,80);
            RooErfPow2Pdf* model_pdf  = new RooErfPow2Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        // For mlvj fit ->Erf*Pow can replace Erf*Exp for sr
        if( model == "ErfPow2_v1_sr"){
            std::cout<< "########### Erf*Pow2 Pdf  for mlvj fit in the SR  ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),4.,2.,8.);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(),-0.5,-2.,0.);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(), 490,440,520);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),50.,30.,80.);
            RooErfPow2Pdf* model_pdf = new RooErfPow2Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
             
            return model_pdf ;
	}
        // For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if( model == "ErfPowExp_v1"){
            std::cout<< "########### Erf*Pow*Exp Pdf  for mlvj fit   ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),11,5,20);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(), 0,-2,2);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(), 470,420,520);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),40,30,50);
            RooErfPowExpPdf* model_pdf  = new RooErfPowExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        // For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if( model == "ErfPowExp_v1_sr"){
            std::cout<< "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),6,2,15);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(), -1,-3,2);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(), 490,440,520);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),50,30,70);
            RooErfPowExpPdf* model_pdf = new RooErfPowExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}

        // For mlvj fit ->Erf*Pow*Exp can replace Erf*Exp 
        if( model == "ErfPowExp_v1_0"){//difference inital value
            std::cout<< "########### Erf*Pow*Exp Pdf for mlvj fit in SR  ############"<<std::endl;
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPowExp"+label+"_"+channel+spectrum).c_str(),20,15,40);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPowExp"+label+"_"+channel+spectrum).c_str(), 1.6,0.5,5);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPowExp"+label+"_"+channel+spectrum).c_str(), 470,420,520);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPowExp"+label+"_"+channel+spectrum).c_str(),47,30,60);
            RooErfPowExpPdf* model_pdf  = new RooErfPowExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        // Keys 
        if( model == "Keys"){
            std::cout<< "########### Erf*Pow*Exp Pdf for Keys  ############"<<std::endl;
            RooDataSet* rdataset  = (RooDataSet*) workspace->data("rdataset_signal_region_mlvj");
            RooKeysPdf* model_pdf = new RooKeysPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rdataset);
            
            return model_pdf ;
	}

        if( model == "GausChebychev"){

            RooRealVar* rrv_mean_gaus  = new RooRealVar(("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean_gaus"+label+"_"+channel+spectrum).c_str(),83.5, 80., 86.5);
            RooRealVar* rrv_sigma_gaus = new RooRealVar(("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma_gaus"+label+"_"+channel+spectrum).c_str(),7.6,3.5,30.);
            RooGaussian* gaus = new RooGaussian(("gaus"+label+"_"+channel+spectrum).c_str(),("gaus"+label+"_"+channel+spectrum).c_str(), *rrv_x,*rrv_mean_gaus,*rrv_sigma_gaus);
                
            RooRealVar* rrv_p0_cheb = new RooRealVar(("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p0_cheb"+label+"_"+channel+spectrum).c_str(), -0.3, -30., 30.);
            RooRealVar* rrv_p1_cheb = new RooRealVar(("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(),("rrv_p1_cheb"+label+"_"+channel+spectrum).c_str(), -0.1, -30., 30.);
            RooChebychev* cheb      = new RooChebychev(("cheb"+label+"_"+channel+spectrum).c_str(),("cheb"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0_cheb,*rrv_p1_cheb));

            RooRealVar* rrv_frac    = new RooRealVar(("rrv_frac"+label+"_"+channel+spectrum).c_str(),("rrv_frac"+label+"_"+channel+spectrum).c_str(),0.5, 0., 1.);
            RooAddPdf* model_pdf    = new RooAddPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),RooArgList(*gaus,*cheb),RooArgList(*rrv_frac),1);
            
            return model_pdf ;
	}
	// 2Gaus per ttbar fit
        if( model == "2Gaus_ttbar"){

            double mean1_tmp = 8.3141e+01;     
            double deltamean_tmp = 6.9129e+00;  
            double sigma1_tmp = 7.5145e+00;     
            double scalesigma_tmp = 3.6819e+00; 
            double frac_tmp = 6.7125e-01;       

            if( TString(label).Contains("herwig") and not TString(label).Contains("data")){

             mean1_tmp = 8.3141e+01; 
             deltamean_tmp = 6.9129e+00; 
             sigma1_tmp = 7.5145e+00; 
             scalesigma_tmp = 3.6819e+00; 
             frac_tmp = 6.7125e-01; 
	    }
           
            RooRealVar* rrv_mean1_gaus = NULL;
            RooRealVar* rrv_sigma1_gaus = NULL;
            // Constrain the peak and the width among electron and muon channel+spectrum to be the same in case of simultaneous fit
            if( channel+spectrum=="el"){                
	      if( workspace->var(("rrv_mean1_gaus%s_mu"+label).c_str()) and workspace->var(("rrv_sigma1_gaus%s_mu"+label).c_str())){
		    rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus%s_mu"+label).c_str());
                    rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus%s_mu"+label).c_str());
	      }
              else{
                    rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4);
	      }
              if(channel+spectrum=="mu" or channel+spectrum == "em"){
                if( workspace->var(("rrv_mean1_gaus%s_el"+label).c_str()) and workspace->var(("rrv_sigma1_gaus%s_el"+label).c_str())){
		  rrv_mean1_gaus  = workspace->var(("rrv_mean1_gaus%s_el"+label).c_str());
		  rrv_sigma1_gaus = workspace->var(("rrv_sigma1_gaus%s_el"+label).c_str());
		}
                else{
                    rrv_mean1_gaus  = new RooRealVar(("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_mean1_gaus"+label+"_"+channel+spectrum).c_str(),mean1_tmp, mean1_tmp-4, mean1_tmp+4);
                    rrv_sigma1_gaus = new RooRealVar(("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_sigma1_gaus"+label+"_"+channel+spectrum).c_str(),sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
		}
	      }
	    }

            RooGaussian* gaus1 = new RooGaussian(("gaus1"+label+"_"+channel+spectrum).c_str(),("gaus1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_mean1_gaus,*rrv_sigma1_gaus);

            RooRealVar* rrv_deltamean_gaus = new RooRealVar(("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_deltamean_gaus"+label+"_"+channel+spectrum).c_str(),deltamean_tmp);//,deltamean_tmp-deltamean_tmp_err*5, deltamean_tmp+deltamean_tmp_err*5);
            RooFormulaVar* rrv_mean2_gaus = new RooFormulaVar(("rrv_mean2_gaus"+label+"_"+channel+spectrum).c_str(),"@0+@1",RooArgList(*rrv_mean1_gaus,*rrv_deltamean_gaus));
            RooRealVar* rrv_scalesigma_gaus = new RooRealVar(("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),("rrv_scalesigma_gaus"+label+"_"+channel+spectrum).c_str(),scalesigma_tmp);//#,scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);

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

            double c0_tmp = -2.9893e-02 ;    
            double offset_tmp = 7.9350e+01 ; double offset_tmp_err = 9.35e+00;
            double width_tmp = 3.3083e+01 ;  

            if( TString(label).Contains("herwig") and not TString(label).Contains("data")){
              c0_tmp = -2.9357e-02 ; 
              offset_tmp = 7.9350e+01 ; 
              width_tmp = 3.3216e+01 ; 
	    }
                                                        
            RooRealVar* rrv_c_ErfExp      = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
            RooRealVar* rrv_offset_ErfExp = new RooRealVar(("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp"+label+"_"+channel+spectrum).c_str(), offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
            RooRealVar* rrv_width_ErfExp  = new RooRealVar(("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp"+label+"_"+channel+spectrum).c_str(), width_tmp,width_tmp-10, width_tmp+10);
                
            RooErfExpPdf* model_pdf = new RooErfExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp,*rrv_offset_ErfExp,*rrv_width_ErfExp);

            //Constrint the parameter within one sigma on the background part of the non resonant component in the pass sample
            //addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
            //addConstraint(rrv_offset_ErfExp,offset_tmp, offset_tmp_err,ConstraintsList);
            //addConstraint(rrv_width_ErfExp,width_tmp, width_tmp_err,ConstraintsList);
            
            return model_pdf ;
	}
        // Erf*Exp for failing events
        if( model == "ErfExp_ttbar_failtau2tau1cut"){

            double c0_tmp = -1.0143e-01 ;    
            double offset_tmp = 2.7718e+02 ; 
            double width_tmp = 7.1891e+01 ;  

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
            //addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
            
            return model_pdf ;
	}
        // extreme fail cut
        if( model == "Exp_ttbar_extremefailtau2tau1cut"){

            double c0_tmp = -3.0278e-02 ; 

            RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
            RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp);
            // constraint within one sigma
            //addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
            
            return model_pdf ;
	}

        if( model == "Exp_bkg_extremefailtau2tau1cut"){
           double c0_tmp = -4.2105e-02 ; 
           RooRealVar* rrv_c_ErfExp = new RooRealVar(("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp"+label+"_"+channel+spectrum).c_str(),c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
           RooExponential* model_pdf = new RooExponential(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_ErfExp);
           // constraint within one sigma
           //addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
            
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
	    
        // For mlvj fit -> Pow function [x/sqrt{s}]^{-[0]-[1]*log(x/sqrt(s)-[2]*log(x/sqrt(s))*log(x/sqrt(s))}                                             
        if( model == "Pow3"){

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_Pow3"+label+"_"+channel+spectrum).c_str(),("rrv_c0_Pow3"+label+"_"+channel+spectrum).c_str(), 5, 0, 20);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_Pow3"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Pow3"+label+"_"+channel+spectrum).c_str(), 0, -5 , 5);
            RooRealVar* rrv_c2 = new RooRealVar(("rrv_c2_Pow3"+label+"_"+channel+spectrum).c_str(),("rrv_c2_Pow3"+label+"_"+channel+spectrum).c_str(), 0, -5 , 5);

            RooPow3Pdf* model_pdf = new RooPow3Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_c2);
            
            return model_pdf ;
	}
	    // For mlvj fit -> 2Pow function){ x^{-[0]}+[1]*x^{-[1]}                                             
        if( model == "2Pow"){

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c0_2Pow"+label+"_"+channel+spectrum).c_str(), 5, 0, 20);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c1_2Pow"+label+"_"+channel+spectrum).c_str(), 0, -5 , 5);
            RooRealVar* rrv_frac = new RooRealVar(("rrv_frac_2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_frac_2Pow"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);
            TString formula;
            formula.Form("TMath::Power(%s,-%s)+%s*TMath::Power(%s,-%s)",rrv_x->GetName(),rrv_c0->GetName(),rrv_frac->GetName(),rrv_x->GetName(),rrv_c1->GetName());
  
            RooGenericPdf* model_pdf = new RooGenericPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_c0,*rrv_c1,*rrv_frac));
            
            return model_pdf ;

	}
        /// polynomial functions   
        if( model == "Chebychev_v2"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);

            RooChebychev* model_pdf = new RooChebychev(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1));
            
            return model_pdf ;
	}
        if( model == "Chebychev_v3"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
 
            RooChebychev* model_pdf = new RooChebychev(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));
            
            return model_pdf ;
	}
	
        if( model == "Chebychev_v4"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
 
            RooChebychev* model_pdf = new RooChebychev(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));
            
            return model_pdf ;
	}

        if( model == "Bernstein_v3"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), 0.3,  0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), 0.01, -1,1);
            RooRealVar*rrv_p2      = new RooRealVar(("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(), 0.01, 0,1);
 
            RooBernstein* model_pdf = new RooBernstein(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));
            
            return model_pdf ;
	}

        if( model == "Bernstein_v4"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
 
            RooBernstein* model_pdf = new RooBernstein(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));

            
            return model_pdf ;
	}
        if( model == "Bernstein_v5"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p0_Pol"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p1_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p2_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p3_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p4      = new RooRealVar(("rrv_p4_Pol"+label+"_"+channel+spectrum).c_str(),("rrv_p4_Pol"+label+"_"+channel+spectrum).c_str(), 0.01,0,1);
            RooBernstein* model_pdf = new RooBernstein(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3,*rrv_p4));
            
            return model_pdf ;

	}
	
	if( model == "ErfExpTail"){//different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExpTail"+label+"_"+channel+spectrum).c_str(),450.,400.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExpTail"+label+"_"+channel+spectrum).c_str(),70.,35.,100.);

            RooRealVar* rrv_s_ExpTail = NULL ;
            RooRealVar* rrv_a_ExpTail = NULL ;

            if(TString(label).Contains("signal_region")){
             rrv_s_ExpTail = new RooRealVar(("rrv_s_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_s_ErfExpTail"+label+"_"+channel+spectrum).c_str(), 500,0.,1e6);
             rrv_a_ExpTail = new RooRealVar(("rrv_a_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_a_ErfExpTail"+label+"_"+channel+spectrum).c_str(), -1e-1,-1.e2,1.);
	    }
	    else{
             rrv_s_ExpTail = new RooRealVar(("rrv_s_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_s_ErfExpTail"+label+"_"+channel+spectrum).c_str(), 500, 0.,1e6);
             rrv_a_ExpTail = new RooRealVar(("rrv_a_ErfExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_a_ErfExpTail"+label+"_"+channel+spectrum).c_str(), -1e-1,-1.e2,1.);
	    }

            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(),rrv_width_erf->GetName());
            RooGenericPdf* erf       = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf));

            RooExpTailPdf* expTail   = new RooExpTailPdf(("expTail"+label+"_"+channel+spectrum).c_str(),("expTail"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_s_ExpTail,*rrv_a_ExpTail);

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*expTail);
            
            return model_pdf ;
	}          
	if( model == "ErfExp_v3" ){ //different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),450.,400.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),60.,35.,100.);
            RooRealVar* rrv_s_ExpTail = NULL ;
            RooRealVar* rrv_a_ExpTail = NULL ;
            RooRealVar* rrv_c_ExpTail = NULL ;

            if(TString(label).Contains("sb_lo")){

             rrv_s_ExpTail = new RooRealVar(("rrv_s_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_s_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -0.004,-0.2,0.);
             rrv_a_ExpTail = new RooRealVar(("rrv_a_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_a_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -10,-50.,-0.001);
             rrv_c_ExpTail = new RooRealVar(("rrv_c_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -2.e-6,-0.001,0.001);
            }
            else{

             rrv_s_ExpTail = new RooRealVar(("rrv_s_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_s_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -0.001,-0.2,0.);
             rrv_a_ExpTail = new RooRealVar(("rrv_a_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_a_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -35,-50.,-0.001);
             rrv_c_ExpTail = new RooRealVar(("rrv_c_ErfExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfExp_v3"+label+"_"+channel+spectrum).c_str(), -2.e-6,-0.001,0.001);
	    }

            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(), rrv_width_erf->GetName());
	    RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf) );
             
	    formula.Form("TMath::Exp(%s+%s*%s+%s*%s*%s)",rrv_a_ExpTail->GetName(),rrv_s_ExpTail->GetName(), rrv_x->GetName(), rrv_c_ExpTail->GetName(),rrv_x->GetName(),rrv_x->GetName());
            RooGenericPdf* exp4 = new RooGenericPdf(("exp4"+label+"_"+channel+spectrum).c_str(),("exp4"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_s_ExpTail,*rrv_a_ExpTail,*rrv_c_ExpTail));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*exp4);

            
            return model_pdf ;
	}
	
        if(model == "Erf2Exp"){ //different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_Erf2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_Erf2Exp"+label+"_"+channel+spectrum).c_str(),450.,400.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_Erf2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_width_Erf2Exp"+label+"_"+channel+spectrum).c_str(),60.,40.,100.);

            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(),rrv_width_erf->GetName());
            RooGenericPdf* erf       = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf));

            RooRealVar* rrv_c1_Exp = new RooRealVar(("rrv_c1_Erf2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Erf2Exp"+label+"_"+channel+spectrum).c_str(),-0.006,-0.1,0.);
            RooRealVar* rrv_c2_Exp = new RooRealVar(("rrv_c2_Erf2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c2_Erf2Exp"+label+"_"+channel+spectrum).c_str(),-5e-4,-0.1,0.);

            RooRealVar* rrv_frac = new RooRealVar(("rrv_frac_Erf2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_frac_Erf2Exp"+label+"_"+channel+spectrum).c_str(),0.005,0.,0.05);

            formula.Form("TMath::Exp(%s*%s)+%s*TMath::Exp(%s*%s)",rrv_x->GetName(),rrv_c1_Exp->GetName(),rrv_frac->GetName(),rrv_x->GetName(),rrv_c2_Exp->GetName());
            RooGenericPdf* Exp = new RooGenericPdf(("2exp"+label+"_"+channel+spectrum).c_str(),("2exp"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_c1_Exp,*rrv_c2_Exp,*rrv_frac));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*Exp);
            
            return model_pdf ;
	}
	
        // ##################################################################
        if( model == "ErfPow_v1"){ //can replace erf*exp                                                                                                                           
  	    RooRealVar* rrv_c      = new RooRealVar(("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_c_ErfPow"+label+"_"+channel+spectrum).c_str(), -5.,-10.,0.);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow"+label+"_"+channel+spectrum).c_str(),50,20,100);

            RooErfPowPdf* model_pdf  = new RooErfPowPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}

        if( model == "ErfPow2_v1"){ //can replace erf*exp                                                                                                                           

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPow2"+label+"_"+channel+spectrum).c_str(),14,1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPow2"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(), 450,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),60,35,100);

            RooErfPow2Pdf* model_pdf  = new RooErfPow2Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        if( model == "ErfPow3_v1"){ //can replace erf*exp                                                                                                                           

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_ErfPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c0_ErfPow3"+label+"_"+channel+spectrum).c_str(),14,-1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_ErfPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c1_ErfPow3"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_c2 = new RooRealVar(("rrv_c2_ErfPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c2_ErfPow3"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_ErfPow3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow3"+label+"_"+channel+spectrum).c_str(), 450,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_ErfPow3"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow3"+label+"_"+channel+spectrum).c_str(),60,35,100);

            RooErfPow3Pdf* model_pdf  = new RooErfPow3Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_c2,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
	
        if( model == "Erf2Pow"){//can replace erf*exp                                                                                                                           
	
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfPow2"+label+"_"+channel+spectrum).c_str(), 450,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfPow2"+label+"_"+channel+spectrum).c_str(),60,35,100);
	
            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(),rrv_width->GetName()); 
            RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
       
            RooRealVar* rrv_c0 = NULL ;
            RooRealVar* rrv_c1 = NULL ;
            RooRealVar* rrv_frac = NULL ;
	
            if( TString(label).Contains("signal_region")){
             rrv_c0 = new RooRealVar(("rrv_c0_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c0_Erf2Pow"+label+"_"+channel+spectrum).c_str(),11,8.,13);
             rrv_c1 = new RooRealVar(("rrv_c1_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Erf2Pow"+label+"_"+channel+spectrum).c_str(),3,2,7);
             rrv_frac = new RooRealVar(("rrv_frac_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_frac_Erf2Pow"+label+"_"+channel+spectrum).c_str(),0.3,0.1,0.4);
	    }
	    else{
             rrv_c0 = new RooRealVar(("rrv_c0_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c0_Erf2Pow"+label+"_"+channel+spectrum).c_str(),3, 2.,7);
             rrv_c1 = new RooRealVar(("rrv_c1_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Erf2Pow"+label+"_"+channel+spectrum).c_str(),3, 2.,7);
             rrv_frac = new RooRealVar(("rrv_frac_Erf2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_frac_Erf2Pow"+label+"_"+channel+spectrum).c_str(),0.2,0.1,0.3);
	    }
		 
            formula.Form("TMath::Power(%s,-%s)+%s*TMath::Power(%s,-%s)",rrv_x->GetName(),rrv_c0->GetName(),rrv_frac->GetName(),rrv_x->GetName(),rrv_c1->GetName());
            RooGenericPdf* Pow = new RooGenericPdf(("2pow"+label+"_"+channel+spectrum).c_str(),("2pow"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_c0,*rrv_c1,*rrv_frac));

            RooProdPdf* model_pdf  = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*Pow);
            
            return model_pdf ;
	}
	  				
        //#################################################################

        if( model == "ErfChebychev_v2"){//can replace erf*exp                                                                                                           

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
            RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_offset,*rrv_width) );

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfChebychev_v2"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
 
            RooChebychev* pol = new RooChebychev(("Chebychev_v2"+label+"_"+channel+spectrum).c_str(),("Chebychev_v2"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);
            
            return model_pdf ;
	}
	
        if( model == "ErfChebychev_v3"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p2_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfChebychev_v3"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width) );
 
            RooChebychev *pol = new RooChebychev(("Chebychev_v3"+label+"_"+channel+spectrum).c_str(),("Chebychev_v3"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);
            
            return model_pdf ;
	}
	
        if( model == "ErfChebychev_v4"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p2_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p3_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfChebychev_v4"+label+"_"+channel+spectrum).c_str(),50,20,100);
            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooChebychev* pol = new RooChebychev(("Chebychev_v4"+label+"_"+channel+spectrum).c_str(),("Chebychev_v4"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);
            
            return model_pdf ;
       }

	//#################################################################

        if( model == "ErfBernstein_v3"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(), 1, -1., 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(), -0.01,-2.,2.);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p2_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(), 0.001,-1.,1.);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfBernstein_v3"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooBernstein* pol = new RooBernstein(("Bernstein_v3"+label+"_"+channel+spectrum).c_str(),("Bernstein_v3"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);
            
            return model_pdf ;
	}

        if( model == "ErfBernstein_v4"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p2_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p3_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfBernstein_v4"+label+"_"+channel+spectrum).c_str(),50,20,100);

 
            TString formula ;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* erf    = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_offset,*rrv_width));
            RooBernstein* pol     = new RooBernstein(("Bernstein_v4"+label+"_"+channel+spectrum).c_str(),("Bernstein_v4"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);

            
            return model_pdf ;
	}

        if( model == "ErfBernstein_v5"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p0_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p1_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p2_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p3_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,-1,1);
            RooRealVar* rrv_p4      = new RooRealVar(("rrv_p4_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p4_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_offset_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_width_ErfBernstein_v5"+label+"_"+channel+spectrum).c_str(),50,20,100);

 
            TString formula;
            formula.Form("(1.+TMath::Erf((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* erf = new RooGenericPdf(("erf"+label+"_"+channel+spectrum).c_str(),("erf"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width) );
            RooBernstein* pol = new RooBernstein(("Bernstein_v5"+label+"_"+channel+spectrum).c_str(),("Bernstein_v5"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3,*rrv_p4));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*erf,*pol);
            
            return model_pdf ;
	}
	
        //########################################################
        //############### Turn On + Falling part #################
        //########################################################
	
        if( model == "AtanExp_v1" ){ //different init-value and range                                                                                                            

            RooRealVar* rrv_c_AtanExp      = new RooRealVar(("rrv_c_AtanExp"+label+"_"+channel+spectrum).c_str(),("rrv_c_AtanExp"+label+"_"+channel+spectrum).c_str(),-0.006,-0.1,0.);
            RooRealVar* rrv_offset_AtanExp = new RooRealVar(("rrv_offset_AtanExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanExp"+label+"_"+channel+spectrum).c_str(),450.,350.,600.);
            RooRealVar* rrv_width_AtanExp  = new RooRealVar(("rrv_width_AtanExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanExp"+label+"_"+channel+spectrum).c_str(),60.,15.,100.);

            RooAtanExpPdf* model_pdf  = new RooAtanExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c_AtanExp,*rrv_offset_AtanExp,*rrv_width_AtanExp);
            
            return model_pdf ;
	}

        if( model == "AtanExpTail" ){ //different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_AtanExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanExpTail"+label+"_"+channel+spectrum).c_str(),450.,350.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_AtanExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanExpTail"+label+"_"+channel+spectrum).c_str(),70.,15.,100.);

            RooRealVar* rrv_s_ExpTail = new RooRealVar(("rrv_s_AtanExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_s_AtanExpTail"+label+"_"+channel+spectrum).c_str(), 250,-1.e6,1e6);
            RooRealVar* rrv_a_ExpTail = new RooRealVar(("rrv_a_AtanExpTail"+label+"_"+channel+spectrum).c_str(),("rrv_a_AtanExpTail"+label+"_"+channel+spectrum).c_str(), 5e-1,-1.e2,1e6);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(), rrv_width_erf->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf));

            RooExpTailPdf* expTail   = new RooExpTailPdf(("expTail"+label+"_"+channel+spectrum).c_str(),("expTail"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_s_ExpTail,*rrv_a_ExpTail);

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*expTail);
            
            return model_pdf ;
	}
        if( model == "AtanExp_v3" ){ //different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),450.,350.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),70.,15.,100.);

            RooRealVar* rrv_s_ExpTail = new RooRealVar(("rrv_s_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_s_AtanExp_v3"+label+"_"+channel+spectrum).c_str(), -0.006,-0.2,0.);

            RooRealVar* rrv_a_ExpTail = new RooRealVar(("rrv_a_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_a_AtanExp_v3"+label+"_"+channel+spectrum).c_str(), -0.1,-50.,-0.001);
            RooRealVar* rrv_c_ExpTail = new RooRealVar(("rrv_c_AtanExp_v3"+label+"_"+channel+spectrum).c_str(),("rrv_c_AtanExp_v3"+label+"_"+channel+spectrum).c_str(), -0.00001,-0.001,0.);

            TString formula; 
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(), rrv_width_erf->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf));

            formula.Form("TMath::Exp(%s+%s*%s+%s*%s*%s)",rrv_a_ExpTail->GetName(),rrv_s_ExpTail->GetName(), rrv_x->GetName(), rrv_c_ExpTail->GetName(),rrv_x->GetName(),rrv_x->GetName());
	    RooGenericPdf*exp4 = new RooGenericPdf(("exp4"+label+"_"+channel+spectrum).c_str(),("exp4"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_s_ExpTail,*rrv_a_ExpTail,*rrv_c_ExpTail));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*exp4);
            
            return model_pdf ;
	}

        if( model == "Atan2Exp" ){ //different init-value and range                                                                                                             
   
            RooRealVar* rrv_offset_erf = new RooRealVar(("rrv_offset_Atan2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_Atan2Exp"+label+"_"+channel+spectrum).c_str(),450.,350.,600.);
            RooRealVar* rrv_width_erf  = new RooRealVar(("rrv_width_Atan2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_width_Atan2Exp"+label+"_"+channel+spectrum).c_str(),70.,15.,100.);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset_erf->GetName(), rrv_width_erf->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset_erf,*rrv_width_erf));
 
            RooRealVar* rrv_c1_Exp = new RooRealVar(("rrv_c1_Atan2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Atan2Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,0.);
            RooRealVar* rrv_c2_Exp = new RooRealVar(("rrv_c2_Atan2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_c2_Atan2Exp"+label+"_"+channel+spectrum).c_str(),-0.05,-0.1,0.);

            RooExponential* exp1 = new RooExponential(("exp1"+label+"_"+channel+spectrum).c_str(),("exp1"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c1_Exp);
            RooExponential* exp2 = new RooExponential(("exp2"+label+"_"+channel+spectrum).c_str(),("exp2"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c2_Exp);

            RooRealVar* rrv_frac = new RooRealVar(("rrv_frac_Atan2Exp"+label+"_"+channel+spectrum).c_str(),("rrv_frac_Atan2Exp"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);

            RooAddPdf* Exp = new RooAddPdf(("2Exp"+label+"_"+channel+spectrum).c_str(),("2Exp"+label+"_"+channel+spectrum).c_str(),RooArgList(*exp1,*exp2),RooArgList(*rrv_frac));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*Exp);
            
            return model_pdf ;
	}
        //##################################################################

        if( model == "AtanPow_v1"){//can replace erf*exp                                                                                 
                                          
            RooRealVar* rrv_c      = new RooRealVar(("rrv_c_AtanPow"+label+"_"+channel+spectrum).c_str(),("rrv_c_AtanPow"+label+"_"+channel+spectrum).c_str(), -5.,-10.,0.);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_AtanPow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanPow"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_AtanPow"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanPow"+label+"_"+channel+spectrum).c_str(),50,20,100);

            RooAtanPowPdf* model_pdf  = new RooAtanPowPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        if( model == "AtanPow2_v1"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_AtanPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c0_AtanPow2"+label+"_"+channel+spectrum).c_str(),14,1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_AtanPow2"+label+"_"+channel+spectrum).c_str(),("rrv_c1_AtanPow2"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_AtanPow2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanPow2"+label+"_"+channel+spectrum).c_str(), 600,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_AtanPow2"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanPow2"+label+"_"+channel+spectrum).c_str(),60,10,100);

            RooAtanPow2Pdf* model_pdf  = new RooAtanPow2Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        if( model == "AtanPow3_v1"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_AtanPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c0_AtanPow3"+label+"_"+channel+spectrum).c_str(),14,1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_AtanPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c1_AtanPow3"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_c2 = new RooRealVar(("rrv_c2_AtanPow3"+label+"_"+channel+spectrum).c_str(),("rrv_c2_AtanPow3"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_AtanPow3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanPow3"+label+"_"+channel+spectrum).c_str(), 600,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_AtanPow3"+label+"_"+channel+spectrum).c_str(), ("rrv_width_AtanPow3"+label+"_"+channel+spectrum).c_str(),60,10,100);

            RooAtanPow3Pdf* model_pdf  = new RooAtanPow3Pdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_c2,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        if( model == "Atan2Pow"){//can replace erf*exp                                                                                                                           

            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_Atan2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_offset_Atan2Pow"+label+"_"+channel+spectrum).c_str(), 600,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_Atan2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_width_Atan2Pow"+label+"_"+channel+spectrum).c_str(),60,10,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_Atan2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c0_Atan2Pow"+label+"_"+channel+spectrum).c_str(),14,1,30);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_Atan2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_c1_Atan2Pow"+label+"_"+channel+spectrum).c_str(), 5,-5,10);
            RooRealVar* rrv_frac = new RooRealVar(("rrv_frac_Atan2Pow"+label+"_"+channel+spectrum).c_str(),("rrv_frac_Atan2Pow"+label+"_"+channel+spectrum).c_str(),0.5,0.,1.);

            formula.Form("TMath::Pow(%s,-%s)+%s*TMath::Pow(%s,-%s)",rrv_x->GetName(),rrv_c0->GetName(),rrv_frac->GetName(),rrv_x->GetName(),rrv_c1->GetName());  
	    RooGenericPdf* Pow = new RooGenericPdf(("2pow"+label+"_"+channel+spectrum).c_str(),("2pow"+label+"_"+channel+spectrum).c_str(),formula.Data(),RooArgList(*rrv_x,*rrv_c0,*rrv_c1,*rrv_frac));

            RooProdPdf* model_pdf  = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*Pow);
            
            return model_pdf ;
	}
        //###################################################################

        if( model == "AtanPowExp_v1"){//can replace erf*exp                                                                                                                         

            RooRealVar* rrv_c0 = new RooRealVar(("rrv_c0_AtanPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c0_AtanPowExp"+label+"_"+channel+spectrum).c_str(),13,5,40);
            RooRealVar* rrv_c1 = new RooRealVar(("rrv_c1_AtanPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_c1_AtanPowExp"+label+"_"+channel+spectrum).c_str(), 2,0,4);
            RooRealVar* rrv_offset = new RooRealVar(("rrv_offset_AtanPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanPowExp"+label+"_"+channel+spectrum).c_str(), 450,400,600);
            RooRealVar* rrv_width  = new RooRealVar(("rrv_width_AtanPowExp"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanPowExp"+label+"_"+channel+spectrum).c_str(),50,15,150);

            RooAtanPowExpPdf* model_pdf  = new RooAtanPowExpPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*rrv_x,*rrv_c0,*rrv_c1,*rrv_offset,*rrv_width);
            
            return model_pdf ;
	}
        //#################################################################

        if( model == "AtanChebychev_v2"){//#can replace erf*exp                                                                                                           

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanChebychev_v2"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
 
            RooChebychev* pol = new RooChebychev(("Chebychev_v2"+label+"_"+channel+spectrum).c_str(),("Chebychev_v2"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}
	
        if( model == "AtanChebychev_v3"){//can replace erf*exp                                                                                               

	    RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p2_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanChebychev_v3"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooChebychev* pol = new RooChebychev(("Chebychev_v3"+label+"_"+channel+spectrum).c_str(),("Chebychev_v3"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}
        if( model == "AtanChebychev_v4"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.3, -3, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p2_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p3_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(), -0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanChebychev_v4"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooChebychev *pol = new RooChebychev(("Chebychev_v4"+label+"_"+channel+spectrum).c_str(),("Chebychev_v4"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}
        //#################################################################

        if( model == "AtanBernstein_v3"){//#can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_p2_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooBernstein* pol = new RooBernstein(("Bernstein_v3"+label+"_"+channel+spectrum).c_str(),("Bernstein_v3"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}
        if( model == "AtanBernstein_v4"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);

            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p2_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_p3_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(), 0.001,-1,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanBernstein_v4"+label+"_"+channel+spectrum).c_str(),50,20,100);


            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooBernstein* pol = new RooBernstein(("Bernstein_v4"+label+"_"+channel+spectrum).c_str(),("Bernstein_v4"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}
        if( model == "AtanBernstein_v5"){//can replace erf*exp                                                                                               

            RooRealVar* rrv_p0      = new RooRealVar(("rrv_p0_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p0_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.3, 0, 3);
            RooRealVar* rrv_p1      = new RooRealVar(("rrv_p1_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p1_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.01,-1,1);
            RooRealVar* rrv_p2      = new RooRealVar(("rrv_p2_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p2_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);
            RooRealVar* rrv_p3      = new RooRealVar(("rrv_p3_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p3_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,-1,1);
            RooRealVar* rrv_p4      = new RooRealVar(("rrv_p4_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(),("rrv_p4_AtanBernstein_v5"+label+"_"+channel+spectrum).c_str(), 0.001,0,1);

            RooRealVar* rrv_offset  = new RooRealVar(("rrv_offset_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_offset_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(), 450,350,600);
            RooRealVar* rrv_width   = new RooRealVar(("rrv_width_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),("rrv_width_AtanBernstein_v3"+label+"_"+channel+spectrum).c_str(),50,20,100);

            TString formula;
            formula.Form("(1.+TMath::ATan((%s-%s)/%s))/2.",rrv_x->GetName(),rrv_offset->GetName(), rrv_width->GetName());
	    RooGenericPdf* atan       = new RooGenericPdf(("atan"+label+"_"+channel+spectrum).c_str(),("atan"+label+"_"+channel+spectrum).c_str(),formula.Data(), RooArgList(*rrv_x,*rrv_offset,*rrv_width));
 
            RooBernstein* pol = new RooBernstein(("Bernstein_v5"+label+"_"+channel+spectrum).c_str(),("Bernstein_v5"+label+"_"+channel+spectrum).c_str(),*rrv_x,RooArgList(*rrv_p0,*rrv_p1,*rrv_p2,*rrv_p3,*rrv_p4));

            RooProdPdf* model_pdf = new RooProdPdf(("model_pdf"+label+"_"+channel+spectrum).c_str(),("model_pdf"+label+"_"+channel+spectrum).c_str(),*atan,*pol);
            
            return model_pdf ;
	}

   return NULL ;

  }
      
	

