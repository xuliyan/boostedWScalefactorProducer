#include "FitUtils.h"

void fit_mj_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & model,const std::string & channel,const std::string & wtagger_label,const int & additionalinformation){

  std::cout<<"############### Fit mj single MC sample"<<fileName<<" "<<label<<"  "<<model<<" ##################"<<std::endl;
  //import variable and dataset
  RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");
  RooDataSet* rdataset_mj     = NULL;
  if( TString(workspace->GetName()).Contains("4bias"))
    rdataset_mj = (RooDataSet*) workspace->data(("rdataset4bias"+label+"_"+channel+"_mj").c_str());
  else if ( TString(workspace->GetName()).Contains("4fit"))
    rdataset_mj = (RooDataSet*) workspace->data(("rdataset4fit"+label+"_"+channel+"_mj").c_str());
  rdataset_mj->Print();

  //## make the extended model
  std::vector<std::string>* constraint_list = new std::vector<std::string>(); 
  RooExtendPdf* model_pdf = NULL ; 
  if(additionalinformation == 1)
    model_pdf = MakeExtendedModel(workspace,label+model,model,"_mj",channel,wtagger_label,constraint_list);
  else
    model_pdf = MakeExtendedModel(workspace,label,model,"_mj",channel,wtagger_label,constraint_list);

  RooFitResult* rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::Extended(kTRUE),RooFit::SumW2Error(kTRUE),RooFit::Minimizer("Minuit2"));
  rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult->Print();

  std::cout<<"####################################################"<<std::endl;
  std::cout<<"######## Normalization Factor in mJ ################"<<std::endl;
  std::cout<<"####################################################"<<std::endl;

  //##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV
  RooArgSet* parameters_list = model_pdf->getParameters(*rdataset_mj);
  TIter par = parameters_list->createIterator();
  par.Reset();
  
  workspace->import(*model_pdf);
  workspace->import(*rfresult);


  //## Plot the result
  RooPlot* mplot = rrv_mass_j->frame(RooFit::Title((label+" fitted by "+model).c_str()), RooFit::Bins(int(rrv_mass_j->getBins())));
  rdataset_mj->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible());

  //## draw the error band for an extend pdf
  draw_error_band_extendPdf(rdataset_mj,model_pdf,rfresult,mplot,2,"L");

  //## draw the function
  model_pdf->plotOn(mplot,RooFit::Name("model_mc")); // remove RooFit.VLines() in order to get right pull in the 1st bin

  //## re-draw the dataset
  rdataset_mj->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),RooFit::Name("data"));

  //## Get the pull
  RooPlot* mplot_pull = get_ratio(rrv_mass_j,rdataset_mj,model_pdf,rfresult,0,1);
  mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.2);

  //## CALCULATE CHI2
  RooDataHist* datahist = rdataset_mj->binnedClone((std::string(rdataset_mj->GetName())+"_binnedClone").c_str(),(std::string(rdataset_mj->GetName())+"_binnedClone").c_str());
  int Nbin = int(rrv_mass_j->getBins()); 
  RooArgList rresult_param = rfresult->floatParsFinal();        
  int nparameters =  rresult_param.getSize();                                         
  RooAbsReal* ChiSquare = model_pdf->createChi2(*datahist,RooFit::Extended(kTRUE),RooFit::DataError(RooAbsData::SumW2));
  float chi_over_ndf= ChiSquare->getVal()/(Nbin - nparameters);

  //## Add Chisquare to mplot_pull
  TString Name ; Name.Form("#chi^{2}/ndf = %0.2f ",float(chi_over_ndf));
  TLatex* cs = new TLatex(0.75,0.8,Name.Data());
  cs->SetNDC();
  cs->SetTextSize(0.12);
  cs->AppendPad("same");
  mplot_pull->addObject(cs);

  TString command; 
  command.Form("mkdir -p plots_%s_%s_g1/mj_fitting/",channel.c_str(),wtagger_label.c_str());
  system(command.Data());

  command.Form("plots_%s_%s_g1/mj_fitting/",channel.c_str(),wtagger_label.c_str());
  draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),label+fileName,model,channel,0,0,GetLumi());
     
  RooPlot* mplot_sys = rrv_mass_j->frame(RooFit::Title("plot sys mj"),RooFit::Bins(int(rrv_mass_j->getBins())));

  if(label == "_STop" or label == "_VV" or label == "_WJets0" or label == "_WW_EWK" or label == "_TTbar"){
       
    rdataset_mj->plotOn(mplot_sys,RooFit::Name("MC Events"),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible());
    model_pdf->plotOn(mplot_sys,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
    draw_error_band_extendPdf(rdataset_mj,model_pdf,rfresult,mplot_sys,kBlack,"F");
    rdataset_mj->plotOn(mplot_sys,RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
    
    if(label == "_WJets0" and workspace->pdf(("model_WJets01_"+model+"_"+channel+"_mj").c_str()))
      workspace->pdf(("model_WJets01"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("alt shape"), RooFit::LineColor(kOrange+1),RooFit::Normalization(workspace->var(("rrv_number_WJets01"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));
                  
               
    TLegend* leg = legend4Plot(mplot_sys,0,0.,0.06,0.16,0.,1,channel);
    mplot_sys->addObject(leg);
    rdataset_mj->plotOn(mplot_sys,RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
    model_pdf->plotOn(mplot_sys,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
    mplot_sys->GetYaxis()->SetRangeUser(0,mplot_sys->GetMaximum()*1.2);
        
    command.Form("mkdir -p plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
    system(command.Data());
       
    command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
    draw_canvas_with_pull(mplot_sys,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),"m_j_extended"+label+fileName,model,channel,0,0,GetLumi());
  }

  
  workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setVal(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());
  workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());

    
}


//## Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted
void fit_mlvj_model_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & in_range, const std::string & mlvj_model, const std::string & channel, const std::string & wtagger_label, const int & logy, const int & ismc, const std::string & label_origin, const std::string & FTest_output_txt){

  std::cout<<"############### Fit mlvj single MC sample "<<fileName<<" "<<label<<" "<<mlvj_model<<" "<<in_range<<""<<std::endl;
  //## imporparam_generatedt variable and dataset         
  RooRealVar* rrv_mass_lvj = workspace->var("rrv_mass_lvj");

  int original_binning = rrv_mass_lvj->getBins();

  RooDataSet* rdataset     = NULL;
  if( TString(workspace->GetName()).Contains("4bias"))
    rdataset = (RooDataSet*) workspace->data(("rdataset4bias"+label+in_range+"_"+channel+"_mlvj").c_str());
  else if ( TString(workspace->GetName()).Contains("4fit"))
    rdataset = (RooDataSet*) workspace->data(("rdataset4fit"+label+in_range+"_"+channel+"_mlvj").c_str());

  std::vector<std::string>* constraintlist = new std::vector<std::string>();

  //## make the extended pdf model
  RooExtendPdf* model = MakeExtendedModel(workspace,label+in_range+mlvj_model,mlvj_model,"_mlvj",channel,wtagger_label,constraintlist,ismc);
  //make the fit
  RooFitResult* rfresult = model->fitTo(*rdataset,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE));
  rfresult = model->fitTo(*rdataset,RooFit::Save(1), RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2") );
  rfresult->Print();
  rrv_mass_lvj->setBins(original_binning);
  //## set the name of the result of the fit and put it in the workspace
  rfresult->SetName(("rfresult"+label+in_range+"_"+channel+"_mlvj").c_str());

  workspace->import(*model);
  workspace->import(*rfresult);
      
  //## plot the result
  RooPlot* mplot = rrv_mass_lvj->frame(RooFit::Title(("M_{lvj"+in_range+"} fitted by "+mlvj_model).c_str()), RooFit::Bins(int(rrv_mass_lvj->getBins())));
  rdataset->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible());
  //## plot the error band but don't store the canvas (only plotted without -b option
  if(mlvj_model != "CBBW_v1") draw_error_band_extendPdf(rdataset,model,rfresult,mplot,2,"L");
  model->plotOn(mplot,RooFit::Name("model_mc"));//#, RooFit.VLines()); in order to have the right pull
  rdataset->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Name("data"));

  //## get the pull
  RooPlot* mplot_pull = get_ratio(rrv_mass_lvj,rdataset,model,rfresult,0,1);
  RooArgSet* parameters_list = model->getParameters(*rdataset);
        
  //##CALCULATE CHI2                                                                                                                                                    
  RooDataHist* datahist   = rdataset->binnedClone((std::string(rdataset->GetName())+"_binnedClone").c_str(),(std::string(rdataset->GetName())+"_binnedClone").c_str());
  TH1F* histo_data = (TH1F*) datahist->createHistogram("histo_data",*rrv_mass_lvj);
  histo_data->SetName("histo_data");
  TH1F* histo_func = (TH1F*) model->createHistogram("histo_func",*rrv_mass_lvj);
  histo_func->SetName("histo_func");
        
  int Nbin     = int(rrv_mass_lvj->getBins());
  RooArgSet rresult_param = rfresult->floatParsFinal();
  int nparameters   = rresult_param.getSize();
  RooAbsReal* ChiSquare = model->createChi2(*datahist,RooFit::Extended(kTRUE),RooFit::DataError(RooAbsData::Poisson));
  float chi_over_ndf  = ChiSquare->getVal()/(Nbin-nparameters);
        
  RooHist* residHist = mplot->residHist("data","model_mc");
  float residual = 0. ;
  for(int  iPoint = 0 ; iPoint < residHist->GetN() ; iPoint++){
    double x = 0. , y = 0. ;
    residHist->GetPoint(iPoint,x,y); 
    residual = residual + y*y ;
  }
        
  mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.2);
     
  if( FTest_output_txt != ""){
    std::ofstream file_out_FTest ;  
    file_out_FTest.open(FTest_output_txt.c_str(),fstream::app);
    std::cout<<" model pdf ----> dump file ftest "<<std::endl;
    file_out_FTest << " ###################### \n";
    file_out_FTest << " Model pdf "+std::string(model->GetName())+" \n";

    TString Line ; Line.Form(" Chi2 Chi2var %f \n",chi_over_ndf);
    file_out_FTest << std::string(Line.Data()); 
    Line.Form(" Residual %f  Nbin %d  nparameters %d \n",residual,Nbin,nparameters);
    file_out_FTest << std::string(Line.Data());
    file_out_FTest.close();
  }

  //##Add Chisquare to mplot_pull   
  TString command; command.Form("#chi^{2}/ndf = %0.2f ",float(chi_over_ndf));                                                                                
  TLatex* cs2 = new TLatex(0.25,0.8,command.Data());
  cs2->SetNDC();
  cs2->SetTextSize(0.12);
  cs2->AppendPad("same");
  mplot_pull->addObject(cs2);

  command.Form("mkdir -p plots_%s_%s_g1/mlvj_fitting/",channel.c_str(),wtagger_label.c_str());
  system(command.Data());

  command.Form("plots_%s_%s_g1/mlvj_fitting/",channel.c_str(),wtagger_label.c_str()); 
  draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),label+fileName,mlvj_model+in_range,"em",0,logy,GetLumi());

  RooPlot* mplot_sys = rrv_mass_lvj->frame( RooFit::Bins(int(rrv_mass_lvj->getBins())));

  if(not TString(label).Contains("_jes") and not TString(label).Contains("_jer")){

    rdataset->plotOn(mplot_sys, RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );
    model->plotOn(mplot_sys,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
    if(mlvj_model != "CBBW_v1") draw_error_band_extendPdf(rdataset, model, rfresult,mplot_sys,1,"F");

    if(workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_up"), RooFit::LineColor(kRed),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

    if(workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_dn"), RooFit::LineColor(kBlue),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

    if(workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer"), RooFit::LineColor(kAzure+1),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

    if(workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_up"), RooFit::LineColor(kGreen+1),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

    if(workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_dn"), RooFit::LineColor(6),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

    if(label == "_WJets0" and workspace->pdf(("model_WJets01"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
      workspace->pdf(("model_WJets01"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("alt shape"), RooFit::LineColor(kOrange+1),RooFit::Normalization(workspace->var(("rrv_number_WJets01"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()/(rdataset->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal())));

  }

  TLegend* legend = legend4Plot(mplot_sys,0,-0.07,0.02,0.01,0.,1,channel);
  mplot_sys->addObject(legend);
  rdataset->plotOn(mplot_sys,RooFit::Name("MC Events"),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));
  model->plotOn(mplot_sys,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
  mplot_sys->GetYaxis()->SetRangeUser(0,mplot_sys->GetMaximum()*1.2);

  command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
  draw_canvas_with_pull(mplot_sys,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),"m_lvj_extended"+label+in_range,mlvj_model,channel,0,logy,GetLumi());

      
  std::cout<<"rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj"<<std::endl;
  std::cout<<"rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj"<<std::endl;
  workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setVal(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal());
  workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setError(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal() );
  workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->Print();
      

}

void get_mj_normalization_insignalregion(RooWorkspace* workspace, const std::string & label, const std::string & model_name, const std::string & channel){

  std::cout<<"################## get mj normalization "<<label<<" ################## "<<std::endl;
  RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");
  RooAbsPdf*  model      = workspace->pdf(("model"+label+model_name+"_"+channel+"_mj").c_str());
  std::cout<<" model "<<"model"+label+"_"+channel+"_mj"<<std::endl;
                                
  RooAbsReal* fullInt   = model->createIntegral(*rrv_mass_j,*rrv_mass_j);
  RooAbsReal* sb_loInt  = model->createIntegral(*rrv_mass_j,*rrv_mass_j,("sb_lo"));
  RooAbsReal* signalInt = model->createIntegral(*rrv_mass_j,*rrv_mass_j,("signal_region"));
  RooAbsReal* sb_hiInt  = model->createIntegral(*rrv_mass_j,*rrv_mass_j,("sb_hi"));

  double fullInt_val   = fullInt->getVal();
  double sb_loInt_val  = sb_loInt->getVal()/fullInt_val;
  double sb_hiInt_val  = sb_hiInt->getVal()/fullInt_val;
  double signalInt_val = signalInt->getVal()/fullInt_val;

  std::cout<<"########### Events Number in MC Dataset: #############"<<std::endl;
  workspace->var(("rrv_number_dataset_sb_lo"+label+"_"+channel+"_mj").c_str())->Print();
  workspace->var(("rrv_number_dataset_signal_region"+label+"_"+channel+"_mj").c_str())->Print();
  workspace->var(("rrv_number_dataset_sb_hi"+label+"_"+channel+"_mj").c_str())->Print();
            
  std::cout<<"########### Events Number get from fit: ##############"<<std::endl;
  RooRealVar* rrv_tmp = workspace->var(("rrv_number"+label+model_name+"_"+channel+"_mj").c_str());
  rrv_tmp->Print();
  
  std::cout<< "Events Number in sideband_low : "<<rrv_tmp->getVal()*sb_loInt_val<<std::endl;
  std::cout<< "Events Number in Signal Region: "<<rrv_tmp->getVal()*signalInt_val<<std::endl;
  std::cout<< "Events Number in sideband_high: "<<rrv_tmp->getVal()*sb_hiInt_val<<std::endl;
  std::cout<< "Total Number in sidebands     : "<<rrv_tmp->getVal()*(sb_loInt_val+sb_hiInt_val)<<std::endl;
  std::cout<< "Ratio signal_region/sidebands : "<<(signalInt_val/(sb_loInt_val+sb_hiInt_val))<<std::endl;

  return ;
}

//##### Function that calculate the normalization inside the mlvj signal region (mass window around the resonance in order to fill datacards)
void get_mlvj_normalization_insignalregion(RooWorkspace* workspace,const std::string & label,  const std::string & model_name, const std::string & mlvj_region,const std::string & channel, const int & isAlpha){

  std::cout<< "############### get mlvj normalization inside SR "<<label<<" "<<model_name<<" ##################"<<std::endl;

  RooRealVar* rrv_mass_lvj = workspace->var("rrv_mass_lvj");
  RooAbsPdf* model ;
  
  if(isAlpha == 0) 
    model = workspace->pdf(("model"+label+mlvj_region+model_name+"_"+channel+"_mlvj").c_str());  
  else
    model = workspace->pdf(("model_pdf"+label+mlvj_region+model_name+"_"+channel+"_after_correct_mlvj").c_str());
  
  std::cout<<"model_pdf"+label+mlvj_region+model_name+"_"+channel+"_after_correct_mlvj"<<std::endl;

  RooAbsReal* fullInt   = model->createIntegral(*rrv_mass_lvj,*rrv_mass_lvj);
  RooAbsReal* signalInt = model->createIntegral(*rrv_mass_lvj,*rrv_mass_lvj,("signal_region"));

  double fullInt_val = fullInt->getVal();
  double signalInt_val = signalInt->getVal()/fullInt_val;

  //## integal in the signal region
  std::cout<< "######### integral in SR: "<<label<<" signalInt ="<<signalInt_val<<std::endl;
  std::cout<< "####### Events Number in MC Dataset: "<<std::endl;
  workspace->var(("rrv_number_dataset_signal_region"+label+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_dataset_AllRange"+label+"_"+channel+"_mlvj").c_str())->Print();

  std::cout<< "########## Events Number get from fit : "<<std::endl;
  RooRealVar* rrv_tmp = workspace->var(("rrv_number"+label+"_signal_region"+model_name+"_"+channel+"_mlvj").c_str());
  std::cout<< "Events Number in Signal Region from fitting: "<<rrv_tmp->getVal()*signalInt_val<<std::endl;

  //#### store the info in the output file
  if( not workspace->var(("rrv_number_fitting_signal_region"+label+"_"+channel+"_mlvj").c_str())){
    RooRealVar* rrv_number_fitting_signal_region_mlvj = new RooRealVar(("rrv_number_fitting_signal_region"+label+"_"+channel+"_mlvj").c_str(),("rrv_number_fitting_signal_region"+label+"_"+channel+"_mlvj").c_str(), rrv_tmp->getVal()*signalInt_val);
    workspace->import(*rrv_number_fitting_signal_region_mlvj);
  }
  else{
    workspace->var(("rrv_number_fitting_signal_region"+label+"_"+channel+"_mlvj").c_str())->setVal(rrv_tmp->getVal()*signalInt_val);
    workspace->var(("rrv_number_fitting_signal_region"+label+"_"+channel+"_mlvj").c_str())->Print();
  }

  return;

}


// method to do fail and pass fit for the scale factor evaluation                                                                                                                      
void ScaleFactorTTbarControlSampleFit(RooWorkspace* workspace, std::map<std::string,std::string > mj_shape, std::map<std::string,int> color_palet, std::vector<std::string>* constraintlist_data, std::vector<std::string>* constraintlist_MC, const std::string & label, const std::string & channel, const std::string & wtagger, const double & ca8_ungroomed_pt_min, const double & ca8_ungroomed_pt_max){

  std::cout<<"############################## Pass: dataset #################################### "<<std::endl;

  RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");
  RooDataSet* rdataset_data_mj  = (RooDataSet*) workspace->data(("rdataset_data"+ label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_TTbar_mj = (RooDataSet*) workspace->data(("rdataset_TTbar"+ label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_STop_mj  = (RooDataSet*) workspace->data(("rdataset_STop"+  label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_VV_mj    = (RooDataSet*) workspace->data(("rdataset_VV"+    label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_WJets_mj = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+channel+"_mj").c_str());

  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj);
  
  std::cout<<"############################### Pass: Histograms for plotting reasons  ################################### "<<std::endl;
  RooAbsPdf* model_histpdf_TTbar = workspace->pdf((std::string(rdataset_TTbar_mj->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_STop  = workspace->pdf((std::string(rdataset_STop_mj->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_VV    = workspace->pdf((std::string(rdataset_VV_mj->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_WJets = workspace->pdf((std::string(rdataset_WJets_mj->GetName())+"_histpdf").c_str());


  RooRealVar* number_TTbar = new RooRealVar(("rrv_number_TTbar"+label+"_"+channel).c_str(),("rrv_number_TTbar"+label+"_"+channel).c_str(),rdataset_TTbar_mj->sumEntries());
  RooRealVar* number_STop  = new RooRealVar(("rrv_number_STop"+label+"_"+channel).c_str(),("rrv_number_STop"+label+"_"+channel).c_str(),rdataset_STop_mj->sumEntries());
  RooRealVar* number_VV    = new RooRealVar(("rrv_number_VV"+label+"_"+channel).c_str(),("rrv_number_VV"+label+"_"+channel).c_str(),rdataset_VV_mj->sumEntries());
  RooRealVar* number_WJets = new RooRealVar(("rrv_number_WJets"+label+"_"+channel).c_str(),("rrv_number_WJets"+label+"_"+channel).c_str(),rdataset_WJets_mj->sumEntries());

  RooAddPdf* model_TTbar_STop_VV_WJets = new RooAddPdf(("model_TTbar_STop_VV_WJets"+label+"_"+channel).c_str(),("model_TTbar_STop_VV_WJets"+label+"_"+channel).c_str(),RooArgList(*model_histpdf_TTbar,*model_histpdf_STop,*model_histpdf_VV,*model_histpdf_WJets),RooArgList(*number_TTbar,*number_STop,*number_VV,*number_WJets));
  workspace->import(*model_TTbar_STop_VV_WJets);
  
  std::cout<<"############################## Fail: dataset #################################### "<<std::endl;
  RooDataSet* rdataset_data_mj_fail  = (RooDataSet*) workspace->data(("rdataset_data"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_TTbar_mj_fail = (RooDataSet*) workspace->data(("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_STop_mj_fail  = (RooDataSet*) workspace->data(("rdataset_STop"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_VV_mj_fail    = (RooDataSet*) workspace->data(("rdataset_VV"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_WJets_mj_fail = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());

  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj_fail);

  std::cout<<"############################### Fail: model ################################### "<<std::endl;
  RooAbsPdf* model_histpdf_TTbar_fail = workspace->pdf((std::string(rdataset_TTbar_mj_fail->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_STop_fail  = workspace->pdf((std::string(rdataset_STop_mj_fail->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_VV_fail    = workspace->pdf((std::string(rdataset_VV_mj_fail->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_WJets_fail = workspace->pdf((std::string(rdataset_WJets_mj_fail->GetName())+"_histpdf").c_str());

  RooRealVar* number_TTbar_fail  = new RooRealVar(("rrv_number_TTbar_fail"+label+"_"+channel).c_str(),("rrv_number_TTbar_fail"+label+"_"+channel).c_str(),rdataset_TTbar_mj_fail->sumEntries());
  RooRealVar* number_STop_fail   = new RooRealVar(("rrv_number_STop_fail"+label +"_"+channel).c_str(),("rrv_number_STop_fail"+label+"_"+channel).c_str(),rdataset_STop_mj_fail->sumEntries());
  RooRealVar* number_VV_fail     = new RooRealVar(("rrv_number_VV_fail"+label   +"_"+channel).c_str(),("rrv_number_VV_fail"+label+"_"+channel).c_str(),rdataset_VV_mj_fail->sumEntries());
  RooRealVar* number_WJets_fail  = new RooRealVar(("rrv_number_WJets_fail"+label+"_"+channel).c_str(),("rrv_number_WJets_fail"+label+"_"+channel).c_str(),rdataset_WJets_mj_fail->sumEntries());

  RooAddPdf* model_TTbar_STop_VV_WJets_fail = new RooAddPdf(("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel).c_str(),("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel).c_str(),RooArgList(*model_histpdf_TTbar_fail,*model_histpdf_STop_fail,*model_histpdf_VV_fail,*model_histpdf_WJets_fail), RooArgList(*number_TTbar_fail,*number_STop_fail,*number_VV_fail,*number_WJets_fail));
  workspace->import(*model_TTbar_STop_VV_WJets_fail);
  
  /// inclusive scale factor from total integral                                                                                                                                   
  double scale_number_TTbar_STop_VV_WJets      = (rdataset_TTbar_mj->sumEntries()+rdataset_STop_mj->sumEntries()+rdataset_VV_mj->sumEntries()+rdataset_WJets_mj->sumEntries())/(rdataset_data_mj->sumEntries()+rdataset_data_mj_fail->sumEntries());
  double scale_number_TTbar_STop_VV_WJets_fail = (rdataset_TTbar_mj_fail->sumEntries()+rdataset_STop_mj_fail->sumEntries()+rdataset_VV_mj_fail->sumEntries()+rdataset_WJets_mj_fail->sumEntries())/( rdataset_data_mj->sumEntries()+rdataset_data_mj_fail->sumEntries());

  RooRealVar* rrv_scale_number_TTbar_STop_VV_WJets = new RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets"+label).c_str(),("rrv_scale_number_TTbar_STop_VV_WJets"+label).c_str(),scale_number_TTbar_STop_VV_WJets);
  RooRealVar* rrv_scale_number_TTbar_STop_VV_WJets_fail = new RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label).c_str(),("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label).c_str(),scale_number_TTbar_STop_VV_WJets_fail);
  workspace->import(*rrv_scale_number_TTbar_STop_VV_WJets);
  workspace->import(*rrv_scale_number_TTbar_STop_VV_WJets_fail);


  /// All the shape parameters and normalization are fixed                                                                                                                          
  std::cout<<"############################### Pass: Single MC model ################################### "<<std::endl;
  RooAbsPdf* model_STop  = get_STop_mj_Model(workspace,"_STop"+label,mj_shape["STop"],channel);
  RooAbsPdf* model_VV    = get_VV_mj_Model(workspace,"_VV"+label,mj_shape["VV"],channel);
  RooAbsPdf* model_WJets = get_General_mj_Model(workspace,"_WJets0"+label,mj_shape["WJets0"],channel);

  std::cout<<"############################### Fail: Single MC model ################################### "<<std::endl;
  RooAbsPdf* model_STop_fail  = get_STop_mj_Model(workspace,"_STop"+label+"_failtau2tau1cut",mj_shape["STop_fail"],channel);
  RooAbsPdf* model_VV_fail    = get_VV_mj_Model(workspace,"_VV"+label+"_failtau2tau1cut",mj_shape["VV_fail"],channel);
  RooAbsPdf* model_WJets_fail = get_General_mj_Model(workspace,"_WJets0"+label+"_failtau2tau1cut",mj_shape["WJets0_fail"],channel);

  /// Model for unmatched events passing and failing the cut --> ErfExp                                                                                                              
  std::cout<<"############################### Pass: TTbar in Data Bkg ################################### "<<std::endl;
  RooAbsPdf* model_bkg_data = MakeExtendedModel(workspace,"_bkg_data"+label,mj_shape["bkg_data"],"_mj",channel,wtagger,constraintlist_data);
  std::cout<<"############################### Fail: TTbar in Data Bkg ################################### "<<std::endl;
  RooAbsPdf* model_bkg_data_fail = MakeExtendedModel(workspace,"_bkg_data"+label+"_failtau2tau1cut",mj_shape["bkg_data_fail"],"_mj",channel,wtagger,constraintlist_data);
  // Model for matched events passing and failing the cut --> 2Gaus_ttbar and GausChebychev_ttbar_failtau2tau1cut                                                                   
  std::cout<<"############################### Pass: Data ttbar resonant component ################################### "<<std::endl;
  RooAbsPdf* model_ttbar_data = MakeModelTTbarControlSample(workspace,"_ttbar_data"+label,mj_shape["signal_data"],"_mj",channel,wtagger,label,constraintlist_data);    
  std::cout<<"############################### Fail: Data ttbar resonant component ################################### "<<std::endl;
  RooAbsPdf* model_ttbar_data_fail = MakeModelTTbarControlSample(workspace,"_ttbar_data"+label+"_failtau2tau1cut",mj_shape["signal_data_fail"],"_mj",channel,wtagger,label,constraintlist_data);


  std::cout<<"############################### Total Data Pdf Fail ################################### "<<std::endl;
  RooAbsPdf* model_data_fail = new RooAddPdf(("model_data"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str(),("model_data+"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str(),RooArgList(*model_ttbar_data_fail,*model_bkg_data_fail,*model_STop_fail,*model_VV_fail,*model_WJets_fail));
  std::cout<<"############################### Total Data Pdf Pass ################################### "<<std::endl;
  RooAbsPdf* model_data      = new RooAddPdf(("model_data"+label+"_"+channel).c_str(),("model_data"+label+"_"+channel).c_str(), RooArgList(*model_ttbar_data,*model_bkg_data,*model_STop,*model_VV,*model_WJets));

  workspace->import(*model_data);
  workspace->import(*model_data_fail);
  
  std::cout<<"############################### Data Event Category ################################### "<<std::endl;
  RooCategory* category_p_f = workspace->cat(("category_p_f"+label+"_"+channel).c_str());

  ///### Define the simultaneous fit                                                                                                                                                    
  RooSimultaneous* simPdf_data = new RooSimultaneous(("simPdf_data"+label+"_"+channel).c_str(),("simPdf_data"+label+"_"+channel).c_str(),*category_p_f);
  simPdf_data->addPdf(*model_data,"pass");
  simPdf_data->addPdf(*model_data_fail,"fail");
  simPdf_data->Print();

  RooDataSet* combData_p_f_data = (RooDataSet*) workspace->data(("combData_p_f_data"+label+"_"+channel).c_str());
  combData_p_f_data->Print();

  std::cout<<" ############################### Data Total Fit ################################### "<<std::endl;
  RooArgSet* pdfconstrainslist_data = new RooArgSet(("pdfconstrainslist_data"+label+"_"+channel).c_str());
  if(constraintlist_data!=NULL){
    for(unsigned int i = 0; i < constraintlist_data->size(); i++){
      workspace->pdf(constraintlist_data->at(i).c_str())->Print();
      pdfconstrainslist_data->add(*workspace->pdf(constraintlist_data->at(i).c_str()));
    }
  }
  std::cout<<"Fit ################ "<<std::endl;       
  RooFitResult* rfresult_data = NULL;
  rfresult_data = simPdf_data->fitTo(*combData_p_f_data,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data),RooFit::Minimizer("Minuit2"));
  rfresult_data = simPdf_data->fitTo(*combData_p_f_data,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data),RooFit::Minimizer("Minuit2"));
  rfresult_data->Print();
  
  /// fit TotalMC                                                                                                                                                                       

  std::cout<< " ############################### Pass: MC TTbar Bkg ################################### "<<std::endl;
  RooAbsPdf* model_bkg_TotalMC = MakeExtendedModel(workspace,"_bkg_TotalMC"+label,mj_shape["bkg_mc"],"_mj",channel,wtagger,constraintlist_MC);
  std::cout<< " ############################### Fail: MC TTbar Bkg ################################### "<<std::endl;
  RooAbsPdf* model_bkg_TotalMC_fail = MakeExtendedModel(workspace,"_bkg_TotalMC"+label+"_failtau2tau1cut",mj_shape["bkg_mc_fail"],"_mj",channel,wtagger,constraintlist_MC);
  std::cout<< " ############################### Pass: MC ttbar Resonant ################################### "<<std::endl;
  RooAbsPdf* model_ttbar_TotalMC  = MakeModelTTbarControlSample(workspace,"_ttbar_TotalMC"+label,mj_shape["signal_mc"],"_mj",channel,wtagger,label,constraintlist_MC ); 
  std::cout<< " ############################### Fail: MC ttbar Resonant ################################### "<<std::endl;
  RooAbsPdf* model_ttbar_TotalMC_fail  = MakeModelTTbarControlSample(workspace,"_ttbar_TotalMC"+label+"_failtau2tau1cut",mj_shape["signal_mc_fail"],"_mj",channel,wtagger,label,constraintlist_MC);                                                                                                                      
  std::cout<< " ############################### Fail: MC total PDF ################################### "<<std::endl;
  RooAddPdf* model_TotalMC_fail = new RooAddPdf(("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str(),("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str(),RooArgList(*model_ttbar_TotalMC_fail,*model_bkg_TotalMC_fail,*model_STop_fail,*model_VV_fail,*model_WJets_fail));
  std::cout<< " ############################### Pass: MC total PDF ################################### "<<std::endl;
  RooAddPdf* model_TotalMC = new RooAddPdf(("model_TotalMC"+label+"_"+channel).c_str(),("model_TotalMC"+label+"_"+channel).c_str(),RooArgList(*model_ttbar_TotalMC,*model_bkg_TotalMC,*model_STop,*model_VV,*model_WJets));
 
  workspace->import(*model_TotalMC_fail);
  workspace->import(*model_TotalMC);


  std::cout<< " ############################### Data MC Category ################################### "<<std::endl;
  RooSimultaneous* simPdf_TotalMC = new RooSimultaneous(("simPdf_TotalMC"+label+"_"+channel).c_str(),("simPdf_TotalMC"+label+"_"+channel).c_str(),*category_p_f);
  simPdf_TotalMC->addPdf(*model_TotalMC,"pass");
  simPdf_TotalMC->addPdf(*model_TotalMC_fail,"fail");

  RooDataSet* combData_p_f_TotalMC = (RooDataSet*) workspace->data(("combData_p_f_TotalMC"+label+"_"+channel).c_str());

  std::cout<< " ############################### MC Total Fit ################################### "<<std::endl;
  RooFitResult* rfresult_TotalMC = NULL ;
  RooArgSet* pdfconstrainslist_TotalMC = new RooArgSet(("pdfconstrainslist_TotalMC"+label+"_"+channel).c_str());
  if(constraintlist_MC!=NULL){
    for(unsigned int i = 0; i < constraintlist_MC->size() ; i++){
      workspace->pdf(constraintlist_MC->at(i).c_str())->Print();
      pdfconstrainslist_TotalMC->add(*workspace->pdf(constraintlist_MC->at(i).c_str()));
    }
  }

  rfresult_TotalMC = simPdf_TotalMC->fitTo(*combData_p_f_TotalMC,RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC), RooFit::Minimizer("Minuit2"), RooFit::SumW2Error(kTRUE));
  rfresult_TotalMC = simPdf_TotalMC->fitTo(*combData_p_f_TotalMC,RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC), RooFit::Minimizer("Minuit2"), RooFit::SumW2Error(kTRUE));
  rfresult_TotalMC->Print();
  
         
  RooPlot* xframe_data             = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));
  RooPlot* xframe_data_fail        = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));
  // RooPlot* xframe_data_extremefail = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));

  TString cut ;
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::Poisson), RooFit::XErrorSize(0) );

  //#############################################################  
  ///plot total pass mc normalizing it to data, passing + failing                                                                                                                   
  //#############################################################  

  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
  
  cut.Form("%s,%s,%s",model_histpdf_STop->GetName(),model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("STop"),RooFit::Components(cut.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("VV"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s",model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("WJets"),RooFit::Components(cut.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

  //pass plots
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());  
  cut.Form("%s,%s,%s",model_histpdf_STop->GetName(), model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());       
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("STop_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV->GetName(),model_histpdf_WJets->GetName());      
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("VV_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s",model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("WJets_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

  // plot data again                                                                                                                                                                 
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data"), RooFit::Cut(cut.Data()), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::Poisson), RooFit::XErrorSize(0));

  // plot mc fit function                                                                                                                                                            
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_TotalMC->plotOn(xframe_data,RooFit::Name("TotalMC_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible() , RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );

  simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("MC fit"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed)), RooFit::DataError(RooAbsData::SumW2);
  //  cut.Form("%s,%s,%s,%s",model_bkg_TotalMC->GetName(),model_STop->GetName(),model_VV->GetName(),model_WJets->GetName());
  // simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  // plot data fit function
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());                                                    
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  simPdf_data->plotOn(xframe_data,RooFit::Name("total data fit"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  cut.Form("%s,%s,%s,%s",model_bkg_data->GetName(),model_STop->GetName(),model_VV->GetName(),model_WJets->GetName());
  simPdf_data->plotOn(xframe_data,RooFit::Name("data fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed));

  ///plot total fail mc normalizing it to data, passing + failing                                                                                                                   
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data_fail,RooFit::Name("data_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  if(TString(label).Contains("herwig"))
    model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar_herwig"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
  else
    model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());

  cut.Form("%s,%s,%s",model_histpdf_STop_fail->GetName(),model_histpdf_VV_fail->GetName(),model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("STop"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());  
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("VV"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s", model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("WJets"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

  //solid line
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s,%s",model_histpdf_STop_fail->GetName(), model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());       
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("STop_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("VV_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s", model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("WJets_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

  
  //fail plots -> plot MC fit                                                                                                                                          
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_TotalMC->plotOn(xframe_data_fail,RooFit::Name("TotalMC_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));

  simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("total MC fit"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed));
  //  cut.Form("%s,%s,%s,%s",model_bkg_TotalMC_fail->GetName(),model_STop_fail->GetName(),model_VV_fail->GetName(),model_WJets_fail->GetName());
  //  simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  //fail plots -> plot data fit                                                                                                                                                       
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data_fail,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  simPdf_data->plotOn(xframe_data_fail,RooFit::Name("total data fit"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data_fit_invisible"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  cut.Form("%s,%s,%s,%s",model_bkg_data_fail->GetName(),model_STop_fail->GetName(),model_VV_fail->GetName(),model_WJets_fail->GetName());
  simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()),RooFit::LineColor(kRed));

  //signal window                                                                                                                                                                     
  TLine* lowerLine = new TLine(65.,0.,65,xframe_data->GetMaximum()*0.7); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kBlack); lowerLine->SetLineStyle(9);
  TLine* upperLine = new TLine(105.,0.,105.,xframe_data->GetMaximum()*0.7); upperLine->SetLineWidth(2); upperLine->SetLineColor(kBlack); upperLine->SetLineStyle(9);
  xframe_data->addObject(lowerLine); xframe_data->addObject(upperLine);
  lowerLine = new TLine(65,0.,65,xframe_data_fail->GetMaximum()*0.7); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kBlack); lowerLine->SetLineStyle(9);
  upperLine = new TLine(105.,0.,105.,xframe_data_fail->GetMaximum()*0.7); upperLine->SetLineWidth(2); upperLine->SetLineColor(kBlack); upperLine->SetLineStyle(9);
  xframe_data_fail->addObject(lowerLine); xframe_data_fail->addObject(upperLine);

  //make the legend 
  TLegend* leg_data = legend4Plot(xframe_data,0,-0.2,0.07,0.04,0.,1,channel);                                                                                          
  xframe_data->addObject(leg_data);
  TLegend* leg_data_fail = legend4Plot(xframe_data_fail,0,-0.2,0.07,0.04,0.,1,channel);                                           
  xframe_data_fail->addObject(leg_data_fail);

  //add mean and width                                                                                                                                                                
  RooRealVar* rrv_mean_gaus_data     = workspace->var(("rrv_mean1_gaus_ttbar_data"    +label+"_"+channel).c_str());
  RooRealVar* rrv_sigma_gaus_data    = workspace->var(("rrv_sigma1_gaus_ttbar_data"   +label+"_"+channel).c_str());
  RooRealVar* rrv_mean_gaus_TotalMC  = workspace->var(("rrv_mean1_gaus_ttbar_TotalMC" +label+"_"+channel).c_str());
  RooRealVar* rrv_sigma_gaus_TotalMC = workspace->var(("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+channel).c_str());

  if(rrv_mean_gaus_TotalMC){
    std::cout<<" IN LOOP rrv_mean_gaus_TotalMC!!!!!!!!!"<<std::endl;
    cut.Form("Mean_{MC } = %3.1f #pm %2.1f",rrv_mean_gaus_TotalMC->getVal(),rrv_mean_gaus_TotalMC->getError());
    TLatex* tl_MC_mean  = new TLatex(0.25 ,0.62, cut.Data());
    cut.Form("Sigma_{MC }= %2.1f #pm %2.1f",rrv_sigma_gaus_TotalMC->getVal(),rrv_sigma_gaus_TotalMC->getError());
    TLatex* tl_MC_sigma = new TLatex(0.25 ,0.57, cut.Data());
    tl_MC_mean->SetNDC(); tl_MC_sigma->SetNDC();
    tl_MC_mean->SetTextSize(0.03);
    tl_MC_sigma->SetTextSize(0.03);
    //xframe_data.addObject(tl_MC_mean);                                                                                                                                           
    //xframe_data.addObject(tl_MC_sigma);                                                                                                                                          
  }
  if(rrv_mean_gaus_data){
    cut.Form("Mean_{data} = %3.1f #pm %2.1f",rrv_mean_gaus_data->getVal(), rrv_mean_gaus_data->getError());
    TLatex* tl_data_mean  = new TLatex(0.25 ,0.52, cut.Data());
    cut.Form("Sigma_{data}= %2.1f #pm %2.1f",rrv_sigma_gaus_data->getVal(), rrv_sigma_gaus_data->getError());
    TLatex* tl_data_sigma = new TLatex(0.25 ,0.47, cut.Data());
    tl_data_mean->SetNDC(); tl_data_sigma->SetNDC();
    tl_data_mean->SetTextSize(0.03);
    tl_data_sigma->SetTextSize(0.03);
    //xframe_data.addObject(tl_data_mean);                                                                                                                                         
    //xframe_data.addObject(tl_data_sigma);                                                                                                                                        
  }

  xframe_data->GetYaxis()->SetRangeUser(0.,xframe_data->GetMaximum()*1.3);
  xframe_data_fail->GetYaxis()->SetRangeUser(0.,xframe_data_fail->GetMaximum()*1.3);
  TString nameDir;  nameDir.Form("plots_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/",channel.c_str(),wtagger.c_str(),wtagger.c_str());
  TString namePlot; namePlot.Form("control%s_%s_%s_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
  draw_canvas(xframe_data,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);    
  namePlot.Form("control%s_%s_%s_fail_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
  draw_canvas(xframe_data_fail,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);

  ShowParam_Pdf(simPdf_data,new RooArgSet(*rrv_mass_j,*category_p_f));
  ShowParam_Pdf(simPdf_TotalMC,new RooArgSet(*rrv_mass_j,*category_p_f));

  rfresult_TotalMC->covarianceMatrix().Print();
  rfresult_data->covarianceMatrix().Print();
  rfresult_TotalMC->Print();
  rfresult_data->Print();

  ShowParam_Pdf(simPdf_data,new RooArgSet(*rrv_mass_j,*category_p_f));

  // //Sf for LP
 //  std::cout<<"########################## Extreme Fail Analysis ############################## "<<std::endl;
  // RooDataSet* rdataset_data_mj_extremefail    = (RooDataSet*) workspace->data(("rdataset_data"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
  // RooDataSet* rdataset_TotalMC_mj_extremefail = (RooDataSet*) workspace->data(("rdataset_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
  // RooDataSet* rdataset_TTbar_mj_extremefail   = (RooDataSet*) workspace->data(("rdataset_TTbar"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
  // RooDataSet* rdataset_STop_mj_extremefail    = (RooDataSet*) workspace->data(("rdataset_STop"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
  // RooDataSet* rdataset_VV_mj_extremefail      = (RooDataSet*) workspace->data(("rdataset_VV"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
  // RooDataSet* rdataset_WJets_mj_extremefail   = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
 //
 //  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_extremefail);
 //  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj_extremefail);
 //  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj_extremefail);
 //  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj_extremefail);
 //
 //  RooAbsPdf* model_histpdf_TTbar_extremefail = workspace->pdf((std::string(rdataset_TTbar_mj_extremefail->GetName())+"_histpdf").c_str());
 //  RooAbsPdf* model_histpdf_STop_extremefail  = workspace->pdf((std::string(rdataset_STop_mj_extremefail->GetName())+"_histpdf").c_str());
 //  RooAbsPdf* model_histpdf_VV_extremefail    = workspace->pdf((std::string(rdataset_VV_mj_extremefail->GetName())+"_histpdf").c_str());
 //  RooAbsPdf* model_histpdf_WJets_extremefail = workspace->pdf((std::string(rdataset_WJets_mj_extremefail->GetName())+"_histpdf").c_str());
 //
  // RooRealVar* number_TTbar_extremefail = new RooRealVar(("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_TTbar_mj_extremefail->sumEntries());
  // RooRealVar* number_STop_extremefail  = new RooRealVar(("rrv_number_STop"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_STop"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_STop_mj_extremefail->sumEntries());
  // RooRealVar*  number_VV_extremefail   = new RooRealVar(("rrv_number_VV"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_VV"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_VV_mj_extremefail->sumEntries());
  // RooRealVar* number_WJets_extremefail = new RooRealVar(("rrv_number_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_WJets_mj_extremefail->sumEntries());

 //  RooAddPdf* model_TTbar_STop_VV_WJets_extremefail = new RooAddPdf(("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail_"+channel).c_str(),RooArgList(*model_histpdf_TTbar_extremefail,*model_histpdf_STop_extremefail,*model_histpdf_VV_extremefail,*model_histpdf_WJets_extremefail), RooArgList(*number_TTbar_extremefail,*number_STop_extremefail,*number_VV_extremefail,*number_WJets_extremefail) );
 //  workspace->import(*model_TTbar_STop_VV_WJets_extremefail);
 //
 //  double scale_number_TTbar_STop_VV_WJets_extremefail = (rdataset_TTbar_mj_extremefail->sumEntries()+rdataset_STop_mj_extremefail->sumEntries()+rdataset_VV_mj_extremefail->sumEntries()+rdataset_WJets_mj_extremefail->sumEntries() )/(rdataset_data_mj_extremefail->sumEntries()) ;
 //
 //  RooRealVar* rrv_scale_number_TTbar_STop_VV_WJets_extremefail = new RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail").c_str(),("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail").c_str(),scale_number_TTbar_STop_VV_WJets_extremefail);
 //  workspace->import(*rrv_scale_number_TTbar_STop_VV_WJets_extremefail);
 //
 //  RooAbsPdf* model_STop_extremefail  = get_STop_mj_Model(workspace,"_STop"+label+"_"+"extremefailtau2tau1cut",mj_shape["STop_extremefail"],channel);
 //  RooAbsPdf* model_VV_extremefail    = get_VV_mj_Model(workspace,"_VV"+label+"_"+"extremefailtau2tau1cut",mj_shape["VV_extremefail"],channel);
 //  RooAbsPdf* model_WJets_extremefail = get_General_mj_Model(workspace,"_WJets0"+label+"_"+"extremefailtau2tau1cut",mj_shape["WJets0_extremefail"],channel);
 //
 //
 //  std::vector<std::string>* constrainslist_data_LP = new std::vector<std::string>();
 //  RooAbsPdf* model_ttbar_data_extremefail = MakeExtendedModel(workspace,"_ttbar_data"+label+"_"+"extremefailtau2tau1cut",mj_shape["data_extremefail"],"_mj",channel,wtagger,constrainslist_data_LP);
 //  RooAbsPdf* model_bkg_data_extremefail   = MakeExtendedModel(workspace,"_bkg_data"+label+"_"+"extremefailtau2tau1cut",mj_shape["data_bkg_extremefail"],"_mj",channel,wtagger,constrainslist_data_LP);
 //
 //  //make the model for the extreme fail fit on data
 //  RooAddPdf* model_data_extremefail = new RooAddPdf(("model_data"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(),("model_data"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(), RooArgList(*model_ttbar_data_extremefail,*model_bkg_data_extremefail,*model_STop_extremefail,*model_VV_extremefail,*model_WJets_extremefail));
 //  workspace->import(*model_data_extremefail);
 //
 //  //constraint the parameter in the extreme fail fit
 //  RooArgSet* pdfconstrainslist_data_LP = new RooArgSet(("pdfconstrainslist_data_LP"+label+"_"+channel).c_str());
 //  if(constrainslist_data_LP!=NULL){
 //    for(unsigned int i = 0; i < constrainslist_data_LP->size(); i++){
 //      workspace->pdf(constrainslist_data_LP->at(i).c_str())->Print();
 //      pdfconstrainslist_data_LP->add(*workspace->pdf(constrainslist_data_LP->at(i).c_str()));
 //    }
 //  }
 //
 //  RooFitResult* rfresult_data_LP = model_data_extremefail->fitTo(*rdataset_data_mj_extremefail,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data_LP),RooFit::Minimizer("Minuit2"));
 //  rfresult_data_LP = model_data_extremefail->fitTo(*rdataset_data_mj_extremefail,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data_LP),RooFit::Minimizer("Minuit2"));
 //  rfresult_data_LP->Print();
 //
 //  cut.Form("rrv_number_bkg_data%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
 //  RooRealVar* rrv_number_bkg_data_extremefailtau2tau1cut   = workspace->var(cut.Data());
 //  cut.Form("rrv_number_ttbar_data%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
 //  RooRealVar* rrv_number_ttbar_data_extremefailtau2tau1cut = workspace->var(cut.Data());
 //  rrv_number_ttbar_data_extremefailtau2tau1cut->setError((rrv_number_ttbar_data_extremefailtau2tau1cut->getVal()+rrv_number_bkg_data_extremefailtau2tau1cut->getVal())/2. );
 //  rrv_number_bkg_data_extremefailtau2tau1cut->setError((rrv_number_ttbar_data_extremefailtau2tau1cut->getVal()+rrv_number_bkg_data_extremefailtau2tau1cut->getVal())/2. );
 //
 //  model_STop_extremefail->Print();
 //  model_VV_extremefail->Print();
 //  model_WJets_extremefail->Print();
 //  rdataset_data_mj_extremefail->Print();
 //
 //  //extremefail fit of the mc
 //  std::vector<std::string>* constrainslist_TotalMC_LP = new std::vector<std::string>();
 //  RooAbsPdf*  model_ttbar_TotalMC_extremefail = MakeExtendedModel(workspace,"_ttbar_TotalMC"+label+"_"+"extremefailtau2tau1cut",mj_shape["mc_extremefail"],"_mj",channel,wtagger,constrainslist_TotalMC_LP);
 //  RooAbsPdf*  model_bkg_TotalMC_extremefail = MakeExtendedModel(workspace,"_bkg_TotalMC"+label+"_"+"extremefailtau2tau1cut",mj_shape["mc_bkg_extremefail"],"_mj",channel,wtagger,constrainslist_TotalMC_LP);
 //
 //  RooAddPdf* model_TotalMC_extremefail = new RooAddPdf(("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(),("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(), RooArgList(*model_ttbar_TotalMC_extremefail,*model_bkg_TotalMC_extremefail,*model_STop_extremefail,*model_VV_extremefail,*model_WJets_extremefail));
 //  workspace->import(*model_TotalMC_extremefail);
 //
 //  RooArgSet* pdfconstrainslist_TotalMC_LP = new RooArgSet(("pdfconstrainslist_TotalMC_LP_"+channel).c_str());
 //  if(constrainslist_TotalMC_LP!=NULL){
 //    for(unsigned int i =0; i < constrainslist_TotalMC_LP->size(); i++){
 //      workspace->pdf(constrainslist_TotalMC_LP->at(i).c_str())->Print();
 //      pdfconstrainslist_TotalMC_LP->add(*workspace->pdf(constrainslist_TotalMC_LP->at(i).c_str()));
 //    }
 //  }
 //
 //  RooFitResult* rfresult_TotalMC_LP = model_TotalMC_extremefail->fitTo(*rdataset_TotalMC_mj_extremefail, RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC_LP),RooFit::Minimizer("Minuit2"));
 //  rfresult_TotalMC_LP = model_TotalMC_extremefail->fitTo(*rdataset_TotalMC_mj_extremefail, RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC_LP),RooFit::Minimizer("Minuit2"));
 //  rfresult_TotalMC_LP->Print();
 //  cut.Form("rrv_number_bkg_TotalMC%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
 //  RooRealVar* rrv_number_bkg_TotalMC_extremefailtau2tau1cut   = workspace->var(cut.Data());
 //  cut.Form("rrv_number_ttbar_TotalMC%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
 //  RooRealVar* rrv_number_ttbar_TotalMC_extremefailtau2tau1cut = workspace->var(cut.Data());
 //  rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->setError((rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut->getVal())/2.);
 //  rrv_number_bkg_TotalMC_extremefailtau2tau1cut->setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut->getVal())/2. );
 //  model_STop_extremefail->Print();
 //  model_VV_extremefail->Print();
 //  model_WJets_extremefail->Print();
 //  rdataset_TotalMC_mj_extremefail->Print();
 //
 //  //LUCA
 //  if(TString(label).Contains("herwig"))
 //    model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("TTbar_herwig"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
 //  else
 //    model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("TTbar"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
 //
 //  cut.Form("%s,%s,%s",model_histpdf_STop_extremefail->GetName(),model_histpdf_VV_extremefail->GetName(),model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("STop"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());
 //  cut.Form("%s,%s",model_histpdf_VV_extremefail->GetName(), model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("VV"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());
 //  cut.Form("%s", model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("WJets"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());
 //
 //  //solid line
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
 //  cut.Form("%s,%s,%s",model_histpdf_STop_extremefail->GetName(), model_histpdf_VV_extremefail->GetName(), model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("STop_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
 //  cut.Form("%s,%s",model_histpdf_VV_extremefail->GetName(), model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("VV_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
 //  cut.Form("%s", model_histpdf_WJets_extremefail->GetName());
 //  model_TTbar_STop_VV_WJets_extremefail->plotOn(xframe_data_extremefail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_extremefail),RooFit::Name("WJets_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
 //  ///////////LUCA
 //
 //  rdataset_data_mj_extremefail->plotOn(xframe_data_extremefail, RooFit::Name("data"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
 //  rdataset_TotalMC_mj_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("TotalMC_invisible"), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
 //
 //  model_TotalMC_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("total MC fit"),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed));
 //  //    cut.Form("%s,%s,%s,%s",model_bkg_TotalMC_extremefail->GetName(),model_STop_extremefail->GetName(),model_VV_extremefail->GetName(),model_WJets_extremefail->GetName());
 //  // model_TotalMC_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("MC fit bkg"),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));
 //
 //
 //  rdataset_data_mj_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data_invisible"), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );
 //
 //  model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("total data fit"),RooFit::NormRange("controlsample_fitting_range"));
 //  model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data_fit_invisible"),RooFit::NormRange("controlsample_fitting_range"));
 //  cut.Form("%s,%s,%s,%s",model_bkg_data_extremefail->GetName(),model_STop_extremefail->GetName(),model_VV_extremefail->GetName(),model_WJets_extremefail->GetName());
 //  //    model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data fit bkg"),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed));
 //
 //  xframe_data_extremefail->GetYaxis()->SetRangeUser(0.,xframe_data_extremefail->GetMaximum()*1.3);
 //  TLegend* leg_data_extremefail = legend4Plot(xframe_data_extremefail,0,-0.2,0.07,0.04,0.,1,channel);
 //  xframe_data_extremefail->addObject(leg_data_extremefail);
 //
 //  nameDir.Form("plots_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/",channel.c_str(),wtagger.c_str(),wtagger.c_str());
 //  namePlot.Form("control%s_%s_%s_extremefail_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
 //  draw_canvas(xframe_data_extremefail,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);
}


///-------------------------------------------
void DrawScaleFactorTTbarControlSample(RooWorkspace* workspace, std::map<std::string,int> color_palet, const std::string & label, const std::string & channel, const std::string & wtagger,const double & ca8_ungroomed_pt_min, const double & ca8_ungroomed_pt_max){

  RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");

  RooDataSet* rdataset_data_mj  = (RooDataSet*) workspace->data(("rdataset_data"+label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_TTbar_mj = (RooDataSet*) workspace->data(("rdataset_TTbar"+label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_STop_mj  = (RooDataSet*) workspace->data(("rdataset_STop"+label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_VV_mj    = (RooDataSet*) workspace->data(("rdataset_VV"+label+"_"+channel+"_mj").c_str());
  RooDataSet* rdataset_WJets_mj = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+channel+"_mj").c_str());

  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_TTbar_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_STop_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_VV_mj);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_WJets_mj);

  RooAbsPdf* model_histpdf_STop  = workspace->pdf((std::string(rdataset_STop_mj->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_VV    = workspace->pdf((std::string(rdataset_VV_mj->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_WJets = workspace->pdf((std::string(rdataset_WJets_mj->GetName())+"_histpdf").c_str());

  RooAbsPdf* model_TTbar_STop_VV_WJets = workspace->pdf(("model_TTbar_STop_VV_WJets"+label+"_"+channel).c_str());

  //dataset fail tau2tau1 cut
  RooDataSet* rdataset_data_mj_fail   = (RooDataSet*) workspace->data(("rdataset_data"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_TTbar_mj_fail  = (RooDataSet*) workspace->data(("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_STop_mj_fail   = (RooDataSet*) workspace->data(("rdataset_STop"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_VV_mj_fail     = (RooDataSet*) workspace->data(("rdataset_VV"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
  RooDataSet* rdataset_WJets_mj_fail  = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+channel+"_mj").c_str());
 
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_TTbar_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_STop_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_VV_mj_fail);
  change_dataset_to_histpdf(workspace,rrv_mass_j, rdataset_WJets_mj_fail);

  RooAbsPdf* model_histpdf_STop_fail  = workspace->pdf((std::string(rdataset_STop_mj_fail->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_VV_fail    = workspace->pdf((std::string(rdataset_VV_mj_fail->GetName())+"_histpdf").c_str());
  RooAbsPdf* model_histpdf_WJets_fail = workspace->pdf((std::string(rdataset_WJets_mj_fail->GetName())+"_histpdf").c_str());

  RooAbsPdf* model_TTbar_STop_VV_WJets_fail = workspace->pdf(("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel).c_str());

  double scale_number_TTbar_STop_VV_WJets = (rdataset_TTbar_mj->sumEntries()+rdataset_STop_mj->sumEntries()+rdataset_VV_mj->sumEntries() +rdataset_WJets_mj->sumEntries() )/( rdataset_data_mj->sumEntries()+rdataset_data_mj_fail->sumEntries() );
  double scale_number_TTbar_STop_VV_WJets_fail = (rdataset_TTbar_mj_fail->sumEntries()+rdataset_STop_mj_fail->sumEntries()+rdataset_VV_mj_fail->sumEntries() +rdataset_WJets_mj_fail->sumEntries() )/( rdataset_data_mj->sumEntries()+rdataset_data_mj_fail->sumEntries() );

  RooAbsPdf* model_data_fail = workspace->pdf(("model_data"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str());
  RooAbsPdf* model_data = workspace->pdf(("model_data"+label+"_"+channel).c_str());

  RooCategory* category_p_f = workspace->cat(("category_p_f"+label+"_"+channel).c_str());

  RooSimultaneous* simPdf_data = new RooSimultaneous(("simPdf_data"+label+"_"+channel).c_str(),("simPdf_data"+label+"_"+channel).c_str(),*category_p_f);
  simPdf_data->addPdf(*model_data,"pass");
  simPdf_data->addPdf(*model_data_fail,"fail");
  RooDataSet* combData_p_f_data = (RooDataSet*) workspace->data(("combData_p_f_data"+label+"_"+channel).c_str());

  simPdf_data->Print();
  combData_p_f_data->Print();

  RooAbsPdf* model_TotalMC_fail = workspace->pdf(("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel).c_str());
  RooAbsPdf* model_TotalMC      = workspace->pdf(("model_TotalMC"+label+"_"+channel).c_str());

  RooSimultaneous* simPdf_TotalMC = new RooSimultaneous(("simPdf_TotalMC"+label+"_"+channel).c_str(),("simPdf_TotalMC"+label+"_"+channel).c_str(),*category_p_f);
  simPdf_TotalMC->addPdf(*model_TotalMC,"pass");
  simPdf_TotalMC->addPdf(*model_TotalMC_fail,"fail");
  RooDataSet* combData_p_f_TotalMC = (RooDataSet*) workspace->data(("combData_p_f_TotalMC"+label+"_"+channel).c_str());

  RooPlot* xframe_data = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));
  RooPlot* xframe_data_fail = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));


  TString cut ;
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );

  //#############################################################  
  ///plot total pass mc normalizing it to data, passing + failing                                                                                                                   
  //#############################################################  

  if(TString(label).Contains("herwig"))    
    model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("TTbar_herwig"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
  else
    model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
  
  cut.Form("%s,%s,%s",model_histpdf_STop->GetName(),model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("STop"),RooFit::Components(cut.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("VV"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s",model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("WJets"),RooFit::Components(cut.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

  //pass plots
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());  
  cut.Form("%s,%s,%s",model_histpdf_STop->GetName(), model_histpdf_VV->GetName(), model_histpdf_WJets->GetName());       
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("STop_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV->GetName(),model_histpdf_WJets->GetName());      
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("VV_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s",model_histpdf_WJets->GetName());
  model_TTbar_STop_VV_WJets->plotOn(xframe_data,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets),RooFit::Name("WJets_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

  // plot data again                                                                                                                                                                 
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data"), RooFit::Cut(cut.Data()), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));

  // plot mc fit function                                                                                                                                                            
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_TotalMC->plotOn(xframe_data,RooFit::Name("TotalMC_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible() , RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );

  simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("total MC fit"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed));     
  // cut.Form("model_bkg_TotalMC_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
  // simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  // plot data fit function
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());                                                    
  combData_p_f_data->plotOn(xframe_data,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  simPdf_data->plotOn(xframe_data,RooFit::Name("data fit"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  //    cut.Form("model_bkg_data_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
  //  simPdf_data->plotOn(xframe_data,RooFit::Name("data fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed));

  ///plot total fail mc normalizing it to data, passing + failing                                                                                                                   
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data_fail,RooFit::Name("data_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  if(TString(label).Contains("herwig"))
    model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar_herwig"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());
  else
    model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar"), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"+label]), RooFit::LineColor(kBlack), RooFit::VLines());

  cut.Form("%s,%s,%s",model_histpdf_STop_fail->GetName(),model_histpdf_VV_fail->GetName(),model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("STop"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());  
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("VV"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());
  cut.Form("%s", model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("WJets"),RooFit::Components(cut.Data()), RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

  //solid line
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s,%s",model_histpdf_STop_fail->GetName(), model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());       
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("STop_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s,%s",model_histpdf_VV_fail->GetName(), model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("VV_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
  cut.Form("%s", model_histpdf_WJets_fail->GetName());
  model_TTbar_STop_VV_WJets_fail->plotOn(xframe_data_fail,RooFit::Normalization(scale_number_TTbar_STop_VV_WJets_fail),RooFit::Name("WJets_line_invisible"),RooFit::Components(cut.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

  
  //fail plots -> plot MC fit                                                                                                                                          
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_TotalMC->plotOn(xframe_data_fail,RooFit::Name("TotalMC_invisible"), RooFit::Cut(cut.Data()), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));

  simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("MC fit"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed));

  //    cut.Form("model_bkg_TotalMC_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
  //    simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

  //fail plots -> plot data fit                                                                                                                                                       
  cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
  combData_p_f_data->plotOn(xframe_data_fail,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

  simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data fit"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data_fit_invisible"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
  //    cut.Form("model_bkg_data_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
  // simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()),RooFit::LineColor(kRed));
      
  //signal window
  TLine* lowerLine = new TLine(65,0.,65,xframe_data->GetMaximum()*0.7); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kBlack); lowerLine->SetLineStyle(9);
  TLine* upperLine = new TLine(105,0.,105,xframe_data->GetMaximum()*0.7); upperLine->SetLineWidth(2); upperLine->SetLineColor(kBlack); upperLine->SetLineStyle(9);
  xframe_data->addObject(lowerLine); xframe_data->addObject(upperLine);
  lowerLine = new TLine(65,0.,65,xframe_data_fail->GetMaximum()*0.7); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kBlack); lowerLine->SetLineStyle(9);
  upperLine = new TLine(105,0.,105,xframe_data_fail->GetMaximum()*0.7); upperLine->SetLineWidth(2); upperLine->SetLineColor(kBlack); upperLine->SetLineStyle(9);
  xframe_data_fail->addObject(lowerLine); xframe_data_fail->addObject(upperLine);

  TLegend* leg_data = legend4Plot(xframe_data,0,-0.2,0.07,0.04,0.,1,channel);                                                                                          
  xframe_data->addObject(leg_data);
  TLegend* leg_data_fail = legend4Plot(xframe_data_fail,0,-0.2,0.07,0.04,0.,1,channel);                                           
  xframe_data_fail->addObject(leg_data_fail);

  std::string tmp_channel = "el";
  if(workspace->var(("rrv_mean1_gaus_ttbar_data"+label+"_"+"el").c_str())) tmp_channel = "el";
  else tmp_channel = "mu";

  RooRealVar* rrv_mean_gaus_data     = workspace->var(("rrv_mean1_gaus_ttbar_data"+label+"_"+tmp_channel).c_str());
  RooRealVar* rrv_sigma_gaus_data    = workspace->var(("rrv_sigma1_gaus_ttbar_data"+label+"_"+tmp_channel).c_str());
  RooRealVar* rrv_mean_gaus_TotalMC  = workspace->var(("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel).c_str());
  RooRealVar* rrv_sigma_gaus_TotalMC = workspace->var(("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel).c_str());

  if(rrv_mean_gaus_TotalMC){
    TString latex ; latex.Form("Mean_{MC } = %3.1f #pm %2.1f",rrv_mean_gaus_TotalMC->getVal(),rrv_mean_gaus_TotalMC->getError()); 
    TLatex* tl_MC_mean  = new TLatex(0.25 ,0.62,latex.Data());
    latex.Form("Sigma_{MC }= %2.1f #pm %2.1f",rrv_sigma_gaus_TotalMC->getVal(),rrv_sigma_gaus_TotalMC->getError());
    TLatex* tl_MC_sigma = new TLatex(0.25 ,0.57,latex.Data());
    tl_MC_mean->SetNDC(); tl_MC_sigma->SetNDC();
    tl_MC_mean->SetTextSize(0.03);
    tl_MC_sigma->SetTextSize(0.03);
    // xframe_data.addObject(tl_MC_mean);
    // xframe_data.addObject(tl_MC_sigma);

  }
  if(rrv_mean_gaus_data){
    TString latex ; latex.Form("Mean_{data} = %3.1f #pm %2.1f",rrv_mean_gaus_data->getVal(),rrv_mean_gaus_data->getError()); 
    TLatex* tl_data_mean  = new TLatex(0.25 ,0.62,latex.Data());
    latex.Form("Sigma_{data}= %2.1f #pm %2.1f",rrv_sigma_gaus_data->getVal(),rrv_sigma_gaus_data->getError());
    TLatex* tl_data_sigma = new TLatex(0.25 ,0.57,latex.Data());
    tl_data_mean->SetNDC(); tl_data_sigma->SetNDC();
    tl_data_mean->SetTextSize(0.03);
    tl_data_sigma->SetTextSize(0.03);
    // xframe_data.addObject(tl_data_mean); xframe_data.addObject(tl_data_sigma);
  }

  xframe_data->GetYaxis()->SetRangeUser(1e-2,xframe_data->GetMaximum()*1.3);
  xframe_data_fail->GetYaxis()->SetRangeUser(1e-2,xframe_data_fail->GetMaximum()*1.3);

  TString nameDir;  nameDir.Form("plots_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/",channel.c_str(),wtagger.c_str(),wtagger.c_str());
  TString namePlot; namePlot.Form("control%s_%s_%s_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
  draw_canvas(xframe_data,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);    
  namePlot.Form("control%s_%s_%s_fail_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
  draw_canvas(xframe_data_fail,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);

}
