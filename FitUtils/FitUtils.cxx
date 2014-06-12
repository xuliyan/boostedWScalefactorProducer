#include "FitUtils.h"

void fit_genHMass(RooWorkspace* workspace, const std::string & label, const std::string & spectrum, const std::string & model, const std::string & channel, const std::string & wtagger, const int & ismc){

  std::cout<<"############## Fit genHMass "<<"  "<<label<<"  "<<model<<"  ##################"<<std::endl;

  RooRealVar* rrv_mass_gen_WW = workspace->var("rrv_mass_gen_WW");
  RooDataSet* rdataset = (RooDataSet*) workspace->data(("rdataset4fit"+label+"_genHMass_"+channel).c_str());
  std::vector<std::string>* constraint_list = new std::vector<std::string> ();

  RooAbsPdf* model_pdf = MakeExtendedModel(workspace,label+model,model,spectrum,channel,wtagger,constraint_list);
  rdataset->Print();

  RooFitResult* rfresult = model_pdf->fitTo(*rdataset,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE));
  rfresult = model_pdf->fitTo(*rdataset,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE),RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  model_pdf->Print();

  rfresult->SetName(("rfresult"+label+"_"+channel+"_mlvj").c_str());
  workspace->import(*rfresult);
  workspace->import(*model_pdf);
 
  RooPlot* mplot = rrv_mass_gen_WW->frame(RooFit::Title("genMHmass_{lvj} fitted by BWRUN"), RooFit::Bins(int(rrv_mass_gen_WW->getBins())));
  rdataset->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));
  model_pdf->plotOn(mplot);

  draw_canvas(mplot,("plots_"+channel+"_"+wtagger+"_g1/m_lvj_fitting/").c_str(),"genHiggsMass"+label,channel,GetLumi(),0,1,0);
}


void fit_mj_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & model,
		      const std::string & channel,const std::string & wtagger_label,const int & additionalinformation, const int & deco, const int & pseudodata){

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

     SystematicUncertaintyHiggs_01jetBin higgsUncertainty;

     RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
     if( not pseudodata) {
      while(param){
	  if(TString(label).Contains("VV") or TString(label).Contains("WW_EWK") or TString(label).Contains("STop") or TString(label).Contains("TTbar")){
	    if(TString(param->GetName()).Contains("rrv_mean1_gaus")){
               param->setRange(param->getMin()+higgsUncertainty.mean_mj_shift,param->getMax()+higgsUncertainty.mean_mj_shift);
               param->setVal(param->getVal()+higgsUncertainty.mean_mj_shift);

	    }
            if(TString(param->GetName()).Contains("rrv_deltamean_gaus")){
                    param->setRange(param->getMin()-higgsUncertainty.mean_mj_shift,param->getMax()-higgsUncertainty.mean_mj_shift);
                    param->setVal(param->getVal()-higgsUncertainty.mean_mj_shift);
	    }
            if(TString(param->GetName()).Contains("rrv_sigma1_gaus")){
                    param->setVal(param->getVal()*higgsUncertainty.sigma_mj_smear);
                    param->setRange(param->getMin()*higgsUncertainty.sigma_mj_smear,param->getMax()*higgsUncertainty.sigma_mj_smear);
	    }
            if(TString(param->GetName()).Contains("rrv_scalesigma_gaus")){
                    param->setRange(param->getMin()/higgsUncertainty.sigma_mj_smear,param->getMax()/higgsUncertainty.sigma_mj_smear);
                    param->setVal(param->getVal()/higgsUncertainty.sigma_mj_smear);
	    }
	  }
          param = dynamic_cast<RooRealVar*>(par.Next());

      }
     }
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
     RooPlot* mplot_pull = get_ratio(rrv_mass_j,rdataset_mj,model_pdf,rfresult,1,1);
     mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.2);

     //## CALCULATE CHI2
     RooDataHist* datahist = rdataset_mj->binnedClone((std::string(rdataset_mj->GetName())+"_binnedClone").c_str(),(std::string(rdataset_mj->GetName())+"_binnedClone").c_str());
     int Nbin = int(rrv_mass_j->getBins()); 
     RooArgList rresult_param = rfresult->floatParsFinal();        
     int nparameters =  rresult_param.getSize();                                         
     RooAbsReal* ChiSquare = model_pdf->createChi2(*datahist,RooFit::Extended(kTRUE),RooFit::DataError(RooAbsData::Poisson));
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
           
       std::cout<<"model"+label+"massvbf_jes_up"+model+"_"+channel+"_mj"<<"   "<<"rrv_scale_to_lumi"+label+"_"+channel<<std::endl;   
       if(workspace->pdf(("model"+label+"massvbf_jes_up"+model+"_"+channel+"_mj").c_str())) 
	  workspace->pdf(("model"+label+"massvbf_jes_up"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_up"), RooFit::LineColor(kRed),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jes_up"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));

       std::cout<<"model"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj"<<"   "<<"rrv_scale_to_lumi"+label+"_"+channel<<std::endl;   
       if(workspace->pdf(("model"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj").c_str())) 
          workspace->pdf(("model"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_dn"), RooFit::LineColor(kBlue),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));

       std::cout<<"model"+label+"massvbf_jer"+model+"_"+channel+"_mj"<<"   "<<"rrv_scale_to_lumi"+label+"_"+channel<<std::endl;   
       if(workspace->pdf(("model"+label+"massvbf_jer"+model+"_"+channel+"_mj").c_str()))
          workspace->pdf(("model"+label+"massvbf_jer"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("jer"), RooFit::LineColor(kAzure+1),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));
  
       std::cout<<"model+"+label+"massvbf_jer_up"+model+"_"+channel+"_mj"<<"   "<<"rrv_scale_to_lumi"+label+"_"+channel<<std::endl;   
       if(workspace->pdf(("model"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str()))
          workspace->pdf(("model"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_up"), RooFit::LineColor(kGreen+1),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));

       if(workspace->pdf(("model"+label+"massvbf_jer_dn"+model+"_"+channel+"_mj").c_str()))
          workspace->pdf(("model"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_dn"), RooFit::LineColor(6),RooFit::Normalization(workspace->var(("rrv_number"+label+"massvbf_jer_dn"+model+"_"+channel+"_mj").c_str())->getVal()/(rdataset_mj->sumEntries()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal())));

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

   if(deco){
      
      //############### comparison between only decorrelated shapes
      std::cout<<"################### Decorrelated mj single mc shape ################"<<std::endl;
      RooAbsPdf* model_pdf_temp = workspace->pdf(("model_pdf"+label+model+"_"+channel+"_mj").c_str()); 
      model_pdf_temp->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE));
      RooFitResult* rfresult_pdf = model_pdf_temp->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Minimizer("Minuit2"));
      RooWorkspace* wsfit_tmp = new RooWorkspace(("wsfit_tmp"+label+model+"_"+channel+"_mj").c_str());
      PdfDiagonalizer* Deco   = new PdfDiagonalizer(("Deco"+label+model+"_"+channel+"_"+wtagger_label+"_mj").c_str(),wsfit_tmp,*rfresult_pdf);
      std::cout<<"##################### diagonalize "<<std::endl;;
      RooAbsPdf* model_deco = Deco->diagonalize(*model_pdf_temp);
      std::cout<<"##################### original parameters "<<std::endl;
      model_pdf_temp->getParameters(rdataset_mj)->Print("v");
      std::cout<<"##################### original decorrelated parameters "<<std::endl;
      model_deco->getParameters(rdataset_mj)->Print("v");
      std::cout<<"##################### original pdf "<<std::endl;

      RooPlot* mplot_deco = rrv_mass_j->frame( RooFit::Title("plot deco mj"),RooFit::Bins(int(rrv_mass_j->getBins())));
      rdataset_mj->plotOn(mplot_deco, RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
      model_deco->plotOn(mplot_deco,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
      RooRealVar* rrv_number_dataset = new RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset_mj->sumEntries());
      rrv_number_dataset->setError(0.);
      draw_error_band(rdataset_mj,model_pdf_temp,rrv_number_dataset,rfresult_pdf,mplot_deco,1,"F");

      if(workspace->pdf(("model"+label+"massvbf_jes_up"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model"+label+"massvbf_jes_up"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("jes_up"), RooFit::LineColor(kRed));

      if(workspace->pdf(("model"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model"+label+"massvbf_jes_dn"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("jes_dn"), RooFit::LineColor(kBlue));

      if(workspace->pdf(("model"+label+"massvbf_jer"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model"+label+"massvbf_jer"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("jer"), RooFit::LineColor(kAzure+1));

      if(workspace->pdf(("model"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model"+label+"massvbf_jer_up"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("jer_up"), RooFit::LineColor(kGreen+1));

      if(workspace->pdf(("model"+label+"massvbf_jer_dn"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model"+label+"massvbf_jer_dn"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("jer_dn"), RooFit::LineColor(6));

      if(label == "_WJets0" and workspace->pdf(("model_WJets01"+model+"_"+channel+"_mj").c_str()))
         workspace->pdf(("model_WJets01"+model+"_"+channel+"_mj").c_str())->plotOn(mplot_deco,RooFit::Name("alt shape"), RooFit::LineColor(kOrange+1));

      TLegend* legend = legend4Plot(mplot_deco,0,0.,0.06,0.16,0.,1,channel);

      mplot_deco->addObject(legend);
      rdataset_mj->plotOn(mplot_deco,RooFit::Name("MC Events"),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));
      model_deco->plotOn(mplot_deco,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
      mplot_deco->GetYaxis()->SetRangeUser(0,mplot_deco->GetMaximum()*1.2);

      command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
      draw_canvas_with_pull(mplot_deco,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),"m_j_shape"+label+fileName,model,channel,0,0,GetLumi());
     
     }

  
     workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setVal(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());
     workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());

    
}


//## Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
void fit_mlvj_model_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & in_range, const std::string & mlvj_model, const std::string & channel, const std::string & wtagger_label, const int & deco, const int & show_constant_parameter, const int & logy, const int & ismc, const std::string & label_origin, const std::string & FTest_output_txt){

      std::cout<<"############### Fit mlvj single MC sample "<<fileName<<" "<<label<<" "<<mlvj_model<<" "<<in_range<<" deco "<<deco<<" ##################"<<std::endl;
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
      RooPlot* mplot_pull = get_ratio(rrv_mass_lvj,rdataset,model,rfresult,1,1);
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

      if(deco){
       std::cout<<"################### Decorrelated mlvj single mc shape ################"<<std::endl;
       RooAbsPdf* model_pdf = workspace->pdf(("model_pdf"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str());
       std::cout<<"model_pdf"+label+in_range+mlvj_model+"_"+channel+"_mlvj"<<std::endl; 
       model_pdf->fitTo(*rdataset,RooFit::Save(1),RooFit::SumW2Error(kTRUE));
       RooFitResult* rfresult_pdf = model_pdf->fitTo(*rdataset, RooFit::Save(1), RooFit::SumW2Error(kTRUE), RooFit::Minimizer("Minuit2"));

       RooWorkspace* wsfit_tmp = new RooWorkspace(("wsfit_tmp"+label+in_range+"_"+channel+"_mlvj").c_str());
       PdfDiagonalizer* Deco   = new PdfDiagonalizer(("Deco"+label+in_range+mlvj_model+"_"+channel+"_"+wtagger_label+"_mlvj").c_str(),wsfit_tmp,*rfresult_pdf); 
       std::cout<<"##################### diagonalize "<<std::endl;;
       RooAbsPdf* model_pdf_deco = Deco->diagonalize(*model_pdf); 
       std::cout<<"##################### workspace for decorrelation "<<std::endl;
       wsfit_tmp->Print("v");
       std::cout<<"##################### original parameters "<<std::endl;
       model_pdf->getParameters(rdataset)->Print("v"); 
       std::cout<<"##################### original decorrelated parameters "<<std::endl;
       model_pdf_deco->getParameters(rdataset)->Print("v");
	     
       if(label == "_TTbar"){
	  
        RooPlot* mplot_deco = rrv_mass_lvj->frame( RooFit::Bins(int(rrv_mass_lvj->getBins())));
	rdataset->plotOn(mplot_deco, RooFit::Name("Powheg Sample"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );
	model_pdf_deco->plotOn(mplot_deco,RooFit::Name("TTbar_Powheg"),RooFit::LineColor(kBlack));
	RooRealVar* rrv_number_dataset = new  RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset->sumEntries());
	rrv_number_dataset->setError(0.); 
	if(mlvj_model != "CBBW_v1") draw_error_band(rdataset,model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,1,"F"); 

	if(workspace->pdf(("model"+label+"_mcanlo"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"_mcanlo"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_mcanlo"), RooFit::LineColor(kBlue));

	if(workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_jes_up"), RooFit::LineColor(kRed));

	if(workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_jes_dn"), RooFit::LineColor(kAzure+1));

	if(workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_jer"), RooFit::LineColor(kGreen+1));

	if(workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_jer_up"), RooFit::LineColor(6));

	if(workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
           workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("TTbar_jer_dn"), RooFit::LineColor(kRed));

        TLegend* leg = legend4Plot(mplot_deco,0,-0.07,0.02,0.01,0.,1,channel);
	mplot_deco->addObject(leg);
	rdataset->plotOn(mplot_deco, RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));

	model_pdf->plotOn(mplot_deco,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
	mplot_deco->GetYaxis()->SetRangeUser(0,mplot_deco->GetMaximum()*1.2);

        command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
        draw_canvas_with_pull(mplot_deco,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),"m_lvj_shape"+label+in_range,mlvj_model,channel,0,logy,GetLumi());
           
       }

       else{

  	 RooPlot* mplot_deco = rrv_mass_lvj->frame( RooFit::Bins(int(rrv_mass_lvj->getBins())));
         rdataset->plotOn(mplot_deco, RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );
         model_pdf->plotOn(mplot_deco,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
	 RooRealVar* rrv_number_dataset = new RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset->sumEntries());
         rrv_number_dataset->setError(0.);
 	 if(mlvj_model != "CBBW_v1") draw_error_band(rdataset,model_pdf,rrv_number_dataset,rfresult_pdf,mplot_deco,1,"F"); 
           
	 if(workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())) 
            workspace->pdf(("model"+label+"massvbf_jes_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("jes_up"), RooFit::LineColor(kRed));

	 if(workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())) 
            workspace->pdf(("model"+label+"massvbf_jes_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("jes_dn"), RooFit::LineColor(kBlue));

	 if(workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())) 
            workspace->pdf(("model"+label+"massvbf_jer"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("jer"), RooFit::LineColor(kAzure+1));

	 if(workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())) 
            workspace->pdf(("model"+label+"massvbf_jer_up"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("jer_up"), RooFit::LineColor(kGreen+1));

	 if(workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())) 
            workspace->pdf(("model"+label+"massvbf_jer_dn"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("jer_dn"), RooFit::LineColor(6));

         if(label == "_WJets0" and workspace->pdf(("model_WJets01"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str()))
            workspace->pdf(("model_WJets01"+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->plotOn(mplot_deco,RooFit::Name("alt shape"), RooFit::LineColor(kOrange+1));

         TLegend* leg = legend4Plot(mplot_deco,0,-0.07,0.02,0.01,0.,1,channel);
         mplot_deco->addObject(leg);
         rdataset->plotOn(mplot_deco, RooFit::Name("MC Events"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );
         model_pdf->plotOn(mplot_deco,RooFit::Name("Nominal MC"),RooFit::LineColor(kBlack));
         mplot_deco->GetYaxis()->SetRangeUser(0,mplot_deco->GetMaximum()*1.2);

         command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger_label.c_str());
         draw_canvas_with_pull(mplot_deco,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),"m_lvj_shape"+label+in_range,mlvj_model,channel,0,logy,GetLumi());
       }


       //## Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
       std::cout<<"rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj"<<std::endl;
       std::cout<<"rrv_scale_to_lumi"+label+"_"+channel+in_range+"_mlvj"<<std::endl;;
        
       workspace->import(*model_pdf_deco);

      }
      
      std::cout<<"rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj"<<std::endl;
      std::cout<<"rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj"<<std::endl;
      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setVal(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal());
      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setError(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label_origin+"_"+channel+in_range+"_mlvj").c_str())->getVal() );
      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->Print();
      

}

//#### make the mj sideband fit on data ti get the Wjets normaliztion
void fit_WJetsNormalization_in_Mj_signal_region(RooWorkspace* workspace,  std::map<std::string,int> color_palet, std::map<std::string,std::string> mj_shape, const std::string & label, const std::string & mass_scale, const std::string & model, const std::string & channel, const std::string wtagger, const int & ttbarcontrolregion, const int & pseudodata, const float & mj_signal_min, const float & mj_signal_max, const std::string & jetBin){

   std::cout<<"############### Fit mj Normalization: "<<label<<" ##################"<<std::endl;
   RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");
   TString Name ; 
   Name.Form("rdataset_data_%s_mj",channel.c_str());
   RooDataSet* rdataset_data_mj = (RooDataSet*) workspace->data(Name.Data());

   //## Fix TTbar, VV and STop
   RooAbsPdf* model_WJets  = NULL;
   RooAbsPdf* model_TTbar  = NULL;
   RooAbsPdf* model_STop   = NULL;
   RooAbsPdf* model_VV     = NULL;
   RooAbsPdf* model_WW_EWK = NULL;

   std::cout<<" shape "<<mj_shape["TTbar"]<<" "<<mj_shape["STop"]<<" "<<mj_shape["VV"]<<"  "<<mj_shape[label]<<std::endl;

   if(ttbarcontrolregion){
         model_WJets  = get_WJets_mj_Model(workspace,"_WJets0"+mass_scale,mj_shape["WJets0"],channel,1,jetBin);
         model_STop   = get_STop_mj_Model(workspace,"_STop"+mass_scale,mj_shape["STop"],channel);
         model_VV     = get_VV_mj_Model(workspace,"_VV"+mass_scale,mj_shape["VV"],channel);
         if( jetBin == "_2jet")   model_WW_EWK = get_WW_EWK_mj_Model(workspace,"_WW_EWK"+mass_scale,mj_shape["WW_EWK"],channel);
	 model_TTbar  = get_TTbar_mj_Model(workspace,label,model,channel,0);
   }
   else{
         model_TTbar  = get_TTbar_mj_Model(workspace,"_TTbar"+mass_scale,mj_shape["TTbar"],channel);
         model_STop   = get_STop_mj_Model(workspace,"_STop"+mass_scale,mj_shape["STop"],channel);
         model_VV     = get_VV_mj_Model(workspace,"_VV"+mass_scale,mj_shape["VV"],channel);
         if( jetBin == "_2jet")   model_WW_EWK = get_WW_EWK_mj_Model(workspace,"_WW_EWK"+mass_scale,mj_shape["WW_EWK"],channel);
	 model_WJets  = get_WJets_mj_Model(workspace,label,model,channel,0,jetBin);
  }

  //## Total Pdf and fit only in sideband
  RooAddPdf* model_data = NULL ;
  if( jetBin =="_2jet") 
   model_data = new RooAddPdf(("model_data"+mass_scale+"_"+channel+"_mj").c_str(),("model_data"+mass_scale+"_"+channel+"_mj").c_str(),RooArgList(*model_WJets,*model_VV,*model_WW_EWK,*model_TTbar,*model_STop));
  else
   model_data = new RooAddPdf(("model_data"+mass_scale+"_"+channel+"_mj").c_str(),("model_data"+mass_scale+"_"+channel+"_mj").c_str(),RooArgList(*model_WJets,*model_VV,*model_TTbar,*model_STop));
                                  
  RooFitResult* rfresult = model_data->fitTo(*rdataset_data_mj,RooFit::Save(1),RooFit::Range("sb_lo,sb_hi") ,RooFit::Extended(kTRUE),RooFit::NumCPU(4));
  rfresult = model_data->fitTo(*rdataset_data_mj,RooFit::Save(1),RooFit::Range("sb_lo,sb_hi"),RooFit::Extended(kTRUE), RooFit::NumCPU(4),RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();
  workspace->import(*model_data);

  //## Total numver of event --> full propagation of error due to all the background sources coming from the fit
  RooRealVar* rrv_number_data_mj = NULL ;

  if (ttbarcontrolregion){
     
    if( jetBin =="_2jet"){
    
     rrv_number_data_mj = new RooRealVar(("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),
                                   workspace->var(("rrv_number+"+label+model+"_"+channel+" _mj").c_str())->getVal()+ //##TTbar
                                   workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getVal()+  //##STop
                                   workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getVal()+    //## VV
                                   workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getVal()+ //## WW_EWK
                                   workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getVal());  //## WJets

     rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getError()+       
                                            workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getError()));         
    }
    else{
    
    rrv_number_data_mj = new RooRealVar(("rrv_number_data_+"+channel+"_mj").c_str(),("rrv_number_data_"+channel+"_mj").c_str(),
                                   workspace->var(("rrv_number+"+label+model+"_"+channel+" _mj").c_str())->getVal()+ //##TTbar
				   workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()*
                                   workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getVal()+    //## VV
                                   workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getVal());  //## WJets

    rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_WJets0"+mass_scale+mj_shape["WJets0"]+"_"+channel+"_mj").c_str())->getError()));         

    }

  workspace->import(*rrv_number_data_mj);

  std::cout<<"TTbar  events: "<<workspace->var(("rrv_number_"+label+model+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"STop   events: "<<workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"VV     events: "<<workspace->var(("rrv_number_VV"+mass_scale+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  if(jetBin == "_2jet") std::cout<<"WW_EWK events: "<<workspace->var(("rrv_number_WW_EWK"+mass_scale+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"WJets  events: "<<workspace->var(("rrv_number_WJets0"+mass_scale+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"Data   events: "<<workspace->var(("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str())->getVal()<<std::endl;

  }
  else{
         

    if( jetBin =="_2jet"){
    std::cout<<"rrv_number"+label+model+"_"+channel+" _mj"<<std::endl;
    std::cout<<"rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj"<<std::endl;
    std::cout<<"rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj"<<std::endl;
    std::cout<<"rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj"<<std::endl;
    std::cout<<"rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj"<<std::endl;

     rrv_number_data_mj = new RooRealVar(("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),
                                         workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getVal()+ //## TTbar
                                         workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getVal()+  //## STop
                                         workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getVal()+    //## VV
                                         workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getVal()+ //## WW_EWK
                                         workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal());  //## WJets

     rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getError()+       
                                              workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
					      workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()));         
    }

    else{
     rrv_number_data_mj = new RooRealVar(("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str(),
                                          workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getVal()+ //## TTbar
                                          workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getVal()+  //## STop
                                          workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getVal()+    //## VV
                                          workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal());  //## WJets

     rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()*
                                              workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getError()+
                                              workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
					      workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()));         
    }

   workspace->import(*rrv_number_data_mj);

   std::cout<< "TTbar  events: ",workspace->var(("rrv_number_TTbar"+mass_scale+mj_shape["TTbar"]+"_"+channel+"_mj").c_str())->getVal();
   std::cout<< "STop   events: ",workspace->var(("rrv_number_STop"+mass_scale+mj_shape["STop"]+"_"+channel+"_mj").c_str())->getVal();
   std::cout<< "VV     events: ",workspace->var(("rrv_number_VV"+mass_scale+mj_shape["VV"]+"_"+channel+"_mj").c_str())->getVal();
   if(jetBin == "_2jet") std::cout<< "WW_EWK events: ",workspace->var(("rrv_number_WW_EWK"+mass_scale+mj_shape["WW_EWK"]+"_"+channel+"_mj").c_str())->getVal();
   std::cout<< "WJets  events: ",workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal();
   std::cout<< "Data   events: ",workspace->var(("rrv_number_data"+mass_scale+"_"+channel+"_mj").c_str())->getVal();
  }

  //## draw the plot for the default WJets Shape
  RooPlot* mplot = rrv_mass_j->frame(RooFit::Title("rrv_mass_j"), RooFit::Bins(int(rrv_mass_j->getBins())));
  rdataset_data_mj->plotOn(mplot,RooFit::Name("data_invisible"), RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::Poisson), RooFit::XErrorSize(0),RooFit::Invisible());
  //## plot solid style
  if (ttbarcontrolregion){
    // plot solid style 
    if(jetBin == "_2jet"){
      Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar"),RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["_TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
   
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WW_EWK"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // plot "dashed" style area
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WW_EWK_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
           
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     /// solid line
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

          
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // dash line
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
    }

    else{

      Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar"),RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["_TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
   
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // plot "dashed" style area
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     /// solid line
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
          
     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());


     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // dash line
     Name.Form("model_WJets0%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model_WJets0%s%s_%s_mj,model_STop%s%s_%s_mj,model_VV%s%s_%s_mj,model%s_%s_mj",mass_scale.c_str(),mj_shape["WJets0"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
    }

  }

  else{

    if(jetBin == "_2jet"){

      Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WW_EWK"),RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
   
     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // plot "dashed" style area
     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WW_EWK_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
           
     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     /// solid line
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());


     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());


     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // dash line
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj,model_WW_EWK%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["WW_EWK"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
    }

    else {

      Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

      model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     // plot "dashed" style area
     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
           
     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("STop_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

     model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
     /// solid line
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());


     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());


     // dash line
     Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

     Name.Form("model%s_%s_mj,model_STop%s%s_%s_mj,model_TTbar%s%s_%s_mj,model_VV%s%s_%s_mj",(label+model).c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["STop"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["TTbar"].c_str(),channel.c_str(),mass_scale.c_str(),mj_shape["VV"].c_str(),channel.c_str());
     model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) , RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()), RooFit::LineStyle(kDashed), RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    }

  }                                                                  
        
  /// draw the error band using the sum of all the entries component MC + fit and the total error == Normalization for the fixed MC, shape + normalization for W+jets
  draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,color_palet["Uncertainty"],"F");
  model_data->plotOn(mplot,RooFit::Name("model_mc"),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"),RooFit::Invisible());

  if (pseudodata == 1) rdataset_data_mj->plotOn( mplot ,RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Name("data"));
  else GetDataPoissonInterval(rdataset_data_mj,rrv_mass_j,mplot);
                                       		
  // Get the pull and plot it
  RooPlot* mplot_pull = get_ratio(rrv_mass_j,rdataset_data_mj,model_data,rfresult,1,1);
  
  // signal window zone with vertical lines
  TLine* lowerLine = new TLine(mj_signal_min,0.,mj_signal_min,mplot->GetMaximum()*0.9); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kGray+2); lowerLine->SetLineStyle(9);
  TLine* upperLine = new TLine(mj_signal_max,0.,mj_signal_max,mplot->GetMaximum()*0.9); upperLine->SetLineWidth(2); upperLine->SetLineColor(kGray+2); upperLine->SetLineStyle(9);
  mplot->addObject(lowerLine);
  mplot->addObject(upperLine);

  // legend of the plot
  TLegend* leg = legend4Plot(mplot,0,-0.2,0.07,0.04,0.,1,channel);
  mplot->addObject(leg);
  mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.5);

  // CALCULATE CHI2
  RooDataHist* datahist = rdataset_data_mj->binnedClone((std::string(rdataset_data_mj->GetName())+"_binnedClone").c_str(),(std::string(rdataset_data_mj->GetName())+"_binnedClone").c_str());
  int Nbin = int(rrv_mass_j->getBins()); 
  RooArgList rresult_param = rfresult->floatParsFinal();        
  int nparameters =  rresult_param.getSize();                                         
  RooAbsReal* ChiSquare = model_data->createChi2(*datahist,RooFit::Extended(kTRUE),RooFit::DataError(RooAbsData::Poisson));
  float chi_over_ndf = ChiSquare->getVal()/(Nbin - nparameters);

  // Add Chisquare to mplot_pull
  Name.Form("#chi^{2}/ndf = %0.2f ",float(chi_over_ndf));
  TLatex* cs  = new TLatex(0.75,0.8,Name.Data());
  cs->SetNDC();
  cs->SetTextSize(0.12);
  cs->AppendPad("same");
  mplot_pull->addObject(cs);

  RooArgSet* parameters_list = model_data->getParameters(rdataset_data_mj);
   
  TString command; command.Form("mkdir -p plots_%s_%s_g1/mj_fitting_%s/",channel.c_str(),wtagger.c_str(),model.c_str()); 
  system(command.Data());

  command.Form("plots_%s_%s_g1/mj_fitting_%s/",channel.c_str(),wtagger.c_str(),model.c_str()); 

  Name.Form("m_j_sideband%s",(label+model).c_str());
  draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),std::string(Name.Data()),model,channel,0,1,GetLumi());

  //## print the different sample normalizations

  get_mj_normalization_insignalregion(workspace,"_data"+mass_scale,"",channel);
  get_mj_normalization_insignalregion(workspace,"_TTbar"+mass_scale,mj_shape["TTbar"],channel);
  get_mj_normalization_insignalregion(workspace,"_STop"+mass_scale,mj_shape["STop"],channel);
  get_mj_normalization_insignalregion(workspace,"_VV"+mass_scale,mj_shape["VV"],channel);
  if(jetBin == "_2jet") get_mj_normalization_insignalregion(workspace,"_WW_EWK"+mass_scale,mj_shape["WW_EWK"],channel);
  get_mj_normalization_insignalregion(workspace,label,model,channel);
           
                
  // to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data
  if(ttbarcontrolregion){

    RooAbsReal* fullInt   = model_TTbar->createIntegral(*rrv_mass_j,*rrv_mass_j);
    RooAbsReal* signalInt = model_TTbar->createIntegral(*rrv_mass_j,*rrv_mass_j,("signal_region"));
    double fullInt_val   = fullInt->getVal();
    double signalInt_val = signalInt->getVal()/fullInt_val;
    // take the value from the fit (normalization) and multiply it from the ratio of the integrals
    RooRealVar* rrv_number_TTbar_in_mj_signal_region_from_fitting = new RooRealVar(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*signalInt_val);

    // Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
    rrv_number_TTbar_in_mj_signal_region_from_fitting->setError(Calc_error_extendPdf(rdataset_data_mj,dynamic_cast<RooExtendPdf*>(model_TTbar),rfresult,"signal_region"));
    std::cout<< "########## error on the normaliztion due to shape + norm = "<<rrv_number_TTbar_in_mj_signal_region_from_fitting->getError()<<std::endl;;
    workspace->import(*rrv_number_TTbar_in_mj_signal_region_from_fitting);
    rrv_number_TTbar_in_mj_signal_region_from_fitting->Print();

    RooRealVar* rrv_TTbar  = workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str()); //nominal parametrization for Wjets                        
    rrv_TTbar->Print();
  }
  else{
         RooAbsReal* fullInt   = model_WJets->createIntegral(*rrv_mass_j,*rrv_mass_j);
         RooAbsReal* signalInt = model_WJets->createIntegral(*rrv_mass_j,*rrv_mass_j,("signal_region"));
         double fullInt_val   = fullInt->getVal();
         double signalInt_val = signalInt->getVal()/fullInt_val;
         //take the value from the fit (normalization) and multiply it from the ratio of the integrals
         RooRealVar* rrv_number_WJets_in_mj_signal_region_from_fitting = new RooRealVar(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*signalInt_val);

         // Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
         rrv_number_WJets_in_mj_signal_region_from_fitting->setError( Calc_error_extendPdf(rdataset_data_mj, dynamic_cast<RooExtendPdf*>(model_WJets), rfresult,"signal_region") );
	 std::cout<< "########## error on the normaliztion due to shape + norm = "<<rrv_number_WJets_in_mj_signal_region_from_fitting->getError()<<std::endl;
         workspace->import(*rrv_number_WJets_in_mj_signal_region_from_fitting);
         rrv_number_WJets_in_mj_signal_region_from_fitting->Print();

         RooRealVar* rrv_WJets0  = workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str()); //nominal parametrization for Wjets                        
         rrv_WJets0->Print();
   }

}


//##### Counting of the events of each component in the signal region taking the lavel for the model
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


//Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape                                                                           
void fit_mlvj_in_Mj_sideband(RooWorkspace* workspace, std::map<std::string,int> color_palet, std::map<std::string,std::string> mlvj_shape, const std::string & label, const std::string & mass_scale, const std::string & mlvj_region, const std::string & mlvj_model, const std::string & channel, const std::string & wtagger, const int & ttbarcontrolregion, const int & logy, const int & pseudodata, const std::string & jetBin, const std::string & FTest_output_txt){

  std::cout<<"############### Fit mlvj in mj sideband: "<<label<<" "<<mlvj_region<<" "<<mlvj_model<<" ##################"<<std::endl;

  RooRealVar* rrv_mass_lvj = workspace->var("rrv_mass_lvj");
  RooDataSet* rdataset_data_mlvj = (RooDataSet*) workspace->data(("rdataset_data"+mlvj_region+"_"+channel+"_mlvj").c_str());

  // get and fix the minor component shapes in the sb low                                                                                                                            
  RooAbsPdf* model_VV_backgrounds     = dynamic_cast<RooAbsPdf*>(get_VV_mlvj_Model(workspace,"_VV"+mass_scale,mlvj_region,mlvj_shape["VV"],channel));
  RooAbsPdf* model_STop_backgrounds   = dynamic_cast<RooAbsPdf*>(get_STop_mlvj_Model(workspace,"_STop"+mass_scale,mlvj_region,mlvj_shape["STop"],channel));
  RooAbsPdf* model_WW_EWK_backgrounds = NULL ;
  if(jetBin == "_2jet") model_WW_EWK_backgrounds = dynamic_cast<RooAbsPdf*>(get_WW_EWK_mlvj_Model(workspace,"_WW_EWK"+mass_scale,mlvj_region,mlvj_shape["WW_EWK"],channel));

  RooAbsPdf* model_TTbar_backgrounds = NULL ; 
  RooAbsPdf* model_WJets_backgrounds = NULL ; 

  if(ttbarcontrolregion == 0)
    model_TTbar_backgrounds  = dynamic_cast<RooAbsPdf*>(get_TTbar_mlvj_Model(workspace,"_TTbar"+mass_scale,mlvj_region,mlvj_shape["TTbar"],channel));
  else
    model_WJets_backgrounds  = dynamic_cast<RooAbsPdf*>(get_WJets_mlvj_Model(workspace,"_WJets0"+mass_scale,mlvj_region,mlvj_shape["WJets0"],channel));

  workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->Print();
  if(jetBin == "_2jet") workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->Print();

  //Make the Pdf for the WJets                                                                                                                                                     
  std::vector<std::string>*   constraint = NULL;
  RooRealVar*   number_WJets_sb_lo = NULL; 
  RooRealVar*   number_TTbar_signal_region = NULL;
  RooExtendPdf* model_WJets = NULL; 
  RooExtendPdf* model_TTbar = NULL ;
  RooAddPdf*    model_data  = NULL ;

  RooAbsPdf*  model_pdf_WJets = NULL ; 
  RooAbsPdf*  model_pdf_TTbar = NULL ; 

  if(ttbarcontrolregion == 0){

    model_pdf_WJets = MakeGeneralPdf(workspace,label+mlvj_region+mlvj_model+std::string("_from_fitting"),mlvj_model,"_mlvj",wtagger,channel,constraint);
    model_pdf_WJets->Print();
    //inititalize the value to what was fitted with the mc in the sideband
    if (TString(label).Contains("WJets01"))                                                                                                          
     number_WJets_sb_lo = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+mlvj_shape["WJets01"]+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));
    else if (TString(label).Contains("WJets0"))
     number_WJets_sb_lo = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+mlvj_shape["WJets0"]+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));
    else if (TString(label).Contains("WJets1"))
     number_WJets_sb_lo = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+mlvj_shape["WJets1"]+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));

    model_WJets = new RooExtendPdf(("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),*model_pdf_WJets,*number_WJets_sb_lo);

    // Add the other bkg component fixed to the total model --> in the extended way
    if(jetBin == "_2jet") 
     model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_WJets,*model_VV_backgrounds,*model_TTbar_backgrounds,*model_STop_backgrounds,*model_WW_EWK_backgrounds));
    else
     model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_WJets,*model_VV_backgrounds,*model_TTbar_backgrounds,*model_STop_backgrounds));

    model_data->Print();
    
  }
  else{

         model_pdf_TTbar = MakeGeneralPdf(workspace,label+mlvj_region+mlvj_model+"_from_fitting",mlvj_model,"_mlvj",wtagger,channel,constraint);
         model_pdf_TTbar->Print();
         // inititalize the value to what was fitted with the mc in the sideband
         number_TTbar_signal_region = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));
    
         model_TTbar     = new  RooExtendPdf(("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),*model_pdf_TTbar,*number_TTbar_signal_region);

         number_TTbar_signal_region->Print();

         // Add the other bkg component fixed to the total model --> in the extended way
         if(jetBin == "_2jet") 
          model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_TTbar,*model_VV_backgrounds,*model_WJets_backgrounds,*model_STop_backgrounds,*model_WW_EWK_backgrounds));
         else
          model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_TTbar,*model_VV_backgrounds,*model_WJets_backgrounds,*model_STop_backgrounds));
        
  }
  
  RooFitResult* rfresult = model_data->fitTo(*rdataset_data_mlvj, RooFit::Save(1) ,RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE));
  rfresult = model_data->fitTo(*rdataset_data_mlvj, RooFit::Save(1) ,RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();
  workspace->import(*model_data);

  RooRealVar* rrv_number_data_sb_lo_mlvj  = NULL ;
  
  if(ttbarcontrolregion == 0){
    model_WJets->getParameters(*rdataset_data_mlvj)->Print("v");

    // data in the sideband plus error from fit        
    if(jetBin == "_2jet"){
     rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                ("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

     rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()));

    }
    else{
      rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                ("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

      rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_TTbar"+mass_scale+mlvj_region+mlvj_shape["TTbar"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()));


    }

    rrv_number_data_sb_lo_mlvj->Print();
    workspace->import(*rrv_number_data_sb_lo_mlvj);
  }

  else{

    model_TTbar->getParameters(*rdataset_data_mlvj)->Print("v");


    // data in the sideband plus error from fit        
    if(jetBin == "_2jet"){
      rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                ("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

     rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_WW_EWK"+mass_scale+mlvj_region+mlvj_shape["WW_EWK"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()));

    }
    else{
     rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                ("rrv_number_data"+mass_scale+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getVal()+
                                                workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

     rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_WJets0"+mass_scale+mlvj_region+mlvj_shape["WJets0"]+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_STop"+mass_scale+mlvj_region+mlvj_shape["STop"]+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_VV"+mass_scale+mlvj_region+mlvj_shape["VV"]+"_"+channel+"_mlvj").c_str())->getError()));

    }
    rrv_number_data_sb_lo_mlvj->Print();
    workspace->import(*rrv_number_data_sb_lo_mlvj);

  }

  RooPlot* mplot = NULL; 
  //plot for WJets default + default shape
  if(TString(label).Contains("_WJets")){

      mplot = rrv_mass_lvj->frame(RooFit::Title("M_lvj fitted in M_j sideband"),RooFit::Bins(int(rrv_mass_lvj->getBins())));

      if(jetBin != "_2jet"){
       rdataset_data_mlvj->plotOn(mplot,RooFit::Invisible(),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible(),RooFit::Name("data_invisible"));

       model_data->plotOn(mplot, RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("WJets"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"),RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components(std::string(model_STop_backgrounds->GetName()).c_str()), RooFit::Name("STop"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::VLines());

      //solid line
      model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("WJets_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName()).c_str())), RooFit::Name("STop_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      }
      else{

       rdataset_data_mlvj->plotOn(mplot,RooFit::Invisible(),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible(),RooFit::Name("data_invisible"));

       model_data->plotOn(mplot, RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WJets"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

       model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"),RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components(std::string(model_STop_backgrounds->GetName()).c_str()), RooFit::Name("STop"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::VLines());

      //solid line
      model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WJets_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName()).c_str())), RooFit::Name("STop_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      }
  }

  //plot for WJets default + default shape
  if(label == "_TTbar"){

            mplot = rrv_mass_lvj->frame(RooFit::Title("M_lvj fitted in TTbar Control Region "), RooFit::Bins(int(rrv_mass_lvj->getBins())));

            if(jetBin == "_2jet"){
             rdataset_data_mlvj->plotOn(mplot,RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Invisible(), RooFit::Name("data_invisible") );

             model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar->GetName())+","+std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("STop"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());


             //solid line
             model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar->GetName())+","+std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("WJets_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("STop_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());
	    }
            else{

             rdataset_data_mlvj->plotOn(mplot,RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Invisible(), RooFit::Name("data_invisible") );

             model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar->GetName())+","+std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("STop"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack), RooFit::VLines());


             //solid line
             model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar->GetName())+","+std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()),RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("WJets_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

             model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("STop_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

	    }

  }

   //draw the error band
   model_data->plotOn( mplot , RooFit::VLines(), RooFit::Invisible());
   model_data->plotOn( mplot , RooFit::Invisible(), RooFit::Name("model_mc"));
   if(pseudodata == 1) rdataset_data_mlvj->plotOn( mplot , RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Name("data"));
   else GetDataPoissonInterval(rdataset_data_mlvj,rrv_mass_lvj,mplot);

   draw_error_band(rdataset_data_mlvj,model_data,rrv_number_data_sb_lo_mlvj,rfresult,mplot,color_palet["Uncertainty"],"F");

   mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.2);
            
   //Add the legend to the plot
   TLegend* leg = legend4Plot(mplot,0,0.,0.06,0.16,0.,1,channel);
   mplot->addObject(leg);

   //### get the pull plot and store the canvas
   RooPlot* mplot_pull = get_ratio(rrv_mass_lvj,rdataset_data_mlvj,model_data,rfresult,1,1);
   RooArgSet* parameters_list = model_data->getParameters(rdataset_data_mlvj);

   //CALCULATE CHI2                                                                                                                                               
   RooDataHist* datahist   = rdataset_data_mlvj->binnedClone((std::string(rdataset_data_mlvj->GetName())+"_binnedClone").c_str(),(std::string(rdataset_data_mlvj->GetName())+"_binnedClone").c_str());
   TH1F* histo_data = (TH1F*) datahist->createHistogram("histo_data",*rrv_mass_lvj) ;
   histo_data->SetName("histo_data");
   TH1F* histo_func = (TH1F*) model_data->createHistogram("histo_func",*rrv_mass_lvj) ;
   histo_func->SetName("histo_func");

   int Nbin     = int(rrv_mass_lvj->getBins());
   RooArgList rresult_param = rfresult->floatParsFinal();
   int nparameters   = rresult_param.getSize();
   RooAbsReal* ChiSquare = model_data->createChi2(*datahist,RooFit::Extended(kTRUE),RooFit::DataError(RooAbsData::Poisson));
   float chi_over_ndf  = ChiSquare->getVal()/(Nbin-nparameters);

   RooHist* residHist = mplot->residHist("data","model_mc");
   double residual = 0. ;
   for(int iPoint = 0 ; iPoint < residHist->GetN(); iPoint++){
         double x = 0; double y = 0;
         residHist->GetPoint(iPoint,x,y); 
         residual = residual + y*y ;
   }
   mplot->GetYaxis()->SetRangeUser(0,mplot->GetMaximum()*1.2);

   if( FTest_output_txt != ""){
       std::ofstream file_out_FTest ;  
       file_out_FTest.open(FTest_output_txt.c_str(),fstream::app);
       file_out_FTest << " ###################### \n";
       file_out_FTest << " Model pdf "+std::string(model_data->GetName())+" \n";

       TString Line ; Line.Form(" Chi2 Chi2var %f \n",chi_over_ndf);
       file_out_FTest << std::string(Line.Data()); 
       Line.Form(" Residual %f  Nbin %d  nparameters %d \n",residual,Nbin,nparameters);
       file_out_FTest << std::string(Line.Data());
       file_out_FTest.close();
   }

   //Add Chisquare to mplot_pull                                                                                                                                             
   TString Name ; Name.Form("#chi^{2}/ndf = %0.2f ",chi_over_ndf);
   TLatex* cs2 = new TLatex(0.75,0.8,Name.Data());
   cs2->SetNDC();
   cs2->SetTextSize(0.12);
   cs2->AppendPad("same");
   mplot_pull->addObject(cs2);
       
   TString command; command.Form("mkdir -p plots_%s_%s_g1/mlvj_fitting_%s/",channel.c_str(),wtagger.c_str(),mlvj_model.c_str()); 
   system(command.Data());

   command.Form("plots_%s_%s_g1/mlvj_fitting_%s/",channel.c_str(),wtagger.c_str(),mlvj_model.c_str()); 

   Name.Form("m_lvj_sb_lo%s_%s",label.c_str(),mlvj_model.c_str());
   draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),std::string(Name.Data()),mlvj_model,channel,0,logy,GetLumi());


   //#### Decorrelate the parameters in order to have a proper shape in the workspace
   RooWorkspace* wsfit_tmp = new RooWorkspace(("wsfit_tmp+"+label+mlvj_model+"_sb_lo_from_fitting_mlvj").c_str());
   PdfDiagonalizer* Deco   = new PdfDiagonalizer(("Deco"+label+"_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str(),wsfit_tmp,*rfresult);
   std::cout<<"#################### diagonalize data sideband fit "<<std::endl;
   RooAbsPdf* model_pdf_WJets_deco = Deco->diagonalize(*model_pdf_WJets);
   std::cout<<"#################### print parameters "<<std::endl;
   model_pdf_WJets_deco->Print("v");
   model_pdf_WJets_deco->getParameters(rdataset_data_mlvj)->Print("");
   workspace->import(*model_pdf_WJets_deco);

   if(not TString(label).Contains("_jes") and not TString(label).Contains("_jer") and not TString(label).Contains("WJets01") and not TString(label).Contains("WJets1")){

     RooPlot* mplot_sys = rrv_mass_lvj->frame( RooFit::Bins(int(rrv_mass_lvj->getBins())));           

     rdataset_data_mlvj->plotOn(mplot_sys, RooFit::Name("Data"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),RooFit::MarkerColor(0),RooFit::LineColor(0) );
     model_pdf_WJets_deco->plotOn(mplot_sys,RooFit::Name("Nominal"),RooFit::LineColor(kBlack));
     RooRealVar* rrv_number_dataset = new RooRealVar("rrv_number_dataset","rrv_number_dataset",rdataset_data_mlvj->sumEntries());
     rrv_number_dataset->setError(0.); //## only shape uncertainty
                            
     draw_error_band(rdataset_data_mlvj,model_pdf_WJets,rrv_number_dataset,rfresult,mplot,color_palet["Uncertainty"],"F");

     if(workspace->pdf(("model_pdf"+label+"massvbf_jes_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jes_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
        workspace->pdf(("model_pdf"+label+"massvbf_jes_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jes_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_up"), RooFit::LineColor(kBlue));

     if(workspace->pdf(("model_pdf"+label+"massvbf_jes_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jes_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
        workspace->pdf(("model_pdf"+label+"massvbf_jes_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jes_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jes_dn"), RooFit::LineColor(kBlue));

     if(workspace->pdf(("model_pdf"+label+"massvbf_jer_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
        workspace->pdf(("model_pdf"+label+"massvbf_jer_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_up_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_up"), RooFit::LineColor(kBlue));

     if(workspace->pdf(("model_pdf"+label+"massvbf_jer_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
        workspace->pdf(("model_pdf"+label+"massvbf_jer_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_dn_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer_dn"), RooFit::LineColor(kBlue));

     if(workspace->pdf(("model_pdf"+label+"massvbf_jer_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_sb_lo_"+mlvj_model+"from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
        workspace->pdf(("model_pdf"+label+"massvbf_jer_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_mlvj_Deco"+label+"massvbf_jer_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("jer"), RooFit::LineColor(kBlue));

     if(label == "_WJets0" and workspace->pdf(("model_pdf_WJets01_sb_lo"+mlvj_shape["WJets01"]+"_from_fitting_"+channel+"_mlvj_Deco_WJets01_sb_lo"+mlvj_shape["WJets01"]+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str()))
      workspace->pdf(("model_pdf_WJets01_sb_lo"+mlvj_shape["WJets01"]+"_from_fitting_"+channel+"_mlvj_Deco_WJets01_sb_lo"+mlvj_shape["WJets01"]+"_from_fitting_"+channel+"_"+wtagger+"_mlvj").c_str())->plotOn(mplot_sys,RooFit::Name("alt shape"), RooFit::LineColor(kOrange+1));
                           
     TLegend* legend = legend4Plot(mplot_sys,0,-0.07,0.02,0.01,0.,1,channel);
     mplot_sys->addObject(legend);
     command.Form("plots_%s_%s_g1/other/",channel.c_str(),wtagger.c_str());
     Name.Form("mlvj_sb_lo%s_%s_%s",label.c_str(),channel.c_str(),wtagger.c_str());
     draw_canvas_with_pull(mplot_sys,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),std::string(Name.Data()),mlvj_model,channel,0,logy,GetLumi());

   }
   
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


//method to get the alpha function to extrapolate the wjets in the signal region
void get_WJets_mlvj_correction_sb_lo_to_signal_region(RooWorkspace* workspace,const std::string & label, const std::string & mlvj_model, const std::string & mlvj_alternate_model, const std::string & mass_spectrum,const std::string & workspace_label, const std::string & channel, const std::string & wtagger, const int & narrow_factor){

  std::cout<<" ############# get the extrapolation function alpha from MC : "<<label<<"   "<<mlvj_model<<"   "<<mass_spectrum<<" ###############"<<std::endl;
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadTickY(0);
  gStyle->cd();

  //take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here
  RooRealVar* rrv_x = workspace->var("rrv_mass_lvj");
  RooDataSet* rdataset_WJets_sb_lo_mlvj = (RooDataSet*) workspace->data(("rdataset"+workspace_label+label+"_sb_lo_"+channel+mass_spectrum).c_str());
  RooDataSet* rdataset_WJets_signal_region_mlvj = (RooDataSet*) workspace->data(("rdataset"+workspace_label+label+"_signal_region_"+channel+mass_spectrum).c_str());
  rrv_x->Print();
  rdataset_WJets_sb_lo_mlvj->Print();
  rdataset_WJets_signal_region_mlvj->Print();
  
  //### create a frame for the next plots
  RooPlot* mplot = rrv_x->frame(RooFit::Title("correlation_pdf"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor))) ;
  mplot->GetYaxis()->SetTitle("arbitrary units");
  
  RooAbsPdf* correct_factor_pdf = NULL ;

  //Model used for Higgs analysis --> parameters in the SR has to be fitted, not yet done in order to take into account correlations between mj and mlvj
  if(mlvj_model == "ErfExp_v1"){

    RooRealVar* rrv_c_sb       = workspace->var(("rrv_c_ErfExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_offset_sb  = workspace->var(("rrv_offset_ErfExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());            

    RooRealVar* rrv_width_sb   = workspace->var(("rrv_width_ErfExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
            
    RooRealVar* rrv_delta_c      = new RooRealVar(("rrv_delta_c_ErfExp"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_ErfExp"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal() -4*rrv_c_sb->getError(), workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());

    RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfExp"+label+"_"+mlvj_model+"_"+channel+mass_spectrum).c_str(),("rrv_delta_offset_ErfExp"+label+"_"+mlvj_model+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_offset_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),workspace->var(("rrv_offset_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()-4*rrv_offset_sb->getError(), workspace->var(("rrv_offset_ErfExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+rrv_offset_sb->getError());

    RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfExp"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),("rrv_delta_width_ErfExp"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region"+mlvj_model+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region"+mlvj_model+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal() -4*rrv_width_sb->getError(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region"+mlvj_model+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError());

    RooFormulaVar* rrv_c_sr      = new RooFormulaVar(("rrv_c_sr"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),      "@0+@1",RooArgList(*rrv_c_sb,*rrv_delta_c));
    RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
    RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+mlvj_model+"_"+channel+mass_spectrum).c_str(),  "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

    correct_factor_pdf = new RooAlpha("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c_sb,*rrv_offset_sb,*rrv_width_sb,rrv_x->getMin(),rrv_x->getMax());
   
  }

  if(mlvj_model == "ErfPow_v1"){

   RooRealVar* rrv_c_sb      = workspace->var(("rrv_c_ErfPow"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPow"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_width_sb  = workspace->var(("rrv_width_ErfPow"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());

   RooRealVar* rrv_delta_c      = new RooRealVar(("rrv_delta_c_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_c_sb->getError(),100*rrv_c_sb->getError());
   RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_offset_sb->getError(),100*rrv_offset_sb->getError());
   RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_width_sb->getError(),100*rrv_width_sb->getError());

   RooFormulaVar* rrv_c_sr      = new RooFormulaVar(("rrv_c_sr"+label+"_"+channel+mass_spectrum).c_str(),      "@0+@1",RooArgList(*rrv_c_sb,*rrv_delta_c));
   RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
   RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(),  "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

   correct_factor_pdf = new RooAlpha4ErfPowPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "ErfPow2_v1"){

   RooRealVar* rrv_c0_sb     = workspace->var(("rrv_c0_ErfPow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_c1_sb     = workspace->var(("rrv_c1_ErfPow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());

   RooRealVar* rrv_width_sb = NULL;
   if(mass_spectrum == "_mlvj_relaxed") rrv_width_sb   = workspace->var(("rrv_width_ErfPow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
   else{
     rrv_width_sb   = workspace->var(("rrv_width_ErfPow2"+label+"_signal_region"+channel+mass_spectrum).c_str());
     rrv_width_sb->setConstant(kTRUE);
   }
   
   RooRealVar* rrv_delta_c0  = new RooRealVar(("rrv_delta_c0_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal(),workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal() -4*rrv_c0_sb->getError(), workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()+4*rrv_c0_sb->getError());

   RooRealVar* rrv_delta_c1     = new RooRealVar(("rrv_delta_c1_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal(),workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal() -4*rrv_c1_sb->getError(), workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()+4*rrv_c1_sb->getError());

   RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_offset_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal() -4*rrv_offset_sb->getError(), workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+4*rrv_offset_sb->getError());            

   RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_width_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_width_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),workspace->var(("rrv_width_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal() -4*rrv_width_sb->getError(), workspace->var(("rrv_width_ErfPow2"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError());
                
   RooFormulaVar* rrv_c0_sr     = new RooFormulaVar(("rrv_c0_sr"+label+"_"+channel+mass_spectrum).c_str(),"@0+@1",RooArgList(*rrv_c0_sb,*rrv_delta_c0));
   RooFormulaVar* rrv_c1_sr     = new RooFormulaVar(("rrv_c1_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_c1_sb,*rrv_delta_c1));
   
   RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
   RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

   correct_factor_pdf = new RooAlpha4ErfPow2Pdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c0_sr,*rrv_c1_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c0_sb,*rrv_c1_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "ErfPowExp_v1"){ //take initial value from what was already fitted in the SR

     RooRealVar* rrv_c0_sb     = workspace->var(("rrv_c0_ErfPowExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
     RooRealVar* rrv_c1_sb     = workspace->var(("rrv_c1_ErfPowExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
     RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPowExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());

     RooRealVar* rrv_width_sb = NULL ;
     if(mass_spectrum == "_mlvj_relaxed")
        rrv_width_sb      = workspace->var(("rrv_width_ErfPowExp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
     else{
           rrv_width_sb   = workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str());
           rrv_width_sb->setConstant(kTRUE);               
     }

     RooRealVar* rrv_delta_c0  = new RooRealVar(("rrv_delta_c0_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                          workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal(),
			  workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()-4*rrv_c0_sb->getError(),
                          workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()+4*rrv_c0_sb->getError());

     RooRealVar* rrv_delta_c1 = new RooRealVar(("rrv_delta_c1_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                          workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal(),
                          workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()-4*rrv_c1_sb->getError(),
			  workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()+4*rrv_c1_sb->getError());

     RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                          ("rrv_delta_offset_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                          workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),
			  workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()-4*rrv_offset_sb->getError(),workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+4*rrv_offset_sb->getError());

    RooRealVar*  rrv_delta_width = new RooRealVar(("rrv_delta_width_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                         ("rrv_delta_width_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                         workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),
                                         workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()-4*rrv_width_sb->getError(),workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError() );

    RooFormulaVar* rrv_c0_sr     = new RooFormulaVar(("rrv_c0_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",    RooArgList(*rrv_c0_sb,*rrv_delta_c0));
    RooFormulaVar* rrv_c1_sr     = new RooFormulaVar(("rrv_c1_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",    RooArgList(*rrv_c1_sb,*rrv_delta_c1));
    RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
    RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1", RooArgList(*rrv_width_sb,*rrv_delta_width));

    correct_factor_pdf = new RooAlpha4ErfPowExpPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c0_sr,*rrv_c1_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c0_sb,*rrv_c1_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "Exp"){
    RooRealVar* rrv_c_sb    = workspace->var(("rrv_c_Exp"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());    
    rrv_c_sb->Print();
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_Exp"+label+"_"+mlvj_model+"_"+channel+mass_spectrum).c_str(),
                            ("rrv_delta_c_Exp"+label+"_"+mlvj_model+"_"+channel+mass_spectrum).c_str(),
                            workspace->var(("rrv_c_Exp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),
                            workspace->var(("rrv_c_Exp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()-4*rrv_c_sb->getError(),
			    workspace->var(("rrv_c_Exp"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());
    rrv_delta_c->Print();

    correct_factor_pdf = new RooExponential("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c);
    correct_factor_pdf->Print();
  }

  if(mlvj_model == "Pow"){

    RooRealVar* rrv_c_sb    = workspace->var(("rrv_c_Pow"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_Pow"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_Pow"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c_sb->getError(),100*rrv_c_sb->getError());
    correct_factor_pdf = new RooPowPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c);
   
  }

  if(mlvj_model == "ExpN"){
            
    RooRealVar* rrv_c_sb  = workspace->var(("rrv_c_ExpN"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_n_sb  = workspace->var(("rrv_n_ExpN"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_ExpN"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_ExpN"+label+"_"+channel+mass_spectrum).c_str(),
                                  workspace->var(("rrv_c_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),
                                  workspace->var(("rrv_c_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()-4*rrv_c_sb->getError(),
				  workspace->var(("rrv_c_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());

    RooRealVar* rrv_delta_n = new RooRealVar(("rrv_delta_n_ExpN"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_n_ExpN"+label+"_"+channel+mass_spectrum).c_str(),
                                  workspace->var(("rrv_n_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal(),
                                  workspace->var(("rrv_n_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal()-4*rrv_n_sb->getError(),
				  workspace->var(("rrv_n_ExpN"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal()+4*rrv_n_sb->getError());

    correct_factor_pdf = new RooExpNPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c,*rrv_delta_n);

  }

  if (mlvj_model == "ExpTail"){

      RooRealVar* rrv_s_sb = workspace->var(("rrv_s_ExpTail"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
      RooRealVar* rrv_a_sb = workspace->var(("rrv_a_ExpTail"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());

      RooRealVar* rrv_delta_s = new RooRealVar(("rrv_delta_s_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_s_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),
                               workspace->var(("rrv_s_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal(),
                               workspace->var(("rrv_s_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal()-4*rrv_s_sb->getError(),
			       workspace->var(("rrv_s_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal()+4*rrv_s_sb->getError());
      RooRealVar* rrv_delta_a = new RooRealVar(("rrv_delta_a_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_a_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),
                               workspace->var(("rrv_a_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal(),
                               workspace->var(("rrv_a_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal()-4*rrv_a_sb->getError(),
			       workspace->var(("rrv_a_ExpTail"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal()+4*rrv_a_sb->getError());

      RooFormulaVar* rrv_a_sr = new RooFormulaVar(("rrv_a_ExpTail_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_a_sb,*rrv_delta_a ) );
      RooFormulaVar* rrv_s_sr = new RooFormulaVar(("rrv_s_ExpTail_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_s_sb,*rrv_delta_s ) );

      correct_factor_pdf = new RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_s_sr,*rrv_a_sr,*rrv_s_sb,*rrv_a_sb);

  }
        
  if(mlvj_model == "Pow2"){

    RooRealVar* rrv_c0_sb    = workspace->var(("rrv_c0_Pow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_c1_sb    = workspace->var(("rrv_c1_Pow2"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());

    RooRealVar* rrv_delta_c0 = new RooRealVar(("rrv_delta_c0_Pow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_Pow2"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c0_sb->getError(),100*rrv_c0_sb->getError());
    RooRealVar* rrv_delta_c1 = new RooRealVar(("rrv_delta_c1_Pow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_Pow2"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c1_sb->getError(),100*rrv_c1_sb->getError());

    correct_factor_pdf = new RooPow2Pdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c0,*rrv_delta_c1);

  }

  //define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr

  RooCategory* data_category = new RooCategory("data_category","data_category");
  data_category->defineType("sideband");
  data_category->defineType("signal_region");
  data_category->Print();
  RooDataSet* combData4fit = NULL;

  combData4fit = (RooDataSet*) workspace->data(("combData"+workspace_label+label+"_"+channel).c_str());           
  combData4fit->Print();
 
  RooAbsPdf*  model_pdf_sb_lo_WJets    = workspace->pdf(("model_pdf"+label+"_sb_lo"+mlvj_model+"_"+channel+mass_spectrum).c_str());
  RooProdPdf* model_pdf_signal_region_WJets = new RooProdPdf(("model_pdf"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str(),
                                                            ("model_pdf"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str(),*model_pdf_sb_lo_WJets,*correct_factor_pdf);
  model_pdf_sb_lo_WJets->Print();
  model_pdf_signal_region_WJets->Print();
  
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",*data_category);
  simPdf->addPdf(*model_pdf_sb_lo_WJets,"sideband");
  simPdf->addPdf(*model_pdf_signal_region_WJets,"signal_region");
  simPdf->Print();

  RooFitResult* rfresult = simPdf->fitTo(*combData4fit,RooFit::Save(kTRUE), RooFit::SumW2Error(kTRUE));
  rfresult = simPdf->fitTo(*combData4fit,RooFit::Save(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();
  
  // Decorrelate the parameters in the alpha shape
  RooWorkspace* wsfit_tmp = new RooWorkspace(("wsfit_tmp"+label+mlvj_model+"_"+channel+"_sim"+mass_spectrum).c_str());
  std::cout<<"############### diagonalizer alpha "<<std::endl;;
  PdfDiagonalizer* Deco      = new PdfDiagonalizer(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str(),wsfit_tmp,*rfresult);
  RooAbsPdf* correct_factor_pdf_deco = Deco->diagonalize(*correct_factor_pdf);
  correct_factor_pdf_deco->Print();
       
  correct_factor_pdf_deco->getParameters(rdataset_WJets_signal_region_mlvj)->Print("v");
  workspace->import(*correct_factor_pdf_deco);
  
  //in case of default Wjets with default shape
  if(label == "_WJets0"){

    //only mc plots in the SB region
    RooPlot* mplot_sb_lo = rrv_x->frame(RooFit::Title("WJets sb low"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor)));

    //plot just W+Jets sb distribution
    rdataset_WJets_sb_lo_mlvj->plotOn(mplot_sb_lo,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),RooFit::Name("data"));
    model_pdf_sb_lo_WJets->plotOn(mplot_sb_lo ,RooFit::Name("model_mc"));
    RooPlot* mplot_pull_sideband = get_ratio(rrv_x,rdataset_WJets_sb_lo_mlvj,model_pdf_sb_lo_WJets,rfresult,0,1);

    RooArgSet* parameters_list     = model_pdf_sb_lo_WJets->getParameters(rdataset_WJets_sb_lo_mlvj);
    mplot_sb_lo->GetYaxis()->SetRangeUser(0,mplot_sb_lo->GetMaximum()*1.2);

    std::string NameDir  = "plots_"+channel+"_"+wtagger+"_g1";
    std::string NamePlot = "m_lvj"+label+"_sb_lo_sim";

    NameDir  = NameDir+"/other/";
    draw_canvas_with_pull(mplot_sb_lo,mplot_pull_sideband,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,channel,0,1,GetLumi());                

    //only W+jets mc plots in the SR region
    RooPlot* mplot_signal_region = rrv_x->frame(RooFit::Title("WJets sr"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor)));

    rdataset_WJets_signal_region_mlvj->plotOn(mplot_signal_region,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0), RooFit::Name("data"));
    model_pdf_signal_region_WJets->plotOn(mplot_signal_region, RooFit::Name("model_mc"));
    RooPlot* mplot_pull_signal_region = get_pull(rrv_x,mplot_signal_region,rdataset_WJets_signal_region_mlvj,model_pdf_signal_region_WJets,rfresult,"data","model_mc",0,1);


      
    parameters_list = model_pdf_signal_region_WJets->getParameters(rdataset_WJets_signal_region_mlvj);
    mplot_signal_region->GetYaxis()->SetRangeUser(0,mplot_signal_region->GetMaximum()*1.2);

    NamePlot = "m_lvj"+label+"_signal_regionsim";
    draw_canvas_with_pull(mplot_signal_region,mplot_pull_signal_region,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,channel,0,1,GetLumi());
     
    //### Total plot shape in sb_lo, sr and alpha
    model_pdf_sb_lo_WJets->plotOn(mplot,RooFit::Name("Sideband"),RooFit::LineStyle(10));
    model_pdf_signal_region_WJets->plotOn(mplot, RooFit::LineColor(kRed) ,RooFit::LineStyle(8), RooFit::Name("Signal Region"));
    correct_factor_pdf_deco->plotOn(mplot, RooFit::LineColor(kBlack),RooFit::Name("#alpha") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange), RooFit::LineStyle(3),RooFit::Name("#alpha: Alternate PS") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+mlvj_alternate_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+mlvj_alternate_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kMagenta), RooFit::LineStyle(7),RooFit::Name("#alpha: Alternate Function") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_up_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_up_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kRed), RooFit::LineStyle(3),RooFit::Name("#alpha: jes up") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_dn_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_dn_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange+1), RooFit::LineStyle(7),RooFit::Name("#alpha: jes dn") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kAzure), RooFit::LineStyle(3),RooFit::Name("#alpha: jer") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_up_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_up_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kYellow+1), RooFit::LineStyle(7),RooFit::Name("#alpha: jer up") );

    if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_dn_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_dn_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kGray), RooFit::LineStyle(3),RooFit::Name("#alpha: jer dn") );
     

  }

  if(label == "_WJets0" or label == "_WJets1"){ //### draw error band ar 1 and 2 sigma using the decorrelated shape
  
    RooArgList* paras = new RooArgList();
    std::cout<<" originate plots .-------------------------- "<<std::endl; 

    // Make a list of paramters as a function of the model after decorrelation
    if(mlvj_model == "ErfExp_v1" || mlvj_model == "ErfPow_v1" || mlvj_model == "2Exp"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig4").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig5").c_str()));
    }
    else if(mlvj_model == "ErfPow2_v1" || mlvj_model == "ErfPowExp_v1"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig4").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig5").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig6").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig7").c_str()));
    }
    else if(mlvj_model == "Exp" || mlvj_model == "Pow"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
    }
    else if(mlvj_model == "ExpN" || mlvj_model == "ExpTail" || mlvj_model == "Pow2"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
    }
       	      
     draw_error_band_shape_Decor(("correct_factor_pdf_Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str(),"rrv_mass_lvj",paras,workspace,double(2.),mplot,kGreen+2,"F",3001,"#alpha #pm",20,400);
     draw_error_band_shape_Decor(("correct_factor_pdf_Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str(),"rrv_mass_lvj",paras,workspace,double(1.),mplot,kGray+3,"F",3001,"#alpha #pm",20,400);
  }
    
  
  //### plot on the same canvas
  correct_factor_pdf_deco->plotOn(mplot, RooFit::LineColor(kBlack),RooFit::Name("#alpha_invisible"));
     
  if(label == "_WJets0"){ //## add also the plot of alternate ps and function on the canvas
     if(workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kMagenta), RooFit::LineStyle(3),RooFit::Name("#alpha_invisible: Alternate PS") );
     if(workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+mlvj_alternate_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()))
       workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+mlvj_alternate_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange), RooFit::LineStyle(7),RooFit::Name("#alpha_invisible: Alternate Function"));
   }
     
   //### Add the legend
   TLegend* leg = legend4Plot(mplot,-1,-0.07,0.02,0.01,0.,1,channel);
   mplot->addObject(leg);

   //## set the Y axis in arbitrary unit
   float tmp_y_max = 0.28;
   mplot->GetYaxis()->SetRangeUser(0.,tmp_y_max);

   //#### Draw another axis with the real value of alpha
   model_pdf_sb_lo_WJets->getVal(RooArgSet(*rrv_x));
   model_pdf_signal_region_WJets->getVal(RooArgSet(*rrv_x));
   correct_factor_pdf_deco->getVal(RooArgSet(*rrv_x));
   double tmp_alpha_ratio = (model_pdf_signal_region_WJets->getVal(RooArgSet(*rrv_x))/model_pdf_sb_lo_WJets->getVal(RooArgSet(*rrv_x)));
   double tmp_alpha_pdf   = correct_factor_pdf_deco->getVal(RooArgSet(*rrv_x))*mplot->getFitRangeBinW(); // value of the pdf in each point
   double tmp_alpha_scale = tmp_alpha_ratio/tmp_alpha_pdf;

   //#add alpha scale axis
   TGaxis* axis_alpha = new TGaxis( rrv_x->getMax(), 0, rrv_x->getMax(), tmp_y_max, 0, tmp_y_max*tmp_alpha_scale, 510, "+L");
   axis_alpha->SetTitle("#alpha");
   axis_alpha->SetTitleOffset(0.65);
   axis_alpha->SetTitleSize(0.05);
   axis_alpha->SetLabelSize(0.045);
   axis_alpha->SetTitleFont(42);
   axis_alpha->SetLabelFont(42);
   mplot->addObject(axis_alpha);

   std::string NameDir  = "plots_"+channel+"_"+wtagger+"_g1";
   std::string NamePlot = "correction_pdf"+label+"_"+mlvj_model+"_M_lvj_signal_regionto_sideband";

   if(mass_spectrum == "_mlvj_relaxed"){
      NameDir  = NameDir+"/other_relaxed/";
      draw_canvas(mplot,NameDir,NamePlot,channel,GetLumi(),0,1,0);
   }
   else{
      NameDir  = NameDir+"/other/";
      draw_canvas(mplot,NameDir,NamePlot,channel,GetLumi(),0,1,0);
   }
   

   std::cout<<"###################### corrector factor pdf "<<std::endl;      
   correct_factor_pdf_deco->getParameters(rdataset_WJets_sb_lo_mlvj)->Print("v");
      
   std::cout<<"###################### decorrelated sb pdf "<<std::endl;      
   RooAbsPdf* model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco = workspace->pdf(("model_pdf"+label+"_sb_lo"+mlvj_model+"_from_fitting_"+channel+mass_spectrum+"_Deco"+label+"_sb_lo"+mlvj_model+"_from_fitting_"+channel+"_"+wtagger+mass_spectrum).c_str());           
   model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco->Print("v");
    
   //### Wjets shape in the SR correctedfunction * sb
   std::cout<<"###################### Prod pdf "<<std::endl;      
   RooProdPdf* model_pdf_WJets_signal_regionafter_correct_mlvj = new RooProdPdf(("model_pdf"+label+"_signal_region"+mlvj_model+"_"+channel+"_after_correct"+mass_spectrum).c_str(),("model_pdf"+label+"_signal_region"+mlvj_model+"_"+channel+"_after_correct"+mass_spectrum).c_str(),*model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco,*workspace->pdf(("correct_factor_pdf_Deco"+label+"_sim_"+mlvj_model+"_"+channel+"_"+wtagger+mass_spectrum).c_str()));
   model_pdf_WJets_signal_regionafter_correct_mlvj->Print();
   //### fix the parmaters and import in the workspace
   workspace->import(*model_pdf_WJets_signal_regionafter_correct_mlvj);

   //##### calculate the normalization and alpha for limit datacard

   if(!TString(workspace_label.c_str()).Contains("4bias")){
    workspace->var(("rrv_number"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->Print();
    workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->Print();
    workspace->var(("rrv_number"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->setVal(workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->getVal());
    workspace->var(("rrv_number"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->setError(workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->getError());
    workspace->var(("rrv_number"+label+"_signal_region"+mlvj_model+"_"+channel+mass_spectrum).c_str())->setConstant(kTRUE);
   }
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

  RooRealVar* number_TTbar_fail  = new RooRealVar(("rrv_number_TTbar_fail"+label+"_"+
channel).c_str(),("rrv_number_TTbar_fail"+label+"_"+channel).c_str(),rdataset_TTbar_mj_fail->sumEntries());
  RooRealVar* number_STop_fail   = new RooRealVar(("rrv_number_STop_fail"+label+"_"+channel).c_str(),("rrv_number_STop_fail"+label+"_"+channel).c_str(),rdataset_STop_mj_fail->sumEntries());
  RooRealVar* number_VV_fail     = new RooRealVar(("rrv_number_VV_fail"+label+"_"+channel).c_str(),("rrv_number_VV_fail"+label+"_"+channel).c_str(),rdataset_VV_mj_fail->sumEntries());
  RooRealVar* number_WJets_fail  = new RooRealVar(("rrv_number_WJets_fail"+label+"_"+channel).c_str(),("rrv_number_WJets_fail"+label+"_"+channel).c_str(),rdataset_WJets_mj_fail->sumEntries());

  RooAddPdf* model_TTbar_STop_VV_WJets_fail = new RooAddPdf(("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel).c_str(),("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel).c_str(),RooArgList(*model_histpdf_TTbar_fail,*model_histpdf_STop_fail,*model_histpdf_VV_fail,*model_histpdf_WJets_fail), RooArgList(*number_TTbar_fail,*number_STop_fail,*number_VV_fail,*number_WJets_fail));
  workspace->import(*model_TTbar_STop_VV_WJets_fail);
  
  /// inclusive scale factor from total integral                                                                                                                                   
  double scale_number_TTbar_STop_VV_WJets = (rdataset_TTbar_mj->sumEntries()+rdataset_STop_mj->sumEntries()+rdataset_VV_mj->sumEntries()+rdataset_WJets_mj->sumEntries())/(rdataset_data_mj->sumEntries()+rdataset_data_mj_fail->sumEntries());
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

  std::cout<<" model data "<<model_data<<" model data fail "<<model_data_fail<<std::endl;

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
  rfresult_data = simPdf_data->fitTo(*combData_p_f_data,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data),RooFit::Minimizer("Minuit2"), RooFit::SumW2Error(kTRUE));
  rfresult_data = simPdf_data->fitTo(*combData_p_f_data,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data),RooFit::Minimizer("Minuit2"), RooFit::SumW2Error(kTRUE));
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

  rfresult_TotalMC = simPdf_TotalMC->fitTo(*combData_p_f_TotalMC,RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC), RooFit::Minimizer("Minuit2"));
  rfresult_TotalMC = simPdf_TotalMC->fitTo(*combData_p_f_TotalMC,RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC), RooFit::Minimizer("Minuit2"));
  rfresult_TotalMC->Print();
  
         
  RooPlot* xframe_data             = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));
  RooPlot* xframe_data_fail        = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));
  RooPlot* xframe_data_extremefail = rrv_mass_j->frame( RooFit::Bins(int(rrv_mass_j->getBins())));

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
  cut.Form("%s,%s,%s,%s",model_bkg_TotalMC->GetName(),model_STop->GetName(),model_VV->GetName(),model_WJets->GetName());
  simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

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
  cut.Form("%s,%s,%s,%s",model_bkg_TotalMC_fail->GetName(),model_STop_fail->GetName(),model_VV_fail->GetName(),model_WJets_fail->GetName());
  simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

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
  RooRealVar* rrv_mean_gaus_data     = workspace->var(("rrv_mean1_gaus_ttbar_data"+label+"_"+channel).c_str());
  RooRealVar* rrv_sigma_gaus_data    = workspace->var(("rrv_sigma1_gaus_ttbar_data"+label+"_"+channel).c_str());
  RooRealVar* rrv_mean_gaus_TotalMC  = workspace->var(("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+channel).c_str());
  RooRealVar* rrv_sigma_gaus_TotalMC = workspace->var(("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+channel).c_str());

  if(rrv_mean_gaus_TotalMC){
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

    //Sf for LP                                                                                                                                                                        
    std::cout<<"########################## Extreme Fail Analysis ############################## "<<std::endl;
    RooDataSet* rdataset_data_mj_extremefail    = (RooDataSet*) workspace->data(("rdataset_data"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
    RooDataSet* rdataset_TotalMC_mj_extremefail = (RooDataSet*) workspace->data(("rdataset_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
    RooDataSet* rdataset_TTbar_mj_extremefail   = (RooDataSet*) workspace->data(("rdataset_TTbar"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
    RooDataSet* rdataset_STop_mj_extremefail    = (RooDataSet*) workspace->data(("rdataset_STop"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
    RooDataSet* rdataset_VV_mj_extremefail      = (RooDataSet*) workspace->data(("rdataset_VV"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());
    RooDataSet* rdataset_WJets_mj_extremefail   = (RooDataSet*) workspace->data(("rdataset_WJets0"+label+"_"+"extremefailtau2tau1cut_"+channel+"_mj").c_str());

    change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_extremefail);
    change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj_extremefail);
    change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj_extremefail);
    change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj_extremefail);

    RooAbsPdf* model_histpdf_TTbar_extremefail = workspace->pdf((std::string(rdataset_TTbar_mj_extremefail->GetName())+"_histpdf").c_str());
    RooAbsPdf* model_histpdf_STop_extremefail  = workspace->pdf((std::string(rdataset_STop_mj_extremefail->GetName())+"_histpdf").c_str());
    RooAbsPdf* model_histpdf_VV_extremefail    = workspace->pdf((std::string(rdataset_VV_mj_extremefail->GetName())+"_histpdf").c_str());
    RooAbsPdf* model_histpdf_WJets_extremefail = workspace->pdf((std::string(rdataset_WJets_mj_extremefail->GetName())+"_histpdf").c_str());

    RooRealVar* number_TTbar_extremefail = new RooRealVar(("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_TTbar"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_TTbar_mj_extremefail->sumEntries());
    RooRealVar* number_STop_extremefail  = new RooRealVar(("rrv_number_STop"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_STop"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_STop_mj_extremefail->sumEntries());
    RooRealVar*  number_VV_extremefail   = new RooRealVar(("rrv_number_VV"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_VV"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_VV_mj_extremefail->sumEntries());
    RooRealVar* number_WJets_extremefail = new RooRealVar(("rrv_number_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),("rrv_number_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),rdataset_WJets_mj_extremefail->sumEntries());

    RooAddPdf* model_TTbar_STop_VV_WJets_extremefail = new RooAddPdf(("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail"+"_"+channel).c_str(),("model_TTbar_STop_VV_WJets"+label+"_"+"extremefail_"+channel).c_str(),RooArgList(*model_histpdf_TTbar_extremefail,*model_histpdf_STop_extremefail,*model_histpdf_VV_extremefail,*model_histpdf_WJets_extremefail), RooArgList(*number_TTbar_extremefail,*number_STop_extremefail,*number_VV_extremefail,*number_WJets_extremefail) );
    workspace->import(*model_TTbar_STop_VV_WJets_extremefail);

    double scale_number_TTbar_STop_VV_WJets_extremefail = (rdataset_TTbar_mj_extremefail->sumEntries()+rdataset_STop_mj_extremefail->sumEntries()+rdataset_VV_mj_extremefail->sumEntries()+rdataset_WJets_mj_extremefail->sumEntries() )/(rdataset_data_mj_extremefail->sumEntries()) ;

    RooRealVar* rrv_scale_number_TTbar_STop_VV_WJets_extremefail = new RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail").c_str(),("rrv_scale_number_TTbar_STop_VV_WJets"+label+"_"+"extremefail").c_str(),scale_number_TTbar_STop_VV_WJets_extremefail);
    workspace->import(*rrv_scale_number_TTbar_STop_VV_WJets_extremefail);
    
    RooAbsPdf* model_STop_extremefail  = get_STop_mj_Model(workspace,"_STop"+label+"_"+"extremefailtau2tau1cut",mj_shape["STop_extremefail"],channel);
    RooAbsPdf* model_VV_extremefail    = get_VV_mj_Model(workspace,"_VV"+label+"_"+"extremefailtau2tau1cut",mj_shape["VV_extremefail"],channel);
    RooAbsPdf* model_WJets_extremefail = get_General_mj_Model(workspace,"_WJets0"+label+"_"+"extremefailtau2tau1cut",mj_shape["WJets0_extremefail"],channel);


    std::vector<std::string>* constrainslist_data_LP = new std::vector<std::string>();
    RooAbsPdf* model_ttbar_data_extremefail = MakeExtendedModel(workspace,"_ttbar_data"+label+"_"+"extremefailtau2tau1cut",mj_shape["data_extremefail"],"_mj",channel,wtagger,constrainslist_data_LP);
    RooAbsPdf* model_bkg_data_extremefail   = MakeExtendedModel(workspace,"_bkg_data"+label+"_"+"extremefailtau2tau1cut",mj_shape["data_bkg_extremefail"],"_mj",channel,wtagger,constrainslist_data_LP);

    //make the model for the extreme fail fit on data                                                                                                                                 
    RooAddPdf* model_data_extremefail = new RooAddPdf(("model_data"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(),("model_data"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(), RooArgList(*model_ttbar_data_extremefail,*model_bkg_data_extremefail,*model_STop_extremefail,*model_VV_extremefail,*model_WJets_extremefail));
    workspace->import(*model_data_extremefail);

    //constraint the parameter in the extreme fail fit                                                                                                                                
    RooArgSet* pdfconstrainslist_data_LP = new RooArgSet(("pdfconstrainslist_data_LP"+label+"_"+channel).c_str());
    if(constrainslist_data_LP!=NULL){
     for(unsigned int i = 0; i < constrainslist_data_LP->size(); i++){
       workspace->pdf(constrainslist_data_LP->at(i).c_str())->Print();
       pdfconstrainslist_data_LP->add(*workspace->pdf(constrainslist_data_LP->at(i).c_str()));
     }
    }
   
    RooFitResult* rfresult_data_LP = model_data_extremefail->fitTo(*rdataset_data_mj_extremefail,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data_LP),RooFit::Minimizer("Minuit2"));
    rfresult_data_LP = model_data_extremefail->fitTo(*rdataset_data_mj_extremefail,RooFit::Save(kTRUE),RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_data_LP),RooFit::Minimizer("Minuit2"));
    rfresult_data_LP->Print();

    cut.Form("rrv_number_bkg_data%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
    RooRealVar* rrv_number_bkg_data_extremefailtau2tau1cut   = workspace->var(cut.Data());
    cut.Form("rrv_number_ttbar_data%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
    RooRealVar* rrv_number_ttbar_data_extremefailtau2tau1cut = workspace->var(cut.Data());
    rrv_number_ttbar_data_extremefailtau2tau1cut->setError((rrv_number_ttbar_data_extremefailtau2tau1cut->getVal()+rrv_number_bkg_data_extremefailtau2tau1cut->getVal())/2. );
    rrv_number_bkg_data_extremefailtau2tau1cut->setError((rrv_number_ttbar_data_extremefailtau2tau1cut->getVal()+rrv_number_bkg_data_extremefailtau2tau1cut->getVal())/2. );

    model_STop_extremefail->Print();
    model_VV_extremefail->Print();
    model_WJets_extremefail->Print();
    rdataset_data_mj_extremefail->Print();

    //extremefail fit of the mc                                                                                                                                                       
    std::vector<std::string>* constrainslist_TotalMC_LP = new std::vector<std::string>();
    RooAbsPdf*  model_ttbar_TotalMC_extremefail = MakeExtendedModel(workspace,"_ttbar_TotalMC"+label+"_"+"extremefailtau2tau1cut",mj_shape["mc_extremefail"],"_mj",channel,wtagger,constrainslist_TotalMC_LP);
    RooAbsPdf*  model_bkg_TotalMC_extremefail = MakeExtendedModel(workspace,"_bkg_TotalMC"+label+"_"+"extremefailtau2tau1cut",mj_shape["mc_bkg_extremefail"],"_mj",channel,wtagger,constrainslist_TotalMC_LP);

    RooAddPdf* model_TotalMC_extremefail = new RooAddPdf(("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(),("model_TotalMC"+label+"_"+"extremefailtau2tau1cut_"+channel).c_str(), RooArgList(*model_ttbar_TotalMC_extremefail,*model_bkg_TotalMC_extremefail,*model_STop_extremefail,*model_VV_extremefail,*model_WJets_extremefail));
    workspace->import(*model_TotalMC_extremefail);

    RooArgSet* pdfconstrainslist_TotalMC_LP = new RooArgSet(("pdfconstrainslist_TotalMC_LP_"+channel).c_str());
    if(constrainslist_TotalMC_LP!=NULL){
     for(unsigned int i =0; i < constrainslist_TotalMC_LP->size(); i++){
      workspace->pdf(constrainslist_TotalMC_LP->at(i).c_str())->Print();
      pdfconstrainslist_TotalMC_LP->add(*workspace->pdf(constrainslist_TotalMC_LP->at(i).c_str()));
     }
    }

    RooFitResult* rfresult_TotalMC_LP = model_TotalMC_extremefail->fitTo(*rdataset_TotalMC_mj_extremefail, RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC_LP),RooFit::Minimizer("Minuit2"));
    rfresult_TotalMC_LP = model_TotalMC_extremefail->fitTo(*rdataset_TotalMC_mj_extremefail, RooFit::Save(kTRUE), RooFit::Range("controlsample_fitting_range"),RooFit::ExternalConstraints(*pdfconstrainslist_TotalMC_LP),RooFit::Minimizer("Minuit2"));
    rfresult_TotalMC_LP->Print();
    cut.Form("rrv_number_bkg_TotalMC%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
    RooRealVar* rrv_number_bkg_TotalMC_extremefailtau2tau1cut   = workspace->var(cut.Data());
    cut.Form("rrv_number_ttbar_TotalMC%s_extremefailtau2tau1cut_%s_mj",label.c_str(),channel.c_str());
    RooRealVar* rrv_number_ttbar_TotalMC_extremefailtau2tau1cut = workspace->var(cut.Data());
    rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->setError((rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut->getVal())/2.);
    rrv_number_bkg_TotalMC_extremefailtau2tau1cut->setError( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut->getVal()+rrv_number_bkg_TotalMC_extremefailtau2tau1cut->getVal())/2. );
    model_STop_extremefail->Print();
    model_VV_extremefail->Print();
    model_WJets_extremefail->Print();
    rdataset_TotalMC_mj_extremefail->Print();


    rdataset_data_mj_extremefail->plotOn(xframe_data_extremefail, RooFit::Name("data"), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));
    rdataset_TotalMC_mj_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("TotalMC_invisible"), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0));

    model_TotalMC_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("total MC fit"),RooFit::NormRange("controlsample_fitting_range"), RooFit::LineStyle(kDashed));
    cut.Form("%s,%s,%s,%s",model_bkg_TotalMC_extremefail->GetName(),model_STop_extremefail->GetName(),model_VV_extremefail->GetName(),model_WJets_extremefail->GetName());
    model_TotalMC_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("MC fit bkg"),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));


    rdataset_data_mj_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data_invisible"), RooFit::MarkerColor(kWhite), RooFit::LineColor(kWhite), RooFit::Invisible(), RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0) );

    model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("total data fit"),RooFit::NormRange("controlsample_fitting_range"));
    model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data_fit_invisible"),RooFit::NormRange("controlsample_fitting_range"));
    cut.Form("%s,%s,%s,%s",model_bkg_data_extremefail->GetName(),model_STop_extremefail->GetName(),model_VV_extremefail->GetName(),model_WJets_extremefail->GetName());
    model_data_extremefail->plotOn(xframe_data_extremefail,RooFit::Name("data fit bkg"),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed));

    xframe_data_extremefail->GetYaxis()->SetRangeUser(0.,xframe_data_extremefail->GetMaximum()*1.3);
    TLegend* leg_data_extremefail = legend4Plot(xframe_data_extremefail,0,-0.2,0.07,0.04,0.,1,channel);                                           
    xframe_data_extremefail->addObject(leg_data_extremefail);

    nameDir.Form("plots_%s_%s/m_j_fitting_TTbar_controlsample_wtaggercut%s_nearbyJets_v3/",channel.c_str(),wtagger.c_str(),wtagger.c_str());
    namePlot.Form("control%s_%s_%s_extremefail_pTbin_%d_%d",label.c_str(),wtagger.c_str(),channel.c_str(),int(ca8_ungroomed_pt_min),int(ca8_ungroomed_pt_max));
    draw_canvas(xframe_data_extremefail,std::string(nameDir),std::string(namePlot),channel,GetLumi(),0,1,0);    
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
    cut.Form("model_bkg_TotalMC_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    simPdf_TotalMC->plotOn(xframe_data,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    // plot data fit function
    cut.Form("category_p_f%s_%s==category_p_f%s_%s::pass",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());                                                    
    combData_p_f_data->plotOn(xframe_data,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

    simPdf_data->plotOn(xframe_data,RooFit::Name("total data fit"),RooFit::Slice(*category_p_f,"pass"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
    cut.Form("model_bkg_data_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
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

    cut.Form("model_bkg_TotalMC_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    simPdf_TotalMC->plotOn(xframe_data_fail,RooFit::Name("MC fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed));

    //fail plots -> plot data fit                                                                                                                                                       
    cut.Form("category_p_f%s_%s==category_p_f%s_%s::fail",label.c_str(),channel.c_str(),label.c_str(),channel.c_str());
    combData_p_f_data->plotOn(xframe_data_fail,RooFit::Name("data_invisible"),RooFit::Cut(cut.Data()),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

    simPdf_data->plotOn(xframe_data_fail,RooFit::Name("total data fit"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
    simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data_fit_invisible"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"));
    cut.Form("model_bkg_data_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    simPdf_data->plotOn(xframe_data_fail,RooFit::Name("data fit bkg"),RooFit::Slice(*category_p_f,"fail"),RooFit::ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),RooFit::NormRange("controlsample_fitting_range"), RooFit::Components(cut.Data()),RooFit::LineColor(kRed));
      
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

///////////////////////////

/*** exponential for I126***/
double exponential(double* x, double* par)
{
  //[0] = N   normalization factor                                                                                                                              
  //[1] = lambda   slope                                                                                                                                        
  double xx = x[0];
  double N = par[0];
  double lambda = par[1];

  //std::cout << "N: " << N << "   lambda: " << lambda << std::endl;                                                                                            
  return N * exp(-1.*lambda*xx);
}


////////////////////////
//---crystal ball: for the signal -----
double crystalBallLowHigh (double* x, double* par) {
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha on the right-hand side
  //[4] = n
  //[5] = alpha2 on the left-hand side
  //[6] = n2

 double xx = x[0];
 double mean = par[1];
 double sigma = fabs (par[2]);
 double alpha = par[3];
 double n = par[4];
 double alpha2 = par[5];
 double n2 = par[6];

 if( (xx-mean)/sigma > fabs(alpha) ) {
  double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
  double B = n/fabs(alpha) - fabs(alpha);

  return par[0] * A * pow(B + (xx-mean)/sigma, -1.*n);
 }

 else if( (xx-mean)/sigma < -1.*fabs(alpha2) ) {
  double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
  double B = n2/fabs(alpha2) - fabs(alpha2);

  return par[0] * A * pow(B - (xx-mean)/sigma, -1.*n2);
 }

 else {
  return par[0] * exp(-1. * (xx-mean)*(xx-mean) / (2*sigma*sigma) );
 }

}

// ---- crystallball+Exp: for S+I ----

double doubleGausCrystalBallLowHighPlusExp (double* x, double* par) {
  //[0] = N
  //[1] = mean
  //[2] = sigma
  //[3] = alpha
  //[4] = n
  //[5] = alpha2
  //[6] = n2

  //[7] = R = ratio between exponential and CB
  //[8] = tau = tau falling of exponential

 double xx = x[0];

// double mean = par[1] ; // mean
// double sigmaP = par[2] ; // sigma of the positive side of the gaussian
// double sigmaN = par[3] ; // sigma of the negative side of the gaussian
// double alpha = par[4] ; // junction point on the positive side of the gaussian
// double n = par[5] ; // power of the power law on the positive side of the gaussian
// double alpha2 = par[6] ; // junction point on the negative side of the gaussian
// double n2 = par[7] ; // power of the power law on the negative side of the gaussian

 double mean = par[1] ; // mean
  double sigmaP = par[2] ; // sigma of the positive side of the gaussian | they are the same!!!
  double sigmaN = par[2] ; // sigma of the negative side of the gaussian |
 double alpha = par[3] ; // junction point on the positive side of the gaussian
 double n = par[4] ; // power of the power law on the positive side of the gaussian
  double alpha2 = par[5] ; // junction point on the negative side of the gaussian
  double n2 = par[6] ; // power of the power law on the negative side of the gaussian

 double R = par[7] ;
 double tau = par[8] ;

 
 if ((xx-mean)/sigmaP > fabs(alpha)) {
  double A = pow(n/fabs(alpha), n) * exp(-0.5 * alpha*alpha);
  double B = n/fabs(alpha) - fabs(alpha);

  return par[0] * ( (1-R)*(A * pow(B + (xx-mean)/sigmaP, -1.*n)) + R * exp(-xx/tau));
 }

  else if ((xx-mean)/sigmaN < -1.*fabs(alpha2)) {
  double A = pow(n2/fabs(alpha2), n2) * exp(-0.5 * alpha2*alpha2);
  double B = n2/fabs(alpha2) - fabs(alpha2);

  return par[0] * ( (1-R)*(A * pow(B - (xx-mean)/sigmaN, -1.*n2)) + R * exp(-xx/tau));
  }

 else if ((xx-mean) > 0) {
   return par[0] * ( (1-R)*exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaP*sigmaP) ) + R * exp(-xx/tau));
  }

 else {
   return par[0] * ( (1-R)*exp(-1. * (xx-mean)*(xx-mean) / (2*sigmaN*sigmaN) ) + R * exp(-xx/tau));
 }

}

//---- division of CBLowHighPlusExp with CBLowHigh ----
Double_t CrystalBallLowHighPlusExpDividedByCrystalBallLowHigh(Double_t *x,Double_t *par) {

 double den = crystalBallLowHigh (x, par + 9) ; // SM signal 
 if (den == 0) return -1. ;
 double num = doubleGausCrystalBallLowHighPlusExp (x, par) ; // signal + I800 +I126  
 double I126 = exponential (x,par+16); //interference from H126

 double alpha = par[18];
 double beta = par[19];
 double zeta = par[20];

 double S_par[7]={par[9],par[10],beta*par[11],par[12],par[13],par[14],par[15]};
 double S = crystalBallLowHigh (x, S_par) ; //SM signal with width rescaled for c',BRnew

 double I = (num - beta*S- sqrt(1-beta)*I126)/sqrt(beta); //I of the heavy-Higgs

 double w = (alpha*S + sqrt(alpha)*I + sqrt(zeta)*I126) / den; //weight

  if (w<0)   return 0;
  else       return w;

}


Double_t getIntWght(std::string wFile , double realMass, double Hmass , double cprime, double BRnew ) {

    TString *readfile = new TString (wFile); //file with the values of the all parameters
    TFile* SI = new TFile(readfile->Data());
    Double_t fill_param[18]; // 9 + 7 + 2 = 18

    TGraph2D* variables_S[7];
    TGraph2D* variables_SI[9];
    TGraph2D* variables_I126[2];
    TF1* crystal_Icorr_qqH;

    double alpha = cprime * (1-BRnew);
    double beta = cprime / (1-BRnew);
    double zeta = 1-cprime;

  TString parameters_normal [9] = {"Norm","Mean_CB","Sigma_CB","alphaR_CB","nR_CB","alphaL_CB","nL_CB","R","Tau"};
  TString parameters_I126 [2] = {"N_exp","Tau_exp"};

    for (int i=0; i<9; i++) {
     TString *name = new TString (parameters_normal[i]);
     name->Append("_SI.txt");
     variables_SI[i] = (TGraph2D*)SI->Get(name->Data());
    }
    for (int i=0; i<7; i++) {
     TString *name = new TString (parameters_normal[i]);
     name->Append("_S.txt");
     variables_S[i] = (TGraph2D*)SI->Get(name->Data());
    }
    for (int i=0; i<2; i++) {
      TString *name = new TString (parameters_I126[i]);
      name->Append("_I126.txt");
      variables_I126[i] = (TGraph2D*)SI->Get(name->Data());
    }



    crystal_Icorr_qqH = new TF1("crystal_Icorr_qqH",CrystalBallLowHighPlusExpDividedByCrystalBallLowHigh,0,3000,16+2+3);

    for (int iVar = 0; iVar<9; iVar++) {
      if (parameters_normal[iVar].Contains("Norm")) {
           crystal_Icorr_qqH->SetParameter(iVar, exp(variables_SI[iVar]->Interpolate(realMass, beta)));
      }
      else {
           crystal_Icorr_qqH->SetParameter(iVar, variables_SI[iVar]->Interpolate(realMass, beta));
      }
    }

    for (int iVar = 0; iVar<7; iVar++) {
      if (parameters_normal[iVar].Contains("Norm")) {
     	   crystal_Icorr_qqH->SetParameter(iVar+9, exp(variables_S[iVar]->Interpolate(realMass, beta)));
      }
      else {
     	   crystal_Icorr_qqH->SetParameter(iVar+9, variables_S[iVar]->Interpolate(realMass, beta));
      }
    }

    for (int iVar = 0; iVar<2; iVar++) {
     	   crystal_Icorr_qqH->SetParameter(iVar+16, variables_I126[iVar]->Interpolate(realMass, beta));
    }

    crystal_Icorr_qqH->SetParameter(6+9+3, alpha);
    crystal_Icorr_qqH->SetParameter(6+9+4, beta);
    crystal_Icorr_qqH->SetParameter(6+9+5, zeta);

    Double_t wInt;
    wInt = 1.;
    wInt = crystal_Icorr_qqH->Eval(Hmass);
    delete SI;

   return wInt;

}
