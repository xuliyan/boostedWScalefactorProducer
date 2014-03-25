#include "FitUtils.h"

void fit_mj_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & model,
const std::string & channel,const std::string & wtagger_label,const int & additionalinformation){

     std::cout<<"############### Fit mj single MC sample"<<fileName<<" "<<label<<"  "<<model<<" ##################"<<std::endl;
     //import variable and dataset
     RooRealVar* rrv_mass_j = workspace->var("rrv_mass_j");
     RooDataSet* rdataset_mj = (RooDataSet*) workspace->data(("rdataset4bias"+label+"_"+channel+"_mj").c_str());
     rdataset_mj->Print();
       
     //## make the extended model
     RooArgList* constraint_list = new RooArgList(); 
     RooExtendPdf* model_pdf = NULL ; 
     if(additionalinformation == 1)
         model_pdf = MakeExtendedModel(workspace,label+model,model,"_mj",channel,wtagger_label,constraint_list);
     else
         model_pdf = MakeExtendedModel(workspace,label,model,"_mj",channel,wtagger_label,constraint_list);

     RooFitResult* rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::Extended(kTRUE),RooFit::SumW2Error(kTRUE));
     rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2"));
     rfresult = model_pdf->fitTo(*rdataset_mj,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2"));
     rfresult->Print();

     //## Plot the result
     RooPlot* mplot = rrv_mass_j->frame(RooFit::Title((label+" fitted by "+model).c_str()), RooFit::Bins(int(rrv_mass_j->getBins())));
     rdataset_mj->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));

     //## draw the error band for an extend pdf
     draw_error_band_extendPdf(rdataset_mj,model_pdf,rfresult,mplot,2,"L");

     //## re-draw the dataset
     rdataset_mj->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0),RooFit::Name("data"));
     //## draw the function
     model_pdf->plotOn(mplot,RooFit::Name("model_mc")); // remove RooFit.VLines() in order to get right pull in the 1st bin

     //## Get the pull
     RooPlot* mplot_pull = get_pull(rrv_mass_j,mplot,rdataset_mj,model_pdf,rfresult,"data","model_mc",0,1);
     mplot->GetYaxis()->SetRangeUser(1e-5,mplot->GetMaximum()*1.2);

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

      RooArgSet* parameters_list = model_pdf->getParameters(*rdataset_mj);

      TString command; 
      command.Form("mkdir -p plots_%s_%s_g1/mj_fitting/",channel.c_str(),wtagger_label.c_str());
      system(command.Data());

      command.Form("plots_%s_%s_g1/mj_fitting/",channel.c_str(),wtagger_label.c_str());
      draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),label+fileName,model,"em",0,1,GetLumi());

      workspace->import(*model_pdf);
      workspace->import(*rfresult);

      std::cout<<"####################################################"<<std::endl;
      std::cout<<"######## Normalization Factor in mJ ################"<<std::endl;
      std::cout<<"####################################################"<<std::endl;

      //normalize the number of total events to lumi --> correct the number to scale to the lumi
      if(additionalinformation == 1){
         workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setVal(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());
         workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());

         workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->Print();
         
         if(TString(label).Contains("ggH")){
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setVal( workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal());
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError());
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->Print();
         }
         if(TString(label).Contains("vbfH")){
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setVal( workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal());
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError());
            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->Print();
         }
      }
      else{
         workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setVal(workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());
         workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel).c_str())->getVal());

         workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->Print();
         
         if(TString(label).Contains("ggH")){
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setVal( workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getVal());
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getError());
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->Print();
         }
         if(TString(label).Contains("vbfH")){
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setVal( workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getVal());
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->setError(workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->getError());
            workspace->var(("rrv_number"+label+"_"+channel+"_mj").c_str())->Print();
        } 
      }
      //##### apply the correction of the mean and sigma from the ttbar control sample to the STop, TTbar and VV
      TIter par = parameters_list->createIterator();
      par.Reset();

      SystematicUncertaintyHiggs higgsUncertainty;

      RooRealVar* param = dynamic_cast<RooRealVar*>(par.Next());
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


//## Define the Extended Pdf for and mlvj fit giving: label, fit model name, list constraint, range to be fitted and do the decorrelation
void fit_mlvj_model_single_MC(RooWorkspace* workspace, const std::string & fileName, const std::string & label, const std::string & in_range, const std::string & mlvj_model, const std::string & channel, const std::string & wtagger_label, const int & deco, const int & show_constant_parameter, const int & logy, const int & ismc){

        std::cout<<"############### Fit mlvj single MC sample "<<fileName<<" "<<label<<" "<<mlvj_model<<" "<<in_range<<" ##################"<<std::endl;
        //## imporparam_generatedt variable and dataset
        RooRealVar* rrv_mass_lvj = workspace->var("rrv_mass_lvj");
        RooDataSet* rdataset     = (RooDataSet*) workspace->data(("rdataset4bias"+label+in_range+"_"+channel+"_mlvj").c_str());
        RooArgList* constraintlist = new RooArgList();

        //## make the extended pdf model
        RooExtendPdf* model = MakeExtendedModel(workspace,label+in_range+mlvj_model,mlvj_model,"_mlvj",channel,wtagger_label,constraintlist,ismc);
        //make the fit
        RooFitResult* rfresult = model->fitTo(*rdataset,RooFit::Save(1),RooFit::SumW2Error(kTRUE),RooFit::Extended(kTRUE));
        rfresult = model->fitTo(*rdataset,RooFit::Save(1), RooFit::SumW2Error(kTRUE) ,RooFit::Extended(kTRUE), RooFit::Minimizer("Minuit2") );
        rfresult->Print();

        //## set the name of the result of the fit and put it in the workspace
        rfresult->SetName(("rfresult"+label+in_range+"_"+channel+"_mlvj").c_str());

        //## plot the result
        RooPlot* mplot = rrv_mass_lvj->frame(RooFit::Title(("M_{lvj"+in_range+"} fitted by "+mlvj_model).c_str()), RooFit::Bins(int(rrv_mass_lvj->getBins())));
        rdataset->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0));
        //## plot the error band but don't store the canvas (only plotted without -b option
        draw_error_band_extendPdf(rdataset,model,rfresult,mplot,2,"L");
        model->plotOn(mplot,RooFit::Name("model_mc"));//#, RooFit.VLines()); in order to have the right pull
        rdataset->plotOn(mplot,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Name("data"));

        //## get the pull
        RooPlot* mplot_pull = get_pull(rrv_mass_lvj,mplot,rdataset,model,rfresult,"data","model_mc",0,1);
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
        
        mplot->GetYaxis()->SetRangeUser(1e-2,mplot->GetMaximum()*1.2);
	/*
        if options.shapetest == 1:
         file_out_FTest.write(" ###################### \n");
         file_out_FTest.write(" Model pdf %s"%(model.GetName()));
         file_out_FTest.write(" Chi2 Chi2var %0.2f "%(chi_over_ndf));
         file_out_FTest.write(" Residual %0.2f   Nbin %0.2f nparameters %0.2f \n"%(residual,Nbin,nparameters));
	*/

        //##Add Chisquare to mplot_pull   
        TString command; command.Form("#chi^{2}/ndf = %0.2f ",float(chi_over_ndf));                                                                                                               TLatex* cs2 = new TLatex(0.75,0.8,command.Data());
        cs2->SetNDC();
        cs2->SetTextSize(0.12);
        cs2->AppendPad("same");
        mplot_pull->addObject(cs2);

      command.Form("mkdir -p plots_%s_%s_g1/mlvj_fitting/",channel.c_str(),wtagger_label.c_str());
      system(command.Data());

      command.Form("plots_%s_%s_g1/mlvj_fitting/",channel.c_str(),wtagger_label.c_str()); 
      draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),label+fileName,mlvj_model,"em",0,logy,GetLumi());

      //## Number of the event in the dataset and lumi scale factor --> set the proper number for bkg extraction or for signal region
      std::cout<<"rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj"<<std::endl;
      std::cout<<"rrv_scale_to_lumi"+label+"_"+channel+in_range+"_mlvj"<<std::endl;;
        
      workspace->import(*model);
      workspace->import(*rfresult);

      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setVal(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel+in_range+"_mlvj").c_str())->getVal());
      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->setError(workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*workspace->var(("rrv_scale_to_lumi"+label+"_"+channel+in_range+"_mlvj").c_str())->getVal() );
      workspace->var(("rrv_number"+label+in_range+mlvj_model+"_"+channel+"_mlvj").c_str())->Print();

}

//#### make the mj sideband fit on data ti get the Wjets normaliztion
void fit_WJetsNormalization_in_Mj_signal_region(RooWorkspace* workspace,  std::map<std::string,int> color_palet, const std::string & label, const std::string & model, const std::string & channel, const std::string wtagger, const int & ttbarcontrolregion, const int & pseudodata, const float & mj_signal_min, const float & mj_signal_max){

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

   if(ttbarcontrolregion){
         model_WJets  = get_WJets_mj_Model(workspace,"_WJets0","",channel);
         model_STop   = get_STop_mj_Model(workspace,"_STop","",channel);
         model_VV     = get_VV_mj_Model(workspace,"_VV","",channel);
         model_WW_EWK = get_WW_EWK_mj_Model(workspace,"_WW_EWK","",channel);
	 model_TTbar  = get_TTbar_mj_Model(workspace,label,model,channel,0);
   }
   else{
         model_TTbar  = get_TTbar_mj_Model(workspace,"_TTbar","",channel);
         model_STop   = get_STop_mj_Model(workspace,"_STop","",channel);
         model_VV     = get_VV_mj_Model(workspace,"_VV","",channel);
         model_WW_EWK = get_WW_EWK_mj_Model(workspace,"_WW_EWK","",channel);
	 model_WJets  = get_WJets_mj_Model(workspace,label,model,channel,0);
  }

  //## Total Pdf and fit only in sideband
  RooAddPdf* model_data = new RooAddPdf(("model_data_"+channel+"_mj").c_str(),("model_data_"+channel+"_mj").c_str(),RooArgList(*model_WJets,*model_VV,*model_WW_EWK,*model_TTbar,*model_STop));
                                    
  RooFitResult* rfresult = model_data->fitTo(*rdataset_data_mj,RooFit::Save(1),RooFit::Range("sb_lo,sb_hi") ,RooFit::Extended(kTRUE),RooFit::NumCPU(4));
  rfresult = model_data->fitTo(*rdataset_data_mj,RooFit::Save(1),RooFit::Range("sb_lo,sb_hi"),RooFit::Extended(kTRUE), RooFit::NumCPU(4),RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();
  workspace->import(*model_data);

  //## Total numver of event --> full propagation of error due to all the background sources coming from the fit
  RooRealVar* rrv_number_data_mj = NULL ;

  if (ttbarcontrolregion){
         
   rrv_number_data_mj = new RooRealVar(("rrv_number_data_+"+channel+"_mj").c_str(),("rrv_number_data_"+channel+"_mj").c_str(),
                                   workspace->var(("rrv_number+"+label+model+"_"+channel+" _mj").c_str())->getVal()+ //##TTbar
                                   workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getVal()+  //##STop
                                   workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getVal()+    //## VV
                                   workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getVal()+ //## WW_EWK
                                   workspace->var(("rrv_number_WJets0_"+channel+"_mj").c_str())->getVal());  //## WJets

   rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getError()+
                                            workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getError()+       
                                            workspace->var(("rrv_number_WJets0_"+channel+"_mj").c_str())->getError()*
                                            workspace->var(("rrv_number_WJets0_"+channel+"_mj").c_str())->getError()));         

  workspace->import(*rrv_number_data_mj);

  std::cout<<"TTbar  events: "<<workspace->var(("rrv_number_"+label+model+"_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"STop   events: "<<workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"VV     events: "<<workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"WW_EWK events: "<<workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"WJets  events: "<<workspace->var(("rrv_number_WJets0_"+channel+"_mj").c_str())->getVal()<<std::endl;
  std::cout<<"Data   events: "<<workspace->var(("rrv_number_data_"+channel+"_mj").c_str())->getVal()<<std::endl;

  }
  else{
         
    rrv_number_data_mj = new RooRealVar(("rrv_number_data_"+channel+"_mj").c_str(),("rrv_number_data_"+channel+"_mj").c_str(),
                                        workspace->var(("rrv_number_TTbar_"+channel+"_mj").c_str())->getVal()+ //## TTbar
                                        workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getVal()+  //## STop
                                        workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getVal()+    //## VV
                                        workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getVal()+ //## WW_EWK
                                        workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal());  //## WJets

    rrv_number_data_mj->setError(TMath::Sqrt(workspace->var(("rrv_number_TTbar_"+channel+"_mj").c_str())->getError()*
                                             workspace->var(("rrv_number_TTbar_"+channel+"_mj").c_str())->getError()+
                                             workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getError()*
                                             workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getError()+
                                             workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getError()*
                                             workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getError()+
                                             workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getError()*
                                             workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getError()+       
                                             workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()*
					     workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getError()));         
   workspace->import(*rrv_number_data_mj);

   std::cout<< "TTbar  events: ",workspace->var(("rrv_number_TTbar_"+channel+"_mj").c_str())->getVal();
   std::cout<< "STop   events: ",workspace->var(("rrv_number_STop_"+channel+"_mj").c_str())->getVal();
   std::cout<< "VV     events: ",workspace->var(("rrv_number_VV_"+channel+"_mj").c_str())->getVal();
   std::cout<< "WW_EWK events: ",workspace->var(("rrv_number_WW_EWK_"+channel+"_mj").c_str())->getVal();
   std::cout<< "WJets  events: ",workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal();
   std::cout<< "Data   events: ",workspace->var(("rrv_number_data_"+channel+"_mj").c_str())->getVal();
   }

  //## draw the plot for the default WJets Shape
  RooPlot* mplot = rrv_mass_j->frame(RooFit::Title("rrv_mass_j"), RooFit::Bins(int(rrv_mass_j->getBins())));
  rdataset_data_mj->plotOn(mplot,RooFit::Name("data_invisible"), RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::Poisson), RooFit::XErrorSize(0),RooFit::Invisible());
  //## plot solid style
  if (ttbarcontrolregion){

    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("TTbar"),RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["_TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
   
    Name.Form("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WW_EWK"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_VV_%s_mj,model_STop_%s_mj",channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj",channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj",channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    // plot "dashed" style area
    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    //plot "dashed" style area
    Name.Form("model_STop_%s_mj,model_WJets0_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WW_EWK_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
           
    //plot "dashed" style area
    Name.Form("model_WJets0_%s_mj,model_VV_%s_mj,model_STop_%s_mj",channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    // plot "dashed" style area
    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj",channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("STop_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
    //plot "dashed" style area
    Name.Form("model_WJets0_%s_mj",channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
    /// solid line
    Name.Form("model_WJets0_%s_mj",channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj",channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj",channel.c_str(),channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model_WJets0_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj,model%s_%s_mj",channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),(label+model).c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
  }

  else{

    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WW_EWK"),RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
   
    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("VV"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("TTbar"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("STop"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WJets"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    // plot "dashed" style area
    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WW_EWK_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    //plot "dashed" style area
    Name.Form("model%s_%s_mj,model_STop_%s_mj,model_TTbar_%s_mj,model_VV_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("VV_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
           
    //plot "dashed" style area
    Name.Form("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    // plot "dashed" style area
    Name.Form("model%s_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("TTbar_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
    //plot "dashed" style area
    Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());

    model_data->plotOn(mplot,RooFit::Name("WJets_invisible"), RooFit::Components(Name.Data()),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack),RooFit::FillStyle(3003),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
	 
    /// solid line
    Name.Form("model%s_%s_mj",(label+model).c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2) ,RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj,model_VV_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());

    Name.Form("model%s_%s_mj,model_TTbar_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WW_EWK_%s_mj",(label+model).c_str(),channel.c_str(),channel.c_str(),channel.c_str(),channel.c_str(),(label+model).c_str());
    model_data->plotOn( mplot,RooFit::Name("_invisible"), RooFit::Components(Name.Data()), RooFit::LineColor(kBlack), RooFit::LineWidth(2),RooFit::NormRange("sb_lo,sb_hi"), RooFit::VLines());
                                                                     

  }
        
  /// draw the error band using the sum of all the entries component MC + fit and the total error == Normalization for the fixed MC, shape + normalization for W+jets
  draw_error_band(rdataset_data_mj, model_data, rrv_number_data_mj,rfresult,mplot,color_palet["Uncertainty"],"F");
  model_data->plotOn(mplot,RooFit::Name("model_mc"),RooFit::Range(rrv_mass_j->getMin(),rrv_mass_j->getMax()),RooFit::NormRange("sb_lo,sb_hi"),RooFit::Invisible());
                        
  if (pseudodata == 1) rdataset_data_mj->plotOn( mplot ,RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Name("data"));               
  else GetDataPoissonInterval(rdataset_data_mj,rrv_mass_j,mplot);
		
  // Get the pull and plot it
  RooPlot* mplot_pull = get_pull(rrv_mass_j,mplot,rdataset_data_mj,model_data,rfresult,"data","model_mc",1,1);
  
  // signal window zone with vertical lines
  TLine* lowerLine = new TLine(mj_signal_min,0.,mj_signal_min,mplot->GetMaximum()*0.9); lowerLine->SetLineWidth(2); lowerLine->SetLineColor(kGray+2); lowerLine->SetLineStyle(9);
  TLine* upperLine = new TLine(mj_signal_max,0.,mj_signal_max,mplot->GetMaximum()*0.9); upperLine->SetLineWidth(2); upperLine->SetLineColor(kGray+2); upperLine->SetLineStyle(9);
  mplot->addObject(lowerLine);
  mplot->addObject(upperLine);

  // legend of the plot
  TLegend* leg = legend4Plot(mplot,0,-0.2,0.07,0.04,0.,1,"em");
  mplot->addObject(leg);
  mplot->GetYaxis()->SetRangeUser(1e-2,mplot->GetMaximum()*1.5);

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

  Name.Form("m_j_sideband%s",label.c_str());
  draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),std::string(Name.Data()),model,"em",0,1,GetLumi());

                
  // to calculate the WJets's normalization and error in M_J signal_region. The error must contain the shape error: model_WJets have new parameters fitting data
  if(ttbarcontrolregion){

    RooAbsReal* fullInt   = model_TTbar->createIntegral(RooArgSet(*rrv_mass_j),RooArgSet(*rrv_mass_j) );
    RooAbsReal* signalInt = model_TTbar->createIntegral(RooArgSet(*rrv_mass_j),RooArgSet(*rrv_mass_j),("signal_region"));
    double fullInt_val   = fullInt->getVal();
    double signalInt_val = signalInt->getVal()/fullInt_val;
    // take the value from the fit (normalization) and multiply it from the ratio of the integrals
    RooRealVar* rrv_number_TTbar_in_mj_signal_region_from_fitting = new RooRealVar(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*signalInt_val);

    // Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
    rrv_number_TTbar_in_mj_signal_region_from_fitting->setError(Calc_error_extendPdf(rdataset_data_mj,dynamic_cast<RooExtendPdf*>(model_TTbar),rfresult,"signal_region"));
    std::cout<< "########## error on the normaliztion due to shape + norm = "<<rrv_number_TTbar_in_mj_signal_region_from_fitting->getError()<<std::endl;;
    workspace->import(*rrv_number_TTbar_in_mj_signal_region_from_fitting);
    rrv_number_TTbar_in_mj_signal_region_from_fitting->Print();

    RooRealVar* rrv_TTbar  = workspace->var(("rrv_number+"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str()); //nominal parametrization for Wjets                        
    rrv_TTbar->Print();

    RooRealVar* rrv_STop    = workspace->var(("rrv_number_dataset_signal_region_STop_"+channel+"_mj").c_str());
    rrv_STop->Print();

    RooRealVar* rrv_WJets0   = workspace->var(("rrv_number_dataset_signal_region_WJtes0_"+channel+"_mj").c_str());
    rrv_WJets0->Print();

    RooRealVar* rrv_TTbar0   = workspace->var(("rrv_number_dataset_signal_region_TTbar_+"+channel+"_mj").c_str());
    rrv_TTbar0->Print();

    RooRealVar* rrv_VV      = workspace->var(("rrv_number_dataset_signal_region_VV_"+channel+"_mj").c_str());
    rrv_VV->Print();

    RooRealVar* rrv_WW_EWK  = workspace->var(("rrv_number_dataset_signal_region_WW_EWK_"+channel+"_mj").c_str());
    rrv_WW_EWK->Print();

  }
  else{
         RooAbsReal* fullInt   = model_WJets->createIntegral(RooArgSet(*rrv_mass_j),RooFit::NormSet(RooArgSet(*rrv_mass_j)));
         RooAbsReal* signalInt = model_WJets->createIntegral(RooArgSet(*rrv_mass_j),RooFit::NormSet(RooArgSet(*rrv_mass_j)),RooFit::Range("signal_region"));
         double fullInt_val   = fullInt->getVal();
         double signalInt_val = signalInt->getVal()/fullInt_val;
         //take the value from the fit (normalization) and multiply it from the ratio of the integrals
         RooRealVar* rrv_number_WJets_in_mj_signal_region_from_fitting = new RooRealVar(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str(),workspace->var(("rrv_number"+label+model+"_"+channel+"_mj").c_str())->getVal()*signalInt_val);

         // Error on the normalization --> from a dedicated function taking into account shape uncertainty on the parameters that are floating in the fit)
         rrv_number_WJets_in_mj_signal_region_from_fitting->setError( Calc_error_extendPdf(rdataset_data_mj, dynamic_cast<RooExtendPdf*>(model_WJets), rfresult,"signal_region") );
	 std::cout<< "########## error on the normaliztion due to shape + norm = "<<rrv_number_WJets_in_mj_signal_region_from_fitting->getError()<<std::endl;
         workspace->import(*rrv_number_WJets_in_mj_signal_region_from_fitting);
         rrv_number_WJets_in_mj_signal_region_from_fitting->Print();

         RooRealVar* rrv_WJets0  = workspace->var(("rrv_number+"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str()); //nominal parametrization for Wjets                        
         rrv_WJets0->Print();

         RooRealVar* rrv_STop    = workspace->var(("rrv_number_dataset_signal_region_STop_"+channel+"_mj").c_str());
         rrv_STop->Print();

         RooRealVar* rrv_TTbar   = workspace->var(("rrv_number_dataset_signal_region_TTbar_"+channel+"_mj").c_str());
         rrv_TTbar->Print();

         RooRealVar* rrv_WJets   = workspace->var(("rrv_number_dataset_signal_region_WJets0_+"+channel+"_mj").c_str());
         rrv_WJets->Print();

         RooRealVar* rrv_VV      = workspace->var(("rrv_number_dataset_signal_region_VV_"+channel+"_mj").c_str());
         rrv_VV->Print();

         RooRealVar* rrv_WW_EWK  = workspace->var(("rrv_number_dataset_signal_region_WW_EWK_"+channel+"_mj").c_str());
         rrv_WW_EWK->Print();

   }


}

//Method to fit data mlvj shape in the sideband -> first step for the background extraction of the shape                                                                           
void fit_mlvj_in_Mj_sideband(RooWorkspace* workspace, std::map<std::string,int> color_palet, const std::string & label, const std::string & mlvj_region, const std::string & mlvj_model, const std::string & channel, const std::string & wtagger, const std::string & fgen, const int & ttbarcontrolregion, const int & logy, const int & pseudodata){

  std::cout<<"############### Fit mlvj in mj sideband: "<<label<<" "<<mlvj_region<<" "<<mlvj_model<<" ##################"<<std::endl;

  RooRealVar* rrv_mass_lvj = workspace->var("rrv_mass_lvj");
  RooDataSet* rdataset_data_mlvj = (RooDataSet*) workspace->data(("rdataset_data+"+mlvj_region+"_"+channel+"_mlvj").c_str());

  // get and fix the minor component shapes in the sb low                                                                                                                            
  RooAbsPdf* model_VV_backgrounds     = get_VV_mlvj_Model(workspace,"_VV",mlvj_region,fgen,channel);
  RooAbsPdf* model_STop_backgrounds   = get_STop_mlvj_Model(workspace,"_STop",mlvj_region,fgen,channel);
  RooAbsPdf* model_WW_EWK_backgrounds = get_WW_EWK_mlvj_Model(workspace,"_WW_EWK",mlvj_region,fgen,channel);

  RooAbsPdf* model_TTbar_backgrounds = NULL ; 
  RooAbsPdf* model_WJets_backgrounds = NULL ; 

  if(ttbarcontrolregion == 0)
     model_TTbar_backgrounds  = get_TTbar_mlvj_Model(workspace,"_TTbar",mlvj_region,fgen,channel);
  else
     model_WJets_backgrounds  = get_WJets_mlvj_Model(workspace,"_WJets0",mlvj_region,fgen,channel);

  workspace->var(("rrv_number_TTbar"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_WJets0"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_STop"+mlvj_region+fgen,+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->Print();
  workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->Print();

  //Make the Pdf for the WJets                                                                                                                                                     
  RooArgList*   constraint = NULL;
  RooRealVar*   number_WJets_sb_lo = NULL; 
  RooRealVar*   number_TTbar_signal_region = NULL;
  RooExtendPdf* model_WJets = NULL; 
  RooExtendPdf* model_TTbar = NULL ;
  RooAddPdf*    model_data  = NULL ;

  if(ttbarcontrolregion == 0){

    RooAbsPdf* model_pdf_WJets = MakeGeneralPdf(workspace,label+mlvj_region+mlvj_model+std::string("_from_fitting"),mlvj_model,"_mlvj",wtagger,channel,constraint);
    model_pdf_WJets->Print();
    //inititalize the value to what was fitted with the mc in the sideband                                                                                                          
    number_WJets_sb_lo = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));

    model_WJets = new RooExtendPdf(("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),*model_pdf_WJets,*number_WJets_sb_lo);
    number_WJets_sb_lo->Print();


    // Add the other bkg component fixed to the total model --> in the extended way
    model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_WJets,*model_VV_backgrounds,*model_TTbar_backgrounds,*model_STop_backgrounds,*model_WW_EWK_backgrounds));
  }
  else{

         RooAbsPdf* model_pdf_TTbar = MakeGeneralPdf(workspace,label+mlvj_region+mlvj_model+"_from_fitting",mlvj_model,"_mlvj",wtagger,channel,constraint);
         model_pdf_TTbar->Print();
         // inititalize the value to what was fitted with the mc in the sideband
         number_TTbar_signal_region = dynamic_cast<RooRealVar*>( workspace->var(("rrv_number"+label+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->clone(("rrv_number"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str()));
    
         model_TTbar        = new  RooExtendPdf(("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),("model"+label+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str(),*model_pdf_TTbar,*number_TTbar_signal_region);
         number_TTbar_signal_region->Print();

         // Add the other bkg component fixed to the total model --> in the extended way
         model_data = new RooAddPdf(("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),("model_data"+label+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),RooArgList(*model_TTbar,*model_VV_backgrounds,*model_WJets_backgrounds,*model_STop_backgrounds,*model_WW_EWK_backgrounds));
        
  }
 
  RooFitResult* rfresult = model_data->fitTo(*rdataset_data_mlvj, RooFit::Save(1) ,RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE));
  rfresult = model_data->fitTo(*rdataset_data_mlvj, RooFit::Save(1) ,RooFit::Extended(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();
  workspace->import(*model_data);
  
  if(ttbarcontrolregion == 0){
    model_WJets->getParameters(*rdataset_data_mlvj)->Print("v");

    // data in the sideband plus error from fit        
    RooRealVar* rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data_"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                            ("rrv_number_data_"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                            workspace->var(("rrv_number_TTbar"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_WJets0"+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

    rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number_WJets0"+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_WJets0"+mlvj_region+mlvj_model+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_TTbar"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_TTbar"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()+
                                                    workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                    workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()));

    rrv_number_data_sb_lo_mlvj->Print();
    workspace->import(*rrv_number_data_sb_lo_mlvj);
  }

  else{

    model_TTbar->getParameters(*rdataset_data_mlvj)->Print("v");

    // data in the sideband plus error from fit        
    RooRealVar* rrv_number_data_sb_lo_mlvj = new RooRealVar(("rrv_number_data_"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                            ("rrv_number_data_"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str(),
                                                            workspace->var(("rrv_number_TTbar"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getVal()+
                                                            workspace->var(("rrv_number_WJets0"+mlvj_region+fgen+"_from_fitting_"+channel+"_mlvj").c_str())->getVal() );

    rrv_number_data_sb_lo_mlvj->setError(TMath::Sqrt(workspace->var(("rrv_number_WJets0"+mlvj_region+fgen+"_from_fitting_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_WJets0"+mlvj_region+fgen+"_from_fitting_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_TTbar"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_TTbar"+mlvj_region+mlvj_model+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_STop"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_WW_EWK"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()+
                                                     workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()*
                                                     workspace->var(("rrv_number_VV"+mlvj_region+fgen+"_"+channel+"_mlvj").c_str())->getError()));

    rrv_number_data_sb_lo_mlvj->Print();
    workspace->import(*rrv_number_data_sb_lo_mlvj);

  }

  RooPlot* mplot = NULL; 
  //plot for WJets default + default shape
  if(label == "_WJets0"){

      mplot = rrv_mass_lvj->frame(RooFit::Title("M_lvj fitted in M_j sideband"),RooFit::Bins(int(rrv_mass_lvj->getBins())));

      rdataset_data_mlvj->plotOn(mplot,RooFit::Invisible(),RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0),RooFit::Invisible(),RooFit::Name("data_invisible"));

      model_data->plotOn(mplot, RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WJets"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WJets"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["WW_EWK"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["VV"]), RooFit::LineColor(kBlack), RooFit::VLines());

      model_data->plotOn(mplot, RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar"),RooFit::DrawOption("F"),RooFit::FillColor(color_palet["TTbar"]), RooFit::LineColor(kBlack),RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components(std::string(model_STop_backgrounds->GetName()).c_str()), RooFit::Name("STop"),RooFit::DrawOption("F"), RooFit::FillColor(color_palet["STop"]), RooFit::LineColor(kBlack),RooFit::VLines());

      //solid line
      model_data->plotOn(mplot,RooFit::Components((std::string(model_WJets->GetName())+","+std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WJets_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())+","+std::string(model_WW_EWK_backgrounds->GetName())).c_str()), RooFit::Name("WW_EWK_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())+","+std::string(model_VV_backgrounds->GetName())).c_str()), RooFit::Name("VV_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_TTbar_backgrounds->GetName())+","+std::string(model_STop_backgrounds->GetName())).c_str()), RooFit::Name("TTbar_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

      model_data->plotOn(mplot,RooFit::Components((std::string(model_STop_backgrounds->GetName()).c_str())), RooFit::Name("STop_line_invisible"), RooFit::LineColor(kBlack), RooFit::LineWidth(2), RooFit::VLines());

  }
  //plot for WJets default + default shape
  if(label == "_TTbar"){

            mplot = rrv_mass_lvj->frame(RooFit::Title("M_lvj fitted in TTbar Control Region "), RooFit::Bins(int(rrv_mass_lvj->getBins())));

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

   //draw the error band
   draw_error_band(rdataset_data_mlvj, model_data,workspace->var(("rrv_number_data_"+mlvj_region+mlvj_region+"_"+channel+"_mlvj").c_str()),rfresult,mplot,color_palet["Uncertainty"],"F");
   model_data->plotOn( mplot , RooFit::VLines(), RooFit::Invisible());
   model_data->plotOn( mplot , RooFit::Invisible(), RooFit::Name("model_mc"));

   if(pseudodata == 1) rdataset_data_mlvj->plotOn( mplot , RooFit::MarkerSize(1.5), RooFit::DataError(RooAbsData::SumW2), RooFit::XErrorSize(0), RooFit::Name("data"));     
   else GetDataPoissonInterval(rdataset_data_mlvj,rrv_mass_lvj,mplot);

   mplot->GetYaxis()->SetRangeUser(1e-2,mplot->GetMaximum()*1.2);
            
   //Add the legend to the plot
   TLegend* leg = legend4Plot(mplot,0,0.,0.06,0.16,0.,1,"em");
   mplot->addObject(leg);

   //### get the pull plot and store the canvas
   RooPlot* mplot_pull = get_pull(rrv_mass_lvj,mplot,rdataset_data_mlvj,model_data,rfresult,"data","model_mc",1,1);
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
   mplot->GetYaxis()->SetRangeUser(1e-2,mplot->GetMaximum()*1.2);

   /*if shapetest == 1:
     file_out_FTest.write(" ###################### \n");
     file_out_FTest.write(" Model pdf %s"%(model_data->GetName()));
     file_out_FTest.write(" Chi2 Chi2var %0.2f "%(chi_over_ndf));
     file_out_FTest.write(" Residual %0.2f   Nbin %0.2f nparameters %0.2f \n"%(residual,Nbin,nparameters));
   */

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
   draw_canvas_with_pull(mplot,mplot_pull,new RooArgList(*parameters_list),std::string(command.Data()),std::string(Name.Data()),mlvj_model,"em",0,logy,GetLumi());

}  


//method to get the alpha function to extrapolate the wjets in the signal region
void get_WJets_mlvj_correction_sb_lo_to_signal_region(RooWorkspace* workspace,const std::string & label, const std::string & mlvj_model, const std::string & mass_spectrum,const std::string & workspace_label, const std::string & channel, const std::string & wtagger, const int & narrow_factor){

  std::cout<<" ############# get the extrapolation function alpha from MC : "<<label<<"   "<<mlvj_model<<"   "<<mass_spectrum<<" ###############"<<std::endl;
  gStyle->SetPadRightMargin(0.08);
  gStyle->SetPadTickY(0);
  gStyle->cd();

  //take input var and datasets from 4fit collection --> mc not scaled to lumi --> just a shape here
  RooRealVar* rrv_x = workspace->var("rrv_mass_lvj");
  RooDataSet* rdataset_WJets_sb_lo_mlvj = (RooDataSet*) workspace->data(("rdataset"+workspace_label+label+"_sb_lo_"+channel+mass_spectrum).c_str());
  RooDataSet* rdataset_WJets_signal_region_mlvj = (RooDataSet*) workspace->data(("rdataset"+workspace_label+label+"_signal_region_"+channel+mass_spectrum).c_str());

  //### create a frame for the next plots
  RooPlot* mplot = rrv_x->frame(RooFit::Title("correlation_pdf"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor))) ;
  mplot->GetYaxis()->SetTitle("arbitrary units");

  RooAbsPdf* correct_factor_pdf = NULL ;

  //Model used for Higgs analysis --> parameters in the SR has to be fitted, not yet done in order to take into account correlations between mj and mlvj
  if(mlvj_model == "ErfExp_v1"){

    RooRealVar* rrv_c_sb       = workspace->var(("rrv_c_ErfExp"+label+"_sb_lo"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_offset_sb  = workspace->var(("rrv_offset_ErfExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());            

    RooRealVar* rrv_width_sb = NULL ;
    if(mass_spectrum == "_mlvj_relaxed") rrv_width_sb   = workspace->var(("rrv_width_ErfExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    else{
          rrv_width_sb   = workspace->var(("rrv_width_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str());
          rrv_width_sb->setConstant(kTRUE);               
    }
            
    RooRealVar* rrv_delta_c      = new RooRealVar(("rrv_delta_c_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal() -4*rrv_c_sb->getError(), workspace->var(("rrv_c_ErfExp"+label+"_signal_region"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());

    RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_offset_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_offset_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),workspace->var(("rrv_offset_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()-4*rrv_offset_sb->getError(), workspace->var(("rrv_offset_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+rrv_offset_sb->getError());

    RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_width_ErfExp"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal() -4*rrv_width_sb->getError(),workspace->var(("rrv_width_ErfExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError());

    RooFormulaVar* rrv_c_sr      = new RooFormulaVar(("rrv_c_sr"+label+"_"+channel+mass_spectrum).c_str(),      "@0+@1",RooArgList(*rrv_c_sb,*rrv_delta_c));
    RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
    RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(),  "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

    correct_factor_pdf = new RooAlpha("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c_sb,*rrv_offset_sb,*rrv_width_sb,rrv_x->getMin(),rrv_x->getMax());
   
  }

  if(mlvj_model == "ErfPow_v1"){

   RooRealVar* rrv_c_sb      = workspace->var(("rrv_c_ErfPow"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPow"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_width_sb  = workspace->var(("rrv_width_ErfPow"+label+"_sb_lo_"+channel+mass_spectrum).c_str());

   RooRealVar* rrv_delta_c      = new RooRealVar(("rrv_delta_c_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_c_sb->getError(),100*rrv_c_sb->getError());
   RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_offset_sb->getError(),100*rrv_offset_sb->getError());
   RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfPow"+label+"_"+channel+mass_spectrum).c_str(),"",0.,-100*rrv_width_sb->getError(),100*rrv_width_sb->getError());

   RooFormulaVar* rrv_c_sr      = new RooFormulaVar(("rrv_c_sr"+label+"_"+channel+mass_spectrum).c_str(),      "@0+@1",RooArgList(*rrv_c_sb,*rrv_delta_c));
   RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
   RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(),  "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

   correct_factor_pdf = new RooAlpha4ErfPowPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "ErfPow2_v1"){

   RooRealVar* rrv_c0_sb     = workspace->var(("rrv_c0_ErfPow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_c1_sb     = workspace->var(("rrv_c1_ErfPow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
   RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());

   RooRealVar* rrv_width_sb = NULL;
   if(mass_spectrum == "_mlvj_relaxed") rrv_width_sb   = workspace->var(("rrv_width_ErfPow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
   else{
     rrv_width_sb   = workspace->var(("rrv_width_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str());
     rrv_width_sb->setConstant(kTRUE);
   }
   
   RooRealVar* rrv_delta_c0  = new RooRealVar(("rrv_delta_c0_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal(),workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal() -4*rrv_c0_sb->getError(), workspace->var(("rrv_c0_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()+4*rrv_c0_sb->getError());

   RooRealVar* rrv_delta_c1     = new RooRealVar(("rrv_delta_c1_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal(),workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal() -4*rrv_c1_sb->getError(), workspace->var(("rrv_c1_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()+4*rrv_c1_sb->getError());

   RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_offset_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal() -4*rrv_offset_sb->getError(), workspace->var(("rrv_offset_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+4*rrv_offset_sb->getError());            

   RooRealVar* rrv_delta_width  = new RooRealVar(("rrv_delta_width_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_width_ErfPow2"+label+"_"+channel+mass_spectrum).c_str(),workspace->var(("rrv_width_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),workspace->var(("rrv_width_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal() -4*rrv_width_sb->getError(), workspace->var(("rrv_width_ErfPow2"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError());
                
   RooFormulaVar* rrv_c0_sr     = new RooFormulaVar(("rrv_c0_sr"+label+"_"+channel+mass_spectrum).c_str(),"@0+@1",RooArgList(*rrv_c0_sb,*rrv_delta_c0));
   RooFormulaVar* rrv_c1_sr     = new RooFormulaVar(("rrv_c1_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_c1_sb,*rrv_delta_c1));
   
   RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
   RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_width_sb,*rrv_delta_width));

   correct_factor_pdf = new RooAlpha4ErfPow2Pdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c0_sr,*rrv_c1_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c0_sb,*rrv_c1_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "ErfPowExp_v1"){ //take initial value from what was already fitted in the SR

     RooRealVar* rrv_c0_sb     = workspace->var(("rrv_c0_ErfPowExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
     RooRealVar* rrv_c1_sb     = workspace->var(("rrv_c1_ErfPowExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
     RooRealVar* rrv_offset_sb = workspace->var(("rrv_offset_ErfPowExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());

     RooRealVar* rrv_width_sb = NULL ;
     if(mass_spectrum == "_mlvj_relaxed")
        rrv_width_sb      = workspace->var(("rrv_width_ErfPowExp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
     else{
           rrv_width_sb   = workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str());
           rrv_width_sb->setConstant(kTRUE);               
     }

     RooRealVar* rrv_delta_c0  = new RooRealVar(("rrv_delta_c0_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                        workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal(),
					workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()-4*rrv_c0_sb->getError(),
                                        workspace->var(("rrv_c0_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c0_sb->getVal()+4*rrv_c0_sb->getError());

     RooRealVar* rrv_delta_c1 = new RooRealVar(("rrv_delta_c1_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                       workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal(),
                                       workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()-4*rrv_c1_sb->getError(),
				       workspace->var(("rrv_c1_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c1_sb->getVal()+4*rrv_c1_sb->getError());

     RooRealVar* rrv_delta_offset = new RooRealVar(("rrv_delta_offset_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                                   ("rrv_delta_offset_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                         workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal(),
					 workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()-4*rrv_offset_sb->getError(),workspace->var(("rrv_offset_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_offset_sb->getVal()+4*rrv_offset_sb->getError());

    RooRealVar*  rrv_delta_width = new RooRealVar(("rrv_delta_width_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_width_ErfPowExp"+label+"_"+channel+mass_spectrum).c_str(),
                                         workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal(),
                                         workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()-4*rrv_width_sb->getError(),workspace->var(("rrv_width_ErfPowExp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_width_sb->getVal()+4*rrv_width_sb->getError() );

    RooFormulaVar* rrv_c0_sr     = new RooFormulaVar(("rrv_c0_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",    RooArgList(*rrv_c0_sb,*rrv_delta_c0));
    RooFormulaVar* rrv_c1_sr     = new RooFormulaVar(("rrv_c1_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",    RooArgList(*rrv_c1_sb,*rrv_delta_c1));
    RooFormulaVar* rrv_offset_sr = new RooFormulaVar(("rrv_offset_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_offset_sb,*rrv_delta_offset));
    RooFormulaVar* rrv_width_sr  = new RooFormulaVar(("rrv_width_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1", RooArgList(*rrv_width_sb,*rrv_delta_width));

    correct_factor_pdf = new RooAlpha4ErfPowExpPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_c0_sr,*rrv_c1_sr,*rrv_offset_sr,*rrv_width_sr,*rrv_c0_sb,*rrv_c1_sb,*rrv_offset_sb,*rrv_width_sb);

  }

  if(mlvj_model == "Exp"){
 
    RooRealVar* rrv_c_sb    = workspace->var(("rrv_c_Exp"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_Exp"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_Exp"+label+"_"+channel+mass_spectrum).c_str(),
                                      workspace->var(("rrv_c_Exp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),
                                      workspace->var(("rrv_c_Exp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()-4*rrv_c_sb->getError(),
				      workspace->var(("rrv_c_Exp"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());

    correct_factor_pdf = new RooExponential("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c);

  }

  if(mlvj_model == "Pow"){

    RooRealVar* rrv_c_sb    = workspace->var(("rrv_c_Pow"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_Pow"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_Pow"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c_sb->getError(),100*rrv_c_sb->getError());
    correct_factor_pdf = new RooPowPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c);
   
  }

  if(mlvj_model == "ExpN"){
            
    RooRealVar* rrv_c_sb  = workspace->var(("rrv_c_ExpN"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_n_sb  = workspace->var(("rrv_n_ExpN"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_delta_c = new RooRealVar(("rrv_delta_c_ExpN"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c_ExpN"+label+"_"+channel+mass_spectrum).c_str(),
                                      workspace->var(("rrv_c_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal(),
                                      workspace->var(("rrv_c_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()-4*rrv_c_sb->getError(),
				      workspace->var(("rrv_c_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_c_sb->getVal()+4*rrv_c_sb->getError());
    RooRealVar* rrv_delta_n = new RooRealVar(("rrv_delta_n_ExpN"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_n_ExpN"+label+"_"+channel+mass_spectrum).c_str(),
                                      workspace->var(("rrv_n_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal(),
                                      workspace->var(("rrv_n_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal()-4*rrv_n_sb->getError(),
				      workspace->var(("rrv_n_ExpN"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_n_sb->getVal()+4*rrv_n_sb->getError());

    correct_factor_pdf = new RooExpNPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c,*rrv_delta_n);

  }

  if (mlvj_model == "ExpTail"){

      RooRealVar* rrv_s_sb = workspace->var(("rrv_s_ExpTail"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
      RooRealVar* rrv_a_sb = workspace->var(("rrv_a_ExpTail"+label+"_sb_lo_"+channel+mass_spectrum).c_str());

      RooRealVar* rrv_delta_s = new RooRealVar(("rrv_delta_s_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_s_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),
                                      workspace->var(("rrv_s_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal(),
                                      workspace->var(("rrv_s_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal()-4*rrv_s_sb->getError(),
				      workspace->var(("rrv_s_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_s_sb->getVal()+4*rrv_s_sb->getError());
      RooRealVar* rrv_delta_a = new RooRealVar(("rrv_delta_a_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_a_ExpTail"+label+"_"+channel+mass_spectrum).c_str(),
                                      workspace->var(("rrv_a_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal(),
                                      workspace->var(("rrv_a_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal()-4*rrv_a_sb->getError(),
				      workspace->var(("rrv_a_ExpTail"+label+"_signal_region_"+channel+mass_spectrum).c_str())->getVal()-rrv_a_sb->getVal()+4*rrv_a_sb->getError());

      RooFormulaVar* rrv_a_sr = new RooFormulaVar(("rrv_a_ExpTail_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_a_sb,*rrv_delta_a ) );
      RooFormulaVar* rrv_s_sr = new RooFormulaVar(("rrv_s_ExpTail_sr"+label+"_"+channel+mass_spectrum).c_str(), "@0+@1",RooArgList(*rrv_s_sb,*rrv_delta_s ) );

      correct_factor_pdf = new RooAlpha4ExpTailPdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_s_sr,*rrv_a_sr,*rrv_s_sb,*rrv_a_sb);

  }
        
  if(mlvj_model == "Pow2"){

    RooRealVar* rrv_c0_sb    = workspace->var(("rrv_c0_Pow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
    RooRealVar* rrv_c1_sb    = workspace->var(("rrv_c1_Pow2"+label+"_sb_lo_"+channel+mass_spectrum).c_str());

    RooRealVar* rrv_delta_c0 = new RooRealVar(("rrv_delta_c0_Pow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c0_Pow2"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c0_sb->getError(),100*rrv_c0_sb->getError());
    RooRealVar* rrv_delta_c1 = new RooRealVar(("rrv_delta_c1_Pow2"+label+"_"+channel+mass_spectrum).c_str(),("rrv_delta_c1_Pow2"+label+"_"+channel+mass_spectrum).c_str(),0., -100*rrv_c1_sb->getError(),100*rrv_c1_sb->getError());

    correct_factor_pdf = new RooPow2Pdf("correct_factor_pdf","correct_factor_pdf",*rrv_x,*rrv_delta_c0,*rrv_delta_c1);

  }

  //define the category and do the simultaneous fit taking the combined dataset of events in mlvj sb and sr

  RooCategory* data_category = new RooCategory("data_category","data_category");
  data_category->defineType("sideband");
  data_category->defineType("signal_region");

  RooDataSet* combData4fit = NULL;

  if(mass_spectrum == "_mlvj_relaxed")  combData4fit = (RooDataSet*) workspace->data(("combData4fit"+label+"_"+channel+"_relaxed").c_str());
  else combData4fit = (RooDataSet*) workspace->data(("combData4fit"+label+"_"+channel).c_str());           

  RooAbsPdf* model_pdf_sb_lo_WJets    = workspace->pdf(("model_pdf"+label+"_sb_lo_"+channel+mass_spectrum).c_str());
  RooProdPdf* model_pdf_signal_region_WJets = new RooProdPdf(("model_pdf+"+label+"_signal_region_"+channel+mass_spectrum).c_str(),("model_pdf"+label+"_signal_region_"+channel+mass_spectrum).c_str(),*model_pdf_sb_lo_WJets,*correct_factor_pdf);

  
  RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",*data_category);
  simPdf->addPdf(*model_pdf_sb_lo_WJets,"sideband");
  simPdf->addPdf(*model_pdf_signal_region_WJets,"signal_region");
  RooFitResult* rfresult = simPdf->fitTo(*combData4fit,RooFit::Save(kTRUE), RooFit::SumW2Error(kTRUE));
  rfresult = simPdf->fitTo(*combData4fit,RooFit::Save(kTRUE), RooFit::SumW2Error(kTRUE), RooFit::Minimizer("Minuit2"));
  rfresult->Print();
  rfresult->covarianceMatrix().Print();

  // Decorrelate the parameters in the alpha shape
  RooWorkspace* wsfit_tmp = new RooWorkspace(("wsfit_tmp"+label+"_sim"+mass_spectrum).c_str());
  std::cout<<"############### diagonalizer alpha "<<std::endl;;
  PdfDiagonalizer* Deco      = new PdfDiagonalizer(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum).c_str(),wsfit_tmp,*rfresult);
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
    model_pdf_sb_lo_WJets->plotOn(mplot_sb_lo,RooFit::Name("model_mc"));
    RooPlot* mplot_pull_sideband = get_pull(rrv_x,mplot,rdataset_WJets_sb_lo_mlvj,model_pdf_sb_lo_WJets,rfresult,"data","model_mc",1,1);

    RooArgSet* parameters_list     = model_pdf_sb_lo_WJets->getParameters(rdataset_WJets_sb_lo_mlvj);
    mplot_sb_lo->GetYaxis()->SetRangeUser(1e-2,mplot_sb_lo->GetMaximum()*1.2);

    std::string NameDir  = "plots_"+channel+"_"+wtagger+"_g1";
    std::string NamePlot = "m_lvj"+label+"_sb_lo_sim";


    if(mass_spectrum == "_mlvj_relaxed"){
      NameDir  = NameDir+"/other_relaxed/";
      draw_canvas_with_pull(mplot_sb_lo,mplot_pull_sideband,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,"em",1,1,GetLumi());
    }
    else{
      NameDir  = NameDir+"/other/";
      draw_canvas_with_pull(mplot_sb_lo,mplot_pull_sideband,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,"em",1,1,GetLumi());                
    }
    //only W+jets mc plots in the SR region
    RooPlot* mplot_signal_region = rrv_x->frame(RooFit::Title("WJets sr"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor)));

    rdataset_WJets_signal_region_mlvj->plotOn(mplot_signal_region,RooFit::MarkerSize(1.5),RooFit::DataError(RooAbsData::SumW2),RooFit::XErrorSize(0), RooFit::Name("data"));
    model_pdf_signal_region_WJets->plotOn(mplot_signal_region, RooFit::Name("model_mc"));
    RooPlot* mplot_pull_signal_region = get_pull(rrv_x,mplot_signal_region,rdataset_WJets_signal_region_mlvj,model_pdf_signal_region_WJets,rfresult,"data","model_mc",1,1);

    parameters_list = model_pdf_signal_region_WJets->getParameters(rdataset_WJets_signal_region_mlvj);
    mplot_signal_region->GetYaxis()->SetRangeUser(1e-2,mplot_signal_region->GetMaximum()*1.2);

    if(mass_spectrum == "_mlvj_relaxed"){
       NamePlot = "m_lvj"+label+"_signal_region_sim";
       draw_canvas_with_pull(mplot_signal_region,mplot_pull_signal_region,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,"em",1,1,GetLumi());
    }
    else{
       NamePlot = "m_lvj"+label+"_signal_region_sim";
       draw_canvas_with_pull(mplot_signal_region,mplot_pull_signal_region,new RooArgList(*parameters_list),NameDir,NamePlot,mlvj_model,"em",1,1,GetLumi());
    }    
     
  }
  
     //plot also what is get from other source if available : alternate PS and shape: 1 PS and 01 is shape or fitting function
     if(label == "_WJets0"){

       //### Total plot shape in sb_lo, sr and alpha
       model_pdf_sb_lo_WJets->plotOn(mplot,RooFit::Name("Sideband"),RooFit::LineStyle(10));
       model_pdf_signal_region_WJets->plotOn(mplot, RooFit::LineColor(kRed) ,RooFit::LineStyle(8), RooFit::Name("Signal Region"));
       correct_factor_pdf_deco->plotOn(mplot, RooFit::LineColor(kBlack),RooFit::Name("#alpha") );

       if(workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+channel+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange), RooFit::LineStyle(3),RooFit::Name("#alpha: Alternate PS") );

       if(workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kMagenta), RooFit::LineStyle(7),RooFit::Name("#alpha: Alternate Function") );

      if(workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_up_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_up_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kRed), RooFit::LineStyle(3),RooFit::Name("#alpha: jes up") );

      if( workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_dn_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jes_dn_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange+1), RooFit::LineStyle(7),RooFit::Name("#alpha: jes dn") );

      if( workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kAzure), RooFit::LineStyle(3),RooFit::Name("#alpha: jer") );

      if( workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_up_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_up_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kYellow+1), RooFit::LineStyle(7),RooFit::Name("#alpha: jer up") );

      if( workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_dn_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets0massvbf_jer_dn_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kGray), RooFit::LineStyle(3),RooFit::Name("#alpha: jer dn") );
     }

     RooArgList* paras = new RooArgList();
     // Make a list of paramters as a function of the model after decorrelation
     if(mlvj_model == "ErfExp_v1" or mlvj_model == "ErfPow_v1" or mlvj_model == "2Exp"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig4").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig5").c_str()));
     }
     if(mlvj_model == "ErfPow2_v1" or mlvj_model == "ErfPowExp_v1"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig4").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig5").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig6").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig7").c_str()));
     }
     if(mlvj_model == "Exp" or mlvj_model == "Pow"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
     }
     if(mlvj_model == "ExpN" or mlvj_model == "ExpTail" or mlvj_model == "Pow2"){
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig0").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig1").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig2").c_str()));
	     paras->add(*workspace->var(("Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum+"_eig3").c_str()));
     }
     
     if(label == "_WJets0" or label == "_WJets1"){ //### draw error band ar 1 and 2 sigma using the decorrelated shape
	      draw_error_band_shape_Decor(("correct_factor_pdf_Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum).c_str(),"rrv_mass_lvj",paras,workspace,double(1.),mplot,kGray+3,"F",3001,"#alpha #pm",20,400);
	      draw_error_band_shape_Decor(("correct_factor_pdf_Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum).c_str(),"rrv_mass_lvj",paras,workspace,double(2.),mplot,kGreen+2,"F",3001,"#alpha #pm",20,400);
	      draw_error_band_shape_Decor(("correct_factor_pdf_Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum).c_str(),"rrv_mass_lvj",paras,workspace,double(2.),mplot,kGray+3,"F",3001,"#alpha_invisible #pm",20,400);
     }

     //### plot on the same canvas
     correct_factor_pdf_deco->plotOn(mplot, RooFit::LineColor(kBlack),RooFit::Name("#alpha_invisible"));

     if(label == "_WJets0"){ //## add also the plot of alternate ps and function on the canvas
	      if(workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+label+"_"+wtagger+mass_spectrum).c_str()))
                workspace->pdf(("correct_factor_pdf_Deco_WJets1_sim_"+label+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kMagenta), RooFit::LineStyle(3),RooFit::Name("#alpha_invisible: Alternate PS") );
	      if(workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+channel+"_"+wtagger+mass_spectrum).c_str()))
		 workspace->pdf(("correct_factor_pdf_Deco_WJets01_sim_"+channel+"_"+wtagger+mass_spectrum).c_str())->plotOn(mplot, RooFit::LineColor(kOrange), RooFit::LineStyle(7),RooFit::Name("#alpha_invisible: Alternate Function"));
      }

      //### Add the legend
      TLegend* leg = legend4Plot(mplot,0,-0.07,0.02,0.01,0.,1,channel);
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
      std::string NamePlot = "correction_pdf"+label+"_"+mlvj_model+"_M_lvj_signal_region_to_sideband";

      if(mass_spectrum == "_mlvj_relaxed"){
        NameDir  = NameDir+"/other_relaxed/";
	draw_canvas(mplot,NameDir,NamePlot,channel,GetLumi(),0,1,0);
      }
      else{
        NameDir  = NameDir+"/other/";
	draw_canvas(mplot,NameDir,NamePlot,channel,GetLumi(),0,1,0);
      }
      correct_factor_pdf_deco->getParameters(rdataset_WJets_sb_lo_mlvj)->Print("v");
      workspace->Print();

      RooAbsPdf* model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco = workspace->pdf(("model_pdf"+label+"_sb_lo_from_fitting_"+channel+mass_spectrum+"_Deco_"+label+"_sb_lo_from_fitting_"+channel+"_"+wtagger+mass_spectrum).c_str());           
      model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco->Print("v");
    
      //### Wjets shape in the SR correctedfunction * sb
      RooProdPdf* model_pdf_WJets_signal_region_after_correct_mlvj = new RooProdPdf(("model_pdf"+label+"_signal_region_"+channel+"_after_correct"+mass_spectrum).c_str(),("model_pdf"+label+"_signal_region_"+channel+"_after_correct"+mass_spectrum).c_str(),*model_pdf_WJets_sb_lo_from_fitting_mlvj_Deco,*workspace->pdf(("correct_factor_pdf_Deco"+label+"_sim_"+channel+"_"+wtagger+mass_spectrum).c_str()));
      model_pdf_WJets_signal_region_after_correct_mlvj->Print();
      //### fix the parmaters and import in the workspace
      workspace->import(*model_pdf_WJets_signal_region_after_correct_mlvj);

      //##### calculate the normalization and alpha for limit datacard
      workspace->var(("rrv_number"+label+"_signal_region_"+channel+mass_spectrum).c_str())->Print();
      workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->Print();
      workspace->var(("rrv_number"+label+"_signal_region_"+channel+mass_spectrum).c_str())->setVal(workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->getVal());
      workspace->var(("rrv_number"+label+"_signal_region_"+channel+mass_spectrum).c_str())->setError(workspace->var(("rrv_number"+label+"_in_mj_signal_region_from_fitting_"+channel).c_str())->getError());
      workspace->var(("rrv_number"+label+"_signal_region_"+channel+mass_spectrum).c_str())->setConstant(kTRUE);
     
}
