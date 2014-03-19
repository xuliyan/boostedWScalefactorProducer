#include "PlotUtils.h"

void GetDataPoissonInterval(const RooAbsData* data, RooRealVar* rrv_x, RooPlot* mplot){
 
  TString Title ; Title.Form("%s_binnedClone",data->GetName());
  RooDataHist* datahist   = new RooDataHist(Title.Data(),Title.Data(),*rrv_x,*((RooDataSet*)data));
  TH1* data_histo         = datahist->createHistogram("histo_data",*rrv_x) ;
  RooHist* data_plot  = new RooHist(*data_histo);
  data_plot->SetMarkerStyle(20);
  data_plot->SetMarkerSize(1.5);
  data_plot->SetName("data");

  double alpha = 1 - 0.6827;
  for( int iPoint = 0 ; iPoint < data_plot->GetN(); iPoint++){

   double N = data_plot->GetY()[iPoint];
   double L , U ;
   if(N==0) L = 0;
   else L = (ROOT::Math::gamma_quantile(alpha/2,N,1.));
   U =  ROOT::Math::gamma_quantile_c(alpha/2,N+1,1);
   data_plot->SetPointEYlow(iPoint,N-L);
   data_plot->SetPointEYhigh(iPoint,U-N);
   data_plot->SetPointEXlow(iPoint,0);
   data_plot->SetPointEXhigh(iPoint,0);
  }

   mplot->addPlotable(data_plot,"PE");

}



// in order to get the pull
RooPlot* get_pull(RooRealVar* rrv_x, RooPlot* mplot_orig, RooDataSet* rdataset, RooAbsPdf* model, RooFitResult* rfresult, const std::string & dataname, const std::string & modelname, const int & makeBand, const int & narrow_factor){

  std::cout<<"############### draw the pull plot ########################"<<std::endl;
  RooHist* hpull = mplot_orig->pullHist(dataname.c_str(),modelname.c_str());
  
  double x = 0.; double y = 0.;

  double err_y_low = 0.; double err_y_high = 0. ;
  double err_x_low = 0.; double err_x_high = 0. ;
  for( int ipoint = 0 ;  ipoint < hpull->GetN(); ipoint++){
    hpull->GetPoint(ipoint,x,y);
    err_x_low = hpull->GetErrorXlow(ipoint);
    err_x_high = hpull->GetErrorXhigh(ipoint);
    err_y_low = hpull->GetErrorYlow(ipoint);
    err_y_high = hpull->GetErrorYhigh(ipoint);
    if(y == 0){
      hpull->SetPoint(ipoint,x,10);
      hpull->SetPointError(ipoint,err_x_low,err_x_high,err_y_low,err_y_high);
    }
  }
        
  RooPlot* mplot_pull = rrv_x->frame(RooFit::Title("Pull Distribution"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor)));
  
  TLine* medianLine = new TLine(rrv_x->getMin(),0.,rrv_x->getMax(),0); 
  medianLine->SetLineWidth(2); 
  medianLine->SetLineColor(kRed);

  if(makeBand) draw_error_band_extendPdf_pull(rdataset,model,rfresult,mplot_pull);

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
  mplot_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma_{data}}");
  mplot_pull->GetYaxis()->CenterTitle();

  return mplot_pull;

}


RooPlot* get_pull_ws(RooRealVar* rrv_x, RooPlot* mplot_orig, TGraphAsymmErrors* plot_graph, const std::string & dataname, const std::string & modelname, const int & narrow_factor){

 std::cout<<"############### draw the pull plot ########################"<<std::endl;
 RooHist* hpull = mplot_orig->pullHist(dataname.c_str(),modelname.c_str());
 double x = 0. ; double y = 0. ;
 double err_y_low = 0. ; double err_y_high = 0. ;
 double err_x_low = 0. ; double err_x_high = 0. ;
        
 for( int ipoint = 0; ipoint < hpull->GetN() ; ipoint++){ 
          
  hpull->GetPoint(ipoint,x,y);
  err_x_low  = hpull->GetErrorXlow(ipoint);
  err_x_high = hpull->GetErrorXhigh(ipoint);
  err_y_low  = hpull->GetErrorYlow(ipoint);
  err_y_high = hpull->GetErrorYhigh(ipoint);
  if(y == 0){
   hpull->SetPoint(ipoint,x,10);
   hpull->SetPointError(ipoint,err_x_low,err_x_high,err_y_low,err_y_high);
   }    
 }
 
 RooPlot* mplot_pull = rrv_x->frame(RooFit::Title("Pull Distribution"), RooFit::Bins(int(rrv_x->getBins()/narrow_factor)));
 TLine* medianLine = new TLine(rrv_x->getMin(),0.,rrv_x->getMax(),0);
 medianLine->SetLineWidth(2); 
 medianLine->SetLineColor(kRed);

 mplot_pull->addObject(plot_graph,"E3");        
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
 mplot_pull->GetYaxis()->SetTitle("#frac{Data-Fit}{#sigma_{data}}");
 mplot_pull->GetYaxis()->CenterTitle();

 return mplot_pull;

}

// method for draw the banner
TLatex* banner4Plot(const std::string & channel, const float & lumi, const int & iswithpull){

 std::cout<<"############### draw the banner ########################"<<std::endl;
 TString bannerName;
 if(channel=="mu") 
  bannerName.Form("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow #mu #nu",lumi);
 else if(channel=="el") 
  bannerName.Form("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow e #nu",lumi);
 else if(channel=="em")  
  bannerName.Form("CMS Preliminary, %.1f fb^{-1} at #sqrt{s} = 8 TeV, W#rightarrow l #nu",lumi);
 
 TLatex* banner = NULL ;
 
 if(iswithpull){
    banner = new TLatex(0.3,0.96,bannerName.Data());
    banner->SetNDC(); banner->SetTextSize(0.04);
 }
 else{
    banner = new TLatex(0.22,0.96,bannerName.Data());
    banner->SetNDC(); banner->SetTextSize(0.04);
 }                                                                                                         
      
 return banner;

}

// in order to make the legend                                                                                                                                                      
TLegend* legend4Plot(RooPlot* plot, const int & left, const double & x_offset_low, const double & y_offset_low,  const double & x_offset_high, const double & y_offset_high, const int & TwoCoulum, const std::string & channel ){
        
  std::cout<<"############### draw the legend ########################"<<std::endl;
  TLegend* theLeg = NULL ;
  if(left ==-1){  
    theLeg = new TLegend(0.65+x_offset_low,0.58+y_offset_low,0.93+x_offset_low,0.87+y_offset_low,"","NDC");
    theLeg->SetName("theLegend");
    theLeg->SetLineColor(0);
    theLeg->SetTextFont(42);
    theLeg->SetTextSize(.04);
  }
  else{
    theLeg = new TLegend(0.41+x_offset_low,0.61+y_offset_low,0.76+x_offset_high,0.93+y_offset_high,"","NDC");
    theLeg->SetName("theLegend");
    if(TwoCoulum) theLeg->SetNColumns(2);
  }
  
  theLeg->SetFillColor(0);
  theLeg->SetFillStyle(0);
  theLeg->SetBorderSize(0);
  theLeg->SetLineColor(0);
  theLeg->SetLineWidth(0);
  theLeg->SetLineStyle(0);
  theLeg->SetTextSize(0.040);
  theLeg->SetTextFont(42);

  int entryCnt = 0;
  std::string objName_before = "";
  TObject* objName_signal_graviton = NULL;
  std::string objNameLeg_signal_graviton = "";
  std::string legHeader ; 
  
  if(channel == "mu")      legHeader="(#mu#nu, 1JHP)";
  else if(channel == "el") legHeader="(e#nu, 1JHP)";
  else if(channel == "em") legHeader="(l#nu, 1JHP)";

  for( int obj = 0 ; obj < int(plot->numItems()); obj++){
          
   std::string objName = plot->nameOf(obj);
   if(objName == "errorband") objName = "Uncertainty";
   
   if( not ( ( (plot->getInvisible(objName.c_str())) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisible") or TString(objName).Contains("TLine") or objName ==objName_before )){
     
    TObject*  theObj = plot->getObject(obj);
    std::string objTitle = objName;
    std::string drawoption= plot->getDrawOptions(objName.c_str()).Data();
    if(drawoption == "P") drawoption = "PE";
    if(TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma")){ 
     objName_before = objName; 
     continue ;
    }
    else if(TString(objName).Contains("Graph")){
     objName_before = objName; 
     continue ;
    }
    else if(TString(objName) == "data"){
     theLeg->AddEntry(theObj,std::string("CMS Data "+legHeader).c_str(),"PE");  
     objName_before = objName;
    }
    else{
     objName_before = objName;
     continue ;
    }
   }
  }
  
  entryCnt = 0;
  objName_before = "";
  objName_signal_graviton = NULL;
  objNameLeg_signal_graviton = "";
  
  for( int obj = 0 ; obj < int(plot->numItems()) ; obj++ ){
    std::string objName = plot->nameOf(obj);
    if( objName == "errorband" ) objName = "Uncertainty";
    if(not ( ( (plot->getInvisible(objName.c_str())) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisible") or TString(objName).Contains("TLine") or objName == objName_before )){
       TObject* theObj = plot->getObject(obj);
       std::string objTitle = objName;
       std::string drawoption = plot->getDrawOptions(objName.c_str()).Data();
       if(drawoption == "P") drawoption = "PE";
       if(TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma")){
        objName_before = objName;
        continue ;
       }    
       else if(TString(objName).Contains("Graph")){
         objName_before = objName;
         continue ;
       }
       else if(TString(objName) == "WJets"){
        theLeg->AddEntry(theObj, "W+jets","F");
        objName_before = objName;
       }
       else{
         objName_before = objName; 
         continue ;
       }
     }
   }   
        
   entryCnt = 0;
   objName_before = "";
   objName_signal_graviton = NULL;
   objNameLeg_signal_graviton = "";

   
   for( int obj = 0 ; obj < int(plot->numItems()) ; obj++){
            
      std::string objName = plot->nameOf(obj);
      std::cout<<" Loop objName "<<objName<<std::endl;
      if(objName == "errorband") objName = "Uncertainty";
      if(not(( (plot->getInvisible(objName.c_str())) and (not TString(objName).Contains("Uncertainty")) ) or TString(objName).Contains("invisible") or TString(objName).Contains("TLine") or objName == objName_before)){
        std::cout<<" Skip invisible and double counting "<<objName<<std::endl;
         TObject* theObj = plot->getObject(obj);
         std::string objTitle = objName;
         std::string drawoption = plot->getDrawOptions(objName.c_str()).Data();
         if(drawoption == "P") drawoption = "PE";
         if(TString(objName).Contains("Uncertainty") or TString(objName).Contains("sigma"))
             theLeg->AddEntry(theObj, objName.c_str(),"F");
         else if(TString(objName).Contains("Graph")){
              if(not (objName_before == "Graph" or objName_before == "Uncertainty"))
               theLeg->AddEntry(theObj, "Uncertainty","F");
	 }
         else{
               if(TString(objName) == "STop") theLeg->AddEntry(theObj, "Single Top","F");
               else if(TString(objName) == "TTbar") theLeg->AddEntry(theObj, "t#bar{t}","F");
               else if(TString(objName) == "VV")    theLeg->AddEntry(theObj, "WW/WZ","F");
               else if(TString(objName) == "data") { objName_before = objName; entryCnt = entryCnt+1; continue;}
               else if(TString(objName) == "WJets"){ objName_before = objName; entryCnt = entryCnt+1; continue;}
               else if(TString(objName).Contains("vbfH")) 
                 theLeg->AddEntry(theObj,(TString(objName).ReplaceAll("vbfH","qqH")).Data() ,"L");
               else if(TString(objName).Contains("Uncertainty"))
                 theLeg->AddEntry(theObj,objTitle.c_str(),drawoption.c_str());
               else if(TString(objName).Contains("Bulk")){
               if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M600") or TString(objName).Contains("BulkG_WW_lvjj_c0p2_M600")){
                 objName_signal_graviton = theObj ;
                 objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.6 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M700")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.7 TeV #tilde{k}=0.5 (#times100)";
                }
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M800")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.8 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M900")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=0.9 TeV #tilde{k}=0.5 (#times100)";
                }   
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1000")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1100")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.1 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1200")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.2 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1300")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.3 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1400")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.4 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1500")){
                    objName_signal_graviton = theObj ;
                    objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.5 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1600") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1600")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.6 TeV #tilde{k}=0.5 (#times100)";
                }
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1700") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1700")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.7 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1800") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1800")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.8 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M1900") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M1900")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=1.9 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2000") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2000")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2 TeV #tilde{k}=0.5 (#times100)";
                }   
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2100") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2100")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.1 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2200") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2200")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.2 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2300") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2300")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.3 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2400") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2400")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.4 TeV #tilde{k}=0.5 (#times100)";
                }    
                if(TString(objName).Contains("BulkG_WW_inclusive_c0p2_M2500") or  TString(objName).Contains("BulkG_WW_lvjj_c0p2_M2500")){
                   objName_signal_graviton = theObj ;
                   objNameLeg_signal_graviton = "Bulk G* M_{G*}=2.5 TeV #tilde{k}=0.5 (#times100)";
		}
                else theLeg->AddEntry(theObj,objTitle.c_str(),drawoption.c_str());
	       }
               entryCnt = entryCnt+1;
	 }
      }
      objName_before = objName;
   }
  
   if(objName_signal_graviton != NULL){
    if(std::string(objName_signal_graviton->GetName()) != "") 
     theLeg->AddEntry(objName_signal_graviton,TString(objNameLeg_signal_graviton).Data(),"L");
   }

   return theLeg;
}

void draw_canvas(RooPlot* in_obj, const std::string & in_directory, const TString & in_file_name, const std::string & channel, const float & lumi, const int & in_range, const int & logy, const int & frompull){

  std::cout<<"############### draw the canvas without pull ########################"<<std::endl;
  TCanvas cMassFit ("cMassFit","cMassFit", 600,600);

  if(frompull and logy)
     in_obj->GetYaxis()->SetRangeUser(1e-2,in_obj->GetMaximum()/100);
  else if(not frompull and logy)
     in_obj->GetYaxis()->SetRangeUser(0.00001,in_obj->GetMaximum());
   

  if(in_range){
    TH2F h2("h2","",100,400,1400,4,0.00001,4);
    h2.Draw();
    in_obj->Draw("same");
  }  
  else in_obj->Draw();
    
  in_obj->GetXaxis()->SetTitleSize(0.045);
  in_obj->GetXaxis()->SetTitleOffset(1.15);
  in_obj->GetXaxis()->SetLabelSize(0.04);

  in_obj->GetYaxis()->SetTitleSize(0.045);
  in_obj->GetYaxis()->SetTitleOffset(1.40);
  in_obj->GetYaxis()->SetLabelSize(0.04);

  TLatex* banner = banner4Plot(channel,lumi,0);
  banner->Draw();
        
  TString Directory(in_directory);
  if(not Directory.EndsWith("/")) Directory = Directory.Append("/");
  system(std::string("mkdir -p "+Directory).c_str());

  TString rlt_file(Directory.Data()+in_file_name);
  if(rlt_file.EndsWith(".root")) rlt_file.ReplaceAll(".root","_rlt_without_pull_and_paramters.png");
  else{
         rlt_file.ReplaceAll(".root","");
         rlt_file = rlt_file.Append(".png");
  }
  
  cMassFit.SaveAs(rlt_file.Data());

  rlt_file.ReplaceAll(".png",".pdf");
  cMassFit.SaveAs(rlt_file.Data());

  rlt_file.ReplaceAll(".pdf",".root");
  cMassFit.SaveAs(rlt_file.Data());

  if(logy){
      in_obj->GetYaxis()->SetRangeUser(1e-2,in_obj->GetMaximum()*100);
      cMassFit.SetLogy() ;
      cMassFit.Update();
      rlt_file.ReplaceAll(".root","_log.root");
      cMassFit.SaveAs(rlt_file.Data());
      rlt_file.ReplaceAll(".root",".pdf");
      cMassFit.SaveAs(rlt_file.Data());
      rlt_file.ReplaceAll(".pdf",".png");
      cMassFit.SaveAs(rlt_file.Data());
  }

}


// draw canvas with plots with pull
void draw_canvas_with_pull(RooPlot* mplot, RooPlot* mplot_pull, RooArgList* parameters_list, const std::string & in_directory, const std::string & in_file_name, const std::string & in_model_name, const std::string & channel, const int & show_parameter, const int & logy, const float & lumi){

  std::cout<<"############### draw the canvas with pull ########################"<<std::endl;
  mplot->GetXaxis()->SetTitleOffset(1.1);
  mplot->GetYaxis()->SetTitleOffset(1.3);
  mplot->GetXaxis()->SetTitleSize(0.055);
  mplot->GetYaxis()->SetTitleSize(0.055);
  mplot->GetXaxis()->SetLabelSize(0.045);
  mplot->GetYaxis()->SetLabelSize(0.045);
  mplot_pull->GetXaxis()->SetLabelSize(0.14);
  mplot_pull->GetYaxis()->SetLabelSize(0.14);
  mplot_pull->GetYaxis()->SetTitleSize(0.15);
  mplot_pull->GetYaxis()->SetNdivisions(205);

  TCanvas cMassFit ("cMassFit_Pull","cMassFit_Pull", 600,600);
  TIter par_first = parameters_list->createIterator();
  par_first.Reset();
  TObject* param_first = par_first.Next();
 
  TPad* pad1 = NULL ;
  TPad* pad2 = NULL ;
  TPad* pad3 = NULL ;
  
  if(param_first and show_parameter != 0){
         
   pad1 = new TPad("pad1","pad1",0.,0. ,0.8,0.24);
   pad2 = new TPad("pad2","pad2",0.,0.24,0.8,1. );
   pad3 = new TPad("pad3","pad3",0.8,0.,1,1);
   pad1->Draw();
   pad2->Draw();
   pad3->Draw();
  }
  else{
   pad1 = new TPad("pad1","pad1",0.,0. ,0.99,0.24);
   pad2 = new TPad("pad2","pad2",0.,0.24,0.99,1. );
   pad1->Draw();
   pad2->Draw();
  }       
   
  pad2->cd();
  mplot->Draw();
  TLatex* banner = banner4Plot(channel,lumi,1);
  banner->Draw();

  pad1->cd();
  mplot_pull->Draw();

  if(param_first and show_parameter != 0){

   pad3->cd();
   TLatex latex;
   latex.SetTextSize(0.1);
   TIter par = parameters_list->createIterator();
   par.Reset();
   TObject* param = par.Next();
   int i = 0;
   while(param){
     RooRealVar* PAR = dynamic_cast<RooRealVar*>(param);
     if((not PAR->isConstant() ) or show_parameter){
       param->Print();
       int icolor = 1;
       if(PAR->isConstant()) icolor = 2;
       TString Name ; Name = Form("#color[%d]{%s}",icolor,param->GetName());
       latex.DrawLatex(0,0.9-i*0.04,Name);
       Name = Form("#color[%d]{%4.3f +/- %2.1f}",icolor,PAR->getVal(),PAR->getError());
       latex.DrawLatex(0,0.9-i*0.04-0.02,Name);
       i=i+1;
       param = par.Next();
    }
   }          
  }      
  // create the directory where store the plots
  TString Directory(in_directory);
  if(not Directory.EndsWith("/")) Directory = Form("%s/",Directory.Data());
  system(("mkdir -p "+std::string(Directory)).c_str());
  
  TString rlt_file;  rlt_file.Form("%s%s",Directory.Data(),in_file_name.c_str());
  if(rlt_file.EndsWith(".root")){
    TString(in_model_name).ReplaceAll(".root","");
    rlt_file.ReplaceAll(".root","_"+in_model_name+"_with_pull.png");
  }   
  else{    
    TString(in_model_name).ReplaceAll(".root","");
    rlt_file.ReplaceAll(".root","");
    rlt_file = rlt_file.Append("_"+in_model_name+"_with_pull.png");
  }
  
  cMassFit.SaveAs(rlt_file.Data());
  rlt_file.ReplaceAll(".png",".pdf");
  cMassFit.SaveAs(rlt_file.Data());
  rlt_file.ReplaceAll(".pdf",".root");
  cMassFit.SaveAs(rlt_file.Data());
  TString string_file_name (in_file_name);
  if(string_file_name.EndsWith(".root"))
    string_file_name.ReplaceAll(".root","_"+in_model_name);
  else{
     string_file_name.ReplaceAll(".root","");
     string_file_name.Append("_"+in_model_name);
  }
  
  if(logy){
     mplot->GetYaxis()->SetRangeUser(1e-2,mplot->GetMaximum()*100);
     pad2->SetLogy() ;
     pad2->Update();
     cMassFit.Update();
     rlt_file.ReplaceAll(".root","_log.root");
     cMassFit.SaveAs(rlt_file.Data());
     rlt_file.ReplaceAll(".root",".pdf");
     cMassFit.SaveAs(rlt_file.Data());
     rlt_file.ReplaceAll(".pdf",".png");
     cMassFit.SaveAs(rlt_file.Data());
  }
 
  draw_canvas(mplot,in_directory,string_file_name,channel,lumi,0,logy,1);

}

// set tdr style function
void setTDRStyle(){

 
 //For the canvas:
 gStyle->SetCanvasBorderMode(0);
 gStyle->SetCanvasColor(kWhite);
 gStyle->SetCanvasDefH(600); 
 gStyle->SetCanvasDefW(600); 
 gStyle->SetCanvasDefX(0); 
 gStyle->SetCanvasDefY(0);
      
 //For the Pad:
 gStyle->SetPadBorderMode(0);
 gStyle->SetPadColor(kWhite);
 gStyle->SetPadGridX(kFALSE);
 gStyle->SetPadGridY(kFALSE);
 gStyle->SetGridColor(0);
 gStyle->SetGridStyle(3);
 gStyle->SetGridWidth(1);
      
 //For the frame:
 gStyle->SetFrameBorderMode(0);
 gStyle->SetFrameBorderSize(1);
 gStyle->SetFrameFillColor(0);
 gStyle->SetFrameFillStyle(0);
 gStyle->SetFrameLineColor(1);
 gStyle->SetFrameLineStyle(1);
 gStyle->SetFrameLineWidth(1);
      
 //For the histo:
 gStyle->SetHistLineColor(1);
 gStyle->SetHistLineStyle(0);
 gStyle->SetHistLineWidth(1);
 gStyle->SetEndErrorSize(2);
 gStyle->SetErrorX(0.);
 gStyle->SetMarkerStyle(20);
      
 //For the fit/function:
 gStyle->SetOptFit(1);
 gStyle->SetFitFormat("5.4g");
 gStyle->SetFuncColor(2);
 gStyle->SetFuncStyle(1);
 gStyle->SetFuncWidth(1);
      
 //For the date:
 gStyle->SetOptDate(0);
      
 //For the statistics box:
 gStyle->SetOptFile(0);
 gStyle->SetOptStat(0);
 gStyle->SetStatColor(kWhite);
 gStyle->SetStatFont(42);
 gStyle->SetStatFontSize(0.025);
 gStyle->SetStatTextColor(1);
 gStyle->SetStatFormat("6.4g");
 gStyle->SetStatBorderSize(1);
 gStyle->SetStatH(0.1);
 gStyle->SetStatW(0.15);
      
  //Margins:
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.06);
      
  //For the Global title:
  gStyle->SetOptTitle(0);
  gStyle->SetTitleFont(42);
  gStyle->SetTitleColor(1);
  gStyle->SetTitleTextColor(1);
  gStyle->SetTitleFillColor(10);
  gStyle->SetTitleFontSize(0.05);
      
  //For the axis titles:
  gStyle->SetTitleColor(1, "XYZ");
  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetTitleSize(0.03, "XYZ");
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(1.5);
      
  //For the axis labels:
  gStyle->SetLabelColor(1, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetLabelOffset(0.007, "XYZ");
  gStyle->SetLabelSize(0.03, "XYZ");
      
  //For the axis:
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(kTRUE);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1); 
  gStyle->SetPadTickY(1);
      
  //Change for log plots:
  gStyle->SetOptLogx(0);
  gStyle->SetOptLogy(0);
  gStyle->SetOptLogz(0);
      
  //Postscript options:
  gStyle->SetPaperSize(20.,20.);
  gStyle->cd();
}

float GetLumi(const std::string & channel){
 
  if(channel=="el") return 19.3;
  else if(channel=="mu") return 19.3;
  else if(channel=="em") return 19.3;

  return -1 ;
}