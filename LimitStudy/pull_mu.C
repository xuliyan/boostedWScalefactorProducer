#include <iostream>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <istream>
#include <sstream>

#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLatex.h"
#include "TMath.h"

#include "TF1.h"

int main (){

TFile* file = TFile::Open("EXP_TIGHT_600.root");

TTree* tree_s = (TTree*)file->Get("limit");

 Double_t limit, limit_temp, rInput=6;
 Double_t limitErr, err_dw_temp, err_up_temp;
 Float_t quantileExpected;

 tree_s->SetBranchAddress("limit",&limit);
 tree_s->SetBranchAddress("limitErr",&limitErr);
 tree_s->SetBranchAddress("quantileExpected",&quantileExpected);

 TH1F* pullDistribution = new TH1F ("pull","mu pull",15,-3,3);
 TH1F* pullnormDistribution = new TH1F ("pullnorm","mu pull norm",15,-3,3);
 TH1F* deltaDistribution = new TH1F ("delta","delta mu",15,-1*rInput,rInput);
 TH1F* sigmaDistribution = new TH1F ("sigma","sigma",15,0,rInput);

for( int iEntry = 0; iEntry < tree_s->GetEntries(); iEntry ++){

        tree_s->GetEntry(iEntry);
        if(quantileExpected==0.5){ limit_temp = limit ; continue ; }
        if(quantileExpected>0.15 && quantileExpected<0.17 ) {err_dw_temp = limit ;  continue ; }
        if(quantileExpected>0.83 && quantileExpected<0.85) {err_up_temp = limit ;  continue ; }
        if(quantileExpected==-1){
             deltaDistribution->Fill(limit_temp-rInput);
             sigmaDistribution->Fill((err_up_temp-err_dw_temp)/2);
             pullDistribution->Fill(2*(limit_temp-rInput)/(err_up_temp-err_dw_temp));
	     //	     std::cout<<limit_temp<<" "<<err_up_temp<<" "<<err_dw_temp<<std::endl;
                        }
              }
//fit and draw sigma distribution

 TF1* fs = new TF1("gaus","gaus",0,20);
 sigmaDistribution->Fit(fs,"R");

 TCanvas* c2 = new TCanvas ("c2","c2");
 c2->cd();
 c2->SetGrid();

 sigmaDistribution->GetYaxis()->SetTitle("N");
 sigmaDistribution->GetXaxis()->SetTitle("sigma");
 sigmaDistribution->SetMarkerStyle(20);
 sigmaDistribution->SetMarkerSize(0.9);
 fs->SetLineColor(4);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1111);

 sigmaDistribution->Draw("pe");
 fs->Draw("same");
 c2->Print("sigma.png",".png");
 c2->Close();



//fit and draw residual distribution

 TF1* fd = new TF1("gaus","gaus",-20,20);
 deltaDistribution->Fit(fd,"R");

 TCanvas* c3 = new TCanvas ("c3","c3");
 c3->cd();
 c3->SetGrid();

 deltaDistribution->GetYaxis()->SetTitle("N");
 deltaDistribution->GetXaxis()->SetTitle("residual");
 deltaDistribution->SetMarkerStyle(20);
 deltaDistribution->SetMarkerSize(0.9);
 fd->SetLineColor(4);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1111);

 deltaDistribution->Draw("pe");
 fd->Draw("same");
 c3->Print("residual.png",".png");
 c3->Close();

 //fit and draw pull distribution

 TF1* f = new TF1("gaus","gaus",-3,3);
 pullDistribution->Fit(f,"R");

 TCanvas* c1 = new TCanvas ("c1","c1");
 c1->cd();
 c1->SetGrid();

 pullDistribution->GetYaxis()->SetTitle("N");
 pullDistribution->GetXaxis()->SetTitle("pull");
 pullDistribution->SetMarkerStyle(20);
 pullDistribution->SetMarkerSize(0.9);
 f->SetLineColor(4);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1111);

 pullDistribution->Draw("pe");
 f->Draw("same");
 c1->Print("pull.png",".png");
 c1->Close();



for( int iEntry = 0; iEntry < tree_s->GetEntries(); iEntry ++){

        tree_s->GetEntry(iEntry);
        if(quantileExpected==0.5){ limit_temp = limit ; continue ; }
        if(quantileExpected>0.15 && quantileExpected<0.17 ) {err_dw_temp = limit ;  continue ; }
        if(quantileExpected>0.83 && quantileExpected<0.85) {err_up_temp = limit ;  continue ; }
        if(quantileExpected==-1){
	  pullnormDistribution->Fill((limit_temp-rInput)/(fs->GetParameter(1)));
                        }
              }


 TF1* fn = new TF1("gaus","gaus",-3,3);
 pullnormDistribution->Fit(fn,"R");

 TCanvas* c4 = new TCanvas ("c4","c4");
 c4->cd();
 c4->SetGrid();

 pullnormDistribution->GetYaxis()->SetTitle("N");
 pullnormDistribution->GetXaxis()->SetTitle("pull norm");
 pullnormDistribution->SetMarkerStyle(20);
 pullnormDistribution->SetMarkerSize(0.9);
 fn->SetLineColor(4);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(1111);

 pullnormDistribution->Draw("pe");
 fn->Draw("same");
 c4->Print("pull_norm.png",".png");
 c4->Close();

return(0);
}
