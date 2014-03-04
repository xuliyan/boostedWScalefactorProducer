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

TFile* file = TFile::Open("cards_em/higgsCombinehwwlvj_ggH600_em_2jet_10_00_unbin.MaxLikelihoodFit.mH600.-282267362.root");

TTree* tree_s = (TTree*)file->Get("limit");

 Double_t limit, limit_temp, rInput=10;
 Double_t limitErr, err_dw_temp, err_up_temp;
 Float_t quantileExpected;

 tree_s->SetBranchAddress("limit",&limit);
 tree_s->SetBranchAddress("limitErr",&limitErr);
 tree_s->SetBranchAddress("quantileExpected",&quantileExpected);

 TH1F* pullDistribution = new TH1F ("pull","mu",15,-1,1);
 TH1F* LimitDistribution = new TH1F ("limit2","limit2",30,-3,3);

for( int iEntry = 0; iEntry < tree_s->GetEntries(); iEntry ++){

        tree_s->GetEntry(iEntry);
        if(quantileExpected==0.5){ limit_temp = limit ; continue ; }
        if(quantileExpected>0.15 && quantileExpected<0.17 ) {err_dw_temp = limit ;  continue ; }
        if(quantileExpected>0.83 && quantileExpected<0.85) {err_up_temp = limit ;  continue ; }
        if(quantileExpected==-1){
             pullDistribution->Fill((limit_temp-rInput)/(err_up_temp-err_dw_temp));
             LimitDistribution->Fill(limit_temp);
                        }
              }

 TF1* f = new TF1("gaus","gaus",-1,1);
 pullDistribution->Fit(f,"R");

 TCanvas* c1 = new TCanvas ("c1","c1");
 c1->cd();
 c1->SetGrid();

 pullDistribution->GetYaxis()->SetTitle("N");
 pullDistribution->GetXaxis()->SetTitle("pull");
 pullDistribution->SetMarkerStyle(20);
 pullDistribution->SetMarkerSize(0.9);
 f->SetLineColor(4);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(1111);

 pullDistribution->Draw("pe");
 f->Draw("same");
 c1->Print("600TL.png",".png");
 c1->Close();

return(0);
}
