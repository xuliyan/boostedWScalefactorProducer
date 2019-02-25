import ROOT,sys
from WTopScalefactorProducer.Fitter.tdrstyle import *
from WTopScalefactorProducer.Fitter.CMS_lumi import *
from WTopScalefactorProducer.Skimmer.getGenEv import getGenEv
from WTopScalefactorProducer.Skimmer.xsec import getXsec, getNev, getSF
setTDRStyle()
from time import sleep

ROOT.gROOT.SetBatch(True)
lumi = 59970. #*0.024

CMS_lumi.lumi_13TeV = "%.1f fb^{-1} (2018)" % lumi
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV"
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4 #iPeriod = 0 for simulation-only plots

dir = "/scratch/zucchett/Ntuple/WSF/"
plotdir = "plots/"
#cut = "(SelectedJet_tau21_ddt_retune<0.57)"
#cut = "(abs(dr_LepJet)>1.5708&&abs(dphi_MetJet)>2.&&abs(dphi_WJet)>2&&Wlep_type==0&&SelectedJet_tau21<0.40)"
#cut = "(1==1)"
cut = "(HLT_Mu50&&nMuon>0&&Muon_pt[0]>55.&&abs(Muon_eta[0])<2.4&&Muon_pfRelIso03_all[0]<0.01&&Muon_highPtId>1&&FatJet_jetId>1&&abs(dr_LepJet)>1.5708&&abs(dphi_MetJet)>1.5708&&maxAK4CSV>0.8484)" #&&abs(dphi_WJet)>2&&Wlep_type==0
vars = ["SelectedJet_softDrop_mass","SelectedJet_tau21", "SelectedJet_tau21_ddt", "SelectedJet_tau21_ddt_retune","FatJet_pt[0]","FatJet_eta[0]","FatJet_phi[0]","FatJet_tau1[0]","FatJet_tau2[0]","FatJet_tau3[0]","FatJet_mass[0]","FatJet_msoftdrop[0]","MET_sumEt[0]","Muon_pt[0]","Muon_eta[0]","Muon_phi[0]","Muon_pfRelIso03_chg[0]","maxAK4CSV","nFatJet", "nJet", "nMuon","W_pt","fabs(dphi_WJet)","fabs(dphi_MetJet)","fabs(dphi_LepJet)","dr_LepJet"] #Pileup_nPU
vars = ["PV_npvs", "MET_pt"]
#vars = ["SelectedJet_softDrop_mass","SelectedJet_tau21", "SelectedJet_tau21_ddt", "SelectedJet_tau21_ddt_retune"]

#Data infile
datas   = ["SingleMuon.root"]

#MC infiles
bkgs = []
STs   = ["ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root", "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root", "ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root", "ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.root"]
TTs   = ["TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.root", "TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.root"]
WJs   = ["WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8.root", "WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.root"]
VVs   = ["WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", "WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", "ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8.root"] #"WWTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8.root", 
bkgs.append(STs)
bkgs.append(VVs)
bkgs.append(WJs)
bkgs.append(TTs)

#For drawing
legs=["Single Top (2017)","VV (2017)","W+jets (2017)", "TT (2017)"]
fillcolor = [432,600,632,417]


#def setTreeWeight(filename):
#  print "For = ", filename
#  file = ROOT.TFile(dir+filename, 'UPDATE')
#  treeE = file.Get('Events')
#  event = treeE.GetEntry(0)
#  xSec = treeE.crossection
#  genEv = getGenEv(file.GetName())
#  LEQ = float(xSec*lumi/genEv)
#  print "xs  = " , xSec; print "N   =" ,genEv; print "Rescaling tree by = ", LEQ
#  treeE.SetWeight(LEQ)
#  treeE.AutoSave()
#  print "Tree weight is now: " ,treeE.GetWeight()

def setTreeWeight(filename):
  print "For = ", filename
  file = ROOT.TFile(dir+filename, 'UPDATE')
  treeE = file.Get('Events')
  event = treeE.GetEntry(0)
  xSec = getXsec(filename)
  xSF  = getSF(filename)
  genEv = getNev(filename)
  LEQ = float(xSec)*float(xSF)*lumi/float(genEv)
  treeE.SetWeight(LEQ)
  treeE.AutoSave()
  print "Tree weight is now: " ,treeE.GetWeight()
  
def getCanvas():
  H_ref = 600
  W_ref = 800
  W = W_ref
  H  = H_ref
  T = 0.08*H_ref
  B = 0.12*H_ref 
  L = 0.12*W_ref
  R = 0.04*W_ref
  canvas = ROOT.TCanvas("c1","c1",50,50,W,H)
  canvas.SetFillColor(0)
  canvas.SetBorderMode(0)
  canvas.SetFrameFillStyle(0)
  canvas.SetFrameBorderMode(0)
  canvas.SetLeftMargin( L/W )
  canvas.SetRightMargin( R/W )
  canvas.SetTopMargin( T/H )
  canvas.SetBottomMargin( B/H )
  canvas.SetTickx(0)
  canvas.SetTicky(0)
  ROOT.gStyle.SetOptStat(0)
  ROOT.gStyle.SetOptTitle(0)
  legend = ROOT.TLegend(0.62,0.7,0.92,0.9,"","brNDC")
  legend.SetBorderSize(0)
  legend.SetLineColor(1)
  legend.SetLineStyle(1)
  legend.SetLineWidth(1)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetTextFont(42)
  addInfo = ROOT.TPaveText(0.73010112,0.2566292,0.8202143,0.5523546,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  return canvas, legend, addInfo
  	
def drawTH1(id,tree,var,cuts,bins,min,max,fillcolor,titlex = "",units = "",drawStyle = "HIST",lumi="%s"%lumi):
  h = ROOT.TH1D("tmpTH1","",bins,min,max)
  h.Sumw2()
  h.SetFillColor(fillcolor)
  if units=="": h.GetXaxis().SetTitle(titlex)
  else: h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
  tree.Draw(var+">>tmpTH1","("+cuts+")","goff")
#  corrString='1'
#  if id.find("data")==-1:
#    corrString = corrString+"*(genWeight)"
#  else:
#    lumi = "1"
#  tree.Draw(var+">>tmpTH1","("+cuts+")*"+lumi+"*("+corrString+")","goff")
#  if id.find("data")==-1:
#      tree.Draw(var+">>tmpTH1","("+cuts+")","goff")
#  else:
#      tree.Draw(var+">>tmpTH1","("+cuts+")*"+lumi+"*(genWeight)","goff")
  return h

def doCP(cutL,postfix=""):
  for var in vars:
    name = var
    canvas,legend,pave = getCanvas()
    minx = 200.
    maxx = 2000.
    bins = 36
    unit = "GeV"
    if var.find("MET")!=-1: minx=40.; maxx=3000.; bins=60; 
    if var.find("softdrop")!=-1 or var.find("mass")!=-1: minx=0.; maxx=130.; bins=16; 
    if var.find("eta")!=-1: minx=-3; maxx=3; bins=25; unit = "";
    if var.find("phi")!=-1: minx=-3.15; maxx=3.15; bins=50; unit = "";
    if var.find("tau")!=-1: minx=0.; maxx=1.0; bins=20; unit = "";
    if var.find("ddt")!=-1: minx=-0.2; maxx=1.2; bins=20; unit = "";
    if var.find("Muon_pt")!=-1: minx=50.; maxx=1000.0; bins=100; unit = "GeV";
    if var.find("RelIso")!=-1: minx=0.; maxx=0.15; bins=100; unit = "";
    if var.find("CSV")!=-1: minx=0.; maxx=1.; bins=100; unit = "";
    if var.find("nJet")!=-1: minx=-0.5; maxx=19.5; bins=20; unit = "";
    if var.find("nFatJet")!=-1: minx=-0.5; maxx=9.5; bins=10; unit = "";
    if var.find("nMuon")!=-1: minx=-0.5; maxx=4.5; bins=5; unit = "";
    if var.find("W_pt")!=-1: minx=0.; maxx=1000.; bins=100; unit = "";
    if var.find("dphi_")!=-1: minx=0; maxx=3.15; bins=30; unit = "";
    if var.find("dr_")!=-1: minx=0; maxx=5; bins=25; unit = "";
    if var.find("npvs")!=-1: minx=0; maxx=80; bins=80; unit = "";
    if var.find("MET_pt")!=-1: minx=0; maxx=500; bins=25; unit = "";
    treeD = ROOT.TChain("Events")
    for file in datas:
      print "Using file: ", ROOT.TString(dir+file)
      fileIn_name = ROOT.TString(dir+file)  
      treeD.Add(fileIn_name.Data())
    cutT = cutL if not var in cutL else ()
    
    cutsData ='*'.join([cutL]) #,"(passedMETfilters==1)"])
    cutsData ='*'.join([cutL])
    datahist = drawTH1("data",treeD,var,cutsData,bins,minx,maxx,1,var.replace("_", " ").replace("[0]", "").replace("FatJet", "AK8 jet"),unit,"HIST","1")
    datahist.SetName("data")
    legend.AddEntry(datahist,"Data (2018)","LEP")
    hists=[]
    stack = ROOT.THStack("stack","")
    dataint		= datahist.Integral(0, datahist.GetNbinsX()+1)
    backint 	= [None]*len(bkgs)
    for i,bg in enumerate(bkgs):
      tree = ROOT.TChain("Events")
      name = bg[0]
      print "Name is: ", name
      for file in bg:
        print "Using file: ", ROOT.TString(dir+file)
        setTreeWeight(file)
        fileIn_name = ROOT.TString(dir+file)  
        tree.Add(fileIn_name.Data())
      hist = drawTH1(str(i),tree,var,cutL,bins,minx,maxx,fillcolor[i],var.replace("_", " ").replace("[0]", "").replace("FatJet", "AK8 jet"),unit,"HIST")
      legend.AddEntry(hist,legs[i],"F")
      hist.SetFillColor(fillcolor[i])
      backint[i] = hist.Integral(0, hist.GetNbinsX()+1)
#      if name.find("TT")!=-1:
#        ttint = hist.Integral()
#        scale = (datahist.Integral()-totalMinoInt)/ttint
#        hist.Scale(scale)
#      else: totalMinoInt += hist.Integral()
      stack.Add(hist)
      hists.append(hist)
    
    
#    print "DATA/MC" ,scale
    canvas.cd()
    datahist.GetYaxis().SetRangeUser(0, datahist.GetMaximum()*1.6);
    if var in ["FatJet_pt[0]","Muon_pt[0]","Muon_pfRelIso03_chg[0]","Muon_pfRelIso03_all[0]"]:
    	datahist.GetYaxis().SetRangeUser(0.1, datahist.GetMaximum()*1000);
    	canvas.SetLogy()
    datahist.Draw("ME")
    stack.Draw("HIST SAME")
    datahist.Draw("ME same")
    legend.Draw("SAME")
    CMS_lumi(canvas, iPeriod, iPos)
    canvas.Update()
    canvas.SaveAs(plotdir+var.replace("[0]","").replace("abs(","").replace("(","").replace(")","")+postfix+".png")
    canvas.SaveAs(plotdir+var.replace("[0]","").replace("abs(","").replace("(","").replace(")","")+postfix+".pdf")
    
    print "--- Summary ----"
    print "Data (2018)", "\t", dataint
    print "-"*20
    for i,bg in enumerate(bkgs):
      print legs[i], "\t", backint[i]
    print "-"*20
    print "Ratio data/bkg", "\t", dataint/sum(backint)
    print "\n"
    b = 2
    print "Scale factor for", legs[i], (dataint - sum([x for i, x in enumerate(backint) if i!=b] )) / backint[b]
    
    sleep(10)

if __name__ == "__main__":
	doCP(cut)
