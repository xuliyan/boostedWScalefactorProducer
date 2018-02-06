import ROOT,sys
from WTopScalefactorProducer.Fitter.tdrstyle import *
from WTopScalefactorProducer.Fitter.CMS_lumi import *
from WTopScalefactorProducer.Skimmer.genEv import getGenEvents
setTDRStyle()
from time import sleep

ROOT.gROOT.SetBatch(True)

CMS_lumi.lumi_13TeV = "41.4 fb^{-1}(2017)"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

lumi = 41368 #RunD+B:  9037773962.937 ub #Run D only: 4235371340.858 ub
dir = "/scratch/thaarres/NANO_06feb/"
cutL = "(MET_sumEt[0]>40&&FatJet_pt[0]>200&&abs(FatJet_eta[0])<2.4&&abs(dr_LepJet)>1.5708&&abs(dphi_MetJet)>2.&&abs(dphi_WJet)>2.&&W_pt[0]>200&&Muon_pt[0]>65)"
vars = ["Muon_eta[0]","Muon_pt[0]","FatJet_pt[0]","FatJet_eta[0]","FatJet_tau1[0]","FatJet_tau2[0]","FatJet_tau3[0]","FatJet_mass[0]","FatJet_msoftdrop[0]","FatJet_tau21","FatJet_tau32","FatJet_softDrop_mass","MET_sumEt[0]"] #Pileup_nPU
bkgs = []
QCDs  = ["QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8.root"]
STs   = ["ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.root" , "ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root" , "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root","ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root","ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.root"]
STs   = ["ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.root" ,"ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root" , "ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.root"]
VVs   = ["WW_TuneCP5_13TeV-pythia8.root" , "WZ_TuneCP5_13TeV-pythia8.root", "ZZ_TuneCP5_13TeV-pythia8.root"]
WJs   = ["WJetsToLNu_HT-70To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root" , "WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"]
TTs   = ["TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root"]

bkgs.append(STs)
bkgs.append(VVs)
bkgs.append(QCDs)
bkgs.append(WJs)
bkgs.append(TTs)
legs=["Single Top (2017)","VV (2017)","QCD (2016)","W+jets (2016)", "TT (2016)"]
genEvs = getGenEvents()

datas = ["SingleMuon_Run2017B-17Nov2017-v1.root","SingleMuon_Run2017C-17Nov2017-v1.root","SingleMuon_Run2017D-17Nov2017-v1.root","SingleMuon_Run2017E-17Nov2017-v1.root","SingleMuon_Run2017F-17Nov2017-v1.root"]
fillcolor = [432,600,619,632,417]

def setTreeWeight(filename):
  print "For = ", filename
  file = ROOT.TFile(dir+filename, 'UPDATE')
  treeE = file.Get('Events')
  event = treeE.GetEntry(0)
  if filename.find("WJetsToLNu_HT-70To100")==-1: xSec = treeE.crossection
  else: xSec = 1270.
  treeR = file.Get('Runs')
  event = treeR.GetEntry(0)
  # neV = treeR.genEventCount
  neV = genEvs[filename]
  weight = xSec/neV
  if filename.find("QCD")!=-1: weight = weight*0.77
  print "XS  = " , xSec; print "N   =" ,neV; print "Rescaling tree by = ", weight
  treeE.SetWeight(weight)
  treeE.AutoSave()
  print "Tree weight is now: " ,treeE.GetWeight()
  return 
  
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
  	
def drawTH1(id,tree,var,cuts,bins,min,max,fillcolor,titlex = "",units = "",drawStyle = "HIST",lumi=lumi):
	h = ROOT.TH1D("tmpTH1","",bins,min,max)
	h.Sumw2()
	h.SetFillColor(fillcolor)
	if units=="":
	    h.GetXaxis().SetTitle(titlex)
	else:
	    h.GetXaxis().SetTitle(titlex+ " ["+units+"]")
	corrections = "1."
	if id.find("data")==-1: corrections = "genWeight"
	tree.Draw(var+">>tmpTH1","("+cuts+")*%s"%lumi+"*("+corrections+")","goff")
	return h

def doCP(postfix=""):
  for var in vars:
    name = var
    canvas,legend,pave = getCanvas()
    minx = 200.
    maxx = 2000.
    bins = 36
    unit = "GeV"
    if var.find("MET")!=-1: minx=40.; maxx=3000.; bins=60; 
    elif var.find("softdrop")!=-1 or var.find("mass")!=-1: minx=40.; maxx=130.; bins=18; 
    elif var.find("eta")!=-1 or var.find("phi")!=-1: minx=-2.5; maxx=2.5; bins=25; unit = "";
    elif var.find("tau")!=-1: minx=0.; maxx=1.0; bins=20; unit = "";
    elif var.find("Muon_pt")!=-1: minx=50.; maxx=1000.0; bins=100; unit = "GeV";
    treeD = ROOT.TChain("Events")
    for file in datas:
      print "Using file: ", ROOT.TString(dir+file)
      fileIn_name = ROOT.TString(dir+file)  
      treeD.Add(fileIn_name.Data())
    cutsData ='*'.join([cutL,"(passedMETfilters==1)"])
    print cutsData
    datahist = drawTH1("data",treeD,var,cutsData,bins,minx,maxx,1,var.replace("_", " ").replace("[0]", "").replace("FatJet", "AK8 jet"),unit,"HIST","1")
    datahist.SetName("data")
    legend.AddEntry(datahist,"Data (2017)","LEP")
    hists=[]
    stack = ROOT.THStack("stack","")
    ttint			= 0
    totalMinoInt 	= 0
    for i,bg in enumerate(bkgs):
      tree = ROOT.TChain("Events")
      name = bg[0]
      print "Name is: ", name
      for file in bg:
        print "Using file: ", ROOT.TString(dir+file)
        setTreeWeight(file)
        fileIn_name = ROOT.TString(dir+file)  
        tree.Add(fileIn_name.Data())
      hist = drawTH1(str(i),tree,var,cutL,bins,minx,maxx,fillcolor[i],var.replace("_", " ").replace("[0]", "").replace("FatJet", "AK8 jet"),unit,"HIST",lumi)
      legend.AddEntry(hist,legs[i],"F")
      hist.SetFillColor(fillcolor[i])
      if name.find("TT")!=-1:
        ttint = hist.Integral()
        hist.Scale(0.807476390014)
      else: totalMinoInt += hist.Integral()
      stack.Add(hist)
      hists.append(hist)
    
    scale = datahist.Integral()-totalMinoInt
    print "DATA/MC" ,scale/ttint
    canvas.cd()
    datahist.GetYaxis().SetRangeUser(0, datahist.GetMaximum()*1.6);
    if var.find("Jet_pt")!=-1:
    	datahist.GetYaxis().SetRangeUser(0.1, datahist.GetMaximum()*1000);
    	canvas.SetLogy()
    datahist.Draw("ME")
    stack.Draw("HISTsame")
    datahist.Draw("MEsame")
    legend.Draw("SAME")
    CMS_lumi(canvas, iPeriod, iPos)
    canvas.Update()
    canvas.SaveAs(var.replace("[0]","").replace("fabs(","").replace(")","")+postfix+".png")
    # sleep(1000)
    
    # ["jetAK8_softDrop_mass","jetAK8_softDrop_mass_unCorr","jetAK8_gen_softDrop_mass","jetAK8_tau1",  "jetAK8_tau2",  "jetAK8_tau3",  "jetAK8_tau21","jetAK8_tau32",
    # "jetAK8_highestSubJetCSV", "jetAK8_pt","jetAK8_eta","jetAK8_gen_pt","jetAK8_csv","lep_pt","Wlep_pt","lep_eta","lep_phi"]


	
	
if __name__ == "__main__":
	doCP()
