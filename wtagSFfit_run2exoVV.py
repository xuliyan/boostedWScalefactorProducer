from optparse import OptionParser
import ROOT, sys
import ROOT as rt
import time
import math

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="em")
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--HP', action="store", type="float",dest="tau2tau1cutHP",default=0.60)
parser.add_option('--LP', action="store", type="float",dest="tau2tau1cutLP",default=0.75)
parser.add_option('--csvMax', action="store",type="float",dest="csvMax",default=100)
parser.add_option('--sample', action="store",type="string",dest="sample",default="")
parser.add_option('--fitTT', action='store_true', dest='fitTT', default=False, help='Only do ttbar fits')
parser.add_option('--lowMass', action='store_true', dest='doLowMass', default=False, help='Use low-mass samples')
parser.add_option('--76X',dest="use76X", default=True, action="store_true", help="Use 76X samples")
parser.add_option('--usePuppiSD',dest="usePuppiSD", default=False, action="store_true", help="Use PUPPI+softdrop")
parser.add_option('--useDDT',dest="useDDT", default=False, action="store_true", help="Use DDT tagger")

(options, args) = parser.parse_args()

ROOT.gSystem.Load(".//PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(".//PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(".//PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(".//PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(".//PDFs/MakePdf_cxx.so")
ROOT.gSystem.Load(".//BiasStudy/BiasUtils_cxx.so")
ROOT.gSystem.Load(".//FitUtils/FitUtils_cxx.so")

from ROOT import RooWorkspace, RooAbsPdf, setTDRStyle, ScaleFactorTTbarControlSampleFit
from ROOT import *

gInterpreter.GenerateDictionary("std::map<std::string,std::string>", "map;string;string")
# gInterpreter.GenerateDictionary("std::vector<std::string>", "vector;string")

RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

if options.noX:
  gROOT.SetBatch(True)

### Fit e+mu together
def getSF():

    print "Getting W-tagging SF for cut " ,options.tau2tau1cutHP
    if options.useDDT: options.usePuppiSD = True
    boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous()

def control_sample(channel="em"):

    print "control_sample"
    
class doFit_wj_and_wlvj_simultaneous:
  
    def __init__(self):
      
      self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_")                           # create workspace
      self.boostedW_fitter_em = doFit_wj_and_wlvj("em", options.sample, 40, 130, self.workspace4fit_) # Define all shapes to be used for Mj, define regions (SB,signal) and input files. 
      self.boostedW_fitter_em.get_datasets_fit_minor_bkg()                                            # Loop over intrees to create datasets om Mj and fit the single MCs.
     
      self.workspace4fit_.Print()
      onlyMCfits = False                                                                              # To be removed. Flag for when testing fits to minor backgrounds only (skip simoultaneous fit) 

      if not options.fitTT and not onlyMCfits:
        
        postfix=""
        if options.usePuppiSD: postfix = "_PuppiSD"
        title = "Pruned jet mass (GeV)"
        if options.usePuppiSD:  title = "PUPPI softdrop jet mass (GeV)"
        
        self.workspace4fit_.data("rdataset_data_em_mj").Print() 
        self.workspace4fit_.data("rdataset_data_failtau2tau1cut_em_mj").Print()
        self.workspace4fit_.data("rdataset_TotalMC_em_mj").Print()
        self.workspace4fit_.data("rdataset_TotalMC_failtau2tau1cut_em_mj").Print() 
 
        #Defining categories
        sample_type = RooCategory("sample_type","sample_type")
        sample_type.defineType("em_pass")
        sample_type.defineType("em_fail")

        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)
 
        #Importing datasets
        rdataset_data_em_mj      = self.workspace4fit_.data("rdataset_data_em_mj")
        rdataset_data_em_mj_fail = self.workspace4fit_.data("rdataset_data_failtau2tau1cut_em_mj")
 
        rrv_mass_j   = self.workspace4fit_.var("rrv_mass_j")

        #Combined dataset
        combData_data = RooDataSet("combData_data","combData_data",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_data_em_mj),RooFit.Import("em_fail",rdataset_data_em_mj_fail) )


        rdataset_TotalMC_em_mj      = self.workspace4fit_.data("rdataset_TotalMC_em_mj")
        rdataset_TotalMC_em_mj_fail = self.workspace4fit_.data("rdataset_TotalMC_failtau2tau1cut_em_mj")
      

        combData_TotalMC = RooDataSet("combData_TotalMC","combData_TotalMC",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_TotalMC_em_mj),RooFit.Import("em_fail",rdataset_TotalMC_em_mj_fail) )
        combData_TotalMC.Print()

        # Import pdf from single fits and define the simultaneous total pdf
        model_data_em      = self.workspace4fit_.pdf("model_data_em")
        model_data_fail_em = self.workspace4fit_.pdf("model_data_failtau2tau1cut_em")
 
        simPdf_data = RooSimultaneous("simPdf_data_em","simPdf_data_em",sample_type)
        simPdf_data.addPdf(model_data_em,"em_pass")
        simPdf_data.addPdf(model_data_fail_em,"em_fail")

       
        constrainslist_data_em = ROOT.std.vector(ROOT.std.string)()
        for i in range(self.boostedW_fitter_em.constrainslist_data.size()):
            constrainslist_data_em.push_back(self.boostedW_fitter_em.constrainslist_data.at(i))
            print self.boostedW_fitter_em.constrainslist_data.at(i)
      

        pdfconstrainslist_data_em = RooArgSet("pdfconstrainslist_data_em")
        for i in range(constrainslist_data_em.size()):
            # self.workspace4fit_.pdf(constrainslist_data_em.at(i)).Print()
            pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em.at(i)) )
     
        # Perform simoultaneous fit to data
        rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em),RooFit.Verbose(kFALSE))
        # rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em),RooFit.Verbose(kFALSE))
        

        frame = rrv_mass_j.frame()
        rdataset_data_em_mj_fail.plotOn(frame,rt.RooFit.DataError(rt.RooAbsData.Poisson),rt.RooFit.Name("rdataset_data_failtau2tau1cut_em_mj"))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error"),RooFit.FillColor(kRed-7),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,rt.RooFit.LineColor(rt.kRed+1),rt.RooFit.Name("model_data_failtau2tau1cut_em"))
        chi2_fail_data = frame.chiSquare("model_data_failtau2tau1cut_em", "rdataset_data_failtau2tau1cut_em_mj")

        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Gauss1"),RooFit.Components("gaus1*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Gauss2"),RooFit.Components("gaus2*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error ErfExp"),RooFit.Components("*ErfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error erfExp"),RooFit.Components("*erfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Cheb"),RooFit.Components("cheb*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Exp*"),RooFit.Components("Exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error exp*"),RooFit.Components("exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "Gaussian 1" ),RooFit.Components("gaus1*"),RooFit.LineStyle(kDashed),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "Gaussian 2" ),RooFit.Components("gaus2*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "ErfExp comp." ),RooFit.Components("*ErfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "erfExp comp." ),RooFit.Components("*erfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "Chebychev" ),RooFit.Components("cheb*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "Exp comp." ),RooFit.Components("Exp*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_fail_em.plotOn(frame,RooFit.Name( "exp comp." ),RooFit.Components("exp*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        rdataset_data_em_mj_fail.plotOn(frame,rt.RooFit.DataError(rt.RooAbsData.Poisson),rt.RooFit.Name("rdataset_data_failtau2tau1cut_em_mj"))



        c1 =rt.TCanvas("c1","",800,800)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetTitleOffset(0.90)
        # frame.GetYaxis().SetLabelSize(0.09)
        frame.SetName("mjjFit")
        frame.GetYaxis().SetTitle("Events")
        frame.GetXaxis().SetTitle(title)
        frame.Draw()

        legend = rt.TLegend(0.6010112,0.7183362,0.8202143,0.919833)
        legend.SetTextSize(0.032)
        legend.SetLineColor(0)
        legend.SetShadowColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetMargin(0.35)
        legend.AddEntry(frame.findObject("rdataset_data_failtau2tau1cut_em_mj"),"CMS data","lpe")
        legend.AddEntry(frame.findObject("model_data_failtau2tau1cut_em"),"Sim. fit","l")
        if frame.findObject("Gaussian 1"):
          legend.AddEntry(frame.findObject("Gaussian 1"),frame.findObject("Gaussian 1").GetName(),"l")
        if frame.findObject("Gaussian 2"):
          legend.AddEntry(frame.findObject("Gaussian 2"),frame.findObject("Gaussian 2").GetName(),"l")
        if frame.findObject("ErfExp comp."):
          legend.AddEntry(frame.findObject("ErfExp comp."),frame.findObject("ErfExp comp.").GetName(),"l")
        if frame.findObject("erfExp comp."):
          legend.AddEntry(frame.findObject("erfExp comp."),frame.findObject("erfExp comp.").GetName(),"l")
        if frame.findObject("Chebychev"):
          legend.AddEntry(frame.findObject("Chebychev"),frame.findObject("Chebychev").GetName(),"l")
        if frame.findObject("Exp comp."):
          legend.AddEntry(frame.findObject("Exp comp."),frame.findObject("Exp comp.").GetName(),"l")
        if frame.findObject("exp comp."):
          legend.AddEntry(frame.findObject("exp comp."),frame.findObject("exp comp.").GetName(),"l")
        legend.Draw("same")

        addInfo = rt.TPaveText(0.2510112,0.2066292,0.4202143,0.3523546,"NDC")
        addInfo.AddText("#chi^{2}/nDOF = %.3f"%chi2_fail_data)
        addInfo.SetFillColor(0)
        addInfo.SetLineColor(0)
        addInfo.SetFillStyle(0)
        addInfo.SetBorderSize(0)
        addInfo.SetTextFont(42)
        addInfo.SetTextSize(0.040)
        addInfo.SetTextAlign(12)
        addInfo.Draw()
        c1.Update()
        c1.SaveAs("plots/DATA-fail_em_HP%.2f%s%s.pdf"%(options.tau2tau1cutHP,options.sample,postfix))


        frame2 = rrv_mass_j.frame()
        rdataset_data_em_mj.plotOn(frame2,rt.RooFit.DataError(rt.RooAbsData.Poisson),rt.RooFit.Name("rdataset_data_em_mj"))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error"),RooFit.FillColor(kRed-7),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,rt.RooFit.LineColor(rt.kRed+1),rt.RooFit.Name("model_data_em"))
        chi2_pass_data = frame2.chiSquare("model_data_em", "rdataset_data_em_mj")
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Gauss1"),RooFit.Components("gaus1*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Gauss2"),RooFit.Components("gaus2*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error ErfExp"),RooFit.Components("*ErfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error erfExp"),RooFit.Components("*erfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Cheb"),RooFit.Components("cheb*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error Exp*"),RooFit.Components("Exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.VisualizeError(rfresult_data,1), RooFit.Name("Fit error exp*"),RooFit.Components("exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "Gaussian 1" ),RooFit.Components("gaus1*"),RooFit.LineStyle(kDashed),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "Gaussian 2" ),RooFit.Components("gaus2*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "ErfExp comp." ),RooFit.Components("*ErfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "erfExp comp." ),RooFit.Components("erfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "Chebychev" ),RooFit.Components("cheb*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "Exp comp." ),RooFit.Components("Exp*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_data_em.plotOn(frame2,RooFit.Name( "exp comp." ),RooFit.Components("exp_*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        rdataset_data_em_mj.plotOn(frame2,rt.RooFit.DataError(rt.RooAbsData.Poisson),rt.RooFit.Name("rdataset_data_em_mj"))


        c2 =rt.TCanvas("c2","",800,800)
        frame2.GetYaxis().SetTitleSize(0.05)
        frame2.GetYaxis().SetTitleOffset(0.90)
        # frame.GetYaxis().SetLabelSize(0.09)
        frame2.SetName("mjjFit")
        frame2.GetYaxis().SetTitle("Events")
        frame2.GetXaxis().SetTitle(title)
        frame2.Draw()



        legend = rt.TLegend(0.6010112,0.7183362,0.8202143,0.919833)
        legend.SetTextSize(0.032)
        legend.SetLineColor(0)
        legend.SetShadowColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetMargin(0.35)
        legend.AddEntry(frame2.findObject("rdataset_data_em_mj"),"CMS data","lpe")
        legend.AddEntry(frame2.findObject("model_data_em"),"Sim. fit","l")
        if frame2.findObject("Gaussian 1"):
          legend.AddEntry(frame2.findObject("Gaussian 1"),frame2.findObject("Gaussian 1").GetName(),"l")
        if frame2.findObject("Gaussian 2"):
          legend.AddEntry(frame2.findObject("Gaussian 2"),frame2.findObject("Gaussian 2").GetName(),"l")
        if frame2.findObject("ErfExp comp."):
          legend.AddEntry(frame2.findObject("ErfExp comp."),frame2.findObject("ErfExp comp.").GetName(),"l")
        if frame2.findObject("erfExp comp."):
          legend.AddEntry(frame2.findObject("erfExp comp."),frame2.findObject("erfExp comp.").GetName(),"l")
        if frame2.findObject("Exp comp."):
          legend.AddEntry(frame2.findObject("Exp comp."),frame2.findObject("Exp comp.").GetName(),"l")
        if frame2.findObject("exp comp."):
          legend.AddEntry(frame2.findObject("exp comp."),frame2.findObject("exp comp.").GetName(),"l")
        legend.Draw("same")
        addInfo = rt.TPaveText(0.2510112,0.2066292,0.4202143,0.3523546,"NDC")
        addInfo.AddText("#chi^{2}/nDOF = %.3f"%chi2_pass_data)
        addInfo.SetFillColor(0)
        addInfo.SetLineColor(0)
        addInfo.SetFillStyle(0)
        addInfo.SetBorderSize(0)
        addInfo.SetTextFont(42)
        addInfo.SetTextSize(0.040)
        addInfo.SetTextAlign(12)
        addInfo.Draw()
        c2.Update()
        c2.SaveAs("plots/DATA-pass_em_HP%.2f%s%s.pdf"%(options.tau2tau1cutHP,options.sample,postfix))



        print "FIT parameters (DATA) :"
        print ""
        print "CHI2 PASS = %.3f    CHI2 FAIL = %.3f" %(chi2_pass_data,chi2_fail_data)
        print ""
        print rfresult_data.Print()
        print ""

  
 

        # fit TotalMC --> define the simultaneous total pdf
        model_TotalMC_em      = self.workspace4fit_.pdf("model_TotalMC_em")
        model_TotalMC_fail_em = self.workspace4fit_.pdf("model_TotalMC_failtau2tau1cut_em")

        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em","simPdf_TotalMC_em",sample_type)
        simPdf_TotalMC.addPdf(model_TotalMC_em,"em_pass")
        simPdf_TotalMC.addPdf(model_TotalMC_fail_em,"em_fail")


        constrainslist_TotalMC_em = ROOT.std.vector(ROOT.std.string)()
        for i in range(self.boostedW_fitter_em.constrainslist_mc.size()):
            constrainslist_TotalMC_em.push_back(self.boostedW_fitter_em.constrainslist_mc.at(i))


        pdfconstrainslist_TotalMC_em = RooArgSet("pdfconstrainslist_TotalMC_em")
        for i in range(constrainslist_TotalMC_em.size()):
            # self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]).Print()
            pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) )

        # Perform simoultaneous fit to MC
        rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em),RooFit.SumW2Error(kTRUE),RooFit.Verbose(kFALSE))
        # rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em),RooFit.SumW2Error(kTRUE),RooFit.Verbose(kFALSE))


        frame3 = rrv_mass_j.frame()

        rdataset_TotalMC_em_mj_fail.plotOn(frame3,rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.Name("rdataset_TotalMC_failtau2tau1cut_em_mj"))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error"),RooFit.FillColor(kRed-7),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,rt.RooFit.LineColor(rt.kRed+1),rt.RooFit.Name("model_TotalMC_failtau2tau1cut_em"))
        chi2_fail_mc = frame3.chiSquare("model_TotalMC_failtau2tau1cut_em", "rdataset_TotalMC_failtau2tau1cut_em_mj")

        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Gauss1"),RooFit.Components("gaus1*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Gauss2"),RooFit.Components("gaus2*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error ErfExp"),RooFit.Components("*ErfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error erfExp"),RooFit.Components("*erfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Cheb"),RooFit.Components("cheb*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Exp*"),RooFit.Components("Exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error exp*"),RooFit.Components("exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "Gaussian 1" ),RooFit.Components("gaus1*"),RooFit.LineStyle(kDashed),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "Gaussian 2" ),RooFit.Components("gaus2*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "ErfExp comp." ),RooFit.Components("*ErfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "erfExp comp." ),RooFit.Components("erfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "Chebychev" ),RooFit.Components("cheb*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "Exp comp." ),RooFit.Components("Exp*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_fail_em.plotOn(frame3,RooFit.Name( "exp comp." ),RooFit.Components("exp_*"),RooFit.LineStyle(9),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        rdataset_TotalMC_em_mj_fail.plotOn(frame3,rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.Name("rdataset_TotalMC_failtau2tau1cut_em_mj"))

        c1 =rt.TCanvas("c1","",800,800)
        frame3.GetYaxis().SetTitleSize(0.05)
        frame3.GetYaxis().SetTitleOffset(0.90)
        # frame.GetYaxis().SetLabelSize(0.09)
        frame3.SetName("mjjFit")
        frame3.GetYaxis().SetTitle("Events")
        frame3.GetXaxis().SetTitle(title)
        frame3.Draw()

        legend = rt.TLegend(0.6510112,0.7183362,0.8202143,0.919833)
        legend.SetTextSize(0.038)
        legend.SetLineColor(0)
        legend.SetShadowColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetMargin(0.35)
        legend.AddEntry(frame3.findObject("rdataset_TotalMC_em_mj_fail"),"Tot. MC","lpe")
        legend.AddEntry(frame3.findObject("model_TotalMC_fail_em"),"Sim. fit","l")
        if frame3.findObject("Gaussian 1"):
          legend.AddEntry(frame3.findObject("Gaussian 1"),frame3.findObject("Gaussian 1").GetName(),"l")
        if frame3.findObject("Gaussian 2"):
          legend.AddEntry(frame3.findObject("Gaussian 2"),frame3.findObject("Gaussian 2").GetName(),"l")
        if frame3.findObject("ErfExp comp."):
          legend.AddEntry(frame3.findObject("ErfExp comp."),frame3.findObject("ErfExp comp.").GetName(),"l")
        if frame3.findObject("erfExp comp."):
          legend.AddEntry(frame3.findObject("erfExp comp."),frame3.findObject("erfExp comp.").GetName(),"l")
        if frame3.findObject("Chebychev"):
          legend.AddEntry(frame3.findObject("Chebychev"),frame3.findObject("Chebychev").GetName(),"l")
        if frame3.findObject("Exp comp."):
          legend.AddEntry(frame3.findObject("Exp comp."),frame3.findObject("Exp comp.").GetName(),"l")
        if frame3.findObject("exp comp."):
          legend.AddEntry(frame3.findObject("exp comp."),frame3.findObject("exp comp.").GetName(),"l")
        legend.Draw("same")

        addInfo = rt.TPaveText(0.2510112,0.2066292,0.4202143,0.3523546,"NDC")
        addInfo.AddText("#chi^{2}/nDOF = %.3f"%chi2_fail_mc)
        addInfo.SetFillColor(0)
        addInfo.SetLineColor(0)
        addInfo.SetFillStyle(0)
        addInfo.SetBorderSize(0)
        addInfo.SetTextFont(42)
        addInfo.SetTextSize(0.040)
        addInfo.SetTextAlign(12)
        addInfo.Draw()
        c1.Update()
        c1.SaveAs("plots/MC-fail_em_HP%.2f%s%s.pdf"%(options.tau2tau1cutHP,options.sample,postfix))


        frame4 = rrv_mass_j.frame()
        rdataset_TotalMC_em_mj.plotOn(frame4,rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.Name("rdataset_TotalMC_em_mj"))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error"),RooFit.FillColor(kRed-7),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,rt.RooFit.LineColor(rt.kRed+1),rt.RooFit.Name("model_TotalMC_em"))
        chi2_pass_mc = frame4.chiSquare("model_TotalMC_em", "rdataset_TotalMC_em_mj")
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Gauss1"),RooFit.Components("gaus1*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Gauss2"),RooFit.Components("gaus2*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error ErfExp"),RooFit.Components("*ErfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error erfExp"),RooFit.Components("*erfExp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Cheb"),RooFit.Components("cheb*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error Exp*"),RooFit.Components("Exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.VisualizeError(rfresult_TotalMC,1), RooFit.Name("Fit error exp*"),RooFit.Components("exp*"),RooFit.FillColor(kRed-9),RooFit.LineColor(kRed-7),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "Gaussian 1" ),RooFit.Components("gaus1*")  ,RooFit.LineStyle(kDashed),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "Gaussian 2" ),RooFit.Components("gaus2*")  ,RooFit.LineStyle(kSolid),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "ErfExp comp."     ),RooFit.Components("*ErfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "erfExp comp."     ),RooFit.Components("*erfExp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "Cheb. comp."     ),RooFit.Components("cheb*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "Exp. comp."     ),RooFit.Components("Exp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        model_TotalMC_em.plotOn(frame4,RooFit.Name( "exp. comp."     ),RooFit.Components("exp*"),RooFit.LineStyle(kDashDotted),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
        rdataset_TotalMC_em_mj.plotOn(frame4,rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.Name("rdataset_TotalMC_em_mj"))


        c2 =rt.TCanvas("c2","",800,800)
        frame4.GetYaxis().SetTitleSize(0.05)
        frame4.GetYaxis().SetTitleOffset(0.90)
        # frame.GetYaxis().SetLabelSize(0.09)
        frame4.SetName("mjjFit")
        frame4.GetYaxis().SetTitle("Events")
        frame4.GetXaxis().SetTitle(title)
        frame4.Draw()


        legend = rt.TLegend(0.6510112,0.7183362,0.8202143,0.919833)
        legend.SetTextSize(0.038)
        legend.SetLineColor(0)
        legend.SetShadowColor(0)
        legend.SetLineStyle(1)
        legend.SetLineWidth(1)
        legend.SetFillColor(0)
        legend.SetFillStyle(0)
        legend.SetMargin(0.35)
        legend.AddEntry(frame4.findObject("rdataset_TotalMC_em_mj"),"Tot. MC","lpe")
        legend.AddEntry(frame4.findObject("model_TotalMC_em"),"Sim. fit","l")
        if frame4.findObject("Gaussian 1"):
          legend.AddEntry(frame4.findObject("Gaussian 1"),frame4.findObject("Gaussian 1").GetName(),"l")
        if frame4.findObject("Gaussian 2"):
          legend.AddEntry(frame4.findObject("Gaussian 2"),frame4.findObject("Gaussian 2").GetName(),"l")
        if frame4.findObject("ErfExp comp."):
          legend.AddEntry(frame4.findObject("ErfExp comp."),frame4.findObject("ErfExp comp.").GetName(),"l")
        if frame4.findObject("erfExp comp."):
          legend.AddEntry(frame4.findObject("erfExp comp."),frame4.findObject("erfExp comp.").GetName(),"l")
        if frame4.findObject("Exp comp."):
          legend.AddEntry(frame4.findObject("Exp comp."),frame4.findObject("Exp comp.").GetName(),"l")
        if frame4.findObject("exp comp."):
          legend.AddEntry(frame4.findObject("exp comp."),frame4.findObject("exp comp.").GetName(),"l")
        legend.Draw("same")


        addInfo = rt.TPaveText(0.2510112,0.2066292,0.4202143,0.3523546,"NDC")
        addInfo.AddText("#chi^{2}/nDOF = %.3f"%chi2_pass_mc)
        addInfo.SetFillColor(0)
        addInfo.SetLineColor(0)
        addInfo.SetFillStyle(0)
        addInfo.SetBorderSize(0)
        addInfo.SetTextFont(42)
        addInfo.SetTextSize(0.040)
        addInfo.SetTextAlign(12)
        addInfo.Draw()
        c2.Update()
        c2.SaveAs("plots/MC-pass_em_HP%.2f%s%s.pdf"%(options.tau2tau1cutHP,options.sample,postfix))

        print "FIT Par. (MC) :"
        print ""
        print "CHI2 PASS= " ,chi2_pass_mc
        print "CHI2 FAIL= " ,chi2_fail_mc
        print ""
        print rfresult_TotalMC.Print()
        print ""


        # draw the plots
        DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_em.color_palet,"","em",self.boostedW_fitter_em.wtagger_label,self.boostedW_fitter_em.AK8_pt_min,self.boostedW_fitter_em.AK8_pt_max)


        ### Efficiency in data and MC
        rrv_eff_MC_em   = self.workspace4fit_.var("eff_ttbar_TotalMC_em_mj")
        rrv_mean_MC_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC_em_mj")
        rrv_sigma_MC_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC_em_mj")

        rrv_eff_data_em   = self.workspace4fit_.var("eff_ttbar_data_em_mj")
        rrv_mean_data_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data_em_mj")
        rrv_sigma_data_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data_em_mj")


        ## GET HP SCALEFACTOR AND UNCERTIANTIES
        pure_wtagger_sf_em             = rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()
        pure_wtagger_sf_em_err         = pure_wtagger_sf_em * ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal() )**2 + (rrv_eff_MC_em.getError()/rrv_eff_MC_em.getVal())**2 )**0.5

        pure_wtagger_mean_shift_em     = rrv_mean_data_em.getVal()-rrv_mean_MC_em.getVal()
        pure_wtagger_mean_shift_err_em = (rrv_mean_data_em.getError()**2 + rrv_mean_MC_em.getError()**2)**0.5

        pure_wtagger_sigma_enlarge_em     = rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()
        pure_wtagger_sigma_enlarge_err_em = ((rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 + (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_em

        mean_sf_error  = (rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal())   * ( (rrv_mean_data_em.getError()/rrv_mean_data_em.getVal())**2   +  (rrv_mean_MC_em.getError() /rrv_mean_MC_em.getVal())**2    )**0.5
        sigma_sf_error = (rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()) * ( (rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 +  (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2   )**0.5
        eff_sf_error   = (rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal())     * ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal())**2     +  (rrv_eff_MC_em.getError()  /rrv_eff_MC_em.getVal())**2     )**0.5




        ## GET EXTREME FAIL NUMBERS IN ORDER TO COMPUTE LP SF:

        rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj")
        rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj    = self.workspace4fit_.var("rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj")

        rrv_number_total_ttbar_TotalMC_em = self.workspace4fit_.var("rrv_number_ttbar_TotalMC_beforetau2tau1cut_em_mj")
        rrv_number_total_ttbar_data_em    = self.workspace4fit_.var("rrv_number_ttbar_data_beforetau2tau1cut_em_mj")

        eff_MC_em_extremefail   = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal()
        eff_data_em_extremefail = rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()    / rrv_number_total_ttbar_data_em.getVal()

        eff_MC_em_extremefail_error   = eff_MC_em_extremefail   * ( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() )**2 + (rrv_number_total_ttbar_TotalMC_em.getError()/rrv_number_total_ttbar_TotalMC_em.getVal())**2 )**0.5
        eff_data_em_extremefail_error = eff_data_em_extremefail * ( (rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()       )**2 + (rrv_number_total_ttbar_data_em.getError()   /rrv_number_total_ttbar_data_em.getVal()   )**2 )**0.5

        eff_SF_extremefail       = eff_data_em_extremefail/eff_MC_em_extremefail
        eff_SF_extremefail_error = eff_SF_extremefail * ( (eff_data_em_extremefail_error/eff_data_em_extremefail)**2 + (eff_MC_em_extremefail_error/eff_MC_em_extremefail)**2 )**0.5


        # w/ extreme fail
        eff_MC_em_LP   = 1.-rrv_eff_MC_em.getVal()   - eff_MC_em_extremefail
        eff_data_em_LP = 1.-rrv_eff_data_em.getVal() - eff_data_em_extremefail

        eff_MC_em_LP_err   = TMath.Sqrt( rrv_eff_MC_em.getError()  **2 + eff_MC_em_extremefail_error**2 )
        eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2 + eff_data_em_extremefail_error**2 )

        pure_wtagger_sf_em_LP     = eff_data_em_LP / eff_MC_em_LP
        pure_wtagger_sf_em_LP_err = pure_wtagger_sf_em_LP * ( (eff_data_em_LP_err/eff_data_em_LP)**2 + (eff_MC_em_LP_err/eff_MC_em_LP)**2 )**0.5

        # w/o extreme fail
        tmpq_eff_MC_em_LP   = 1. - rrv_eff_MC_em.getVal()
        tmpq_eff_data_em_LP = 1. - rrv_eff_data_em.getVal()

        tmpq_eff_MC_em_LP_err   = TMath.Sqrt( rrv_eff_MC_em.getError()  **2)
        tmpq_eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2)

        pureq_wtagger_sf_em_LP = tmpq_eff_data_em_LP / tmpq_eff_MC_em_LP
        pureq_wtagger_sf_em_LP_err = pureq_wtagger_sf_em_LP*TMath.Sqrt( (tmpq_eff_data_em_LP_err/tmpq_eff_data_em_LP)**2 + (tmpq_eff_MC_em_LP_err/tmpq_eff_MC_em_LP)**2 )


        print "-----------------------------------------------------------------------------------------------------------------------------"
        print "                                     HP                                    "
        print "-----------------------------------------------------------------------------------------------------------------------------"
        print "Pure W-tagging SF            : %0.3f +/- %0.3f" %(pure_wtagger_sf_em, pure_wtagger_sf_em_err)
        print "Pure W-tagging mean shift    : %0.3f +/- %0.3f" %(pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em)
        print "Pure W-tagging sigma enlarge : %0.3f +/- %0.3f" %(pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em)
        print ""
        print "Parameter                 Data                          Simulation                          Data/Simulation"
        print " < m >              %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_mean_data_em.getVal(),rrv_mean_data_em.getError(),rrv_mean_MC_em.getVal(),rrv_mean_MC_em.getError()    , rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal()  ,  mean_sf_error)
        print " #sigma             %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_sigma_data_em.getVal(),rrv_sigma_data_em.getError(),rrv_sigma_MC_em.getVal(),rrv_sigma_MC_em.getError(), rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal(),  sigma_sf_error)
        print ""
        print "HP W-tag eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_eff_data_em.getVal(),rrv_eff_data_em.getError(),rrv_eff_MC_em.getVal(),rrv_eff_MC_em.getError()        , rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()    ,  eff_sf_error)
        print ""
        print "                                EXTREME FAIL                               "
        print "-----------------------------------------------------------------------------------------------------------------------------"
        print ""
        print "Parameter                     Data                          Simulation                          Data/Simulation"
        print "Extreme fail eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_extremefail,eff_data_em_extremefail_error,eff_MC_em_extremefail,eff_MC_em_extremefail_error, eff_SF_extremefail, eff_SF_extremefail_error)
        print ""
        print "                                    LP                                     "
        print "-----------------------------------------------------------------------------------------------------------------------------"
        print "Parameter                                   Data                          Simulation                          Data/Simulation"
        print "LP W-tag eff+SF ( w/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_LP,eff_data_em_LP_err,eff_MC_em_LP,eff_MC_em_LP_err,pure_wtagger_sf_em_LP,pure_wtagger_sf_em_LP_err)
        print "LP W-tag eff+SF (wo/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(tmpq_eff_data_em_LP,tmpq_eff_data_em_LP_err,tmpq_eff_MC_em_LP,tmpq_eff_MC_em_LP_err,pureq_wtagger_sf_em_LP,pureq_wtagger_sf_em_LP_err)
        print ""
        print "-----------------------------------------------------------------------------------------------------------------------------"


        self.boostedW_fitter_em.file_out_ttbar_control.write("\n                                     HP                                    ")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nPure W-tagging SF of          : %0.3f +/- %0.3f" %(pure_wtagger_sf_em, pure_wtagger_sf_em_err))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nPure W-tagging mean shift     : %0.3f +/- %0.3f" %(pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nPure W-tagging sigma enlarge  : %0.3f +/- %0.3f" %(pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nParameter                 Data                          Simulation                          Data/Simulation")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n < m >              %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_mean_data_em.getVal(),rrv_mean_data_em.getError(),rrv_mean_MC_em.getVal(),rrv_mean_MC_em.getError()    , rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal()  ,  mean_sf_error))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n #sigma             %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_sigma_data_em.getVal(),rrv_sigma_data_em.getError(),rrv_sigma_MC_em.getVal(),rrv_sigma_MC_em.getError(), rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal(),  sigma_sf_error))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nHP W-tag eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_eff_data_em.getVal(),rrv_eff_data_em.getError(),rrv_eff_MC_em.getVal(),rrv_eff_MC_em.getError(), rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal(),  eff_sf_error))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n                                EXTREME FAIL                               ")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nParameter                     Data                          Simulation                          Data/Simulation")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nExtreme fail eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_extremefail,eff_data_em_extremefail_error,eff_MC_em_extremefail,eff_MC_em_extremefail_error, eff_SF_extremefail, eff_SF_extremefail_error))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n                                    LP                                     ")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nParameter                                   Data                          Simulation                          Data/Simulation")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nLP W-tag eff+SF ( w/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_LP,eff_data_em_LP_err,eff_MC_em_LP,eff_MC_em_LP_err,pure_wtagger_sf_em_LP,pure_wtagger_sf_em_LP_err))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\nLP W-tag eff+SF (wo/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(tmpq_eff_data_em_LP,tmpq_eff_data_em_LP_err,tmpq_eff_MC_em_LP,tmpq_eff_MC_em_LP_err,pureq_wtagger_sf_em_LP,pureq_wtagger_sf_em_LP_err))
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n")
        self.boostedW_fitter_em.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")

class doFit_wj_and_wlvj:

    ## COnstructor: Input is channel (mu,ele,em), range in mj and a workspace
    def __init__(self, in_channel, in_sample, in_mj_min=40, in_mj_max=130, input_workspace=None):
      
      RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
      RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)
      
      ### set the channel type
      self.channel = in_channel
      
      print "CHANNEL = %s" %in_channel

      ### shapes to be used in mj                                                                                                                                          
      self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string)()
      self.mj_shape["TTbar"]            = "2Gaus_ttbar"
      self.mj_shape["TTbar_fail"]       = "ErfExp_ttbar_failtau2tau1cut"
      
      self.mj_shape["TTbar_realW"]      = "2Gaus_ttbar"
      self.mj_shape["TTbar_realW_fail"] = "2Gaus_ttbar"
      self.mj_shape["TTbar_fakeW"]      = "ErfExp_ttbar"
      self.mj_shape["TTbar_fakeW_fail"] = "ErfExp_ttbar_failtau2tau1cut"
      
      if (options.tau2tau1cutHP==0.60):
        print "Using Tau21 leptonic HP cut of" ,options.tau2tau1cutHP 
        self.mj_shape["STop"]                 = "ErfExpGaus_sp"
        self.mj_shape["STop_fail"]            = "ExpGaus"
        self.mj_shape["STop_extremefail"]     = "Exp"
        self.mj_shape["VV"]                   = "ErfExpGaus_sp"
        self.mj_shape["VV_fail"]              = "ExpGaus"
        self.mj_shape["VV_extremefail"]       = "Exp"
        self.mj_shape["WJets0"]               = "ErfExp"
        self.mj_shape["WJets0_fail"]          = "ErfExp"
        self.mj_shape["WJets0_extremefail"]   = "Exp"
        
        self.mj_shape["bkg_mc_fail"]          = "ErfExp_ttbar_failtau2tau1cut"
        self.mj_shape["bkg_data_fail"]        = "ErfExp_ttbar_failtau2tau1cut"    
        
        self.mj_shape["signal_mc_fail"]       = "GausChebychev_ttbar_failtau2tau1cut"        
        self.mj_shape["signal_data_fail"]     = "GausExp_ttbar_failtau2tau1cut"
        
      elif (options.tau2tau1cutHP==0.45 or options.tau2tau1cutHP==0.55):
        self.mj_shape["STop"]             = "ErfExpGaus_sp"       
        self.mj_shape["STop_fail"]          = "ExpGaus"    
        self.mj_shape["STop_extremefail"]   = "Exp"
        self.mj_shape["VV"]                 = "ExpGaus"
        self.mj_shape["VV_fail"]            = "ExpGaus"
        self.mj_shape["VV_extremefail"]     = "Exp"
        self.mj_shape["WJets0"]             = "ErfExp"
        # self.mj_shape["WJets0"]             = "ErfExp"  #test
        # self.mj_shape["WJets0_fail"]        = "Exp"
        self.mj_shape["WJets0_fail"]        = "ErfExp" #test
        self.mj_shape["WJets0_extremefail"] = "Exp"
        
        self.mj_shape["bkg_mc_fail"]          = "ErfExp_ttbar_failtau2tau1cut"
        self.mj_shape["bkg_data_fail"]        = "ErfExp_ttbar_failtau2tau1cut"    
        
        # self.mj_shape["signal_mc_fail"]       = "GausChebychev_ttbar_failtau2tau1cut"
        # self.mj_shape["signal_data_fail"]     = "GausChebychev_ttbar_failtau2tau1cut"
        
        self.mj_shape["signal_mc_fail"]       = "2Gaus_ttbar"
        self.mj_shape["signal_data_fail"]     = "2Gaus_ttbar"
        # if options.usePuppiSD:
 #          self.mj_shape["signal_mc_fail"]       = "2Gaus_ttbar"
 #          self.mj_shape["signal_data_fail"]     = "2Gaus_ttbar"
        # if options.tau2tau1cutHP==0.55:
        #   self.mj_shape["VV_fail"]              = "ErfExp"

        
      else:
        print "NO CHANNEL IS DEFINED!!! ABORT!!"
        sys.exit(0)
         
      self.mj_shape["bkg_data"]             = "ErfExp_ttbar" 
      self.mj_shape["bkg_mc"]               = "ErfExp_ttbar" 
      
      self.mj_shape["signal_data"]          = "2Gaus_ttbar"
      self.mj_shape["signal_mc"]            = "2Gaus_ttbar"
      
      self.mj_shape["data_extremefail"]     = "Exp_ttbar_extremefailtau2tau1cut"
      self.mj_shape["mc_extremefail"]       = "Exp_ttbar_extremefailtau2tau1cut"
      
      self.mj_shape["data_bkg_extremefail"] = "Exp_bkg_extremefailtau2tau1cut"
      self.mj_shape["mc_bkg_extremefail"]   = "Exp_bkg_extremefailtau2tau1cut"
        

      self.Lumi=2300
      self.BinWidth_mj = 5.
      self.narrow_factor = 1.

      self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor
      nbins_mj         = int( (in_mj_max - in_mj_min) / self.BinWidth_mj )
      in_mj_max        = in_mj_min+nbins_mj*self.BinWidth_mj
      
      jetMass = "pruned jet mass"
      if options.usePuppiSD: jetMass = "PUPPI softdrop jet mass"

      rrv_mass_j = RooRealVar("rrv_mass_j", jetMass ,(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV")
      rrv_mass_j.setBins(nbins_mj)
 
      # Create workspace and import variable
      if input_workspace is None:
          self.workspace4fit_ = RooWorkspace("workspace4fit_","Workspace4fit_")
      else:
          self.workspace4fit_ = input_workspace
      getattr(self.workspace4fit_,"import")(rrv_mass_j)

      # Signal region between 65 and 105 GeV
      self.mj_sideband_lo_min = in_mj_min
      self.mj_sideband_lo_max = 65
      self.mj_signal_min      = 65
      self.mj_signal_max      = 105
      self.mj_sideband_hi_min = 105
      self.mj_sideband_hi_max = in_mj_max
 
      # Setting ranges...
      rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max) # 30-65 GeV
      rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max)   # 65-105 GeV
      rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max) # 105-135 GeV
      rrv_mass_j.setRange("controlsample_fitting_range",40,130) # ---> what is this????
      
      
      postfix = ""
      if options.use76X: 
        postfix ="_76X"
      # Tree directory and input files
      self.file_Directory         = "/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/Wtag/WWTree_%s/"%(self.channel)
      if options.doLowMass: 
        self.file_Directory       = "/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/Wtag/lowMass/WWTree_%s/"%(self.channel)
      self.file_data              = ("ExoDiBosonAnalysis.WWTree_data%s.root") %postfix
      self.file_pseudodata        = ("ExoDiBosonAnalysis.WWTree_pseudodata%s%s.root")%(in_sample,postfix)       
      self.file_WJets0_mc         = ("ExoDiBosonAnalysis.WWTree_WJets%s.root") %postfix
      self.file_VV_mc             = ("ExoDiBosonAnalysis.WWTree_VV%s.root") %postfix          
      self.file_TTbar_mc          = ("ExoDiBosonAnalysis.WWTree_TTbar%s%s.root")%(in_sample,postfix)      
      if options.use76X: 
        self.file_TTbar_mc        = ("ExoDiBosonAnalysis.WWTree_TTbar_powheg_76X.root")
        if in_sample.find("herwig")!=-1: self.file_TTbar_mc        = ("ExoDiBosonAnalysis.WWTree_TTbar_herwig_76X.root")
      self.file_STop_mc           = ("ExoDiBosonAnalysis.WWTree_STop%s.root") %postfix          
      
      
      if options.usePuppiSD:
        self.file_Directory         = "/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/Wtag/WWTree_%s/"%(self.channel)
        self.file_data              = ("ExoDiBosonAnalysis.WWTree_data_76X_PUPPISD.root")
        self.file_pseudodata        = ("ExoDiBosonAnalysis.WWTree_pseudodata_76X_PUPPISD.root")     
        self.file_WJets0_mc         = ("ExoDiBosonAnalysis.WWTree_WJets_76X_PUPPISD.root")
        self.file_VV_mc             = ("ExoDiBosonAnalysis.WWTree_VV_76X_PUPPISD.root")        
        self.file_TTbar_mc          = ("ExoDiBosonAnalysis.WWTree_TTbar_powheg_76X_PUPPISD.root")
        self.file_STop_mc           = ("ExoDiBosonAnalysis.WWTree_STop_76X_PUPPISD.root")
        
        
        
        
        
      # Define Tau21 WP
      self.wtagger_label = options.category;

      if self.wtagger_label == "HP" :
        self.wtagger_cut = options.tau2tau1cutHP
        self.wtagger_cut_min = 0.

      if self.wtagger_label == "LP":
          self.wtagger_cut = options.tau2tau1cutLP 
          self.wtagger_cut_min = options.tau2tau1cutHP

      if self.wtagger_label == "nocut":
          self.wtagger_cut = 10000
      
      if options.usePuppiSD: 
        postfix = postfix + "_PuppiSD"
      if options.useDDT: 
        postfix = postfix + "_DDT"
        print postfix
        print postfix
        print postfix
        print postfix
      if (options.tau2tau1cutHP==0.60):
        self.wtagger_label = self.wtagger_label + "0v60%s%s"%(in_sample,postfix) 
      elif (options.tau2tau1cutHP==0.45):  
        self.wtagger_label = self.wtagger_label + "0v45%s%s"%(in_sample,postfix) 
      elif (options.tau2tau1cutHP==0.55):  
        self.wtagger_label = self.wtagger_label + "0v55%s%s"%(in_sample,postfix)  
       
      self.color_palet = ROOT.std.map(ROOT.std.string, int) ()
      self.color_palet["data"]              = 1
      self.color_palet["WJets"]             = 2
      self.color_palet["VV"]                = 4
      self.color_palet["WW_EWK"]            = 6
      self.color_palet["STop"]              = 7
      self.color_palet["TTbar"]             = 210
      self.color_palet["ggH"]               = 1
      self.color_palet["vbfH"]              = 12
      self.color_palet["Signal"]            = 1
      self.color_palet["Uncertainty"]       = 1
      self.color_palet["Other_Backgrounds"] = 1    
      
      # Cuts (dont really need this cuts because thay are already implemented in ROOT tree)
      self.vpt_cut      = 200   # hadronic and leptonic W cut
      self.mass_lvj_max = 5000. # invariant mass of 3 body max
      self.mass_lvj_min = 0.    # invariant mass of 3 body min
      self.pfMET_cut    = 40.    # missing transverse energy
      self.lpt_cut      = 53.    # lepton pT
      self.AK8_pt_min   = 200
      self.AK8_pt_max   = 5000  
      if self.channel  == "el":
        self.pfMET_cut = 80
        self.lpt_cut = 120
        
      
      # Out .txt file
      self.file_ttbar_control_txt = "WtaggingSF_%s%.2f%s.txt"%(self.channel,options.tau2tau1cutHP,options.sample)
      self.file_out_ttbar_control = open(self.file_ttbar_control_txt,"w")
                                                                                                                                                             
      setTDRStyle()

    def get_datasets_fit_minor_bkg(self):

      rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
      
      print "#########################################"
      print "################ TTbar %s################"%options.sample
      print "#########################################"
      print ""

      
      if options.fitTT:
        self.get_mj_dataset(self.file_TTbar_mc,"_TTbar_realW")
        self.get_mj_dataset(self.file_TTbar_mc,"_TTbar_fakeW")

        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar_realW",self.mj_shape["TTbar_realW"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar_realW_failtau2tau1cut",self.mj_shape["TTbar_realW_fail"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar_fakeW",self.mj_shape["TTbar_fakeW"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_TTbar_mc,"_TTbar_fakeW_failtau2tau1cut",self.mj_shape["TTbar_fakeW_fail"],self.channel,self.wtagger_label)
      else:
        self.get_mj_dataset(self.file_TTbar_mc,"_TTbar")


        # Build single-t fit pass and fail distributions
        print ""
        print ""
        print "##################################################"
        print "############### Single Top DataSet ###############"
        print "##################################################"
        print ""
        print ""

        self.get_mj_dataset(self.file_STop_mc,"_STop")

        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop"                        ,self.mj_shape["STop"],self.channel,self.wtagger_label) #Start value and range of parameters defined in PDFs/MakePDF.cxx
        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_failtau2tau1cut"        ,self.mj_shape["STop_fail"],self.channel,self.wtagger_label)

        # fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_extremefailtau2tau1cut" ,self.mj_shape["STop_extremefail"],self.channel,self.wtagger_label)


        ### Build WJet fit pass and fail distributions
        print "###########################################"
        print "############### WJets Pythia ##############"
        print "###########################################"
        print ""
        print ""

        self.get_mj_dataset(self.file_WJets0_mc,"_WJets0")

        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0",self.mj_shape["WJets0"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_failtau2tau1cut",self.mj_shape["WJets0_fail"],self.channel,self.wtagger_label)
        # fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_extremefailtau2tau1cut",self.mj_shape["WJets0_extremefail"],self.channel,self.wtagger_label)



        # Build VV fit pass and fail distributions
        print "#########################################"
        print "############### VV Pythia ###############"
        print "#########################################"
        print ""
        print ""

        self.get_mj_dataset(self.file_VV_mc,"_VV")

        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_failtau2tau1cut",self.mj_shape["VV_fail"],self.channel,self.wtagger_label)
        # fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_extremefailtau2tau1cut",self.mj_shape["VV_extremefail"],self.channel,self.wtagger_label)

        print "################################################"
        print "############## Pseudo Data Powheg ##############"
        print "################################################"
        print ""
        print ""
        self.get_mj_dataset(self.file_pseudodata,"_TotalMC")

        print "#################################"
        print "############# Data ##############"
        print "#################################"
        print ""
        print ""
        self.get_mj_dataset(self.file_data,"_data")


        # self.print_yields()

        self.constrainslist_data = ROOT.std.vector(ROOT.std.string)()
        self.constrainslist_mc   = ROOT.std.vector(ROOT.std.string)()

        ScaleFactorTTbarControlSampleFit(self.workspace4fit_,self.mj_shape,self.color_palet,self.constrainslist_data,self.constrainslist_mc,"",self.channel,self.wtagger_label,self.AK8_pt_min,self.AK8_pt_max)

        rrv_scale_number                      = self.workspace4fit_.var("rrv_scale_number_TTbar_STop_VV_WJets").getVal()
        rrv_scale_number_fail                 = self.workspace4fit_.var("rrv_scale_number_TTbar_STop_VV_WJets_fail").getVal()

        print " Pass MC / all data = %.3f" %(rrv_scale_number)
        print " Fail MC / all data = %.3f" %(rrv_scale_number_fail)




    def print_yields(self):

        # Print dataset yields in the signal region
        print ""
        print ""
        print ""
        self.workspace4fit_.var("rrv_number_dataset_signal_region_data_"    +self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_VV_"      +self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_WJets0_"  +self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_STop_"    +self.channel+"_mj").Print()
        self.workspace4fit_.var("rrv_number_dataset_signal_region_TTbar_"   +self.channel+"_mj").Print()
        print ""
        print ""
        print ""

        number_dataset_signal_region_data_mj                      = self.workspace4fit_.var("rrv_number_dataset_signal_region_data_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_error2_data_mj               = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data_"+self.channel+"_mj").getVal()
        
        number_dataset_signal_region_TotalMC_mj                   = self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_error2_TotalMC_mj            = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC_"+self.channel+"_mj").getVal()
        
        number_dataset_signal_region_before_cut_data_mj           = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_data_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_before_cut_error2_data_mj    = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_data_"+self.channel+"_mj").getVal()
        
        number_dataset_signal_region_before_cut_TotalMC_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_TotalMC_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_before_cut_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_TotalMC_"+self.channel+"_mj").getVal()
        
        wtagger_eff_MC                                            = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj
        wtagger_eff_data                                          = number_dataset_signal_region_data_mj/number_dataset_signal_region_before_cut_data_mj

        wtagger_eff_reweight                                      = wtagger_eff_data/wtagger_eff_MC
        wtagger_eff_reweight_err                                  = wtagger_eff_reweight*TMath.Sqrt(number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +number_dataset_signal_region_before_cut_error2_data_mj/number_dataset_signal_region_before_cut_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_before_cut_error2_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj)
        
        print ""
        print "Nr. data events in signal_region                  : %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj**.5)
        print ""
        print "Nr. MC events in signal_region                    : %s +/- sqrt(%s)"%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj)
        print ""
        print "Nr. dataevents in signalregion before cut on tau21: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj)
        print ""
        print "Nr. MC events in signalregion before cut on tau21 : %s +/- sqrt(%s) "%(number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj)
        print ""
        print ""
        print ""                                                     
        print "W-tagging efficiency (pre-fit):"
        print "W-tagging eff. MC       = %.3f "%(wtagger_eff_MC)
        print "W-tagging eff. data     = %.3f "%(wtagger_eff_data)
        print "W-tagging SF            = %.3f +/- %.3f"%(wtagger_eff_reweight, wtagger_eff_reweight_err)
        print ""
        print ""
        
        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel))
        self.file_out_ttbar_control.write("Nr. events in signal_region                        : %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj))
        self.file_out_ttbar_control.write("Nr. TotalMC in signal_region                       : %s +/- sqrt(%s) \n"%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj))
        self.file_out_ttbar_control.write("event number of data in signalregion before_cut    : %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj))
        self.file_out_ttbar_control.write("event number of TotalMC in signal_region before_cut: %s +/- sqrt(%s) \n"%(number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj))
        self.file_out_ttbar_control.write("wtagger_eff_MC         = %s       \n"%(wtagger_eff_MC ))
        self.file_out_ttbar_control.write("wtagger_eff_data       = %s       \n"%(wtagger_eff_data ))
        self.file_out_ttbar_control.write("wtagger_eff_reweight   = %s +/- %s\n"%(wtagger_eff_reweight, wtagger_eff_reweight_err))
            
    # Loop over trees
    def get_mj_dataset(self,in_file_name, label, jet_mass="Whadr_pruned"): 
      
      if options.usePuppiSD: jet_mass="Whadr_puppi_softdrop"
      
      print "Using mass variable " ,jet_mass
    
      fileIn_name = TString(self.file_Directory+in_file_name)
      
      print "Using file " ,fileIn_name
      
      fileIn      = TFile(fileIn_name.Data())
      treeIn      = fileIn.Get("tree")
      
      rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
      rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)

      # Mj dataset before tau2tau1 cut : Passed
      rdataset_mj     = RooDataSet("rdataset"     +label+"_"+self.channel+"_mj","rdataset"    +label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rdataset4fit_mj = RooDataSet("rdataset4fit" +label+"_"+self.channel+"_mj","rdataset4fit"+label+"_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rrv_number_pass = RooRealVar("rrv_number_ttbar"+label+"_passtau2tau1cut_em_mj","rrv_number_ttbar"+label+"_passtau2tau1cut_em_mj",0.,10000000.) #LUCA
  
      # Mj dataset before tau2tau1 cut : Total
      rdataset_beforetau2tau1cut_mj     = RooDataSet("rdataset"     +label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset"    +label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rdataset4fit_beforetau2tau1cut_mj = RooDataSet("rdataset4fit" +label+"_beforetau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_beforetau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rrv_number_before = RooRealVar("rrv_number_ttbar"+label+"_beforetau2tau1cut_em_mj","rrv_number_ttbar"+label+"_beforetau2tau1cut_em_mj",0.,10000000.) #LUCA
 
      ### Mj dataset failed tau2tau1 cut :
      rdataset_failtau2tau1cut_mj     = RooDataSet("rdataset"     +label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset"    +label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rdataset4fit_failtau2tau1cut_mj = RooDataSet("rdataset4fit" +label+"_failtau2tau1cut_"+self.channel+"_mj","rdataset4fit"+label+"_failtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rrv_number_fail = RooRealVar("rrv_number_ttbar"+label+"_failtau2tau1cut_em_mj","rrv_number_ttbar"+label+"_failtau2tau1cut_em_mj",0.,10000000.) #LUCA

      ### Mj dataset extreme failed tau2tau1 cut: > 0.75
      rdataset_extremefailtau2tau1cut_mj     = RooDataSet("rdataset"    +label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset"     +label+"_extremefailtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rdataset4fit_extremefailtau2tau1cut_mj = RooDataSet("rdataset4fit"+label+"_extremefailtau2tau1cut_"+self.channel+"_mj","rdataset4fit" +label+"_extremefailtau2tau1cut_"+self.channel+"_mj",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight) )
      rrv_number_extremefail = RooRealVar("rrv_number_ttbar"+label+"_extremefailtau2tau1cut_em_mj","rrv_number_ttbar"+label+"_extremefailtau2tau1cut_em_mj",0.,10000000.) #LUCA
      
      # category_cut = RooCategory("category_cut"+"_"+self.channel,"category_cut"+"_"+self.channel) #---->Think this can be removed!!!!
 #      category_cut.defineType("cut",1)
 #      category_cut.defineType("beforecut",2)
 #      combData4cut = RooDataSet("combData4cut"+"_"+self.channel,"combData4cut"+"_"+self.channel,RooArgSet(rrv_mass_j, category_cut, rrv_weight),RooFit.WeightVar(rrv_weight) )
      
      
      # Define categories
      if self.workspace4fit_.cat("category_p_f"+"_"+self.channel):
        category_p_f = self.workspace4fit_.cat("category_p_f"+"_"+self.channel)
      else:
        category_p_f = RooCategory("category_p_f"+"_"+self.channel,"category_p_f"+"_"+self.channel)
        category_p_f.defineType("pass")
        category_p_f.defineType("fail")
        getattr(self.workspace4fit_,"import")(category_p_f)
      
      combData_p_f = RooDataSet("combData_p_f"+label+"_"+self.channel,"combData_p_f"+label+"_"+self.channel,RooArgSet(rrv_mass_j, category_p_f, rrv_weight),RooFit.WeightVar(rrv_weight))
      
      print "N entries: ", treeIn.GetEntries()
      
      hnum_4region                    = TH1D("hnum_4region"       +label+"_"+self.channel,"hnum_4region"        +label+"_"+self.channel,4, -1.5, 2.5) # m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
      hnum_4region_error2             = TH1D("hnum_4region_error2"+label+"_"+self.channel,"hnum_4region_error2" +label+"_"+self.channel,4, -1.5, 2.5) # m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total

      hnum_4region_before_cut         = TH1D("hnum_4region_before_cut"        +label+"_"+self.channel,"hnum_4region_before_cut"       +label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total
      hnum_4region_before_cut_error2  = TH1D("hnum_4region_before_cut_error2" +label+"_"+self.channel,"hnum_4region_before_cut_error2"+label+"_"+self.channel,4,-1.5,2.5);# m_j -1: sb_lo; 0:signal_region; 1: sb_hi; 2:total

      hnum_2region                    = TH1D("hnum_2region"       +label+"_"+self.channel,"hnum_2region"        +label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total --> There is only 1 and that is SIGNAL REGION?!
      hnum_2region_error2             = TH1D("hnum_2region_error2"+label+"_"+self.channel,"hnum_2region_error2" +label+"_"+self.channel,2,-0.5,1.5);# m_lvj 0: signal_region; 1: total
      
  
      #-------------------------------------------------------------------------------------------
      # Loop over tree entries
      tmp_scale_to_lumi = 1
      for i in range(treeIn.GetEntries()):
          if i % 5000 == 0: print "iEntry: ",i
          treeIn.GetEntry(i)
          
          if TString(label).Contains("realW") and not getattr(treeIn,"Whadr_isW"): 
            continue
          if TString(label).Contains("fakeW") and     getattr(treeIn,"Whadr_isW"): 
            continue
            
          #b-tag veto!
          if getattr(treeIn,"Whadr_csv") > options.csvMax:
            continue
                
          discriminantCut = 0
          wtagger = getattr(treeIn,"Whadr_tau21")
          if options.useDDT:
            wtagger = getattr(treeIn,"Whadr_puppi_tau2")/getattr(treeIn,"Whadr_puppi_tau1")+ (0.063 * math.log( (pow( getattr(treeIn,"Whadr_puppi_softdrop"),2))/getattr(treeIn,"Whadr_puppi_pt") ))
          # print "Puppi tau21 " ,getattr(treeIn,"Whadr_puppi_tau2")/getattr(treeIn,"Whadr_puppi_tau1")
          # print "Puppi softdrop " ,getattr(treeIn,"Whadr_puppi_softdrop")
          # print "Puppi pt " ,getattr(treeIn,"Whadr_puppi_pt")
          # print "Puppi ddt " ,wtagger

          if wtagger <= options.tau2tau1cutHP: # HP
              discriminantCut = 2
          elif wtagger > options.tau2tau1cutHP and wtagger <= options.tau2tau1cutLP: #LP
              discriminantCut = 1
          elif wtagger > options.tau2tau1cutLP: # Extreme fail
              discriminantCut = 0

          tmp_jet_mass = getattr(treeIn, jet_mass);
          
          if i==0: 
            tmp_scale_to_lumi = treeIn.lumiweight*self.Lumi ## weigth for xs and lumi

          if not TString(label).Contains("data"):
            tmp_event_weight     = getattr(treeIn,"weight")*self.Lumi 
            tmp_event_weight4fit = getattr(treeIn,"genweight")*getattr(treeIn,"puweight")*getattr(treeIn,"hltweight") #Missing b-tag weight and Vtg weight!
            tmp_event_weight4fit = tmp_event_weight4fit*treeIn.lumiweight*self.Lumi /tmp_scale_to_lumi
            # tmp_event_weight4fit = tmp_event_weight/tmp_scale_to_lumi
          else:
            tmp_event_weight = 1.
            tmp_event_weight4fit = 1.

          #  HP category
          if discriminantCut == 2  and tmp_jet_mass > rrv_mass_j.getMin() and tmp_jet_mass < rrv_mass_j.getMax():   
             rrv_mass_j.setVal(tmp_jet_mass)
             
             rdataset_mj    .add(RooArgSet(rrv_mass_j), tmp_event_weight)
             rdataset4fit_mj.add(RooArgSet(rrv_mass_j), tmp_event_weight4fit)

             if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                 hnum_4region.Fill(-1,tmp_event_weight )
                 
             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                 # hnum_2region.Fill(1,tmp_event_weight)
                 hnum_4region.Fill(0,tmp_event_weight)
                 hnum_4region_error2.Fill(0,tmp_event_weight*tmp_event_weight)
             if tmp_jet_mass >= self.mj_sideband_hi_min and tmp_jet_mass < self.mj_sideband_hi_max:
                 hnum_4region.Fill(1,tmp_event_weight)

             hnum_4region.Fill(2,tmp_event_weight) 
             
             # category_cut.setLabel("cut");
 #             combData4cut.add(RooArgSet(rrv_mass_j,category_cut),tmp_event_weight4fit)
             category_p_f.setLabel("pass")
             combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight)
          
          # TOTAL category (no Tau21 )
          if discriminantCut == 2 or discriminantCut == 1 or discriminantCut == 0 and (rrv_mass_j.getMin() < tmp_jet_mass < rrv_mass_j.getMax()):   

              rrv_mass_j.setVal(tmp_jet_mass)

              if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass <self.mj_signal_max :
                 hnum_4region_before_cut.Fill(0,tmp_event_weight)
                 hnum_4region_before_cut_error2.Fill(0,tmp_event_weight*tmp_event_weight)

              rdataset_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight)
              rdataset4fit_beforetau2tau1cut_mj.add(RooArgSet(rrv_mass_j),tmp_event_weight4fit)
          
          # 1 minus HP category (LP+extreme fail)   
          if (discriminantCut==1 or discriminantCut==0) and (rrv_mass_j.getMin() < tmp_jet_mass < rrv_mass_j.getMax()):
  
              rrv_mass_j.setVal(tmp_jet_mass)

              rdataset_failtau2tau1cut_mj     .add(RooArgSet(rrv_mass_j), tmp_event_weight)
              rdataset4fit_failtau2tau1cut_mj .add(RooArgSet(rrv_mass_j), tmp_event_weight4fit )
    
              category_p_f.setLabel("fail");
              combData_p_f.add(RooArgSet(rrv_mass_j,category_p_f),tmp_event_weight)
           #-------------------------------------------------------------------------------------------
           # Extreme fail category (Tau21 > LP_max)   
          if discriminantCut==0 and (rrv_mass_j.getMin() < tmp_jet_mass < rrv_mass_j.getMax()):

            rdataset_extremefailtau2tau1cut_mj     .add(RooArgSet(rrv_mass_j), tmp_event_weight)
            rdataset4fit_extremefailtau2tau1cut_mj .add(RooArgSet(rrv_mass_j), tmp_event_weight4fit )

      print label
      print "rrv_scale_to_lumi"+label+"_"                       +self.channel
      rrv_scale_to_lumi                        = RooRealVar("rrv_scale_to_lumi"+label+"_"                       +self.channel,"rrv_scale_to_lumi"+label+"_"                       +self.channel,tmp_scale_to_lumi)
      rrv_scale_to_lumi_failtau2tau1cut        = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"       +self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"       +self.channel,tmp_scale_to_lumi)
      rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,tmp_scale_to_lumi)
      
      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)
      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut)
      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut)
        
      rrv_number_pass.setVal(rdataset_mj.sumEntries())
      rrv_number_pass.setError(TMath.Sqrt(rdataset_mj.sumEntries()))
      # rrv_number_pass.Print()
      rrv_number_before.setVal(rdataset_beforetau2tau1cut_mj.sumEntries()) 
      rrv_number_before.setError(TMath.Sqrt(rdataset_beforetau2tau1cut_mj.sumEntries())) 
      # rrv_number_before.Print()
      rrv_number_fail.setVal(rdataset_failtau2tau1cut_mj.sumEntries())
      rrv_number_fail.setError(TMath.Sqrt(rdataset_failtau2tau1cut_mj.sumEntries()))
      # rrv_number_fail.Print()
      rrv_number_extremefail.setVal(rdataset_extremefailtau2tau1cut_mj.sumEntries())
      rrv_number_extremefail.setError(TMath.Sqrt(rdataset_extremefailtau2tau1cut_mj.sumEntries()))
      # rrv_number_extremefail.Print()
 
      getattr(self.workspace4fit_,"import")(rrv_number_pass)
      getattr(self.workspace4fit_,"import")(rrv_number_before)
      getattr(self.workspace4fit_,"import")(rrv_number_fail)
      getattr(self.workspace4fit_,"import")(rrv_number_extremefail)

      #prepare m_j dataset
      rrv_number_dataset_sb_lo_mj                 = RooRealVar("rrv_number_dataset_sb_lo"               +label+"_"+self.channel+"_mj","rrv_number_dataset_sb_lo"                +label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(1))
      rrv_number_dataset_signal_region_mj         = RooRealVar("rrv_number_dataset_signal_region"       +label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region"        +label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(2))
      rrv_number_dataset_signal_region_error2_mj  = RooRealVar("rrv_number_dataset_signal_region_error2"+label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_error2" +label+"_"+self.channel+"_mj",hnum_4region_error2.GetBinContent(2))

      rrv_number_dataset_signal_region_before_cut_mj        = RooRealVar("rrv_number_dataset_signal_region_before_cut"        +label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut"       +label+"_"+self.channel+"_mj",hnum_4region_before_cut.GetBinContent(2))
      rrv_number_dataset_signal_region_before_cut_error2_mj = RooRealVar("rrv_number_dataset_signal_region_before_cut_error2" +label+"_"+self.channel+"_mj","rrv_number_dataset_signal_region_before_cut_error2"+label+"_"+self.channel+"_mj",hnum_4region_before_cut_error2.GetBinContent(2))
      rrv_number_dataset_sb_hi_mj                           = RooRealVar("rrv_number_dataset_sb_hi"                           +label+"_"+self.channel+"_mj","rrv_number_dataset_sb_hi"                          +label+"_"+self.channel+"_mj",hnum_4region.GetBinContent(3))

      getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_lo_mj)
      getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_mj)
      getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_error2_mj)
      getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_mj)
      getattr(self.workspace4fit_,"import")(rrv_number_dataset_signal_region_before_cut_error2_mj)
      getattr(self.workspace4fit_,"import")(rrv_number_dataset_sb_hi_mj)
      getattr(self.workspace4fit_,"import")(combData_p_f)

      print "N_rdataset_mj: "
      getattr(self.workspace4fit_,"import")(rdataset_mj)
      getattr(self.workspace4fit_,"import")(rdataset4fit_mj)
      getattr(self.workspace4fit_,"import")(rdataset_beforetau2tau1cut_mj)
      getattr(self.workspace4fit_,"import")(rdataset4fit_beforetau2tau1cut_mj)
      getattr(self.workspace4fit_,"import")(rdataset_failtau2tau1cut_mj)
      getattr(self.workspace4fit_,"import")(rdataset4fit_failtau2tau1cut_mj)
      getattr(self.workspace4fit_,"import")(rdataset_extremefailtau2tau1cut_mj)
      getattr(self.workspace4fit_,"import")(rdataset4fit_extremefailtau2tau1cut_mj)

      rdataset_mj.Print()
      rdataset4fit_mj.Print()
      rdataset_failtau2tau1cut_mj.Print()
      rdataset4fit_failtau2tau1cut_mj.Print()
      rdataset_extremefailtau2tau1cut_mj.Print()
      rdataset4fit_extremefailtau2tau1cut_mj.Print()
      rrv_number_dataset_sb_lo_mj.Print()
      rrv_number_dataset_signal_region_mj.Print()
      rrv_number_dataset_signal_region_error2_mj.Print()
      rrv_number_dataset_signal_region_before_cut_mj.Print()
      rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
      rrv_number_dataset_sb_hi_mj.Print()

      rdataset_mj.Print()
      rdataset_beforetau2tau1cut_mj.Print()
      rdataset_failtau2tau1cut_mj.Print()
      rdataset_extremefailtau2tau1cut_mj.Print()
      rrv_number_dataset_signal_region_mj.Print()
      rrv_number_dataset_signal_region_error2_mj.Print()
      rrv_number_dataset_signal_region_before_cut_mj.Print()
      rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
      combData_p_f.Print("v")
      
      
      ### Start  main
if __name__ == '__main__':
  
  channel = options.channel
  
  print 'fitwtagger adding el+mu sample' #I am actually not doing a simoultaneous fit. So..... change this
  getSF()

