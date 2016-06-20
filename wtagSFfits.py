from optparse import OptionParser
import ROOT
import sys
import time
import math
import CMS_lumi, tdrstyle
from ROOT import *

parser = OptionParser()
parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="em")
parser.add_option('--HP', action="store", type="float",dest="tau2tau1cutHP",default=0.45)
parser.add_option('--LP', action="store", type="float",dest="tau2tau1cutLP",default=0.75)
parser.add_option('--sample', action="store",type="string",dest="sample",default="powheg")
parser.add_option('--fitTT', action='store_true', dest='fitTT', default=False, help='Only do ttbar fits')
parser.add_option('--fitMC', action='store_true', dest='fitMC', default=False, help='Only do MC fits')
parser.add_option('--doUnbinned',dest="doUnbinnedFit", default=False, action="store_true", help="Do unbinned fit")
parser.add_option('--76X',dest="use76X", default=False, action="store_true", help="Use 76X samples")
parser.add_option('--useDDT',dest="useDDT", default=False, action="store_true", help="Use DDT tagger")
parser.add_option('--usePuppiSD',dest="usePuppiSD", default=False, action="store_true", help="Use PUPPI+softdrop")

(options, args) = parser.parse_args()

ROOT.gSystem.Load(".//PlotStyle/Util_cxx.so")
ROOT.gSystem.Load(".//PlotStyle/PlotUtils_cxx.so")
ROOT.gSystem.Load(".//PDFs/PdfDiagonalizer_cc.so")
ROOT.gSystem.Load(".//PDFs/HWWLVJRooPdfs_cxx.so")
ROOT.gSystem.Load(".//PDFs/MakePdf_cxx.so")
ROOT.gSystem.Load(".//BiasStudy/BiasUtils_cxx.so")
ROOT.gSystem.Load(".//FitUtils/FitUtils_cxx.so")

tdrstyle.setTDRStyle()
CMS_lumi.lumi_13TeV = "2.3 fb^{-1}"
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
iPos = 11
if( iPos==0 ): CMS_lumi.relPosX = 0.12
iPeriod = 4

gInterpreter.GenerateDictionary("std::map<std::string,std::string>", "map;string;string")
gStyle.SetOptTitle(0)
RooMsgService.instance().setGlobalKillBelow(RooFit.FATAL)

if options.noX: gROOT.SetBatch(True)

def getLegend():
  legend = TLegend(0.6010112,0.7183362,0.8202143,0.919833)
  legend.SetTextSize(0.032)
  legend.SetLineColor(0)
  legend.SetShadowColor(0)
  legend.SetLineStyle(1)
  legend.SetLineWidth(1)
  legend.SetFillColor(0)
  legend.SetFillStyle(0)
  legend.SetMargin(0.35)
  return legend

def getPavetext():
  addInfo = TPaveText(0.2510112,0.2066292,0.4202143,0.3523546,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  return addInfo
    
def drawFrameGetChi2(variable,fitResult,dataset,pdfModel):
    
    wpForPlotting ="%.2f"%options.tau2tau1cutHP
    wpForPlotting = wpForPlotting.replace(".","v")
    postfix=""
    if options.usePuppiSD: postfix = "_PuppiSD"
    if options.useDDT: postfix = "_PuppiSD_DDT"
    title = "Pruned jet mass (GeV)"
    if options.usePuppiSD or options.useDDT:  title = "PUPPI softdrop jet mass (GeV)"
    
    
    frame = variable.frame()
    dataset.plotOn(frame,RooFit.DataError(RooAbsData.Poisson),RooFit.Name(dataset.GetName ()))
    pdfModel.plotOn(frame,RooFit.VisualizeError(fitResult,1), RooFit.Name("Fit error"),RooFit.FillColor(kGray),RooFit.LineColor(kGray)) #,RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    pdfModel.plotOn(frame,RooFit.LineColor(kBlack),RooFit.Name(pdfModel.GetName())) #,RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    chi2 = frame.chiSquare(pdfModel.GetName(), dataset.GetName ())
    pdfModel.plotOn(frame,RooFit.Name( "Gaussian 2" ),RooFit.Components("gaus2*"),RooFit.LineStyle(kSolid),RooFit.LineColor(kBlue+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    pdfModel.plotOn(frame,RooFit.Name( "ErfExp comp." ),RooFit.Components("model_bkg*"),RooFit.LineStyle(9),RooFit.LineColor(kBlue+2),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    pdfModel.plotOn(frame,RooFit.Name( "Gaussian 1" ),RooFit.Components("gaus1*"),RooFit.LineStyle(kDashed),RooFit.LineColor(kRed+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    pdfModel.plotOn(frame,RooFit.Name( "Chebyshev comp." ),RooFit.Components("cheb*"),RooFit.LineStyle(kDashed),RooFit.LineColor(kBlue+3),RooFit.Normalization(1.0,RooAbsReal.RelativeExpected))
    dataset.plotOn(frame,RooFit.DataError(RooAbsData.Poisson),RooFit.Name(dataset.GetName ()))
    c1 =TCanvas("data_fail","",800,800)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetTitleOffset(0.90)
    frame.SetName("mjjFit")
    frame.GetYaxis().SetTitle("Events")
    frame.GetXaxis().SetTitle(title)
    frame.Draw()
    frame.SetMinimum(0.)
    legend = getLegend()
    legend.AddEntry(frame.findObject(dataset.GetName ()),"CMS data","lpe")
    legend.AddEntry(frame.findObject(pdfModel.GetName()),"Sim. fit","l")
    if frame.findObject("Gaussian 1"):
      legend.AddEntry(frame.findObject("Gaussian 1"),frame.findObject("Gaussian 1").GetName(),"l")
    if frame.findObject("Chebyshev comp."):
      legend.AddEntry(frame.findObject("Chebyshev comp."),frame.findObject("Chebyshev comp.").GetName(),"l")
    if frame.findObject("Gaussian 2"):
      legend.AddEntry(frame.findObject("Gaussian 2"),frame.findObject("Gaussian 2").GetName(),"l")
    if frame.findObject("ErfExp comp."):
      legend.AddEntry(frame.findObject("ErfExp comp."),frame.findObject("ErfExp comp.").GetName(),"l")
    legend.Draw("same")
    CMS_lumi.CMS_lumi(c1, iPeriod, iPos)
    addInfo = getPavetext()
    addInfo.AddText("#chi^{2}/nDOF = %.3f"%chi2)
    addInfo.Draw()
    c1.Update()
    cname = pdfModel.GetName()
    c1.SaveAs("plots/%s.pdf"%(cname))
    c1.SaveAs("plots/%s.root"%(cname))
    return chi2
    
def getSF():
    print "Getting W-tagging SF for cut " ,options.tau2tau1cutHP
    if options.useDDT: options.usePuppiSD = True
    boostedW_fitter_sim = doWtagFits()

def doFitsToMatchedTT():
    workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_")
    ttMC_fitter = initialiseFits("em", options.sample, 40, 130, workspace4fit_)

    ttMC_fitter.get_mj_dataset(ttMC_fitter.file_TTbar_mc,"_TTbar_realW")
    ttMC_fitter.get_mj_dataset(ttMC_fitter.file_TTbar_mc,"_TTbar_fakeW")

    fit_mj_single_MC(ttMC_fitter.workspace4fit_,ttMC_fitter.file_TTbar_mc,"_TTbar_realW",ttMC_fitter.mj_shape["TTbar_realW"],ttMC_fitter.channel,ttMC_fitter.wtagger_label)
    fit_mj_single_MC(ttMC_fitter.workspace4fit_,ttMC_fitter.file_TTbar_mc,"_TTbar_realW_failtau2tau1cut",ttMC_fitter.mj_shape["TTbar_realW_fail"],ttMC_fitter.channel,ttMC_fitter.wtagger_label)
    fit_mj_single_MC(ttMC_fitter.workspace4fit_,ttMC_fitter.file_TTbar_mc,"_TTbar_fakeW",ttMC_fitter.mj_shape["TTbar_fakeW"],ttMC_fitter.channel,ttMC_fitter.wtagger_label)
    fit_mj_single_MC(ttMC_fitter.workspace4fit_,ttMC_fitter.file_TTbar_mc,"_TTbar_fakeW_failtau2tau1cut",ttMC_fitter.mj_shape["TTbar_fakeW_fail"],ttMC_fitter.channel,ttMC_fitter.wtagger_label)
    
    print "Finished fitting matched tt MC! Plots can be found in plots_*_MCfits. Printing workspace:"
    workspace4fit_.Print()
        
def doFitsToMC():
    workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_")
    boostedW_fitter_em = initialiseFits("em", options.sample, 40, 130, workspace4fit_)
    boostedW_fitter_em.get_datasets_fit_minor_bkg()
    print "Finished fitting MC! Plots can be found in plots_*_MCfits. Printing workspace:"
    workspace4fit_.Print()
                    
def GetWtagScalefactors(workspace,fitter):

    #------------- Calculate scalefactors -------------
    ### Efficiency in data and MC
    rrv_eff_MC_em   = workspace.var("eff_ttbar_TotalMC_em_mj")
    rrv_mean_MC_em  = workspace.var("rrv_mean1_gaus_ttbar_TotalMC_em_mj")
    rrv_sigma_MC_em = workspace.var("rrv_sigma1_gaus_ttbar_TotalMC_em_mj")

    rrv_eff_data_em   = workspace.var("eff_ttbar_data_em_mj")
    rrv_mean_data_em  = workspace.var("rrv_mean1_gaus_ttbar_data_em_mj")
    rrv_sigma_data_em = workspace.var("rrv_sigma1_gaus_ttbar_data_em_mj")

    ## GET HP SCALEFACTOR AND UNCERTIANTIES
    wtagger_sf_em             = rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()
    wtagger_sf_em_err         = wtagger_sf_em * ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal() )**2 + (rrv_eff_MC_em.getError()/rrv_eff_MC_em.getVal())**2 )**0.5

    wtagger_mean_shift_em     = rrv_mean_data_em.getVal()-rrv_mean_MC_em.getVal()
    wtagger_mean_shift_err_em = (rrv_mean_data_em.getError()**2 + rrv_mean_MC_em.getError()**2)**0.5

    wtagger_sigma_enlarge_em     = rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()
    wtagger_sigma_enlarge_err_em = ((rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 + (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2 )**0.5* wtagger_sigma_enlarge_em

    mean_sf_error  = (rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal())   * ( (rrv_mean_data_em.getError()/rrv_mean_data_em.getVal())**2   +  (rrv_mean_MC_em.getError() /rrv_mean_MC_em.getVal())**2    )**0.5
    sigma_sf_error = (rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()) * ( (rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 +  (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2   )**0.5
    eff_sf_error   = (rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal())     * ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal())**2     +  (rrv_eff_MC_em.getError()  /rrv_eff_MC_em.getVal())**2     )**0.5

    ## GET EXTREME FAIL NUMBERS IN ORDER TO COMPUTE LP SF:
    rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj = workspace.var("rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj")
    rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj    = workspace.var("rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj")

    rrv_number_total_ttbar_TotalMC_em = workspace.var("rrv_number_ttbar_TotalMC_beforetau2tau1cut_em_mj")
    rrv_number_total_ttbar_data_em    = workspace.var("rrv_number_ttbar_data_beforetau2tau1cut_em_mj")

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

    wtagger_sf_em_LP     = eff_data_em_LP / eff_MC_em_LP
    wtagger_sf_em_LP_err = wtagger_sf_em_LP * ( (eff_data_em_LP_err/eff_data_em_LP)**2 + (eff_MC_em_LP_err/eff_MC_em_LP)**2 )**0.5

    # w/o extreme fail
    tmpq_eff_MC_em_LP   = 1. - rrv_eff_MC_em.getVal()
    tmpq_eff_data_em_LP = 1. - rrv_eff_data_em.getVal()

    tmpq_eff_MC_em_LP_err   = TMath.Sqrt( rrv_eff_MC_em.getError()  **2)
    tmpq_eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2)

    pureq_wtagger_sf_em_LP = tmpq_eff_data_em_LP / tmpq_eff_MC_em_LP
    pureq_wtagger_sf_em_LP_err = pureq_wtagger_sf_em_LP*TMath.Sqrt( (tmpq_eff_data_em_LP_err/tmpq_eff_data_em_LP)**2 + (tmpq_eff_MC_em_LP_err/tmpq_eff_MC_em_LP)**2 )


    #------------- Print results -------------
    print "--------------------------------------------------------------------------------------------"
    print "                                             HP                                             "
    print "--------------------------------------------------------------------------------------------"
    print "W-tagging SF      : %0.3f +/- %0.3f (rel. error = %0.1f percent)" %(wtagger_sf_em, wtagger_sf_em_err,(wtagger_sf_em_err/wtagger_sf_em*100))
    print "Mass shift [GeV]  : %0.3f +/- %0.3f" %(wtagger_mean_shift_em, wtagger_mean_shift_err_em)
    print "Mass resolution SF: %0.3f +/- %0.3f" %(wtagger_sigma_enlarge_em, wtagger_sigma_enlarge_err_em)
    print ""
    print "Parameter             Data               Simulation        Data/Simulation"
    print " < m >           %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(rrv_mean_data_em.getVal(),rrv_mean_data_em.getError(),rrv_mean_MC_em.getVal(),rrv_mean_MC_em.getError()    , rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal()  ,  mean_sf_error)
    print " #sigma          %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(rrv_sigma_data_em.getVal(),rrv_sigma_data_em.getError(),rrv_sigma_MC_em.getVal(),rrv_sigma_MC_em.getError(), rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal(),  sigma_sf_error)
    print ""
    print "HP W-tag eff+SF  %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(rrv_eff_data_em.getVal(),rrv_eff_data_em.getError(),rrv_eff_MC_em.getVal(),rrv_eff_MC_em.getError()        , rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()    ,  eff_sf_error)
    print "--------------------------------------------------------------------------------------------"
    print "                                        EXTREME FAIL                                        "
    print "--------------------------------------------------------------------------------------------"
    print "Parameter                     Data            Simulation            Data/Simulation"
    print "Extreme fail eff+SF     %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(eff_data_em_extremefail,eff_data_em_extremefail_error,eff_MC_em_extremefail,eff_MC_em_extremefail_error, eff_SF_extremefail, eff_SF_extremefail_error)
    print "--------------------------------------------------------------------------------------------"
    print "                                             LP                                             "
    print "--------------------------------------------------------------------------------------------"
    print "Parameter                               Data              Simulation         Data/Simulation"
    print "LP W-tag eff+SF ( w/ext fail)      %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(eff_data_em_LP,eff_data_em_LP_err,eff_MC_em_LP,eff_MC_em_LP_err,wtagger_sf_em_LP,wtagger_sf_em_LP_err)
    print "LP W-tag eff+SF (wo/ext fail)      %0.3f +/- %0.3f      %0.3f +/- %0.3f      %0.3f +/- %0.3f" %(tmpq_eff_data_em_LP,tmpq_eff_data_em_LP_err,tmpq_eff_MC_em_LP,tmpq_eff_MC_em_LP_err,pureq_wtagger_sf_em_LP,pureq_wtagger_sf_em_LP_err)
    print "--------------------------------------------------------------------------------------------"
    print "--------------------------------------------------------------------------------------------"
    print ""; print "";

    #Write results to file
    fitter.file_out_ttbar_control.write("\n                                     HP                                    ")
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\nW-tagging SF      : %0.3f +/- %0.3f (rel. error = %.1f percent)" %(wtagger_sf_em, wtagger_sf_em_err, (wtagger_sf_em_err/wtagger_sf_em*100)))
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\nMass shift [GeV]  : %0.3f +/- %0.3f" %(wtagger_mean_shift_em, wtagger_mean_shift_err_em))
    fitter.file_out_ttbar_control.write("\nMass resolution SF: %0.3f +/- %0.3f" %(wtagger_sigma_enlarge_em, wtagger_sigma_enlarge_err_em))
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\nParameter                 Data                          Simulation                          Data/Simulation")
    fitter.file_out_ttbar_control.write("\n < m >              %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_mean_data_em.getVal(),rrv_mean_data_em.getError(),rrv_mean_MC_em.getVal(),rrv_mean_MC_em.getError()    , rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal()  ,  mean_sf_error))
    fitter.file_out_ttbar_control.write("\n #sigma             %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_sigma_data_em.getVal(),rrv_sigma_data_em.getError(),rrv_sigma_MC_em.getVal(),rrv_sigma_MC_em.getError(), rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal(),  sigma_sf_error))
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\nHP W-tag eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_eff_data_em.getVal(),rrv_eff_data_em.getError(),rrv_eff_MC_em.getVal(),rrv_eff_MC_em.getError(), rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal(),  eff_sf_error))
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
    fitter.file_out_ttbar_control.write("\n                                EXTREME FAIL                               ")
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
    fitter.file_out_ttbar_control.write("\nParameter                     Data                          Simulation                          Data/Simulation")
    fitter.file_out_ttbar_control.write("\nExtreme fail eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_extremefail,eff_data_em_extremefail_error,eff_MC_em_extremefail,eff_MC_em_extremefail_error, eff_SF_extremefail, eff_SF_extremefail_error))
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
    fitter.file_out_ttbar_control.write("\n                                    LP                                     ")
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\nParameter                                   Data                          Simulation                          Data/Simulation")
    fitter.file_out_ttbar_control.write("\nLP W-tag eff+SF ( w/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_LP,eff_data_em_LP_err,eff_MC_em_LP,eff_MC_em_LP_err,wtagger_sf_em_LP,wtagger_sf_em_LP_err))
    fitter.file_out_ttbar_control.write("\nLP W-tag eff+SF (wo/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(tmpq_eff_data_em_LP,tmpq_eff_data_em_LP_err,tmpq_eff_MC_em_LP,tmpq_eff_MC_em_LP_err,pureq_wtagger_sf_em_LP,pureq_wtagger_sf_em_LP_err))
    fitter.file_out_ttbar_control.write("\n")
    fitter.file_out_ttbar_control.write("\n-----------------------------------------------------------------------------------------------------------------------------")
     
     
     
class doWtagFits:
    def __init__(self):
        self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_")                           # create workspace
        self.boostedW_fitter_em = initialiseFits("em", options.sample, 40, 130, self.workspace4fit_)    # Define all shapes to be used for Mj, define regions (SB,signal) and input files. 
        self.boostedW_fitter_em.get_datasets_fit_minor_bkg()                                            # Loop over intrees to create datasets om Mj and fit the single MCs.
       
        print "Printing workspace:"; self.workspace4fit_.Print(); print ""
        
        self.boostedW_fitter_em.get_sim_fit_components()     

        #Defining categories
        sample_type = RooCategory("sample_type","sample_type")
        sample_type.defineType("em_pass")
        sample_type.defineType("em_fail")

        #Importing fit variables
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
        rrv_weight = RooRealVar("rrv_weight","rrv_weight",0. ,10000000.)
        
        #-------------IMPORT DATA-------------
        #Importing datasets
        rdataset_data_em_mj      = self.workspace4fit_.data("rdataset_data_em_mj")
        rdataset_data_em_mj_fail = self.workspace4fit_.data("rdataset_data_failtau2tau1cut_em_mj")

        #For binned fit (shorter computing time, more presise when no SumW2Error is used!)
        if not options.doUnbinnedFit:
          #Converting to RooDataHist
          rdatahist_data_em_mj      = RooDataHist(rdataset_data_em_mj.binnedClone())
          rdatahist_data_em_mj_fail = RooDataHist(rdataset_data_em_mj_fail.binnedClone())

          #Converting back to RooDataSet
          rdataset_data_em_mj_2 = rdataset_data_em_mj.emptyClone()
          for i in range(0,rdatahist_data_em_mj.numEntries()):
            rdataset_data_em_mj_2.add(rdatahist_data_em_mj.get(i),rdatahist_data_em_mj.weight())

          rdataset_data_em_mj_fail_2 = rdataset_data_em_mj_fail.emptyClone()
          for i in range(0,rdatahist_data_em_mj_fail.numEntries()):
            rdataset_data_em_mj_fail_2.add(rdatahist_data_em_mj_fail.get(i),rdatahist_data_em_mj_fail.weight())

          #Combined dataset
          combData_data = RooDataSet("combData_data","combData_data",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_data_em_mj_2),RooFit.Import("em_fail",rdataset_data_em_mj_fail_2) )

        #For unbinned fit
        else:
          combData_data = RooDataSet("combData_data","combData_data",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_data_em_mj),RooFit.Import("em_fail",rdataset_data_em_mj_fail) )

        #-------------IMPORT MC-------------
        #Importing MC datasets
        rdataset_TotalMC_em_mj      = self.workspace4fit_.data("rdataset_TotalMC_em_mj")
        rdataset_TotalMC_em_mj_fail = self.workspace4fit_.data("rdataset_TotalMC_failtau2tau1cut_em_mj")

        if not options.doUnbinnedFit:
          #Converting to RooDataHist
          rdatahist_TotalMC_em_mj      = RooDataHist(rdataset_TotalMC_em_mj.binnedClone())
          rdatahist_TotalMC_em_mj_fail = RooDataHist(rdataset_TotalMC_em_mj_fail.binnedClone())

          #Converting back to RooDataSet
          rdataset_TotalMC_em_mj_2 = rdataset_TotalMC_em_mj.emptyClone()
          for i in range(0,rdatahist_TotalMC_em_mj.numEntries()):
            rdataset_TotalMC_em_mj_2.add(rdatahist_TotalMC_em_mj.get(i),rdatahist_TotalMC_em_mj.weight())

          rdataset_TotalMC_em_mj_fail_2 = rdataset_TotalMC_em_mj_fail.emptyClone()
          for i in range(0,rdatahist_TotalMC_em_mj_fail.numEntries()):
            rdataset_TotalMC_em_mj_fail_2.add(rdatahist_TotalMC_em_mj_fail.get(i),rdatahist_TotalMC_em_mj_fail.weight())

          #Combined MC dataset
          combData_TotalMC = RooDataSet("combData_TotalMC","combData_TotalMC",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_TotalMC_em_mj_2),RooFit.Import("em_fail",rdataset_TotalMC_em_mj_fail_2) )

        else:
         combData_TotalMC = RooDataSet("combData_TotalMC","combData_TotalMC",RooArgSet(rrv_mass_j,rrv_weight),RooFit.WeightVar(rrv_weight),RooFit.Index(sample_type),RooFit.Import("em_pass",rdataset_TotalMC_em_mj),RooFit.Import("em_fail",rdataset_TotalMC_em_mj_fail) )

        #-------------Define and perform fit to data-------------
        #Import pdf from single fits and define the simultaneous total pdf
        model_data_em      = self.workspace4fit_.pdf("model_data_em")
        model_data_fail_em = self.workspace4fit_.pdf("model_data_failtau2tau1cut_em")

        simPdf_data = RooSimultaneous("simPdf_data_em","simPdf_data_em",sample_type)
        simPdf_data.addPdf(model_data_em,"em_pass")
        simPdf_data.addPdf(model_data_fail_em,"em_fail")

        #Import Gaussian constraints to propagate error to likelihood
        constrainslist_data_em = ROOT.std.vector(ROOT.std.string)()
        for i in range(self.boostedW_fitter_em.constrainslist_data.size()):
            constrainslist_data_em.push_back(self.boostedW_fitter_em.constrainslist_data.at(i))
            print self.boostedW_fitter_em.constrainslist_data.at(i)
        pdfconstrainslist_data_em = RooArgSet("pdfconstrainslist_data_em")
        for i in range(constrainslist_data_em.size()):
          pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em.at(i)) )
          pdfconstrainslist_data_em.Print()

        # Perform simoultaneous fit to data
        if not options.doUnbinnedFit:
          rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_data_em))#, RooFit.SumW2Error(kTRUE))
          rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_data_em))#, RooFit.SumW2Error(kTRUE))
        else:
          rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_data_em), RooFit.SumW2Error(kTRUE))
          rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_data_em), RooFit.SumW2Error(kTRUE))

        #Draw       
        chi2FailData = drawFrameGetChi2(rrv_mass_j,rfresult_data,rdataset_data_em_mj_fail,model_data_fail_em)
        chi2PassData = drawFrameGetChi2(rrv_mass_j,rfresult_data,rdataset_data_em_mj,model_data_em)

        #Print final data fit results
        print "FIT parameters (DATA) :"; print ""
        print "CHI2 PASS = %.3f    CHI2 FAIL = %.3f" %(chi2PassData,chi2FailData)
        print ""; print rfresult_data.Print(); print ""

        #-------------Define and perform fit to MC-------------

        # fit TotalMC --> define the simultaneous total pdf
        model_TotalMC_em      = self.workspace4fit_.pdf("model_TotalMC_em")
        model_TotalMC_fail_em = self.workspace4fit_.pdf("model_TotalMC_failtau2tau1cut_em")
        simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em","simPdf_TotalMC_em",sample_type)
        simPdf_TotalMC.addPdf(model_TotalMC_em,"em_pass")
        simPdf_TotalMC.addPdf(model_TotalMC_fail_em,"em_fail")

        #Import Gaussian constraints  for fixed paramters to propagate error to likelihood
        constrainslist_TotalMC_em = ROOT.std.vector(ROOT.std.string)()
        for i in range(self.boostedW_fitter_em.constrainslist_mc.size()):
            constrainslist_TotalMC_em.push_back(self.boostedW_fitter_em.constrainslist_mc.at(i))
        pdfconstrainslist_TotalMC_em = RooArgSet("pdfconstrainslist_TotalMC_em")
        for i in range(constrainslist_TotalMC_em.size()):
          pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) )

        # Perform simoultaneous fit to MC
        if not options.doUnbinnedFit:
          rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em))#, RooFit.SumW2Error(kTRUE))--> Removing due to unexected behaviour. See https://root.cern.ch/phpBB3/viewtopic.php?t=16917, https://root.cern.ch/phpBB3/viewtopic.php?t=16917
          rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em))#, RooFit.SumW2Error(kTRUE))        
        else:
          rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em), RooFit.SumW2Error(kTRUE))
          rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.Verbose(kFALSE), RooFit.Minimizer("Minuit2"),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em), RooFit.SumW2Error(kTRUE))
          
          
        chi2FailMC = drawFrameGetChi2(rrv_mass_j,rfresult_TotalMC,rdataset_TotalMC_em_mj_fail,model_TotalMC_fail_em)
        chi2PassMC = drawFrameGetChi2(rrv_mass_j,rfresult_TotalMC,rdataset_TotalMC_em_mj,model_TotalMC_em)
        
        #Print final MC fit results
        print "FIT Par. (MC) :"; print ""
        print "CHI2 PASS = %.3f    CHI2 FAIL = %.3f" %(chi2PassMC,chi2FailMC)
        print ""; print rfresult_TotalMC.Print(); print ""
        
        # draw the final fit results
        DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_em.color_palet,"","em",self.boostedW_fitter_em.wtagger_label,self.boostedW_fitter_em.AK8_pt_min,self.boostedW_fitter_em.AK8_pt_max)
       
        # Get W-tagging scalefactor and efficiencies
        GetWtagScalefactors(self.workspace4fit_,self.boostedW_fitter_em)
        
        # Delete workspace
        del self.workspace4fit_

class initialiseFits:

    # Constructor: Input is channel (mu,ele,em), range in mj and a workspace
    def __init__(self, in_channel, in_sample, in_mj_min=40, in_mj_max=130, input_workspace=None):
      
      RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
      RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)
      
      # Set channel 
      self.channel = in_channel
      
      print "CHANNEL = %s" %in_channel
      print "Using Tau21 HP cut of " ,options.tau2tau1cutHP 

      # Map of shapes to be used for the various fits (defined in PDFs/MakePdf.cxx)                                                                                                                                         
      self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string)()
      
      # Fit functions for matched tt MC
      self.mj_shape["TTbar_realW"]      = "GausErfExp_ttbar" #before "2Gaus_ttbar"
      self.mj_shape["TTbar_realW_fail"] = "GausErfExp_ttbar_failtau2tau1cut" #before "GausChebychev_ttbar_failtau2tau1cut"
      self.mj_shape["TTbar_fakeW"]      = "ErfExp_ttbar"
      self.mj_shape["TTbar_fakeW_fail"] = "ErfExp_ttbar_failtau2tau1cut"
      
      # Fit functions for minor backgrounds
      self.mj_shape["VV"]                 = "ExpGaus"
      self.mj_shape["VV_fail"]            = "ExpGaus"
      self.mj_shape["WJets0"]             = "ErfExp"
      self.mj_shape["WJets0_fail"]        = "ErfExp"
      self.mj_shape["STop"]               = "ErfExpGaus_sp"       
      self.mj_shape["STop_fail"]          = "ExpGaus"  
          
      if (options.tau2tau1cutHP==0.60): self.mj_shape["VV"]                 = "ErfExpGaus_sp"
        
      if (options.useDDT):
        self.mj_shape["VV"]                 = "ExpGaus"
        self.mj_shape["VV_fail"]            = "ErfExpGaus_sp"
        self.mj_shape["STop_fail"]          = "ErfExpGaus_sp" 
        self.mj_shape["STop"]               = "ExpGaus"  
        
      # Fit functions used in simultaneous fit of pass and fail categories
      self.mj_shape["bkg_mc_fail"]          = "ErfExp_ttbar_failtau2tau1cut"
      self.mj_shape["bkg_data_fail"]        = "ErfExp_ttbar_failtau2tau1cut"    
      
      self.mj_shape["signal_mc_fail"]       = "GausErfExp_ttbar_failtau2tau1cut" #Before GausChebychev_ttbar_failtau2tau1cut
      self.mj_shape["signal_data_fail"]     = "GausErfExp_ttbar_failtau2tau1cut"
         
      self.mj_shape["bkg_data"]             = "ErfExp_ttbar" 
      self.mj_shape["bkg_mc"]               = "ErfExp_ttbar" 
      
      self.mj_shape["signal_data"]          = "GausErfExp_ttbar" #Before 2Gaus_ttbar
      self.mj_shape["signal_mc"]            = "GausErfExp_ttbar"
      
      if (options.useDDT): 
        self.mj_shape["signal_mc_fail"]       = "GausChebychev_ttbar_failtau2tau1cut" 
        self.mj_shape["signal_data_fail"]     = "GausChebychev_ttbar_failtau2tau1cut"
      
      # #TESTS USING OLD METHOD, EG ADDING TWO ADDITIONAL FIT FUNCTIONS
      # self.mj_shape["signal_data"]          = "2Gaus_ttbar"
      # self.mj_shape["signal_mc"]            = "2Gaus_ttbar"
      # self.mj_shape["signal_mc_fail"]       = "GausChebychev_ttbar_failtau2tau1cut"
      # self.mj_shape["signal_data_fail"]     = "GausChebychev_ttbar_failtau2tau1cut"
      
      #Set lumi  
      self.Lumi=2198. #74
      if options.use76X: self.Lumi=2300. #76
          
      self.BinWidth_mj = 5.
      self.narrow_factor = 1.

      self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor
      nbins_mj         = int( (in_mj_max - in_mj_min) / self.BinWidth_mj )
      in_mj_max        = in_mj_min+nbins_mj*self.BinWidth_mj
      
      jetMass = "pruned jet mass"
      if options.usePuppiSD: jetMass = "PUPPI softdrop jet mass"

      rrv_mass_j = RooRealVar("rrv_mass_j", jetMass ,(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV")
      rrv_mass_j.setBins(nbins_mj)
 
      # Create workspace and import fit variable
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
        

      # Directory and input files
      self.file_Directory         = "$HOME/EXOVVAnalysisRunII/AnalysisOutput/Wtag/PRUNED/WWTree_%s/"%(self.channel)
    
      # if options.usePuppiSD or options.useDDT:
      #   self.file_Directory       = "$HOME/EXOVVAnalysisRunII/AnalysisOutput/Wtag/PUPPISD/WWTree_%s/"%(self.channel)
      #   options.use76X = True

      postfix = ""
      if options.use76X: postfix ="_76X"  
          
      self.file_data              = ("ExoDiBosonAnalysis.WWTree_data%s.root") %postfix        
      self.file_WJets0_mc         = ("ExoDiBosonAnalysis.WWTree_WJets%s.root") %postfix
      self.file_VV_mc             = ("ExoDiBosonAnalysis.WWTree_VV%s.root") %postfix      
      self.file_STop_mc           = ("ExoDiBosonAnalysis.WWTree_STop%s.root") %postfix           
      self.file_TTbar_mc          = ("ExoDiBosonAnalysis.WWTree_TTbar_%s%s.root")%(in_sample,postfix)
      self.file_pseudodata        = ("ExoDiBosonAnalysis.WWTree_pseudodata_%s%s.root")%(in_sample,postfix)      #Important! ROOT tree containing all backgrounds added together (tt+singleT+VV+Wjets). Used for fit to total MC    
             
      # Define Tau21 WP
      self.wtagger_label = "HP"
      self.wtagger_cut = options.tau2tau1cutHP
      self.wtagger_cut_min = 0.
      
      if options.usePuppiSD: 
        postfix = postfix + "_PuppiSD"
      if options.useDDT: 
        postfix = postfix + "_PuppiSD_DDT"    
      
      # Define label used for plots and choosing fit paramters in PDFs/MakePdf.cxx  
      wp = "%.2f" %options.tau2tau1cutHP
      wp = wp.replace(".","v")
      self.wtagger_label = self.wtagger_label + "%s%s%s"%(wp,in_sample,postfix) 

      
      #Color pallett for plots
      self.color_palet = ROOT.std.map(ROOT.std.string, int) ()
      self.color_palet["data"]              = 1
      self.color_palet["WJets"]             = 2
      self.color_palet["VV"]                = 4
      self.color_palet["STop"]              = 7
      self.color_palet["TTbar"]             = 210
      self.color_palet["Signal"]            = 1
      self.color_palet["Uncertainty"]       = 1
      self.color_palet["Other_Backgrounds"] = 1    
      
      # Cuts (dont need these cuts if they are already implemented in ROOT tree)
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
      
      # Out .txt file with final SF numbers
      self.file_ttbar_control_txt = "WtaggingSF.txt"
      self.file_out_ttbar_control = open(self.file_ttbar_control_txt,"w")
                                                                                                                                                             
      setTDRStyle()

    def get_datasets_fit_minor_bkg(self):
        
        rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
  
        # Build single-t fit pass and fail distributions
        print "##################################################"
        print "############### Single Top DataSet ###############"
        print "##################################################"
        print ""

        self.get_mj_dataset(self.file_STop_mc,"_STop")
        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop"                        ,self.mj_shape["STop"],self.channel,self.wtagger_label) #Start value and range of parameters defined in PDFs/MakePDF.cxx
        fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_failtau2tau1cut"        ,self.mj_shape["STop_fail"],self.channel,self.wtagger_label)

        ### Build WJet fit pass and fail distributions
        print "###########################################"
        print "############### WJets Pythia ##############"
        print "###########################################"
        print ""

        self.get_mj_dataset(self.file_WJets0_mc,"_WJets0")
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0",self.mj_shape["WJets0"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_failtau2tau1cut",self.mj_shape["WJets0_fail"],self.channel,self.wtagger_label)


        # Build VV fit pass and fail distributions
        print "#########################################"
        print "############### VV Pythia ###############"
        print "#########################################"
        print ""
        print ""

        self.get_mj_dataset(self.file_VV_mc,"_VV")
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label)
        fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_failtau2tau1cut",self.mj_shape["VV_fail"],self.channel,self.wtagger_label)

        if options.fitMC:
            return
        
        # Get dataset for ttbar    
        print "#########################################"
        print "################ TTbar %s################"%options.sample
        print "#########################################"
        print ""

        self.get_mj_dataset(self.file_TTbar_mc,"_TTbar")    

        # Get dataset used for fit to total MC
        print "################################################"
        print "############## Pseudo Data Powheg ##############"
        print "################################################"
        print ""
        self.get_mj_dataset(self.file_pseudodata,"_TotalMC")

        print "#################################"
        print "############# Data ##############"
        print "#################################"
        print ""
        print ""
        self.get_mj_dataset(self.file_data,"_data")

        # print "Saving workspace in myworkspace.root! To save time when debugging use option --WS myworkspace.root to avoid recreating workspace every time"
        # filename = "myworkspace.root"
        # self.workspace4fit_.writeToFile(filename)

    def get_sim_fit_components(self):
      # self.print_yields()
      self.constrainslist_data = ROOT.std.vector(ROOT.std.string)()
      self.constrainslist_mc   = ROOT.std.vector(ROOT.std.string)()
        
      #Construct pass/fail models (fix minor backgrounds, create sim. fit total PDFS)
      ScaleFactorTTbarControlSampleFit(self.workspace4fit_,self.mj_shape,self.color_palet,self.constrainslist_data,self.constrainslist_mc,"",self.channel,self.wtagger_label,self.AK8_pt_min,self.AK8_pt_max)
     
      #Get data/MC scalefactors
      rrv_scale_number                      = self.workspace4fit_.var("rrv_scale_number_TTbar_STop_VV_WJets").getVal()
      rrv_scale_number_fail                 = self.workspace4fit_.var("rrv_scale_number_TTbar_STop_VV_WJets_fail").getVal()
      
      #Print data/MC scalefactors
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

      if options.usePuppiSD or options.useDDT: 
        jet_mass="Whadr_puppi_softdrop"
      
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
      i = 0
      for i in range(treeIn.GetEntries()):
          if i % 5000 == 0: print "iEntry: ",i
          treeIn.GetEntry(i)
          
          if TString(label).Contains("realW") and not getattr(treeIn,"Whadr_isW_def1"): 
            continue
          if TString(label).Contains("fakeW") and     getattr(treeIn,"Whadr_isW_def1"): 
            continue
          
          wtagger = getattr(treeIn,"Whadr_tau21")
          if options.usePuppiSD:
            wtagger = getattr(treeIn,"Whadr_puppi_tau2")/getattr(treeIn,"Whadr_puppi_tau1")
          if options.useDDT:
            wtagger = getattr(treeIn,"Whadr_puppi_tau2")/getattr(treeIn,"Whadr_puppi_tau1")+ (0.063 * TMath.log( (pow( getattr(treeIn,"Whadr_puppi_softdrop"),2))/getattr(treeIn,"Whadr_puppi_pt") ))        
            
          discriminantCut = 0
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
            tmp_event_weight4fit = getattr(treeIn,"genweight")*getattr(treeIn,"puweight")*getattr(treeIn,"hltweight") #Missing b-tag weight and Vtg weight! ##FROM JEN!
            tmp_event_weight4fit = tmp_event_weight4fit*treeIn.lumiweight*self.Lumi /tmp_scale_to_lumi      ##FROM JEN!
         
         
            # tmp_event_weight4fit = getattr(treeIn,"weight")*self.Lumi/tmp_scale_to_lumi ##FROM QUN
            
            # tmp_event_weight = getattr(treeIn,"weight")*self.Lumi  ##LUCA
            # tmp_event_weight4fit = getattr(treeIn,"genweight")*getattr(treeIn,"puweight")*getattr(treeIn,"hltweight")  ##LUCA
            # tmp_event_weight4fit = tmp_event_weight/tmp_scale_to_lumi
          else:
            tmp_event_weight = 1.
            tmp_event_weight4fit = 1.

          
          if options.fitTT:
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
    channel = options.channel ## ele, mu or ele+mu combined
    if options.fitTT:
        print "Doing fits to matched tt MC. Tree must contain branch with flag for match/no-match to generator level W!"
        doFitsToMatchedTT()
    elif options.fitMC:
        print "Doing fits to MC only"
        doFitsToMC()
    else:
        print 'Getting W-tagging scalefactor for %s sample for n-subjettiness < %.2f' %(channel,options.tau2tau1cutHP) #I am actually not doing a simoultaneous fit. So..... change this
        getSF()

