from optparse import OptionParser
import ROOT, sys


parser = OptionParser()
parser.add_option('-c', '--channel',action="store",type="string",dest="channel",default="mu")
parser.add_option('-s','--simultaneous', action='store_true', dest='fitwtaggersim', default=False, help='Fit W-tagger adding mu+ele channels')
parser.add_option('--category', action="store",type="string",dest="category",default="HP")
parser.add_option('--tau2tau1cutHP', action="store", type="float",dest="tau2tau1cutHP",default=0.60)
parser.add_option('--tau2tau1cutLP', action="store", type="float",dest="tau2tau1cutLP",default=0.75)

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
gInterpreter.GenerateDictionary("std::vector<std::string>", "vector;string")

### Fit e+mu together
def control_sample_simultaneous():

    print "control_sample_simultaneous"
    boostedW_fitter_sim = doFit_wj_and_wlvj_simultaneous()

def control_sample(channel="mu"):

    print "control_sample"
    
class doFit_wj_and_wlvj_simultaneous:
  
    def __init__(self):
      
      self.workspace4fit_ = RooWorkspace("workspace4fit_","workspace4fit_")           # create workspace
      
      self.boostedW_fitter_em = doFit_wj_and_wlvj("em", 40, 130, self.workspace4fit_) # Define all shapes to be used for Mj, define regions (SB,signal) and input files.
                                                              
      self.boostedW_fitter_em.fit_TTbar_controlsample()                               # Loop over intrees to create datasets om Mj and fit the single MCs.
 
      self.workspace4fit_.data("rdataset_data_em_mj").Print() # Then I PRINT the contents of my workspace

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
      combData_data.Print()

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

      label = ""

      constrainslist_data_em = ROOT.std.vector(ROOT.std.string)()
      for i in range(self.boostedW_fitter_em.constrainslist_data.size()):
          constrainslist_data_em.push_back(self.boostedW_fitter_em.constrainslist_data.at(i))

      pdfconstrainslist_data_em = RooArgSet("pdfconstrainslist_data_em"+label)
      for i in range(constrainslist_data_em.size()):
          self.workspace4fit_.pdf(constrainslist_data_em.at(i)).Print()
          pdfconstrainslist_data_em.add(self.workspace4fit_.pdf(constrainslist_data_em.at(i)) )
      pdfconstrainslist_data_em.Print()

      rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em))
      rfresult_data = simPdf_data.fitTo(combData_data,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_data_em))

      # fit TotalMC --> define the simultaneous total pdf
      model_TotalMC_em      = self.workspace4fit_.pdf("model_TotalMC"+label+"_em")
      model_TotalMC_fail_em = self.workspace4fit_.pdf("model_TotalMC"+label+"_failtau2tau1cut_em")

      simPdf_TotalMC = RooSimultaneous("simPdf_TotalMC_em"+label,"simPdf_TotalMC_em"+label,sample_type)
      simPdf_TotalMC.addPdf(model_TotalMC_em,"em_pass")
      simPdf_TotalMC.addPdf(model_TotalMC_fail_em,"em_fail")

      constrainslist_TotalMC_em = ROOT.std.vector(ROOT.std.string)()
      for i in range(self.boostedW_fitter_em.constrainslist_mc.size()):
          constrainslist_TotalMC_em.push_back(self.boostedW_fitter_em.constrainslist_mc.at(i))

      pdfconstrainslist_TotalMC_em = RooArgSet("pdfconstrainslist_TotalMC_em")
      for i in range(constrainslist_TotalMC_em.size()):
          self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]).Print()
          pdfconstrainslist_TotalMC_em.add(self.workspace4fit_.pdf(constrainslist_TotalMC_em[i]) )
      pdfconstrainslist_TotalMC_em.Print()

      rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em))
      rfresult_TotalMC = simPdf_TotalMC.fitTo(combData_TotalMC,RooFit.Save(kTRUE),RooFit.ExternalConstraints(pdfconstrainslist_TotalMC_em))

      ## draw the plots
      DrawScaleFactorTTbarControlSample(self.workspace4fit_,self.boostedW_fitter_em.color_palet,label,"em",self.boostedW_fitter_em.wtagger_label,self.boostedW_fitter_em.AK8_pt_min,self.boostedW_fitter_em.AK8_pt_max)

      rfresult_TotalMC.Print()
      rfresult_data.Print()

      ### Efficiency in data and MC
      rrv_eff_MC_em   = self.workspace4fit_.var("eff_ttbar_TotalMC"+label+"_em_mj")
      rrv_mean_MC_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_TotalMC"+label+"_em_mj")
      rrv_sigma_MC_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_em_mj")

      rrv_eff_data_em   = self.workspace4fit_.var("eff_ttbar_data"+label+"_em_mj")
      rrv_mean_data_em  = self.workspace4fit_.var("rrv_mean1_gaus_ttbar_data"+label+"_em_mj")
      rrv_sigma_data_em = self.workspace4fit_.var("rrv_sigma1_gaus_ttbar_data"+label+"_em_mj")

      # rrv_eff_MC_em.Print()
 #      rrv_eff_MC_em.Print()
 #      rrv_eff_data_em.Print()
 #      rrv_eff_data_em.Print()
 #      rrv_mean_MC_em.Print()
 #      rrv_mean_data_em.Print()
 #      rrv_sigma_MC_em.Print()
 #      rrv_sigma_data_em.Print()

      ## GET HP SCALEFACTOR AND UNCERTIANTIES
      pure_wtagger_sf_em            = rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()
      
      pure_wtagger_mean_shift_em    = rrv_mean_data_em.getVal()-rrv_mean_MC_em.getVal()
      pure_wtagger_sigma_enlarge_em = rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()

      pure_wtagger_sf_em_err = ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal() )**2 + (rrv_eff_MC_em.getError()/rrv_eff_MC_em.getVal())**2 )**0.5* pure_wtagger_sf_em

      pure_wtagger_mean_shift_err_em = (rrv_mean_data_em.getError()**2 + rrv_mean_MC_em.getError()**2)**0.5

      pure_wtagger_sigma_enlarge_err_em = ((rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 + (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2 )**0.5* pure_wtagger_sigma_enlarge_em
      
      mean_sf_error  = (rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal())   * ( (rrv_mean_data_em.getError()/rrv_mean_data_em.getVal())**2   +  (rrv_mean_MC_em.getError() /rrv_mean_MC_em.getVal())**2   )**0.5
      sigma_sf_error = (rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal()) * ( (rrv_sigma_data_em.getError()/rrv_sigma_data_em.getVal())**2 +  (rrv_sigma_MC_em.getError()/rrv_sigma_MC_em.getVal())**2 )**0.5
      eff_sf_error  = (rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal())      * ( (rrv_eff_data_em.getError()/rrv_eff_data_em.getVal())**2     +  (rrv_eff_MC_em.getError()  /rrv_eff_MC_em.getVal())**2     )**0.5

      print "                                     HP                                    "
      print "---------------------------------------------------------------------------"
      print ""
      print "Pure W-tagging SF of em %s         : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em, pure_wtagger_sf_em_err)
      print ""
      print "Pure W-tagging mean shift em %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em)
      print "Pure W-tagging sigma enlarge em %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em)
      print ""
      print "Parameter                 Data                          Simulation                          Data/Simulation"
      print " < m >              %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_mean_data_em.getVal(),rrv_mean_data_em.getError(),rrv_mean_MC_em.getVal(),rrv_mean_MC_em.getError()    , rrv_mean_data_em.getVal()/rrv_mean_MC_em.getVal()  ,  mean_sf_error)
      print " #sigma             %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_sigma_data_em.getVal(),rrv_sigma_data_em.getError(),rrv_sigma_MC_em.getVal(),rrv_sigma_MC_em.getError(), rrv_sigma_data_em.getVal()/rrv_sigma_MC_em.getVal(),  sigma_sf_error)
      print ""
      print ""
      print "HP W-tag eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(rrv_eff_data_em.getVal(),rrv_eff_data_em.getError(),rrv_eff_MC_em.getVal(),rrv_eff_MC_em.getError()        , rrv_eff_data_em.getVal()/rrv_eff_MC_em.getVal()    ,  eff_sf_error)
      print ""
      print "---------------------------------------------------------------------------"
      
      self.boostedW_fitter_em.file_out_ttbar_control.write( "\nPure W-tagger SF of em %s     : %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em, pure_wtagger_sf_em_err))
      self.boostedW_fitter_em.file_out_ttbar_control.write( "\nPure W-tagger mean shift em %s    : %0.3f +/- %0.3f"%(label,pure_wtagger_mean_shift_em, pure_wtagger_mean_shift_err_em))
      self.boostedW_fitter_em.file_out_ttbar_control.write( "\nPure W-tagger sigma enlarge em %s : %0.3f +/- %0.3f"%(label,pure_wtagger_sigma_enlarge_em, pure_wtagger_sigma_enlarge_err_em))


      ## GET EXTREME FAIL NUMBERS IN ORDER TO COMPUTE LP SF:   
      rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_extremefailtau2tau1cut_em_mj")
      rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj    = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_extremefailtau2tau1cut_em_mj")
 
      # rrv_number_ttbar_TotalMC_failtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_failtau2tau1cut_em_mj")
      # rrv_number_ttbar_data_failtau2tau1cut_em_mj    = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_failtau2tau1cut_em_mj")
      # rrv_number_ttbar_TotalMC_beforetau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_beforetau2tau1cut_em_mj")
      # rrv_number_ttbar_data_beforetau2tau1cut_em_mj    = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_beforetau2tau1cut_em_mj")
      # rrv_number_ttbar_TotalMC_passtau2tau1cut_em_mj = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_passtau2tau1cut_em_mj")
      # rrv_number_ttbar_data_passtau2tau1cut_em_mj    = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_passtau2tau1cut_em_mj")

      # rrv_number_total_ttbar_TotalMC_em = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_em_mj")
      # rrv_number_total_ttbar_data_em    = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_em_mj")
      rrv_number_total_ttbar_TotalMC_em = self.workspace4fit_.var("rrv_number_ttbar_TotalMC"+label+"_beforetau2tau1cut_em_mj")
      rrv_number_total_ttbar_data_em    = self.workspace4fit_.var("rrv_number_ttbar_data"+label+"_beforetau2tau1cut_em_mj")
      
      
      
      
      
      # rrv_number_total_ttbar_TotalMC_el = self.workspace4fit_.var("rrv_number_total_ttbar_TotalMC"+label+"_el");
      # rrv_number_total_ttbar_data_el = self.workspace4fit_.var("rrv_number_total_ttbar_data"+label+"_el");
      # print "el TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal());
      # print "el data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal());
      #
      # tmp_eff_MC_el_extremefail = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_TotalMC_el.getVal();
      # tmp_eff_data_el_extremefail= rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() / rrv_number_total_ttbar_data_el.getVal();
      #
      # tmp_eff_MC_el_extremefail_error =tmp_eff_MC_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_TotalMC_el.getError()/rrv_number_total_ttbar_TotalMC_el.getVal() )**2 );
      # tmp_eff_data_el_extremefail_error =tmp_eff_data_el_extremefail* TMath.Sqrt( (rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_el_mj.getVal() )**2+ (rrv_number_total_ttbar_data_el.getError()/rrv_number_total_ttbar_data_el.getVal() )**2 );
      #
      #
      # tmp_eff_MC_el_LP =1. - rrv_eff_MC_el.getVal() - tmp_eff_MC_el_extremefail;
      # tmp_eff_data_el_LP =1. - rrv_eff_data_el.getVal() - tmp_eff_data_el_extremefail;
      # tmp_eff_MC_el_LP_err = TMath.Sqrt( rrv_eff_MC_el.getError()**2 + tmp_eff_MC_el_extremefail_error**2 );
      # tmp_eff_data_el_LP_err = TMath.Sqrt( rrv_eff_data_el.getError()**2 + tmp_eff_data_el_extremefail_error**2 );
      #
      # tmpq_eff_MC_el_LP =1. - rrv_eff_MC_el.getVal() ;
      # tmpq_eff_data_el_LP =1. - rrv_eff_data_el.getVal() ;
      #
      # tmpq_eff_MC_el_LP_err = TMath.Sqrt( rrv_eff_MC_el.getError()**2);# + tmpq_eff_MC_el_extremefail_error**2 );
      # tmpq_eff_data_el_LP_err = TMath.Sqrt( rrv_eff_data_el.getError()**2);# + tmpq_eff_data_el_extremefail_error**2 );
      # print "LP Eff of el data %s: %0.3f +/- %0.3f"%(label,tmp_eff_data_el_LP, tmp_eff_data_el_LP_err);
      # print "LP Eff of el MC %s: %0.3f +/- %0.3f"%(label,tmp_eff_MC_el_LP, tmp_eff_MC_el_LP_err);
      #
      # pure_wtagger_sf_el_LP = tmp_eff_data_el_LP / tmp_eff_MC_el_LP;
      #  pure_wtagger_sf_el_LP_err = pure_wtagger_sf_el_LP*TMath.Sqrt( (tmp_eff_data_el_LP_err/tmp_eff_data_el_LP)**2 + (tmp_eff_MC_el_LP_err/tmp_eff_MC_el_LP)**2 );
      #
      #  pureq_wtagger_sf_el_LP = tmpq_eff_data_el_LP / tmpq_eff_MC_el_LP;
      #  pureq_wtagger_sf_el_LP_err = pureq_wtagger_sf_el_LP*TMath.Sqrt( (tmpq_eff_data_el_LP_err/tmpq_eff_data_el_LP)**2 + (tmpq_eff_MC_el_LP_err/tmpq_eff_MC_el_LP)**2 );







      # print ""
      # print "em TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal())
      # print "em data Eff of extremefail %s   : %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()    / rrv_number_total_ttbar_data_em.getVal())
      # print ""
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data Eff of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_data_em.getVal()))
      #
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC number of extremefail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data number of extremefail %s: %0.6f"%(label, rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()))
      #
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC number of fail %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_failtau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data number of fail %s: %0.6f"%(label, rrv_number_ttbar_data_failtau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC number of before %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_beforetau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data number of before %s: %0.6f"%(label, rrv_number_ttbar_data_beforetau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC number of pass %s: %0.6f"%(label, rrv_number_ttbar_TotalMC_passtau2tau1cut_em_mj.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data number of pass %s: %0.6f"%(label, rrv_number_ttbar_data_passtau2tau1cut_em_mj.getVal()))
      #
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel TotalMC number of total %s: %0.6f"%(label, rrv_number_total_ttbar_TotalMC_em.getVal()))
      # self.boostedW_fitter_em.file_out_ttbar_control.write("\nel data number of total %s: %0.6f"%(label, rrv_number_total_ttbar_data_em.getVal()))

      eff_MC_em_extremefail   = rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() / rrv_number_total_ttbar_TotalMC_em.getVal()
      eff_data_em_extremefail = rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()    / rrv_number_total_ttbar_data_em.getVal()

      eff_MC_em_extremefail_error   = eff_MC_em_extremefail   * ( (rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_TotalMC_extremefailtau2tau1cut_em_mj.getVal() )**2 + (rrv_number_total_ttbar_TotalMC_em.getError()/rrv_number_total_ttbar_TotalMC_em.getVal())**2 )**0.5
      eff_data_em_extremefail_error = eff_data_em_extremefail * ( (rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getError()/rrv_number_ttbar_data_extremefailtau2tau1cut_em_mj.getVal()       )**2 + (rrv_number_total_ttbar_data_em.getError()   /rrv_number_total_ttbar_data_em.getVal()   )**2 )**0.5

      eff_SF_extremefail       = eff_data_em_extremefail/eff_MC_em_extremefail
      eff_SF_extremefail_error = eff_SF_extremefail * ( (eff_data_em_extremefail_error/eff_data_em_extremefail)**2 + (eff_MC_em_extremefail_error/eff_MC_em_extremefail)**2 )**0.5
      
      # print ""
      # print ""
      # print ""
      # print "eff_MC_em_extremefail_error   %s: %f"%(label,eff_MC_em_extremefail_error)
      # print "eff_data_em_extremefail_error %s: %f"%(label,eff_data_em_extremefail_error)
      # print ""
      # print ""
      # print ""
      
      
      print "                                EXTREME FAIL                               "
      print "---------------------------------------------------------------------------"
      print ""
      print "Parameter                     Data                          Simulation                          Data/Simulation"
      print "Extreme fail eff+SF     %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_extremefail,eff_data_em_extremefail_error,eff_MC_em_extremefail,eff_MC_em_extremefail_error, eff_SF_extremefail, eff_SF_extremefail_error)
      print ""
      print "---------------------------------------------------------------------------"
      
      

      self.boostedW_fitter_em.file_out_ttbar_control.write("\neff_MC_em_extremefail_error %s: %f"%(label,eff_MC_em_extremefail_error))
      self.boostedW_fitter_em.file_out_ttbar_control.write("\neff_data_em_extremefail_error %s: %f"%(label,eff_data_em_extremefail_error))

      self.boostedW_fitter_em.file_out_ttbar_control.write("\nTotalMC Eff of pass %s: %0.6f"%(label, rrv_eff_MC_em.getVal()))
      self.boostedW_fitter_em.file_out_ttbar_control.write("\ndata Eff of pass %s: %0.6f"%(label, rrv_eff_data_em.getVal()))

      self.boostedW_fitter_em.file_out_ttbar_control.write("\nerr TotalMC Eff of pass %s: %0.6f"%(label, rrv_eff_MC_em.getError()))
      self.boostedW_fitter_em.file_out_ttbar_control.write("\nerr data Eff of pass %s: %0.6f"%(label, rrv_eff_data_em.getError()))


      ## GET LP SCALEFACTOR AND UNCERTIANTIES
      eff_MC_em_LP   = 1.-rrv_eff_MC_em.getVal()   - eff_MC_em_extremefail
      eff_data_em_LP = 1.-rrv_eff_data_em.getVal() - eff_data_em_extremefail

      eff_MC_em_LP_err   = TMath.Sqrt( rrv_eff_MC_em.getError()  **2 + eff_MC_em_extremefail_error**2 )
      eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2 + eff_data_em_extremefail_error**2 )
      
      
      
      tmpq_eff_MC_em_LP   = 1. - rrv_eff_MC_em.getVal()
      tmpq_eff_data_em_LP = 1. - rrv_eff_data_em.getVal()

      tmpq_eff_MC_em_LP_err   = TMath.Sqrt( rrv_eff_MC_em.getError()  **2)
      tmpq_eff_data_em_LP_err = TMath.Sqrt( rrv_eff_data_em.getError()**2)
      
      pureq_wtagger_sf_em_LP = tmpq_eff_data_em_LP / tmpq_eff_MC_em_LP
      pureq_wtagger_sf_em_LP_err = pureq_wtagger_sf_em_LP*TMath.Sqrt( (tmpq_eff_data_em_LP_err/tmpq_eff_data_em_LP)**2 + (tmpq_eff_MC_em_LP_err/tmpq_eff_MC_em_LP)**2 )
            
      # print ""
      # print ""
      # print ""
      # print "LP Eff of em data %s: %0.3f +/- %0.3f"%(label,eff_data_em_LP, eff_data_em_LP_err)
      # print "LP Eff of em MC %s: %0.3f +/- %0.3f"%(label,eff_MC_em_LP, eff_MC_em_LP_err)
      # print ""
      # print ""
      # print ""

      self.boostedW_fitter_em.file_out_ttbar_control.write("\nLP Eff of em data %s: %f +/- %f"%(label,eff_data_em_LP, eff_data_em_LP_err))
      self.boostedW_fitter_em.file_out_ttbar_control.write("\nLP Eff of em MC %s: %f +/- %f"  %(label,eff_MC_em_LP, eff_MC_em_LP_err))

      pure_wtagger_sf_em_LP     = eff_data_em_LP / eff_MC_em_LP
      pure_wtagger_sf_em_LP_err = pure_wtagger_sf_em_LP * ( (eff_data_em_LP_err/eff_data_em_LP)**2 + (eff_MC_em_LP_err/eff_MC_em_LP)**2 )**0.5

      print ""
      print ""
      print ""
      print "Pure W-tagger LP SF of em %s: %0.3f +/- %0.3f"%(label,pure_wtagger_sf_em_LP, pure_wtagger_sf_em_LP_err)
      print ""
      print ""
      print ""
      
      print "                                    LP                                     "
      print "---------------------------------------------------------------------------"
      print ""
      print "Parameter                                   Data                          Simulation                          Data/Simulation"
      print "LP W-tag eff+SF ( w/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(eff_data_em_LP,eff_data_em_LP_err,eff_MC_em_LP,eff_MC_em_LP_err,pure_wtagger_sf_em_LP,pure_wtagger_sf_em_LP_err)
      print "LP W-tag eff+SF (wo/ext fail)         %0.3f +/- %0.3f                  %0.3f +/- %0.3f                        %0.3f +/- %0.3f" %(tmpq_eff_data_em_LP,tmpq_eff_data_em_LP_err,tmpq_eff_MC_em_LP,tmpq_eff_MC_em_LP_err,pureq_wtagger_sf_em_LP,pureq_wtagger_sf_em_LP_err)
      print ""
      print "---------------------------------------------------------------------------"
      
      

      self.boostedW_fitter_em.file_out_ttbar_control.write("\nPure W-tagger LP SF of em %s: %f +/- %f"%(label,pure_wtagger_sf_em_LP, pure_wtagger_sf_em_LP_err))

class doFit_wj_and_wlvj:

    ## COnstructor: Input is channel (mu,ele,em), range in mj and a workspace
    def __init__(self, in_channel, in_mj_min=40, in_mj_max=130, input_workspace=None):
      
      RooAbsPdf.defaultIntegratorConfig().setEpsRel(1e-9)
      RooAbsPdf.defaultIntegratorConfig().setEpsAbs(1e-9)
      
      ### set the channel type
      self.channel = in_channel
      
      print "CHANNEL = %s" %in_channel

      ### shapes to be used in mj                                                                                                                                          
      self.mj_shape = ROOT.std.map(ROOT.std.string,ROOT.std.string)()
      self.mj_shape["TTbar"]  = "2Gaus_ErfExp"
      
      if (options.tau2tau1cutHP==0.60):
        self.mj_shape["STop"]               = "ExpGaus"
        self.mj_shape["STop_fail"]          = "ExpGaus"
        self.mj_shape["STop_extremefail"]   = "Exp"
        self.mj_shape["VV"]                 = "ExpGaus"
        self.mj_shape["VV_fail"]            = "ExpGaus"
        self.mj_shape["VV_extremefail"]     = "Exp"
        self.mj_shape["WJets0"]               = "ErfExp"
        self.mj_shape["WJets0_fail"]          = "Exp"
        self.mj_shape["WJets0_extremefail"]   = "Exp"

      elif (options.tau2tau1cutHP==0.45):
        if self.channel == "mu":
          self.mj_shape["STop"]               = "ErfExpGaus_sp" # Qun: ErfExpGaus_sp for MU and "ExpGaus" for EM/ELE
        else:
          self.mj_shape["STop"]               = "ExpGaus"       # Qun: "ExpGaus" for EM/ELE
        self.mj_shape["STop_fail"]          = "ExpGaus"       # Qun: "ExpGaus" for EM/ELE
        self.mj_shape["STop_extremefail"]   = "Exp"           # Qun: "Exp"
        self.mj_shape["VV"]                 = "Gaus"          # Qun: "ExpGaus"
        self.mj_shape["VV_fail"]            = "Exp"           # Qun: "ExpGaus"
        self.mj_shape["VV_extremefail"]     = "Exp"           # Qun: "Exp"
        self.mj_shape["WJets0"]               = "ExpGaus"        # Qun: "ExpGaus"
        self.mj_shape["WJets0_fail"]          = "Exp"           # Qun: "ExpGaus"
        self.mj_shape["WJets0_extremefail"]   = "Exp"           # Qun: "Exp"
      else:
        print "NO CHANNEL IS DEFINED!!! ABORT!!"
        sys.exit(0)
        
      self.mj_shape["bkg_data"]         = "ErfExp_ttbar"
      self.mj_shape["bkg_data_fail"]    = "ErfExp_ttbar_failtau2tau1cut"
      self.mj_shape["signal_data"]      = "2Gaus_ttbar" 
      self.mj_shape["signal_data_fail"] = "GausChebychev_ttbar_failtau2tau1cut"
      self.mj_shape["bkg_mc"]           = "ErfExp_ttbar"
      self.mj_shape["bkg_mc_fail"]      = "ErfExp_ttbar_failtau2tau1cut"
      self.mj_shape["signal_mc"]        = "2Gaus_ttbar" 
      self.mj_shape["signal_mc_fail"]   = "GausChebychev_ttbar_failtau2tau1cut"
      self.mj_shape["data_extremefail"]     = "Exp_ttbar_extremefailtau2tau1cut"
      self.mj_shape["data_bkg_extremefail"] = "Exp_bkg_extremefailtau2tau1cut"
      self.mj_shape["mc_extremefail"]       = "Exp_ttbar_extremefailtau2tau1cut"
      self.mj_shape["mc_bkg_extremefail"]   = "Exp_bkg_extremefailtau2tau1cut"
        

      self.Lumi=2100

      # Mj binning (for plotting)
      self.BinWidth_mj = 5.
      
      # Higgs-Combination-Tools generates binned sample, need bin width narrower than 5 as above.
      self.narrow_factor = 1.

      # Set range
      self.BinWidth_mj = self.BinWidth_mj/self.narrow_factor
      nbins_mj         = int( (in_mj_max - in_mj_min) / self.BinWidth_mj )
      in_mj_max        = in_mj_min+nbins_mj*self.BinWidth_mj

      # Declare RooRealVar
      rrv_mass_j = RooRealVar("rrv_mass_j","pruned jet mass",(in_mj_min+in_mj_max)/2.,in_mj_min,in_mj_max,"GeV")
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
      self.mj_sideband_hi_min = 135
      self.mj_sideband_hi_max = in_mj_max
 
      # Setting ranges...
      rrv_mass_j.setRange("sb_lo",self.mj_sideband_lo_min,self.mj_sideband_lo_max) # 30-65 GeV
      rrv_mass_j.setRange("signal_region",self.mj_signal_min,self.mj_signal_max)   # 65-105 GeV
      rrv_mass_j.setRange("sb_hi",self.mj_sideband_hi_min,self.mj_sideband_hi_max) # 105-135 GeV
      rrv_mass_j.setRange("controlsample_fitting_range",40,130) # ---> what is this????

      # Tree directory and input files
      self.file_Directory         = "/shome/thaarres/EXOVVAnalysisRunII/AnalysisOutput/Wtag/WWTree_%s/"%(self.channel)
      self.file_data              = ("ExoDiBosonAnalysis.WWTree_data.root")
      self.file_pseudodata        = ("ExoDiBosonAnalysis.WWTree_pseudodata.root")
      self.file_WJets0_mc         = ("ExoDiBosonAnalysis.WWTree_WJets.root")
      self.file_VV_mc             = ("ExoDiBosonAnalysis.WWTree_VV.root")            
      self.file_TTbar_mc          = ("ExoDiBosonAnalysis.WWTree_TTbar.root")        
      self.file_STop_mc           = ("ExoDiBosonAnalysis.WWTree_STop.root")          
      
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
      self.AK8_pt_min = 200
      self.AK8_pt_max = 5000  
      if self.channel == "el":
        self.pfMET_cut = 80
        self.lpt_cut = 120
        
      
      # Out .txt file
      self.file_ttbar_control_txt = "WtaggingSF_%s.txt"%(self.channel)
      self.file_out_ttbar_control = open(self.file_ttbar_control_txt,"w")
                                                                                                                                                             
      setTDRStyle()

    def fit_TTbar_controlsample(self):
  
      print "Fit ttbar control sample"
      rrv_mass_j = self.workspace4fit_.var("rrv_mass_j")
      
      ### Build single-t fit pass and fail distributions
      print "##################################################"
      print "############### Single Top DataSet ###############"
      print "##################################################"
    
      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_STop_mc,"_STop")

      fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop",self.mj_shape["STop"],self.channel,self.wtagger_label,1) #Start value and range of parameters defined in PDFs/MakePDF.cxx
      fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_failtau2tau1cut",self.mj_shape["STop_fail"],self.channel,self.wtagger_label,1)
      fit_mj_single_MC(self.workspace4fit_,self.file_STop_mc,"_STop_extremefailtau2tau1cut",self.mj_shape["STop_extremefail"],self.channel,self.wtagger_label,1)
 

      ### Build WJet fit pass and fail distributions
      print "###########################################"
      print "############### WJets Pythia ##############"
      print "###########################################"

      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_WJets0_mc,"_WJets0")

      fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0",self.mj_shape["WJets0"],self.channel,self.wtagger_label,1)
      fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_failtau2tau1cut",self.mj_shape["WJets0_fail"],self.channel,self.wtagger_label,1)
      fit_mj_single_MC(self.workspace4fit_,self.file_WJets0_mc,"_WJets0_extremefailtau2tau1cut",self.mj_shape["WJets0_extremefail"],self.channel,self.wtagger_label,1)



      ### Build VV fit pass and fail distributions
      print "#########################################"
      print "############### VV Pythia ###############"
      print "#########################################"

      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_VV_mc,"_VV")

      fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV",self.mj_shape["VV"],self.channel,self.wtagger_label,1)
      fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_failtau2tau1cut",self.mj_shape["VV_fail"],self.channel,self.wtagger_label,1)
      fit_mj_single_MC(self.workspace4fit_,self.file_VV_mc,"_VV_extremefailtau2tau1cut",self.mj_shape["VV_extremefail"],self.channel,self.wtagger_label,1)


      print "#########################################"
      print "############# TTbar Powheg ##############"
      print "#########################################"
      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_TTbar_mc,"_TTbar")

      print "################################################"
      print "############## Pseudo Data Powheg ##############"
      print "################################################"
      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_pseudodata,"_TotalMC")

      print "#################################"
      print "############# Data ##############"
      print "#################################"
      self.get_mj_and_mlvj_dataset_TTbar_controlsample(self.file_data,"_data")
 

      self.fit_mj_TTbar_controlsample(self.file_data);

      self.constrainslist_data = ROOT.std.vector(ROOT.std.string)()
      self.constrainslist_mc   = ROOT.std.vector(ROOT.std.string)()
  
      ScaleFactorTTbarControlSampleFit(self.workspace4fit_,self.mj_shape,self.color_palet,self.constrainslist_data,self.constrainslist_mc,"",self.channel,self.wtagger_label,self.AK8_pt_min,self.AK8_pt_max)


    def fit_mj_TTbar_controlsample(self,in_file_name):

        ##### Print number of events passing cut and before cut + the efficiency for the W-tagging -> dataset yields in the signal region
        
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

        number_dataset_signal_region_data_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_data_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_data_"+self.channel+"_mj").getVal()
        print ""
        print "Nr. data events in signal_region: %s +/- sqrt(%s)"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj)
        print ""

        self.file_out_ttbar_control.write("%s channel SF: \n"%(self.channel));
        self.file_out_ttbar_control.write("Nr. events in signal_region: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_data_mj, number_dataset_signal_region_error2_data_mj))

        number_dataset_signal_region_TotalMC_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_TotalMC_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_error2_TotalMC_"+self.channel+"_mj").getVal()
        print ""
        print "Nr. MC events in signal_region: %s +/- sqrt(%s) "%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj)
        print ""

        self.file_out_ttbar_control.write("Nr.  TotalMC in signal_region: %s +/- sqrt(%s) \n"%(number_dataset_signal_region_TotalMC_mj, number_dataset_signal_region_error2_TotalMC_mj))


        number_dataset_signal_region_before_cut_data_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_data_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_before_cut_error2_data_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_data_"+self.channel+"_mj").getVal()
        print ""
        print "Nr dataevents in signalregion before cut on tau21: %s +/- sqrt(%s)"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj)
        print ""
        self.file_out_ttbar_control.write("event number of data in signal_region before_cut: %s +/- sqrt(%s)\n"%(number_dataset_signal_region_before_cut_data_mj, number_dataset_signal_region_before_cut_error2_data_mj))

        number_dataset_signal_region_before_cut_TotalMC_mj        = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_TotalMC_"+self.channel+"_mj").getVal()
        number_dataset_signal_region_before_cut_error2_TotalMC_mj = self.workspace4fit_.var("rrv_number_dataset_signal_region_before_cut_error2_TotalMC_"+self.channel+"_mj").getVal()
        print ""
        print "Nr MC events in signal_region before cut on tau21: %s +/- sqrt(%s) "%(number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj)
        print ""
        self.file_out_ttbar_control.write("event number of TotalMC in signal_region before_cut: %s +/- sqrt(%s) \n"%(number_dataset_signal_region_before_cut_TotalMC_mj, number_dataset_signal_region_before_cut_error2_TotalMC_mj))
                                                             
        # wtagger_eff reweight: only reweight the efficiency difference between MC and data
        wtagger_eff_MC   = number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj
        wtagger_eff_data = number_dataset_signal_region_data_mj/number_dataset_signal_region_before_cut_data_mj

        wtagger_eff_reweight     = wtagger_eff_data/wtagger_eff_MC
        wtagger_eff_reweight_err = wtagger_eff_reweight*TMath.Sqrt(number_dataset_signal_region_error2_data_mj/number_dataset_signal_region_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_error2_TotalMC_mj/number_dataset_signal_region_TotalMC_mj/number_dataset_signal_region_TotalMC_mj +number_dataset_signal_region_before_cut_error2_data_mj/number_dataset_signal_region_before_cut_data_mj/number_dataset_signal_region_data_mj + number_dataset_signal_region_before_cut_error2_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj/number_dataset_signal_region_before_cut_TotalMC_mj)
        
        print "W-tagging efficiency for the %s channel:"%(self.channel)
        print "W-tagging eff. MC       = %s "%(wtagger_eff_MC)
        print "W-tagging eff. data     = %s "%(wtagger_eff_data)
        print "EFF. DATA/EFF. MC       = %s +/- %s"%(wtagger_eff_reweight, wtagger_eff_reweight_err)

        self.file_out_ttbar_control.write("wtagger_eff_MC       = %s       \n"%(wtagger_eff_MC ))
        self.file_out_ttbar_control.write("wtagger_eff_data     = %s       \n"%(wtagger_eff_data ))
        self.file_out_ttbar_control.write("wtagger_eff_reweight = %s +/- %s\n"%(wtagger_eff_reweight, wtagger_eff_reweight_err))
            
    # Loop over trees
    def get_mj_and_mlvj_dataset_TTbar_controlsample(self,in_file_name, label, jet_mass="Mjpruned"): 
    
      fileIn_name = TString(self.file_Directory+in_file_name)
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
      for i in range(treeIn.GetEntries()):
          if i % 5000 == 0: print "iEntry: ",i
          treeIn.GetEntry(i)
              
          discriminantCut = 0
          wtagger = getattr(treeIn,"tau21")

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
            # tmp_event_weight4fit = getattr(treeIn,"genweight")*getattr(treeIn,"puweight")*getattr(treeIn,"hltweight") #Missing b-tag weight and Vtg weight!
            # tmp_event_weight4fit = tmp_event_weight4fit*getattr(treeIn,"lumiweight")*self.Lumi/tmp_scale_to_lumi      #this is essentially the same as tmp_event_weight/tmp_scale_to_lumi?
            tmp_event_weight4fit = tmp_event_weight*self.Lumi/tmp_scale_to_lumi
          else:
            tmp_event_weight = 1.
            tmp_event_weight4fit = 1.

          #  HP category
          if discriminantCut == 2  and (rrv_mass_j.getMin() < tmp_jet_mass < rrv_mass_j.getMax()):   
          
             rrv_mass_j.setVal(tmp_jet_mass)
             
             rdataset_mj    .add(RooArgSet(rrv_mass_j), tmp_event_weight)
             rdataset4fit_mj.add(RooArgSet(rrv_mass_j), tmp_event_weight4fit)
  
             if tmp_jet_mass >= self.mj_sideband_lo_min and tmp_jet_mass < self.mj_sideband_lo_max:
                 hnum_4region.Fill(-1,tmp_event_weight )
             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max:
                 hnum_2region.Fill(1,tmp_event_weight)
             if tmp_jet_mass >= self.mj_signal_min and tmp_jet_mass < self.mj_signal_max :
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
      


      
      
      rrv_scale_to_lumi                        = RooRealVar("rrv_scale_to_lumi"+label+"_"                       +self.channel,"rrv_scale_to_lumi"+label+"_"                       +self.channel,rdataset_mj.sumEntries()/rdataset4fit_mj.sumEntries()) 
      rrv_scale_to_lumi_failtau2tau1cut        = RooRealVar("rrv_scale_to_lumi"+label+"_failtau2tau1cut_"       +self.channel,"rrv_scale_to_lumi"+label+"_failtau2tau1cut_"       +self.channel,rdataset_failtau2tau1cut_mj.sumEntries()/rdataset4fit_failtau2tau1cut_mj.sumEntries())
      rrv_scale_to_lumi_extremefailtau2tau1cut = RooRealVar("rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,"rrv_scale_to_lumi"+label+"_extremefailtau2tau1cut_"+self.channel,rdataset_extremefailtau2tau1cut_mj.sumEntries()/rdataset4fit_extremefailtau2tau1cut_mj.sumEntries()) 
        
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
      
      

      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi)
      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_failtau2tau1cut)
      getattr(self.workspace4fit_,"import")(rrv_scale_to_lumi_extremefailtau2tau1cut)
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

      # rdataset_mj.Print()
      # rdataset4fit_mj.Print()
      # rdataset_failtau2tau1cut_mj.Print()
      # rdataset4fit_failtau2tau1cut_mj.Print()
      # rdataset_extremefailtau2tau1cut_mj.Print()
      # rdataset4fit_extremefailtau2tau1cut_mj.Print()
      # rrv_number_dataset_sb_lo_mj.Print()
      # rrv_number_dataset_signal_region_mj.Print()
      # rrv_number_dataset_signal_region_error2_mj.Print()
      # rrv_number_dataset_signal_region_before_cut_mj.Print()
      # rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
      # rrv_number_dataset_sb_hi_mj.Print()
      #
      # rdataset_mj.Print()
      # rdataset_beforetau2tau1cut_mj.Print()
      # rdataset_failtau2tau1cut_mj.Print()
      # rdataset_extremefailtau2tau1cut_mj.Print()
      # rrv_number_dataset_signal_region_mj.Print()
      # rrv_number_dataset_signal_region_error2_mj.Print()
      # rrv_number_dataset_signal_region_before_cut_mj.Print()
      # rrv_number_dataset_signal_region_before_cut_error2_mj.Print()
      # combData_p_f.Print("v")

    def make_Pdf(self, label, in_model_name, mass_spectrum="_mj", ConstraintsList=[]):

      if TString(mass_spectrum).Contains("_mj"): 
        rrv_x = self.workspace4fit_.var("rrv_mass_j")
        
      if in_model_name == "Gaus":
        rrv_mean_gaus  = RooRealVar("rrv_mean_gaus"+label+"_"+self.channel,"rrv_mean_gaus"+label+"_"+self.channel,84,78,92);
        rrv_sigma_gaus = RooRealVar("rrv_sigma_gaus"+label+"_"+self.channel,"rrv_sigma_gaus"+label+"_"+self.channel,7,1,15);
        model_pdf = RooGaussian("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_mean_gaus,rrv_sigma_gaus);

      if in_model_name == "ErfExp":
        
        rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.026,-0.1, -1.e-4)
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,60.,5.,120)
        rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,1.,100.)
        model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)



      if in_model_name == "Exp":
        rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+self.channel,"rrv_c_Exp"+label+"_"+self.channel, -0.030, -0.2, 0.05)
        model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_Exp)
      if in_model_name == "ErfExpGaus_sp":
        rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,-0.04,-0.2,0.);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel,30.,10,300.);
        rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,80,40,100);
        erfExp = ROOT.RooErfExpPdf("erfExp"+label+"_"+self.channel,"erfExp"+label+"_"+self.channel,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp);

        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,0.,40);
        gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.5,0.,1.);
        model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfExp,gaus),RooArgList(rrv_high))


      if in_model_name == "ExpGaus":
        
        rrv_c_Exp = RooRealVar("rrv_c_ExpGaus"+label+"_"+self.channel,"rrv_c_ExpGaus"+label+"_"+self.channel,-0.05,-0.5,0.5)
        exp = ROOT.RooExponential("ExpGaus"+label+"_"+self.channel,"ExpGaus"+label+"_"+self.channel,rrv_x,rrv_c_Exp)
        rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,84,60,120)
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,7,4,40)
        gaus = RooGaussian("gaus"+label+"_"+self.channel,"gaus"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
        rrv_high = RooRealVar("rrv_high"+label+"_"+self.channel,"rrv_high"+label+"_"+self.channel,0.5,0.,1.);
        model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))
      if in_model_name == "2Gaus_ttbar":
        mean1_tmp = 7.6717e+01; mean1_tmp_err = 1.84612e+00
        deltamean_tmp = 1.6880e+00; deltamean_tmp_err = 7.28832e-01
        sigma1_tmp = 7.3633e+00; sigma1_tmp_err = 9.88917e-02
        scalesigma_tmp = 3.794e+00; scalesigma_tmp_err = 9.91214e-01
        frac_tmp = 6.24209e-01; frac_tmp_err = 8.42986e-02

        ## constrain the peak and the width to be the same in case of simultaneous fit
        if self.channel=="el":
            if self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label)):
                rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label));
                rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label));
            else:
                rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-25, mean1_tmp+25);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-20,sigma1_tmp+20);
        if self.channel=="mu":
            if self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label)):
                rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label));
                rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label));
            else:
                rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-25, mean1_tmp+25);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-20,sigma1_tmp+20 );

        gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);
        rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#,deltamean_tmp-deltamean_tmp_err*5, deltamean_tmp+deltamean_tmp_err*5);
        rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#,scalesigma_tmp-scalesigma_tmp_err*4, scalesigma_tmp+scalesigma_tmp_err*4);
        rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);
        rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp ,frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);
        model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)


      if in_model_name == "2Gaus_ttbar_a":
        
        mean1_tmp = 8.0717e+01; mean1_tmp_err = 1.84612e+00;
        deltamean_tmp = 1.6880e+00; deltamean_tmp_err = 7.28832e-01;
        sigma1_tmp = 7.3633e+00; sigma1_tmp_err = 9.88917e-02;
        scalesigma_tmp = 3.794e+00; scalesigma_tmp_err = 9.91214e-01;
        frac_tmp = 6.24209e-01; frac_tmp_err = 8.42986e-02;

        ## constrain the peak and the width among electron and muon channel to be the same in case of simultaneous fit
        if self.channel=="el":
            if self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label)):
                rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_mu"%(label));
                rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_mu"%(label));
            else:
                rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-10, mean1_tmp+10);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-15,sigma1_tmp+15);
        if self.channel=="mu":
            if self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label)) and self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label)):
                rrv_mean1_gaus = self.workspace4fit_.var("rrv_mean1_gaus%s_el"%(label));
                rrv_sigma1_gaus = self.workspace4fit_.var("rrv_sigma1_gaus%s_el"%(label));
            else:
                rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-10, mean1_tmp+10);
                rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-15,sigma1_tmp+15 );

        gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        p0_tmp = -4.73757e-01; p0_tmp_err = 2.42754e-02;
        p1_tmp = 9.57687e-02; p1_tmp_err = 2.93041e-02;
        frac_tmp = 2.90629e-01; frac_tmp_err = 1.86392e-02;


        if TString(label).Contains("herwig") and not TString(label).Contains("data") and self.channel == "mu":
          p0_tmp = -6.6824e-01; p0_tmp_err = 5.04e-02
          p1_tmp = -5.3365e-02; p1_tmp_err = 6.74e-02
          frac_tmp = 2.0421e-01; frac_tmp_err = 2.48e-02
          
        rrv_p0_cheb = RooRealVar("rrv_p0_cheb"+label+"_"+self.channel,"rrv_p0_cheb"+label+"_"+self.channel,p0_tmp)
        rrv_p1_cheb = RooRealVar("rrv_p1_cheb"+label+"_"+self.channel,"rrv_p1_cheb"+label+"_"+self.channel,p1_tmp)
        cheb = RooChebychev("cheb"+label+"_"+self.channel,"cheb"+label+"_"+self.channel, rrv_x, RooArgList(rrv_p0_cheb, rrv_p1_cheb) )

        rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp)
        model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,cheb),RooArgList(rrv_frac),1)# qun

      if in_model_name == "GausChebychev_ttbar_failtau2tau1cut":
        
        p0_tmp = -4.73757e-01 ; p0_tmp_err = 2.42754e-02
        p1_tmp = 9.57687e-02  ; p1_tmp_err = 2.93041e-02
        frac_tmp = 2.90629e-01; frac_tmp_err = 1.86392e-02

        if TString(label).Contains("herwig") and not TString(label).Contains("data") and self.channel == "mu":
          p0_tmp = -6.6824e-01; p0_tmp_err = 5.04e-02
          p1_tmp = -5.3365e-02; p1_tmp_err = 6.74e-02
          frac_tmp = 2.0421e-01; frac_tmp_err = 2.48e-02;

        ## take the same gaussian used in the pass sample
        if TString(label).Contains("data") and not TString(label).Contains("herwig"):
              gaus = self.workspace4fit_.pdf("gaus1_ttbar_data_"+self.channel);
        if TString(label).Contains("TotalMC") and not TString(label).Contains("herwig"):
              gaus = self.workspace4fit_.pdf("gaus1_ttbar_TotalMC_"+self.channel);

        if TString(label).Contains("data") and TString(label).Contains("herwig"):
              gaus = self.workspace4fit_.pdf("gaus1_ttbar_data_herwig_"+self.channel);
        if TString(label).Contains("TotalMC") and TString(label).Contains("herwig"):
              gaus = self.workspace4fit_.pdf("gaus1_ttbar_TotalMC_herwig_"+self.channel);

        rrv_p0_cheb = RooRealVar("rrv_p0_cheb"+label+"_"+self.channel,"rrv_p0_cheb"+label+"_"+self.channel,p0_tmp)#,p0_tmp-p0_tmp_err*4,p0_tmp+p0_tmp_err*4);
        rrv_p1_cheb = RooRealVar("rrv_p1_cheb"+label+"_"+self.channel,"rrv_p1_cheb"+label+"_"+self.channel,p1_tmp)#,p1_tmp-p1_tmp_err*4,p1_tmp+p1_tmp_err*4);
        cheb = RooChebychev("cheb"+label+"_"+self.channel,"cheb"+label+"_"+self.channel, rrv_x, RooArgList(rrv_p0_cheb, rrv_p1_cheb) );

        rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel,frac_tmp)#,frac_tmp-frac_tmp_err*4,frac_tmp+frac_tmp_err*4);
        model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus,cheb),RooArgList(rrv_frac),1)
        
        model_pdf = RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

      if in_model_name == "ErfExp_ttbar":
        
        c0_tmp = -4.18408e-02;  c0_tmp_err = 3.06900e-03
        offset_tmp = 8.27928e+01; offset_tmp_err = 3.20278e+00
        width_tmp = 2.94111e+01; width_tmp_err = 8.81193e-01
        
        if TString(label).Contains("herwig") and not TString(label).Contains("data"):
          
          c0_tmp = -2.9357e-02 ; c0_tmp_err = 6.83e-03;
          offset_tmp = 7.9350e+01 ; offset_tmp_err = 9.35e+00;
          width_tmp = 3.3216e+01 ; width_tmp_err = 2.97e+00;

        rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 );
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp,width_tmp-10, width_tmp+10);

        model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        ## Constrint the parameter within one sigma on the background part of the non resonant component in the pass sample
        self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);
        self.addConstraint(rrv_offset_ErfExp,offset_tmp, offset_tmp_err,ConstraintsList);
        self.addConstraint(rrv_width_ErfExp,width_tmp, width_tmp_err,ConstraintsList);

      if in_model_name == "ErfExp_ttbar_failtau2tau1cut":
        
        c0_tmp = -5.44643e-02 ; c0_tmp_err = 1.43399e-03;
        offset_tmp = 6.594e+01; offset_tmp_err = 3.40396e+00;
        width_tmp = 3.90e+01  ; width_tmp_err = 1.92162e+00;
        
        if TString(label).Contains("herwig") and not TString(label).Contains("data"):
          c0_tmp = -1.0141e-01 ; c0_tmp_err = 1.46e-02;
          offset_tmp = 2.6730e+02 ; offset_tmp_err = 4.92e+01;
          width_tmp = 7.1505e+01 ; width_tmp_err = 4.69e+00;


        rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp,c0_tmp-4e-1, c0_tmp+4e-1);
        rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel,offset_tmp, offset_tmp-offset_tmp_err*20,offset_tmp+offset_tmp_err*20);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-width_tmp_err*20, width_tmp+width_tmp_err*20);
        model_pdf = ROOT.RooErfExpPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);
        ## constrain just the slope of the exponential
        self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);

      if in_model_name == "Exp_ttbar_extremefailtau2tau1cut":
        
        c0_tmp = -3.0278e-02 ; c0_tmp_err = 5.16e-03
        rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
        model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp);
        ## constraint within one sigma
        self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);

      if in_model_name == "Exp_bkg_extremefailtau2tau1cut":
        c0_tmp = -4.2105e-02 ; c0_tmp_err = 2.61e-03;
        rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-1, c0_tmp+4e-1 );
        model_pdf = ROOT.RooExponential("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp);
        ## constraint within one sigma
        self.addConstraint(rrv_c_ErfExp,c0_tmp, c0_tmp_err,ConstraintsList);



      if in_model_name == "2Gaus_ErfExp":
        
        mean1_tmp =8.3145e+01; mean1_tmp_err =6.94e-01;
        deltamean_tmp =6.6321e+00; deltamean_tmp_err =1.21e+00;
        sigma1_tmp =7.5097e+00; sigma1_tmp_err =2.01e-01;
        scalesigma_tmp=3.8707e+00; scalesigma_tmp_err=2.20e-01;
        frac_tmp =6.4728e-01; frac_tmp_err =2.03e-02;

        if self.wtagger_cut==0.60:
          mean1_tmp =7.6273e+01; mean1_tmp_err =6.94e-01;
          deltamean_tmp =6.9129e+00; deltamean_tmp_err =1.24e+00;
          sigma1_tmp =7.9599e+00; sigma1_tmp_err =9.84e-01;
          scalesigma_tmp=3.6819e+00; scalesigma_tmp_err=2.11e-01;
          frac_tmp =6.7125e-01; frac_tmp_err =2.09e-02;

        rrv_mean1_gaus = RooRealVar("rrv_mean1_gaus"+label+"_"+self.channel,"rrv_mean1_gaus"+label+"_"+self.channel,mean1_tmp, mean1_tmp-15, mean1_tmp+15);
        rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+self.channel,"rrv_sigma1_gaus"+label+"_"+self.channel,sigma1_tmp, sigma1_tmp-4,sigma1_tmp+4 );
        gaus1 = RooGaussian("gaus1"+label+"_"+self.channel,"gaus1"+label+"_"+self.channel, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus);

        rrv_deltamean_gaus = RooRealVar("rrv_deltamean_gaus"+label+"_"+self.channel,"rrv_deltamean_gaus"+label+"_"+self.channel,deltamean_tmp)#, deltamean_tmp, deltamean_tmp);
        rrv_mean2_gaus = RooFormulaVar("rrv_mean2_gaus"+label+"_"+self.channel,"@0+@1",RooArgList(rrv_mean1_gaus, rrv_deltamean_gaus));
        rrv_scalesigma_gaus = RooRealVar("rrv_scalesigma_gaus"+label+"_"+self.channel,"rrv_scalesigma_gaus"+label+"_"+self.channel,scalesigma_tmp)#, scalesigma_tmp, scalesigma_tmp);
        rrv_sigma2_gaus = RooFormulaVar("rrv_sigma2_gaus"+label+"_"+self.channel,"@0*@1", RooArgList(rrv_sigma1_gaus,rrv_scalesigma_gaus));
        gaus2 = RooGaussian("gaus2"+label+"_"+self.channel,"gaus2"+label+"_"+self.channel, rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus);

        rrv_frac_2gaus = RooRealVar("rrv_frac_2gaus"+label+"_"+self.channel,"rrv_frac_2gaus"+label+"_"+self.channel,frac_tmp);#, frac_tmp-frac_tmp_err*4, frac_tmp+frac_tmp_err*4);

        c0_tmp = -2.8628e-02 ; c0_tmp_err = 6.08e-03;
        offset_tmp= 7.6259e+01 ; offset_tmp_err = 9.17e+00;
        width_tmp = 3.4207e+01 ; width_tmp_err = 3.18e+00;

        if self.wtagger_cut==0.60:
          c0_tmp = -2.9893e-02 ; c0_tmp_err = 6.83e-03;
          offset_tmp= 7.9350e+01 ; offset_tmp_err = 9.35e+00;
          width_tmp = 3.3083e+01 ; width_tmp_err = 2.97e+00;

        rrv_c_ErfExp = RooRealVar("rrv_c_ErfExp"+label+"_"+self.channel,"rrv_c_ErfExp"+label+"_"+self.channel,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 );
        rrv_offset_ErfExp= RooRealVar("rrv_offset_ErfExp"+label+"_"+self.channel,"rrv_offset_ErfExp"+label+"_"+self.channel, offset_tmp)#, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4);
        rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+self.channel,"rrv_width_ErfExp"+label+"_"+self.channel, width_tmp, width_tmp-10, width_tmp+10);
        erfexp = ROOT.RooErfExpPdf("erfexp"+label+"_"+self.channel+mass_spectrum,"erfexp"+label+"_"+self.channel+mass_spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp);

        rrv_frac = RooRealVar("rrv_frac"+label+"_"+self.channel,"rrv_frac"+label+"_"+self.channel, 0.5,0.,1.);
        model_pdf =RooAddPdf("model_pdf"+label+"_"+self.channel+mass_spectrum,"model_pdf"+label+"_"+self.channel+mass_spectrum,RooArgList(erfexp, gaus1,gaus2),RooArgList(rrv_frac, rrv_frac_2gaus),1)
      
      getattr(self.workspace4fit_,"import")(model_pdf)
      return self.workspace4fit_.pdf("model_pdf"+label+"_"+self.channel+mass_spectrum)



### Start  main
if __name__ == '__main__':
  
  channel = options.channel
   
  if options.fitwtaggersim:
    print 'fitwtagger adding el+mu sample' #I am actually not doing a simoultaneous fit. So..... change this
    control_sample_simultaneous()

  elif not options.fitwtaggersim:
    print 'fitwtagger for %s sample'%(channel)
    control_sample(channel)