import os,sys
import ROOT as rt

def getPavetext():
  addInfo = rt.TPaveText(0.73010112,0.2566292,0.8202143,0.5523546,"NDC")
  addInfo.SetFillColor(0)
  addInfo.SetLineColor(0)
  addInfo.SetFillStyle(0)
  addInfo.SetBorderSize(0)
  addInfo.SetTextFont(42)
  addInfo.SetTextSize(0.040)
  addInfo.SetTextAlign(12)
  return addInfo
  
def getLegend(mplot,channel,x1=0.62,x2=0.7,y1=0.92,y2=0.9):
    theLeg = rt.TLegend(x1,x2,y1,y2,"","brNDC")
    for obj in range(0,int(mplot.numItems())):
		
		objName = mplot.nameOf(obj);
		
		if objName == "errorband":
			 objName = "Uncertainty"
		if not objName.find("Uncertainty")!=-1 or objName.find("invisible")!=-1 or objName.find("TLine")!=-1:
		   theObj = mplot.getObject(obj)
		   objTitle = objName
		   drawoption = mplot.getDrawOptions(objName)
		   label = drawoption
		   if drawoption == "P": drawoption = "pe"
		   if label: label = "pe";
		if objName.find("Uncertainty")!=-1 or objName.find("sigma")!=-1: theLeg.AddEntry(theObj, objName,label)
		elif objName.find("Graph")!=-1: theLeg.AddEntry(theObj, "Uncertainty",label)
		else:
			if     objName == "STop": theLeg.AddEntry(theObj, "Single Top",label);
			elif objName.find("_TTbar_realW_failtau2tau1cut")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (W-jet)","ple");
			elif objName.find("_TTbar_realW")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (W-jet)","ple");
			elif objName.find("_TTbar_fakeW_failtau2tau1cut")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (non-W)","ple");
			elif objName.find("_TTbar_fakeW")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (non-W)","ple");
			elif objName.find("TTbar_realW")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (merged)","f");
			elif objName.find("TTbar_fakeW")!=-1: theLeg.AddEntry(theObj, "t#bar{t} (unmerged)","f");
			elif objName.find("TTbar")!=-1: theLeg.AddEntry(theObj, "t#bar{t}",label);
			elif objName.find("VV")!=-1:    theLeg.AddEntry(theObj, "WW/WZ/ZZ",label);
			elif objName.find("data")!=-1:  theLeg.AddEntry(theObj, "data",label);
			elif objName.find("WJets")!=-1: theLeg.AddEntry(theObj, "W+jets",label);
	
def fit_mj_single_MC(workspace,fileName,label, model,channel, wtagger_label):

	rrv_mass_j  = workspace.var("rrv_mass_j")
	rdataset_mj = workspace.data("rdataset4fit"+label+"_"+channel+"_mj")

	constraint_list = []
	from python.makepdf import MakeExtendedModel
	MakeExtendedModel(workspace,label,model,"_mj",channel,wtagger_label,constraint_list)
	model_pdf = workspace.pdf("model"+label+"_"+channel+"_mj")
	
	rfresult = model_pdf.fitTo(rdataset_mj,rt.RooFit.Save(1),rt.RooFit.SumW2Error(rt.kTRUE),rt.RooFit.Extended(rt.kTRUE), rt.RooFit.Minimizer("Minuit2"),rt.RooFit.Verbose(rt.kFALSE))
	rfresult = model_pdf.fitTo(rdataset_mj,rt.RooFit.Save(1),rt.RooFit.SumW2Error(rt.kTRUE),rt.RooFit.Extended(rt.kTRUE), rt.RooFit.Minimizer("Minuit2"),rt.RooFit.Verbose(rt.kFALSE))
	rfresult = model_pdf.fitTo(rdataset_mj,rt.RooFit.Save(1),rt.RooFit.SumW2Error(rt.kTRUE),rt.RooFit.Extended(rt.kTRUE), rt.RooFit.Minimizer("Minuit2"))

	getattr(workspace,'import')(rfresult)


	mplot = rrv_mass_j.frame(rt.RooFit.Title((label+" fitted by "+model)), rt.RooFit.Bins(int(rrv_mass_j.getBins())))
	mplot.GetYaxis().SetRangeUser(0,mplot.GetMaximum()*1.2)
	rdataset_mj.plotOn(mplot,rt.RooFit.MarkerSize(1.5),rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.XErrorSize(0),rt.RooFit.Invisible())
  
	model_pdf.plotOn(mplot,rt.RooFit.Name( "Gaussian comp. 1" ),rt.RooFit.Components("gaus1*"),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.LineColor(rt.kRed+3))
	model_pdf.plotOn(mplot,rt.RooFit.Name( "Gaussian comp. 2" ),rt.RooFit.Components("gaus2*"),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.LineColor(rt.kRed+4))
	model_pdf.plotOn(mplot,rt.RooFit.Name( "ErfExp comp." ),rt.RooFit.Components("erfExp*"),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.LineColor(rt.kRed+3))
	model_pdf.plotOn(mplot,rt.RooFit.Name( "Exp. comp." ),rt.RooFit.Components("exp*"),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.LineColor(rt.kRed+2))
	model_pdf.plotOn(mplot,rt.RooFit.Name( "Chebychev comp." ),rt.RooFit.Components("cheb*"),rt.RooFit.LineStyle(rt.kDashed),rt.RooFit.LineColor(rt.kRed+3))
	model_pdf.plotOn(mplot,rt.RooFit.Name((model) )) 

	rdataset_mj.plotOn(mplot,rt.RooFit.MarkerSize(1.),rt.RooFit.Name( (label) ), rt.RooFit.DataError(rt.RooAbsData.SumW2), rt.RooFit.XErrorSize(1))
	leg1 = getLegend(mplot,channel)
	mplot.addObject(leg1)

	#TODO RooPlot* mplot_pull = get_ratio(rrv_mass_j,rdataset_mj,model_pdf,rfresult,0,1)
	mplot.GetYaxis().SetRangeUser(0,mplot.GetMaximum()*1.2)
	mplot.GetYaxis().SetTitle("MC events / 5 GeV")
  
	
	datahist = rdataset_mj.binnedClone(rdataset_mj.GetName()+"_binnedClone", rdataset_mj.GetName()+"_binnedClone")
	Nbin = int(rrv_mass_j.getBins())
	rresult_param = rfresult.floatParsFinal()
	nparameters =  rresult_param.getSize()
	ChiSquare = model_pdf.createChi2(datahist,rt.RooFit.Extended(rt.kTRUE),rt.RooFit.DataError(rt.RooAbsData.Poisson))
	chi_over_ndf= ChiSquare.getVal()/(Nbin - nparameters)

	
	Name  = "#chi^{2}/ndf = %0.2f " %float(chi_over_ndf)
	cs = getPavetext()
	cs.AppendPad("same")
	mplot.addObject(cs)
   
	try: os.stat("plots/MCfits/") 
	except: os.makedirs("plots/MCfits/")

	# TODO # draw_canvas_with_pull(mplot,mplot_pull,RooArgList(*parameters_list),"plots/MCfits/",label+fileName,model,channel,0,0,GetLumi())
  
	workspace.var("rrv_number"+label+"_"+channel+"_mj").setVal(workspace.var("rrv_number"+label+"_"+channel+"_mj").getVal()*workspace.var("rrv_scale_to_lumi"+label+"_"+channel).getVal())
	workspace.var("rrv_number"+label+"_"+channel+"_mj").setError(workspace.var("rrv_number"+label+"_"+channel+"_mj").getError()*workspace.var("rrv_scale_to_lumi"+label+"_"+channel).getVal())
                                                                                                          
def ScaleFactorTTbarControlSampleFit(workspace, mj_shape, color_palet, constraintlist_data, constraintlist_MC, channel,wtagger):
  
  data      = workspace.data(("rdataset_data_beforetau2tau1cut_"  +channel+"_mj")).sumEntries()
  tt        = workspace.data(("rdataset_TTbar_beforetau2tau1cut_" +channel+"_mj")).sumEntries()
  # #// tt       = workspace.data(("rdataset_TTbar_realW_"+channel+"_mj")).sumEntries()
 #  #// tt            +=  workspace.data(("rdataset_TTbar_fakeW_"+channel+"_mj")).sumEntries()
  minorBKG  = workspace.data(("rdataset_WJets0_beforetau2tau1cut_"+channel+"_mj")).sumEntries()
  minorBKG += workspace.data(("rdataset_VV_beforetau2tau1cut_"    +channel+"_mj")).sumEntries()
  minorBKG += workspace.data(("rdataset_STop_beforetau2tau1cut_"  +channel+"_mj")).sumEntries()

  ttScalefactor = (data-minorBKG)/tt
  
  scalefactor = RooRealVar("tt_scalefactor","tt_scalefactor",ttScalefactor)
  getattr(workspace,'import')(scalefactor)

  label = ""
  rrv_mass_j = workspace.var("rrv_mass_j")
  rdataset_data_mj  =  workspace.data(("rdataset_data"+ label+"_"+channel+"_mj"))
  rdataset_TTbar_mj =  workspace.data(("rdataset_TTbar"+ label+"_"+channel+"_mj"))
  rdataset_STop_mj  =  workspace.data(("rdataset_STop"+  label+"_"+channel+"_mj"))
  rdataset_VV_mj    =  workspace.data(("rdataset_VV"+    label+"_"+channel+"_mj"))
  rdataset_WJets_mj =  workspace.data(("rdataset_WJets0"+label+"_"+channel+"_mj"))
 
  rdataset_TTbar_mj_merged   =  workspace.data(("rdataset_TTbar_realW"+ label+"_"+channel+"_mj"))
  rdataset_TTbar_mj_unmerged =  workspace.data(("rdataset_TTbar_fakeW"+ label+"_"+channel+"_mj"))
  
  from python.makepdf import change_dataset_to_histpdf
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj)

  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_merged)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_unmerged)
  
  model_histpdf_TTbar = workspace.pdf(((rdataset_TTbar_mj.GetName())+"_histpdf"))
  model_histpdf_STop  = workspace.pdf(((rdataset_STop_mj.GetName())+"_histpdf"))
  model_histpdf_VV    = workspace.pdf(((rdataset_VV_mj.GetName())+"_histpdf"))
  model_histpdf_WJets = workspace.pdf(((rdataset_WJets_mj.GetName())+"_histpdf"))

  model_histpdf_TTbar_merged   = workspace.pdf(((rdataset_TTbar_mj_merged.GetName())+"_histpdf"))
  model_histpdf_TTbar_unmerged = workspace.pdf(((rdataset_TTbar_mj_unmerged.GetName())+"_histpdf"))
   

  number_TTbar = RooRealVar(("rrv_number_TTbar"+label+"_"+channel),("rrv_number_TTbar"+label+"_"+channel),rdataset_TTbar_mj.sumEntries()*ttScalefactor)
  number_STop  = RooRealVar(("rrv_number_STop"+label+"_"+channel),("rrv_number_STop"+label+"_"+channel),rdataset_STop_mj.sumEntries())
  number_VV    = RooRealVar(("rrv_number_VV"+label+"_"+channel),("rrv_number_VV"+label+"_"+channel),rdataset_VV_mj.sumEntries())
  number_WJets = RooRealVar(("rrv_number_WJets"+label+"_"+channel),("rrv_number_WJets"+label+"_"+channel),rdataset_WJets_mj.sumEntries())

  number_TTbar_merged   = RooRealVar(("rrv_number_TTbar_realW"+label+"_"+channel),("rrv_number_TTbar_realW"+label+"_"+channel),rdataset_TTbar_mj_merged.sumEntries()*ttScalefactor)
  number_TTbar_unmerged = RooRealVar(("rrv_number_TTbar_fakeW"+label+"_"+channel),("rrv_number_TTbar_fakeW"+label+"_"+channel),rdataset_TTbar_mj_unmerged.sumEntries()*ttScalefactor)
  model_TTbar_STop_VV_WJets = RooAddPdf(("model_TTbar_STop_VV_WJets"+label+"_"+channel),("model_TTbar_STop_VV_WJets"+label+"_"+channel),RooArgList(model_histpdf_TTbar_merged,model_histpdf_TTbar_unmerged,model_histpdf_STop,model_histpdf_VV,model_histpdf_WJets),RooArgList(number_TTbar_merged,number_TTbar_unmerged,number_STop,number_VV,number_WJets))
  
  getattr(workspace,'import')(model_TTbar_STop_VV_WJets)

  rdataset_data_mj_fail  =  workspace.data(("rdataset_data"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_TTbar_mj_fail =  workspace.data(("rdataset_TTbar"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_STop_mj_fail  =  workspace.data(("rdataset_STop"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_VV_mj_fail    =  workspace.data(("rdataset_VV"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_WJets_mj_fail =  workspace.data(("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_TTbar_mj_fail_merged   =  workspace.data(("rdataset_TTbar_realW"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))
  rdataset_TTbar_mj_fail_unmerged =  workspace.data(("rdataset_TTbar_fakeW"+label+"_"+"failtau2tau1cut_"+channel+"_mj"))

  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_fail)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_STop_mj_fail)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_VV_mj_fail)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_WJets_mj_fail)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_fail_merged)
  change_dataset_to_histpdf(workspace,rrv_mass_j,rdataset_TTbar_mj_fail_unmerged)

  model_histpdf_TTbar_fail = workspace.pdf(((rdataset_TTbar_mj_fail.GetName())+"_histpdf"))
  model_histpdf_STop_fail  = workspace.pdf(((rdataset_STop_mj_fail.GetName())+"_histpdf"))
  model_histpdf_VV_fail    = workspace.pdf(((rdataset_VV_mj_fail.GetName())+"_histpdf"))
  model_histpdf_WJets_fail = workspace.pdf(((rdataset_WJets_mj_fail.GetName())+"_histpdf"))

  model_histpdf_TTbar_fail_merged   = workspace.pdf(((rdataset_TTbar_mj_fail_merged.GetName())+"_histpdf"))
  model_histpdf_TTbar_fail_unmerged = workspace.pdf(((rdataset_TTbar_mj_fail_unmerged.GetName())+"_histpdf"))

  number_TTbar_fail  = RooRealVar(("rrv_number_TTbar_fail"+label+"_"+channel),("rrv_number_TTbar_fail"+label+"_"+channel),rdataset_TTbar_mj_fail.sumEntries()*ttScalefactor)
  number_STop_fail   = RooRealVar(("rrv_number_STop_fail"+label +"_"+channel),("rrv_number_STop_fail"+label+"_"+channel),rdataset_STop_mj_fail.sumEntries())
  number_VV_fail     = RooRealVar(("rrv_number_VV_fail"+label   +"_"+channel),("rrv_number_VV_fail"+label+"_"+channel),rdataset_VV_mj_fail.sumEntries())
  number_WJets_fail  = RooRealVar(("rrv_number_WJets_fail"+label+"_"+channel),("rrv_number_WJets_fail"+label+"_"+channel),rdataset_WJets_mj_fail.sumEntries())

  number_TTbar_fail_merged    = RooRealVar(("rrv_number_TTbar_fail_realW"+label+"_"+channel),("rrv_number_TTbar_fail_realW"+label+"_"+channel),rdataset_TTbar_mj_fail_merged.sumEntries()*ttScalefactor)
  number_TTbar_fail_unmerged  = RooRealVar(("rrv_number_TTbar_fail_fakeW"+label+"_"+channel),("rrv_number_TTbar_fail_fakeW"+label+"_"+channel),rdataset_TTbar_mj_fail_unmerged.sumEntries()*ttScalefactor)
  
  model_TTbar_STop_VV_WJets_fail = RooAddPdf(("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel),("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel),RooArgList(model_histpdf_TTbar_fail_merged,model_histpdf_TTbar_fail_unmerged,model_histpdf_STop_fail,model_histpdf_VV_fail,model_histpdf_WJets_fail), RooArgList(number_TTbar_fail_merged,number_TTbar_fail_unmerged,number_STop_fail,number_VV_fail,number_WJets_fail))
  getattr(workspace,'import')(model_TTbar_STop_VV_WJets_fail)
  
  
  scale_number_TTbar_STop_VV_WJets      = (rdataset_TTbar_mj_merged.sumEntries()*ttScalefactor+rdataset_TTbar_mj_unmerged.sumEntries()*ttScalefactor+rdataset_STop_mj.sumEntries()+rdataset_VV_mj.sumEntries()+rdataset_WJets_mj.sumEntries())/(rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries())
  scale_number_TTbar_STop_VV_WJets_fail = (rdataset_TTbar_mj_fail_merged.sumEntries()*ttScalefactor+rdataset_TTbar_mj_fail_unmerged.sumEntries()*ttScalefactor+rdataset_STop_mj_fail.sumEntries()+rdataset_VV_mj_fail.sumEntries()+rdataset_WJets_mj_fail.sumEntries())/( rdataset_data_mj.sumEntries()+rdataset_data_mj_fail.sumEntries())

  rrv_scale_number_TTbar_STop_VV_WJets      = RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets"+label),("rrv_scale_number_TTbar_STop_VV_WJets"+label),scale_number_TTbar_STop_VV_WJets)
  rrv_scale_number_TTbar_STop_VV_WJets_fail = RooRealVar(("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label),("rrv_scale_number_TTbar_STop_VV_WJets_fail"+label),scale_number_TTbar_STop_VV_WJets_fail)
  getattr(workspace,'import')(rrv_scale_number_TTbar_STop_VV_WJets)
  getattr(workspace,'import')(rrv_scale_number_TTbar_STop_VV_WJets_fail)


  #/// Fix all shape parameters and normalization                                                                                                         
  model_STop       = get_General_mj_Model(workspace,"_STop"+label                      ,mj_shape["STop"]       ,channel)
  model_VV         = get_General_mj_Model(workspace,"_VV"+label                        ,mj_shape["VV"]         ,channel)
  model_WJets      = get_General_mj_Model(workspace,"_WJets0"+label                    ,mj_shape["WJets0"]     ,channel)
  model_STop_fail  = get_General_mj_Model(workspace,"_STop"+label+"_failtau2tau1cut"   ,mj_shape["STop_fail"]  ,channel)
  model_VV_fail    = get_General_mj_Model(workspace,"_VV"+label+"_failtau2tau1cut"     ,mj_shape["VV_fail"]    ,channel)
  model_WJets_fail = get_General_mj_Model(workspace,"_WJets0"+label+"_failtau2tau1cut" ,mj_shape["WJets0_fail"],channel)


  model_bkg_data         = MakeExtendedModel(workspace,"_bkg_data"+label,mj_shape["bkg_data"],"_mj",channel,wtagger,constraintlist_data)
  model_bkg_data_fail    = MakeExtendedModel(workspace,"_bkg_data"+label+"_failtau2tau1cut",mj_shape["bkg_data_fail"],"_mj",channel,wtagger,constraintlist_data)                                                      
  model_ttbar_data       = MakeModelTTbarControlSample(workspace,"_ttbar_data"+label,mj_shape["signal_data"],"_mj",channel,wtagger,label,constraintlist_data)    
  model_ttbar_data_fail  = MakeModelTTbarControlSample(workspace,"_ttbar_data"+label+"_failtau2tau1cut",mj_shape["signal_data_fail"],"_mj",channel,wtagger,label,constraintlist_data)


  model_data_fail = RooAddPdf(("model_data"+label+"_"+"failtau2tau1cut"+"_"+channel),("model_data+"+label+"_"+"failtau2tau1cut"+"_"+channel),RooArgList(model_ttbar_data_fail,model_bkg_data_fail,model_STop_fail,model_VV_fail,model_WJets_fail))
  model_data      = RooAddPdf(("model_data"+label+"_"+channel),("model_data"+label+"_"+channel), RooArgList(model_ttbar_data,model_bkg_data,model_STop,model_VV,model_WJets))

  getattr(workspace,'import')(model_data)
  getattr(workspace,'import')(model_data_fail)

  
  model_bkg_TotalMC        = MakeExtendedModel(workspace,"_bkg_TotalMC"+label,mj_shape["bkg_mc"],"_mj",channel,wtagger,constraintlist_MC)
  model_bkg_TotalMC_fail   = MakeExtendedModel(workspace,"_bkg_TotalMC"+label+"_failtau2tau1cut",mj_shape["bkg_mc_fail"],"_mj",channel,wtagger,constraintlist_MC)
  model_ttbar_TotalMC      = MakeModelTTbarControlSample(workspace,"_ttbar_TotalMC"+label,mj_shape["signal_mc"],"_mj",channel,wtagger,label,constraintlist_MC ) 
  model_ttbar_TotalMC_fail = MakeModelTTbarControlSample(workspace,"_ttbar_TotalMC"+label+"_failtau2tau1cut",mj_shape["signal_mc_fail"],"_mj",channel,wtagger,label,constraintlist_MC)                                                                                                                      
  model_TotalMC_fail       = RooAddPdf(("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel),("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel),RooArgList(model_ttbar_TotalMC_fail,model_bkg_TotalMC_fail,model_STop_fail,model_VV_fail,model_WJets_fail))
  model_TotalMC            = RooAddPdf(("model_TotalMC"+label+"_"+channel),("model_TotalMC"+label+"_"+channel),RooArgList(model_ttbar_TotalMC,model_bkg_TotalMC,model_STop,model_VV,model_WJets))
 
  getattr(workspace,'import')(model_TotalMC_fail)
  getattr(workspace,'import')(model_TotalMC)

# #///-------------------------------------------
# def DrawScaleFactorTTbarControlSample(RooWorkspace* workspace, std.map<,int> color_palet, const  & label, const  & channel, const  & wtagger,const double & ca8_ungroomed_pt_min, const double & ca8_ungroomed_pt_max, const  & sample){
#
#
#  #// ttSF = workspace.var("tt_scalefactor").getValV()
# ttSF =1
#   std.cout<< "Using tt scalefactor of " << ttSF << std.endl
#
#   rrv_mass_j = workspace.var("rrv_mass_j")
#
#   model_histpdf_STop                = workspace.pdf(("rdataset_STop"+label+"_"+channel+"_mj_histpdf"))
#   model_histpdf_VV                  = workspace.pdf(("rdataset_VV"+label+"_"+channel+"_mj_histpdf"))
#   model_histpdf_WJets               = workspace.pdf(("rdataset_WJets0"+label+"_"+channel+"_mj_histpdf"))
#   model_histpdf_TTbar_merged        = workspace.pdf(("rdataset_TTbar_realW"+label+"_"+channel+"_mj_histpdf"))
#   model_histpdf_TTbar_unmerged      = workspace.pdf(("rdataset_TTbar_fakeW"+label+"_"+channel+"_mj_histpdf"))
#   model_TTbar_STop_VV_WJets         = workspace.pdf(("model_TTbar_STop_VV_WJets"+label+"_"+channel))
#
#
#   model_histpdf_STop_fail           = workspace.pdf(("rdataset_STop"+label+"_"+"failtau2tau1cut_"+channel+"_mj_histpdf"))
#   model_histpdf_VV_fail             = workspace.pdf(("rdataset_VV"+label+"_"+"failtau2tau1cut_"+channel+"_mj_histpdf"))
#   model_histpdf_WJets_fail          = workspace.pdf(("rdataset_WJets0"+label+"_"+"failtau2tau1cut_"+channel+"_mj_histpdf"))
#   model_histpdf_TTbar_merged_fail   = workspace.pdf(("rdataset_TTbar_realW"+label+"_"+"failtau2tau1cut_"+channel+"_mj_histpdf"))
#   model_histpdf_TTbar_unmerged_fail = workspace.pdf(("rdataset_TTbar_fakeW"+label+"_"+"failtau2tau1cut_"+channel+"_mj_histpdf"))
#   model_TTbar_STop_VV_WJets_fail    = workspace.pdf(("model_TTbar_STop_VV_WJets_fail"+label+"_"+channel))
#
#
#   double scale_number_TTbar_STop_VV_WJets      = workspace.var("rrv_scale_number_TTbar_STop_VV_WJets").getValV()
#   double scale_number_TTbar_STop_VV_WJets_fail = workspace.var("rrv_scale_number_TTbar_STop_VV_WJets_fail").getValV()
#
#
#   model_data_fail = workspace.pdf(("model_data"+label+"_"+"failtau2tau1cut"+"_"+channel))
#   model_data      = workspace.pdf(("model_data"+label+"_"+channel))
#
#   RooCategory* category_p_f = workspace.cat(("category_p_f"+label+"_"+channel))
#
#   RooSimultaneous* simPdf_data = RooSimultaneous(("simPdf_data"+label+"_"+channel),("simPdf_data"+label+"_"+channel),*category_p_f)
#   simPdf_data.addPdf(*model_data,"pass")
#   simPdf_data.addPdf(*model_data_fail,"fail")
#   combData_p_f_data =  workspace.data(("combData_p_f_data"+label+"_"+channel))
#
#   simPdf_data.Print()
#   combData_p_f_data.Print()
#
#   model_TotalMC_fail = workspace.pdf(("model_TotalMC"+label+"_"+"failtau2tau1cut"+"_"+channel))
#   model_TotalMC      = workspace.pdf(("model_TotalMC"+label+"_"+channel))
#
#   RooSimultaneous* simPdf_TotalMC = RooSimultaneous(("simPdf_TotalMC"+label+"_"+channel),("simPdf_TotalMC"+label+"_"+channel),*category_p_f)
#   simPdf_TotalMC.addPdf(*model_TotalMC,"pass")
#   simPdf_TotalMC.addPdf(*model_TotalMC_fail,"fail")
#   combData_p_f_TotalMC =  workspace.data(("combData_p_f_TotalMC"+label+"_"+channel))
#
#   RooPlot* xframe_data = rrv_mass_j.frame( rt.RooFit.Bins(int(rrv_mass_j.getBins())))
#   RooPlot* xframe_data_fail = rrv_mass_j.frame( rt.RooFit.Bins(int(rrv_mass_j.getBins())))
#   xframe_data.GetYaxis().SetTitle(" Events / (5 GeV)")
#   xframe_data_fail.GetYaxis().SetTitle(" Events / (5 GeV)")
#   xframe_data     .GetYaxis().SetTitleOffset(1.39)
#   xframe_data_fail.GetYaxis().SetTitleOffset(1.39)
#
#   #//#############################################################
#   #///Plot total pass MC normalizing it to data
#   #//#############################################################
#
#
#   TString cut
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.pass",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data,rt.RooFit.Name("data_invisible"),rt.RooFit.Cut(cut.Data()),rt.RooFit.MarkerSize(1.5), rt.RooFit.DataError(rt.RooAbsData.SumW2), rt.RooFit.XErrorSize(0) )
#
#   #// Draw filled histograms
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar"), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s,%s",model_histpdf_TTbar_merged.GetName(),model_histpdf_TTbar_unmerged.GetName(),model_histpdf_STop.GetName(),model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar_realW"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar_realW"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s",model_histpdf_TTbar_unmerged.GetName(),model_histpdf_STop.GetName(),model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar_fakeW"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar_fakeW"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#
#   cut.Form("%s,%s,%s",model_histpdf_STop.GetName(),model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("STop"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["STop"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s,%s",model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("VV"),rt.RooFit.Components(cut.Data()), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["VV"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s",model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("WJets"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["WJets"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#
#   #// Draw black lines
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar_line_invisible"), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s,%s",model_histpdf_TTbar_merged.GetName(),model_histpdf_TTbar_unmerged.GetName(),model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar_realW_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s",model_histpdf_TTbar_unmerged.GetName(),model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("TTbar_fakeW_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s",model_histpdf_STop.GetName(), model_histpdf_VV.GetName(), model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("STop_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s",model_histpdf_VV.GetName(),model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("VV_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s",model_histpdf_WJets.GetName())
#   model_TTbar_STop_VV_WJets.plotOn(xframe_data,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets),rt.RooFit.Name("WJets_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#
#   #// plot data again
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.pass",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data,rt.RooFit.Name("data"), rt.RooFit.Cut(cut.Data()), rt.RooFit.MarkerSize(1.5), rt.RooFit.DataError(rt.RooAbsData.Poisson), rt.RooFit.XErrorSize(0))
#
#   #// plot mc fit function
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.pass",label,channel,label,channel)
#   simPdf_TotalMC.plotOn(xframe_data,rt.RooFit.Name("MC fit"),rt.RooFit.Slice(*category_p_f,"pass"), rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.LineStyle(kSolid), rt.RooFit.LineColor(rt.kRed))
#   cut.Form("model_bkg_TotalMC_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel,channel,channel,channel)
#   simPdf_TotalMC.plotOn(xframe_data,rt.RooFit.Name("mc fit bkg_invisible"),rt.RooFit.Slice(*category_p_f,"pass"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(rt.kRed), rt.RooFit.LineStyle(rt.kDashed))
#
#   #// plot data fit function
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.pass",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data,rt.RooFit.Name("data_invisible"),rt.RooFit.Cut(cut.Data()),rt.RooFit.MarkerSize(1.5),rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.XErrorSize(0))
#
#   simPdf_data.plotOn(xframe_data,rt.RooFit.Name("Data fit"),rt.RooFit.Slice(*category_p_f,"pass"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.LineStyle(kSolid), rt.RooFit.LineColor(kBlue))
#   cut.Form("model_bkg_data_%s_mj,model_STop_%s_mj,model_VV_%s_mj,model_WJets0_%s_mj",channel,channel,channel,channel)
#   simPdf_data.plotOn(xframe_data,rt.RooFit.Name("dat fit bkg_invisible"),rt.RooFit.Slice(*category_p_f,"pass"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.Components(cut.Data()), rt.RooFit.LineStyle(rt.kDashed), rt.RooFit.LineColor(kBlue))
#
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.fail",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data_fail,rt.RooFit.Name("data_invisible"), rt.RooFit.Cut(cut.Data()), rt.RooFit.MarkerSize(1.5), rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.XErrorSize(0))
#
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar"), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#
#   cut.Form("%s,%s,%s,%s,%s",model_histpdf_TTbar_merged_fail.GetName(),model_histpdf_TTbar_unmerged_fail.GetName(),model_histpdf_STop_fail.GetName(),model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar_realW"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar_realW"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s",model_histpdf_TTbar_unmerged_fail.GetName(),model_histpdf_STop_fail.GetName(),model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar_fakeW"),rt.RooFit.Components(cut.Data()),rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["TTbar_fakeW"+label]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#
#
#   cut.Form("%s,%s,%s",model_histpdf_STop_fail.GetName(),model_histpdf_VV_fail.GetName(),model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("STop"),rt.RooFit.Components(cut.Data()), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["STop"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s,%s",model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("VV"),rt.RooFit.Components(cut.Data()), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["VV"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#   cut.Form("%s", model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("WJets"),rt.RooFit.Components(cut.Data()), rt.RooFit.DrawOption("F"), rt.RooFit.FillColor(color_palet["WJets"]), rt.RooFit.LineColor(kBlack), rt.RooFit.VLines())
#
#   #//solid line
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar_line_invisible"), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s,%s",model_histpdf_TTbar_merged_fail.GetName(),model_histpdf_TTbar_unmerged_fail.GetName(),model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar_realW_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s,%s,%s",model_histpdf_TTbar_unmerged_fail.GetName(),model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("TTbar_fakeW_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#
#
#   cut.Form("%s,%s,%s",model_histpdf_STop_fail.GetName(), model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("STop_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s,%s",model_histpdf_VV_fail.GetName(), model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("VV_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#   cut.Form("%s", model_histpdf_WJets_fail.GetName())
#   model_TTbar_STop_VV_WJets_fail.plotOn(xframe_data_fail,rt.RooFit.Normalization(scale_number_TTbar_STop_VV_WJets_fail),rt.RooFit.Name("WJets_line_invisible"),rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(kBlack), rt.RooFit.LineWidth(2), rt.RooFit.VLines())
#
#
#   #//fail plots . plot MC fit
#
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.fail",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data_fail,rt.RooFit.Name("data_invisible"),rt.RooFit.Cut(cut.Data()),rt.RooFit.MarkerSize(1.5),rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.XErrorSize(0))
#   simPdf_TotalMC.plotOn(xframe_data_fail,rt.RooFit.Name("MC fit")    ,rt.RooFit.Slice(*category_p_f,"fail"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.LineStyle(kSolid), rt.RooFit.LineColor(rt.kRed))
#   cut.Form("model_bkg_TotalMC_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel,channel,channel,channel)
#   simPdf_TotalMC.plotOn(xframe_data_fail,rt.RooFit.Name("MC fit bkg"),rt.RooFit.Slice(*category_p_f,"fail"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_TotalMC),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.Components(cut.Data()), rt.RooFit.LineColor(rt.kRed), rt.RooFit.LineStyle(rt.kDashed))
#
#   #//fail plots . plot data fit
#   cut.Form("category_p_f%s_%s==category_p_f%s_%s.fail",label,channel,label,channel)
#   combData_p_f_data.plotOn(xframe_data_fail,rt.RooFit.Name("data_invisible"),rt.RooFit.Cut(cut.Data()),rt.RooFit.MarkerSize(1.5),rt.RooFit.DataError(rt.RooAbsData.SumW2),rt.RooFit.XErrorSize(0))
#
#   simPdf_data.plotOn(xframe_data_fail,rt.RooFit.Name("Data fit"),rt.RooFit.Slice(*category_p_f,"fail"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),rt.RooFit.NormRange("controlsample_fitting_range"),rt.RooFit.LineStyle(kSolid), rt.RooFit.LineColor(kBlue))
#   simPdf_data.plotOn(xframe_data_fail,rt.RooFit.Name("data_fit_invisible"),rt.RooFit.Slice(*category_p_f,"fail"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),rt.RooFit.NormRange("controlsample_fitting_range"))
#   cut.Form("model_bkg_data_failtau2tau1cut_%s_mj,model_STop_failtau2tau1cut_%s_mj,model_VV_failtau2tau1cut_%s_mj,model_WJets0_failtau2tau1cut_%s_mj",channel,channel,channel,channel)
#   simPdf_data.plotOn(xframe_data_fail,rt.RooFit.Name("data fit bkg"),rt.RooFit.Slice(*category_p_f,"fail"),rt.RooFit.ProjWData(RooArgSet(*category_p_f),*combData_p_f_data),rt.RooFit.NormRange("controlsample_fitting_range"), rt.RooFit.Components(cut.Data()),rt.RooFit.LineStyle(rt.kDashed), rt.RooFit.LineColor(kBlue))
#
#   #// #//signal window
# #//   TLine* lowerLine = TLine(65,0.,65,xframe_data.GetMaximum()*0.7) lowerLine.SetLineWidth(2) lowerLine.SetLineColor(kBlack) lowerLine.SetLineStyle(9)
# #//   TLine* upperLine = TLine(105,0.,105,xframe_data.GetMaximum()*0.7) upperLine.SetLineWidth(2) upperLine.SetLineColor(kBlack) upperLine.SetLineStyle(9)
# #//   xframe_data.addObject(lowerLine) xframe_data.addObject(upperLine)
# #//   lowerLine = TLine(65,0.,65,xframe_data_fail.GetMaximum()*0.7) lowerLine.SetLineWidth(2) lowerLine.SetLineColor(kBlack) lowerLine.SetLineStyle(9)
# #//   upperLine = TLine(105,0.,105,xframe_data_fail.GetMaximum()*0.7) upperLine.SetLineWidth(2) upperLine.SetLineColor(kBlack) upperLine.SetLineStyle(9)
# #//   xframe_data_fail.addObject(lowerLine) xframe_data_fail.addObject(upperLine)
#
#   #// TLegend* leg_data = legend4Plot(xframe_data,0,-0.2,0.07,0.04,0.,1,channel)
#   #// xframe_data.addObject(leg_data)
#   #// TLegend* leg_data_fail = legend4Plot(xframe_data_fail,0,-0.2,0.07,0.04,0.,1,channel)
#   #// xframe_data_fail.addObject(leg_data)
#
#   TLegend*  theLeg = TLegend(0.3885213,0.6640827,0.9774937,0.8992248,"","NDC")
#   theLeg.SetName("theLegend")
#   #// theLeg.SetNColumns(3)
#   theLeg.SetNColumns(2)
#   theLeg.SetFillColor(0)
#   theLeg.SetFillStyle(0)
#   theLeg.SetBorderSize(0)
#   theLeg.SetLineColor(0)
#   theLeg.SetLineWidth(0)
#   theLeg.SetLineStyle(0)
#   theLeg.SetTextSize(0.039)
#   theLeg.SetTextFont(42)
#   theLeg.AddEntry(xframe_data.findObject("data")       ,"CMS data"            ,"PLE")
#   theLeg.AddEntry(xframe_data.findObject("TTbar_fakeW"),"t#bar{t} (unmerged)" ,"F")
#   theLeg.AddEntry(xframe_data.findObject("Data fit")   ,"Data fit"            ,"L")
#   theLeg.AddEntry(xframe_data.findObject("STop")       ,"Single top"          ,"F")
#   theLeg.AddEntry(xframe_data.findObject("MC fit")     ,"MC fit"              ,"L")
#   theLeg.AddEntry(xframe_data.findObject("WJets")      ,"W+jets"              ,"F")
#   theLeg.AddEntry(xframe_data.findObject("TTbar_realW"),"t#bar{t} (merged)"   ,"F")
#   theLeg.AddEntry(xframe_data.findObject("VV")         ,"WW/WZ/ZZ"            ,"F")
#
#   xframe_data.addObject(theLeg)
#   xframe_data_fail.addObject(theLeg)
#
#
#    tmp_channel = "el"
#   if(workspace.var(("rrv_mean1_gaus_ttbar_data"+label+"_"+"el"))) tmp_channel = "el"
#   else tmp_channel = "mu"
#
#   rrv_mean_gaus_data     = workspace.var(("rrv_mean1_gaus_ttbar_data"+label+"_"+tmp_channel))
#   rrv_sigma_gaus_data    = workspace.var(("rrv_sigma1_gaus_ttbar_data"+label+"_"+tmp_channel))
#   rrv_mean_gaus_TotalMC  = workspace.var(("rrv_mean1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel))
#   rrv_sigma_gaus_TotalMC = workspace.var(("rrv_sigma1_gaus_ttbar_TotalMC"+label+"_"+tmp_channel))
#
#   if(rrv_mean_gaus_TotalMC){
#     TString latex  latex.Form("Mean_{MC } = %3.1f #pm %2.1f",rrv_mean_gaus_TotalMC.getVal(),rrv_mean_gaus_TotalMC.getError())
#     TLatex* tl_MC_mean  = TLatex(0.25 ,0.62,latex.Data())
#     latex.Form("Sigma_{MC }= %2.1f #pm %2.1f",rrv_sigma_gaus_TotalMC.getVal(),rrv_sigma_gaus_TotalMC.getError())
#     TLatex* tl_MC_sigma = TLatex(0.25 ,0.57,latex.Data())
#     tl_MC_mean.SetNDC() tl_MC_sigma.SetNDC()
#     tl_MC_mean.SetTextSize(0.03)
#     tl_MC_sigma.SetTextSize(0.03)
#     #// xframe_data.addObject(tl_MC_mean)
#     #// xframe_data.addObject(tl_MC_sigma)
#
#   }
#   if(rrv_mean_gaus_data){
#     TString latex  latex.Form("Mean_{data} = %3.1f #pm %2.1f",rrv_mean_gaus_data.getVal(),rrv_mean_gaus_data.getError())
#     TLatex* tl_data_mean  = TLatex(0.25 ,0.62,latex.Data())
#     latex.Form("Sigma_{data}= %2.1f #pm %2.1f",rrv_sigma_gaus_data.getVal(),rrv_sigma_gaus_data.getError())
#     TLatex* tl_data_sigma = TLatex(0.25 ,0.57,latex.Data())
#     tl_data_mean.SetNDC() tl_data_sigma.SetNDC()
#     tl_data_mean.SetTextSize(0.03)
#     tl_data_sigma.SetTextSize(0.03)
#     #// xframe_data.addObject(tl_data_mean) xframe_data.addObject(tl_data_sigma)
#   }
#
#   xframe_data.GetYaxis().SetRangeUser(1e-2,xframe_data.GetMaximum()*1.4)
#   xframe_data_fail.GetYaxis().SetRangeUser(1e-2,xframe_data_fail.GetMaximum()*1.4)
#
#   TString nameDir  nameDir.Form("plots/TotalFit/",wtagger)
#   TString namePlot namePlot.Form("TotalFit")
#   draw_canvas(xframe_data,(nameDir),(namePlot),channel,GetLumi(),0,0,0)
#   namePlot.Form("TotalFit_fail")
#   draw_canvas(xframe_data_fail,(nameDir),(namePlot),channel,GetLumi(),0,0,0)
#
# }