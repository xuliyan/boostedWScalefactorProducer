import ROOT
from ROOT import RooHistPdf,RooRealVar,RooExponential,RooGaussian,RooExtendPdf,RooAddPdf,RooArgList,RooArgSet,RooFormulaVar

ROOT.gSystem.Load("PDFs/HWWLVJRooPdfs_cxx.so")

def addConstraint(workspace,rrv_x, x_mean, x_sigma, ConstraintsList):
	rrv_x_mean = RooRealVar(rrv_x.GetName()+"_mean",rrv_x.GetName()+"_mean",x_mean )
	rrv_x_sigma = RooRealVar(rrv_x.GetName()+"_sigma",rrv_x.GetName()+"_sigma",x_sigma )
	constrainpdf_x = RooGaussian("constrainpdf_"+rrv_x.GetName(),"constrainpdf_"+rrv_x.GetName(),rrv_x, rrv_x_mean, rrv_x_sigma)
	ConstraintsList.append(constrainpdf_x.GetName())
	getattr(workspace,"import")(constrainpdf_x)
	return constrainpdf_x
			
def change_dataset_to_histpdf(workspace,x,dataset):
  datahist = dataset.binnedClone(dataset.GetName()+"_binnedClone",dataset.GetName()+"_binnedClone")
  histpdf  = RooHistPdf(dataset.GetName()+"_histpdf",dataset.GetName()+"_histpdf",RooArgSet(x),datahist)
  getattr(workspace,'import')(histpdf)

def MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label, channel,constraint,ismc=1):
	
	print "Making general PDF for ","model_pdf"+label+"_"+channel+spectrum 
	rrv_x = workspace.var("rrv_mass_j")

	if model == "ErfExp":
		rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+channel+spectrum,"rrv_c_ErfExp"+label+"_"+channel+spectrum,-0.026)#,-0.05, 0.05)
  		rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum,"rrv_offset_ErfExp"+label+"_"+channel+spectrum,41.,0.,100)
  		rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+channel+spectrum,"rrv_width_ErfExp"+label+"_"+channel+spectrum,30.,1.,100.)
  		model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)

	if model == "ExpGaus":
		
		if label.find("_STop_failtau2tau1cut")!=-1:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.03,-0.5,0.5)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,84,60,120) #Too narrow limits here often lead to error!! eg max 80
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,7,4,60)
		
		
		elif label.find("fail")!=-1:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.05,-0.5,0.5)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,80,70.,95.)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,7,6.,20.)
			rrv_high        = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.,0.6,1.)
        
		else:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.01,-1.,0.)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,84 ,75,89) 
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,11,4,15)
			rrv_high        = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.7,0.,1.)
		  
		exp		= RooExponential("exp"+label+"_"+channel+spectrum,"exp"+label+"_"+channel+spectrum,rrv_x,rrv_c_Exp)
		gaus		= RooGaussian("gaus"+label+"_"+channel+spectrum,"gaus"+label+"_"+channel+spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
		model_pdf	= RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))
     
	if model == "ErfExp_ttbar" or model == "ErfExp_ttbar_failtau2tau1cut":
		
		if model.find("failtau2tau1cut")!=-1:
			c0_tmp     = -3.5160e-02
			offset_tmp = 6.0935e+01 
			width_tmp  = 3.9583e+01 
			c0_tmp_err = 1.46e-02 
			offset_tmp_err = 4.92e+01 
			width_tmp_err = 4.69e+00
			
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp )
			
		else:	
			c0_tmp         =  -2.1046e-02      
			offset_tmp     =  80.1       
			width_tmp      =  29.5  
			c0_tmp_err     = 6.83e-03
			offset_tmp_err = 9.3
			width_tmp_err  = 2.9
			
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp, offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp , width_tmp-10, width_tmp+10)
 			 
		rrv_c_ErfExp  = RooRealVar("rrv_c_ErfExp"  +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp    , c0_tmp-4e-2, c0_tmp+4e-2 )
		model_pdf     = ROOT.RooErfExpPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum ,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)
		gaus1         = addConstraint(workspace,rrv_c_ErfExp,rrv_c_ErfExp.getVal(),c0_tmp_err,constraint)
      
		if model.find("failtau2tau1cut")==-1:
			gaus2 = addConstraint(workspace,rrv_offset_ErfExp,rrv_offset_ErfExp.getVal(),offset_tmp_err,constraint)
			gaus3 = addConstraint(workspace,rrv_width_ErfExp,rrv_width_ErfExp.getVal(),width_tmp_err,constraint)
	
	if model == "GausErfExp_ttbar" or model =="GausErfExp_ttbar_failtau2tau1cut":
		
		if model.find("fail")!=-1:
			c0_tmp      = -0.05
			offset_tmp  = 199. 
			offset_tmp_err = 20
			width_tmp   = 64.8 
			frac_tmp    = 0.15 
			mean1_tmp   = 80.0 
			sigma1_tmp  = 9.41
			if label.find("data")!=-1: 
				gaus1 = workspace.pdf("gaus1_ttbar_data_"+channel+spectrum)
			else:					   
				gaus1 = workspace.pdf("gaus1_ttbar_TotalMC_"+channel+spectrum)
      
		else:
			c0_tmp     = -2.1
			offset_tmp =  8.0
			offset_tmp_err = 20
			width_tmp  =  2.9
			frac_tmp    = 0.80
			mean1_tmp  = 80.8
			sigma1_tmp =  8.4
			rangeMean  =  5. 
			rangeWidth =  5. 
			
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum ,mean1_tmp, mean1_tmp-rangeMean, mean1_tmp+rangeMean)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,sigma1_tmp, sigma1_tmp-rangeWidth,sigma1_tmp+rangeWidth )
			gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
			
		rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"     +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 )
		rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err*4,offset_tmp+offset_tmp_err*4)
		rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-10, width_tmp+10)
		
		rrv_frac = RooRealVar("rrv_frac"+label+"_"+channel+spectrum ,"rrv_frac"+label+"_"+channel+spectrum ,frac_tmp,0.,1.)

		rrv_c_ErfExp.setConstant(ROOT.kTRUE)
		rrv_offset_ErfExp.setConstant(ROOT.kTRUE)
		rrv_width_ErfExp.setConstant(ROOT.kTRUE)
		
		erfExp    = ROOT.RooErfExpPdf("erfExp"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)
		model_pdf = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,gaus1,erfExp,rrv_frac)

	if model == "Exp" :
		rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.030, -2., 0.05)
		model_pdf = RooExponential("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_Exp)		
    
	if model == "Gaus":
		rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus" +label+"_"+channel+spectrum ,"rrv_mean1_gaus" +label+"_"+channel+spectrum,80,40,100)
		rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,7,0.,15)  
		model_pdf        = RooGaussian("gaus"+label+"_"+channel+spectrum,"gaus"+label+"_"+channel+spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
    
	if model == "ErfExpGaus_sp":
		
		rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"    +channel+spectrum,"rrv_c_ErfExp"    +label+"_"+channel+spectrum,-0.04,-0.2,0.)
		rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+channel+spectrum,"rrv_width_ErfExp"+label+"_"+channel+spectrum,30.,0.,300.)
		rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"  +channel+spectrum,"rrv_mean1_gaus"  +label+"_"+channel+spectrum,80,60,100)
		rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_" +channel+spectrum,"rrv_sigma1_gaus" +label+"_"+channel+spectrum,7,10.,40)
		erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+channel+spectrum,"erfExp"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp)
		gaus             = RooGaussian ("gaus"+label+"_"+channel+spectrum  ,"gaus"+label+"_"+channel+spectrum  , rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
		rrv_high   = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.3,0.0,0.8)
		model_pdf  = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,erfExp,gaus,rrv_high)

	getattr(workspace,'import')(model_pdf)
	return workspace.pdf("model_pdf"+label+"_"+channel+spectrum)
			
def MakeExtendedModel(workspace, label, model,spectrum, channel, wtagger_label,constraint):
	
	print "Making extended model with name: " , "model"+label+"_"+channel+spectrum
	rrv_number = RooRealVar("rrv_number"+label+"_"+channel+spectrum,"rrv_number"+label+"_"+channel+spectrum,500.,0.,1e5)
	model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label,channel,constraint)

	model_extended = RooExtendPdf("model"+label+"_"+channel+spectrum,"model"+label+"_"+channel+spectrum,model_pdf, rrv_number)
	
	getattr(workspace,'import')(rrv_number)
	getattr(workspace,'import')(model_extended)
	return workspace.pdf("model"+label+"_"+channel+spectrum)

def fixParameters(workspace,label,channel,fix=1):
	name = "rdataset%s_%s_mj"%(label,channel)
	rdataset_General_mj = workspace.data(name)
	model_General = workspace.pdf("model"+label+"_"+channel+"_mj")
	rdataset_General_mj.Print()
	model_General.Print()
	
	parameters_General = model_General.getParameters(rdataset_General_mj)
	par=parameters_General.createIterator()
	par.Reset()
	param=par.Next()
	while (param):
		param.setConstant(ROOT.kTRUE)
		param=par.Next()
	print "Failing for " ,"model"+label+"_"+channel+"_mj"
	return workspace.pdf("model"+label+"_"+channel+"_mj")
			
def makeTTbarModel(workspace,label, model,channel, wtagger, constraint=[], spectrum="_mj"):
	 
	 info =""
	 if label.find("_ttbar_data")!=-1 and label.find("fail")==-1:
		 rrv_number_total = RooRealVar("rrv_number_total_ttbar_data"+info+"_"+channel,"rrv_number_total_ttbar_data"+info+"_"+channel,500,0.,1e7)
		 eff_ttbar = RooRealVar("eff_ttbar_data"+info+"_"+channel,"eff_ttbar_data"+info+"_"+channel,0.7,0.3,0.9)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total))

	 elif label.find("_ttbar_data")!=-1 and label.find("fail")!=-1:
		 rrv_number_total = workspace.var("rrv_number_total_ttbar_data"+info+"_"+channel)
		 eff_ttbar        = workspace.var("eff_ttbar_data"+info+"_"+channel)
		 rrv_number       = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total))
		 
	 elif label.find("_ttbar_TotalMC")!=-1 and label.find("fail")==-1:
		 rrv_number_total = RooRealVar("rrv_number_total_ttbar_TotalMC"+info+"_"+channel,"rrv_number_total_ttbar_TotalMC"+info+"_"+channel,500,0.,1e7)
		 eff_ttbar = RooRealVar("eff_ttbar_TotalMC"+info+"_"+channel,"eff_ttbar_TotalMC"+info+"_"+channel,0.7,0.3,0.9)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total))
		 
	 elif label.find("_ttbar_TotalMC")!=-1 and label.find("fail")!=-1:
		 rrv_number_total = workspace.var("rrv_number_total_ttbar_TotalMC"+info+"_"+channel)
		 eff_ttbar = workspace.var("eff_ttbar_TotalMC"+info+"_"+channel)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total)) 
		 
	 model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger,channel,constraint)
	 model = RooExtendPdf("model"+label+"_"+channel+spectrum,"model"+label+"_"+channel+spectrum, model_pdf, rrv_number)
	 getattr(workspace,"import")(model)
	 return workspace.pdf("model"+label+"_"+channel+spectrum)
	