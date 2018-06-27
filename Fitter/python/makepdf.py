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

def MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label, channel,constraint,peak="W"):
	
	print "Making general PDF for ","model_pdf"+label+"_"+channel+spectrum 
	rrv_x = workspace.var("rrv_mass_j")
	gaus_means  = 8.2653e+01
	gaussigmas   = 7
	if peak == "t":
	    gaus_means  = 180
	    gaussigmas   = 18
    
	if model == "ErfExp":
		rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"+label+"_"+channel+spectrum,"rrv_c_ErfExp"+label+"_"+channel+spectrum,-0.026,-0.05, 0.05)
  		rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum,"rrv_offset_ErfExp"+label+"_"+channel+spectrum,41.,0.,100)
  		rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp"+label+"_"+channel+spectrum,"rrv_width_ErfExp"+label+"_"+channel+spectrum,30.,1.,100.)
  		model_pdf         = ROOT.RooErfExpPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)

	if model == "ExpGaus":
		
		if label.find("_STop_failtau2tau1cut")!=-1:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.03,-0.5,0.5)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,gaus_means,gaus_means*.8,gaus_means*1.2) #Too narrow limits here often lead to error!! eg max 80
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,gaussigmas,gaussigmas*.5,gaussigmas*1.5)
		
		
		elif label.find("fail")!=-1:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.05,-0.5,0.5)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,gaus_means,gaus_means*.8,gaus_means*1.2)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,gaussigmas,gaussigmas*.5,gaussigmas*1.5)
			rrv_high        = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.3,0.3,1.)
        
		else:
			rrv_c_Exp       = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.01,-1.,0.)
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum,"rrv_mean1_gaus"+label+"_"+channel+spectrum,gaus_means,gaus_means*.8,gaus_means*1.2) 
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,gaussigmas,gaussigmas*.5,gaussigmas*1.5)
			rrv_high        = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.7,0.,1.)
		  
		exp		= RooExponential("exp"+label+"_"+channel+spectrum,"exp"+label+"_"+channel+spectrum,rrv_x,rrv_c_Exp)
		gaus		= RooGaussian("gaus"+label+"_"+channel+spectrum,"gaus"+label+"_"+channel+spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
		model_pdf	= RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,RooArgList(exp,gaus),RooArgList(rrv_high))
     
	if model == "ErfExp_ttbar" or model == "ErfExp_ttbar_ddt" or model == "ErfExp_ttbar_failtau2tau1cut" or model == "ErfExp_ttbar_failtau2tau1cut_fitMC" or model == "ErfExp_ttbar_failtau2tau1cut_ddt":
		print "PRINTING MODEL " ,model
        if peak == "Wt":
            c0_tmp     =  -3.5239e-02 ; c0_tmp_err     = 4.52e-03;
            offset_tmp =  9.3695e+01  ; offset_tmp_err = 3.53e+00;
            width_tmp  =  3.2006e+01  ; width_tmp_err  = 1.84e+00;
            rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,93, 0,200) 
            rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,32,30,100)#,width_tmp-10, width_tmp+10)
            rrv_c_ErfExp  = RooRealVar("rrv_c_ErfExp"  +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,-0.02    , -2, 2 )
            
        
        elif model.find("failtau2tau1cut")!=-1:
        	c0_tmp = -2.6259e-02    ; c0_tmp_err = 1.46e-02;
        	offset_tmp = 9.1220e+01 ; offset_tmp_err = 4.92;
        	width_tmp = 3.1356e+01  ; width_tmp_err = 4.69e+00;
        	rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp, offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err) #Can be fixed to avoid downwoards slope towards zero
        	rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp,width_tmp-10, width_tmp+10)
        	rrv_c_ErfExp  = RooRealVar("rrv_c_ErfExp"  +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp, c0_tmp-4e-2, c0_tmp+4e-2 )
        	if model.find("_ddt")!=-1:	 
        	    rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp, 0.,100.) #Can be fixed to avoid downwoards slope towards zero
        	    rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp,width_tmp-10, width_tmp+10)
        	
        else:	
        	c0_tmp     =  -3.5239e-02 ; c0_tmp_err     = 4.52e-03;
        	offset_tmp =  9.3695e+01  ; offset_tmp_err = 3.53e+00;
        	width_tmp  =  3.2006e+01  ; width_tmp_err  = 1.84e+00;
        	
        	rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp, offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err)
        	rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-10, width_tmp+10)
        	rrv_c_ErfExp  = RooRealVar("rrv_c_ErfExp"  +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp    , c0_tmp-4e-2, c0_tmp+4e-2 )	 
        	rrv_c_ErfExp.Print()
        	print "PRINTED"
        
        
        model_pdf     = ROOT.RooErfExpPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum ,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)
        gaus1         = addConstraint(workspace,rrv_c_ErfExp,rrv_c_ErfExp.getVal(),c0_tmp_err,constraint)
        
        if model.find("failtau2tau1cut")==-1:
        	gaus2 = addConstraint(workspace,rrv_offset_ErfExp,rrv_offset_ErfExp.getVal(),offset_tmp_err,constraint)
        	gaus3 = addConstraint(workspace,rrv_width_ErfExp,rrv_width_ErfExp.getVal(),width_tmp_err,constraint)
        
	if model == "GausErfExp_ttbar" or model == "GausErfExp_ttbar_fitMC" or model =="GausErfExp_ttbar_failtau2tau1cut" or model =="GausErfExp_ttbar_failtau2tau1cut_fitMC":

		frac_tmp = 1.; #Use only Gaussian component to model real W-jets, ErfExp set to zero by setting fraction to 1!
		
		if model.find("fitMC")!=-1:
			mean1_tmp = 8.2653e+01;
			sigma1_tmp = 7.5932e+00;
			rangeMean = 8. ;
			rangeWidth = 5. ;
			c0_tmp     = -2.7180e-02 ; c0_tmp_err     = 6.83e-03;
			offset_tmp =  8.6888e+01 ; offset_tmp_err = 19.35e+00;
			width_tmp  =  3.6383e+01 ; width_tmp_err  = 1.34e+00;
				
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum   ,gaus_means,gaus_means*.8,gaus_means*1.2)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,gaussigmas,gaussigmas*.5,gaussigmas*2.5 )
			gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)     
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, 30., 100.)
		    
		        
		elif model.find("fail")!=-1:
			c0_tmp     = -7.3903e-03 ; c0_tmp_err = 1.46e-02;
			offset_tmp = 8.9205e+01  ; offset_tmp_err = 4.92e+01;
			width_tmp  = 7.6305e+01  ; width_tmp_err = 14.69;
			mean1_tmp  = 8.2278e+01;
			sigma1_tmp = 9.6909e+00;
            #frac_tmp = 5.8066e-01
			if label.find("data")!=-1:
			  gaus1 = workspace.pdf("gaus1_ttbar_data_"+channel+spectrum)
			else:
			  gaus1 = workspace.pdf("gaus1_ttbar_TotalMC_"+channel+spectrum)
            # rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum   ,gaus_means,gaus_means*.8,gaus_means*1.2)
            # rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,gaussigmas,gaussigmas*.5,gaussigmas*1.5 )
            # gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
            
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+2*offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-width_tmp_err, width_tmp+2*width_tmp_err)
		
		else:
			mean1_tmp = 8.4666e+01;
			sigma1_tmp = 8.5006e+00;
			rangeMean = 8. ;
			rangeWidth = 5. ;
			c0_tmp     = -5.8699e-02 ; c0_tmp_err     = 6.83e-03;
			offset_tmp =  1.0199e+02 ; offset_tmp_err = 9.35e+00;
			width_tmp  =  3.0009e+01 ; width_tmp_err  = 1.34e+00;
				
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum   ,90,gaus_means*.8,gaus_means*1.2)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,5,0,6)
			gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)     
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-10, width_tmp+10)
			
		
		
		if model.find("fitMC")!=-1:
		    rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"     +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp,c0_tmp-1, c0_tmp+1 )
		    
		    rrv_frac = RooRealVar("rrv_frac"+label+"_"+channel+spectrum ,"rrv_frac"+label+"_"+channel+spectrum ,frac_tmp,0.,1.)
		
		    
		else:
		    rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"     +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 )
		    
		    rrv_frac = RooRealVar("rrv_frac"+label+"_"+channel+spectrum ,"rrv_frac"+label+"_"+channel+spectrum ,frac_tmp)
		
		    rrv_frac.setConstant(ROOT.kTRUE)  
		    rrv_c_ErfExp.setConstant(ROOT.kTRUE)#Force to constant, if not RooFit will try to fit them despite frac==0 for ErfExp
		    rrv_offset_ErfExp.setConstant(ROOT.kTRUE)
		    rrv_width_ErfExp.setConstant(ROOT.kTRUE)
		
		erfExp    = ROOT.RooErfExpPdf("erfExp"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)
		model_pdf = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,RooArgList(gaus1,erfExp),RooArgList(rrv_frac),1)

	if model == "Gaus2ErfExp_ttbar" or model == "Gaus2ErfExp_ttbar_fitMC" or model =="Gaus2ErfExp_ttbar_failtau2tau1cut" or model =="Gaus2ErfExp_ttbar_failtau2tau1cut_fitMC":
		
		frac_tmp = 1.0; #Use only Gaussian component to model real W-jets, ErfExp set to zero by setting fraction to 1!
		
		if model.find("fitMC")!=-1:
			mean1_tmp = 8.2653e+01;
			sigma1_tmp = 7.5932e+00;
			rangeMean = 8. ;
			rangeWidth = 5. ;
			c0_tmp     = -2.7180e-02 ; c0_tmp_err     = 6.83e-03;
			offset_tmp =  8.6888e+01 ; offset_tmp_err = 19.35e+00;
			width_tmp  =  3.6383e+01 ; width_tmp_err  = 1.34e+00;
				
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum   ,80.,75.,85.)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,7.6,5.,10. )
			gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)     
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, 30., 100.)
            
			rrv_mean2_gaus  = RooRealVar("rrv_mean2_gaus"+label+"_"+channel+spectrum ,"rrv_mean2_gaus"+label+"_"+channel+spectrum   ,170,150,180)
			rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+channel+spectrum ,"rrv_sigma2_gaus"+label+"_"+channel+spectrum ,13,10.,20. )
			gaus2 = RooGaussian("gaus2"+label+"_"+channel+spectrum ,"gaus2"+label+"_"+channel+spectrum ,rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus)
		    
		        
		elif model.find("fail")!=-1:
			c0_tmp     = -7.3903e-03 ; c0_tmp_err = 1.46e-02;
			offset_tmp = 8.9205e+01  ; offset_tmp_err = 4.92e+01;
			width_tmp  = 7.6305e+01  ; width_tmp_err = 14.69;
			mean1_tmp  = 8.2278e+01;
			sigma1_tmp = 9.6909e+00;
            #frac_tmp = 5.8066e-01
			if label.find("data")!=-1:
			  gaus1 = workspace.pdf("gaus1_ttbar_data_"+channel+spectrum)
			else:
			  gaus1 = workspace.pdf("gaus1_ttbar_TotalMC_"+channel+spectrum)
			
			rrv_mean2_gaus  = RooRealVar("rrv_mean2_gaus"+label+"_"+channel+spectrum ,"rrv_mean2_gaus"+label+"_"+channel+spectrum   ,gaus_means+100,(gaus_means+100)*.8,(gaus_means+100)*1.2)
			rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+channel+spectrum ,"rrv_sigma2_gaus"+label+"_"+channel+spectrum ,gaussigmas+10,(gaussigmas+10)*.5,(gaussigmas+10)*1.5 )
			gaus2 = RooGaussian("gaus2"+label+"_"+channel+spectrum ,"gaus2"+label+"_"+channel+spectrum ,rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus)            
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+2*offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-width_tmp_err, width_tmp+2*width_tmp_err)
		
		else:
			mean1_tmp = 8.4666e+01;
			sigma1_tmp = 8.5006e+00;
			rangeMean = 8. ;
			rangeWidth = 5. ;
			c0_tmp     = -5.8699e-02 ; c0_tmp_err     = 6.83e-03;
			offset_tmp =  1.0199e+02 ; offset_tmp_err = 9.35e+00;
			width_tmp  =  3.0009e+01 ; width_tmp_err  = 1.34e+00;
				
			rrv_mean1_gaus  = RooRealVar("rrv_mean1_gaus"+label+"_"+channel+spectrum ,"rrv_mean1_gaus"+label+"_"+channel+spectrum   ,gaus_means,gaus_means*.8,gaus_means*1.2)
			rrv_sigma1_gaus = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum ,gaussigmas,gaussigmas*.5,gaussigmas*1.5)
			gaus1 = RooGaussian("gaus1"+label+"_"+channel+spectrum ,"gaus1"+label+"_"+channel+spectrum ,rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)     
			rrv_offset_ErfExp = RooRealVar("rrv_offset_ErfExp"+label+"_"+channel+spectrum ,"rrv_offset_ErfExp"+label+"_"+channel+spectrum ,offset_tmp,offset_tmp-offset_tmp_err,offset_tmp+offset_tmp_err)
			rrv_width_ErfExp  = RooRealVar("rrv_width_ErfExp" +label+"_"+channel+spectrum ,"rrv_width_ErfExp" +label+"_"+channel+spectrum ,width_tmp, width_tmp-10, width_tmp+10)
			
			rrv_mean2_gaus  = RooRealVar("rrv_mean2_gaus"+label+"_"+channel+spectrum ,"rrv_mean2_gaus"+label+"_"+channel+spectrum   ,gaus_means+100,(gaus_means+100)*.8,(gaus_means+100)*1.2)
			rrv_sigma2_gaus = RooRealVar("rrv_sigma2_gaus"+label+"_"+channel+spectrum ,"rrv_sigma2_gaus"+label+"_"+channel+spectrum ,gaussigmas+10,(gaussigmas+10)*.5,(gaussigmas+10)*1.5 )
			gaus2 = RooGaussian("gaus2"+label+"_"+channel+spectrum ,"gaus2"+label+"_"+channel+spectrum ,rrv_x,rrv_mean2_gaus,rrv_sigma2_gaus)
			
		
		
		if model.find("fitMC")!=-1:
		    rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"     +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp,c0_tmp-1, c0_tmp+1 )
		    
		    rrv_frac = RooRealVar("rrv_frac"+label+"_"+channel+spectrum ,"rrv_frac"+label+"_"+channel+spectrum ,frac_tmp,0.,1.)
		
		    
		else:
		    rrv_c_ErfExp      = RooRealVar("rrv_c_ErfExp"     +label+"_"+channel+spectrum ,"rrv_c_ErfExp"     +label+"_"+channel+spectrum ,c0_tmp,c0_tmp-4e-2, c0_tmp+4e-2 )
		    
		    rrv_frac = RooRealVar("rrv_frac"+label+"_"+channel+spectrum ,"rrv_frac"+label+"_"+channel+spectrum ,frac_tmp)
		
		    rrv_frac.setConstant(ROOT.kTRUE)  
		    rrv_c_ErfExp.setConstant(ROOT.kTRUE)#Force to constant, if not RooFit will try to fit them despite frac==0 for ErfExp
		    rrv_offset_ErfExp.setConstant(ROOT.kTRUE)
		    rrv_width_ErfExp.setConstant(ROOT.kTRUE)
		
		erfExp    = ROOT.RooErfExpPdf("erfExp"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_offset_ErfExp,rrv_width_ErfExp)
        # model_pdf = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,RooArgList(gaus1,erfExp),RooArgList(rrv_frac),1)
		model_pdf = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,RooArgList(gaus1,gaus2),RooArgList(rrv_frac),1)

	if model == "Exp" :
		rrv_c_Exp = RooRealVar("rrv_c_Exp"+label+"_"+channel+spectrum,"rrv_c_Exp"+label+"_"+channel+spectrum,-0.030, -2., 0.05)
		model_pdf = RooExponential("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,rrv_x,rrv_c_Exp)		
    
	if model == "Gaus":
		rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus" +label+"_"+channel+spectrum ,"rrv_mean1_gaus" +label+"_"+channel+spectrum,gaus_means,gaus_means*.8,gaus_means*1.2)
		rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_"+channel+spectrum ,"rrv_sigma1_gaus"+label+"_"+channel+spectrum,gaussigmas,gaussigmas*.5,gaussigmas*1.5)  
		model_pdf        = RooGaussian("gaus"+label+"_"+channel+spectrum,"gaus"+label+"_"+channel+spectrum, rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
    
	if model == "ErfExpGaus_sp":
		
		rrv_c_ErfExp     = RooRealVar("rrv_c_ErfExp"+label+"_"    +channel+spectrum,"rrv_c_ErfExp"    +label+"_"+channel+spectrum,-0.04,-1.,1.)
		rrv_width_ErfExp = RooRealVar("rrv_width_ErfExp"+label+"_"+channel+spectrum,"rrv_width_ErfExp"+label+"_"+channel+spectrum,30.,0.,400.)
		rrv_mean1_gaus   = RooRealVar("rrv_mean1_gaus"+label+"_"  +channel+spectrum,"rrv_mean1_gaus"  +label+"_"+channel+spectrum,gaus_means,gaus_means*.8,gaus_means*1.2)
		rrv_sigma1_gaus  = RooRealVar("rrv_sigma1_gaus"+label+"_" +channel+spectrum,"rrv_sigma1_gaus" +label+"_"+channel+spectrum,gaussigmas,gaussigmas*.5,gaussigmas*1.5)
		erfExp           = ROOT.RooErfExpPdf("erfExp"+label+"_"+channel+spectrum,"erfExp"+label+"_"+channel+spectrum,rrv_x,rrv_c_ErfExp,rrv_mean1_gaus,rrv_width_ErfExp)
		gaus             = RooGaussian ("gaus"+label+"_"+channel+spectrum  ,"gaus"+label+"_"+channel+spectrum  , rrv_x,rrv_mean1_gaus,rrv_sigma1_gaus)
		rrv_high   = RooRealVar("rrv_high"+label+"_"+channel+spectrum,"rrv_high"+label+"_"+channel+spectrum,0.3,0.0,0.99)
		model_pdf  = RooAddPdf("model_pdf"+label+"_"+channel+spectrum,"model_pdf"+label+"_"+channel+spectrum,erfExp,gaus,rrv_high)

	getattr(workspace,'import')(model_pdf)
	return workspace.pdf("model_pdf"+label+"_"+channel+spectrum)
			
def MakeExtendedModel(workspace, label, model,spectrum, channel, wtagger_label,constraint,peak="W"):
	
	print "Making extended model with name: " , "model"+label+"_"+channel+spectrum
	rrv_number = RooRealVar("rrv_number"+label+"_"+channel+spectrum,"rrv_number"+label+"_"+channel+spectrum,500.,0.,1e5)
	model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger_label,channel,constraint,peak)

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
			
def makeTTbarModel(workspace,label, model,channel, wtagger, constraint=[],peak="W", spectrum="_mj"):
	 
	 info =""
	 if label.find("_ttbar_data")!=-1 and label.find("fail")==-1:
		 rrv_number_total = RooRealVar("rrv_number_total_ttbar_data"+info+"_"+channel,"rrv_number_total_ttbar_data"+info+"_"+channel,500,0.,1e10)
		 eff_ttbar = RooRealVar("eff_ttbar_data"+info+"_"+channel,"eff_ttbar_data"+info+"_"+channel,0.7,0.2,1.0)
		 if peak == "Wt":eff_ttbar = RooRealVar("eff_ttbar_data"+info+"_"+channel,"eff_ttbar_data"+info+"_"+channel,0.7,0.0,1.0)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total))

	 elif label.find("_ttbar_data")!=-1 and label.find("fail")!=-1:
		 rrv_number_total = workspace.var("rrv_number_total_ttbar_data"+info+"_"+channel)
		 eff_ttbar        = workspace.var("eff_ttbar_data"+info+"_"+channel)
		 rrv_number       = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total))
		 
	 elif label.find("_ttbar_TotalMC")!=-1 and label.find("fail")==-1:
		 rrv_number_total = RooRealVar("rrv_number_total_ttbar_TotalMC"+info+"_"+channel,"rrv_number_total_ttbar_TotalMC"+info+"_"+channel,500,0.,1e10)
		 eff_ttbar = RooRealVar("eff_ttbar_TotalMC"+info+"_"+channel,"eff_ttbar_TotalMC"+info+"_"+channel,0.7,0.2,1.0)
		 if peak == "Wt": eff_ttbar = RooRealVar("eff_ttbar_TotalMC"+info+"_"+channel,"eff_ttbar_TotalMC"+info+"_"+channel,0.7,0.0,1.0)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "@0*@1", RooArgList(eff_ttbar,rrv_number_total))
		 
	 elif label.find("_ttbar_TotalMC")!=-1 and label.find("fail")!=-1:
		 rrv_number_total = workspace.var("rrv_number_total_ttbar_TotalMC"+info+"_"+channel)
		 eff_ttbar = workspace.var("eff_ttbar_TotalMC"+info+"_"+channel)
		 rrv_number = RooFormulaVar("rrv_number"+label+"_"+channel+spectrum, "(1-@0)*@1", RooArgList(eff_ttbar,rrv_number_total)) 
		 
	 model_pdf = MakeGeneralPdf(workspace,label,model,spectrum,wtagger,channel,constraint,peak)
	 model = RooExtendPdf("model"+label+"_"+channel+spectrum,"model"+label+"_"+channel+spectrum, model_pdf, rrv_number)
	 getattr(workspace,"import")(model)
	 return workspace.pdf("model"+label+"_"+channel+spectrum)
	