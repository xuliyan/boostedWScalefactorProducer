from ROOT import RooAbsPdf


def ErfExp( x,  c,  offset,  width):
	if width<1e-2: width = 1e-2
	if c==0: c=-1e-7
	return TMath.Exp(c*x)*(1.+TMath.Erf((x-offset)/width))/2.

def ErfExp( x,  x_min,  x_max,  c,  offset,  width):
	if width<1e-2: width = 1e-2
	if c==0: c=-1e-7
	minTerm = (TMath.Exp(c*c*width*width/4+c*offset) * TMath.Erf((2*x_min-c*width*width-2*offset)/2/width) - TMath.Exp(c*x_min) * TMath.Erf((x_min-offset)/width) - TMath.Exp(c*x_min))/-2/c
	maxTerm = (TMath.Exp(c*c*width*width/4+c*offset) * TMath.Erf((2*x_max-c*width*width-2*offset)/2/width) - TMath.Exp(c*x_max) * TMath.Erf((x_max-offset)/width) - TMath.Exp(c*x_max))/-2/c
	integral=(maxTerm-minTerm)
	return TMath.Exp(c*x)*(1.+TMath.Erf((x-offset)/width))/2./integral
	
def Exp( x,  c):
	return TMath.Exp(c*x)

def Exp( x,  x_min,  x_max,  c):
	integral = 0
	if c==0. : integral=x_max-x_min
	else: integral= ( TMath.Exp(c*x_max)-TMath.Exp(c*x_min) ) / c
	return TMath.Exp(c*x)/integral


def RooErfExpPdf(name, title, x, c, offset, width):
	width_tmp=width 
	if width<1e-2 : width_tmp=1e-2
	return ErfExp(x,c,offset,width_tmp)