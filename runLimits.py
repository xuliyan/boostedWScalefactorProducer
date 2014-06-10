#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen
from ROOT import gROOT, gStyle, gSystem, TLatex, TGaxis, TPaveText, TH2D, TColor, gPad

ROOT.gStyle.SetPadRightMargin(0.16);

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

#### basic allowed methods
parser.add_option('--makeCards',     action='store_true', dest='makeCards',         default=False, help='options to produce datacards and workspaces via fitting analysis')
parser.add_option('--computeLimits', action='store_true', dest='computeLimits',     default=False, help='basic option in order to compute asymptotic limits')
parser.add_option('--plotLimits',    action='store_true', dest='plotLimits',        default=False, help='basic option to plot asymptotic and pvalue')
parser.add_option('--computePvalue', action='store_true', dest='computePvalue',     default=False, help='basic option to plot asymptotic and pvalue')
parser.add_option('--biasStudy',     action='store_true', dest='biasStudy',         default=False, help='basic option to perform bias study with our own tool')
parser.add_option('--maximumLikelihoodFit', action='store_true', dest='maximumLikelihoodFit', default=False, help='basic option to run max Likelihood fit inside combination tool')
parser.add_option('--generateOnly',  action='store_true', dest='generateOnly', default=False, help='basic option to run the only generation with combiner')
parser.add_option('--makeLikelihoodScan', action='store_true', dest='makeLikelihoodScan', default=False, help='basic option to run the likelihood scan')

##### submit jobs to condor, lxbatch and hercules 
parser.add_option('--batchMode',      action='store_true', dest='batchMode',      default=False, help='to run jobs on condor fnal')
parser.add_option('--lxbatchCern',    action='store_true', dest='lxbatchCern',    default=False, help='run jobs on lxbatch system at cern, default is condor fnal')
parser.add_option('--herculesMilano', action='store_true', dest='herculesMilano', default=False, help='run jobs on hercules system at Milano, default is condor fnal')

##### other basci options for all the methods 
parser.add_option('--datacardDIR', action="store", type="string", dest="datacardDIR", default="")
parser.add_option('--channel',     action="store", type="string", dest="channel",     default="em")
parser.add_option('--queque',      action="store", type="string", dest="queque",      default="")
parser.add_option('--pseudodata',  action="store", type="int",    dest="pseudodata",  default=0)
parser.add_option('--systematics', action="store", type="int",    dest="systematics", default=1)
parser.add_option('--massPoint',   action="store", type="int",    dest="massPoint",   default=-1)
parser.add_option('--closuretest', action="store", type="int",    dest="closuretest", default=0)
parser.add_option('--cPrime',      action="store", type="int",    dest="cPrime",      default=-1)
parser.add_option('--brNew',       action="store", type="int",    dest="brNew",       default=-1)
parser.add_option('--odir',        action="store", type="string", dest="odir",        default=".")
parser.add_option('--sigChannel',  action="store", type="string", dest="sigChannel",  default="")
parser.add_option('--jetBin',      action="store", type="string", dest="jetBin",      default="")
parser.add_option('--skipJetSystematics',   action="store", type="int", dest="skipJetSystematics",    default=0)
parser.add_option('--turnOnAnalysis',       action="store", type="int",  dest="turnOnAnalysis",       default=0)
parser.add_option('--injectSingalStrenght', action="store", type=float,  dest="injectSingalStrenght", default=0., help='inject a singal in the toy generation')

parser.add_option('--higgsCombination', action="store", type="int", dest="higgsCombination", default=0)
parser.add_option('--interferenceModel', action="store", type="string", dest="interferenceModel", default="3")


###### options for Bias test in the combination tool
parser.add_option('--nToys',        action="store", type="int",    dest="nToys",       default=0)
parser.add_option('--crossedToys',  action="store", type="int",    dest="crossedToys", default=0)
parser.add_option('--inputGeneratedDataset', action="store", type="string",    dest="inputGeneratedDataset", default="")
parser.add_option('--outputTree',   action="store", type="int",    dest="outputTree",  default=0)


##### options specific for the bias tool 
parser.add_option('--shapetest',           action="store", type="int",    dest="shapetest",           default=0)
parser.add_option('--ttbarcontrolregion',  action="store", type="int",    dest="ttbarcontrolregion",  default=0)
parser.add_option('--mlvjregion',          action="store", type="string", dest="mlvjregion",          default="_sb_lo")
parser.add_option('--fitjetmass',          action="store", type="int",    dest="fitjetmass",          default=0)
parser.add_option('--onlybackgroundfit',   action="store", type="int",    dest="onlybackgroundfit",   default=0, help='run only background fit')
parser.add_option('--inflatejobstatistic', action="store", type="int",    dest="inflatejobstatistic", default=1, help='enlarge the generated statistics in the fit')
parser.add_option('--scalesignalwidth',    action="store", type="int",    dest="scalesignalwidth",    default=1, help='reduce the signal width by a factor x')

##### final plot options
parser.add_option('--makeSMLimitPlot',       action="store", type="int",    dest="makeSMLimitPlot",        default=0)
parser.add_option('--makeBSMLimitPlotMass',  action="store", type="int",    dest="makeBSMLimitPlotMass",   default=0)
parser.add_option('--makeBSMLimitPlotBRnew', action="store", type="int",    dest="makeBSMLimitPlotBRnew",  default=0)
parser.add_option('--makeBSMLimitPlot2D',    action="store", type="int",    dest="makeBSMLimitPlot2D",     default=0)
parser.add_option('--blindObservedLine',     action="store", type="int",    dest="blindObservedLine",      default=1)
parser.add_option('--plotPValue',            action="store", type="int",    dest="plotPValue",             default=0)
parser.add_option('--plotSignalStrenght',    action="store", type="int",    dest="plotSignalStrenght",     default=0)
parser.add_option('--plotLikelihoodScan',    action="store", type="int",    dest="plotLikelihoodScan",     default=0)

(options, args) = parser.parse_args()

#############################################
######### Get Some Global Variables #########
#############################################

mass  =  [ 600, 700, 800, 900,1000] ## define the mass point considered for the analysis
ccmlo =  [ 550, 600, 700, 750, 800]
ccmhi =  [ 700, 850, 950,1100,1200]
mjlo  =  [ 40, 40, 40, 40, 40]      ## min mJ cut
mjhi  =  [ 130, 130, 130, 130, 130] ## max mJ cut

################## options turnOn Analysis

if options.turnOnAnalysis :

 mlo      =  [ 400, 400, 400, 400, 400] ## min mlvj cut 
 mhi      =  [1500,1500,1500,1500,1500] ## max mlvj cut

else: 

 mlo      =  [ 550, 550, 550, 550, 550] ## min mlvj cut
 mhi      =  [1500,1500,1500,1500,1500] ## max mlvj cut


################## options for makeCards
 
if options.makeCards and options.turnOnAnalysis:
 
 shapeAlt =  ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"] ## basic shape
 shape    =  ["ErfPow2_v1", "ErfPow2_v1","ErfPow2_v1","ErfPow2_v1","ErfPow2_v1"] ## alternate one

elif not options.turnOnAnalysis and options.makeCards:

 shape      =  ["Exp","Exp","Exp","Exp","Exp"] ## basic shape
 shapeAlt   =  ["Pow","Pow","Pow","Pow","Pow"] ## alternate one

################## options for bias Study

if options.biasStudy:

 if not options.turnOnAnalysis and not options.fitjetmass:

#  shape_gen = ["Exp","Exp","Exp","Exp","Exp"]    
  shape_fit = ["Exp","Exp","Exp","Exp","Exp"]
  shape_gen = ["Pow2","Pow2","Pow2","Pow2","Pow2"]    
#  shape_fit = ["Pow2","Pow2","Pow2","Pow2","Pow2"]
#  shape_gen = ["Pow","Pow","Pow","Pow","Pow"]    
#  shape_fit = ["Pow","Pow","Pow","Pow","Pow"]

 elif options.turnOnAnalysis and not options.fitjetmass:

  shape_gen = ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"];    
  shape_fit = ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"];
#  shape_gen = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1"];    
#  shape_fit = ["ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1","ErfPowExp_v1"];
#  shape_gen = ["ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1"];    
#  shape_fit = ["ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1","ErfPow_v1"];

 elif options.fitjetmass:

#  shape_gen = ["ErfExp","ErfExp","ErfExp","ErfExp","ErfExp"];    
#  shape_fit = ["ErfExp","ErfExp","ErfExp","ErfExp","ErfExp"];
#  shape_gen = ["User1","User1","User1","User1","User1"];    
  shape_fit = ["User1","User1","User1","User1","User1"];
  shape_gen = ["ErfPow","ErfPow","ErfPow","ErfPow","ErfPow"];    
#  shape_fit = ["ErfPow","ErfPow","ErfPow","ErfPow","ErfPow"];

 isMC      = [0,0,0,0,0];

cprime = [01,02,03,05,07]
BRnew  = [00,01,02,03,04,05]
#cprime = [10]
#BRnew = [00]

########################################
###### Make Asymptotic Limit Plot ######
########################################

def getAsymLimits(file):
        
 f = ROOT.TFile(file);
 t = f.Get("limit"); 
 entries = t.GetEntries();
    
 lims = [0,0,0,0,0,0];

 limit_entries = 0;  
    
 for i in range(entries):
        
  t.GetEntry(i);
  t_quantileExpected = t.quantileExpected;
  t_limit = t.limit;
        
  #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;        
  if t_quantileExpected == -1.: lims[0] += t_limit; 
  elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] += t_limit;
  elif t_quantileExpected >= 0.15 and t_quantileExpected <= 0.17: lims[2] += t_limit;            
  elif t_quantileExpected == 0.5: lims[3] += t_limit; limit_entries += 1 ;            
  elif t_quantileExpected >= 0.83 and t_quantileExpected <= 0.85: lims[4] += t_limit;
  elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] += t_limit;
  else: print "Unknown quantile!"

 lims[0] = lims[0]/limit_entries ;
 lims[1] = lims[1]/limit_entries ;
 lims[2] = lims[2]/limit_entries ;
 lims[3] = lims[3]/limit_entries ;
 lims[4] = lims[4]/limit_entries ;
 lims[5] = lims[5]/limit_entries ;
    
 return lims;

########################################
###### Submit batch job for cards ######
########################################

def submitBatchJob( command, fn ):
    
 currentDir = os.getcwd();
    
 # create a dummy bash/csh
 outScript = open(fn+".sh","w");
 
 if not options.lxbatchCern and not options.herculesMilano :
  outScript.write('#!/bin/bash');
  outScript.write("\n"+'date');
  outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
  outScript.write("\n"+'echo "condor dir: " ${_CONDOR_SCRATCH_DIR}');    
  outScript.write("\n"+'cd '+currentDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'cd -');    
  outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
  outScript.write("\n"+'echo ${PATH}');
  outScript.write("\n"+'ls'); 
  outScript.write("\n"+command);  
  outScript.write("\n"+'tar -cvzf outputFrom_'+fn+'.tar.gz *');    
  outScript.close();

  condorScript = open("condor_"+fn,"w");
  condorScript.write('universe = vanilla')
  condorScript.write("\n"+"Executable = "+fn+".sh")
  condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
  condorScript.write("\n"+'Should_Transfer_Files = YES')
  condorScript.write("\n"+'Transfer_Input_Files = doFit_class_higgs.py, BiasStudy/do_fitBias_higgs.py, AutoDict_std__map_std__string_std__string__cxx.so')    
  condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
  condorScript.write("\n"+'Output = out_$(Cluster).stdout')
  condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
  condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
  condorScript.write("\n"+'Log    = out_$(Cluster).log')
  condorScript.write("\n"+'Notification    = Error')
  condorScript.write("\n"+'Queue 1')
  condorScript.close();

  os.system("condor_submit "+"condor_"+fn);

 elif options.lxbatchCern and not options.herculesMilano:
  outScript.write('#!/bin/bash');
  outScript.write("\n"+'cd '+currentDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'export PATH=${PATH}:'+currentDir);
  outScript.write("\n"+'echo ${PATH}');
  outScript.write("\n"+'ls');  
  outScript.write("\n"+command);
  outScript.write("\n"+'rm *.out');  
  outScript.close();
         
  os.system("chmod 777 "+currentDir+"/"+fn+".sh");
  if options.queque != "":
   os.system("bsub -q "+options.queque+" -cwd "+currentDir+" "+fn+".sh");
  else:
   os.system("bsub -q 8nh -cwd "+currentDir+" "+fn+".sh");
      
 elif not options.lxbatchCern and options.herculesMilano:

  outScript.write('#!/bin/bash');
  outScript.write("\n"+'cd '+currentDir);
  outScript.write("\n"+'eval `scram runtime -sh`');
  outScript.write("\n"+'cd -');
  outScript.write("\n"+'cp '+currentDir+'/BiasStudy/do_fitBias_higgs.py ./');
  outScript.write("\n"+'cp '+currentDir+'/doFit_class_higgs.py ./');
  outScript.write("\n"+'ls');  
  outScript.write("\n"+"unbuffer "+command+" > /gwteray/users/brianza/output"+fn+".txt");
  outScript.close();
         
  os.system("chmod 777 "+currentDir+"/"+fn+".sh");
  
  if options.queque != "":
   os.system("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh");
  else:
   os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
      

##########################################
###### Submit batch job for combine ######
##########################################

def submitBatchJobCombine( command, fn, mass, cprime, BRnew ):
    
    
    if options.sigChannel !="": 
     SIGCH = options.jetBin+"_"+options.sigChannel;
    else:
     SIGCH = options.jetBin;
    currentDir = os.getcwd();
    
    # create a dummy bash/csh
    outScript = open(fn+".sh","w");

    file1 = "hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass,SIGCH,cprime,BRnew);
    file2 = "hwwlvj_ggH%03d_mu%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);
    file3 = "hwwlvj_ggH%03d_el%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);

    if not options.lxbatchCern and not options.herculesMilano :   
     outScript.write('#!/bin/bash');
     outScript.write("\n"+'date');
     outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
     outScript.write("\n"+'cd '+currentDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\n"+'ls');    
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
      outScript.write("\n "+"mv higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
      
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();
    
     # link a condor script to your shell script
     condorScript = open("subCondor_"+fn,"w");
     condorScript.write('universe = vanilla')
     condorScript.write("\n"+"Executable = "+fn+".sh")
     condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
     condorScript.write("\n"+'Should_Transfer_Files = YES')
     condorScript.write("\n"+'Transfer_Input_Files = '+file1+', '+file2+', '+file3)    
     condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
     condorScript.write("\n"+'Output = out_$(Cluster).stdout')
     condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
     condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
     condorScript.write("\n"+'Log    = out_$(Cluster).log')
     condorScript.write("\n"+'Notification    = Error')
     condorScript.write("\n"+'Queue 1')
     condorScript.close();

     # submit the condor job     
     os.system("condor_submit "+"subCondor_"+fn)
 
    elif options.lxbatchCern and not options.herculesMilano: 

     outScript.write('#!/bin/bash');
     outScript.write("\n"+'date');
     outScript.write("\n"+'cd '+currentDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\n"+'ls');    
    
     file1 = "hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt"%(mass,options.channel,SIGCH,cprime,BRnew);
     file2 = "hwwlvj_ggH%03d_mu%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);
     file3 = "hwwlvj_ggH%03d_el%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);
     file4 = "hwwlvj_ggH%03d_em%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);

     outScript.write("\n"+'cp '+currentDir+'/'+file1+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file2+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file3+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file4+' ./');
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
      outScript.write("\n "+"mv higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();

     os.system("chmod 777 "+currentDir+"/"+fn+".sh");
     if options.queque!="" :
      os.system("bsub -q "+options.queque+" -cwd "+currentDir+" "+fn+".sh");
     else: 
      os.system("bsub -q 8nh -cwd "+currentDir+" "+fn+".sh");

    elif not options.lxbatchCern and options.herculesMilano: 

     outScript.write('#!/bin/bash');
     outScript.write("\n"+'date');
     outScript.write("\n"+'cd '+currentDir);
     outScript.write("\n"+'eval `scram runtime -sh`');
     outScript.write("\n"+'cd -');
     outScript.write("\n"+'ls');    
    
     file1 = "hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt"%(mass,options.channel,SIGCH,cprime,BRnew);
     file2 = "hwwlvj_ggH%03d_mu%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);
     file3 = "hwwlvj_ggH%03d_el%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);
     file4 = "hwwlvj_ggH%03d_em%s_%02d_%02d_workspace.root"%(mass,SIGCH,cprime,BRnew);

     outScript.write("\n"+'cp '+currentDir+'/'+file1+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file2+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file3+' ./');
     outScript.write("\n"+'cp '+currentDir+'/'+file4+' ./');
     outScript.write("\n"+command);
     if options.inputGeneratedDataset != "" and options.generateOnly:
      outScript.write("\n "+"mv higgsCombine* "+currentDir+"/"+options.inputGeneratedDataset);
     outScript.write("\n "+"rm rootstats* ");
     outScript.close();

     os.system("chmod 777 "+currentDir+"/"+fn+".sh");

     if options.queque != "":
      os.system("qsub -V -d "+currentDir+" -q "+options.queque+" "+currentDir+"/"+fn+".sh");
     else:
      os.system("qsub -V -d "+currentDir+" -q longcms "+currentDir+"/"+fn+".sh");
         

#####################################
###### definition of the style ######
#####################################

def setStyle():

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetPadBottomMargin(0.12);
  gStyle.SetPadLeftMargin(0.12);
  gStyle.SetCanvasColor(ROOT.kWhite);
  gStyle.SetCanvasDefH(600); #Height of canvas                                                                                                                                            
  gStyle.SetCanvasDefW(600); #Width of canvas                                                                                                                                             
  gStyle.SetCanvasDefX(0);   #POsition on screen
  gStyle.SetCanvasDefY(0);
  gStyle.SetPadTopMargin(0.05);
  gStyle.SetPadBottomMargin(0.15);#0.13);
  gStyle.SetPadLeftMargin(0.15);#0.16);
  gStyle.SetPadRightMargin(0.05);#0.02);                                                                                                                                                  
  # For the Pad:                                                                                                                                                                        
  gStyle.SetPadBorderMode(0);
  gStyle.SetPadColor(ROOT.kWhite);
  gStyle.SetPadGridX(ROOT.kFALSE);
  gStyle.SetPadGridY(ROOT.kFALSE);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);

  # For the frame:                                                                                                                                                                       
  gStyle.SetFrameBorderMode(0);
  gStyle.SetFrameBorderSize(1);
  gStyle.SetFrameFillColor(0);
  gStyle.SetFrameFillStyle(0);
  gStyle.SetFrameLineColor(1);
  gStyle.SetFrameLineStyle(1);
  gStyle.SetFrameLineWidth(1);

  gStyle.SetAxisColor(1, "XYZ");
  gStyle.SetStripDecimals(ROOT.kTRUE);
  gStyle.SetTickLength(0.03, "XYZ");
  gStyle.SetNdivisions(505, "XYZ");
  gStyle.SetPadTickX(1);  # To get tick marks on the opposite side of the frame                                                                                                        
  gStyle.SetPadTickY(1);
  gStyle.SetGridColor(0);
  gStyle.SetGridStyle(3);
  gStyle.SetGridWidth(1);

  gStyle.SetTitleColor(1, "XYZ");
  gStyle.SetTitleFont(42, "XYZ");
  gStyle.SetTitleSize(0.05, "XYZ");
  gStyle.SetTitleXOffset(1.15);#0.9);                                                                                                                                                  
  gStyle.SetTitleYOffset(1.3); # => 1.15 if exponents                                                                                                                                   
  gStyle.SetLabelColor(1, "XYZ");
  gStyle.SetLabelFont(42, "XYZ");
  gStyle.SetLabelOffset(0.007, "XYZ");
  gStyle.SetLabelSize(0.045, "XYZ");

  gStyle.SetPadBorderMode(0);
  gStyle.SetFrameBorderMode(0);
  gStyle.SetTitleTextColor(1);
  gStyle.SetTitleFillColor(10);
  gStyle.SetTitleFontSize(0.05);

  gStyle.SetOptStat(0);
  gStyle.SetOptTitle(0)
  gStyle.SetOptFit(1)

####################################################
### Get PValue from combine -M ProfileLikelihood ###
####################################################

def getPValueFromCard(file,observed):

 f = ROOT.TFile(file);
 t = f.Get("limit");
 entries = t.GetEntries();

 lims = 0;

 for i in range(entries):
  t.GetEntry(i);
  lims += t.limit ;

 return lims/entries;

##############################
#### Make SM Limits Plots ####  
##############################  

def makeSMLimitPlot(SIGCH,cprime = 10, brnew = 00):

    nPoints = len(mass);
    
    xbins     = array('f', []); xbins_env = array('f', []);
    ybins_exp = array('f', []); ybins_obs = array('f', []);
    ybins_1s  = array('f', []); ybins_2s  = array('f', []);

    setStyle();
     
    for i in range(len(mass)):
	curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
	curAsymLimits = getAsymLimits(curFile);
        xbins.append( mass[i] );
	xbins_env.append( mass[i] );
	ybins_exp.append( curAsymLimits[3] );
        ybins_obs.append( curAsymLimits[0] );
        ybins_2s.append( curAsymLimits[1] );
	ybins_1s.append( curAsymLimits[2] );

    for i in range( len(mass)-1, -1, -1 ):
	curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
	curAsymLimits = getAsymLimits(curFile);
        xbins_env.append( mass[i] );
        ybins_2s.append( curAsymLimits[5] );
	ybins_1s.append( curAsymLimits[4] );
                                                                    
    
    curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
    curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
    curGraph_1s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_1s);
    curGraph_2s  = ROOT.TGraphAsymmErrors(nPoints*2,xbins_env,ybins_2s);


    curGraph_obs.SetMarkerStyle(20);
    curGraph_obs.SetLineWidth(3);
    curGraph_obs.SetLineStyle(1);
    curGraph_obs.SetMarkerSize(1.6);
    curGraph_exp.SetMarkerSize(1.3);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_exp.SetLineStyle(2);
    curGraph_exp.SetLineWidth(3);
    curGraph_exp.SetMarkerSize(2);
    curGraph_exp.SetMarkerStyle(24);
    curGraph_exp.SetMarkerColor(ROOT.kBlack);

    curGraph_1s.SetFillColor(ROOT.kGreen);
    curGraph_1s.SetFillStyle(1001);
    curGraph_1s.SetLineStyle(ROOT.kDashed);
    curGraph_1s.SetLineWidth(3);

    curGraph_2s.SetFillColor(ROOT.kYellow);
    curGraph_2s.SetFillStyle(1001);
    curGraph_2s.SetLineStyle(ROOT.kDashed);
    curGraph_2s.SetLineWidth(3);
                               

    oneLine = ROOT.TF1("oneLine","1",599,1001);
    oneLine.SetLineColor(ROOT.kRed);
    oneLine.SetLineWidth(3);

    setStyle();
    
    can_SM = ROOT.TCanvas("can_SM","can_SM",600,650);
    hrl_SM = can_SM.DrawFrame(599,0.0,1001,ROOT.TMath.MaxElement(curGraph_2s.GetN(),curGraph_2s.GetY())*1.2);

    hrl_SM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_SM.GetYaxis().SetTitleOffset(1.35);
    hrl_SM.GetYaxis().SetTitleSize(0.045);
    hrl_SM.GetYaxis().SetTitleFont(42);

    hrl_SM.GetXaxis().SetTitle("M_{H} (GeV)");
    hrl_SM.GetXaxis().SetTitleSize(0.045);
    hrl_SM.GetXaxis().SetTitleFont(42);

    hrl_SM.GetYaxis().SetNdivisions(505);
    can_SM.SetGridx(1);
    can_SM.SetGridy(1);
                   
    curGraph_2s.Draw("F");
    curGraph_1s.Draw("Fsame");
    if not options.blindObservedLine : curGraph_obs.Draw("PCsame");
    curGraph_exp.Draw("Csame");

    leg2 = ROOT.TLegend(0.3,0.65,0.75,0.9);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);

    leg2.AddEntry(curGraph_exp,"Asympt. CL_{S} Expected","L")
    leg2.AddEntry(curGraph_1s, "Asympt. CL_{S} Expected #pm 1#sigma","LF")
    leg2.AddEntry(curGraph_2s, "Asympt. CL_{S} Expected #pm 2#sigma","LF")
                                       
    if not options.blindObservedLine:     leg2.AddEntry(curGraph_obs,"Asympt. CL_{S} Observed","LP")

    can_SM.Update();
    can_SM.RedrawAxis();
    can_SM.RedrawAxis("g");
    can_SM.Update();

    leg2.Draw();

    banner = TPaveText( 0.145, 0.953, 0.56, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.038);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.038);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    label_sqrt.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    
    can_SM.SaveAs("limitFigs/SMLim_%s_%02d_%02d.png"%(options.channel,cprime,brnew));
    can_SM.SaveAs("limitFigs/SMLim_%s_%02d_%02d.pdf"%(options.channel,cprime,brnew));

##############################
#### Make SM PValue Plots ####  
##############################  

def makeSMPValuePlot(SIGCH,cprime = 10,brnew = 00):

    nPoints = len(mass);
    
    xbins_obs     = array('f', []);
    xbins_exp     = array('f', []);

    ybins_obs     = array('f', []);
    ybins_exp     = array('f', []);

    setStyle();

    for i in range(len(mass)):
	curFile_obs = "higgsCombinehwwlvj_pval_obs_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
	curFile_exp = "higgsCombinehwwlvj_pval_exp_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
        xbins_obs.append(mass[i]); 
        xbins_exp.append(mass[i]); 
        print "higgsCombinehwwlvj_pval_obs_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
        ybins_obs.append(getPValueFromCard(curFile_obs,1));
        print "higgsCombinehwwlvj_pval_exp_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]); 
        ybins_exp.append(getPValueFromCard(curFile_exp,0));

    gr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins_obs,ybins_obs);
    gr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins_exp,ybins_exp);
                    
    gr_obs.SetLineColor(1); gr_obs.SetMarkerColor(1); gr_obs.SetMarkerStyle(20); gr_obs.SetLineWidth(3);gr_obs.SetMarkerSize(1.6);
    gr_exp.SetLineColor(2); gr_exp.SetMarkerColor(2); gr_exp.SetMarkerStyle(20); gr_exp.SetLineWidth(3);gr_exp.SetMarkerSize(1.6);
    gr_exp.SetLineStyle(8);

    oneSLine = ROOT.TF1("oneSLine","1.58655253931457074e-01",mass[0],mass[len(mass)-1]);
    oneSLine.SetLineColor(ROOT.kRed); oneSLine.SetLineWidth(2); oneSLine.SetLineStyle(2);
    twoSLine = ROOT.TF1("twoSLine","2.27501319481792155e-02",mass[0],mass[len(mass)-1]);
    twoSLine.SetLineColor(ROOT.kRed); twoSLine.SetLineWidth(2); twoSLine.SetLineStyle(2);
    threeSLine = ROOT.TF1("threeSLine","1.34989803163009588e-03",mass[0],mass[len(mass)-1]);
    threeSLine.SetLineColor(ROOT.kRed); threeSLine.SetLineWidth(2); threeSLine.SetLineStyle(2);
    fourSLine = ROOT.TF1("fourSLine","3.16712418331199785e-05",mass[0],mass[len(mass)-1]);
    fourSLine.SetLineColor(ROOT.kRed); fourSLine.SetLineWidth(2); fourSLine.SetLineStyle(2);
    
    banner = TLatex(0.32,0.955,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    ban1s = TLatex(950,1.58655253931457074e-01,("1 #sigma"));
    ban1s.SetTextSize(0.028); ban1s.SetTextColor(1)
    ban2s = TLatex(950,2.27501319481792155e-02,("2 #sigma"));
    ban2s.SetTextSize(0.028); ban2s.SetTextColor(1)
    ban3s = TLatex(950,1.34989803163009588e-03,("3 #sigma"));
    ban3s.SetTextSize(0.028); ban3s.SetTextColor(1);
    ban4s = TLatex(950,3.16712418331199785e-05,("4 #sigma"));
    ban4s.SetTextSize(0.028); ban4s.SetTextColor(1)

    leg2 = ROOT.TLegend(0.25,0.2,0.6,0.35);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    if options.jetBin == "":
     leg2.AddEntry( gr_obs, "obs signif, %s "%options.channel, "pl" );
     leg2.AddEntry( gr_exp, "exp signif, SM Higgs", "pl" );
    else:
     leg2.AddEntry( gr_obs, "obs signif, VBF e+mu", "pl" );
     leg2.AddEntry( gr_exp, "exp signif, SM Higgs", "pl" );

    zeroLine = ROOT.TF1("zeroLine","0.5",mass[0],mass[len(mass)-1]);
    zeroLine.SetLineColor(ROOT.kBlue); oneSLine.SetLineWidth(2);

    can = ROOT.TCanvas("can","can",600,650);
    hrl = can.DrawFrame(mass[0],1e-6,mass[len(mass)-1],1.);
    hrl.GetYaxis().SetTitle("p-value");
    hrl.GetXaxis().SetTitle("m_{H} (GeV)");
    can.SetGrid();
    ROOT.gPad.SetLogy();
    gr_obs.Draw("PL");
    gr_exp.Draw("PLsame");
    zeroLine.Draw("same");
    oneSLine.Draw("same");
    twoSLine.Draw("same");
    threeSLine.Draw("same");
    fourSLine.Draw("same");
    banner.Draw();
    leg2.Draw();
    ban1s.Draw();
    ban2s.Draw();
    ban3s.Draw();
    ban4s.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can.SaveAs("limitFigs/SMPvals_%s_%02d_%02d.pdf"%(options.channel,cprime,brnew),"pdf");
    can.SaveAs("limitFigs/SMPvals_%s_%02d_%02d.png"%(options.channel,cprime,brnew),"png");

####################################
#### Make Signal Strenght Plots ####  
####################################  

def makeSignalStrenghtPlot(SIGCH,cprime = 10,brnew = 00):

    nPoints = len(mass);

    xbins_mu         = array('f', []);
    xbins_mu_err_up  = array('f', []);
    xbins_mu_err_dn  = array('f', []);
    ybins_mu         = array('f', []);
    ybins_mu_err_up  = array('f', []);
    ybins_mu_err_dn  = array('f', []);
    ybins_mu_err_up_2s  = array('f', []);
    ybins_mu_err_dn_2s  = array('f', []);
    ybins_mu_ratio   = array('f', []);
    xbins_obs     = array('f', [0.]);
    ybins_obs     = array('f', [0.]);

    setStyle();
        
    for i in range(len(mass)):

     curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.MaxLikelihoodFit.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
 
     f = ROOT.TFile(curFile);
     t = f.Get("limit");
     entries = t.GetEntries();

     mu = 0;
     effect_entry = 0 ;
     mu_err_up = 0 ; mu_err_dn = 0 ;
     mu_err_up_2s = 0 ; mu_err_dn_2s = 0 ;
     
     mu_temp = 0 ;  
     for ientry in range(t.GetEntries()):
       t.GetEntry(ientry);
       if t.quantileExpected == 0.5 : mu += t.limit ;  effect_entry += 1; mu_temp = t.limit ; continue ; 
       if t.quantileExpected > 0.15 and t.quantileExpected < 0.17 : mu_err_dn += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected > 0.83 and t.quantileExpected < 0.85 : mu_err_up += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected >= 0.024 and t.quantileExpected <= 0.026: mu_err_up_2s += abs(t.limit-mu_temp) ; continue ;
       if t.quantileExpected >= 0.974 and t.quantileExpected <= 0.976: mu_err_dn_2s += abs(t.limit-mu_temp) ; continue ;
                         
     mu = mu/effect_entry;
     mu_err_dn = mu_err_dn/effect_entry;
     mu_err_up = mu_err_up/effect_entry;
     mu_err_dn_2s = mu_err_dn_2s/effect_entry;
     mu_err_up_2s = mu_err_up_2s/effect_entry;
     
     xbins_mu.append(mass[i]); 
     xbins_mu_err_up.append(0.); 
     xbins_mu_err_dn.append(0.); 

     ybins_mu.append(mu);
     ybins_mu_err_up.append(mu_err_up);
     ybins_mu_err_dn.append(mu_err_dn);
     ybins_mu_err_up_2s.append(mu_err_up_2s);
     ybins_mu_err_dn_2s.append(mu_err_dn_2s);

     if mu >= 0 : ybins_mu_ratio.append(mu/mu_err_dn);
     else:        ybins_mu_ratio.append(mu/mu_err_up);

     if options.plotPValue == 1:
	curFile_obs = "higgsCombinehwwlvj_pval_obs_ggH%03d_%s%s_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
        xbins_obs.append(mass[i]); 
        ybins_obs.append(getPValueFromCard(curFile_obs,1));

                    
    gr_mu_1 = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu);
    gr_mu_1.SetLineColor(1); gr_mu_1.SetMarkerColor(1); gr_mu_1.SetMarkerStyle(20); gr_mu_1.SetLineWidth(3);gr_mu_1.SetMarkerSize(1.6);

    gr_mu_ratio = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu_ratio);
    gr_mu_ratio.SetLineColor(1); gr_mu_ratio.SetMarkerColor(1); gr_mu_ratio.SetMarkerStyle(20); gr_mu_ratio.SetLineWidth(3);gr_mu_ratio.SetMarkerSize(1.6);

    gr_mu_1_1s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn,ybins_mu_err_up);
    gr_mu_1_1s.SetLineColor(1); gr_mu_1_1s.SetMarkerColor(1); gr_mu_1_1s.SetMarkerStyle(20); gr_mu_1_1s.SetLineWidth(5); gr_mu_1_1s.SetMarkerSize(1.6);

    gr_mu_1_2s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,ybins_mu,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn_2s,ybins_mu_err_up_2s);
    gr_mu_1_2s.SetLineColor(ROOT.kBlue); gr_mu_1_2s.SetMarkerColor(ROOT.kBlue); gr_mu_1_2s.SetMarkerStyle(20); gr_mu_1_2s.SetLineWidth(3); gr_mu_1_2s.SetMarkerSize(1.6);

    gr_mu_zero_1s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,xbins_mu_err_up,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn,ybins_mu_err_up);
    gr_mu_zero_2s = ROOT.TGraphAsymmErrors(nPoints,xbins_mu,xbins_mu_err_up,xbins_mu_err_dn,xbins_mu_err_up,ybins_mu_err_dn_2s,ybins_mu_err_up_2s);
    
    gr_mu_zero_1s.SetFillStyle(3001);
    gr_mu_zero_1s.SetFillColor(ROOT.kRed);
    gr_mu_zero_2s.SetFillStyle(3002);
    gr_mu_zero_2s.SetFillColor(210);

    gr_obs = ROOT.TGraphAsymmErrors(nPoints+1,xbins_obs,ybins_obs);                    
    gr_obs.SetLineColor(ROOT.kBlue); gr_obs.SetMarkerColor(ROOT.kBlue); gr_obs.SetMarkerStyle(28); gr_obs.SetLineWidth(3);gr_obs.SetMarkerSize(1.6);

    oneSLine0 = ROOT.TF1("oneSLine","0",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine0.SetLineColor(ROOT.kRed); oneSLine0.SetLineWidth(2);
    oneSLine1 = ROOT.TF1("oneSLine","1",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine1.SetLineColor(ROOT.kRed); oneSLine1.SetLineWidth(2); oneSLine1.SetLineStyle(2);
    oneSLine2 = ROOT.TF1("oneSLine","2",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine2.SetLineColor(ROOT.kRed); oneSLine2.SetLineWidth(2); oneSLine2.SetLineStyle(2);
    oneSLine3 = ROOT.TF1("oneSLine","3",mass[0]-50,mass[len(mass)-1]+50);
    oneSLine3.SetLineColor(ROOT.kRed); oneSLine3.SetLineWidth(2); oneSLine3.SetLineStyle(2);
    
    banner = TLatex(0.32,0.955,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    ban1s = TLatex(950,1.,("#mu SM"));
    ban1s.SetTextSize(0.028); ban1s.SetTextColor(1)

    leg2 = ROOT.TLegend(0.25,0.2,0.6,0.35);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    if options.jetBin == "":
     leg2.AddEntry( gr_mu_1, "signal strenght %s "%options.channel, "pl" );
    else:
     leg2.AddEntry( gr_mu_1, "signal strenght, VBF e+#mu", "pl" );

    can = ROOT.TCanvas("can","can",600,650);
    hrl = can.DrawFrame(mass[0],-5,mass[len(mass)-1],5);
    hrl.GetYaxis().SetTitle("signal strenght");
    hrl.GetXaxis().SetTitle("m_{H} (GeV)");
    can.SetGrid();
    gr_mu_1.Draw("P");
    gr_mu_zero_2s.Draw("e3same");
    gr_mu_zero_1s.Draw("e3same"); 
    gr_mu_1.Draw("Psame");
    oneSLine0.Draw("same");
   
    banner.Draw();
    leg2.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can.SaveAs("limitFigs/signal_strenght_band_%s%s_%02d_%02d.pdf"%(options.channel,SIGCH,cprime,brnew),"pdf");
    can.SaveAs("limitFigs/signal_strenght_band_%s%s_%02d_%02d.png"%(options.channel,SIGCH,cprime,brnew),"png");

    can_2 = ROOT.TCanvas("can_2","can_2",600,650);
    hrl_2 = can_2.DrawFrame(mass[0]-50,-5,mass[len(mass)-1]+50,5);
    hrl_2.GetYaxis().SetTitle("signal strenght");
    hrl_2.GetXaxis().SetTitle("m_{H} (GeV)");
    can_2.SetGrid();
    gr_mu_1_2s.Draw("P");
    gr_mu_1_1s.Draw("P"); 
    oneSLine0.Draw("same");
    oneSLine1.Draw("same");
    oneSLine2.Draw("same");
    oneSLine3.Draw("same");
   
    banner.Draw();
    leg2.Draw();

    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can_2.SaveAs("limitFigs/signal_strenght_%s%s_%02d_%02d.pdf"%(options.channel,SIGCH,cprime,brnew),"pdf");
    can_2.SaveAs("limitFigs/signal_strenght_%s%s_%02d_%02d.png"%(options.channel,SIGCH,cprime,brnew),"png");


    can_ratio = ROOT.TCanvas("can_ratio","can_ratio",600,650);
    hrl_ratio = can_ratio.DrawFrame(mass[0]-50,-5,mass[len(mass)-1]+50,5);
    hrl_ratio.GetYaxis().SetTitle("#mu/#sigma_{#mu}");
    hrl_ratio.GetXaxis().SetTitle("m_{H} (GeV)");
    can_ratio.SetGrid();
    leg2.Clear()
    if options.jetBin == "":
     leg2.AddEntry( gr_mu_ratio, "#sigma_{#mu}/#mu %s "%options.channel, "pl" );
    else:
     leg2.AddEntry( gr_mu_ratio, "#sigma_{#mu}/#mu, VBF e+#mu", "pl" );
    gr_mu_ratio.Draw("P");
    if options.plotPValue == 1 : gr_obs.Draw("P"); leg2.AddEntry( gr_obs, "local significance", "pl" );
    
    oneSLine0.Draw("same");
    oneSLine1.Draw("same");
    oneSLine2.Draw("same");
    oneSLine3.Draw("same");
    banner.Draw();
    leg2.Draw();
    
    os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
    can_ratio.SaveAs("limitFigs/signal_strenght_ratio_%s%s_%02d_%02d.pdf"%(options.channel,SIGCH,cprime,brnew),"pdf");
    can_ratio.SaveAs("limitFigs/signal_strenght_ratio_%s%s_%02d_%02d.png"%(options.channel,SIGCH,cprime,brnew),"png");

#################################
### make likelihood scan plot ###
#################################    

def makeLikelihoodScanPlot(SIGCH,cprime,brnew):

    nPoints = len(mass);
    
    setStyle();

    banner = TLatex(0.32,0.955,("CMS Preliminary, 19.3 fb^{-1} at #sqrt{s}=8TeV"));
    banner.SetNDC(); banner.SetTextSize(0.035);

    leg2 = ROOT.TLegend(0.55,0.65,0.8,0.85);
    leg2.SetFillStyle(0);
    leg2.SetBorderSize(1);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    
    for i in range(len(mass)):

	curFile = "higgsCombinehwwlvj_LikelihoodScan_ggH%03d_%s%s_%02d_%02d_unbin.MultiDimFit.mH%02d.root"%(mass[i],options.channel,SIGCH,cprime,brnew,mass[i]);
        f = ROOT.TFile(curFile);
        t = f.Get("limit");
        
        xbins_mu = array('f',[]); 
        ybins_mu = array('f',[]); 

        for ientry in range(t.GetEntries()):
         if ientry == 0 : continue ;   
         t.GetEntry(ientry);   
         xbins_mu.append(t.r); 
         ybins_mu.append(2*t.deltaNLL); 

        gr_mu = ROOT.TGraph(t.GetEntries()-1,xbins_mu,ybins_mu);
        gr_mu.Sort();
                                    
        gr_mu.SetLineWidth(2);
        gr_mu.SetLineColor(2);
        gr_mu.SetMarkerStyle(20);

        if options.jetBin == "":
         leg2.AddEntry(gr_mu, "#mu Likelihood scan %s "%options.channel, "pl" );
        else:
         leg2.AddEntry(gr_mu, "#mu Likelihood scan, VBF e+#mu", "pl" );

        can = ROOT.TCanvas("can_%d"%(mass[i]),"can_%d"%(mass[i]),600,650);
        gr_mu.GetYaxis().SetTitle("-2#Delta ln(L)");
        gr_mu.GetXaxis().SetTitle("signal strenght");
        can.SetGrid();
        gr_mu.Draw("AC");
        banner.Draw();
        leg2.Draw();

        os.system("mkdir -p %s/limitFigs/"%(os.getcwd()));
        can.SaveAs("limitFigs/LikelihoodScan_%03d_%s%s_%02d_%02d.pdf"%(mass[i],options.channel,SIGCH,cprime,brnew),"pdf");
        can.SaveAs("limitFigs/LikelihoodScan_%03d_%s%s_%02d_%02d.png"%(mass[i],options.channel,SIGCH,cprime,brnew),"png");
        leg2.Clear();
    
                                           
##############################################
### Make the BSM Limit plot vs mass and c' ###
##############################################    
def makeBSMLimitPlotMass(SIGCH):

    print "module ===> makeBSMLimits_vsMass";
    curcolors = [1,2,210,4,6,12,7];
    brnew = 00; nPoints = len(mass);

    massCS  = [];
    if SIGCH == "" or SIGCH == "_2jet":
        massCS.append((0.5230+0.09688));
        massCS.append((0.2290+0.06330));
        massCS.append((0.1097+0.04365));
        massCS.append((0.0571+0.03164));
        massCS.append((0.0320+0.02399));
    elif SIGCH == "_ggH":
        massCS.append((0.5230));
        massCS.append((0.2290));
        massCS.append((0.1097));
        massCS.append((0.0571));
        massCS.append((0.0320));
    elif SIGCH == "_vbfH":
        massCS.append((0.09688));
        massCS.append((0.06330));
        massCS.append((0.04365));
        massCS.append((0.03164));
        massCS.append((0.02399));
    else:
        print "problem!"
        
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];

    gridMax    = -999;
    gridMaxSig = -999;

    tGraphs_exp = [];
    tGraphs_obs = [];
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];
    tGraphs_csXbr_th = [];
    
    for j in range(len(cprime)):

        xbins           = array('f', []); ybins_exp       = array('f', []);
        ybins_obs       = array('f', []); ybins_csXbr_exp = array('f', []);
        ybins_csXbr_obs = array('f', []); ybins_csXbr_th  = array('f', []);

        for i in range(len(mass)):
         if (cprime[j]*0.1)<=(1-brnew*0.1):
            curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[i],options.channel,SIGCH,cprime[j],brnew,mass[i]);

            curAsymLimits = getAsymLimits(curFile);
            xbins.append( mass[i] );
            ybins_exp.append( curAsymLimits[3] );
            ybins_obs.append( curAsymLimits[0] );
            ybins_csXbr_exp.append(curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            ybins_csXbr_obs.append(curAsymLimits[0]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            ybins_csXbr_th.append(1.*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );

            if gridMax < curAsymLimits[3]:
                gridMax = curAsymLimits[3];

            cscur = ( curAsymLimits[3]*massCS[i]*cprime[j]*0.1*(1-brnew*0.1)*massBRWW[i] );
            if gridMaxSig < cscur:
                gridMaxSig = cscur;

        curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_obs);
        curGraph_csXbr_th  = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_th);
        curGraph_csXbr_exp.SetLineStyle(10);
        curGraph_csXbr_exp.SetLineWidth(2);
        curGraph_csXbr_obs.SetLineWidth(2);
        curGraph_csXbr_exp.SetMarkerSize(2);
        curGraph_csXbr_obs.SetMarkerSize(2);
        curGraph_csXbr_th.SetLineWidth(2);

        tGraphs_exp.append( curGraph_exp );
        tGraphs_obs.append( curGraph_obs );
        tGraphs_csXbr_exp.append( curGraph_csXbr_exp );
        tGraphs_csXbr_obs.append( curGraph_csXbr_obs );
        tGraphs_csXbr_th.append( curGraph_csXbr_th );


    setStyle();

    can_BSM = ROOT.TCanvas("can_BSM","can_BSM",630,600);
    hrl_BSM = can_BSM.DrawFrame(599,0.0,1001,gridMax*1.5);

    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.35);
    hrl_BSM.GetYaxis().SetTitleSize(0.045);
    hrl_BSM.GetYaxis().SetTitleFont(42);

    hrl_BSM.GetXaxis().SetTitle("M_{H} (GeV)");
    hrl_BSM.GetXaxis().SetTitleSize(0.045);
    hrl_BSM.GetXaxis().SetTitleFont(42);

    hrl_BSM.GetYaxis().SetNdivisions(505);
    can_BSM.SetGridx(1);
    can_BSM.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(cprime)):
        tGraphs_exp[k].SetLineStyle(1);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetMarkerColor(curcolors[k]);
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);
        tGraphs_exp[k].Draw("PL");
        if not options.blindObservedLine: tGraphs_obs[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f"%(float((cprime[k])/10.));
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L");
        tmplabel = "obs., C'^{ 2} = %1.1f"%(float((cprime[k])/10.));
        if not options.blindObservedLine: leg2.AddEntry(tGraphs_obs[k],tmplabel,"L");


    banner2 = TLatex(0.17,0.91,("BR_{new} = 0"));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can_BSM.Update();
    can_BSM.RedrawAxis();
    can_BSM.RedrawAxis("g");
    can_BSM.Update();
    
    banner = TPaveText( 0.145, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.036);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.036);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    label_sqrt.Draw();

    leg2.Draw();
    banner2.Draw();
    can_BSM.SaveAs("limitFigs/BSMLim%s_%s.png"%(SIGCH,options.channel));
    can_BSM.SaveAs("limitFigs/BSMLim%s_%s.pdf"%(SIGCH,options.channel));

    can_BSMsig = ROOT.TCanvas("can_BSMsig","can_BSMsig",630,600);
    hrl_BSMsig = can_BSMsig.DrawFrame(599,0.001,1001,gridMaxSig*15);

    hrl_BSMsig.GetYaxis().SetTitle("#sigma #times BR_{WW}");
    hrl_BSMsig.GetYaxis().SetTitleOffset(1.35);
    hrl_BSMsig.GetYaxis().SetTitleSize(0.045);
    hrl_BSMsig.GetYaxis().SetTitleFont(42);

    hrl_BSMsig.GetXaxis().SetTitle("M_{H} (GeV)");
    hrl_BSMsig.GetXaxis().SetTitleSize(0.045);
    hrl_BSMsig.GetXaxis().SetTitleFont(42);

    hrl_BSMsig.GetYaxis().SetNdivisions(505);
    can_BSMsig.SetGridx(1);
    can_BSMsig.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.70,0.75,0.88);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    can_BSMsig.Update();
    can_BSMsig.RedrawAxis();
    can_BSMsig.RedrawAxis("g");
    can_BSMsig.Update();
    ROOT.gPad.SetLogy();
    
    for k in range(len(cprime)):
        tGraphs_csXbr_exp[k].SetLineStyle(1);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_th[k].SetLineStyle(2);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);
        tGraphs_csXbr_exp[k].Draw("PL");
        if not options.blindObservedLine : tGraphs_csXbr_obs[k].Draw("PL");
        tGraphs_csXbr_th[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")
        tmplabel = "obs., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
        if not options.blindObservedLine: leg2.AddEntry(tGraphs_csXbr_obs[k],tmplabel,"L");
        tmplabel = "th., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) )
        leg2.AddEntry(tGraphs_csXbr_th[k],tmplabel,"L");

    leg2.Draw();
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();
    
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_%s_Sigma.png"%(SIGCH,options.channel));
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_%s_Sigma.pdf"%(SIGCH,options.channel));


###############################################
### Make the BSM Limit plot vs c' and brNew ###
###############################################
    
def makeBSMLimitPlotBRnew(SIGCH,mass):

    print "module ===> makeBSMLimits_vsBRnew";

    curcolors = [1,2,210,4,6,12,7];
    BRnew_x = [0,0.1,0.2,0.3,0.4,0.5];
    massindex = {600:0,700:1,800:2,900:3,1000:4}
    nPoints = len(BRnew);
    massXS  = [];

    if SIGCH == "" or SIGCH == "_2jet":
	massXS.append((0.5230+0.09688));
        massXS.append((0.2288+0.06330));
        massXS.append((0.1095+0.04365));
        massXS.append((0.05684+0.03164));
        massXS.append((0.03163+0.02399));
    elif SIGCH == "_ggH":
        massXS.append((0.5230));
        massXS.append((0.2288));
        massXS.append((0.1095));
        massXS.append((0.05684));
        massXS.append((0.03163));
    elif SIGCH == "_vbfH":
        massXS.append((0.09688));
        massXS.append((0.06330));
        massXS.append((0.04365));
        massXS.append((0.03164));
        massXS.append((0.02399));
    else:
        print "problem!"
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];

    gridMax = -999;
    gridMaxSig = -999;

    tGraphs_exp = [];
    tGraphs_obs = [];
    tGraphs_csXbr_exp = [];
    tGraphs_csXbr_obs = [];
    tGraphs_csXbr_th = [];

    for j in range(len(cprime)):

        xbins           = array('f', []);
        ybins_exp       = array('f', []);
        ybins_obs       = array('f', []);
        ybins_csXbr_exp = array('f', []);
        ybins_csXbr_obs = array('f', []);
        ybins_csXbr_th  = array('f', []);

        for i in range(len(BRnew)):
         if (cprime[j]*0.1)<=(1-BRnew[i]*0.1):
            curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass,options.channel,SIGCH,cprime[j],BRnew[i],mass);
            curAsymLimits = getAsymLimits(curFile);
            xbins.append(BRnew_x[i]);
            ybins_exp.append( curAsymLimits[3] );
            ybins_obs.append( curAsymLimits[0] );
            ybins_csXbr_exp.append( curAsymLimits[3]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            ybins_csXbr_obs.append( curAsymLimits[0]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            ybins_csXbr_th.append( 1.*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );

            if gridMax < curAsymLimits[3]: gridMax = curAsymLimits[3];
            cscur = ( curAsymLimits[3]*massXS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]] );
            if gridMaxSig < cscur: gridMaxSig = cscur;

        curGraph_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_exp);
        curGraph_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_obs);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_exp);
        curGraph_csXbr_obs = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_obs);
        curGraph_csXbr_th  = ROOT.TGraphAsymmErrors(nPoints,xbins,ybins_csXbr_th);
        curGraph_obs.SetMarkerStyle(20);
        curGraph_obs.SetLineWidth(3);
        curGraph_obs.SetLineStyle(1);
        curGraph_obs.SetMarkerSize(1);
        curGraph_exp.SetMarkerSize(1.3);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);
        curGraph_exp.SetLineStyle(2);
        curGraph_exp.SetLineWidth(3);
        curGraph_exp.SetMarkerSize(1);
        curGraph_exp.SetMarkerStyle(24);
        curGraph_exp.SetMarkerColor(ROOT.kBlack);

        curGraph_csXbr_exp.SetMarkerSize(1);
        curGraph_csXbr_exp.SetMarkerStyle(24);
        curGraph_csXbr_exp.SetMarkerColor(ROOT.kBlack);

        tGraphs_exp.append(curGraph_exp);
        tGraphs_obs.append(curGraph_obs);
        tGraphs_csXbr_exp.append(curGraph_csXbr_exp);
        tGraphs_csXbr_obs.append(curGraph_csXbr_obs);
        tGraphs_csXbr_th.append(curGraph_csXbr_th);

    setStyle();
    can_BSM = ROOT.TCanvas("can_BSM_BR","can_BSM",630,600);
    hrl_BSM = can_BSM.DrawFrame(0.0,0.0,0.5,gridMax*1.5);

    hrl_BSM.GetYaxis().SetTitle("#mu = #sigma_{95%} / #sigma_{SM}");
    hrl_BSM.GetYaxis().SetTitleOffset(1.35);
    hrl_BSM.GetYaxis().SetTitleSize(0.045);
    hrl_BSM.GetYaxis().SetTitleFont(42);

    hrl_BSM.GetXaxis().SetTitle("BR_{new}");
    hrl_BSM.GetXaxis().SetTitleSize(0.045);
    hrl_BSM.GetXaxis().SetTitleFont(42);

    hrl_BSM.GetYaxis().SetNdivisions(505);
    can_BSM.SetGridx(1);
    can_BSM.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(cprime)): 
        tGraphs_exp[k].SetLineStyle(1);
        tGraphs_exp[k].SetLineColor(curcolors[k]);
        tGraphs_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_obs[k].SetLineColor(curcolors[k]);
        tGraphs_obs[k].SetMarkerColor(curcolors[k]);
        tGraphs_exp[k].SetLineWidth(2);
        tGraphs_obs[k].SetLineWidth(2);
        tGraphs_exp[k].Draw("PL");
        if not options.blindObservedLine: tGraphs_obs[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f    "%( float((cprime[k])/10.) )
        leg2.AddEntry(tGraphs_exp[k],tmplabel,"L")


    can_BSM.Update();
    can_BSM.RedrawAxis();
    can_BSM.RedrawAxis("g");
    can_BSM.Update();

    banner = TPaveText( 0.145, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.038);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.5,0.953,0.96,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.038);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");

    banner2 = TLatex(0.17,0.91,("Higgs Mass, %i GeV"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    leg2.Draw();
    banner.Draw();
    banner2.Draw();
    label_sqrt.Draw();

    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.png"%(SIGCH,mass));
    can_BSM.SaveAs("limitFigs/BSMLim%s_Mu_vsBRnew_%i.pdf"%(SIGCH,mass));

    can_BSMsig = ROOT.TCanvas("can_BSM_Sig","can_BSM",630,600);
    hrl_BSMSig = can_BSMsig.DrawFrame(0.,0.0,0.5,gridMaxSig*1.5);

    hrl_BSMSig.GetYaxis().SetTitle("#sigma #times BR_{WW}");
    hrl_BSMSig.GetYaxis().SetTitleOffset(1.35);
    hrl_BSMSig.GetYaxis().SetTitleSize(0.045);
    hrl_BSMSig.GetYaxis().SetTitleFont(42);

    hrl_BSMSig.GetXaxis().SetTitle("BR_{new}");
    hrl_BSMSig.GetXaxis().SetTitleSize(0.045);
    hrl_BSMSig.GetXaxis().SetTitleFont(42);

    hrl_BSMSig.GetYaxis().SetNdivisions(505);
    can_BSMsig.SetGridx(1);
    can_BSMsig.SetGridy(1);

    leg2 = ROOT.TLegend(0.25,0.65,0.75,0.85);
    leg2.SetFillColor(0);
    leg2.SetShadowColor(0);
    leg2.SetTextFont(42);
    leg2.SetTextSize(0.028);
    leg2.SetNColumns(2);

    for k in range(len(cprime)):
        tGraphs_csXbr_exp[k].SetLineStyle(1);
        tGraphs_csXbr_exp[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_obs[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_th[k].SetLineStyle(2);
        tGraphs_csXbr_th[k].SetLineColor(curcolors[k]);
        tGraphs_csXbr_exp[k].SetLineWidth(2);
        tGraphs_csXbr_obs[k].SetLineWidth(2);
        tGraphs_csXbr_exp[k].SetMarkerColor(curcolors[k]);
        tGraphs_csXbr_exp[k].Draw("PL");
        if not options.blindObservedLine: tGraphs_csXbr_obs[k].Draw("PL");

        tmplabel = "exp., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) );
        leg2.AddEntry(tGraphs_csXbr_exp[k],tmplabel,"L")
        tmplabel = "obs., C'^{ 2} = %1.1f"%( float((cprime[k])/10.) );

    leg2.Draw();
    banner.Draw();
    banner2.Draw();
    label_sqrt.Draw();
     
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.png"%(SIGCH,mass));
    can_BSMsig.SaveAs("limitFigs/BSMLim%s_Sigma_vsBRnew_%i.pdf"%(SIGCH,mass));

#############################################
### Make the BSM 2D Scane vs c' and brNew ###
#############################################

def makeBSMLimitPlot2D( SIGCH, mass ):


    stylePath = os.getenv("ROOTStyle");
#    gROOT.ProcessLine(".x "+stylePath+"/rootPalette.C");
    
    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];            
    print "module ===> makeBSMLimits_2D";

    massindex = {600:0,700:1,800:2,900:3,1000:4}
    mass_XS  = [];
    BRnew_y  = array('d',[0.,0.15,0.25,0.35,0.45,0.55]);
    cprime_x = array('d',[0.,0.15,0.25,0.35,0.55,0.75]);

    if SIGCH == "" or SIGCH == "_2jet":
        mass_XS.append((0.5230+0.09688));
        mass_XS.append((0.2288+0.06330));
        mass_XS.append((0.1095+0.04365));
        mass_XS.append((0.05684+0.03164));
        mass_XS.append((0.03163+0.02399));
    elif SIGCH == "_ggH":
        mass_XS.append((0.5230));
        mass_XS.append((0.2288));
        mass_XS.append((0.1095));
        mass_XS.append((0.05684));
        mass_XS.append((0.03163));
    elif SIGCH == "_vbfH":
        mass_XS.append((0.09688));
        mass_XS.append((0.06330));
        mass_XS.append((0.04365));
        mass_XS.append((0.03164));
        mass_XS.append((0.02399));
    else:
        print "problem!"

    setStyle();

    h2d_exp = TH2D("h2d_exp","",len(cprime_x)-1,cprime_x,len(BRnew_y)-1,BRnew_y);
    h2d_obs = TH2D("h2d_obs","",len(cprime_x)-1,cprime_x,len(BRnew_y)-1,BRnew_y);
    h2d_csXbr_exp = TH2D("h2d_csXbr_exp","",len(cprime_x)-1,cprime_x,len(BRnew_y)-1,BRnew_y);
    h2d_csXbr_obs = TH2D("h2d_csXbr_obs","",len(cprime_x)-1,cprime_x,len(BRnew_y)-1,BRnew_y);

    for j in range(len(cprime)):
        for i in range(len(BRnew)):
         if (cprime[j]*0.1)<=(1-BRnew[i]*0.1):
            curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass,options.channel,SIGCH,cprime[j],BRnew[i],mass); 
            curAsymLimits = getAsymLimits(curFile);
            h2d_exp.SetBinContent(j+1,i+1,curAsymLimits[3]);
            h2d_obs.SetBinContent(j+1,i+1,curAsymLimits[0]);
            h2d_csXbr_exp.SetBinContent(j+1,i+1,curAsymLimits[3]*mass_XS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]]);
            h2d_csXbr_obs.SetBinContent(j+1,i+1,curAsymLimits[0]*mass_XS[massindex[mass]]*cprime[j]*0.1*(1-BRnew[i]*0.1)*massBRWW[massindex[mass]]);
            

    h2d_obs.GetXaxis().SetTitle("C^{'2}");
    h2d_csXbr_exp.GetXaxis().SetTitle("C^{'2}");
    h2d_csXbr_obs.GetXaxis().SetTitle("C^{'2}");

    h2d_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_obs.GetYaxis().SetTitle("BR_{new}");
    h2d_csXbr_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_csXbr_obs.GetYaxis().SetTitle("BR_{new}");

    h2d_csXbr_obs.GetZaxis().SetLimits(0,h2d_csXbr_obs.GetMaximum());
    h2d_csXbr_exp.GetZaxis().SetLimits(0,h2d_csXbr_exp.GetMaximum());
    h2d_exp.GetZaxis().SetLimits(0,h2d_exp.GetMaximum());
    h2d_obs.GetZaxis().SetLimits(0,h2d_obs.GetMaximum());


    h2d_exp.SetLineStyle(0);
    h2d_exp.SetMarkerStyle(20);
    h2d_exp.GetXaxis().SetTitle("C^{'2}");
    h2d_exp.GetXaxis().SetNdivisions(504);
    h2d_exp.GetXaxis().SetLabelFont(42);
    h2d_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_exp.GetXaxis().SetLabelSize(0.045);
    h2d_exp.GetXaxis().SetTitleSize(0.045);
    h2d_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_exp.GetXaxis().SetTitleFont(42);
    h2d_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_exp.GetYaxis().SetNdivisions(9);
    h2d_exp.GetYaxis().SetLabelFont(42);
    h2d_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_exp.GetYaxis().SetLabelSize(0.045);
    h2d_exp.GetYaxis().SetTitleSize(0.045);
    h2d_exp.GetYaxis().SetTitleOffset(1.25);
    h2d_exp.GetYaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitle("signal strenght excluded at 95% C.L");
    h2d_exp.GetZaxis().SetLabelFont(42);
    h2d_exp.GetZaxis().SetLabelSize(0.03);
    h2d_exp.GetZaxis().SetTitleOffset(0.95);
    h2d_exp.GetZaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitleSize(0.05);

    h2d_obs.SetLineStyle(0);
    h2d_obs.SetMarkerStyle(20);
    h2d_obs.GetXaxis().SetTitle("C^{'2}");
    h2d_obs.GetXaxis().SetNdivisions(504);
    h2d_obs.GetXaxis().SetLabelFont(42);
    h2d_obs.GetXaxis().SetLabelOffset(0.007);
    h2d_obs.GetXaxis().SetLabelSize(0.045);
    h2d_obs.GetXaxis().SetTitleSize(0.045);
    h2d_obs.GetXaxis().SetTitleOffset(1.02);
    h2d_obs.GetXaxis().SetTitleFont(42);
    h2d_obs.GetYaxis().SetTitle("BR_{new}");
    h2d_obs.GetYaxis().SetNdivisions(9);
    h2d_obs.GetYaxis().SetLabelFont(42);
    h2d_obs.GetYaxis().SetLabelOffset(0.007);
    h2d_obs.GetYaxis().SetLabelSize(0.045);
    h2d_obs.GetYaxis().SetTitleSize(0.045);
    h2d_obs.GetYaxis().SetTitleOffset(1.25);
    h2d_obs.GetYaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitle("signal strenght excluded at 95% C.L");
    h2d_obs.GetZaxis().SetLabelFont(42);
    h2d_obs.GetZaxis().SetLabelSize(0.03);
    h2d_obs.GetZaxis().SetTitleOffset(0.95);
    h2d_obs.GetZaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitleSize(0.05);


    banner = TPaveText( 0.125, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.033);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.45,0.953,0.89,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.033);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
    banner2 = TLatex(0.17,0.91,("Higgs mass, %i GeV"%(mass)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can1_BSM2D = ROOT.TCanvas("can1_BSM2D","can1_BSM2D",1,1,600,600);
    can1_BSM2D.SetHighLightColor(2);
    can1_BSM2D.SetFillColor(0);
    can1_BSM2D.SetBorderMode(0);
    can1_BSM2D.SetBorderSize(2);
    can1_BSM2D.SetGridx();
    can1_BSM2D.SetGridy();
    can1_BSM2D.SetTickx(1);
    can1_BSM2D.SetTicky(1);
    can1_BSM2D.SetLeftMargin(0.15);
    can1_BSM2D.SetRightMargin(0.17);
    can1_BSM2D.SetTopMargin(0.05);
    can1_BSM2D.SetBottomMargin(0.1);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);

    h2d_exp.Draw("colz");
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();

    gPad.Update();
    palette = h2d_exp.GetListOfFunctions().FindObject("palette"); 
    myx1=palette.GetX1NDC() +0.01;
    myx2=palette.GetX2NDC() +0.01;
    palette.SetX1NDC(myx1);
    palette.SetX2NDC(myx2);
                                                                                             
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.png"%(SIGCH,mass));
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_%i.pdf"%(SIGCH,mass));

    if not options.blindObservedLine:
        can2_BSM2D = ROOT.TCanvas("can2_BSM2D","can2_BSM2D",1,1,600,600);
        can2_BSM2D.SetHighLightColor(2);
        can2_BSM2D.SetFillColor(0);
        can2_BSM2D.SetBorderMode(0);
        can2_BSM2D.SetBorderSize(2);
        can2_BSM2D.SetGridx();
        can2_BSM2D.SetGridy();
        can2_BSM2D.SetTickx(1);
        can2_BSM2D.SetTicky(1);
        can2_BSM2D.SetLeftMargin(0.15);
        can2_BSM2D.SetRightMargin(0.17);
        can2_BSM2D.SetTopMargin(0.05);
        can2_BSM2D.SetBottomMargin(0.1);
        can2_BSM2D.SetFrameFillStyle(0);
        can2_BSM2D.SetFrameBorderMode(0);
        can2_BSM2D.SetFrameFillStyle(0);
        can2_BSM2D.SetFrameBorderMode(0);
 
        h2d_obs.Draw("colz");
        banner.Draw();
        label_sqrt.Draw();
        banner2.Draw();

        gPad.Update();
        palette = h2d_obs.GetListOfFunctions().FindObject("palette"); 
        myx1=palette.GetX1NDC() +0.01;
        myx2=palette.GetX2NDC() +0.01;
        palette.SetX1NDC(myx1);
        palette.SetX2NDC(myx2);
        
        can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.png"%(SIGCH,mass));
        can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_%i.pdf"%(SIGCH,mass));

    can3_BSM2D = ROOT.TCanvas("can3_BSM2D","can3_BSM2D",1,1,630,576);
    can3_BSM2D.SetHighLightColor(2);
    can3_BSM2D.SetFillColor(0);
    can3_BSM2D.SetBorderMode(0);
    can3_BSM2D.SetBorderSize(2);
    can3_BSM2D.SetGridx();
    can3_BSM2D.SetGridy();
    can3_BSM2D.SetTickx(1);
    can3_BSM2D.SetTicky(1);
    can3_BSM2D.SetLeftMargin(0.15);
    can3_BSM2D.SetRightMargin(0.17);
    can3_BSM2D.SetTopMargin(0.05);
    can3_BSM2D.SetBottomMargin(0.15);
    can3_BSM2D.SetFrameFillStyle(0);
    can3_BSM2D.SetFrameBorderMode(0);
    can3_BSM2D.SetFrameFillStyle(0);
    can3_BSM2D.SetFrameBorderMode(0);

    h2d_csXbr_exp.SetLineStyle(0);
    h2d_csXbr_exp.SetMarkerStyle(20);
    h2d_csXbr_exp.GetXaxis().SetTitle("C^{'}");
    h2d_csXbr_exp.GetXaxis().SetNdivisions(504);
    h2d_csXbr_exp.GetXaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_csXbr_exp.GetXaxis().SetLabelSize(0.045);
    h2d_csXbr_exp.GetXaxis().SetTitleSize(0.045);
    h2d_csXbr_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_csXbr_exp.GetXaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetYaxis().SetTitle("BR_{new}");
    h2d_csXbr_exp.GetYaxis().SetNdivisions(9);
    h2d_csXbr_exp.GetYaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_csXbr_exp.GetYaxis().SetLabelSize(0.045);
    h2d_csXbr_exp.GetYaxis().SetTitleSize(0.045);
    h2d_csXbr_exp.GetYaxis().SetTitleOffset(1.25);
    h2d_csXbr_exp.GetYaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetZaxis().SetTitle("#sigma #times BR_{WW}  excluded at 95% C.L");
    h2d_csXbr_exp.GetZaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetZaxis().SetLabelSize(0.03);
    h2d_csXbr_exp.GetZaxis().SetTitleOffset(1);
    h2d_csXbr_exp.GetZaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetZaxis().SetTitleSize(0.05);

    h2d_csXbr_exp.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_%i.png"%(SIGCH,mass));
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_%i.pdf"%(SIGCH,mass));

    if not options.blindObservedLine:
        can4_BSM2D = ROOT.TCanvas("can4_BSM2D","can4_BSM2D",630,600);
        can4_BSM2D.SetHighLightColor(2);
        can4_BSM2D.SetFillColor(0);
        can4_BSM2D.SetBorderMode(0);
        can4_BSM2D.SetBorderSize(2);
        can4_BSM2D.SetGridx();
        can4_BSM2D.SetGridy();
        can4_BSM2D.SetTickx(1);
        can4_BSM2D.SetTicky(1);
        can4_BSM2D.SetLeftMargin(0.15);
        can4_BSM2D.SetRightMargin(0.17);
        can4_BSM2D.SetTopMargin(0.05);
        can4_BSM2D.SetBottomMargin(0.15);
        can4_BSM2D.SetFrameFillStyle(0);
        can4_BSM2D.SetFrameBorderMode(0);
        can4_BSM2D.SetFrameFillStyle(0);
        can4_BSM2D.SetFrameBorderMode(0);
        h2d_csXbr_obs.Draw("colz");
        banner.Draw();
        banner2.Draw();
        can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_%i.png"%(SIGCH,mass));
        can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_%i.pdf"%(SIGCH,mass));


#############################################
### Make the BSM 2D Scane vs c' and brNew ###
#############################################

def makeBSMLimitPlot2DBRnew( SIGCH, brNew):


    stylePath = os.getenv("ROOTStyle");
#    gROOT.ProcessLine(".x "+stylePath+"/rootPalette.C");

    massBRWW = [5.58E-01,5.77E-01,5.94E-01,6.09E-01,6.21E-01];            

    print "module ===> makeBSMLimits_2DMass";

    massindex = {600:0,700:1,800:2,900:3,1000:4}
    mass_XS  = [];
    cprime_x = array('d',[0.1,0.2,0.3,0.5,0.75]);
    mass_x   = array('d',[550,650,750,850,950,1050]);

    setStyle();

    if SIGCH == "" or SIGCH == "_2jet":
        mass_XS.append((0.5230+0.09688));
        mass_XS.append((0.2288+0.06330));
        mass_XS.append((0.1095+0.04365));
        mass_XS.append((0.05684+0.03164));
        mass_XS.append((0.03163+0.02399));
    elif SIGCH == "_ggH":
        mass_XS.append((0.5230));
        mass_XS.append((0.2288));
        mass_XS.append((0.1095));
        mass_XS.append((0.05684));
        mass_XS.append((0.03163));
    elif SIGCH == "_vbfH":
        mass_XS.append((0.09688));
        mass_XS.append((0.06330));
        mass_XS.append((0.04365));
        mass_XS.append((0.03164));
        mass_XS.append((0.02399));
    else:
        print "problem!"
        
    h2d_exp = TH2D("h2d_exp_brNew","",len(mass_x)-1,mass_x,len(cprime_x)-1,cprime_x);
    h2d_obs = TH2D("h2d_obs_brNew","",len(mass_x)-1,mass_x,len(cprime_x)-1,cprime_x);
    h2d_csXbr_exp = TH2D("h2d_csXbr_exp_brNew","",len(mass_x)-1,mass_x,len(cprime_x)-1,cprime_x);
    h2d_csXbr_obs = TH2D("h2d_csXbr_obs_brNew","",len(mass_x)-1,mass_x,len(cprime_x)-1,cprime_x);

    for j in range(len(mass)):
        for i in range(len(cprime)):
            curFile = "higgsCombinehwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(mass[j],options.channel,SIGCH,cprime[i],brNew,mass[j]);
            curAsymLimits = getAsymLimits(curFile);
            h2d_exp.SetBinContent(j+1,i+1,curAsymLimits[3]);
            h2d_obs.SetBinContent(j+1,i+1,curAsymLimits[0]);
            h2d_csXbr_exp.SetBinContent(j+1,i+1,curAsymLimits[3]*mass_XS[massindex[mass[j]]]*cprime[i]*0.1*(1-brNew*0.1)*massBRWW[massindex[mass[j]]]);
            h2d_csXbr_obs.SetBinContent(j+1,i+1,curAsymLimits[0]*mass_XS[massindex[mass[j]]]*cprime[i]*0.1*(1-brNew*0.1)*massBRWW[massindex[mass[j]]]);


    h2d_obs.GetYaxis().SetTitle("C^{'2}");
    h2d_csXbr_exp.GetYaxis().SetTitle("C^{'2}");
    h2d_csXbr_obs.GetYaxis().SetTitle("C^{'2}");
    h2d_exp.GetYaxis().SetTitle("C^{'2}");

    h2d_exp.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_obs.GetXaxis().SetTitle("m_{H} (GeV)");
    h2d_csXbr_exp.GetYaxis().SetTitle("C^{'2}");
    h2d_csXbr_obs.GetYaxis().SetTitle("C^{'2}");

    h2d_csXbr_obs.GetZaxis().SetLimits(0,h2d_csXbr_obs.GetMaximum());
    h2d_csXbr_exp.GetZaxis().SetLimits(0,h2d_csXbr_exp.GetMaximum());
    h2d_exp.GetZaxis().SetLimits(0,h2d_exp.GetMaximum());
    h2d_obs.GetZaxis().SetLimits(0,h2d_obs.GetMaximum());


    h2d_exp.SetLineStyle(0);
    h2d_exp.SetMarkerStyle(20);
    h2d_exp.GetXaxis().SetNdivisions(504);
    h2d_exp.GetXaxis().SetLabelFont(42);
    h2d_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_exp.GetXaxis().SetLabelSize(0.045);
    h2d_exp.GetXaxis().SetTitleSize(0.045);
    h2d_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_exp.GetXaxis().SetTitleFont(42);
    h2d_exp.GetYaxis().SetNdivisions(9);
    h2d_exp.GetYaxis().SetLabelFont(42);
    h2d_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_exp.GetYaxis().SetLabelSize(0.045);
    h2d_exp.GetYaxis().SetTitleSize(0.045);
    h2d_exp.GetYaxis().SetTitleOffset(1.25);
    h2d_exp.GetYaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitle("signal strenght excluded at 95% C.L");
    h2d_exp.GetZaxis().SetLabelFont(42);
    h2d_exp.GetZaxis().SetLabelSize(0.03);
    h2d_exp.GetZaxis().SetTitleOffset(0.95);
    h2d_exp.GetZaxis().SetTitleFont(42);
    h2d_exp.GetZaxis().SetTitleSize(0.05);

    h2d_obs.SetLineStyle(0);
    h2d_obs.SetMarkerStyle(20);
    h2d_obs.GetXaxis().SetNdivisions(504);
    h2d_obs.GetXaxis().SetLabelFont(42);
    h2d_obs.GetXaxis().SetLabelOffset(0.007);
    h2d_obs.GetXaxis().SetLabelSize(0.045);
    h2d_obs.GetXaxis().SetTitleSize(0.045);
    h2d_obs.GetXaxis().SetTitleOffset(1.02);
    h2d_obs.GetXaxis().SetTitleFont(42);
    h2d_obs.GetYaxis().SetNdivisions(9);
    h2d_obs.GetYaxis().SetLabelFont(42);
    h2d_obs.GetYaxis().SetLabelOffset(0.007);
    h2d_obs.GetYaxis().SetLabelSize(0.045);
    h2d_obs.GetYaxis().SetTitleSize(0.045);
    h2d_obs.GetYaxis().SetTitleOffset(1.25);
    h2d_obs.GetYaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitle("signal strenght excluded at 95% C.L");
    h2d_obs.GetZaxis().SetLabelFont(42);
    h2d_obs.GetZaxis().SetLabelSize(0.03);
    h2d_obs.GetZaxis().SetTitleOffset(0.95);
    h2d_obs.GetZaxis().SetTitleFont(42);
    h2d_obs.GetZaxis().SetTitleSize(0.05);


    banner = TPaveText( 0.125, 0.953, 0.76, 0.975, "brNDC");
    banner.SetFillColor(ROOT.kWhite);
    banner.SetTextSize(0.033);
    banner.SetTextAlign(11);
    banner.SetTextFont(62);
    banner.SetBorderSize(0);
    leftText = "CMS Preliminary";
    banner.AddText(leftText);
    banner.Draw();

    label_sqrt = TPaveText(0.45,0.953,0.89,0.975, "brNDC");
    label_sqrt.SetFillColor(ROOT.kWhite);
    label_sqrt.SetBorderSize(0);
    label_sqrt.SetTextSize(0.033);
    label_sqrt.SetTextFont(62);
    label_sqrt.SetTextAlign(31); # align right                                                                                                                                         
    label_sqrt.AddText("L = 19.3 fb^{-1} at #sqrt{s} = 8 TeV");
    
    banner2 = TLatex(0.17,0.91,("BR new = %0.1f"%(brNew)));
    banner2.SetNDC(); banner2.SetTextSize(0.028);

    can1_BSM2D = ROOT.TCanvas("can1_BSM2D","can1_BSM2D",1,1,600,600);
    can1_BSM2D.SetHighLightColor(2);
    can1_BSM2D.SetFillColor(0);
    can1_BSM2D.SetBorderMode(0);
    can1_BSM2D.SetBorderSize(2);
    can1_BSM2D.SetGridx();
    can1_BSM2D.SetGridy();
    can1_BSM2D.SetTickx(1);
    can1_BSM2D.SetTicky(1);
    can1_BSM2D.SetLeftMargin(0.15);
    can1_BSM2D.SetRightMargin(0.17);
    can1_BSM2D.SetTopMargin(0.05);
    can1_BSM2D.SetBottomMargin(0.1);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);
    can1_BSM2D.SetFrameFillStyle(0);
    can1_BSM2D.SetFrameBorderMode(0);

    h2d_exp.Draw("colz");
    banner.Draw();
    label_sqrt.Draw();
    banner2.Draw();

    gPad.Update();
    palette = h2d_exp.GetListOfFunctions().FindObject("palette"); 
    myx1=palette.GetX1NDC() +0.01;
    myx2=palette.GetX2NDC() +0.01;
    palette.SetX1NDC(myx1);
    palette.SetX2NDC(myx2);
                                                                                             
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_brNew_%0.1f.png"%(SIGCH,brNew));
    can1_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpMu_brNew_%0.1f.pdf"%(SIGCH,brNew));

    if not options.blindObservedLine:
        can2_BSM2D = ROOT.TCanvas("can2_BSM2D","can2_BSM2D",1,1,600,600);
        can2_BSM2D.SetHighLightColor(2);
        can2_BSM2D.SetFillColor(0);
        can2_BSM2D.SetBorderMode(0);
        can2_BSM2D.SetBorderSize(2);
        can2_BSM2D.SetGridx();
        can2_BSM2D.SetGridy();
        can2_BSM2D.SetTickx(1);
        can2_BSM2D.SetTicky(1);
        can2_BSM2D.SetLeftMargin(0.15);
        can2_BSM2D.SetRightMargin(0.17);
        can2_BSM2D.SetTopMargin(0.05);
        can2_BSM2D.SetBottomMargin(0.1);
        can2_BSM2D.SetFrameFillStyle(0);
        can2_BSM2D.SetFrameBorderMode(0);
        can2_BSM2D.SetFrameFillStyle(0);
        can2_BSM2D.SetFrameBorderMode(0);

        h2d_obs.Draw("colz");
        banner.Draw();
        banner2.Draw();
        can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_brNew_%0.1f.png"%(SIGCH,brNew));
        can2_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsMu_brNew_%0.1f.pdf"%(SIGCH,brNew));

    can3_BSM2D = ROOT.TCanvas("can3_BSM2D","can3_BSM2D",1,1,630,576);
    can3_BSM2D.SetHighLightColor(2);
    can3_BSM2D.SetFillColor(0);
    can3_BSM2D.SetBorderMode(0);
    can3_BSM2D.SetBorderSize(2);
    can3_BSM2D.SetGridx();
    can3_BSM2D.SetGridy();
    can3_BSM2D.SetTickx(1);
    can3_BSM2D.SetTicky(1);
    can3_BSM2D.SetLeftMargin(0.15);
    can3_BSM2D.SetRightMargin(0.17);
    can3_BSM2D.SetTopMargin(0.05);
    can3_BSM2D.SetBottomMargin(0.15);
    can3_BSM2D.SetFrameFillStyle(0);
    can3_BSM2D.SetFrameBorderMode(0);
    can3_BSM2D.SetFrameFillStyle(0);
    can3_BSM2D.SetFrameBorderMode(0);

    h2d_csXbr_exp.SetLineStyle(0);
    h2d_csXbr_exp.SetMarkerStyle(20);
    h2d_csXbr_exp.GetXaxis().SetNdivisions(504);
    h2d_csXbr_exp.GetXaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetXaxis().SetLabelOffset(0.007);
    h2d_csXbr_exp.GetXaxis().SetLabelSize(0.045);
    h2d_csXbr_exp.GetXaxis().SetTitleSize(0.045);
    h2d_csXbr_exp.GetXaxis().SetTitleOffset(1.02);
    h2d_csXbr_exp.GetXaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetYaxis().SetNdivisions(9);
    h2d_csXbr_exp.GetYaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetYaxis().SetLabelOffset(0.007);
    h2d_csXbr_exp.GetYaxis().SetLabelSize(0.045);
    h2d_csXbr_exp.GetYaxis().SetTitleSize(0.045);
    h2d_csXbr_exp.GetYaxis().SetTitleOffset(1.25);
    h2d_csXbr_exp.GetYaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetZaxis().SetTitle("#sigma #times BR_{WW}  excluded at 95% C.L");
    h2d_csXbr_exp.GetZaxis().SetLabelFont(42);
    h2d_csXbr_exp.GetZaxis().SetLabelSize(0.03);
    h2d_csXbr_exp.GetZaxis().SetTitleOffset(1);
    h2d_csXbr_exp.GetZaxis().SetTitleFont(42);
    h2d_csXbr_exp.GetZaxis().SetTitleSize(0.05);

    h2d_csXbr_exp.Draw("colz");
    banner.Draw();
    banner2.Draw();
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_brNew_%0.1f.png"%(SIGCH,brNew));
    can3_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ExpSigma_brNew_%0.1f.pdf"%(SIGCH,brNew));

    if not options.blindObservedLine:
        can4_BSM2D = ROOT.TCanvas("can4_BSM2D","can4_BSM2D",630,600);
        can4_BSM2D.SetHighLightColor(2);
        can4_BSM2D.SetFillColor(0);
        can4_BSM2D.SetBorderMode(0);
        can4_BSM2D.SetBorderSize(2);
        can4_BSM2D.SetGridx();
        can4_BSM2D.SetGridy();
        can4_BSM2D.SetTickx(1);
        can4_BSM2D.SetTicky(1);
        can4_BSM2D.SetLeftMargin(0.15);
        can4_BSM2D.SetRightMargin(0.17);
        can4_BSM2D.SetTopMargin(0.05);
        can4_BSM2D.SetBottomMargin(0.15);
        can4_BSM2D.SetFrameFillStyle(0);
        can4_BSM2D.SetFrameBorderMode(0);
        can4_BSM2D.SetFrameFillStyle(0);
        can4_BSM2D.SetFrameBorderMode(0);
        h2d_csXbr_obs.Draw("colz");
        banner.Draw();
        banner2.Draw();
        can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_brNew_%0.1f.png"%(SIGCH,brNew));
        can4_BSM2D.SaveAs("limitFigs/BSMLim%s_2D_ObsSigma_brNew_%0.1f.pdf"%(SIGCH,brNew));

##################################
########### Main Code ############
##################################    

if __name__ == '__main__':

        
    CHAN = options.channel;
    DIR  = options.datacardDIR;
    if options.sigChannel !="": 
     SIGCH = options.jetBin+"_"+options.sigChannel;
    else:
     SIGCH = options.jetBin;
        
    moreCombineOpts = "";
    print "channel ",CHAN," directiory ",DIR," signal channel ",SIGCH," more options ",moreCombineOpts ;
    
    if options.makeCards and DIR != "":
        if not os.path.isdir(DIR):
            os.system("mkdir "+DIR);
        else: 
            print "Directory "+DIR+" already exists...";

    if options.computeLimits or options.plotLimits: os.chdir(DIR);

    nMasses = len(mass);
    mLo = 0;
    mHi = nMasses;
    nCprimes = len(cprime);
    cpLo = 0;
    cpHi = nCprimes;
    nBRnews = len(BRnew);
    brLo = 0;
    brHi = nBRnews;
    
    if options.massPoint > 0:   
        curIndex = mass.index(options.massPoint);
        mLo = curIndex;
        mHi = curIndex+1;
    if options.cPrime > 0:   
        curIndex = cprime.index(options.cPrime);
        cpLo = curIndex;
        cpHi = curIndex+1;
        nCprimes = 1;
    if options.brNew >= 0:   
        curIndex = BRnew.index(options.brNew);
        brLo = curIndex;
        brHi = curIndex+1;
        nBRnews = 1;    


    if options.injectSingalStrenght == 0 : rMin = -10 ; rMax = 10 ;
    else: rMin = - 50 ; rMax = 50

    #######################################################
    # ======= make the datacards running the fit =======  #
    #######################################################
    
    if options.makeCards:
        if not os.path.isdir("log"): os.system("mkdir log" );
        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";                
                    print "R U N N I N G   F I T S" 
                    print "mass = ",mass[i],", cprime = ",cprime[j],", brnew = ",BRnew[k],", channel: ",options.channel," pseudodata ",options.pseudodata
                    print "--------------------------------------------------";                
                    print "--------------------------------------------------";  
                    
                    time.sleep(0.3);
                    
                    command = "python doFit_class_higgs.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew %01d --inPath %s/ --jetBin %s --channel %s --pseudodata %d --closuretest %d --skipJetSystematics %d --interferenceModel %s"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], os.getcwd(), options.jetBin, options.channel,options.pseudodata,options.closuretest,options.skipJetSystematics,options.interferenceModel);
                    print command ;

                    if options.batchMode :
                        fn = "fitScript_%s_%03d_%02d_%02d_%s"%(options.channel,mass[i],cprime[j],BRnew[k],shape[i]);
                        submitBatchJob( command, fn );
                    if not options.batchMode: 
                        os.system(command);

    #########################################
    # ===================================== #
    #########################################

    if options.computeLimits:

        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";
                    print "analyzing card: hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                    print "--------------------------------------------------";                

                    ###############################################
                    ####### Combine single channel cards ##########
                    ###############################################


                    if options.channel =="em" and options.jetBin != "_2jet" :
                     combineCmmd = "combineCards.py hwwlvj_ggH%03d_el%s_%02d_%02d_unbin.txt hwwlvj_ggH%03d_mu%s_%02d_%02d_unbin.txt > hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k]);
                     print "combineCmmd ",combineCmmd; 
                     os.system(combineCmmd);

                    if options.higgsCombination == 1:
                     options.channel = "combo" ; 
                     combineCmmd = "combineCards.py hwwlvj_ggH%03d_mu%s_%02d_%02d_unbin.txt hwwlvj_ggH%03d_el%s_%02d_%02d_unbin.txt hwwlvj_ggH%03d_em_2jet%s_%02d_%02d_unbin.txt > hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k],mass[i],options.channel,SIGCH,cprime[j],BRnew[k]);
                     print "combineCmmd ",combineCmmd; 
                     os.system(combineCmmd);
                        

                    ########################################
                    ####### Generate only options ##########
                    ########################################

                    os.system("mkdir "+options.inputGeneratedDataset);
                     
                    if options.generateOnly == 1 :
                      if options.outputTree == 0 :
                          for iToy in range(options.nToys):                              
                           runCmmd =  "combine -M GenerateOnly --saveToys -s -1 -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -t 1 --expectSignal=%d "%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.injectSingalStrenght);
                           print "runCmmd ",runCmmd;
                           if options.batchMode:
                              fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                              cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                           else: 
                              os.system(runCmmd);
                              os.system("mv higgsCombine* "+options.inputGeneratedDataset);   
                           continue ;

                      else:
                           runCmmd =  "combine -M GenerateOnly --saveToys -s -1 -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -t %d --expectSignal=%d "%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.nToys,options.injectSingalStrenght);
                           print "runCmmd ",runCmmd;
                           if options.batchMode:
                              fn = "combineScript_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                              cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                           else: 
                              os.system(runCmmd);
                              os.system("mv higgsCombine* "+options.inputGeneratedDataset);   
                           continue ;

                    ######################################
                    ####### Maximum Likelihood Fits ######
                    ######################################   
    
                    elif options.maximumLikelihoodFit == 1:
                        
                       ################################################# 
                       #### just one fit with the defined datacards ####
                       #################################################
                        
                       if options.nToys == 0 and options.crossedToys == 0 : 
                        runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties  -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2  --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);                     
                        print "runCmmd ",runCmmd;
                        if options.batchMode:
                           fn = "combineScript_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                           submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                        else:   
                         os.system(runCmmd);
                         
                       ######################################################## 
                       #### run many toys and fits with the same datacards  ###
                       ########################################################

                       elif options.nToys != 0 and options.crossedToys == 0 :
                          if options.outputTree == 0:  
                           for iToy in range(options.nToys):
                             runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveToys --saveWithUncertainties --toysNoSystematics -s -1 -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin_%d -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -t 1 --expectSignal=%d --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,cprime[j],BRnew[k],iToy,mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.injectSingalStrenght);                     
                             print "runCmmd ",runCmmd;
                             if options.batchMode:
                              fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                              cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                             else: 
                              os.system(runCmmd);
                           continue ;
                          else:
                             runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties  --toysNoSystematics --saveToys -s -1 -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -t %d --expectSignal=%d --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.nToys,options.injectSingalStrenght);                     
                             print "runCmmd ",runCmmd;
                             if options.batchMode:
                              fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],options.nToys);
                              cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                              submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                             else: 
                              os.system(runCmmd);
                             continue ;

                       ##################################
                       ###### Make the crossed toys ##### 
                       ##################################  

                       elif options.nToys != 0 and options.crossedToys == 1 :

                          os.system("ls "+options.inputGeneratedDataset+" | grep root | grep higgsCombine | grep ggH"+str(mass[i])+" > list_temp.txt"); 
                          iToy = 0 ;
                          if options.outputTree == 0:  
                           with open("list_temp.txt") as input_list:
                            for line in input_list:
                             for name in line.split():
                                if iToy >= options.nToys: continue ; 
                                runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin_%d -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -s -1 -t 1  --toysFile %s/%s --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,cprime[j],BRnew[k],iToy,mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.inputGeneratedDataset,name);
                                iToy = iToy + 1 ;
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                                  cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

                          else:
                           with open("list_temp.txt") as input_list:
                              for line in input_list:
                               for name in line.split():                                                                       
                                runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --rMin %d --rMax %d --saveNormalizations --saveWithUncertainties -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -s -1 -t %d --toysFile %s/%s --robustFit=1 --do95=1"%(rMin,rMax,mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts,options.nToys,options.inputGeneratedDataset,name);                     
                                print "runCmmd ",runCmmd;                                
                                if options.batchMode:
                                  fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                                  cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                                  submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                                else: 
                                  os.system(runCmmd);

                          os.system("rm list_temp.txt")
                          continue ;

                    ###############################
                    #### Asymptotic Limit part  ###
                    ###############################
                      
                    elif options.systematics == 0 and not options.computePvalue == 1 and not options.makeLikelihoodScan == 1:
                       runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -S 0"%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);
                       print "runCmmd ",runCmmd ;

                       if options.batchMode:
                        fn = "combineScript_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                        cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                        submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                       else: 
                        os.system(runCmmd);


                    elif not options.computePvalue == 1 and not options.makeLikelihoodScan == 1:

                       #############################################
                       ###### run Asymptotic on the final card ##### 
                       #############################################  

                       if options.nToys == 0: 
                        runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);                                        
                        print "runCmmd ",runCmmd;

                        if options.batchMode:
                         fn = "combineScript_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                         submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                        else: 
                         os.system(runCmmd);

                       else:

                       #############################################
                       ###### run Asymptotic making many toys  ##### 
                       #############################################  

                        for iToy in range(options.nToys):   
                         runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin_%d -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 -t 1 -s -1"%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],iToy,mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);                                        
                         print "runCmmd ",runCmmd;

                         if options.batchMode:
                          fn = "combineScript_%s_%03d%s_%02d_%02d_iToy%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                          submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                         else: 
                          os.system(runCmmd);
                           
                      
                    elif options.computePvalue == 1 and not options.makeLikelihoodScan == 1: 

                       ##################################################
                       ###### run the observed and expected pvalue  ##### 
                       ##################################################  

                        runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n hwwlvj_pval_obs_ggH%03d_%s%s_%02d_%02d_unbin -m %03d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2"%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);
                        print "runCmmd ",runCmmd;

                        if options.batchMode:
                         fn = "combineScript_ProfileLikelihood_obs_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                         submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                        else:
                         os.system(runCmmd);

                        for iToy in range(options.nToys):
                         runCmmd = "combine -M ProfileLikelihood --signif --pvalue -n hwwlvj_pval_exp_ggH%03d_%s%s_%02d_%02d_unbin_%d -m %03d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt %s -v 2 --toysFreq -t 1 --expectSignal=1 -s -1 "%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],iToy,mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);
                         print "runCmmd ",runCmmd;

                         if options.batchMode:
                           fn = "combineScript_ProfileLikelihood_exp_%s_%03d%s_%02d_%02d_%d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k],iToy);
                           submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                         else:
                          os.system(runCmmd);

                    elif options.makeLikelihoodScan == 1:
                        
                         runCmmd = "combine -M MultiDimFit -n hwwlvj_LikelihoodScan_ggH%03d_%s%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_%s%s_%02d_%02d_unbin.txt  --algo=grid --points=150 --setPhysicsModelParameterRanges r=-1,5 %s"%(mass[i],options.channel,SIGCH,cprime[j],BRnew[k],mass[i],mass[i],options.channel,SIGCH,cprime[j],BRnew[k],moreCombineOpts);
                         print "runCmmd ",runCmmd;
 
                         if options.batchMode:
                          fn = "combineScript_LikelihoodScan_%s_%03d%s_%02d_%02d"%(options.channel,mass[i],SIGCH,cprime[j],BRnew[k]);
                          submitBatchJobCombine(runCmmd, fn, mass[i], cprime[j], BRnew[k]);
                         else:
                          os.system(runCmmd);
                                                                                                                               
    ####################################################
    # =================== Bias Analysis ============== #
    ####################################################
    
    if options.biasStudy:

        for i in range(mLo,mHi): 
            print "--------------------------------------------------";                
            print "--------------------------------------------------";                
            print "B I A S  S T U D Y   F I T S" 
            print "mass = ",mass[i]," channel: ",options.channel," pseudodata ",options.pseudodata
            print "--------------------------------------------------";                
            print "--------------------------------------------------";  

            if options.jetBin == "_2jet" :
             command = "python BiasStudy/do_fitBias_higgs.py ggH%d %d %d %d %d -b --pseudodata %d --fgen %s --fres %s --nexp %d --isMC %d --storeplot %d --channel %s --inPath %s --ttbarcontrolregion %d --fitjetmass %d --mlvjregion %s --onlybackgroundfit %d --inflatejobstatistic %d --scalesignalwidth %0.1f --injectSingalStrenght %0.1f --jetBin %s"%(mass[i],mlo[i],mhi[i],mjlo[i],mjhi[i],options.pseudodata,shape_gen[i],shape_fit[i],options.nToys,isMC[i],1,options.channel,os.getcwd(),options.ttbarcontrolregion,options.fitjetmass,options.mlvjregion,options.onlybackgroundfit,options.inflatejobstatistic,options.scalesignalwidth,options.injectSingalStrenght,options.jetBin);
            else:
             command = "python BiasStudy/do_fitBias_higgs.py ggH%d %d %d %d %d -b --pseudodata %d --fgen %s --fres %s --nexp %d --isMC %d --storeplot %d --channel %s --inPath %s --ttbarcontrolregion %d --fitjetmass %d --mlvjregion %s --onlybackgroundfit %d --inflatejobstatistic %d --scalesignalwidth %0.1f --injectSingalStrenght %0.1f "%(mass[i],mlo[i],mhi[i],mjlo[i],mjhi[i],options.pseudodata,shape_gen[i],shape_fit[i],options.nToys,isMC[i],1,options.channel,os.getcwd(),options.ttbarcontrolregion,options.fitjetmass,options.mlvjregion,options.onlybackgroundfit,options.inflatejobstatistic,options.scalesignalwidth,options.injectSingalStrenght);
                
            print command ;
            if options.batchMode:
             suffix = options.channel;
             if options.jetBin == "_2jet" : suffix = suffix+"_2jet";
             if options.ttbarcontrolregion : suffix = suffix+"_ttbar";
             if options.fitjetmass         : suffix = suffix+"_jetmass";
             if options.turnOnAnalysis     : suffix = suffix+"_turnOn";
             if options.scalesignalwidth !=1 : suffix = suffix+ ("_width_%0.1f")%(options.scalesignalwidth);
             if options.injectSingalStrenght !=0 : suffix = suffix+"_SB";
             else :suffix = suffix+"_B";
             if options.onlybackgroundfit  : suffix = suffix+"_B";
             else: suffix = suffix+"_SB";
             fn = "biasScript_ggH%03d_%s_%s%s"%(mass[i],shape_gen[i],shape_fit[i],suffix);
             submitBatchJob( command, fn );
            else: 
             os.system(command);

    #############################################################
    # =================== Plot of the Limit  ================== #
    #############################################################
    
    if options.plotLimits:

      if options.makeSMLimitPlot == 1:
          makeSMLimitPlot(SIGCH,10,00);
          makeSMLimitPlot(SIGCH,01,00);
          makeSMLimitPlot(SIGCH,01,01);
          makeSMLimitPlot(SIGCH,01,02);
          makeSMLimitPlot(SIGCH,01,03);
          makeSMLimitPlot(SIGCH,01,04);
          makeSMLimitPlot(SIGCH,01,05);
          makeSMLimitPlot(SIGCH,02,00);
          makeSMLimitPlot(SIGCH,02,01);
          makeSMLimitPlot(SIGCH,02,02);
          makeSMLimitPlot(SIGCH,02,03);
          makeSMLimitPlot(SIGCH,02,04);
          makeSMLimitPlot(SIGCH,02,05);
          makeSMLimitPlot(SIGCH,03,00);
          makeSMLimitPlot(SIGCH,03,01);
          makeSMLimitPlot(SIGCH,03,02);
          makeSMLimitPlot(SIGCH,03,03);
          makeSMLimitPlot(SIGCH,03,04);
          makeSMLimitPlot(SIGCH,03,05);
          makeSMLimitPlot(SIGCH,05,00);
          makeSMLimitPlot(SIGCH,05,01);
          makeSMLimitPlot(SIGCH,05,02);
          makeSMLimitPlot(SIGCH,05,03);
          makeSMLimitPlot(SIGCH,05,04);
          makeSMLimitPlot(SIGCH,05,05);
          makeSMLimitPlot(SIGCH,07,00);
          makeSMLimitPlot(SIGCH,07,01);
          makeSMLimitPlot(SIGCH,07,02);
          makeSMLimitPlot(SIGCH,07,03);


          if options.plotPValue == 1:
              makeSMPValuePlot(SIGCH,10,00);
              makeSMPValuePlot(SIGCH,10,01);
              makeSMPValuePlot(SIGCH,10,02);
              makeSMPValuePlot(SIGCH,10,03);
              makeSMPValuePlot(SIGCH,10,04);
              makeSMPValuePlot(SIGCH,10,05);

              makeSMPValuePlot(SIGCH,01,00);
              makeSMPValuePlot(SIGCH,01,01);
              makeSMPValuePlot(SIGCH,01,02);
              makeSMPValuePlot(SIGCH,01,03);
              makeSMPValuePlot(SIGCH,01,04);
              makeSMPValuePlot(SIGCH,01,05);
              makeSMPValuePlot(SIGCH,02,00);
              makeSMPValuePlot(SIGCH,02,01);
              makeSMPValuePlot(SIGCH,02,02);
              makeSMPValuePlot(SIGCH,02,03);
              makeSMPValuePlot(SIGCH,02,04);
              makeSMPValuePlot(SIGCH,02,05);
              makeSMPValuePlot(SIGCH,03,00);
              makeSMPValuePlot(SIGCH,03,01);
              makeSMPValuePlot(SIGCH,03,02);
              makeSMPValuePlot(SIGCH,03,03);
              makeSMPValuePlot(SIGCH,03,04);
              makeSMPValuePlot(SIGCH,03,05);
              makeSMPValuePlot(SIGCH,05,00);
              makeSMPValuePlot(SIGCH,05,01);
              makeSMPValuePlot(SIGCH,05,02);
              makeSMPValuePlot(SIGCH,05,03);
              makeSMPValuePlot(SIGCH,05,04);
              makeSMPValuePlot(SIGCH,05,05);
              makeSMPValuePlot(SIGCH,07,00);
              makeSMPValuePlot(SIGCH,07,01);
              makeSMPValuePlot(SIGCH,07,02);
              makeSMPValuePlot(SIGCH,07,03);

      if options.makeBSMLimitPlotMass == 1:
          makeBSMLimitPlotMass(SIGCH);
                 
      if options.makeBSMLimitPlotBRnew == 1:
          
          makeBSMLimitPlotBRnew(SIGCH,600);
          makeBSMLimitPlotBRnew(SIGCH,700);
          makeBSMLimitPlotBRnew(SIGCH,800);
          makeBSMLimitPlotBRnew(SIGCH,900);
          makeBSMLimitPlotBRnew(SIGCH,1000);

      if options.makeBSMLimitPlot2D == 1:
          makeBSMLimitPlot2D(SIGCH,600);
          makeBSMLimitPlot2D(SIGCH,700);
          makeBSMLimitPlot2D(SIGCH,800);
          makeBSMLimitPlot2D(SIGCH,900);
          makeBSMLimitPlot2D(SIGCH,1000);
      
          makeBSMLimitPlot2DBRnew(SIGCH,0);
          makeBSMLimitPlot2DBRnew(SIGCH,0.1);
          makeBSMLimitPlot2DBRnew(SIGCH,0.2);
          makeBSMLimitPlot2DBRnew(SIGCH,0.3);
          makeBSMLimitPlot2DBRnew(SIGCH,0.4);
          makeBSMLimitPlot2DBRnew(SIGCH,0.5);

      if options.plotSignalStrenght == 1:
          makeSignalStrenghtPlot(SIGCH,10,00);

      if options.plotLikelihoodScan == 1:
          makeLikelihoodScanPlot(SIGCH,10,00);
          
