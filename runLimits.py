#! /usr/bin/env python
import os
import glob
import math
import array
import sys
import time
import subprocess
import ROOT

from optparse import OptionParser
from subprocess import Popen
from ROOT import gROOT, gStyle, gSystem, TLatex

ROOT.gStyle.SetPadRightMargin(0.16);

############################################
#            Job steering                  #
############################################

parser = OptionParser()

parser.add_option('-b', action='store_true', dest='noX', default=False, help='no X11 windows')

parser.add_option('--makeCards',     action='store_true', dest='makeCards',         default=False, help='no X11 windows')
parser.add_option('--computeLimits', action='store_true', dest='computeLimits', default=False, help='no X11 windows')
parser.add_option('--plotLimits',    action='store_true', dest='plotLimits',       default=False, help='no X11 windows')
parser.add_option('--biasStudy',     action='store_true', dest='biasStudy',        default=False, help='no X11 windows')

# submit jobs to condor 
parser.add_option('--batchMode', action='store_true', dest='batchMode', default=False, help='no X11 windows')


parser.add_option('--channel',    action="store",type = "string", dest="channel",    default="em")
parser.add_option('--pseudodata', action="store",type = int,      dest="pseudodata", default = 0)
parser.add_option('--systematics',action="store",type = "int",    dest="systematics",default=1)
parser.add_option('--saveParam',action="store",type="int",dest="saveParam",default=1)
parser.add_option('--massPoint',action="store",type="int",dest="massPoint",default=-1)
parser.add_option('--cPrime',action="store",type="int",dest="cPrime",default=-1)
parser.add_option('--brNew',action="store",type="int",dest="brNew",default=-1)
parser.add_option('--odir',action="store",type="string",dest="odir",default=".")
parser.add_option('--sigChannel',action="store",type="string",dest="sigChannel",default="")

parser.add_option('--turnOnAnalysis', action="store",type="int",   dest="turnOnAnalysis",default=0)
parser.add_option('--shapetest',      action="store",type="int",   dest="shapetest",default=0)


(options, args) = parser.parse_args()

#############################################
######### Get Some Global Variables #########
#############################################

mass  =  [ 600, 700, 800, 900,1000]
ccmlo =  [ 550, 600, 700, 750, 800]
ccmhi =  [ 700, 850, 950,1100,1200]
mjlo  =  [ 40, 40, 40, 40, 40]
mjhi  =  [ 130, 130, 130, 130, 130]

if options.turnOnAnalysis:

 mlo      =  [ 400, 400, 400, 400, 400]
 mhi      =  [1500,1500,1500,1500,1500]
 shape    =  ["ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1","ErfExp_v1"]
 shapeAlt =  ["ErfPow2_v1", "ErfPow2_v1","ErfPow2_v1","ErfPow2_v1","ErfPow2_v1"]

else:

 mlo      =  [ 550, 550, 550, 550, 550]
 mhi      =  [1500,1500,1500,1500,1500]
 shape    =  ["Exp","Exp","Exp","Exp","Exp"]
 shapeAlt =  ["Pow","Pow","Pow","Pow","Pow"]

if options.biasStudy:

 shape_gen = ["Exp","Exp","Exp","Exp","Exp"]    
 shape_fit = ["Exp","Exp","Exp","Exp","Exp"]
 nexp      = [1000,1000,1000,1000,1000]
 isMC      = [0,0,0,0,0]

BRnew  = [0];
cprime = [10];
  
##############################################

def getAsymLimits(file):
        
 f = ROOT.TFile(file);
 t = f.Get("limit");
 entries = t.GetEntries();
    
 lims = [0,0,0,0,0,0];
    
 for i in range(entries):
        
  t.GetEntry(i);
  t_quantileExpected = t.quantileExpected;
  t_limit = t.limit;
        
  #print "limit: ", t_limit, ", quantileExpected: ",t_quantileExpected;
        
  if t_quantileExpected == -1.: lims[0] = t_limit;
  elif t_quantileExpected >= 0.024 and t_quantileExpected <= 0.026: lims[1] = t_limit;
  elif t_quantileExpected >= 0.15 and t_quantileExpected <= 0.17: lims[2] = t_limit;            
  elif t_quantileExpected == 0.5: lims[3] = t_limit;            
  elif t_quantileExpected >= 0.83 and t_quantileExpected <= 0.85: lims[4] = t_limit;
  elif t_quantileExpected >= 0.974 and t_quantileExpected <= 0.976: lims[5] = t_limit;
  else: print "Unknown quantile!"
    
  return lims;

######################################

def submitBatchJob( command, fn ):
    
 currentDir = os.getcwd();
    
 # create a dummy bash/csh
 outScript = open(fn+".sh","w");
    
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
    
 # link a condor script to your shell script
 condorScript = open("condor_"+fn,"w");
 condorScript.write('universe = vanilla')
 condorScript.write("\n"+"Executable = "+fn+".sh")
 condorScript.write("\n"+'Requirements = Memory >= 199 &&OpSys == "LINUX"&& (Arch != "DUMMY" )&& Disk > 1000000')
 condorScript.write("\n"+'Should_Transfer_Files = YES')
 condorScript.write("\n"+'Transfer_Input_Files = doFit_class.py, BiasStudy/do_fitBias_vbf.py')    
 condorScript.write("\n"+'WhenToTransferOutput  = ON_EXIT_OR_EVICT')
 condorScript.write("\n"+'Output = out_$(Cluster).stdout')
 condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
 condorScript.write("\n"+'Error  = out_$(Cluster).stderr')
 condorScript.write("\n"+'Log    = out_$(Cluster).log')
 condorScript.write("\n"+'Notification    = Error')
 condorScript.write("\n"+'Queue 1')
 condorScript.close();
    
 # submit the condor job 
 os.system("condor_submit "+"condor_"+fn)

# ----------------------------------------

def submitBatchJobCombine( command, fn, mass, cprime, BRnew ):
    
    
    SIGCH = "_"+options.sigChannel;    
    currentDir = os.getcwd();
    
    # create a dummy bash/csh
    outScript = open(fn+".sh","w");
    
    outScript.write('#!/bin/bash');
    outScript.write("\n"+'date');
    outScript.write("\n"+'source /uscmst1/prod/sw/cms/bashrc prod');
    outScript.write("\n"+'cd '+currentDir);
    outScript.write("\n"+'eval `scram runtime -sh`');
    outScript.write("\n"+'cd -');
    outScript.write("\n"+'ls');    
    outScript.write("\n"+command);
    outScript.close();
    
    file1 = "hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass,SIGCH,cprime,BRnew);
    file2 = "hwwlvj_ggH%03d_mu_%02d_%02d_workspace.root"%(mass,cprime,BRnew);
    file3 = "hwwlvj_ggH%03d_el_%02d_%02d_workspace.root"%(mass,cprime,BRnew);    

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

############################################################

if __name__ == '__main__':

    
    CHAN = options.channel;
    DIR  = "cards_"+CHAN;
    SIGCH = "_"+options.sigChannel;    
    moreCombineOpts = "";
    
    ###############    
    if options.makeCards:
        if not os.path.isdir(DIR): os.system("mkdir "+DIR);
        else: 
            print "Directory "+DIR+" already exists...";
            #sys.exit(0);

    if options.computeLimits or options.plotLimits: os.chdir(DIR);

    # put in functionality to test just one mass point or just one cprime
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

    # ======= make the datacards running the fit =======  #

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
                    
                    #command = "python doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew %01d --odir %s > %s/log/log_%s_ggH%03d_%02d_%02d_%02d_%02d_%02d_%02d_%s_%s_cprime_%02d_BRnew_%02d "%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], options.odir, options.odir, CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k]);
                    command = "python doFit_class.py %s ggH%03d %02d %02d %02d %02d %02d %02d %s %s -b -m --cprime %01d --BRnew %01d --inPath %s"%(CHAN, mass[i], ccmlo[i], ccmhi[i], mjlo[i], mjhi[i], mlo[i], mhi[i], shape[i], shapeAlt[i], cprime[j], BRnew[k], os.getcwd());
                    print command ;

                    unbinnedCard = options.odir+"/cards_%s/hwwlvj_ggH%03d_%s_%02d_%02d_unbin.txt"%(options.channel,mass[i],options.channel,cprime[j],BRnew[k]);
                    fileExists = os.path.isfile(unbinnedCard)
                    print "fileExists: ",fileExists,", cards: ", unbinnedCard
                    if options.batchMode and not fileExists:
                        fn = "fitScript_%s_%03d_%02d_%02d"%(options.channel,mass[i],cprime[j],BRnew[k]);
                        submitBatchJob( command, fn );
                    if not options.batchMode: 
                        os.system(command);

    # =====================================

    if options.computeLimits:

        for i in range(mLo,mHi):
            for j in range(cpLo,cpHi):
                for k in range(brLo,brHi):

                    print "--------------------------------------------------";
                    print "creating card: hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                    print "--------------------------------------------------";                

                    combineCmmd = "combineCards.py hwwlvj_ggH%03d_el%s_%02d_%02d_unbin.txt hwwlvj_ggH%03d_mu%s_%02d_%02d_unbin.txt > hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k],mass[i],SIGCH,cprime[j],BRnew[k]);
                    os.system(combineCmmd);

                    
                    print "running card: hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt"%(mass[i],SIGCH,cprime[j],BRnew[k]);

                    if options.saveParam == 1:
                       runCmmd =  "combine -M MaxLikelihoodFit --minimizerAlgo Minuit2 --minimizerStrategy 2 --saveWorkspace -n hwwlvj_ggH%03d_em%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt %s -v 2 -t 0 --expectSignal=0 "%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],mass[i],SIGCH,cprime[j],BRnew[k],moreCombineOpts);             

                    elif options.systematics == 0:
                       runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_em%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt %s -v 2 -S 0"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],mass[i],SIGCH,cprime[j],BRnew[k],moreCombineOpts);                    

                    else: 
                       runCmmd = "combine -M Asymptotic --minimizerAlgo Minuit2 --minosAlgo stepping -n hwwlvj_ggH%03d_em%s_%02d_%02d_unbin -m %03d -d hwwlvj_ggH%03d_em%s_%02d_%02d_unbin.txt %s -v 2"%(mass[i],SIGCH,cprime[j],BRnew[k],mass[i],mass[i],SIGCH,cprime[j],BRnew[k],moreCombineOpts);                                        
                    
                    if options.batchMode:
                        fn = "combineScript_%03d%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                        cardStem = "hwwlvj_ggH%03d_em%s_%02d_%02d"%(mass[i],SIGCH,cprime[j],BRnew[k]);
                        submitBatchJobCombine( runCmmd, fn, mass[i], cprime[j], BRnew[k] );
                    else: 
                        os.system(runCmmd);


    # =================== Bias Analysis ============== #
    if options.biasStudy:

        for i in range(mLo,mHi): 
            print "--------------------------------------------------";                
            print "--------------------------------------------------";                
            print "B I A S  S T U D Y   F I T S" 
            print "mass = ",mass[i]," channel: ",options.channel," pseudodata ",options.pseudodata
            print "--------------------------------------------------";                
            print "--------------------------------------------------";  

            runCmmd = "python do_fitBias_vbf.py ggH%03d %03d %03d %03d %03d -b --pseudodata %d --fgen %s --fres %s --nexp %d --isMC %d --storeplot %d --channel %s"%(mass[i],mlo[i],mhi[i],mjlo[i],mjhi[i],options.pseudodata,shape_gen[i],shape_fit[i],nexp[i],isMC[i],1,options.channel);             
            print runCmmd ;
            if options.batchMode:
              fn = "biasScript_ggH%03d_%s_%s"%(mass[i],shape_gen[i],shape_fit[i]);
              submitBatchJob( command, fn );
            else: 
              os.system(runCmmd);

    # =================== Plot of the Limit  ================== #

    if options.plotLimits:

       nGraphs = nCprimes*2 + 2;
       tGraphs = [];
       nPoints = len(mass);

       for j in range(cpHi-1,cpLo-1,-1):
            
            xbins     = array('d', [])
            ybins_exp = array('d', [])
            ybins_obs = array('d', [])            
            xbins_env = array('d', [])
            ybins_1s  = array('d', [])                        
            ybins_2s  = array('d', [])                        
            
            for i in range(mLo,mHi):

                curFile = "higgsCombinehwwlvj_ggH%03d_em_2jet_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
                curAsymLimits = getAsymLimits(curFile);
                
                xbins.append( mass[i] );
                ybins_exp.append( curAsymLimits[3] );                
                ybins_obs.append( curAsymLimits[0] );                                
                xbins_env.append( mass[i] );                                
                ybins_2s.append( curAsymLimits[1] );                                
                ybins_1s.append( curAsymLimits[2] );                                                

                print "mass: ", mass[i], ", cprime: ", cprime[j], " -- obs: ", curAsymLimits[0], ", exp: ", curAsymLimits[3];

            for i in range(mHi-1,mLo-1,-1):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_2jet_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
                curAsymLimits = getAsymLimits(curFile);
                
                xbins_env.append( mass[i] );                                
                ybins_2s.append( curAsymLimits[5] );                                
                ybins_1s.append( curAsymLimits[4] );                                                

            curGraph_exp = ROOT.TGraph(nPoints,xbins,ybins_exp);
            curGraph_obs = ROOT.TGraph(nPoints,xbins,ybins_obs);
            curGraph_1s = ROOT.TGraph(nPoints*2,xbins_env,ybins_1s);
            curGraph_2s = ROOT.TGraph(nPoints*2,xbins_env,ybins_2s);

            curGraph_exp.SetLineStyle(2);
            curGraph_exp.SetLineWidth(2);
            curGraph_obs.SetLineWidth(2);                    
            curGraph_exp.SetMarkerSize(2);
            curGraph_obs.SetMarkerSize(2);                    
            curGraph_1s.SetFillColor(ROOT.kGreen);
            curGraph_2s.SetFillColor(ROOT.kYellow);

            tGraphs.append(curGraph_exp);
            tGraphs.append(curGraph_obs);
            tGraphs.append(curGraph_1s);
            tGraphs.append(curGraph_2s);
    
       banner = TLatex(0.25,0.92,("CMS Preliminary, 19.2-19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
       banner.SetNDC(); banner.SetTextSize(0.028);
        
       leg = ROOT.TLegend(0.65,0.8,0.85,0.9);
       leg.SetFillStyle(0);
       leg.SetBorderSize(0);
       leg.SetNColumns(2);
     #       tmplabel = "observed, C' = %1.1f"%( float((11.-cprime[k])/10.) )
     #       leg.AddEntry(tGraphs[k*4+1],tmplabel,"L");        
       for k in range(nCprimes):
            print "k: ",k
            if k == 0: leg.AddEntry(tGraphs[0],"expected","L")
                
                
       drawColors = [];
       drawColors.append(ROOT.kRed+2);
       drawColors.append(ROOT.kGreen+2);
       drawColors.append(ROOT.kBlue+2);
       drawColors.append(ROOT.kMagenta+2);
       drawColors.append(ROOT.kCyan+2);
       drawColors.append(ROOT.kGray+2);                          
       drawColors.append(ROOT.kOrange+2);                                          
       drawColors.append(ROOT.kOrange+2);                                          
       drawColors.append(ROOT.kOrange+2);                                          
       drawColors.append(ROOT.kOrange+2);                                          
                
       oneLine = ROOT.TF1("oneLine","1",599,1001);
       oneLine.SetLineColor(ROOT.kRed);
       oneLine.SetLineWidth(3);

       can0 = ROOT.TCanvas("can0","can0",800,800);
       hrl0 = can0.DrawFrame(599,0.5,1001,20.0);
       hrl0.GetYaxis().SetTitle("#mu");
       hrl0.GetXaxis().SetTitle("mass (GeV)");
       can0.SetGrid();

       tGraphs[3].SetFillStyle(3001);
       tGraphs[2].SetFillStyle(3001);        
       tGraphs[3].Draw("f");
       tGraphs[2].Draw("f");
#      tGraphs[1].Draw("PL");
       tGraphs[0].Draw("PL");
       oneLine.Draw("LSAMES");
       leg.Draw();
#      ROOT.gPad.SetLogy();

       if not os.path.isdir("limitFigs"): 
          os.system("mkdir limitFigs" );
       if options.systematics == 1:        
          can0.SaveAs("Limit_with_syst.png");
       else:                               
          can0.SaveAs("Limit_no_syst.png");

       can = ROOT.TCanvas("can","can",800,800);
       hrl = can.DrawFrame(599,0.5,1001,100.0);
       hrl.GetYaxis().SetTitle("#mu' = C' #times #mu");
       hrl.GetXaxis().SetTitle("mass (GeV)");
       can.SetGrid();
       tGraphs[3].Draw("f");
       tGraphs[2].Draw("f");
       tGraphs[1].Draw("PL");
       tGraphs[0].Draw("PL");
            
       oneLine.Draw("LSAMES");

       for k in range(nCprimes):
            tGraphs[k*4].SetLineColor(drawColors[k]);
            tGraphs[k*4+1].SetLineColor(drawColors[k]);
            tGraphs[k*4].SetMarkerColor(drawColors[k]);
            tGraphs[k*4+1].SetMarkerColor(drawColors[k]);
            tGraphs[k*4].Draw("PL");
            tGraphs[k*4+1].Draw("PL");   
                         
       leg.Draw();

       ROOT.gPad.SetLogy();
       can.SaveAs("limitFigs/test_cPrime.eps");

       # ----------
        
       banner = TLatex(0.25,0.92,("CMS Preliminary, 19.2-19.3 fb^{-1} at #sqrt{s}=8TeV, e+#mu"));
       banner.SetNDC(); banner.SetTextSize(0.028);
        
       cprimeHist = [0.05,0.15,0.25,0.35,0.45,0.6,0.8,1.2];
       massHist = [550,650,750,850,950,1050];
       limitsHist_exp = ROOT.TH2F("limitsHist_exp","",nPoints,550,1050,nCprimes,0.05,1.05);
       limitsHist_obs = ROOT.TH2F("limitsHist_obs","",nPoints,550,1050,nCprimes,0.05,1.05);

       limitsHist_exp_600_2D = ROOT.TH2F("limitsHist_exp_600_2D","",nCprimes,0.05,1.05,nBRnews,-0.05,0.55);
       limitsHist_obs_600_2D = ROOT.TH2F("limitsHist_obs_600_2D","",nCprimes,0.05,1.05,nBRnews,-0.05,0.55); 
        
       for j in range(cpLo,cpHi):
            for i in range(mLo,mHi):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_2jet_%02d_00_unbin.Asymptotic.mH%03d.root"%(mass[i],cprime[j],mass[i]);
                curAsymLimits = getAsymLimits(curFile);
                limitsHist_exp.SetBinContent(i+1,j+1,curAsymLimits[3]);
                limitsHist_obs.SetBinContent(i+1,j+1,curAsymLimits[0]);          

       for j in range(cpLo,cpHi):
            for i in range(brLo,brHi):
                curFile = "higgsCombinehwwlvj_ggH%03d_em_2jet_%02d_%02d_unbin.Asymptotic.mH%03d.root"%(600,cprime[j],BRnew[i],600);
                curAsymLimits = getAsymLimits(curFile);
                limitsHist_exp_600_2D.SetBinContent(j+1,i+1,curAsymLimits[3]);
                limitsHist_obs_600_2D.SetBinContent(j+1,i+1,curAsymLimits[0]);          


       limitsHist_exp.SetTitle("; mass (GeV); C'^{2}; Expected Upper Limit, #mu'");
       limitsHist_obs.SetTitle("; mass (GeV); C'^{2}; Observed Upper Limit, #mu'");
       limitsHist_exp.GetZaxis().RotateTitle(1);
       limitsHist_obs.GetZaxis().RotateTitle(1);
    
       limitsHist_exp_600_2D.SetTitle("; C'^{2}; BR_{New}; Expected Upper Limit, #mu'");
       limitsHist_obs_600_2D.SetTitle("; C'^{2}; BR_{New}; Observed Upper Limit, #mu'");
       limitsHist_exp_600_2D.GetZaxis().RotateTitle(1);
       limitsHist_obs_600_2D.GetZaxis().RotateTitle(1);

       can2d_exp = ROOT.TCanvas("can2d_exp","can2d_exp",1000,800);
       limitsHist_exp.SetStats(0);
       limitsHist_exp.Draw("colz");
       banner.Draw();
       can2d_exp.SaveAs("limitFigs/test2D_exp.eps");

       can2d_obs = ROOT.TCanvas("can2d_obs","can2d_obs",1000,800);
       limitsHist_obs.SetStats(0);
       limitsHist_obs.Draw("colz");
       banner.Draw();
       can2d_obs.SaveAs("limitFigs/test2D_obs.eps");

       can2dbr_exp = ROOT.TCanvas("can2dbr_exp","can2dbr_exp",1000,800);
       limitsHist_exp_600_2D.SetStats(0);
       limitsHist_exp_600_2D.Draw("colz");
       banner.Draw();
       can2dbr_exp.SaveAs("limitFigs/test2Dbr_exp.eps");

       can2dbr_obs = ROOT.TCanvas("can2dbr_obs","can2dbr_obs",1000,800);
       limitsHist_obs_600_2D.SetStats(0);
       limitsHist_obs_600_2D.Draw("colz");
       banner.Draw();
       can2dbr_obs.SaveAs("limitFigs/test2Dbr_obs.eps");
