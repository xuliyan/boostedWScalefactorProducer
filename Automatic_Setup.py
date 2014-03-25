
#! /usr/bin/env python
import os
import glob
import math
import array
import ROOT
import ntpath
import sys
import subprocess

from optparse import OptionParser
from subprocess import Popen
from ROOT import gROOT, gStyle, gSystem, TLatex



if __name__ == "__main__":

  ROOT.gSystem.AddIncludePath("-I$ROOFITSYS/include");

  inPath = os.getenv("PWD")

  os.chdir(inPath+"/PlotStyle");
  ROOT.gROOT.ProcessLine(".L Util.cxx+");
  ROOT.gSystem.Load("Util_cxx.so");

  ROOT.gROOT.ProcessLine(".L PlotUtils.cxx+");
  ROOT.gSystem.Load("PlotUtils_cxx.so");

 
  os.chdir(inPath+"/PDFs");
  ROOT.gROOT.ProcessLine(".L PdfDiagonalizer.cc+");
  ROOT.gSystem.Load("PdfDiagonalizer_cc.so");

  ROOT.gROOT.ProcessLine(".L PdfDiagonalizer.cc+");
  ROOT.gSystem.Load("PdfDiagonalizer_cc.so");

  ROOT.gROOT.ProcessLine(".L HWWLVJRooPdfs.cxx+");
  ROOT.gSystem.Load("HWWLVJRooPdfs_cxx.so");

  ROOT.gROOT.ProcessLine(".L MakePdf.cxx+");
  ROOT.gSystem.Load("MakePdf_cxx.so");


  os.chdir(inPath+"/BiasStudy");
  ROOT.gROOT.ProcessLine(".L BiasUtils.cxx+");
  ROOT.gSystem.Load("BiasUtils_cxx.so");

  os.chdir(inPath+"/FitUtils");
  ROOT.gROOT.ProcessLine(".L FitUtils.cxx+");
  ROOT.gSystem.Load("FitUtils_cxx.so");
  
  
