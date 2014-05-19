#!/usr/bin/env python
import re
from sys import argv, stdout, stderr, exit
from optparse import OptionParser

# import ROOT with a fix to get batch mode (http://root.cern.ch/phpBB3/viewtopic.php?t=3198)
hasHelp = False
for X in ("-h", "-?", "--help"):
 if X in argv:
    hasHelp = True
    argv.remove(X)
argv.append( '-b-' )

import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
argv.remove( '-b-' )

if hasHelp: argv.append("-h")

## parsing inputs 

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--vtol", "--val-tolerance", dest="vtol", default=0.30, type="float", help="Report nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol", "--sig-tolerance", dest="stol", default=0.10, type="float", help="Report nuisances whose sigma changes by more than this amount")
parser.add_option("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type="float", help="Report severely nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol2", "--sig-tolerance2", dest="stol2", default=0.50, type="float", help="Report severely nuisances whose sigma changes by more than this amount")
parser.add_option("-a", "--all",      dest="all", default=False, action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
parser.add_option("-A", "--absolute", dest="abs", default=False,  action="store_true", help="Report also absolute values of nuisance values and errors, not only the ones normalized to the input sigma")
parser.add_option("-p", "--poi",      dest="poi",    default="r",    type="string",  help="Name of signal strength parameter (default is 'r' as per text2workspace.py)")
parser.add_option("-f", "--format",   dest="format", default="text", type="string",  help="Output format ('text', 'latex', 'twiki'")
parser.add_option("-g", "--histogram", dest="plotfile", default=None, type="string", help="If true, plot the pulls of the nuisances to the given file.")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

## open the input file 
file = ROOT.TFile(args[0])

nuis_histo_list_s  = [];
pulls_s            = [];
pulls_name_s       = [];
errors_s           = [];

nuis_histo_list_b  = [];
pulls_b            = [];
pulls_name_b       = [];
errors_b           = [];

isFlagged = {};
table     = {};

prefit = file.Get("nuisances_prefit")
fit_s  = file.Get("fit_s")
fit_b  = file.Get("fit_b")

tree_s = file.Get("tree_fit_sb")
tree_b = file.Get("tree_fit_b")

if tree_s == None or tree_s.ClassName() != "TTree": raise RuntimeError, "File %s does not contain the output of the signal fit tree 'fit_s'" % args[0]
if tree_b == None or tree_b.ClassName() != "TTree": raise RuntimeError, "File %s does not contain the output of the background fit tree 'fit_b'" % args[0]
if prefit == None or prefit.ClassName() != "RooArgSet": raise RuntimeError, "File %s does not contain the prefit nuisances 'nuisances_prefit'" % args[0]

branches_list_b = tree_b.GetListOfBranches()
branches_list_s = tree_s.GetListOfBranches()

fpf_s = fit_s.floatParsFinal()
fpf_b = fit_b.floatParsFinal()

### Loop on signal tree

for i in range (branches_list_s.GetEntries()):
  hname_s = branches_list_s[i].GetName();    
  row    = [];
  flag   = False;
  nuis_p = prefit.find(hname_s);

  if ( hname_s.find("nll")==-1 and hname_s.find("n_exp")==-1 and hname_s.find("status")==-1 and hname_s.find("In")==-1 and hname_s.find("sigma")==-1):

    if nuis_p == None:
      if not options.abs: continue      
    else: 
      mean_p, sigma_p = (nuis_p.getVal(), nuis_p.getError())
      if options.abs: row += [ "%.2f +/- %.2f" % (nuis_p.getVal(), nuis_p.getError()) ]

    for fit_name, nuis_x in [('b',(branches_list_s.FindObject(hname_s)).GetName()),('s',hname_s)]:

     if fit_name == 's' and fpf_s.find(nuis_x) == None:
        row += [ " n/a " ]
     elif fit_name == 'b' and fpf_b.find(nuis_x) == None:
        row += [ " n/a " ]
     elif fit_name == 's' and fpf_s.find(nuis_x) != None:                 

        histo_nuis_s       = ROOT.TH1F(nuis_x+"_s","",100,-5,5);          
        histo_nuis_sigma_s = ROOT.TH1F(nuis_x+"_sigma_s","",100,-5,5);    
        tree_s.Draw(nuis_x+" >> "+nuis_x+"_s", "" ,"goff");

        if nuis_x != "mu" :
          if branches_list_s.FindObject("sigma_from_fit_"+nuis_x):   
           nuis_x_sigma_s     = (branches_list_s.FindObject("sigma_from_fit_"+nuis_x)).GetName();
           tree_s.Draw(nuis_x_sigma_s+" >> "+histo_nuis_sigma_s.GetName(), "" ,"goff");
          else:
            nuis_x_sigma_s =  fpf_s.find(nuis_x);
            histo_nuis_sigma_s.Fill(nuis_x_sigma_s.getError());
                             
          pulls_name_s.append(nuis_x);            
          row += [ "%+.2f +/- %.2f" % (histo_nuis_s.GetMean(), histo_nuis_sigma_s.GetMean())];
          pulls_s.append(float(histo_nuis_s.GetMean()-nuis_p.getVal())/histo_nuis_sigma_s.GetMean());
          errors_s.append(histo_nuis_sigma_s.GetMean()/nuis_p.getError());

          if options.abs:
           row[-1] += " (%+4.2fsig, %4.2f)" % (pulls_s[len(pulls_s)-1],errors_s[len(errors_s)-1]);
          else:
           row[-1] = " %+4.2f, %4.2f" % (pulls_s[len(pulls_s)-1],errors_s[len(errors_s)-1]);

          if(abs(pulls_s[len(pulls_s)-1]) > options.vtol2 or abs(errors_s[len(errors_s)-1]-1) > options.stol2):
            isFlagged[(nuis_p.GetName(),fit_name)] = 2;
            flag = True;
          elif(abs(pulls_s[len(pulls_s)-1]) > options.vtol or abs(errors_s[len(errors_s)-1]-1) > options.stol):
            if options.all: isFlagged[(nuis_p.GetName(),fit_name)] = 1
            flag = True
          elif options.all:
            flag = True
             
     elif fit_name == 'b' and fpf_b.find(nuis_x) != None:

        histo_nuis_b       = ROOT.TH1F(nuis_x+"_b","",100,-5,5);          
        histo_nuis_sigma_b = ROOT.TH1F(nuis_x+"_sigma_b","",100,-5,5);    
        tree_b.Draw(nuis_x+" >> "+nuis_x+"_b", "" ,"goff");
    
        if nuis_x != "mu":
            
          if branches_list_b.FindObject("sigma_from_fit_"+nuis_x) :       
             nuis_x_sigma_b     = (branches_list_b.FindObject("sigma_from_fit_"+nuis_x)).GetName();
             tree_b.Draw(nuis_x_sigma_b+" >> "+histo_nuis_sigma_b.GetName(), "" ,"goff");
          else:    
             nuis_x_sigma_b =  fpf_b.find(nuis_x);
             histo_nuis_sigma_b.Fill(nuis_x_sigma_b.getError());

          nuis_histo_list_b.append(histo_nuis_b);

          pulls_name_b.append(nuis_x);            
          row += [ "%+.2f +/- %.2f" % (histo_nuis_b.GetMean(),histo_nuis_sigma_b.GetMean())];
          pulls_b.append(float(histo_nuis_b.GetMean()-nuis_p.getVal())/histo_nuis_sigma_b.GetMean());
          errors_b.append(histo_nuis_sigma_b.GetMean()/nuis_p.getError());
          if options.abs:
           row[-1] += " (%+4.2fsig, %4.2f)" % (pulls_b[len(pulls_b)-1],errors_b[len(errors_b)-1]);
          else:
           row[-1] = " %+4.2f, %4.2f" % (pulls_b[len(pulls_b)-1],errors_b[len(errors_b)-1]);
          if(abs(pulls_b[len(pulls_b)-1]) > options.vtol2 or abs(errors_b[len(errors_b)-1]-1) > options.stol2):
            isFlagged[(nuis_p.GetName(),fit_name)] = 2;
            flag = True;
          elif(abs(pulls_b[len(pulls_b)-1]) > options.vtol or abs(errors_b[len(errors_b)-1]-1) > options.stol):
           if options.all: isFlagged[(nuis_p.GetName(),fit_name)] = 1
           flag = True;
          elif options.all:
           flag = True;
    row += [ "%+4.2f"  % fit_s.correlation(nuis_p.GetName(), options.poi) ]            
    if flag or options.all:
        table[nuis_p.GetName()] = row;

### print table
fmtstring = "%-40s %15s %15s %10s"
highlight = "*%s*"
morelight = "!%s!"
pmsub, sigsub = None, None
if options.format == 'text':
    if options.abs:
         fmtstring = "%-40s %15s %30s %30s %10s"
         print fmtstring % ('name', 'pre fit', 'b-only fit', 's+b fit', 'rho')
    else:
         print fmtstring % ('name', 'b-only fit', 's+b fit', 'rho')
elif options.format == 'latex':
  pmsub = (r"(\S+) \+/- (\S+)", r"$\1 \\pm \2$")
  sigsub = ("sig", r"$\\sigma$")
  highlight = "\\textbf{%s}"
  morelight = "{{\\color{red}\\textbf{%s}}}"
  if options.abs:
    fmtstring = "%-40s & %15s & %30s & %30s & %6s \\\\"
    print "\\begin{tabular}{|l|r|r|r|r|} \\hline ";
    print (fmtstring % ('name', 'pre fit', '$b$-only fit', '$s+b$ fit', r'$\rho(\theta, \mu)$')), " \\hline"
  else:
    fmtstring = "%-40s & %15s & %15s & %6s \\\\"
    print "\\begin{tabular}{|l|r|r|r|} \\hline ";
    what = r"\Delta x/\sigma_{\text{in}}$, $\sigma_{\text{out}}/\sigma_{\text{in}}$"
    print fmtstring % ('', '$b$-only fit', '$s+b$ fit', '')
    print (fmtstring % ('name', what, what, r'$\rho(\theta, \mu)$')), " \\hline"
elif options.format == 'twiki':
  pmsub = (r"(\S+) \+/- (\S+)", r"\1 &plusmn; \2")
  sigsub = ("sig", r"&sigma;")
  highlight = "<b>%s</b>"
  morelight = "<b style='color:red;'>%s</b>"
  if options.abs:
    fmtstring = "| <verbatim>%-40s</verbatim> | %-15s | %-30s | %-30s | %-15s |"
    print "| *name* | *pre fit* | *b-only fit* | *s+b fit* | "
  else:
    fmtstring = "| <verbatim>%-40s</verbatim> | %-15s | %-15s | %-15s |"
    print "| *name* | *b-only fit* | *s+b fit* | *corr.* |"

elif options.format == 'html':
  pmsub = (r"(\S+) \+/- (\S+)", r"\1 &plusmn; \2")
  sigsub = ("sig", r"&sigma;")
  highlight = "<b>%s</b>"
  morelight = "<strong>%s</strong>"
  print """
<html><head><title>Comparison of nuisances</title>
<style type="text/css">
  td, th { border-bottom: 1px solid black; padding: 1px 1em; }
  td { font-family: 'Consolas', 'Courier New', courier, monospace; }
  strong { color: red; font-weight: bolder; }
</style>
</head><body style="font-family: 'Verdana', sans-serif; font-size: 10pt;"><h1>Comparison of nuisances</h1>
<table>
"""
  if options.abs:
    print "<tr><th>nuisance</th><th>pre fit</th><th>background fit </th><th>signal fit</th><th>correlation</th></tr>"
    fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td><td> %-30s </td><td> %-30s </td><td> %-15s </td></tr>"
  else:
    what = "&Delta;x/&sigma;<sub>in</sub>, &sigma;<sub>out</sub>/&sigma;<sub>in</sub>";
    print "<tr><th>nuisance</th><th>background fit<br/>%s </th><th>signal fit<br/>%s</th><th>&rho;(&mu;, &theta;)</tr>" % (what,what)
    fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td><td> %-15s </td><td> %-15s </td></tr>"

names = table.keys()
names.sort()
highlighters = { 1:highlight, 2:morelight };

for n in names:
 v = table[n]
 if options.format == "latex": n = n.replace(r"_", r"\_")
 if pmsub != None: v = [ re.sub(pmsub[0], pmsub[1], i) for i in v ]
 if sigsub != None: v = [ re.sub(sigsub[0], sigsub[1], i) for i in v ]
 if (n,'b') in isFlagged: v[-3] = highlighters[isFlagged[(n,'b')]] % v[-3]
 if (n,'s') in isFlagged: v[-2] = highlighters[isFlagged[(n,'s')]] % v[-2]
 if options.abs:
   print fmtstring % (n, v[0], v[1], v[2], v[3])
 else:
   print fmtstring % (n, v[0], v[1], v[2])

if options.format == "latex":
 print " \\hline\n\end{tabular}"
elif options.format == "html":
 print "</table></body></html>"
                                                                                              
if options.plotfile:

 import ROOT
 histo_pull_s_fit = ROOT.TH1F("histo_pull_s_fit","pull S+B fit",len(pulls_s),0,len(pulls_s));
 histo_pull_b_fit = ROOT.TH1F("histo_pull_b_fit","pull B fit",len(pulls_b),0,len(pulls_b));
    
 ROOT.gStyle.SetPadBottomMargin(0.50)
 ROOT.gStyle.SetPadLeftMargin  (0.05)
 ROOT.gStyle.SetPadRightMargin (0.05)
 ROOT.gStyle.SetOptStat(0)
 ROOT.gStyle.SetStatBorderSize(0)
 ROOT.gStyle.SetStatColor(10)
 ROOT.gStyle.SetStatFont(42)
 ROOT.gStyle.SetStatX(0.94)
 ROOT.gStyle.SetStatY(0.91)
 ROOT.gStyle.cd()

                                            
 ROOT.gStyle.SetGridStyle(2)
 gr_s_prefit = ROOT.TGraphErrors(len(pulls_b));
 gr_s_prefit.SetName("gr_s_prefit");
 gr_b_prefit = ROOT.TGraphErrors(len(pulls_b));
 gr_b_prefit.SetName("gr_b_prefit");
 
 for ipull in range(len(pulls_b)):
  histo_pull_b_fit.SetBinContent(ipull+1,pulls_b[ipull]);
  histo_pull_b_fit.SetBinError(ipull+1,errors_b[ipull]);        
  histo_pull_b_fit.GetXaxis().SetBinLabel(ipull+1,pulls_name_b[ipull]);
  gr_b_prefit.SetPoint(ipull,ipull,0);
  gr_b_prefit.SetPointError(ipull,0,1);
 
 gr_b_prefit.SetPoint(len(pulls_b),len(pulls_b),0);
 gr_b_prefit.SetPointError(len(pulls_b),0,1);

 histo_pull_b_fit.GetXaxis().SetLabelSize(0.035);        
 histo_pull_b_fit.GetXaxis().LabelsOption("v");

 canvas = ROOT.TCanvas("Pulls B fit", "Pulls B fit")
 histo_pull_b_fit.GetYaxis().SetTitle("pull")
 histo_pull_b_fit.GetYaxis().SetTitleOffset(0.55)    
 histo_pull_b_fit.SetMarkerStyle(20)
 histo_pull_b_fit.SetMarkerSize(0.5)
 gr_b_prefit.SetFillColor(5);
 gr_b_prefit.SetFillStyle(3001);
 canvas.SetGridy();


 histo_pull_b_fit.SetMaximum(3);
 histo_pull_b_fit.SetMinimum(-3);    
 histo_pull_b_fit.SetLineColor(1);    
 histo_pull_b_fit.SetMarkerColor(1);    
 histo_pull_b_fit.Draw("pe1")
 gr_b_prefit.Draw("e3same");

 line_up = ROOT.TLine(0,1,len(pulls_b),1);
 line_dn = ROOT.TLine(0,-1,len(pulls_b),-1);
 line_up.SetLineColor(ROOT.kRed);
 line_dn.SetLineColor(ROOT.kRed);
 line_up.SetLineWidth(2);
 line_dn.SetLineWidth(2);
 line_up.Draw("same") 
 line_dn.Draw("same") 

 histo_pull_b_fit.Draw("pe1same")    

 canvas.RedrawAxis();
 canvas.RedrawAxis("g");
 canvas.Update();
        
 canvas.SaveAs(options.plotfile+"_B.png","png")
 canvas.SaveAs(options.plotfile+"_B.pdf","pdf")

 for ipull in range(len(pulls_s)):

  histo_pull_s_fit.SetBinContent(ipull+1,pulls_s[ipull]);
  histo_pull_s_fit.SetBinError(ipull+1,errors_s[ipull]);        
  histo_pull_s_fit.GetXaxis().SetBinLabel(ipull+1,pulls_name_s[ipull]);
  gr_s_prefit.SetPoint(ipull,ipull,0);
  gr_s_prefit.SetPointError(ipull,0,1);

 gr_s_prefit.SetPoint(len(pulls_s),len(pulls_s),0);
 gr_s_prefit.SetPointError(len(pulls_s),0,1);

 histo_pull_s_fit.GetXaxis().SetLabelSize(0.035);        
 histo_pull_s_fit.GetXaxis().LabelsOption("v");

 canvas = ROOT.TCanvas("Pulls S+B fit", "Pulls S+B fit")
 histo_pull_s_fit.GetYaxis().SetTitle("pull")
 histo_pull_s_fit.GetYaxis().SetTitleOffset(0.55)    
 histo_pull_s_fit.SetMarkerStyle(20)
 histo_pull_s_fit.SetMarkerSize(0.5)
 gr_s_prefit.SetFillColor(5);
 gr_s_prefit.SetFillStyle(3001);


 canvas.SetGridy();
 histo_pull_s_fit.SetLineColor(1);    
 histo_pull_s_fit.SetMarkerColor(1);    
 histo_pull_s_fit.SetMaximum(3);
 histo_pull_s_fit.SetMinimum(-3);    
 histo_pull_s_fit.Draw("pe1")
 gr_s_prefit.Draw("3same");    
 line_up.Draw("same") 
 line_dn.Draw("same") 
 histo_pull_s_fit.Draw("pe1same")    

 canvas.RedrawAxis();
 canvas.RedrawAxis("g");
 canvas.Update();

 canvas.SaveAs(options.plotfile+"_SB.png","png")
 canvas.SaveAs(options.plotfile+"_SB.pdf","pdf")

