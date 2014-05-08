#python diffNuisances.py file_haddato.root -g NOME_PLOT -t 's'
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

parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
parser.add_option("--vtol", "--val-tolerance", dest="vtol", default=0.00, type="float", help="Report nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol", "--sig-tolerance", dest="stol", default=0.10, type="float", help="Report nuisances whose sigma changes by more than this amount")
parser.add_option("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type="float", help="Report severely nuisances whose value changes by more than this amount of sigmas")
parser.add_option("--stol2", "--sig-tolerance2", dest="stol2", default=0.50, type="float", help="Report severely nuisances whose sigma changes by more than this amount")
parser.add_option("-a", "--all",      dest="all",    default=False,  action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
parser.add_option("-A", "--absolute", dest="abs",    default=True,  action="store_true", help="Report also absolute values of nuisance values and errors, not only the ones normalized to the input sigma")
parser.add_option("-p", "--poi",      dest="poi",    default="r",    type="string",  help="Name of signal strength parameter (default is 'r' as per text2workspace.py)")
parser.add_option("-f", "--format",   dest="format", default="text", type="string",  help="Output format ('text', 'latex', 'twiki'")
parser.add_option("-g", "--histogram", dest="plotfile", default=None, type="string", help="If true, plot the pulls of the nuisances to the given file.")
parser.add_option("-t", "--type", dest="type", default='b', type="string", help="b=background only, s=s+b")

(options, args) = parser.parse_args()
if len(args) == 0:
    parser.print_usage()
    exit(1)

file = ROOT.TFile(args[0])

tree_s = file.Get("tree_fit_sb")
tree_b = file.Get("tree_fit_b")

prefit= file.Get("nuisances_prefit")

branches_list_b = tree_b.GetListOfBranches()
branches_list_s = tree_s.GetListOfBranches()

histo_list_In = [];
histo_list_fin= [];
list_mean_fin= [];
pulls_name= [];
errors=[];

pulls = []
sig_x = []


for i in range (branches_list_s.GetEntries()):

    hname = branches_list_s[i].GetName();    

    if ( hname.find("nll")==-1 and hname.find("n_exp")==-1 and hname.find("status")==-1 and hname.find("In")==-1 and hname.find("sigma")==-1):
      
       if hname.find("mu")==-1: 

            hname_sigma = branches_list_s[i+1].GetName();
            histo_temp = ROOT.TH1F(hname,"",100,-5,5);          
            histo_temp_sigma = ROOT.TH1F(hname_sigma,"",100,-5,5);    
            print hname;
            print hname_sigma;

            if options.type == 's':
                tree_s.Draw( branches_list_s[i].GetName()+" >> "+hname, "" ,"goff")
                tree_s.Draw( branches_list_s[i+1].GetName()+" >> "+hname_sigma, "" ,"goff")

            if options.type == 'b':
                tree_b.Draw( branches_list_s[i].GetName()+" >> "+hname, "" ,"goff")
                tree_b.Draw( branches_list_s[i+1].GetName()+" >> "+hname_sigma, "" ,"goff")            

            histo_list_fin.append(histo_temp);
            pulls_name.append(hname);
            pulls.append(float(histo_temp.GetMean()));
            errors.append(histo_temp_sigma.GetMean());

       else:

            histo_temp = ROOT.TH1F(hname,"",30,0,20);
            if options.type == 's':
                tree_s.Draw( branches_list_s[i].GetName()+" >> "+hname, "" ,"goff")
            if options.type == 'b':
                tree_b.Draw( branches_list_s[i].GetName()+" >> "+hname, "" ,"goff")
                
       canvas_temp = ROOT.TCanvas();
       histo_temp.Draw();
       canvas_temp.SaveAs(hname+".png","png");

if options.plotfile:
    import ROOT
#    ROOT.gROOT.SetStyle("Plain")
#    ROOT.gStyle.SetOptFit(1)
    histogram = ROOT.TH1F("","",len(pulls),0,len(pulls));
    histogram2 = ROOT.TH1F("a","a",len(pulls),0,len(pulls));    
    
    ROOT.gStyle.SetPadBottomMargin(0.50)
    ROOT.gStyle.SetPadLeftMargin  (0.05)
    ROOT.gStyle.SetPadRightMargin (0.05)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat       (   0)
    ROOT.gStyle.SetStatBorderSize(   0)
    ROOT.gStyle.SetStatColor     (  10)
    ROOT.gStyle.SetStatFont      (  42)
    ROOT.gStyle.SetStatX         (0.94)
    ROOT.gStyle.SetStatY         (0.91)
    ROOT.gStyle.cd()


                                            
    ROOT.gStyle.SetGridStyle (2)
    gr = ROOT.TGraphErrors();

#    print len(pulls);
    i=0;
    c=1;
    for pull in pulls:

#      if not pulls_name[i]=="mu":  
#        print pulls_name[i];
#        print errors[i];
        histogram.SetBinContent(i+1,pull);
        histogram.GetXaxis().SetBinLabel(i+1,pulls_name[i]);
        if pulls_name[i].find("Deco")==-1:
            histogram.GetXaxis().SetLabelSize(0.040);
        else:
            histogram.GetXaxis().SetLabelSize(0.035);        
#        histogram.GetXaxis().SetLabelFont(42);
        histogram.GetXaxis().LabelsOption("v");
        gr.SetPoint(i+1,i,0);
        gr.SetPointError(i+1,0,1);
        histogram.SetBinError(i+1,errors[i]);
        
#        histogram2.SetBinContent(i+1,errors[i]);
#        histogram2.GetXaxis().SetBinLabel(i+1,pulls_name[i]);
#        histogram2.GetXaxis().SetLabelSize(0.035);
#        histogram2.GetXaxis().SetLabelFont(42);
#        histogram2.GetXaxis().LabelsOption("v");

        i=i+1;

    gr.SetPoint(i+1,i,0);
    gr.SetPointError(i+1,0,1);

    canvas = ROOT.TCanvas("asdf", "asdf")
    histogram.GetYaxis().SetTitle("pull")
    histogram.GetYaxis().SetTitleOffset(0.55)    
    histogram.SetTitle("Post-fit nuisance pull distribution")
    histogram.SetMarkerStyle(20)
    histogram.SetMarkerSize(0.5)
    #histogram.Fit("gaus")
#    gr.SetLineColor(1);
    gr.SetFillColor(5);
    gr.SetFillStyle(3001);


    histogram.SetMaximum(3);
    histogram.SetMinimum(-3);    
    histogram.Draw("pe1")
#    gr.Draw("same");
    gr.Draw("e3same");    
    histogram.Draw("pe1same")    
   # ROOT.gPad.RedrawAxis();
#    histogram.Draw("pe1")    

    canvas.SetGridy();
    canvas.SaveAs(options.plotfile+".png","png")


    '''
    canvas2 = ROOT.TCanvas("asdf2", "asdf2")
    histogram2.GetYaxis().SetTitle("(k-1)/RMS")
    histogram2.GetYaxis().SetTitleOffset(0.55)    
    histogram2.SetTitle("(k-1)/RMS")
    histogram2.SetMarkerStyle(20)
    histogram2.SetMarkerSize(0.5)

    histogram2.SetMaximum(4);
    histogram2.SetMinimum(-0.1);    
    histogram2.Draw("p")

    canvas2.SetGridy();
    canvas2.SaveAs(options.plotfile+"_k.png","png")
    '''
