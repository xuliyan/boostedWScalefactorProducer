#########################################
### How to run the new runLimits code ###
#########################################

1) Options to be used -> complete list

   General Options:

   --makeCards : in order to run the fitting analysis and produce workspaces and datacards
   --computeLimits: skip the cards production and go to the combine analysis 
   --plotLimits : plot asymptotic limits and p-values 
   --biasStudy  : in order to run bias study for mJ fit, bias study for mWW in the ttbar control region for the ttbar, in the sb of mWW for the w+jets   
   --maximumLikelihoodFit : perform maximumLikelihood fits inside combination tool
   --generateOnly : generate pseudodata toys according to a datacard model

   
   Submission options:

   --batchMode : to submit job on condor cluster
   --lxbatchCern : if also batchMode is specified, submit jobs on lxbatch system at Cern
   --herculesMilano : if also batchMode is specified, submit jobs on hercules system at Milano


   Basic Options:

   --datacardDIR : specify the directory with the datacards to be used .. the directory will set as working dir (os.chdir())
   --channel     : specify the channel label "em" for the merged category, "mu" for muons and "el" for electrons
   --systematics : higgs mass point to be considered
   --cPrime      : c' for the bsm model to be considered
   --brNew       : BR new for the bsm model
   --odir        : output directory for cards
   --sigChannel  : label for signal: nothing ggH + vbfH, otherwise only ggH or vbfH cards
   --jetBin      : define the jet bin, in order to use the correct datacards
   --turnOnAnalysis : if 0 use Exp or Pow by default, otherwise ErfExp or ErfPow --> shapes are defined in static vectors inside the runLimtisCode


   MakeCards specific options:

   --pseudodata  : run the cards analysis on data or pseudodata

   Asymptotic Limit specific options: 

   --systematics : run the asympotic limit without systematics applied in case it is set to 1

   MaximumLikelihood Fit and Generation only options:

   --injectSingalStrengh : amount of signal to be injected wrt what present in the datacard
   --nToys               : numer of toys -t N of combine
   --crossedToys         : if 1 it is saying that you want to use a different fitting model wrt to the generation

   --inputGeneratedDataset :  if you use this options when --generateOnly is true, it tells the directory where to store the generated dataset
                              if you use this options when --maximumLikelihoodFi is true, it tells the directory location where to take the dataset for the fit

   --outputTree : if this is equal to zero, it means that it will create as many jobs as the nToys number
                  if this is equal to one, it means that it will create one single job with a tree or a list of dataset in the output file

   
   Our own bias study options :
  
   --injectSingalStrenght : amount of signal to be injected
   --shapetes             : perform the f-test
   --ttbarcontrolregion   : analyze or not the ttbar control region
   --mlvjregion           : _sb_lo for the lower sideband, _sb_hi for the higher one, _signal_region for the signal region where generate toys
   --fitjetmass           : if 0 look at the mWW in the --mlvjregion, otherwise look at the mJ distribution
   --onlybackgroundfit    : if 1 only B fit in the toys, if 0 do S+B fits
   --inflatejobstatistic  : generate N times the expected yields of event accorting to the prefit model done on the expected yields
   --scalesignalwidth     : change the with of the signal as -> width/scalesignalwidth 



2) Pratical examples for the Bias study:

   -> generate toys with a specific model without signal, foe a mass point :

      python runLimits.py --computeLimits --generateOnly --batchMode --herculesMilano 
                          --datacardDIR <nome datacard> --channel em --massPoint 600 --jetBin _2jet --injectSingalStrenght 0
                           --nToys 1000 --inputGeneratedDataset Generated_SB --outputTree 0
       
     in this case it will generate nToys files, in each just one datasetm in the Generated_SB directory under (inside) the --datacardDIR directory


   -> Make Fit and generation together (not crossed study --> exampl generate with Exp and re-fit with Exp for signal injection test):

      python runLimits.py --computeLimits --maximumLikelihood --batchMode --herculesMilano 
                          --datacardDIR <nome datacard> --channel em --massPoint 600 --jetBin _2jet --injectSingalStrenght 5
                           --nToys 1000 --outputTree 0 --crossedToys 0


   -> Make the Fit on a set of generated dataset, using a different model between generation and fit 

      First of all you must generate toys according with a defined model as explained above.

      Then you have to:

      python runLimits.py --computeLimits --maximumLikelihood --batchMode --herculesMilano --datacardDIR <nome datacard>
                          --channel em --massPoint 600 --jetBin _2jet  --nToys 1000 --inputGeneratedDataset ../<dacard_path>/Generated_SB --outputTree 1 --crossedToys 1



      