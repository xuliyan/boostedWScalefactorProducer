def getXsec(sample):
  if sample.find("QCD_HT100to200"                       ) !=-1 : return 23700000.0;
  elif sample.find("QCD_HT200to300"                       ) !=-1 : return 1547000.0;
  elif sample.find("QCD_HT300to500"                       ) !=-1 : return 322600.0;
  elif sample.find("QCD_HT500to700"                       ) !=-1 : return 29980.0;
  elif sample.find("QCD_HT700to1000"                      ) !=-1 : return 6334.0;
  elif sample.find("QCD_HT1000to1500"                     ) !=-1 : return 1088.0;
  elif sample.find("QCD_HT1500to2000"                     ) !=-1 : return 99.11;
  elif sample.find("QCD_HT2000toInf"                      ) !=-1 : return 20.23;
  elif sample.find("ST_s-channel_4f_leptonDecays"          ) !=-1 : return 3.74*0.3272 ;
  elif sample.find("ST_t-channel_top_4f_inclusiveDecays"   ) !=-1 : return 113.3        ;
  elif sample.find("ST_t-channel_antitop_4f_inclusiveDecays") !=-1 : return 67.91       ;
  elif sample.find("ST_tW_antitop_5f_inclusiveDecays"      ) !=-1 : return 34.97         ;
  elif sample.find("ST_tW_top_5f_inclusiveDecays"          ) !=-1 : return 34.91         ;
  elif sample.find("TT_TuneCH3_13TeV-powheg-herwig7"      ) !=-1 : return 722.8	;
  elif sample.find("TTJets_TuneCP5_13TeV-amcatnloFXFX"    ) !=-1 : return 722.8;
  elif sample.find("TTToSemiLeptonic_TuneCP5"  ) !=-1 : return 687.1*0.3272*0.6741*2.;
  elif sample.find("TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8") !=-1 : return 687.1*0.3272*0.3272;
  elif sample.find("WJetsToLNu_HT-70To100"                 ) !=-1 : return 1292.0*1.22    ;
  elif sample.find("WJetsToLNu_HT-100To200"                ) !=-1 : return 1395.0*1.22    ;
  elif sample.find("WJetsToLNu_HT-200To400"                ) !=-1 : return 407.9*1.22     ;
  elif sample.find("WJetsToLNu_HT-400To600"                ) !=-1 : return 57.48*1.22    ;
  elif sample.find("WJetsToLNu_HT-600To800"                ) !=-1 : return 12.87*1.22   ;
  elif sample.find("WJetsToLNu_HT-800To1200"               ) !=-1 : return 5.366*1.22    ;
  elif sample.find("WJetsToLNu_HT-1200To2500"              ) !=-1 : return 1.074*1.22    ;
  elif sample.find("WJetsToLNu_HT-2500ToInf"               ) !=-1 : return 0.008001*1.22 ;
  elif sample.find("WJetsToLNu_TuneCP5"                    ) !=-1 : return 52940.0;
  elif sample.find("WW_Tune"                               ) !=-1 : return 75.8        ;
  elif sample.find("WZ_Tune"                               ) !=-1 : return 27.6        ;
  elif sample.find("ZZ_Tune"                               ) !=-1 : return 12.14         ;
  elif sample.find("SingleMuon")!=-1  or sample.find("SingleElectron") !=-1 or sample.find("JetHT") !=-1 or sample.find("data") !=-1 : return 1.
  else:
	  print "Cross section not defined for this sample!!"
	  return 0

def getNGenerated(sample):
  if sample.find("QCD_HT100to200"                       ) !=-1 : return 0;
  elif sample.find("QCD_HT200to300"                       ) !=-1 : return 0;
  elif sample.find("QCD_HT300to500"                       ) !=-1 : return 0;
  elif sample.find("QCD_HT500to700"                       ) !=-1 : return 0;
  elif sample.find("QCD_HT700to1000"                      ) !=-1 : return 0;
  elif sample.find("QCD_HT1000to1500"                     ) !=-1 : return 0;
  elif sample.find("QCD_HT1500to2000"                     ) !=-1 : return 0;
  elif sample.find("QCD_HT2000toInf"                      ) !=-1 : return 0;
  elif sample.find("ST_s-channel_4f_leptonDecays"          ) !=-1 : return 9883805;
  elif sample.find("ST_t-channel_top_4f_inclusiveDecays"   ) !=-1 : return 5982064;
  elif sample.find("ST_t-channel_antitop_4f_inclusiveDecays") !=-1 : return 3675910;
  elif sample.find("ST_tW_antitop_5f_inclusiveDecays"      ) !=-1 : return 7977430;
  elif sample.find("ST_tW_top_5f_inclusiveDecays"          ) !=-1 : return 7794186;
  elif sample.find("TT_TuneCH3_13TeV-powheg-herwig7"      ) !=-1 : return 97531483;
  elif sample.find("TTJets_TuneCP5_13TeV-amcatnloFXFX"    ) !=-1 : return 153896982;
  elif sample.find("TTToSemiLeptonic"                     ) !=-1 : return 197670000;
  elif sample.find("TTTo2L2Nu"                            ) !=-1 : return 9000000;
  elif sample.find("WJetsToLNu_HT-70To100"                 ) !=-1 : return 22255124;
  elif sample.find("WJetsToLNu_HT-100To200"                ) !=-1 : return 35862893;
  elif sample.find("WJetsToLNu_HT-200To400"                ) !=-1 : return 21250517;
  elif sample.find("WJetsToLNu_HT-400To600"                ) !=-1 : return 14313274;
  elif sample.find("WJetsToLNu_HT-600To800"                ) !=-1 : return 21709087;
  elif sample.find("WJetsToLNu_HT-800To1200"               ) !=-1 : return 20432728;
  elif sample.find("WJetsToLNu_HT-1200To2500"              ) !=-1 : return 20258624;
  elif sample.find("WJetsToLNu_HT-2500ToInf"               ) !=-1 : return 21495421;
  elif sample.find("WJetsToLNu_TuneCP5"                    ) !=-1 : return 30008250;
  elif sample.find("WW_Tune"                               ) !=-1 : return 7765828;
  elif sample.find("WZ_Tune"                               ) !=-1 : return 3807850;
  elif sample.find("ZZ_Tune"                               ) !=-1 : return 1949768;
  elif sample.find("SingleMuon")!=-1  or sample.find("SingleElectron") !=-1 or sample.find("JetHT") !=-1 or sample.find("data") !=-1 : return 1.
  else:
	  print "N events not defined for this sample!!"
	  return 0

# def getSF(sample):
 # if sample.find("WJetsToLNu_HT-70To100"                 ) !=-1 : return 1.26    ;
 # elif sample.find("WJetsToLNu_HT-100To200"                ) !=-1 : return 1.26    ;
 # elif sample.find("WJetsToLNu_HT-200To400"                ) !=-1 : return 1.48     ;
 # elif sample.find("WJetsToLNu_HT-400To600"                ) !=-1 : return 1.26    ;
 # elif sample.find("WJetsToLNu_HT-600To800"                ) !=-1 : return 1.03   ;
 # elif sample.find("WJetsToLNu_HT-800To1200"               ) !=-1 : return 1.05    ;
 # elif sample.find("WJetsToLNu_HT-1200To2500"              ) !=-1 : return 0.77    ;
 # elif sample.find("WJetsToLNu_HT-2500ToInf"               ) !=-1 : return 0.77 ;
 # else:
	  # return 1.

def getLumi(sample):
  return 41530.
