import os
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-m', '--mass', action = "store",     type = "string", dest = "mass", default = "600")

(options, args) = parser.parse_args()

mass=options.mass;

for i in range(0,980):
  outScript = open("Jobtemp/bjob_"+mass+"_"+str(i)+".sh","w");
  outScript.write('#!/bin/bash');
  outScript.write("\n"+'cd /gwpool/users/brianza/PHD/CMSSW_6_1_1/');
  outScript.write("\n"+'eval `scramv1 runtime -sh`');
  outScript.write("\n"+'export SCRAM_ARCH=slc5_amd64_gcc472');
  outScript.write("\n"+'cd /gwpool/users/brianza/PHD/CMSSW_6_1_1/src/boostedWWAnalysis_SAVENORM/');
  outScript.write("\n"+'python runLimits600.py --channel EXP_TIGHT_NARROW --number '+str(i)+' --massPoint '+mass+' -b');
  outScript.close();

outScript2 = open("lanciatutto.sh","w");
for i in range(0,980):
  outScript2.write("\n"+'qsub -V -d /gwpool/users/brianza/PHD/CMSSW_6_1_1/src/boostedWWAnalysis_SAVENORM/Jobtemp -q longcms /gwpool/users/brianza/PHD/CMSSW_6_1_1/src/boostedWWAnalysis_SAVENORM/Jobtemp/bjob_'+mass+'_'+str(i)+'.sh');
outScript2.close();

#os.system("sh lanciatutto.sh");
#os.system("rm lanciatutto.sh");
