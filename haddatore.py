#! /usr/bin/env python
import os
import glob
import math
from array import array
import sys
import time
import subprocess
import ROOT

mass=[600,700,800,900,1000]
cprime=[01,02,03,05,07,10]
BRnew=[00,01,02,03,04,05]

for i in range(len(mass)):
    for j in range(len(cprime)):
        for k in range(len(BRnew)):

            command="hadd -f cards_combo/higgsCombinehwwlvj_pval_exp_ggH%03d_combo_%02d_%02d_unbin.ProfileLikelihood.mH%03d.root cards_combo/higgsCombinehwwlvj_pval_exp_ggH%03d_combo_%02d_%02d_unbin_*"%(mass[i],cprime[j],BRnew[k],mass[i],mass[i],cprime[j],BRnew[k]);
            os.system(command);

