#!/bin/bash
#export XRD_NETWORKSTACK=IPv4
WD=$PWD

export EOSPATH=/eos/user/m/mhuwiler

RUNDIR="${9}"

CURRDIR="${10}"

mkdir -p "$RUNDIR"



export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch/
source $VO_CMS_SW_DIR/cmsset_default.sh
export SCRAM_ARCH=slc6_amd64_gcc630
scram project CMSSW CMSSW_10_2_6
cd CMSSW_10_2_6/src
eval $(scram runtime -sh)
git clone https://github.com/BoostedScalefactors/WTopScalefactorProducer.git
git clone https://github.com/cms-nanoAOD/nanoAOD-tools.git PhysicsTools/NanoAODTools
scram b -j2

cd $WD

export X509_USER_PROXY=${WD}/x509up #x509up_my

voms-proxy-info -all

#tar zxf /eos/home-m/mhuwiler/software/PhotonID/xgbosandbox.tgz
#cd xgbo
#pip install --user .

#cd $WD

#cp $CURRDIR/src/electron_id.py . 

python CMSSW_10_2_6/src/WTopScalefactorProducer/Skimmer/qsub/qsub_script_SFs.py $@ 


    





retval=$?

echo
echo
echo "Job finished with exit code $?"
echo "Files in ouput folder"


#cp $RUNDIR/results.txt ${CURRDIR}

#cp $RUNDIR/plots.root ${CURRDIR}

#cp -r $WD/electron_id ${CURRDIR} 


ls -ltr
if [[ $retval == 0 ]]; then
    errors=""
    for file in $(find -name '*.root' -or -name '*.xml'); do
        #cp -pv $file /eos/user/m/mhuwiler/data/PhotonId/RunSPVDProduction1correctedfull/output
        if [[ $? != 0 ]]; then
            errors="$errors $file($?)"
        fi
    done
    if [[ -n "$errors" ]]; then
       echo "Errors while staging files"
       echo "$errors"
       exit -2
    fi
fi


#exit $retval

