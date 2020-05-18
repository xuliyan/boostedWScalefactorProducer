#!/bin/bash


cd ${_CONDOR_SCRATCH_DIR}

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh
if [ -r CMSSW_10_2_6/src ] ; then 
 echo release CMSSW_10_2_6 already exists
else
scram p CMSSW CMSSW_10_2_6
fi
cd CMSSW_10_2_6/src
eval `scram runtime -sh`
echo pwd

# Untar transfered files
mv ../../transfer.tgz ./
tar -xvf transfer.tgz

scram b

cd boostedWScalefactorProducer/Skimmer/condorsub


# Use condor process ID to select block of files to process
# $2 means the $2th block containing 2 files each
n_file_start=$(($2*1))
n_file_end=$((($2+1)*1))
echo $n_file_start
echo $n_file_end
n_file_current=0
while IFS= read -r line
do
    [ "$line" == "" ] && continue
    echo $n_file_current
    if [ $n_file_current -lt $n_file_start ]
    then
	((n_file_current++))
    else
	echo "Processing:" $line
	python script_skimming.py $line output
	full_name=(${line//// }) 
	dataset=${full_name[3]}
	filename=${full_name[7]%.*}
	mv output/${filename}_Skim.root ./
	hadd -f6 ${filename}_Skimmed.root ${filename}_Skim.root
	xrdcp -f -p ${filename}_Skimmed.root root://cmseos.fnal.gov//store/user/xuyan/WtaggingSF/${dataset}_v5/${filename}_Skimmed.root
	((n_file_current++))
	[ $n_file_current -ge $n_file_end ] && break
    fi
    
done < fileLists_Fall17_nanoAODv6/$1



