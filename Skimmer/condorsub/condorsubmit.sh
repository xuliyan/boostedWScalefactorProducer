#!/bin/bash



#for i in {"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.txt","ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.txt","ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8.txt","ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt","ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt"}

#for i in {"WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.txt","WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8.txt"}
for i in {"TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt","234"}
#for i in {"WZ_TuneCP5_13TeV-pythia8.txt","ZZ_TuneCP5_13TeV-pythia8.txt","WW_TuneCP5_13TeV-pythia8.txt"}
#for i in {"SingleMuon-Run2017B.txt","SingleMuon-Run2017C.txt","SingleMuon-Run2017D.txt","SingleMuon-Run2017E.txt","SingleMuon-Run2017F.txt"}
#for i in {"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.txt","TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.txt"}

#for i in {"WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.txt","123.txt"}

do
    dataset=${i%%.txt}
    mkdir $dataset
    mkdir $dataset/Condor_error
    mkdir $dataset/Condor_output
    mkdir $dataset/Condor_log
    echo `cat fileLists_Fall17_nanoAODv6/${i} | wc -l`
    n_file=`cat fileLists_Fall17_nanoAODv6/${i} | wc -l`
    #n_job=$((n_file/1 + 1))
	n_job=$((n_file/1))
    cat > Condor_submit_${i%%.txt}.jdl <<EOF
universe                = vanilla
executable              = Skimmer_condor.sh
arguments               = ${i} \$(Process)
should_transfer_files   = YES
transfer_input_files    = transfer.tgz
transfer_output_files   = _condor_stderr, _condor_stdout
when_to_transfer_output = ON_EXIT
output                  = ${dataset}/Condor_output/Condor_job.\$(Cluster).\$(Process).out
error                   = ${dataset}/Condor_error/Condor_job.\$(Cluster).\$(Process).err
log                     = ${dataset}/Condor_log/Condor_job.\$(Cluster).log
request_cpus            = 1
request_memory          = 1 GB
+JobFlavour             = "testmatch"

queue ${n_job}

EOF

    condor_submit Condor_submit_${i%%.txt}.jdl
done
