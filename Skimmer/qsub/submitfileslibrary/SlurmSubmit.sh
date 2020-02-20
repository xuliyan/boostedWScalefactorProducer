#!/bin/bash
# 
#SBATCH -p wn
#SBATCH --account=t3
#SBATCH --time 08:59:00
#SBATCH --mem=6gb       # Job memory request
#SBATCH -e cn-test.err 
#SBATCH -o cn-test.out  # replace default slurm-SLURM_JOB_ID.out
#SBATCH --get-user-env


echo HOME: $HOME 
echo USER: $USER 
echo SLURM_JOB_ID: $SLURM_JOB_ID
echo HOSTNAME: $HOSTNAME

# each worker node has local /scratch space to be used during job run
mkdir -p /scratch/$USER/${SLURM_JOB_ID}
export TMPDIR=/scratch/$USER/${SLURM_JOB_ID}

sleep 10
# here comes a computation

JobList=$1
TaskCmd=$(cat $JobList | sed ''$TaskID'q;d')

echo "Going to execute" $TaskCmd
eval $TaskCmd

Joblist=$1
TaskCmd=$(cat $JobList | sed ''${SLURM_ARRAY_TASK_ID}'q;d')

eval $TaskCmd

# cleaning of temporal working dir when job was completed:
rmdir  -rf /scratch/$USER/${SLURM_JOB_ID}

date
