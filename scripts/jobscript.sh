#!/bin/bash

#SBATCH --job-name=nanoclass
#SBATCH --output=cluster/%x-%u-%A-%a.log
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=...@...
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=58
#SBATCH --time=4000
#SBATCH --mem=460G


start=`date "+%s"`
echo "$SLURM_JOB_NAME started at `date` on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."


## Conda init
__conda_setup="$('/home/$USER/personal/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/$USER/personal/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/$USER/personal/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/$USER/personal/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup


## Make sure to use programs on the amplicomics group share
conda activate /zfs/omics/projects/amplicomics/miniconda3/envs/snakemake-ampl

srun mkdir -p /scratch/$USER/tmp/
export TMPDIR=/scratch/$USER/tmp/

## Run nanoclass
cmd="srun --cores $SLURM_CPUS_ON_NODE snakemake --use-conda --cores $SLURM_CPUS_ON_NODE --nolock --rerun-incomplete"
echo "Running: $cmd"
eval $cmd


## Create report
#cmd="srun snakemake --report report/NanoClass-`date "+%Y%m%dT%H%M"`.zip"
#echo "Running: $cmd"
#eval $cmd


end=`date "+%s"`
runtime=$((end-start))
echo "$SLURM_JOB_NAME finished at `date` in $runtime seconds."
