#!/bin/bash

#SBATCH --job-name=NanoClass2
#SBATCH --output=%x-%u-%A-%a.log
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=...@...
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=58
#SBATCH --time=4000
#SBATCH --mem=460G

# Set LC_ALL and export it
export LC_ALL=en_US.UTF-8

start=`date "+%s"`
echo "$SLURM_JOB_NAME started at `date` on node $SLURM_NODEID using $SLURM_CPUS_ON_NODE cpus."

#ensure that conda is recognized
source ~/.bashrc

## Make sure to use the snakemake env installed in the amplicomics share
conda activate /zfs/omics/projects/amplicomics/miniconda3/envs/snakemake_nanoclass2

srun mkdir -p /scratch/$USER/tmp/
export TMPDIR=/scratch/$USER/tmp/

## Run nanoclass2
cmd="srun --cores $SLURM_CPUS_ON_NODE snakemake \
    -s /zfs/omics/projects/amplicomics/bin/NanoClass2/Snakefile \
    --configfile config.yaml --use-conda \
    --conda-prefix /zfs/omics/projects/amplicomics/bin/NanoClass2/.snakemake/conda \
    --cores $SLURM_CPUS_ON_NODE --nolock --rerun-incomplete"

echo "Running: $cmd"
eval $cmd


## Create report
#cmd="srun snakemake --report report/NanoClass-`date "+%Y%m%dT%H%M"`.zip"
#echo "Running: $cmd"
#eval $cmd


end=`date "+%s"`
runtime=$((end-start))
echo "$SLURM_JOB_NAME finished at `date` in $runtime seconds."

