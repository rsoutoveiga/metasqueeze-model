#!/bin/bash
#SBATCH --ntasks=2
#SBATCH --mem=1024
#SBATCH --job-name=job
#SBATCH --error=./tools/parallel/sbatch_out/job_meta.%J.stderr
#SBATCH --output=./tools/parallel/sbatch_out/job_meta.%J.stdout
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rsoutoveiga@uni-potsdam.de

cd ${SLURM_SUBMIT_DIR}

Rscript tools/parallel/metasqueeze-parallel.R $1 $2 $3
