#!/bin/bash

#SBATCH -J rgnp
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o out/rgnp.%j.out
#SBATCH -e out/rgnp.%j.err
#SBATCH --time=5-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jingma@fredhutch.org
#SBATCH --array=1-20

Rscript "run_sims.R" ${SLURM_ARRAY_TASK_ID}

