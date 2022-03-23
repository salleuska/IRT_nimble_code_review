#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=sally.paganin@berkeley.edu
#SBATCH -o timss_3PLpara_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5
######################

FILES=(models/parametric3PL_long/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=output/posterior_samples \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=output/posterior_samples \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered