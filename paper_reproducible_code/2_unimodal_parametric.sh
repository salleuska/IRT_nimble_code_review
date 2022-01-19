#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=sally.paganin@berkeley.edu
#SBATCH -o unimodal_parametric_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5
######################


FILES=(models/parametric/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=default


Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=centered


# Rscript 1_runNimbleModels.R  \
# --model=models/parametric/parametric_IRT_constrainedAbilities.R \
# --data=data/simulation_unimodal.rds \
# --niter=500 \
# --nburnin=50 \
# --nthin=1 \
# --mode=default

