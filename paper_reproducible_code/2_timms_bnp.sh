#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o timss_bnp_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=0-3
######################

FILES=(models/bnp_long/*.R)

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=/scratch/users/sallypaganin \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=25000 \
--nthin=1 \
--mode=default

Rscript 1_runNimbleModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--dirResults=/scratch/users/sallypaganin \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=25000 \
--nthin=1 \
--mode=centered