#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o stan_bimodal.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#SBATCH --array=1-10
###########################

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal.rds \
--nsamples=10000 \
--nwarmup=5000 \
--mode=rep