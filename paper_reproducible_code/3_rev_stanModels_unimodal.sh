#!/bin/bash
#SBATCH --mail-type=ALL                       
#SBATCH --mail-user=spaganin@hsph.harvard.edu
#SBATCH -o stan_unimodal_%j.out                 # File to which STDERR will be written, including job ID
#SBATCH --cpus-per-task=1
#########################

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal.rds \
--nsamples=4000 \
--nwarmup=4000 

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal_I_30_N_1000.rds \
--nsamples=4000 \
--nwarmup=4000 

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal_I_30_N_5000.rds \
--nsamples=4000 \
--nwarmup=4000 

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal_I_10_N_5000.rds \
--nsamples=4000 \
--nwarmup=4000 

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal_I_10_N_1000.rds \
--nsamples=4000 \
--nwarmup=4000 
