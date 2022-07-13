#!/bin/bash

#############################################
## 2) extract MCMC samples and postprocess them
#############################################

## Extract results for parametric models 


for filename in output/posterior_samples/simulation_bimodal_I_10_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##########
## TERM 1

for filename in output/posterior_samples/simulation_bimodal/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_10_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_30_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_30_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##################################################################

##########
## TERM 2

## Extract results for parametric models 

for filename in output/posterior_samples/simulation_unimodal/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_10_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_10_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_30_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_30_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##################################################################

##########
## TERM 3

## Extract results for parametric models 

for filename in output/posterior_samples/simulation_multimodal/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_multimodal_I_10_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done


for filename in output/posterior_samples/simulation_multimodal_I_10_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_multimodal_I_30_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_multimodal_I_30_N_5000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done
### done sim 2

##################################################################



for filename in output/posterior_samples/simulation_bimodal_I_10_N_1000/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##########
## TERM 1

for filename in output/posterior_samples/simulation_bimodal/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_10_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_10_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_30_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal_I_30_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done


##########
## TERM 2

## Extract results for bnp models 

for filename in output/posterior_samples/simulation_unimodal/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_10_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_10_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_30_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_unimodal_I_30_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##################################################################

##########
## TERM 3

## Extract results for bnp models 

## sim2 - local
for filename in output/posterior_samples/simulation_multimodal/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

## sim2 - server
for filename in output/posterior_samples/simulation_multimodal2_I_10_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

## sim2 - server
for filename in output/posterior_samples/simulation_multimodal2_I_10_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

## sim2 - rerun - wait
for filename in output/posterior_samples/simulation_multimodal2_I_30_N_1000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_multimodal2_I_30_N_5000/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

##################################################################
## Stan models

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan2.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan3.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan4.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan5.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan6.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan7.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan8.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan9.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan10.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan11.rds


Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan2.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan3.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan4.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan5.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan6.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan7.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan8.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan9.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan10.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan11.rds



Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan.rds

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal_I_10_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal_I_10_N_5000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal_I_30_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal_I_30_N_5000/parametric/parametric_IRT_stan.rds

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal_I_10_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal_I_10_N_5000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal_I_30_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal_I_30_N_5000/parametric/parametric_IRT_stan.rds

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal_I_10_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal_I_10_N_5000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal_I_30_N_1000/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal_I_30_N_5000/parametric/parametric_IRT_stan.rds


Rscript 2_extractResults.R --resFileName=output/posterior_samples/data_timss/parametric/parametric_IRT_stan.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/data_health/parametric/parametric_IRT_stan.rds

##################
## TMP for plot
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal/parametric/parametric_IRT_constrainedItem.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal/parametric/parametric_IRT_unconstrained.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_unimodal/parametric/parametric_IRT_stan.rds

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_constrainedItem.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_unconstrained.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_bimodal/parametric/parametric_IRT_stan.rds

Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_constrainedItem.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_unconstrained.rds
Rscript 2_extractResults.R --resFileName=output/posterior_samples/simulation_multimodal/parametric/parametric_IRT_stan.rds
