#!/bin/bash

#############################################
## 2) extract MCMC samples and postprocess them
#############################################

## Extract results for parametric models 

# for filename in output/posterior_samples/simulation_bimodal/parametric/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

for filename in output/posterior_samples/simulation_bimodal_I_10_N_1000/parametric/*.rds; do
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

## Extract results for parametric models 

# for filename in output/posterior_samples/simulation_unimodal/parametric/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

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

##################################################################


