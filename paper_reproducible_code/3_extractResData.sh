#!/bin/bash
##################################################################

## Extract results for parametric models 
for filename in /scratch/users/sallypaganin/data_timss/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename 
done

# for filename in output/posterior_samples/data_health/parametric/*.rds; do
# 	echo $filename
# 	Rscript 2_extractResults.R --resFileName=$filename
# done

############################################################
## Extract results for semiparametric models 

for filename in output/posterior_samples/data_timss/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done


for filename in output/posterior_samples/data_health/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

################