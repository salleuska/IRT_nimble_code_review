##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: June, 22 2021
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## nimble version 0.11.1
##-----------------------------------------##
# This files compute some quantities for plots and tables
dir.create("figures/dataForFigures", recursive = TRUE, showWarnings = FALSE)
##-----------------------------------------#
## Data
# rm(list = ls()# 
library(bayestestR) ## For HDI intervals
# ## Compute simulation MSE
source("R_functions/ggplot_settings.R")
source("R_functions/multimodalDensity.R")


##--------------------------------##
## Unimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_unimodal"
load(paste0("data/",dataName,"_allValues.RData"))
## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				     eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

metricsUnimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))

##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)
 
## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
indices <- seq(10, 45000, by = 10)
bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3])))*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
# indexSample <- seq(1, 2000, length = 50)
# truePerc <- pnorm(etaAbility[indexSample], 0, sd = 1.25)

truePerc <- pnorm(etaAbility, 0, sd = 1.25)
indexSample <- order(truePerc)[round(seq(1, 2000, length = 50))]

## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}



##########################

unimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure, 
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(unimodalRes, file = "figures/dataForFigures/unimodal.rds")

##--------------------------------##
## Bimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_bimodal"
load(paste0("data/",dataName,"_allValues.RData"))

## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				   eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

## Save simulation metrics computed using true values
metricsBimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_unimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_unimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))

save(metricsUnimodal, metricsBimodal, file = "figures/dataForFigures/table3_metricsSimulations.RData")
#--------------------------------------------------------------------------------#
# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2,function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#

niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid + paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
indices <- seq(10, 45000, by = 10)
bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute percentiles simulation - qqplot
##------------------------------------------------------------#
## Take a grid
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
# indexSample <- seq(1, 2000, length = 50)
# truePerc <- 0.5*pnorm(etaAbility[indexSample] , -2, sd = 1.25) +
# 					 0.5*pnorm(etaAbility[indexSample] , 2, sd = 1.25)

truePerc <- 0.5*pnorm(etaAbility , -2, sd = 1.25) +
					 0.5*pnorm(etaAbility , 2, sd = 1.25)

# f.freq <- function(x)
# {
#   tab <- data.frame(table(x))
#   tab$Percent <- tab$Freq*100/length(x)
#   tab$Cum.Percent[1] <- tab$Percent[1]
#   for(i in 2:length(tab[,1]))
#     tab$Cum.Percent[i] <- tab$Cum.Percent[i-1] + tab$Percent[i]
#   tab
# }

# f.freq(scores)

# scores <- apply(Y, 1, sum)
# plot(ecdf(scores))

# plot(etaAbility, scores)
# plot(apply(paraModel$etaSamp, 2, mean), scores)
# plot(apply(bnpModel$etaSamp, 2, mean), scores)


## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}


quantSampBnp <- apply(bnpPerc, 2, function(x) quantile(x, seq(0.1, 0.9, by = 0.1)))
quantBNP <- apply(quantSampBnp, 1, mean)
quantTrue <- quantile(truePerc, seq(0.1, 0.9, by = 0.1))

plot(quantTrue, quantBNP)

plot(truePerc)
points(1:2000, apply(paraPerc, 2, mean), col = 4)
points(1:2000, apply(bnpPerc, 2, mean), col = 3)



plot(sort(truePerc), apply(paraPerc, 2, mean)[order(truePerc)], ylim = c(0,1))

plot(truePerc, apply(bnpPerc, 2, mean), ylim = c(0,1))
abline(0,1)

plot(apply(bnpPerc, 2, mean),apply(paraPerc, 2, mean), ylim = c(0,1))

##-----------------------------------------------##

bimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(bimodalRes, file = "figures/dataForFigures/bimodal.rds")
##------------------------------------------------------------#
##--------------------------------##
## multimodal simulation
##--------------------------------##
## load true values
dataName <- "simulation_multimodal"
#dataName <- "simulation_multimodal"
load(paste0("data/",dataName,"_allValues.RData"))
## Set a grid for density computation
grid <- seq(-8, 8, len = 400) 

trueValues <- list(beta   = beta0,
		   		   lambda = lambda0,
				     eta    = etaAbility)

bestModel <- modelData$model[modelData$data == dataName]


bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


biasParametric <-  sapply(1:3, function(i) paraEstimates[[i]] - trueValues[[i]])
biasBnp        <-  sapply(1:3, function(i) bnpEstimates[[i]] - trueValues[[i]])

metricsMultimodal <- data.frame(parameters =c("Difficulties", "Discrimination", "Abilities"),  
							   MAE_multimodal_para = unlist(lapply(biasParametric, function(x) mean(abs(x)))), 
							   MSE_multimodal_para = unlist(lapply(biasParametric, function(x) mean(x^2))),
							   MAE_multimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(abs(x)))), 
							   MSE_multimodal_bnp  = unlist(lapply(biasBnp, function(x) mean(x^2))))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
# indices <- seq(10, 45000, by = 10)
#bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3])))*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute percentile using a fixed grid
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
# truePerc <- pnorm(etaAbility , 0, sd = 1.25)

grid <- seq(-5, 4, len = 20) 

truePerc <- pMultiModal(grid,
              weights = c(2,4,4), 
              means= c(-2, 0, 3))

paraPerc  <- matrix(0, ncol = length(grid), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- grid*paraModel$scaleShiftEta[i] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}
bnpPerc  <- matrix(0, ncol = length(grid), nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*grid + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

# plot(apply(bnpPerc, 2, mean), truePerc)
# abline(0,1)
# plot(apply(paraPerc, 2, mean), truePerc)

## Plot ##
library(ggplot2)
library(bayestestR)
dfPercentile <- data.frame(ind       = rep(1:length(grid), 2),
                           estimate  = c(apply(paraPerc, 2, mean),
                                         apply(bnpPerc, 2, mean)),
                           CI_low    = c(apply(paraPerc, 2, quantile, 0.025),
                                         apply(bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(paraPerc, 2, quantile, 0.975),
                                         apply(bnpPerc, 2, quantile, 0.975)))


dfPercentile$trueVal <- rep(truePerc, 2)
dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pUniPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title=element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
        	width = 0.8, position = position_dodge(width = 0.6)) + 
        geom_point(aes(x = ind, y = trueVal*100, fill = "True value"), color = "black", size = 1.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank"), 
                                                     shape    = c(16, 16),
                                                     color    = c(paraColor, bnpColor)))) + 
        labs(y = "Percentile", x = "Individual", 
                title = "Multimodal simulation") 	

ggsave(filename = "figures/REV_fixed_grid_percentiles.pdf", plot = pUniPerc,
        width = 20, height = 9 , dpi = 300, units = "cm", device='pdf')

##------------------------------------------------------------#
## Compute individual percentile - us
##------------------------------------------------------------#
load(paste0("data/", dataName,"_allValues.RData"))
# truePerc <- pnorm(etaAbility , 0, sd = 1.25)

indexSample <- round(seq(1, 2000, length = 50))
truePerc <- pMultiModal(etaAbility[indexSample],
              weights = c(2,4,4), 
              means= c(-2, 0, 3))

## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

# plot(apply(etaSamplesPara, 2, mean), etaAbility[indexSample])
## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}


## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}



dfPercentile <- data.frame(ind       = rep(1:length(truePerc), 2),
                           estimate  = c(apply(paraPerc, 2, mean)[order(truePerc)],
                                         apply(bnpPerc, 2, mean)[order(truePerc)]),
                           CI_low    = c(apply(paraPerc, 2, quantile, 0.025)[order(truePerc)],
                                         apply(bnpPerc, 2, quantile, 0.025)[order(truePerc)]),
                           CI_upp    = c(apply(paraPerc, 2, quantile, 0.975)[order(truePerc)],
                                         apply(bnpPerc, 2, quantile, 0.975)[order(truePerc)]))



dfPercentile$trueVal <- rep(truePerc[order(truePerc)], 2)
dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pUniPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title=element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
        	width = 0.8, position = position_dodge(width = 0.6)) + 
        geom_point(aes(x = ind, y = trueVal*100, fill = "True value"), color = "black", size = 1.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank"), 
                                                     shape    = c(16, 16),
                                                     color    = c(paraColor, bnpColor)))) + 
        labs(y = "Percentile", x = "Individual", 
                title = "Multimodal simulation") 	
pUniPerc

ggsave(filename = "figures/REV_multimodal_percentiles.pdf", plot = pUniPerc,
        width = 20, height = 9 , dpi = 300, units = "cm", device='pdf')

pdf("REV_ability_estimate.pdf")
plot(etaAbility[indexSample], apply(etaSamplesPara, 2, mean), col = paraColor, pch = 16, 
	xlab = "true ability", ylab = "estimated ability")
points(etaAbility[indexSample], apply(etaSamplesBnp, 2, mean), col = bnpColor, pch = 16)
abline(0,1)
dev.off()

### test n 3 
## take the original answer and scores

scores <- apply(Y, 1, sum)
x.ranks <- ecdf(scores)
val <- sort(unique(scores))
empirical <- cbind(val,cum.prop=x.ranks(val))

trueAbility <- 


##########################

multimodalRes <- list(truValues = trueValues,
					paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure, 
					truePerc = truePerc[indexSample],
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)

saveRDS(multimodalRes, file = "figures/dataForFigures/multimodal.rds")
##------------------------------------------------------------#
## Data
##------------------------------------------------------------#
##--------------------------------##
## Health data
##--------------------------------##
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-15, 30, len = 800) 

dataName <- "data_health"

bestModel <- modelData$model[modelData$data == dataName]
# bestModel <- "SI_unconstrained"

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#

niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  # rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  rescGrid <- paraModel$scaleShiftEta[i]*grid +  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
niter2 <- length(bnpG0)
indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
	# rescGrid <- bnpModel$scaleShiftEta[i]*grid -  bnpModel$locationShiftEta[i]
	rescGrid <- bnpModel$scaleShiftEta[i]*grid +  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}

##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
indexSample <- seq(1, length(bnpEstimates$eta), length = 50)

## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]

## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

#####

healthRes   <- list(paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)


saveRDS(healthRes, file = "figures/dataForFigures/health.rds")

##--------------------------------##
## Timss data
##--------------------------------##
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-8, 8, len = 800) 

dataName <- "data_timss"
bestModel <- modelData$model[modelData$data == dataName]

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp/bnp_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric/parametric_", bestModel, ".rds"))


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
					  lambda = apply(paraModel$lambdaSamp, 2, mean), 
					  eta    = apply(paraModel$etaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
					  lambda = apply(bnpModel$lambdaSamp, 2, mean), 
					  eta    = apply(bnpModel$etaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))



bnpLow <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
				 lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
			     eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
				   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
			       eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid -  paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp_", bestModel, ".rds"))
niter2 <- length(bnpG0)
indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]

densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
  	rescGrid <- bnpModel$scaleShiftEta[i]*grid -  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#

indexSample <- seq(1, length(bnpEstimates$eta), length = 50)

## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] - paraModel$locationShiftEta[i]
	 paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
	 rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] - bnpModel$locationShiftEta[i]
	 bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

##------------------------------------------------------------#

timssRes   <- list(paraEstimates = paraEstimates, 
					paraLow = paraLow, 
					paraUpper = paraUpper,
					bnpEstimates = bnpEstimates, 
					bnpLow = bnpLow, 
					bnpUpper = bnpUpper,
					grid = grid, 
					densitySamplesPara = densitySamplesPara, 
					densityDPMeasure = densityDPMeasure,
					paraPerc = paraPerc, 
					bnpPerc = bnpPerc)


saveRDS(timssRes, file = "figures/dataForFigures/timss.rds")
##------------------------------------------------------------#
## TIMSS 3PLx
##------------------------------------------------------------#
rm(list = ls());gc()
library(bayestestR) ## For HDI intervals

## Compute simulation MSE
source("R_functions/ggplot_settings.R")

## Set a grid for density computation
grid <- seq(-10, 8, len = 800) 

dataName <- "data_timss"
bestModel <- "IRT_unconstrained"

bnpModel  <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/bnp3PL/bnp3PL_", bestModel, ".rds"))
paraModel <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/parametric3PL/parametric3PL_", bestModel, ".rds"))

paraModel$deltaSamp <- paraModel$otherParSamp[, grep("delta", colnames(paraModel$otherParSamp))]
bnpModel$deltaSamp  <- bnpModel$otherParSamp[, grep("delta", colnames(paraModel$otherParSamp))]


paraEstimates <- list(beta   = apply(paraModel$betaSamp, 2, mean),
                  lambda = apply(paraModel$lambdaSamp, 2, mean), 
                  eta    = apply(paraModel$etaSamp, 2, mean), 
                  delta  = apply(paraModel$deltaSamp, 2, mean))

bnpEstimates  <- list(beta   = apply(bnpModel$betaSamp, 2, mean),
                      lambda = apply(bnpModel$lambdaSamp, 2, mean), 
                      eta    = apply(bnpModel$etaSamp, 2, mean), 
                      delta  = apply(bnpModel$deltaSamp, 2, mean))


# Compute HDI for item and abilities estimates
paraLow <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
                lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
                eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low), 
                delta    = apply(paraModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))

paraUpper <-  list(beta   = apply(paraModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
                   lambda = apply(paraModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                    eta    = apply(paraModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high),
                    delta   = apply(paraModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))


bnpLow <-   list(beta    = apply(bnpModel$betaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low),
                lambda   = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95,  method = "HDI")$CI_low), 
                eta      = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low), 
                delta    = apply(bnpModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_low))



bnpUpper <-   list(beta   = apply(bnpModel$betaSamp, 2, function(x) ci(x,ci = 0.95, method = "HDI")$CI_high),
                   lambda = apply(bnpModel$lambdaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                   eta    = apply(bnpModel$etaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high), 
                   delta  = apply(bnpModel$deltaSamp, 2, function(x) ci(x, ci = 0.95, method = "HDI")$CI_high))

##------------------------------------------------------------#
## Compute abilities distribution (predictive) 
## - using samples from the DP
## - parametric counterpart
##------------------------------------------------------------#
niter <- nrow(paraModel$etaSamp)

## Samples for parametric density
densitySamplesPara <- matrix(0, ncol = length(grid), nrow = niter)
muParaSamples <- paraModel$otherParSamp[ ,"mu"]
s2ParaSamples <- paraModel$otherParSamp[ ,"s2.eta"]

for(i in 1:niter){
  rescGrid <- paraModel$scaleShiftEta[i]*grid + paraModel$locationShiftEta[i]
  
  densitySamplesPara[i, ] <- sapply(rescGrid,
                function(x) dnorm(x, muParaSamples[i], sqrt(s2ParaSamples[i]))*paraModel$scaleShiftEta[i])
}

## Samples BNP density

bnpG0 <- readRDS(paste0("output/posterior_samples_elaborated/", dataName, "/DPG0_bnp3PL_", bestModel, ".rds"))
niter2 <- length(bnpG0)

indices <- seq(10, niter2, by = 10)

bnpG0 <- bnpG0[indices]
densityDPMeasure <- matrix(0, ncol = length(grid), nrow = length(bnpG0))

for(i in seq_len(length(bnpG0))) {  
        rescGrid <- bnpModel$scaleShiftEta[i]*grid +  bnpModel$locationShiftEta[i]

    densityDPMeasure[i, ] <- sapply(rescGrid,
                function(x)(
                  sum(bnpG0[[i]][,1]*
                    dnorm(x, bnpG0[[i]][,2], sqrt(bnpG0[[i]][,3]))
                    )*bnpModel$scaleShiftEta[i]))
}
##------------------------------------------------------------#
## Compute individual percentile
##------------------------------------------------------------#
indexSample <- seq(1, length(bnpEstimates$eta), length = 50)
## take some individuals
etaSamplesPara   <- paraModel$etaSamp[, indexSample]


## Take a grid
paraPerc  <- matrix(0, ncol = dim(etaSamplesPara)[2], nrow = niter)

for(i in 1:niter) {  
         rescaled <- paraModel$scaleShiftEta[i]*etaSamplesPara[i,] + paraModel$locationShiftEta[i]
         paraPerc[i, ] <- sapply(rescaled, function(x)
                                   pnorm(x,
                                   mean = muParaSamples[i],
                                   sd   = sqrt(s2ParaSamples[i]),
                                   lower.tail = TRUE ))
}

## take some individuals
etaSamplesBnp   <- bnpModel$etaSamp[, indexSample]

bnpPerc  <- matrix(0, ncol = dim(etaSamplesBnp)[2], nrow = niter)

for(i in 1:niter) {  
         rescaled <- bnpModel$scaleShiftEta[i]*etaSamplesBnp[i,] + bnpModel$locationShiftEta[i]
         bnpPerc[i, ] <- sapply(rescaled, function(x) sum(bnpG0[[i]][,1] *pnorm(x,
                                   mean = bnpG0[[i]][,2],
                                   sd   = sqrt(bnpG0[[i]][,3]),
                                   lower.tail = TRUE )))
}

##------------------------------------------------------------#
timssRes   <- list(paraEstimates = paraEstimates, 
                                        paraLow = paraLow, 
                                        paraUpper = paraUpper,
                                        bnpEstimates = bnpEstimates, 
                                        bnpLow = bnpLow, 
                                        bnpUpper = bnpUpper,
                                        grid = grid, 
                                        densitySamplesPara = densitySamplesPara, 
                                        densityDPMeasure = densityDPMeasure,
                                        paraPerc = paraPerc, 
                                        bnpPerc = bnpPerc)

saveRDS(timssRes, file = "figures/dataForFigures/timss3PL.rds")
