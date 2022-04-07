####################################
####################################
## Run a model multiple times and saves the multiESS
## using the small simulation

library(nimble)
library(mcmcse)
source("R_functions/rescalingFunctions.R")

args <- R.utils::commandArgs(asValue=TRUE)


MCMCcontrol 		<- list()
MCMCcontrol$niter 	<- as.numeric(args$niter)
MCMCcontrol$nburnin <- as.numeric(args$nburnin)

# ## MCMC settings
# MCMCcontrol 		<- list()
# MCMCcontrol$niter 	<- 10000
# MCMCcontrol$nburnin <- 1000


data 	<- list(y = readRDS("data/simulation_unimodal_I_10_N_1000.rds"))

source("models/parametric/parametric_IRT_constrainedAbilities.R")

scores 		<- apply(data$y, 1, sum)
Sscores 	<- (scores - mean(scores))/sd(scores)
inits$eta 	<- Sscores


model <- nimbleModel(code = code,
					 data 		= data,  
					 constants	= constants,
					 inits 		= inits, 
					 calculate 	= FALSE)


## update monitors
mcmcConf <- configureMCMC(model, monitors = monitors)
mcmcConf$addMonitors("eta")

mcmc <- buildMCMC(mcmcConf)	




compilationTime <- system.time({
    Cmodel <- try(compileNimble(model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble model.")
    }
    Cmcmc <- try(compileNimble(mcmc, project = model))
    if(inherits(Cmodel, 'try-error')) {
      stop("There was a problem compiling the nimble MCMC.")
    }
})


##---------------------------------------------------##
## Run MCMC 
##---------------------------------------------------##
nTimes <- 20
out <- numeric(20)

for(i in 1:nTimes){

	runningTime <- system.time({try(
		res <- runMCMC(Cmcmc, 
					   niter 	= MCMCcontrol$niter, 
					   nburnin  = MCMCcontrol$nburnin, 
					   setSeed = i))
		if(inherits(res, 'try-error')) {
	  		warning(paste0("There was a problem running nimble MCMC."))
	  }
	})


	modelRes <- posteriorRescalingBeta(samples  = res[, -grep("^eta", colnames(res))],
	                                   samples2 = res[, grep("^eta", colnames(res))],
	                                   thinEta  = 1, 
					                  rescale  = TRUE)

	## matrices for ESS evaluations
	onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
	              modelRes[grepl("beta", names(modelRes))][[1]])

	itemsAndAbility <- cbind(onlyItems, modelRes$etaSamp)
	itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]

	try(out[i] <- mcmcse::multiESS(itemsAndAbilityMultiESS, 
	                                            method = "bm", r = 1, 
	                                            adjust = FALSE))
}

save(out, file = paste0("output/multiESS_", MCMCcontrol$niter, ".rds"))

######################
## Plots

# d1 <- data.frame(mESS = readRDS("output/multiESS_10000.rds"))
# d2 <- data.frame(mESS = readRDS("output/multiESS_20000.rds"))
# d4 <- data.frame(mESS = readRDS("output/multiESS_40000.rds"))
# d5 <- data.frame(mESS = readRDS("output/multiESS_50000.rds"))

# d1$Niter <- paste0("N samples = " , 9000, "\n(niter = 10000, 5% burnin)")
# d2$Niter <- paste0("N samples = " , 18000, "\n(niter = 20000, 5% burnin)")
# d4$Niter <- paste0("N samples = " , 36000, "\n(niter = 40000, 5% burnin)")
# d5$Niter <- paste0("N samples = " , 45000, "\n(niter = 50000, 5% burnin)")

# library(ggplot2)

# df <- rbind(d1, d2, d4, d5)
# df$Niter <- factor(df$Niter, levels = levels(as.factor(df$Niter))[c(4,1, 2,3 )])

# bp <- ggplot(df, aes(x=Niter, y=mESS, group=Niter)) + 
#   geom_boxplot() + theme_bw() +
#   theme(axis.text.x=element_blank() )+
#   facet_wrap(~Niter, ncol = 2, scales = "free") + 
#   ggtitle("NIMBLE model - mESS across 20 replications")
# bp
# ggsave(bp, file = "output/multiESSNimble.pdf")

# # boxplot(mESS ~ Niter, data = dz)


# # boxplot(mESS ~ Niter, data = df)

# s1 <- data.frame(mESS = readRDS("output/Stan_multiESS__warmup_5000.rds"))
# s2 <- data.frame(mESS = readRDS("output/Stan_multiESS__warmup_10000.rds"))

# d4$Niter <- paste0("N samples = " , 5000, "\n warmup = 5000")
# d5$Niter <- paste0("N samples = " , 10000 "\n warmup = 10000")


# df <- rbind(d1, d2)
# df$Niter <- factor(df$Niter, levels = levels(as.factor(df$Niter))[c(2,1)])

# bp <- ggplot(df, aes(x=Niter, y=mESS, group=Niter)) + 
#   geom_boxplot() + theme_bw() +
#   theme(axis.text.x=element_blank() )+
#   facet_wrap(~Niter, ncol = 2, scales = "free") + 
#   ggtitle("Stan Model - mESS across 20 replications")

# bp  
# ggsave(bp, file = "output/multiESSStan.pdf")




