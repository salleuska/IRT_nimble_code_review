####################################
#### TO FINISH!! ########
####################################
## Run a model multiple times and saves the multiESS
## using the small simulation

library(nimble)
library(mcmcse)

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


## MCMC settings
MCMCcontrol 		<- list()
MCMCcontrol$niter 	<- 10000
MCMCcontrol$nburnin <- 2000


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
runningTime <- system.time({try(
	res <- runMCMC(Cmcmc, 
				   niter 	= MCMCcontrol$niter, 
				   nburnin  = MCMCcontrol$nburnin,
				   setSeed  = seed))
	if(inherits(res, 'try-error')) {
  		warning(paste0("There was a problem running nimble MCMC."))
  }
})

