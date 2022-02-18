##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: Feb 10, 2022
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
## nimble version ??
##-----------------------------------------#
# args <- R.utils::commandArgs(asValue=TRUE)

args <- list()

args$model <- "models/parametric/parametric_IRT_constrainedAbilities.R"
args$dirResults <- "output/posterior_waic/"
args$data <- "data/simulation_unimodal.rds"
args$niter <- 100
args$nthin <- 10
args$mode <- "default"

## Script options from bash
## --model=
## --dirResults=
## --data=
## --niter=
## --nburnin=
## --nthin=
## --mode=

# install_github("nimble-dev/nimble", 
# 	subdir = "packages/nimble", #subdir
# 	ref = "c6137a569",   #branch
# 	lib.loc = "/R/x86_64-pc-linux-gnu-library/dev") #local library
##-----------------------------------------##
## Set variables 
##-----------------------------------------##
calcWAIC <- TRUE

## results directory
if(is.null(args$dirResults)) dir <- "output/posterior_samples" else dir <- args$dirResults

## filename used for output
filename <- unlist(strsplit(basename(args$model), "[\\.]"))[1]

## MCMC settings
MCMCcontrol 		<- list()
MCMCcontrol$niter 	<- as.numeric(args$niter)
MCMCcontrol$nburnin <- as.numeric(args$nburnin)
## thinning for second set of monitors
if(is.null(args$nthin)) MCMCcontrol$thin2 <- 1 else MCMCcontrol$thin2 <- as.numeric(args$nthin)

## set seed based on slurm task id
task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID")
if(task_id == "") seed <- 1 else seed <- 1 + as.numeric(task_id)

MCMCcontrol$seed <- seed

cat("##--------------------------------##\n")
cat("Model ", filename, "\n")
cat("##--------------------------------##\n")

## load library and functions
library(nimble)
source("R_functions/customSamplers.R")
# source("R_functions/monitorLogProb.R")
logLik_summer <- nimbleFunction(
	name = 'logLik_summer',
	contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
          dataNodes <- model$getNodeNames(dataOnly = TRUE)
    },
    run = function() {
        model[[target]] <<- model$getLogProb(dataNodes)  

        copy(from = model, to = mvSaved, 
    	row = 1, nodes = target, logProb = FALSE)

   }, 
   methods = list( reset = function () {})
)
 
############

## Handle data differently if TIMSS (long format)
if(grepl("timss", args$data)){
	alldata <- readRDS(args$data)
	data 	<- list(y = alldata$y)
} else {
	data 	<- list(y = readRDS(args$data))
}
## read model
source(args$model)

##---------------------------------------##
## Initialization
##---------------------------------------##
## init random effects using standardized raw score

if(grepl("timss", args$data)) {
	scores 		<- as.vector(by(alldata, alldata$id, function(x) sum(x$y)/length(x$y), simplify = T))
	Sscores 	<- (scores - mean(scores))/sd(scores)
	inits$eta 	<- Sscores
} else {
	scores 		<- apply(data$y, 1, sum)
	Sscores 	<- (scores - mean(scores))/sd(scores)
	inits$eta 	<- Sscores
}

## BNP inits for data application
if(grepl("bnp", args$model)) {
	if(grepl("health", args$data)) {
		inits$zi 	<- kmeans(Sscores, 3)$cluster
		inits$a 	<- 1        # gamma prior parameters
		inits$b 	<- 3        # gamma prior parameters
	}

	## hyperamenters for simulated data
	if(grepl("simulation", args$data)) {
		inits$zi 	<- kmeans(Sscores, 4)$cluster
		inits$a 	<- 2
		inits$b 	<- 4
	}
	## hyperamenters for timss data
	if(grepl("timss", args$data)) {
		inits$zi 	<- kmeans(Sscores, 3)$cluster
		constants$M <- 30   	# number of clusters
		inits$a 	<- 1        # gamma prior parameters
		inits$b 	<- 3        # gamma prior parameters
	}
}

##---------------------------------------------------##
## Create model and MCMC configuration
##---------------------------------------------------##


model <- nimbleModel(code 			= code2PL,
										 data 			= data,  
										 constants	= constants,
										 inits 			= inits, 
										 calculate 	= FALSE)

# monitors <- c(monitors, "logProbAll", "logProbSum", "logLik")
monitors <- c(monitors, "logLik")
## Flag for WAIC
if(calcWAIC) {
## conditional WAIC - grouped students
  if(grepl("timss", args$data)) {
		indList <- split(seq_along(alldata$id), alldata$id)
		groups <- sapply(indList, function(x)  paste0('y[', x, ']'))  
  }
  else {
		groups <- paste0('y[', 1:constants$N, ', 1:', constants$I, ']')
  }


	mcmcConf <- configureMCMC(model, monitors = monitors, 
		enableWAIC = TRUE, 
		waicControl = list(dataGroups = groups))

} else {
	mcmcConf <- configureMCMC(model, monitors = monitors)

}

## Add samplers to monitor log posterior probability and likelihood
nodeList = c("beta", "lambda", "eta")
if("gamma" %in% monitors) nodeList[1] <- "gamma"

# mcmcConf$removeSampler("logProbAll", "logProbSum", "logLik")
# mcmcConf$addSampler("logProbSum", type =  "logProb_summer", 
#  				control = list(nodeList = nodeList) )

# mcmcConf$removeSampler("logProbAll")
# mcmcConf$addSampler("logProbAll", type =  "logProb_summer")

mcmcConf$removeSampler("logLik")
mcmcConf$addSampler("logLik", type =  "logLik_summer")

mcmcConf

## sampler configuration changes according to mode
if(args$mode == "centered" ) {

	if(("gamma" %in% monitors) & grepl("constrainedAbilities|unconstrained", filename)){ 

		  mcmcConf$removeSamplers("log_lambda")
		  # mcmcConf$removeSamplers("gamma")
  
 		  mcmcConf$addSampler(type = 'centered',
       			target = c('log_lambda', 'gamma'),
       			control = list(nodesToCenter = 'eta', scale = 0.1, adaptive = TRUE))

	} else {
	  q(save = 'no')
	}
}


mcmcConf$addMonitors2("eta")
mcmcConf$setThin2(MCMCcontrol$thin2)
mcmc <- buildMCMC(mcmcConf)	


##---------------------------------------------------##
## Compile model & MCMC 
##---------------------------------------------------##
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


##---------------------------------------------------##
## Calculate WAIC
##---------------------------------------------------##
# if(calcWAIC) {
# 	modelWAIC <- try(Cmcmc$calculateWAIC(nburnin = MCMCcontrol$nburnin))
# 	if(inherits(res, 'try-error')) {
# 	  	warning(paste0("There was a problem running nimble MCMC."))
# 	  }

# 	  ##  "manual" thinning of etaSamples 
# 	  ## so other scripts do not need to be modified
# 	  res <- list(samples = res)

# 	  indicesEta <- seq(from = 1,  
# 				  to = MCMCcontrol$niter - MCMCcontrol$nburnin, 
# 				  by = MCMCcontrol$thin2)

# 	  etaSampIndex <- grep("^eta", colnames(res$samples))
# 	  logLambdaSampIndex <- grep("^log_lambda", colnames(res$samples))

# 	  res$samples2 <- res$samples[indicesEta, etaSampIndex]
# 	  res$samples <- res$samples[, -c(logLambdaSampIndex, etaSampIndex)]
# }
##---------------------------------------------------##
## Save results, times, settings
##---------------------------------------------------##
results <- list(samples  = res,
				compilationTime  = compilationTime,
				samplingTime     = runningTime*(1 - MCMCcontrol$niter/MCMCcontrol$nburnin),
				runningTime      = runningTime,
				MCMCcontrol      = MCMCcontrol)

if(calcWAIC) results$modelWAIC <- Cmcmc$getWAIC()$WAIC


##---------------------------------------------------##
## directory for output
##---------------------------------------------------##
# modelType       <- unlist(strsplit(basename(args$model), "[\\_\\.]"))[1]
# dataName        <- unlist(strsplit(basename(args$data), "[.]"))[1]

# outDir <- paste0(dir, "/", dataName, "/", modelType, "/")

# dir.create(file.path(outDir), recursive = TRUE, showWarnings = FALSE)


# if(grepl("centered", args$mode)) {
# 	filenameOutput <- paste0(outDir, filename, "_centered.rds")
# } else {
# 	filenameOutput <- paste0(outDir, filename, ".rds")
# }

# saveRDS(results, file = filenameOutput )
