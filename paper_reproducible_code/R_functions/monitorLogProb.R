#############
## Custom sampler to monitor posterior logProbabilities
############# 
library(nimble)
############# 

## sampler to monitor 
logProb_summer <- nimbleFunction(
	name = 'logProb_summer',
	contains = sampler_BASE,
    setup = function(model, mvSaved, target,  control) {
    	if(!is.null(control$nodeList)){
    		nodes <- model$expandNodeNames(control$nodeList)
    	} else {
        ## use all nodes in the model if the user does not provide
        ## a list of nodes to be monitored
  			nodes <- model$getNodeNames()
    	}
      nodes <- nodes[nodes != target]

    },
    run = function() {
        model[[target]] <<- model$getLogProb(nodes)  
        
    # copy(from = model, to = mvSaved, 
    # 	row = 1, nodes = target, logProb = TRUE)
   }, 
   methods = list( reset = function () {} )
)
 
logLik_summer <- nimbleFunction(
	name = 'logLik_summer',
	contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
          dataNodes <- model$getNodeNames(dataOnly = TRUE)
          dataNodes <- dataNodes[dataNodes != target]
    },
    run = function() {
        model[[target]] <<- model$getLogProb(dataNodes)  

        copy(from = model, to = mvSaved, 
    	row = 1, nodes = target, logProb = TRUE)

   }, 
   methods = list( reset = function () {})
)
 

# # monitor logProb_sum
# conf$removeSampler("logProb_sum")
# conf$addSampler("logProb_sum", type =  "logProb_summer", 
# 				control = list(nodeList = c("beta", "lambda", "eta")) )

data <- list(y = readRDS("data/simulation_bimodal.rds"))

###------------------------------------------------ ##
## Parametric 2PL - constraints on abilities ----
##------------------------------------------------##
code2PL <- nimbleCode({
  for(i in 1:I) {
    for(j in 1:N) {
      y[j, i] ~ dbern(pi[j, i])
      logit(pi[j, i]) <-  lambda[i]*(eta[j] - beta[i])
    }
  }  
  
  for(i in 1:I) {
    log(lambda[i]) ~ dnorm(0.5, var = 0.5)   
    beta[i] ~ dnorm(0,  var = 3)
  } 
  
  for(j in 1:N) {
    eta[j] ~ dnorm(0, 1)
  }

  # ## dummy nodes to track log porbability and log likelihood
  logProbAll  ~ dnorm(0,1)
  # logProbSum ~ dnorm(0,1)
  # logLik      ~ dnorm(0,1)

})

constants <- list(I= dim(data$y)[2], N = dim(data$y)[1])

inits <- list(beta   = rnorm(constants$I, 0, 1),
              log_lambda = runif(constants$I, -1, 1))

inits$lambda <- exp(inits$log_lambda)

monitors <- c("beta", "lambda")


model <- nimbleModel(code  = code2PL,
					data 			= data,  
					constants	= constants,
					inits 			= inits, 
					calculate 	= FALSE)

# monitors <- c(monitors, "logProbAll", "logProbSum", "logLik")
monitors <- c(monitors, "logProbAll")

mcmcConf <- configureMCMC(model, monitors = monitors)

# mcmcConf$removeSampler("logLik")
# mcmcConf$addSampler("logLik", type =  "logLik_summer")
mcmcConf$removeSampler("logProbAll")
mcmcConf$addSampler("logProbAll", type =  "logProb_summer")

# mcmcConf$addSampler("logProbSum", type =  "logProb_summer", 
# control = list(nodeList = c("beta", "lambda", "eta")))
mcmcConf

mcmc <- buildMCMC(mcmcConf)	

mcmc$run(1)

Cmodel <- compileNimble(model)
Cmcmc <- compileNimble(mcmc, project = model, 
  showCompilerOutput = TRUE)

# Error in eval(call(".Call", nimbleUserNamespace$sessionSpecificDll$populateValueMapAccessorsFromNodeNames,  :
#   VECTOR_ELT() can only be applied to a 'list', not a 'NULL'

res <- runMCMC(Cmcmc, 
		   niter 	= 1000,  
		   nburnin  = 100, 
		   setSeed  = 32)

colnames(res)
res[, "logProbAll"]


res[, "logProS"]
res[, "logLik"]











