res <- readRDS("output/posterior_samples_elaborated/simulation_unimodal/parametric/parametric_IRT_constrainedAbilities.rds")

str(res)

library(mcmcse)
xx <- cbind(res[grepl("lambda", names(res))][[1]],
            res[grepl("beta", names(res))][[1]])

multiESS(xx) 


res1 <- readRDS("output/posterior_samples/simulation_unimodal/parametric/parametric_IRT_constrainedAbilities.rds")


ess_multi <- multiESS(res1$samples$samples) 
