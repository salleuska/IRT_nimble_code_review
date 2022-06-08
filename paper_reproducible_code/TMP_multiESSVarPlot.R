#######################
library(ggplot2)
#######################

## Plots

d1 <- data.frame(mESS = readRDS("output/multiESS_10000.rds"))
d2 <- data.frame(mESS = readRDS("output/multiESS_20000.rds"))
d4 <- data.frame(mESS = readRDS("output/multiESS_40000.rds"))
d5 <- data.frame(mESS = readRDS("output/multiESS_50000.rds"))

d1$Niter <- paste0("N samples = " , 9000, "\n(niter = 10000, 5% burnin)")
d2$Niter <- paste0("N samples = " , 18000, "\n(niter = 20000, 5% burnin)")
d4$Niter <- paste0("N samples = " , 36000, "\n(niter = 40000, 5% burnin)")
d5$Niter <- paste0("N samples = " , 45000, "\n(niter = 50000, 5% burnin)")

df <- rbind(d1, d2, d4, d5)
df$Niter <- factor(df$Niter, levels = levels(as.factor(df$Niter))[c(4,1, 2,3 )])

bp <- ggplot(df, aes(x=Niter, y=mESS, group=Niter)) + 
  geom_boxplot() + theme_bw() +
  geom_point() + 
  theme(axis.text.x=element_blank() )+
  facet_wrap(~Niter, ncol = 2, scales = "free") + 
  ggtitle("NIMBLE model - mESS across 20 replications")
bp
ggsave(bp, file = "output/multiESSNimble.pdf", width = 8, height = 5)

# boxplot(mESS ~ Niter, data = dz)


# boxplot(mESS ~ Niter, data = df)

########
## Stan
s1 <- data.frame(mESS = readRDS("output/Stan_multiESS_10000_warmup_5000.rds"))
s2 <- data.frame(mESS = readRDS("output/Stan_multiESS__warmup_10000.rds"))
s3 <- data.frame(mESS = readRDS("output/Stan_multiESS_30000_warmup_15000.rds"))
s4 <- data.frame(mESS = readRDS("output/Stan_multiESS__warmup_5000.rds"))

s1$Niter <- paste0("N samples (after warmup) = " , 5000, "\n warmup = 5000")
s2$Niter <- paste0("N samples (after warmup) = " , 10000, "\n warmup = 10000")
s3$Niter <- paste0("N samples (after warmup) = " , 15000, "\n warmup = 15000")
s4$Niter <- paste0("N samples (after warmup) = " , 45000, "\n warmup = 5000")


df <- rbind(s1, s2, s3, s4)
df$Niter <- factor(df$Niter, levels = levels(as.factor(df$Niter))[c(4,1, 2, 3)])

bp <- ggplot(df, aes(x=Niter, y=mESS, group=Niter)) + 
  geom_boxplot() + theme_bw() +
  geom_point() + 
  theme(axis.text.x=element_blank() )+
  facet_wrap(~Niter, ncol = 2, scales = "free") + 
  ggtitle("Stan Model - mESS across 20 replications")

bp  
# ggsave(bp, file = "output/multiESSStan.pdf",  width = 8, height = 5)
##################
## Results for one example (1000, 1000)
## STAN
##################
library(rstan)
res <- readRDS("output/Stan_multiESS_res.rds")

names(res)




multiESSVec <- numeric(length(res$samples))
minESSVec <- numeric(length(res$samples))

multiESSVecInd <- numeric(length(res$samples))
minESSVecInd <- numeric(length(res$samples))

essList <- list()
essListStan <- list()
minPar <- numeric(length(res$samples))



for(i in 1:length(res$samples)){
	cat("sample - ", i , "\n")
	modelRes <- res$samples[[i]]
	onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
	              modelRes[grepl("beta", names(modelRes))][[1]])

	itemsAndAbility <- cbind(onlyItems , modelRes$etaSamp)

	itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]


	multiESSVec[i]<- mcmcse::multiESS(itemsAndAbilityMultiESS, 
	                method = "bm",
	                r = 1, 
	                adjust = FALSE)
	essList[[i]] <- coda::effectiveSize(itemsAndAbility)
	minESSVec[i] <- min(essList[[i]])
	minPar[i] <- which.min(essList[[i]])

	essListStan[[i]] <- apply(itemsAndAbilityMultiESS, 2, ess_bulk)

	multiESSVecInd[i]<- mcmcse::multiESS(itemsAndAbilityMultiESS[sample(1:10000, replace = F), ], 
	                method = "bm",
	                r = 1, 
	                adjust = FALSE)
	minESSVecInd[i] <- min(coda::effectiveSize(itemsAndAbility[sample(1:10000, replace = F), ]))

}

# essListStan[[6]]

out <-list(multiESSVec = multiESSVec,
			 minESSVec = minESSVec, 
			 essList = essList, 
			 multiESSVecInd = multiESSVecInd, 
			 minESSVecInd = minESSVecInd)
save(out, file = "output/infoESSStan.rds")
load("output/infoESSStan.rds")

########
xx <- res$params[[1]]




treeDepthMean <- sapply(1:20, function(x) mean(res$params[[x]][[1]][, "treedepth__"]))
plot(treeDepthMean, out$multiESSVec)
plot(treeDepthMean, out$minESSVec)

#sapply(1:20, function(x) sum(res$params[[x]][[1]][, "treedepth__"] == ))

leapfrogMean <- sapply(1:20, function(x) mean(res$params[[x]][[1]][, "n_leapfrog__"]))
plot(leapfrogMean, out$minESSVec,
	xlab = "average of the leapfrog parameters across 20 runs", 
	ylab = "minESS", 
	main = "" )

leapfrogMean <- sapply(1:20, function(x) mean(res$params[[x]][[1]][, "n_leapfrog__"]))
plot(leapfrogMean, out$minESSVec,
	xlab = "average of the leapfrog parameters across 20 runs", 
	ylab = "minESS", 
	main = "" )

stepSizeMean <- sapply(1:20, function(x) mean(res$params[[x]][[1]][, "stepsize__"]))
plot(stepSizeMean, out$minESSVec,
	xlab = "average of the step-size parameter across 20 runs", 
	ylab = "minESS", 
	main = "Stan - 20 runs")

plot(stepSizeMean, minESSVec)



minPar
colnames(itemsAndAbility[, minPar])
minESSVec

# plot(res[[20]]$betaSamp[,10], type = "l")
# coda::effectiveSize(res[[20]]$betaSamp[,10])


#####

load("output/infoESSStan.rds")
stan <- out

out$minESSVec

hist(out$essList[[1]], breaks = 50, main = "stan - low ESS run", 
	xlim = c(2500, 12000) )
hist(out$essList[[9]], breaks = 50, main = "stan - high ESS run", 
xlim = c(8000, 27000))

pdf(file = "stan_mESSvsminESS.pdf", width = 10, height = 5)
plot(out$multiESSVec, out$minESSVec, 
	xlab = 'mESS', ylab = 'minESS', 
	main = 'Stan runs - mESS vs minESS  (20 runs)', 
	pch = 16)
dev.off()

pdf(file = "stan_mESS.pdf", width = 5, height = 4)
hist(out$multiESSVec, breaks = 10,
 main = "stan - mESS distribution (20 runs)")
dev.off()

pdf(file = "stan_minESS.pdf", width = 5, height = 4)
hist(out$minESSVec, breaks = 10, main = "stan - minESS distribution (20 runs)")
dev.off()


##################
## Results for one example (1000, 1000)
## NIMBLE
##################
# res <- readRDS("output/NIMBLE_multiESS50000res.rds")

# multiESSVec <- numeric(length(res))
# minESSVec <- numeric(length(res))

# multiESSVecInd <- numeric(length(res))
# minESSVecInd <- numeric(length(res))

# essList <- list()

# for(i in 1:length(res)){
# 	cat("sample - ", i , "\n")
# 	modelRes <- res[[i]]
# 	onlyItems <- cbind(modelRes[grepl("lambda", names(modelRes))][[1]],
# 	              modelRes[grepl("beta", names(modelRes))][[1]])

# 	itemsAndAbility <- cbind(onlyItems , modelRes$etaSamp)

# 	itemsAndAbilityMultiESS <- itemsAndAbility[, !grepl("(beta\\[1\\])|(lambda\\[1\\])", colnames(itemsAndAbility))]


# 	multiESSVec[i]<- mcmcse::multiESS(itemsAndAbilityMultiESS, 
# 	                method = "bm",
# 	                r = 1, 
# 	                adjust = FALSE)
# 	minESSVec[i] <- min(coda::effectiveSize(itemsAndAbility))
# 	essList[[i]] <- coda::effectiveSize(itemsAndAbility)

# 	multiESSVecInd[i]<- mcmcse::multiESS(itemsAndAbilityMultiESS[sample(1:10000, replace = F), ], 
# 	                method = "bm",
# 	                r = 1, 
# 	                adjust = FALSE)
# 	minESSVecInd[i] <- min(coda::effectiveSize(itemsAndAbility[sample(1:10000, replace = F), ]))

# }

# out <-list(multiESSVec = multiESSVec,
# 			 minESSVec = minESSVec, 
# 			 essList = essList, 
# 			 multiESSVecInd = multiESSVecInd, 
# 			 minESSVecInd = minESSVecInd)
# save(out, file = "output/infoESSNimble.rds")


load("output/infoESSNimble.rds")
nimble <- out
out$minESSVec

cbind(stanWorst = unlist(lapply(stan$essList, function(x) names(which.min(x)))),
nimbleWorst = unlist(lapply(nimble$essList, function(x) names(which.min(x)))))

hist(stan$essList[[1]], breaks = 50, main = "stan - low ESS run",
xlab = "ESS/second")

hist(nimble$essList[[1]])


which(stan$essList[[1]] < 6000)
which(nimble$essList[[1]] < 4000)

which(stan$essList[[1]] < 6000)
hist(stan$essList[[3]])
which(stan$essList[[3]] < 10000)

# # ## nimble 
# # 1175.458
# # ## stan
# # 1228.449
# # hist(nimble$essList[[1]]/1175.458, breaks = 50, main = "nimble ",
# # xlab = "ESS/second")


# hist(stan$essList[[1]]/1228.449, breaks = 50, main = "stan - low ESS run",
# xlab = "ESS/second")

# hist(stan$essList[[9]]/1228.449, breaks = 50, main = "stan - high ESS run", 
# xlab = "ESS/second")

hist(out$essList[[1]], breaks = 50, main = "nimble" )

pdf(file = "nimble_mESSvsminESS.pdf", width = 10, height = 5)
plot(out$multiESSVec, out$minESSVec, 
	xlab = 'mESS', ylab = 'minESS', 
	main = 'NIMBLE runs - mESS vs minESS  (20 runs)', 
	pch = 16)
dev.off()

pdf(file = "nimble_mESS.pdf", width = 5, height = 4)
hist(out$multiESSVec, breaks = 10,
 main = "nimble - mESS distribution (20 runs)")
dev.off()

pdf(file = "nimble_minESS.pdf", width = 5, height = 4)
hist(out$minESSVec, breaks = 10, main = "nimble - minESS distribution (20 runs)")
dev.off()




