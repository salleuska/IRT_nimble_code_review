##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: March 2022
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## nimble version 0.11.1
##-----------------------------------------#
source("R_functions/ggplot_settings.R")
##-----------------------------------------#
## Set dimensions
plot_width  <- 21 ## a4 paper
plot_height <- plot_width/5*2
unit   <- "cm" 

dir <- "output/posterior_samples/mcmc_time/"
paraFileName <- "parametric_efficiency.txt"
bnpFileName <- "bnp_efficiency.txt"
##-----------------------------------------#
## Plot efficiency results for simulation
##-----------------------------------------#
fileList <- list.files("output/mcmc_time", full.names = TRUE)

##-----------------------------------------#
## Unimodal
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", paraFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]

}     

unimodalDf <- as.data.frame(do.call(rbind, unimodalList))

# unique(unimodalDf$fileName)
# unique(labelData$R_label)

# unimodalDf$ESS_second <- unimodalDf$essCodaLogLik/unimodalDf$runningTime
unimodalDf$ESS_second <- unimodalDf$essCodaLogPostItemsAbility/unimodalDf$runningTime
# unimodalDf$ESS_second <- unimodalDf$multiEssItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
levels(dfParametricEff$Simulation) <- c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")


dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

str(dfParametricEff)
dfParametricEff$ESS

ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='free') +
      ylab("min ESS/second (total time)") + xlab("") + 
#      ylim(c(0,6)) + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "unimodalMultiESS.png", plot = p,
        width = 30, height = 30 , 
        dpi = 300, units = unit, device='png')


p1 <- ggplot(dfParametricEff, 
      aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
            legend.position = "bottom")

ggsave(filename = "unimodalMultiESS2.png", plot = p1,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

p2 <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS, group = Strategy, color = Strategy)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = "bottom")

ggsave(filename = "unimodalMultiESS3.png", plot = p2,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", paraFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

unique(bimodalDf$fileName)
unique(labelData$R_label)

bimodalDf$ESS_second <- bimodalDf$multiEssItemsAbility/bimodalDf$runningTime
# bimodalDf$ESS_second <- bimodalDf$essCodaLogLik/bimodalDf$runningTime
# bimodalDf$ESS_second <- bimodalDf$essCodaLogPostItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
levels(dfParametricEff$Simulation) <- c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='free') +
      ylab("min ESS/second (total time)") + xlab("") + 
#      ylim(c(0,6)) + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "bimodalMultiESS.png", plot = p,
        width = 30, height = 30 , 
        dpi = 300, units = unit, device='png')

p1 <- ggplot(dfParametricEff, 
      aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
            legend.position = "bottom")
p1
ggsave(filename = "bimodalMultiESS2.png", plot = p1,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

p2 <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS, group = Strategy, color = Strategy)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = "bottom")
p2
ggsave(filename = "bimodalMultiESS3.png", plot = p2,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$multiEssItemsAbility/multimodalDf$runningTime
# multimodalDf$ESS_second <- multimodalDf$essCodaLogLik/multimodalDf$runningTime
# multimodalDf$ESS_second <- multimodalDf$essCodaLogPostItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
levels(dfParametricEff$Simulation) <- c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

str(dfParametricEff)
dfParametricEff$ESS
ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='free') +
      ylab("min ESS/second (total time)") + xlab("") + 
#      ylim(c(0,6)) + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p



ggsave(filename = "multimodalMultiESS.png", plot = p,
        width = 30, height = 30 , 
        dpi = 300, units = unit, device='png')

p1 <- ggplot(dfParametricEff, 
      aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
            legend.position = "bottom")
p1
ggsave(filename = "multimodalMultiESS2.png", plot = p1,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

p2 <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS, group = Strategy, color = Strategy)) +
  geom_line() + geom_point() +
  ylab("min ESS/second (total time)") + xlab("") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
      legend.position = "bottom")
p2
ggsave(filename = "multimodalMultiESS3.png", plot = p2,
        width = 20, height = 12 , 
        dpi = 300, units = unit, device='png')

##-----------------------------------------#


# xx <- readRDS("output/posterior_samples_elaborated/simulation_unimodal/parametric/parametric_SI_unconstrained.rds")
# plot(xx$otherParSamp[, "myLogLik"], type = "l")
# plot(xx$otherParSamp[, "myLogProbAll"], type = "l")
# plot(xx$otherParSamp[, "myLogProbSome"], type = "l")

# coda::effectiveSize(xx$otherParSamp[, "myLogProbSome"])
# coda::effectiveSize(xx$otherParSamp[, "myLogProbAll"])
# coda::effectiveSize(xx$otherParSamp[, "myLogLik"])





##-----------------------------------------#
## OLD CODE - NOT RUN
##-----------------------------------------#
## Plot Figure 2a
##-----------------------------------------#
ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='fixed') +
      ylab("min ESS/second (total time)") + xlab("") + 
      ylim(c(0,6)) + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p


ggsave(filename = "figures/fig3a_simulation_efficiencies_multiESS.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Plot Figure 2b
## Comparison with sampling times
##-----------------------------------------#
bimodal$ESS_second2 <- bimodal$multiEssItemsAbility/bimodal$samplingTime
unimodal$ESS_second2 <- unimodal$multiEssItemsAbility/unimodal$samplingTime

dfParametricEffSampling <- dfParametricEff
dfParametricEffSampling$ESS <- c(unimodal$ESS_second2, bimodal$ESS_second2)


ylabel <- 'min ESS/second'
title <- paste0("Minimum effective sample size per second (sampling time)")
p <-  ggplot(dfParametricEffSampling,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black", width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='fixed') +
      ylab("min ESS/second (sampling time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricEffSampling$Strategy))) +
      scale_fill_manual(values = labelData$colors) 
p

ggsave(filename = "figures/fig3b_simulation_efficiencies_sampling_multiESS.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Efficiency fof bnp vs parametric - simulation
##-----------------------------------------#

bnpUnimodal <- read.table(paste0("output/mcmc_time/simulation_unimodal/", bnpFileName), header = T)
bnpBimodal  <- read.table(paste0("output/mcmc_time/simulation_bimodal/", bnpFileName), header = T)

bnpUnimodal$ESS_second <- bnpUnimodal$multiEssItemsAbility/bnpUnimodal$runningTime
bnpBimodal$ESS_second <- bnpBimodal$multiEssItemsAbility/bnpBimodal$runningTime

bnpBimodal$simulation <- "Bimodal simulation"
bnpUnimodal$simulation <- "Unimodal simulation"

bnpUnimodal$labels <- gsub("bnp_", "", bnpUnimodal$fileName)
bnpBimodal$labels <- gsub("bnp_", "", bnpBimodal$fileName)

## match R labels to plot labels
bnpUnimodal$labels <- droplevels(labelData[match(bnpUnimodal$labels, labelData$R_label), ]$plot_label)
bnpBimodal$labels  <- droplevels(labelData[match(bnpBimodal$labels, labelData$R_label), ]$plot_label)

## select parametric models with bnp equivalent and create data frame for plotting
dfParametricBnp <- data.frame(rbind(unimodal[which(unimodal$label %in% levels(bnpUnimodal$labels)), 
							c("labels", "ESS_second", "simulation")], 
						   bimodal[which(bimodal$label %in% levels(bnpBimodal$labels)), 
							c("labels", "ESS_second", "simulation")])) 
dfParametricBnp <- droplevels(dfParametricBnp)
dfParametricBnp$model <- "Parametric"

bnpUnimodal$model <- "Semiparametric"						  
bnpBimodal$model  <- "Semiparametric"						  

dfParametricBnp <- rbind(dfParametricBnp, 
                         bnpUnimodal[, c("labels", "ESS_second", "simulation", "model")],
 				                 bnpBimodal[, c("labels", "ESS_second", "simulation", "model")])


colnames(dfParametricBnp) <- c("Strategy", "ESS", "Simulation", "Model")
dfParametricBnp$Strategy  <- droplevels(factor(dfParametricBnp$Strategy, levels = labelData$plot_label))
## set level order
dfParametricBnp$Simulation <- factor(dfParametricBnp$Simulation, levels = c("Unimodal simulation", "Bimodal simulation"))

colorsParametricBnp <- labelData$colors[match(levels(dfParametricBnp$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 3a
##-----------------------------------------#

ylabel <- 'Effective sample size per second'
title <- paste0("Minimum effective sample size per second (total time)")

p <-  ggplot(dfParametricBnp,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation,ncol=2, scales='fixed') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricBnp$Strategy))) +
      scale_fill_manual(values = colorsParametricBnp) 
p

ggsave(filename = "figures/fig4_simulation_efficiencies_bnp_multiESS.png", plot = p,
        width = plot_width, height = plot_height/2*3 , dpi = 300, units = unit, device='png')
