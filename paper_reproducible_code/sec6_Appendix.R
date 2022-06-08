##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: March 2022
## R version 4.1.0 (2021-05-18) -- "Camp Pontanezen"
## nimble version 0.11.1
##-----------------------------------------#
source("R_functions/ggplot_settings.R")
library(ggstance)
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
## Unimodal - extraSimulation
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
# unimodalDf$ESS_second <- unimodalDf$essCodaLogPostItemsAbility/unimodalDf$runningTime
unimodalDf$ESS_second <- unimodalDf$multiEssItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

# p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
#       geom_bar(position= position_dodge(),stat='identity',colour = "black",
#        width = 0.8) +
#       facet_wrap(~ Simulation, ncol=2, scales='free') +
#       ylab(yLabel) + xlab("") + 
# #      ylim(c(0,6)) + 
#       theme(legend.position = "none") +
#       coord_flip() +
#       scale_fill_manual(values = labelData$colors[-1]) +
#       scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

# p

# ggsave(filename = "unimodalMultiESS.png", plot = p,
#         width = 30, height = 30 , 
#         dpi = 300, units = unit, device='png')


# p1 <- ggplot(dfParametricEff, 
#       aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
#             legend.position = "bottom")

# p1
# ggsave(filename = "unimodalMultiESS2.png", plot = p1,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')


dfParametricEff$Simulation

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Scenario")

pUni

ggsave(filename = "figures/SM_fig1_unimodalMultiESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

# ##################################
# ### MINIMUM ESS

# dfParametricEff$ESS_second <- unimodalDf$essCodaItemsAbility/unimodalDf$runningTime

# yLabel <- paste0("min ESS/second (total time)")
# p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
#       geom_bar(position= position_dodge(),stat='identity',colour = "black",
#        width = 0.8) +
#       facet_wrap(~ Simulation, ncol=2, scales='free') +
#       ylab(yLabel) + xlab("") + 
# #      ylim(c(0,6)) + 
#       theme(legend.position = "none") +
#       coord_flip() +
#       scale_fill_manual(values = labelData$colors) +
#       scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

# p

# ggsave(filename = "unimodalMinESS.png", plot = p,
#         width = 30, height = 30 , 
#         dpi = 300, units = unit, device='png')


# p1 <- ggplot(dfParametricEff, 
#       aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
#             legend.position = "bottom")

# ggsave(filename = "unimodalMinESS2.png", plot = p1,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')

# p2 <- ggplot(dfParametricEff, 
#       aes(x = Simulation, y= ESS, group = Strategy, color = Strategy)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
#       legend.position = "bottom")

# ggsave(filename = "unimodalMinESS3.png", plot = p2,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')


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

bimodalDf$ESS_second <- bimodalDf$multiEssItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

# p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
#       geom_bar(position= position_dodge(),stat='identity',colour = "black",
#        width = 0.8) +
#       facet_wrap(~ Simulation, ncol=2, scales='free') +
#       ylab(yLabel) + xlab("") + 
# #      ylim(c(0,6)) + 
#       theme(legend.position = "none") +
#       coord_flip() +
#       scale_fill_manual(values = labelData$colors[-1]) +
#       scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

# p

# ggsave(filename = "bimodalMultiESS.png", plot = p,
#         width = 30, height = 30 , 
#         dpi = 300, units = unit, device='png')

# p1 <- ggplot(dfParametricEff, 
#       aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
#             legend.position = "bottom")
# p1
# ggsave(filename = "bimodalMultiESS2.png", plot = p1,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Scenario")

pBi

ggsave(filename = "figures/SM_fig2_bimodalMultiESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

### MINIMUM ESS

# dfParametricEff$ESS_second <- bimodalDf$essCodaItemsAbility/bimodalDf$runningTime


# yLabel <- paste0("min ESS/second (total time)")

# p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
#       geom_bar(position= position_dodge(),stat='identity',colour = "black",
#        width = 0.8) +
#       facet_wrap(~ Simulation, ncol=2, scales='free') +
#       ylab(yLabel) + xlab("") + 
# #      ylim(c(0,6)) + 
#       theme(legend.position = "none") +
#       coord_flip() +
#       scale_fill_manual(values = labelData$colors) +
#       scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

# p

# ggsave(filename = "bimodalMinESS.png", plot = p,
#         width = 30, height = 30 , 
#         dpi = 300, units = unit, device='png')

# p1 <- ggplot(dfParametricEff, 
#       aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
#             legend.position = "bottom")
# p1
# ggsave(filename = "bimodalMinESS2.png", plot = p1,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')

# p2 <- ggplot(dfParametricEff, 
#       aes(x = Simulation, y= ESS, group = Strategy, color = Strategy)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5),
#       legend.position = "bottom")
# p2
# ggsave(filename = "bimodalMinESS3.png", plot = p2,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')

##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal2", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$multiEssItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
# dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])


# yLabel <- paste0("mESS/second (total time)")
# p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
#       geom_bar(position= position_dodge(),stat='identity',colour = "black",
#        width = 0.8) +
#       facet_wrap(~ Simulation, ncol=2, scales='free') +
#       ylab(yLabel) + xlab("") + 
# #      ylim(c(0,6)) + 
#       theme(legend.position = "none") +
#       coord_flip() +
#       scale_fill_manual(values = labelData$colors) +
#       scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

# p



# ggsave(filename = "multimodalMultiESS.png", plot = p,
#         width = 30, height = 30 , 
#         dpi = 300, units = unit, device='png')

# p1 <- ggplot(dfParametricEff, 
#       aes(x = Strategy, y= ESS, group = Simulation, color = Simulation)) +
#   geom_line() + geom_point() +
#   ylab(yLabel) + xlab("") + 
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.5), 
#             legend.position = "bottom")
# p1
# ggsave(filename = "multimodalMultiESS2.png", plot = p1,
#         width = 20, height = 12 , 
#         dpi = 300, units = unit, device='png')

pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Scenario")

pMulti

ggsave(filename = "figures/SM_fig3_multimodalMultiESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#


