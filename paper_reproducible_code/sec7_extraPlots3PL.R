##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## last update: June, 22 2021
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
##-----------------------------------------#
## Real data - efficiency comparisons
##-----------------------------------------#
paraTimss  <- read.table(paste0("output/mcmc_time/data_timss/", paraFileName), header = T)
paraTimss2  <- read.table(paste0("output/mcmc_time/data_timss/parametric3PL_efficiency.txt"), header = T)

paraTimss2$ESS_second <- paraTimss2$multiEssItemsAbility/paraTimss2$runningTime
paraTimss$ESS_second  <- paraTimss$multiEssItemsAbility/paraTimss$runningTime

paraTimss$simulation  <- "2PL"
paraTimss2$simulation <- "3PL"

paraTimss2$labels <- gsub("parametric3PL_", "", paraTimss2$fileName)
paraTimss$labels  <- gsub("parametric_", "", paraTimss$fileName)

## match R labels to plot labels
paraTimss2$labels <- labelData[match(paraTimss2$labels, labelData$R_label), ]$plot_label
paraTimss$labels  <- labelData[match(paraTimss$labels, labelData$R_label), ]$plot_label

paraTimss <- droplevels(paraTimss[!is.na(paraTimss$labels), ])

paraTimss2$model <- "Parametric"             
paraTimss$model  <- "Parametric"             

## data frame for plotting
dfParametricEff <- data.frame(rbind(paraTimss2[, c("labels", "ESS_second", "simulation", "model")],
                                    paraTimss[, c("labels", "ESS_second", "simulation","model")]))

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation", "Model")

dfParametricEff$Strategy  <- droplevels(factor(dfParametricEff$Strategy, levels = labelData$plot_label))

colorsParametricData <- labelData$colors[match(levels(dfParametricEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5
##-----------------------------------------#
ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation, ncol=2, scales='fixed') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = colorsParametricData) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "figures/REV_3PL_fig5_data_efficiencies.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Compare bnp efficiencies - data
##-----------------------------------------#

bnpTimss  <- read.table(paste0("output/mcmc_time/data_timss/", bnpFileName), header = T)
bnpTimss2  <- read.table(paste0("output/mcmc_time/data_timss/bnp3PL_efficiency.txt"), header = T)

bnpTimss2$ESS_second <- bnpTimss2$multiEssItemsAbility/bnpTimss2$runningTime
bnpTimss$ESS_second <- bnpTimss$multiEssItemsAbility/bnpTimss$runningTime

paraTimss$simulation  <- "2PL"
paraTimss2$simulation <- "3PL"

bnpTimss2$labels <- gsub("bnp3PL_", "", bnpTimss2$fileName)
bnpTimss$labels <- gsub("bnp_", "", bnpTimss$fileName)

## match R labels to plot labels
bnpTimss2$labels <- labelData[match(bnpTimss2$labels, labelData$R_label), ]$plot_label
bnpTimss$labels  <- labelData[match(bnpTimss$labels, labelData$R_label), ]$plot_label

bnpTimss2$model <- "Semiparametric"             
bnpTimss$model  <- "Semiparametric"             
## data frame for plotting
dfBnpEff <- data.frame(rbind(bnpTimss2[, c("labels", "ESS_second", "simulation", "model")],
                                    bnpTimss[, c("labels", "ESS_second", "simulation","model")]))

colnames(dfBnpEff) <- c("Strategy", "ESS", "Simulation", "Model")

dfBnpEff$Strategy  <- droplevels(factor(dfBnpEff$Strategy, levels = labelData$plot_label))
dfBnpEff <- droplevels(dfBnpEff[-grep("constrained item", dfBnpEff$Strategy), ])

colorsBnp <- labelData$colors[match(levels(dfBnpEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5b
##-----------------------------------------#

ylabel <- 'Effective sample size per second'
title <- paste0("Minimum effective sample size per second (total time)")

p <-  ggplot(dfBnpEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation,ncol=2, scales='free_x') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfBnpEff$Strategy))) +
      scale_fill_manual(values = colorsBnp) 
p

ggsave(filename = "figures/REV_3PL_fig5b_data_efficiencies_bnp.png", plot = p,
        width = plot_width, height = plot_height/6*4 , dpi = 300, units = unit, device='png')
