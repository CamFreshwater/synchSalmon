#*******************************************************************************
# suppSim_3ProdOMs_Figs.R
# Date revised: April 25, 2019
# Inputs: recoverySim.R
# Outputs: pdf plots
# Explainer: Takes output data from simulation run in runModel_mainSim.R and
# generates multipanel figures for three productivity OMs instead of two.
#*******************************************************************************

# Check if required packages are installed and run
listOfPackages <- c("plyr", "here",  
                    "reshape2", "tidyverse", "gsl", 
                    "forcats", "ggpubr", "viridis", 
                    "samSim")

# Load packages
newPackages <- listOfPackages[!(listOfPackages %in% 
                                  installed.packages()[ , "Package"])]
if(length(newPackages)) install.packages(newPackages)
lapply(listOfPackages, require, character.only = TRUE)


# Import data for plotting from runModel_mainSim.R
plotDat <- read.csv(here("outputs/generatedData", "summaryTable_3ProdOMs.csv"),
                    stringsAsFactors = F) %>% 
  mutate(synch = fct_relevel(synch, "high", after = Inf),
         sigma = fct_relevel(sigma, "high", after = Inf),
         om = fct_relevel(om, c("Reference Prod.", "Low Prod.", 
                                "Low Prod. Heavy Tails")))


# Set up operating model names
simPar <- read.csv(here("data/sox/fraserOMInputs_varyCorr.csv"), 
                   stringsAsFactors = F)
simParTrim <- subset(simPar, scenario %in% c("lowSig", "medSig", "highSig",
                                             "lowSigLowA", "medSigLowA",
                                             "highSigLowA",
                                             "lowSigLowStudT", "medSigLowStudT",
                                             "highSigLowStudT"))
scenNames <- unique(simParTrim$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simParTrim$species),
                                                sep = "_"))


# Default plotting variable settings
dotSize = 3.25; lineSize = 0.8; legSize = 0.7; axSize = 10; facetSize = 0.95
colPal <- viridis(length(levels(plotDat$synch)), begin = 0, end = 1)
names(colPal) <- levels(plotDat$synch)

#make dataframe of labels to annotate facets
consVars <- c("medRecRY", "stdRecRY", "ppnMixedOpen", "ppnCUUpper")
consYLabs <- c("Agg. Return\nAbundance", "CU-Specific Std.\nReturn Abundance", 
               "Prop. MUs\nAbove Esc. Goal", "Prop. CUs\nAbove Benchmark")
consLabs <- data.frame(om = rep(factor(unique(plotDat$om),
                                       levels = unique(plotDat$om)),
                                each = 4),
                       lab = c("a)", "d)", "g)", "j)", "b)", "e)", "h)", "k)",
                               "c)", "f)", "i)", "l)"),
                       var = rep(factor(consVars, levels = unique(consVars)),
                                 times = 3)
)

consPlots <- lapply(seq_along(consVars), function(i) {
  temp <- plotDat %>%
    filter(var == consVars[i])
  q <- ggplot(temp, aes(x = sigma, y = medn, ymin = lowQ, ymax = highQ,
                        fill = synch)) +
    labs(x = "Component Variance", y = consYLabs[i], 
         fill = "Sim.\nParameter\nValue") +
    geom_pointrange(shape = 21, fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
    scale_fill_manual(name = "Synchrony", values = colPal,
                      labels = c("low" = "Low",
                                 "med" = "Moderate",
                                 "high" = "High")) +
    geom_text(data = consLabs %>% filter(var == consVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
                            hjust = 1, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed", ncol = 4, nrow = 1)
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 2 | i == 3) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 4) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  return(q)
})

catchVars <- c("medCatch", "ppnYrsHighCatch", "stabilityCatch")
catchYLabs <- c("Catch\nAbundance", "Prop. Years Above\nCatch Threshold",
                "Catch\nStability")
catchLabs <- data.frame(om = rep(factor(unique(plotDat$om),
                                        levels = unique(plotDat$om)),
                                 each = 3),
                        lab = c("a)", "d)", "g)", "b)", "e)", "h)", "c)", "f)",
                                "i)"),
                        var = rep(factor(catchVars, levels = unique(catchVars)),
                                  times = 3)
)
catchPlots <- lapply(seq_along(catchVars), function(i) {
  temp <- plotDat %>%
    filter(var == catchVars[i])
  q <- ggplot(temp, aes(x = sigma, y = medn, ymin = lowQ, ymax = highQ,
                        fill = synch)) +
    labs(x = "Component Variance", y = catchYLabs[i], 
         fill = "Sim.\nParameter\nValue") +
    geom_pointrange(shape = 21, fatten = dotSize, size = lineSize,
                    position = position_dodge(width = 0.65)) +
    scale_x_discrete(labels = c("low" = "Low",
                                "med" = "Reference",
                                "high" = "High")) +
    scale_fill_manual(name = "Synchrony", values = colPal,
                      labels = c("low" = "Low",
                                 "med" = "Moderate",
                                 "high" = "High")) +
    geom_text(data = catchLabs %>% filter(var == catchVars[i]),
              mapping = aes(x = 0.75, y = min(temp$lowQ), label = lab,
                            hjust = 1, vjust = 0),
              show.legend = FALSE, inherit.aes = FALSE) +
    facet_wrap(~om, scales = "fixed", ncol = 4, nrow = 1)
  if (i == 1) {
    q <- q + theme_sleekX(position = "top", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 2) {
    q <- q + theme_sleekX(position = "mid", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  if (i == 3) {
    q <- q + theme_sleekX(position = "bottom", legendSize = legSize,
                          axisSize = axSize, facetSize = facetSize)
  }
  return(q)
})

png(file = paste(here(),"/figs/FigS2_consGroupedPlots_3ProdOMs.png", sep = ""),
    height = 5.5, width = 6, units = "in", res = 600)
ggarrange(consPlots[[1]], consPlots[[2]], consPlots[[3]], consPlots[[4]],
          ncol = 1, nrow = 4, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1,1.2))
dev.off()
png(file = paste(here(),"/figs/FigS3_catchGroupedPlots_3ProdOMs.png", sep = ""),
    height = 4.5, width = 6, units = "in", res = 600)
ggarrange(catchPlots[[1]], catchPlots[[2]], catchPlots[[3]],
          ncol = 1, nrow = 3, common.legend = TRUE, legend = "right",
          align = "v", heights = c(1.1,1,1.2))
dev.off()