#*************************************************************************************
# distributionPlots.R
# Date revised: Oct 7 2018
# Outputs: figure showing distributions used in different productivity OMs for synch
# analysis
#*************************************************************************************

require(here); require(ggplot2); require(viridis); require(sn); require(tidyverse)
source(here("scripts/synchFunctions.R"))
source(here("scripts/synchFunctions.R"))

dat <- data.frame(norm = rst(n = 100000, xi = 0, alpha = 0, nu = Inf, omega = 0.75),
                  sNorm = rst(n = 100000, xi = 0, alpha = log(0.65), nu = Inf, 
                              omega = 0.75),
                  sT = rst(n = 100000, xi = 0, alpha = log(0.65), nu = 3, omega = 0.75)) %>% 
  gather(key = dist, value = values) %>% 
  mutate(dist = as.factor(dist))

dat <- data.frame(norm = rst(n = 100000, xi = 0, alpha = log(0.65), nu = 5, omega = 1),
                  sNorm = rst(n = 100000, xi = 0, alpha = log(0.75), nu = 20, 
                              omega = 1),
                  sT = rst(n = 100000, xi = 0, alpha = log(0.85), nu = Inf, omega = 1)) %>% 
  gather(key = dist, value = values) %>% 
  mutate(dist = as.factor(dist))

colPal <- c("black", "#fd8d3c", "#bd0026")
p <- ggplot(dat, aes(x = values)) +
  geom_density(aes(group = dist, color = dist, fill = dist), position = "identity", 
               alpha = 0.1) +
  scale_x_continuous(limits = c(-5, 5)) +
  labs(x = "Standard Deviations", y = "Probabiliy Density") +
  geom_vline(xintercept = 0, color = "grey30", linetype = 2) +
  scale_color_manual(name = "Distribution", values = colPal,
                     labels = c("norm" = "Normal",
                                "sNorm" = "Skewed Normal",
                                "sT" = "Skewed Student t")) +
  scale_fill_manual(name = "Distribution", values = colPal,
                    labels = c("norm" = "Normal",
                               "sNorm" = "Skewed Normal",
                               "sT" = "Skewed Student t")) +
  theme_sleekX(legendSize = 0.65)

png(here("outputs/distFig.png"), height = 3, width = 5, res = 300, units = "in")
print(p)
dev.off()

  