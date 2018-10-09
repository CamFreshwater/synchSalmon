#*************************************************************************************
# distributionPlots.R
# Date revised: Oct 7 2018
# Outputs: figure showing distributions used in different productivity OMs for synch
# analysis
#*************************************************************************************

require(here); require(ggplot2); require(viridis); require(sn); require(tidyverse)
source(here("scripts/synchFunctions.R"))

dat <- data.frame(norm = rst(n = 100000, xi = 0, alpha = 0, nu = 1000, omega = 0.75),
                  sNorm = rst(n = 100000, xi = 0, alpha = log(0.65), nu = 1000, 
                              omega = 0.75),
                  sT = rst(n = 100000, xi = 0, alpha = log(0.65), nu = 3, omega = 0.75)) %>% 
  gather(key = dist, value = values) %>% 
  mutate(dist = as.factor(dist))

ggplot(dat, aes(x = values)) +
  stat_density(aes(group = dist, color = dist), position = "identity", 
               geom = "line") +
  scale_x_continuous(limits = c(-10, 10)) +
  geom_vline(xintercept = 0) +
  theme_sleekX()

ggplot(dat, aes(x = values)) +
  stat_density(aes(group = dist, color = dist), position = "identity", 
               geom = "line") +
  scale_x_continuous(limits = c(-10, 10)) +
  scale_y_continuous(limits = c(0, 0.02)) +
  theme_sleekX()
