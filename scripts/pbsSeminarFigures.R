#*******************************************************************************
# pbsSeminarFigures.R
# Date revised: April 10, 2019
# Explainer: Tweaks previously generated plots for use in April seminar at
# Pacific Biological Station.
#*******************************************************************************


require(here); require(synchrony); require(zoo); require(ggplot2); 
require(tidyverse); require(viridis); require(ggpubr); require(samSim)

select <- dplyr::select

## Data clean
# SR data
recDat <- read.csv(here("/data/sox/fraserRecDat.csv"), 
                   stringsAsFactors = FALSE) 
recDat <- recDat %>%
  mutate(prod = (rec/ets),
         logProd = log(rec/ets),
         prodZ = NA,
         spwnRetYr = yr + 4)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stkN <- recDat$stk[i]
  if (stkN == 1 | stkN == 2 | stkN == 6 | stkN == 815 | stkN == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}

# Add rolling estimates and lagged EFF abundances for larkin models
recDat$rollProd <- NA
recDat$ets1 <- NA
recDat$ets2 <- NA
recDat$ets3 <- NA
stkIndex <- unique(recDat$stk)
for(j in seq_along(stkIndex)) {
  d <- subset(recDat, stk == stkIndex[j])
  d$rollProd <- rollmeanr(d$logProd, 4, fill = NA)
  for (i in 1:nrow(d)) { #add top-fitting SR model
    d$ets1[i] <- ifelse(i < 1, NA, d$ets[i - 1])
    d$ets2[i] <- ifelse(i < 2, NA, d$ets[i - 2])
    d$ets3[i] <- ifelse(i < 3, NA, d$ets[i - 3])
    d$prodZ <- as.numeric(scale(d$logProd, center = TRUE, scale = TRUE))
  }
  recDat[recDat$stk == stkIndex[j], c("rollProd", "prodZ", "ets1", "ets2", 
                                      "ets3")] <- d[, c("rollProd", "prodZ", 
                                                        "ets1", "ets2", "ets3")]
}
recDatTrim1 <- subset(recDat, !is.na(prod) & !is.na(ets3) & !yr == "2012")

rawDat <- recDatTrim1 %>%
  select(stk, yr, rollProd, prodZ, logProd) %>% 
  # filter(stk %in% longStks$stk) %>% 
  mutate(stk = as.factor(stk))

## Trends in productivity
stkN <- length(unique(rawDat$stk))
plotYrs <- c(min(rawDat$yr), max(rawDat$yr))

png(file = here("figs", "presentationFigs", "prodTrends.png"),
    height = 5, width = 6, units = "in", res = 450)
ggplot(rawDat, aes(x = yr, y = rollProd, col = stk)) + 
  labs(x = "", y = "Productivity (log(Recruits/Spawner))") + 
  ylim(min(rawDat$rollProd),  max(rawDat$rollProd)) +
  geom_line(size = 0.55) +
  xlim(plotYrs) +
  scale_color_manual(values = rep("grey", length.out = stkN), guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.35)  +
  theme_sleekX(axisSize = 13) +
  theme(axis.text.y = element_text(hjust = 0.5))
dev.off()