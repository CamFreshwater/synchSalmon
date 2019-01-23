#*******************************************************************************
# synchRetroAnalysis_suppShortTS.R
# Date revised: Dec. 5 2018
# Inputs: stock-recruit data csv
# Outputs: figure pdfs
# Explainer: Basically equivalent to synchRetroAnalysis but restricted to short
# time series and without a few exploratory analyses. Produced Fig S1 for ms
#*******************************************************************************

require(here); require(zoo); require(ggplot2); 
require(tidyverse); require(ggpubr); require(samSim)

select <- dplyr::select

## Data clean
# SR data
recDat <- read.csv(here("/data/sox/fraserRecDat.csv"), 
                   stringsAsFactors = FALSE) 
recDat <- recDat %>%
  mutate(prod = (rec/ets),
         logProd = log(rec/ets),
         rollProd = NA,
         spwnRetYr = yr + 4)
for (i in 1:nrow(recDat)) { #add top-fitting SR model
  stkN <- recDat$stk[i]
  if (stkN == 1 | stkN == 2 | stkN == 6 | stkN == 8 | stkN == 9) {
    recDat$model[i] <- "larkin"
  }	else {
    recDat$model[i] <- "ricker"
  } 
}

stkIndex <- unique(recDat$stk)
for(j in seq_along(stkIndex)) {
  d <- subset(recDat, stk == stkIndex[j])
  d$rollProd <- rollmeanr(d$logProd, 4, fill = NA)
  recDat[recDat$stk == stkIndex[j], "rollProd"] <- d[, "rollProd"]
}

# Trim and convert to matrix w/ stocks w/ time series from 1975
recDatTrimS <- recDat %>% 
  dplyr::filter(yr > 1975, !stk == "11", !stk == "815") #more stocks, shorter TS

yrs <- unique(recDatTrimS$yr)
stks <- unique(recDatTrimS$stk)
recMatS <- recDatTrimS %>% 
  select(stk, yr, rec) %>% 
  spread(stk, rec) %>% 
  select(-yr) %>% 
  as.matrix()

saveRDS(recMatS, file = here("outputs", "generatedData", "recMatShort.rds"))


# Trim catch data
catchDat <- read.csv(here("/data/sox/fraserCatchDatTrim.csv"), 
                     stringsAsFactors = FALSE)
catchDatShort <- catchDat %>% 
  filter(yr > min(yrs),
         stk %in% stks) %>% 
  group_by(yr) %>% 
  summarize(catch = sum(totalCatch), 
            ret = sum(totalRun)) %>% 
  mutate(rollCatch = rollmeanr(catch, 4, fill = NA),
         rollRet = rollmeanr(ret, 4, fill = NA)) %>% 
  filter(!is.na(rollRet))

# Trim productivity data
rawDat <- recDatTrimS %>%
  select(stk, yr, logProd, rollProd) %>% 
  filter(!is.na(rollProd)) %>% 
  mutate(stk = as.factor(stk))

# Import TMB model outputs
modOut <- readRDS(here("outputs", "generatedData", "tmbSynchEst_suppShort.rds"))
modOut$broodYr <- rep(yrs[12:36], times = 3)

#note that in Sean's output cv_s is "species" and cv_c is "community" following
#Thib and Connolly 2013
retroCVc <- modOut %>% 
  filter(term == "log_cv_s")  
retroPhi <- modOut %>% 
  filter(term == "logit_phi")
retroCVa <- modOut %>% 
  filter(term == "log_cv_c")


#_________________________________________________________________________
## Plot
stkN <- length(unique(rawDat$stk))
labDat <- data.frame(lab = c("a)", "b)", "c)", "d)", "e)", "f)"),
                     var = c("prod", "rec", "catch", "cv", "synch", "aggCV"))
plotYrs <- c(min(rawDat$yr), max(catchDatShort$yr))

rawProdPlot <- ggplot(rawDat, aes(x = yr, y = rollProd, col = stk)) + 
  labs(x = "", y = "log(Recruits/Spawner)") + 
  ylim(min(rawDat$rollProd), max(rawDat$rollProd)) +
  geom_line(size = 0.75) +
  xlim(plotYrs) +
  scale_color_manual(values = rep("grey", length.out = stkN), guide = FALSE) +
  stat_summary(fun.y = mean, colour = "black", geom = "line", size = 1.25)  +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "prod"),
            mapping = aes(x = min(plotYrs), y = max(rawDat$rollProd), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4)
aggRetPlot <- ggplot(catchDatShort, aes(x = yr, y = rollRet)) +
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "rec"),
            mapping = aes(x = min(plotYrs), y = max(catchDatShort$rollRet), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Return")
aggCatchPlot <- ggplot(catchDatShort, aes(x = yr, y = rollCatch)) + 
  geom_line(size = 1.15) +
  theme_sleekX(position = "top", axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  xlim(plotYrs) +
  geom_text(data = labDat %>% filter(var == "catch"),
            mapping = aes(x = min(plotYrs), 
                          y = max(catchDatShort$rollCatch), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Catch") 
compCVPlot <- ggplot(retroCVc, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.15) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "cv"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroCVc$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Component Variability") 
synchPlot <- ggplot(retroPhi, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "synch"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroPhi$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Synchrony")
agCVPlot <- ggplot(retroCVa, aes(x = broodYr, y = est)) + 
  geom_line(size = 1.25) +
  xlim(plotYrs) +
  geom_ribbon(aes(ymin = lwr, ymax = upr), alpha=0.2) +
  theme_sleekX(axisSize = 12) +
  theme(axis.text.y = element_text(angle = 90, vjust = 1, hjust = 0.5)) +
  geom_text(data = labDat %>% filter(var == "aggCV"),
            mapping = aes(x = min(plotYrs), 
                          y = max(retroCVa$upr, na.rm = TRUE), 
                          label = lab, hjust = 0.25, vjust = 0.5), 
            show.legend = FALSE, inherit.aes = FALSE, size = 4) +
  labs(x = "", y = "Aggregate Variability")

png(here("figs/SFig2New_Retro_ShortTS.png"), height = 4.5, width = 6.5,
    units = "in", res = 300)
ggarrange(rawProdPlot, aggRetPlot, aggCatchPlot, compCVPlot, synchPlot, 
          agCVPlot, nrow = 2, ncol = 3, heights = c(1, 1.1))
dev.off()
