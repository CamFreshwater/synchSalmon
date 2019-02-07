# synchSalmon
## Analyses to explore drivers and consequences of component variability and synchrony using a combination of retrospective analyses and closed-loop simulations

-----

*Authors: C. Freshwater*

*Date: **ONGOING***

-----

### Summary
The `synchSalmon` repository contains the data and models necessary to complete the analysis within *MANUSCRIPT TITLE*. The goal of this study was two-fold. First, we conducted a retrospective analysis using time series of stock recruitment data from the Fraser River sockeye salmon stock aggregate to quantify changes in component variability (mean temporal coefficient of variation within a stock) and synchrony (covariance among stocks). Second, we used a closed-loop simulation model to evaluate how different levels of component variability and synchrony influence the probability of achieving a suite of management objectives. The closed-loop simulation model was developed to address generic questions related to rebuilding and is contained within a separate package `samSim`. 

A summary of relevant files and how to run a simulation are provided below. Most functions are provided in `samSim` and contain relatively detailed documentation (and sometimes functioning examples). Details on the operating model (biological dynamics and fishery interactions) and the management procedures (harvest control rule and assessment process) will be provided in the `samSim` vignette at a future date.

*Note:* development of `samSim` is ongoing and current package versions may no longer be directly compatible with the code provided here. Use the version number (Git release?) associated with the published paper to guarantee compatibility. 

-----

### Files
The repository contains the following directories:

#### data
Includes subdirectory:
- *sox* - contains stock-recruitment data necessary to complete retrospective analysis as well input .csv files used as inputs to the `recoverySim()` function

#### figs
Contains figures (main text and supplementary) included in *MANUSCRIPT TITLE*.

#### outputs
Include