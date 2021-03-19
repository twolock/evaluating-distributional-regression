# Evaluating distributional regression for modelling partner age distributions
This repository contains the analysis code for [Wolock et al. (2021)](https://arxiv.org/abs/2103.08341). These scripts were designed to be run in parallel on the Imperial College Research Computing Service cluster. This repository contains three scripts:

1. `R/probability_distribution_analysis.R`: the analysis script for the probability distribution analysis.
2. `R/distreg_analysis.R`: the analysis script for the distributional regression evaluation
3. `R/helpers.R`: the functions necessary to construct a number of custom BRMS families

The model lists in `data/` were used to specify the numerous models fit in this analysis.

**NB:** Due to data use restrictions, we only provide the code here, not the data necessary to replicate the analysis. We provide [a separate repository](https://github.com/twolock/distreg-illustration) that illustrates how to fit distributional sinh-arcsinh models in BRMS.