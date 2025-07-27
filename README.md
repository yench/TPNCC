# TPNCC
R code for "Prediction of transition probabilities in multi-state models with nested case-control data" by Yen Chang, Anastasia Ivanova, Demetrius Albanes, Jason P. Fine, and Yei Eun Shin. The code accompanies Section 4 and Appendices A-C of the paper, and is designed to help readers understand the proposed methods, and to provide simple working examples with simulated data that readers can modify to analyze real data.  

## Competing risks
- cr_analysis.R: a working example of predicting transition probabilities in the competing risks model with weight calibration and/or the proportional baselines model 
- cr_functions.R: R functions called in cr_analysis.R

## Illness-death model
- idm_analysis.R: a working example of predicting transition probabilities in the illness-death model with weight calibration and/or the proportional baselines model
- idm_functions.R: R functions called in idm_analysis.R

.cpp files contain C++ functions sourced by TPNCC_cr_functions.R or TPNCC_idm_functions.R.
