# TPNCC
R code for "Improved prediction of transition probabilities in multi-state models with nested case-control data" by Yen Chang, Anastasia Ivanova, Jason P. Fine, Demetrius Albanes, and Yei Eun Shin. 

## Competing risks
- cr_analysis.R: for predicting transition probabilities in the competing risks model with weight calibration and/or the proportional baselines model 
- cr_functions.R: R functions used in cr_analysis.R

## Illness-death model
- idm_analysis.R: code for predicting transition probabilities in the illness-death model with weight calibration and/or the proportional baselines model
- idm_functions.R: R functions used in idm_analysis.R

.cpp files contain C++ functions sourced by TPNCC_cr_functions.R or TPNCC_idm_functions.R.
