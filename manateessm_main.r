#///////////////////////////////////////////////////
#### Fit Poisson-Gaussian SSM on simulated data ####
#///////////////////////////////////////////////////


### // load libraries ----

library(TMB)
library(splines) # for B-splines bs()

# TODO: adapt and simplify scripts for illustration with simulated data

# # system.time(compile("PoiGauAR1SSM_alt.cpp")) # run only once
# dyn.load(dynlib("PoiGauAR1SSM_alt")) # run for every new session
# # ^ v0.3-1 forked from v0.3: piecewise linear temp fx separate from Zmat
# 















