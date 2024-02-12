#///////////////////////////////////////////////////
#### Fit Poisson-Gaussian SSM on simulated data ####
#///////////////////////////////////////////////////

# rm(list=ls()) # clear global environment

# setwd('.') # set to cloned manateessm GitHub repository


### // load libraries ----

library(TMB)
library(splines) # for B-splines bs()

# TODO: adapt and simplify scripts for illustration with simulated data


### // compile and load model C++ script ----

# TODO

# # system.time(compile("PoiGauAR1SSM_alt.cpp")) # run only once
# dyn.load(dynlib("PoiGauAR1SSM_alt")) # run for every new session
# # ^ v0.3-1 forked from v0.3: piecewise linear temp fx separate from Zmat


### // load data from supplied txt file ----

dat <- read.table('manateessm_data.txt',sep=',',header=T)
str(dat) # 100 rows = 100 months
# ^ month imported as character string by default
dat$ts <- as.Date(paste0(dat$month,'-01')) # fake day 01 for format YYYY-MM-DD
str(dat) # ts = time stamp

nT <- dim(dat)[1] # total number of time points = 100


### // pre-processing: linear time trend ----

cov.timetrend <- (0:(nT-1))/(nT-1)
# ^ whole time window seen as [0,1], with 0 = first month and 1 = last month


### // pre-processing: seasonality as two seasons, winter vs. others ----

ts.month <- as.integer(substr(dat$ts,6,7)) # month as integer, 1-12
# ^ extracts only month from ts
table(ts.month)
# ^ each month appears 8-9 times

cov.seas <- rep(1L,nT) # ini with 1s, 2 seasons dummy code, winter vs. others
for (i in 1:nT){
	if (ts.month[i]%in%c(12,1,2)){ # 12=Dec, 1=Jan, 2=Feb
		cov.seas[i] <- 0L # winter months coded as 0
	}
}
cbind(ts.month, cov.seas)
# ^ winter months set as 0 => winter is the reference in the intercept and the
#   coefficient associated with cov.seas represents the average difference
#   between Mar-Nov and {Jan, Feb, Dec}.


### // pre-processing: piecewise linear effect of temperature ----

thresh.tem <- 24.7
# ^ hard-coded threshold for the two linear pieces

tem.01 <- as.numeric(dat$temperature > thresh.tem)
# ^ months with temperature <= threshold coded as 0, months with temperature >
#   threshold coded as 1

Z.tem <- cbind(
	1-tem.01,                   # intercept tem <= threshold
	dat$temperature*(1-tem.01), # slope tem <= threshold
	tem.01,                     # intercept tem > threshold
	dat$temperature*tem.01      # slope tem > threshold
)
summary(Z.tem)
head(Z.tem)
tail(Z.tem)
# ^ all intercept and slopes, constraints on coefficients hard-coded in cpp


### // setup objectsn for fit ----

# here!!!

