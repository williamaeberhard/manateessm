#///////////////////////////////////////////////////
#### Fit Poisson-Gaussian SSM on simulated data ####
#///////////////////////////////////////////////////

# rm(list=ls()) # clear global environment

# setwd('.') # set to cloned manateessm GitHub repository


### // load libraries ----

library(TMB) # Template Model Builder, required to compile the C++ script


### // compile and load model C++ script ----

compile("PoiGauAR1SSM.cpp") # run only once
dyn.load(dynlib("PoiGauAR1SSM")) # run for every new session


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


### // setup objects for fit ----
add.zerosurv <- 1e-5 # small arbitrary number to avoid log(0)

Z.cov <- cbind(
	rep(1,nT),               # intercept
	cov.timetrend,           # linear time trend
	dat$surveyor_experience, # surveyor experience dummy
	cov.seas                 # seasonality dummy, 2 season, winter in intercept
)

parlist <- list(
	'beta'=rep(0,dim(Z.cov)[2]), # covariates in lin comb detfx
	'betatemp'=c(0,0),           # tem piecewise lin fx, only 2 free param
	'tphi'=0,                    # AR(1) coef
	'logsigma'=0,                # AR(1) Gaussian noise sd
	'X'=rep(0,nT)                # AR(1) process itself
)
# ^ list with all param initial values, incl random effect X

datalist <- list(
	'obs'=dat$manatee_counts,
	'obsind'=as.integer(!is.na(dat$manatee_counts)), # 0/1 if obs available
	'logoffset'=log(dat$number_surveys + add.zerosurv),
	'Zmat'=Z.cov, # covariates in lin comb detfx
	'Zmattemp'=Z.tem, # temp piecewise lin fx, all 4 col
	'tempthresh'=thresh.tem # temperature threshold junction 2 linear pieces
)
# ^ list with all data and fixed inputs



### // fit ----

obj <- MakeADFun(
	data=datalist,
	parameters=parlist,
	random=c('X'),
	DLL="PoiGauAR1SSM",
	silent=T
)

# obj$fn() # objective function = marginal negative log-likelihood
# obj$gr() # gradient wrt param beta, betatemp, tphi, and logsigma

opt <- nlminb(start=obj$par, obj=obj$fn, gr=obj$gr,
							control=list(eval.max=5000,iter.max=5000))
opt$message # ok



### // param estimates ----

rep <- obj$rep() # report
unlist(rep[c('phi','sigma','inisd')])
# ^ AR(1) estimated param. phi rather close to 0

sdrep <- summary(sdreport(obj)) # report with standard errors
summary.beta <- sdrep[dimnames(sdrep)[[1]]=='beta',][1:dim(Z.cov)[2],]
dimnames(summary.beta)[[1]] <- c(
	'Intercept',
	'Linear time trend',
	'Surveyor experience',
	'Seasonality'
)

summary.betatemp <-  sdrep[dimnames(sdrep)[[1]]=='betatemp',][1:2,]
dimnames(summary.betatemp)[[1]] <- c(
	'Temperature slope <= threshold',
	'Temperature slope > threshold'
)

summary.beta
# ^ estimates and std err for effects in Z.cov

summary.betatemp
# ^ estimates and std err for piecewise linear effect of temperature


### // visualize fit and temperature effect ----
yearticks <- as.Date(paste0(unique(substr(dat$ts,1,4)),'-01-01'))
# ^ year tickmarks/vertical lines on plots

### observation time series and fitted values 
plot(dat$ts,dat$manatee_counts/(dat$number_surveys+add.zerosurv),
		 pch=19,cex=0.6,
		 main='Standardized observations and fitted values',
		 ylab='Counts standardized by number of surveys',xlab='')
abline(v=yearticks,lty=3)
lines(dat$ts,rep$fitted/(dat$number_surveys+add.zerosurv),col='red')
legend('topleft',legend=c('Observations','Fitted values'),
			 col=c(1,'red'),pch=c(19,NA),lty=c(NA,1),pt.cex=c(0.6,NA))


### decomposition of fitted values
plot(dat$ts,log(rep$fitted/(dat$number_surveys+add.zerosurv)),col='red',type='l')
lines(dat$ts,rep$zb,col='blue',lty=2)
lines(dat$ts,rep$X,col='limegreen',lty=2)
# ^ red curve (fitted values) = blue (fixed effects) + green (random effect)
# ^ relative contribution: here fixed effects make up the most part


### piecewise linear effect of temperature
gridtem <- seq(min(dat$temperature), max(dat$temperature), length.out=300)
gridtem.01 <- as.numeric(gridtem > thresh.tem)
Z.gridtem <- cbind(
	1-gridtem.01,                   # intercept tem <= threshold
	gridtem*(1-gridtem.01), # slope tem <= threshold
	gridtem.01,                     # intercept tem > threshold
	gridtem*gridtem.01      # slope tem > threshold
)
linpred.tem <- as.numeric(Z.gridtem%*%rep$betatempvec)
# ^ piecewise linear effect of temperature (marginal)

linpred <- vector('list',4)
linpred[[1]] <- as.numeric(c(1,mean(Z.cov[,2]),0,0)%*%rep$beta) + linpred.tem
linpred[[2]] <- as.numeric(c(1,mean(Z.cov[,2]),1,0)%*%rep$beta) + linpred.tem
linpred[[3]] <- as.numeric(c(1,mean(Z.cov[,2]),0,1)%*%rep$beta) + linpred.tem
linpred[[4]] <- as.numeric(c(1,mean(Z.cov[,2]),1,1)%*%rep$beta) + linpred.tem
# ^ 00 = less experienced, winter
# ^ 10 = experienced, winter
# ^ 01 = less experienced, other seasons
# ^ 11 = experienced, other seasons

colvec.month <- c('red','blue')[cov.seas+1]
# ^ symbol coloring by season: red=winter, blue=others
colvec.curves <- rep(c('red','blue'),each=2)
ltyvec.curves <- rep(2:1,2)

plot(dat$temperature, dat$manatee_counts/(dat$number_surveys+add.zerosurv),
		 main='Estimated temperature effect, by season and surveyor experience',
		 xlab='Monthly max temperature in °C',
		 ylab='Counts standardized by number of surveys',
		 col=colvec.month, # coloring by seasons
		 pch=ifelse(as.logical(dat$surveyor_experience),24,25)
)
for (j in 1:4){
	lines(gridtem, exp(linpred[[j]]),col=colvec.curves[j],lty=ltyvec.curves[j])
}
legend('topleft',c('winter, less experienced','winter, experienced',
									 'spring+fall+summer, less experienced',
									 'spring+fall+summer, experienced'),
			 col=colvec.curves,pch=rep(c(25,24),2),lty=ltyvec.curves
)



# end manateessm_main.r
