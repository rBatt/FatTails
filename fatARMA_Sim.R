# This script will generate a set of ARMA models (see Options below), then analyze the 'annual' maxima using GEV, then ends w/ a plot

# The simulations simulate a few short ARMA time series.
# Each short series represents dynamics within a year – the maximum of each short ARMA series is an annual maxima
# GEV is fit to the maxima

# I look across all of the ARMA series that were simulated, and pick the 1 w/ the thinnest and the 1 w/ the fattest tail
# I use these two sets of time series to illustrate the concept of our method:
	# Start with a full time series (all of the short ARMA series concatenated) [Panels A&B, gray lines]
	# Take annual maxima of those time series (max of each short series) [Panels A&B, colored circles]
	# The "parent" distribution is the distribution of the full time series (Panels C&D, gray polygon)
	# The empirical distribution taken from the annual maxima (i.e., maxima of parent) are the blue/ red polygons in Panels C&D
	# To get an estimate of tailedness, fit the GEV to the time series of maxima
	# The distributions of the GEV fit to the 2 example time series are shown in Panel E
	# The red distribution is fat-tailed (xi = 0.34), the blue distribution is thin-tailed (xi = -0.42)


# =================
# = Load packages =
# =================
library("evir")
library("plyr")

# =========================
# = Set working directory =
# =========================
setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

# ================
# = Load scripts =
# ================
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") #also loads GenSA and DEoptim packages
# source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R")

# =============================
# = Define Simulation Options =
# =============================
nPerYear <- 10 # number of observations per "year"
Ps <- 1 # vector of the orders of AR to simulate (e.g., 1:2 would simulate time series that were AR(1), and that were AR(2))
Qs <- 0 # vector of MA orders
chooseDists <- c("normal", "rcoin") # vector of distributions – can be "normal", "cauchy", "lnorm", and "t"
nReps <- 50 # number of reps to do for each combination of P, Q, and Distributionrt;
simPars <- expand.grid(P=Ps, Q=Qs, Distribution=chooseDists, Rep=1:nReps, N=c(nPerYear)) # set up combinations of simulation options


# ===============================================
# = Simulate ARMA time series and calculate GEV =
# ===============================================
set.seed(439) # Set seed for random number generator
simXi0 <- dlply(.data=simPars, .variables=c("P","Q","Distribution","Rep", "N"), .fun=myFatSim) # conduct simulations and GEV analysis
reorgSum <- function(x)x$summary # convenience function to pull out the 'summary' output from simXi0
reorgFull <- function(x)x$fullTS # convenience function to pull out the 'fullTS' output from simXi0
reorgMax <- function(x)x$maxTS # convenience function to pull out the 'maxTS' output from simXi0

simXi <- ldply(simXi0, .fun=reorgSum) # pull out summary
simXiFull <- t(ldply(simXi0, .fun=reorgFull)[,-(1:5)]) # pull out full time series and drop useless columns
row.names(simXiFull) <- NULL
simXiMax <- t(ldply(simXi0, .fun=reorgMax)[,-(1:5)]) # pull out time series maxima
row.names(simXiMax) <- NULL

statSim <- simXi[,"Lambda"]<1 # subset to stationary time series
simXiS <- simXi[statSim,]
simXiS[,"critXi"] <- fattestSig(simXiS[,"Xi"], simXiS[,"xi.se"])

# ===================================================================
# = Grab the thinnest and fattest time series (that are stationary) =
# ===================================================================

fDl <- simXiS[,"Distribution"]==chooseDists[2]
fattestI <- which(fDl & simXiS[,"Xi"]==max(simXiS[,"Xi"])) #which.max(simXiS[,"critXi"])

tDl <- simXiS[,"Distribution"]==chooseDists[1]
xiPse <- (abs(simXiS[,"Xi"])+simXiS[,"xi.se"])
thinnestI <- which(tDl & (xiPse==min(xiPse[tDl]))) # which.min(abs(simXiS[,"Xi"])+simXiS[,"xi.se"])

print(simXiS[c(fattestI,thinnestI),])
flush.console()

ffTS <- simXiFull[,fattestI] # Fattest full time series
fmTS <- simXiMax[,fattestI] # Fattest max time series
tfTS <- simXiFull[,thinnestI] # Thinnest full time series
tmTS <- simXiMax[,thinnestI] # Thinnest max time series




ddply(simXiS, "Distribution", function(x)x[which.max(x[,"Xi"]),])
boxplot(Xi~Distribution, data=simXiS)




