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
source("ARMAFunctions.R") #also loads GenSA and DEoptim packages
source("fatPlot_Functions.R")

# =============================
# = Define Simulation Options =
# =============================
nPerYear <- 10 # number of observations per "year"
Ps <- 1 # vector of the orders of AR to simulate (e.g., 1:2 would simulate time series that were AR(1), and that were AR(2))
Qs <- 0:1 # vector of MA orders
chooseDists <- c("normal", "lnorm") # vector of distributions – can be "normal", "cauchy", "lnorm", and "t"
nReps <- 5 # number of reps to do for each combination of P, Q, and Distribution
simPars <- expand.grid(P=Ps, Q=Qs, Distribution=chooseDists, Rep=1:nReps, N=c(nPerYear)) # set up combinations of simulation options


# ===============================================
# = Simulate ARMA time series and calculate GEV =
# ===============================================
set.seed(449) # Set seed for random number generator
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

# ===================================================================
# = Grab the thinnest and fattest time series (that are stationary) =
# ===================================================================
fattestI <- which.max(simXiS[,"Xi"])
thinnestI <- which.min(simXiS[,"Xi"])

ffTS <- simXiFull[,fattestI] # Fattest full time series
fmTS <- simXiMax[,fattestI] # Fattest max time series
tfTS <- simXiFull[,thinnestI] # Thinnest full time series
tmTS <- simXiMax[,thinnestI] # Thinnest max time series


# ========================
# = Prepare for plotting =
# ========================
# Calculate densities for GEV fit (Panel E)
adat <- c(tfTS[tmTS], ffTS[fmTS]) # adat is 'all data' – both thin and fat maxima time series
ma <- min(adat) # find the smallest value of the maxima
dSeq <- seq(ma-sign(ma)*1.5*ma, max(adat)*1.1, by=0.1) # create a sequence of values over which to calculate the density
fatGEV <- dgev(dSeq, xi=simXiS[fattestI,"Xi"], simXiS[fattestI,"mu"], simXiS[fattestI,"sig"]) # probs for fat time series using GEV fit (package 'evir')
fatGEV[!is.finite(fatGEV)] <- 0 
thinGEV <- dgev(dSeq, xi=simXiS[thinnestI,"Xi"], simXiS[thinnestI,"mu"], simXiS[thinnestI,"sig"]) # probs for thin time series
thinGEV[!is.finite(thinGEV)] <- 0

# Prepare the Xi values to be plotted on figure (Panel E)
tXi <- round(simXiS[thinnestI,"Xi"], 2)
fXi <- round(simXiS[fattestI,"Xi"], 2)
tLabXi <- parse(text=paste("xi", tXi, sep=" = "))
fLabXi <- parse(text=paste("xi", fXi, sep=" = "))

# ===============================
# = Plot Example w/ reversed xy =
# ===============================
# Set up figure space
dev.new(width=3.5, height=5) # open graphical device
cols1 <- rep(rep(1:3, each=3), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(4,5,3), each=3), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2), ncol=7) # create layout matrix
layout(lmat) # define graphical device layout
par(mar=c(2.5,0.0,0.5,0.1), oma=c(0, 2, 0, 0.25), ps=10, cex=1, mgp=c(1, 0.3, 0), tcl=-0.25, family="Times") # set graphical parameters

# Plot thin-tailed time series
ylim1 <- range(tfTS)*c(1, 1.15)
plot(tfTS, type="l", col="gray", xlab="", xaxt="n", ylab="", ylim=ylim1, bty="l")
text(0.025*length(tfTS), y=max(tfTS)*1.05, "A", font=2)
ats1 <- axTicks(1)
labs1 <- ats1/nPerYear
axis(side=1, labels=labs1, at=ats1)
points(tmTS, tfTS[tmTS], col="blue")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot fat-tailed time series
ylim2 <- range(ffTS)*c(1, 1.15)
plot(ffTS, type="l", col="gray", xlab="", ylab="", xaxt="n", ylim=ylim2, bty="l")
text(0.025*length(ffTS), y=max(ffTS)*1.05, "B", font=2)
ats2 <- axTicks(1)
labs2 <- ats2/nPerYear
axis(side=1, labels=labs2, at=ats2)
points(fmTS, ffTS[fmTS], col="red")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot GEV densities
colorPoly(quants=dSeq, dents=list(thinGEV,fatGEV), cols=c("blue", "red"), bty="l")
cm0 <- min(dSeq)
text(x=cm0+sign(cm0)*cm0*0.05, y=max(c(thinGEV,fatGEV))*0.95, "E", font=2)
mtext("density", side=2, line=1.25)
mtext("y", side=1, line=1.25)
cma0 <- max(dSeq)
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.75, bquote(xi~"="~.(tXi)), pos=4, col="blue")
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.58, bquote(xi~"="~.(fXi)), pos=4, col="red")

# Plot empirical densities for thin-tailed
colorDens(vals=list(tfTS, tfTS[tmTS]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1)
cm1 <- min(tfTS)
text(y=cm1+sign(cm1)*cm1*0.15, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
mtext("density", side=1, line=1.25)

# Plot empirical densities for fat-tailed
cm2 <- min(ffTS)
colorDens(vals=list(ffTS, ffTS[fmTS]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim2)
text(y=cm2+sign(cm2)*cm2*0.15, x=0.25*max(density(ffTS)$y, density(ffTS[fmTS])$y), "D", font=2)
mtext("density", side=1, line=1.25)










