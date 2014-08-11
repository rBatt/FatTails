

# ================
# = Load scripts =
# ================
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") #also loads GenSA and DEoptim packages
# source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R")


# =================
# = Tony Function =
# =================
# ======================================================================================
# Series of yearly maxima
# ======================================================================================
yearly.Max <- function(X=NULL, n.per.year=NULL){
  n <- length(X)
  n.year <- n/n.per.year
  
  Y <- array(0,c(n.year,1))
  for(i in 1:n.year){
    Y[i] <- max(X[(1+(i-1)*n.per.year):(i*n.per.year)])
  }
  return(Y)
}

yearly.Max2 <- function(X=NULL, n.per.year=NULL){
  n <- length(X)
  n.year <- n/n.per.year
  
  Y <- array(0,c(n.year,2))
  Y0 <- c()
  x.index <- 1:length(X)
  for(i in 1:n.year){
	t.window <- (1+(i-1)*n.per.year):(i*n.per.year)
	Y0[i] <- x.index[t.window][which.max(X[t.window])]
  }
  Y[,1] <- Y0
  Y[,2] <- X[Y0]
  return(Y)
}

# ===========================
# = First half of Tony code =
# ===========================
B <- function(b,phi) max(1,exp(b*(1+phi)))
S <- function(s,phi) (1+exp(-s*(1+phi)))^(-1)
FF <- function(a,X) (1+a*X)^(-1)

n.per.sample <- 25 # number of time steps per sample
samples.per.year <- 50 # number of samples per year
years <- 30 # number of years

tmax <- n.per.sample * samples.per.year * years
Tmax <- samples.per.year * years

variables <- array(0,c(tmax,4))
samples <- array(0,c(Tmax,3))
samples2 <- array(0,c(Tmax,5), dimnames=list(NULL, c("phi1","phi2","rB","rS","X")))

b <- 0.2
s <- 2
a <- 0.1

sd1 <- 0.25
sd2 <- 0.25

X <- 25
cR <- 1
Time <- 1

X.t <- X
compare <- c()

for(time in 1:tmax){
	
	phi1 <- rnorm(n=1,m=0,sd=sd1)
	phi2 <- rnorm(n=1,m=0,sd=sd2)
	
	#rB <- B(b,phi1 + phi2)
	#rS <- S(s,phi1 - phi2)
	rB <- B(b,phi1)
	rS <- S(s,phi2)
	
	variables[time,] <- c(phi1, phi2, rB, rS)
	
	cR <- cR*rB*rS

	# =================
	# = Added by Ryan =
	# =================
	# lastT <- which((1:time) %% n.per.sample == 0) # what is the last value of X.t that lined up with X?
	# last.XT <- ifelse(length(lastT)==0L, X, X.t[time]) # select the last value of X.t that lined up with X, or, if time < n.per.sample, then just use X
	# X.t[time+1] <- rB*rS*X.t[time]*FF(a, last.XT)^(1/n.per.sample) # Ryan
	X.t[time+1] <- rB*rS*X.t[time]*FF(a, X)^(1/n.per.sample) # Ryan
	# ============
	# = End Ryan =
	# ============
	
	
	if(time %% n.per.sample == 0){
		X <- X * cR * FF(a,X)
		samples[Time,] <- c(phi1, phi2, X)
		samples2[Time,] <- c(phi1, phi2, rB, rS, X)
		
		# =================
		# = Added by Ryan =
		# =================
		compare[Time] <- X.t[time+1] # Ryan
		# ============
		# = End Ryan =
		# ============
		
		cR <- 1
		Time <- Time + 1
	}
}

# =================
# = Added by Ryan =
# =================
# plot(head(samples[,3], 10), type="l"); lines(head(compare, 10), type="l", col="red")
# ============
# = End Ryan =
# ============


# ============================
# = Second half of Tony Code =
# ============================
dev.new()
par(mfrow=c(2,2))

# values of output for year 2
annual.variables <- variables[(n.per.sample * samples.per.year+1):(2*(n.per.sample * samples.per.year)),]
annual.samples <- samples[(samples.per.year+1):(2*samples.per.year),]

# A: within-year environmental variables (phi1 and phi2)
plot((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,1], type='l', col='black', xlab='Time t within a year', ylab='Standardized value')
lines((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,2], col='blue')
points(1:samples.per.year,annual.samples[,1], col='black')
points(1:samples.per.year,annual.samples[,2], col='blue')

# B: within-year population variables (B, S, and X)
plot((1:(n.per.sample * samples.per.year))/n.per.sample, annual.variables[,3], type='l', col='green', xlab='Time t within a year', ylab='Standardized value', ylim=c(0.5,1.5), xlim=c(0, samples.per.year))
lines((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,4], col='red')
# points(1:samples.per.year, 3*annual.samples[,3]/max(annual.samples[,3]), col='black')
par(new=TRUE);plot(1:samples.per.year, annual.samples[,3], col='black', type="p", yaxt="n", xlim=c(0, samples.per.year)); axis(side=4)

# C: among-year environmental variables (phi1 and phi2)
plot(samples[,1],type='l')
Y <- yearly.Max2(samples[,1],n.per.year=samples.per.year)
points(Y)

gev.fit2(xdat=Y[,2], show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

# D: among-year environmental variables (phi1 and phi2)
plot(samples[,3],type='l')
Y <- yearly.Max2(samples[,3],n.per.year=samples.per.year)
points(Y)

gev.fit2(xdat=Y[,2], show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

##############################################
############ Iterate to check    #############
##############################################

nreps <- 100

n.per.sample <- 25 # number of time steps per sample
samples.per.year <- 20 # number of samples per year
years <- 30 # number of years

tmax <- n.per.sample * samples.per.year * years
Tmax <- samples.per.year * years

b <- 0.2
s <- 2
a <- 0.1

sd1 <- 0.2
sd2 <- 0.2

xis <- array(0,c(nreps,2))
for(i in 1:nreps){
	variables <- array(0,c(tmax,4))
	samples <- array(0,c(Tmax,3))
	
	X <- 25
	cR <- 1
	T <- 1
	for(t in 1:tmax){
		
		phi1 <- rnorm(n=1,m=0,sd=sd1)
		phi2 <- rnorm(n=1,m=0,sd=sd2)
		
		#rB <- B(b,phi1 + phi2)
		#rS <- S(s,phi1 - phi2)
		rB <- B(b,phi1)
		rS <- S(s,phi2)
		
		variables[t,] <- c(phi1, phi2, rB, rS)
		
		cR <- cR*rB*rS
		
		if(t %% n.per.sample == 0){
			X <- X * cR * FF(a,X)
			samples[T,] <- c(phi1, phi2, X)
			
			cR <- 1
			T <- T + 1
		}
	}
	
	Y <- yearly.Max(samples[,1],n.per.year=samples.per.year)
	gev.envir <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
	
	Y <- yearly.Max(samples[,3],n.per.year=samples.per.year)
	gev.pop <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
	
	xis[i,] <- c(gev.envir, gev.pop)
	
	show(xis[i,])
}

dev.new()
par(mfrow=c(3,1))
hist(xis[,1],main='Envir', xlim=c(-0.8,0.8), breaks=seq(-0.8,0.8,0.1))
hist(xis[,2],main='Pop', xlim=c(-0.8,0.8), breaks=seq(-0.8,0.8,0.1))
hist(xis[,2] - xis[,1],main='Pop-Envir', xlim=c(-1,1), breaks=seq(-1,1,0.1))

colMeans(xis)
mean(xis[,2] - xis[,1])
cor(xis)


# ===============
# = Ryan Figure =
# ===============
library("evir")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

dev.new(width=3.5, height=7)
cols1 <- rep(rep(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7), each=4), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(8, 8, 8, 9, 9, 9, 0, 0, 10, 10, 10, 10, 0, 0, 11, 11, 11, 11, 7, 7, 7, 7), each=4), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2), ncol=7) # create layout matrix
layout(lmat) # define graphical device layout
# par(mfcol=c(6,2), mar=c(1,2,0.5,0.5), ps=8, cex=1, mgp=c(0.25, 0.0, 0), tcl=0.15)
par(mar=c(0.5,0,0.5,0), oma=c(2, 2, 0.1, 0.1), ps=8, cex=1, mgp=c(0.25, 0.0, 0), tcl=0.15)

	# ==========
	# = part 0 =
	# ==========
# part 0: annual maxima of environmental variables, birth/ survival rates, and population size
samples2.max <- matrix(c(apply(samples2, 2, yearly.Max2, n.per.year=samples.per.year)), nrow=years)
colnames(samples2.max) <- paste(rep(c("phi1","phi2","rB","rS","X"), each=2), rep(c(".i",""), 5), sep="")

	# ============
	# = part 0.5 =
	# ============
# part 0.5: Calculate densities for GEV fit (Panel E)
adat <- c(apply(samples2.max[,c("phi1","phi2","X")], 2, scale)) # adat is 'all data' – both thin and fat maxima time series
ma <- min(adat) # find the smallest value of the maxima
# dSeq <- seq(ma, max(adat)*1.5, by=0.1) # create a sequence of values over which to calculate the density
dSeq <- seq(2, 7, by=0.1) # create a sequence of values over which to calculate the density

phi1.gev <- gev.fit2(scale(samples2.max[,"phi1"]))$mle
phi1GEV <- dgev(dSeq, xi=phi2.gev[3], phi2.gev[1], phi2.gev[2]) # probs for fat time series using GEV fit (package 'evir')
phi1GEV[!is.finite(phi1GEV)] <- 0 

phi2.gev <- gev.fit2(scale(samples2.max[,"phi2"]))$mle
phi2GEV <- dgev(dSeq, xi=phi2.gev[3], phi2.gev[1], phi2.gev[2]) # probs for thin time series
phi2GEV[!is.finite(phi2GEV)] <- 0

X.gev <- gev.fit2(scale(samples2.max[,"X"]))$mle
XGEV <- dgev(dSeq, xi=X.gev[3], X.gev[1], X.gev[2]) # probs for thin time series
XGEV[!is.finite(XGEV)] <- 0

	# ==========
	# = part 1 =
	# ==========
# part 1: time series of all samples of environmental variables, phi1 and phi2
plot(samples2[,"phi1"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi1
points(samples2.max[,"phi1.i"], samples2.max[,"phi1"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1

par(mar=c(1, 0, 0.5, 0))
plot(samples2[,"phi2"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi2
points(samples2.max[,"phi2.i"], samples2.max[,"phi2"], col="red", xpd=TRUE) # add circles for annual maxima of phi2

par(mar=c(1, 0, 1, 0))
plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region
# text(1, 1, labels=quote(atop({phantom() %down% phantom()}, {B(phi1[1])*","~S(phi1[2])})), cex=2)
text(1, 1, labels=quote(phantom() %dbldown% phantom()), cex=3, xpd=TRUE) # add double arrow to visually indicate that part 1 gives rise to part 2

	# ==========
	# = part 2 =
	# ==========
# part 2: time series of all samples of birth and survival rates, B(phi1) and S(phi2)
pt2.ylim <- range(samples2[,c("rB", "rS")])
plot(samples2[,"rB"], type="l", col="gray", xlab="", ylab="", bty="l", ylim=pt2.ylim, xpd=TRUE) # plot full samples of variable phi1
points(samples2.max[,"rB.i"], samples2.max[,"rB"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1

lines(samples2[,"rS"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi2
points(samples2.max[,"rS.i"], samples2.max[,"rS"], col="red", xpd=TRUE) # add circles for annual maxima of phi2

plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region
text(1, 1, labels=quote(phantom() %dbldown% phantom()), cex=3, xpd=TRUE) # add double arrow to visually indicate that part 2 gives rise to part 3

	# ==========
	# = part 3 =
	# ==========
# part 3: time series of all samples of the population size, X
plot(samples2[,"X"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi1
points(samples2.max[,"X.i"], samples2.max[,"X"], col="black", xpd=TRUE) # add circles for annual maxima of phi1


	# ==========
	# = part 4 =
	# ==========
# part 4: Distributions (parent and maxima for parts 1-3)
phi1.xi <- round(phi1.gev[3], 2)
phi2.xi <- round(phi2.gev[3], 2)
X.xi <- round(X.gev[3], 2)
phi1.LabXi <- parse(text=paste("xi", phi1.xi, sep=" = "))
phi2.LabXi <- parse(text=paste("xi", phi2.xi, sep=" = "))
X.LabXi <- parse(text=paste("xi", X.xi, sep=" = "))

par(xaxt="n", yaxt="n")
colorPoly(quants=dSeq, dents=list(phi1GEV,phi2GEV, XGEV), cols=c("blue", "red", "black"), bty="l")
cm0 <- min(dSeq)
text(x=cm0+sign(cm0)*cm0*0.05, y=max(c(phi1GEV,phi2GEV, XGEV))*0.95, "E", font=2)
mtext("density", side=2, line=1.25)
mtext("y", side=1, line=1.25)
cma0 <- max(dSeq)
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.75, bquote(xi~"="~.(phi1.LabXi)), pos=4, col="blue")
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.55, bquote(xi~"="~.(phi2.LabXi)), pos=4, col="red")
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.35, bquote(xi~"="~.(X.LabXi)), pos=4, col="black")





