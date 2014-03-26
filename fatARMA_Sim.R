# Simulate ARMA time series of varying order and noise, determine if their Xi values differ.
#_v0 (21-Jan-2014)
#_v2 (27-Jan-2014) Changing arima.sim() to allow for non-stationary simulations. This is desirable to simulate a wider variety of time series, as all of the stationary ARMA time series were thin-tailed or had bounded tails.
rm(list=ls())

library("plyr")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("/Data/fatARMA.RData") #this is the data file containing the completed ARMA analysis. Note that _v2 is the same as _v1, because tonyARMA_short _v4 and _v5 compute the ARMA the same, but differ in the way they compute the sigma metrics. fatARMA_vX.RData only contains the ARMA fit, not the sigma metrics.
load("/Data/All_Params_TurnExtreme_Fat_Data.RData")
load("finalFrame.RData")

source("FatTails_Functions.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("ARMAFunctions.R") #also loads GenSA and DEoptim packages
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")

# simPars <- expand.grid(P=0:3, Q=0:3, SD=seq(1,40, by=3), Distribution=c("normal", "t"), Rep=1:5)
simPars <- expand.grid(P=0:2, Q=0:1, Distribution=c("normal", "t", "cauchy"), Rep=1:5, N=c(50, 100, 200))

which.quantile <- function (x, probs, na.rm = FALSE){
  if (! na.rm & any (is.na (x)))
  return (rep (NA_integer_, length (probs)))

  o <- order (x)
  n <- sum (! is.na (x))
  o <- o [seq_len (n)]

  nppm <- n * probs - 0.5
  j <- floor(nppm)
  h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
  j <- j + h

  j [j == 0] <- 1
  o[j]
}

# =======================================================
# = Modify arima.sim to allow nonstationary simulations =
# =======================================================
arima.sim2 <- function (model, n, rand.gen = rnorm, innov = rand.gen(n, ...), 
    n.start = NA, start.innov = rand.gen(n.start, ...), ...) 
{
    if (!is.list(model)) 
        stop("'model' must be list")
    if (n <= 0L) 
        stop("'n' must be strictly positive")
    p <- length(model$ar)
    if (p) {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
		# CHANGED removed this stop to permit nonstationary simulation
        # if (minroots <= 1) 
        #     stop("'ar' part of model is not stationary")
    }
    q <- length(model$ma)
    if (is.na(n.start)) 
		#CHANGED added the max(), because the logic breaks if minroots is <= 1
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(max(minroots, 1.1))), 
            0)
    if (n.start < p + q) 
        stop("burn-in 'n.start' must be as long as 'ar + ma'")
    d <- 0
    if (!is.null(ord <- model$order)) {
        if (length(ord) != 3L) 
            stop("'model$order' must be of length 3")
        if (p != ord[1L]) 
            stop("inconsistent specification of 'ar' order")
        if (q != ord[3L]) 
            stop("inconsistent specification of 'ma' order")
        d <- ord[2L]
        if (d != round(d) || d < 0) 
            stop("number of differences must be a positive integer")
    }
    if (!missing(start.innov) && length(start.innov) < n.start) 
        stop(sprintf(ngettext(n.start, "'start.innov' is too short: need %d point", 
            "'start.innov' is too short: need %d points"), n.start), 
            domain = NA)
    x <- ts(c(start.innov[seq_len(n.start)], innov[1L:n]), start = 1 - 
        n.start)
    if (length(model$ma)) {
        x <- filter(x, c(1, model$ma), sides = 1L)
        x[seq_along(model$ma)] <- 0
    }
    if (length(model$ar)) 
        x <- filter(x, model$ar, method = "recursive")
    if (n.start > 0) 
        x <- x[-(seq_len(n.start))]
    if (d > 0) 
        x <- diffinv(x, differences = d)
    as.ts(x)
}




myFatSim <- function(x){
	N <- x[,"N"]
	# AR Coefficients
	if(x[,"P"]>=1){
		arC <- runif(-1, 1, n=(x[,"P"]))
		Lambda <- max(abs(Eig(arC)))
		minRoot <- min(Mod(polyroot(c(1, -arC))))
	}else{
		arC <- NULL
		Lambda <- NA
		minRoot <- NA
	}

	# MA Coefficients
	if(x[,"Q"]>=1){
		maC <- runif(min=-1, max=1, n=(x[,"Q"]))
	}else{
		maC <- NULL
	}
	
	simOrder <- c(x[,"P"], 0, x[,"Q"])
	if(x[,"Distribution"]=="normal"){
		simResid <- rnorm(N, 0, 1)
		simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=N, innov=simResid, start.innov=rep(0,1E4)))
		# simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=50, sd=x[,"SD"]))
	}
	if(x[,"Distribution"]=="t"){
		simResid <- rt(N, 5)
		simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=N, innov=simResid, start.innov=rep(0,1E4)))
		# simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=50, rand.gen=rt, df=x[,"SD"]))
	}
	if(x[,"Distribution"]=="cauchy"){
		simResid <- rcauchy(n=N)
		simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=N, innov=simResid, start.innov=rep(0,1E4)))
		# simTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=50, rand.gen=rt, df=x[,"SD"]))
	}
	
	
	# Xi <- tryCatch({gev.fit(simTS)$mle["sh_0"]}, error=function(cond)NA)
	Xi <- gev.fit2(simTS)$mle["sh_0"]
	Xi_resid <- gev.fit2(simResid)$mle["sh_0"]
	Order <- x[,"P"] + x[,"Q"]
	# data.frame(x, "Xi"=Xi, "Lambda"=Lambda, "Order"=Order, "minRoot"=minRoot)
	
	arC2 <- c("AR1"=NA,"AR2"=NA,"AR3"=NA)
	arC2[0:x[,"P"]] <- arC
	maC2 <- c("MA1"=NA,"MA2"=NA,"MA3"=NA)
	maC2[0:x[,"Q"]] <- maC
	
	adf <- data.frame(x, "Order"=Order, "Lambda"=Lambda, "minRoot"=minRoot, t(arC2), t(maC2), "Xi"=Xi, "xiResid"=Xi_resid)
	list("summary"=adf, "data"=simTS)
}

# simXi <- ddply(.data=simPars, .variables=c("P","Q","SD"), .fun=myFatSim)



simXi0 <- dlply(.data=simPars, .variables=c("P","Q","Distribution","Rep", "N"), .fun=myFatSim)
reorgSum <- function(x)x$summary
reorgData <- function(x)x$data

simXi <- ldply(simXi0, .fun=reorgSum)
simXiData <- t(ldply(simXi0, .fun=reorgData)[,-(1:3)])


# dev.new(height=5, width=7)
# par(mfrow=c(2,3), mar=c(3,3,0.5,0.5), ps=9, tcl=-0.45, mgp=c(1.5, 0.5, 0))
# hist(simXi[,"Xi"], main="", ylab=bquote(xi))
# plot(simXi[simXi["Distribution"]=="normal","SD"], simXi[simXi["Distribution"]=="normal","Xi"], xlab="sd", ylab=bquote(xi))
# plot(simXi[simXi["Distribution"]=="t","SD"], simXi[simXi["Distribution"]=="t","Xi"], xlab="df", ylab=bquote(xi))
# plot(simXi[,"Lambda"], simXi[,"Xi"], ylab=bquote(xi), xlab=bquote(lambda))
# plot(log10(simXi[,"minRoot"]), simXi[,"Xi"])
# fatI <- simXi[,"Xi"] > 0
# plot(log10(simXi[fatI,"minRoot"]), simXi[fatI,"Xi"])

# for(i in 1:2){
# 	tDistrib <- c("normal", "t")[i]
# 	mainName <- c("Normal Errors", "Student's T Errors")[i]
# 	tsimXi <- simXi[simXi[,"Distribution"]==tDistrib,]
# 	dev.new(height=5, width=7)
# 	par(mfrow=c(2,3), mar=c(3,3,0.5,0.5), ps=10, tcl=-0.45, mgp=c(1.5, 0.5, 0), oma=c(0.1, 0.1, 1, 0.1))
# 	hist(tsimXi[,"Xi"], main="", ylab=bquote(xi))
# 	plot(tsimXi[,"SD"], tsimXi[,"Xi"], xlab="sd", ylab=bquote(xi))
# 	plot(tsimXi[,"Lambda"], tsimXi[,"Xi"], ylab=bquote(xi), xlab=bquote(lambda))
# 	plot(log10(tsimXi[,"minRoot"]), tsimXi[,"Xi"])
# 	fatI <- tsimXi[,"Xi"] > 0
# 	plot(log10(tsimXi[fatI,"minRoot"]), tsimXi[fatI,"Xi"])
# 	mtext(mainName, line=0, outer=TRUE)
# }




# 
# 
# plotFats <- function(qxi, distrib, ...){
# 	# ind <- simXi[,"Distribution"]==distrib & !is.na(simXi[,"Xi"])
# 	# ind <- simXi[,"Distribution"]==distrib & !is.na(simXi[,"Xi"] & simXi[,"minRoot"]>10)
# 	ind <- simXi[,"Distribution"]==distrib & !is.na(simXi[,"Xi"]) & (simXi[,"minRoot"]>0.9 | is.na(simXi[,"minRoot"]))
# 	
# 	foci <- c(0,0.01, 0.05, 0.1, 0.5, 0.9, 0.95, 0.99, 1)[qxi]
# 	
# 	# indQuo <- which.quantile(1:sum(ind), (qxi-1)/8) #floor(sum(ind)/9)*qxi
# 	# indQuo <- which.quantile(1:sum(ind), (qxi-1)/8) #floor(sum(ind)/9)*qxi
# 	indQuo <- which.quantile(1:sum(ind), foci) #floor(sum(ind)/9)*qxi
# 	closestI <- order(simXi[ind,"Xi"])[indQuo]
# 	closestXi <- round(simXi[ind,][closestI,"Xi"],2)
# 	plot(simXiData[,ind][,closestI], ..., ylab="", xlab="Time", main=bquote(xi==.(closestXi)), type="l")
# }
# dev.new(height=7, width=7)
# par(mfrow=c(3,3), mar=c(3,3,1.5,0.5), ps=9, tcl=-0.45, mgp=c(1.5, 0.5, 0))
# for(i in c(1:9)){
# 	plotFats(i, distrib="normal")
# }
# 
# 
# 
# plotFats2 <- function(qxi, distrib, ...){
# 	ind <- simXi[,"Distribution"]==distrib
# 	closestI <- which.quantile(simXi[ind,"Xi"], qxi, na.rm=TRUE)
# 	closestXi <- round(simXi[ind,][closestI,"Xi"],2)
# 	plot(simXiData[,ind][,closestI], ..., ylab="", xlab="Time", main=bquote(xi==.(closestXi)), type="l")
# }
# dev.new(height=7, width=7)
# par(mfrow=c(3,3), mar=c(3,3,1.5,0.5), ps=9, tcl=-0.45, mgp=c(1.5, 0.5, 0))
# smalls <- c(0,1E-3,1E-2,1E-1)
# for(i in c(smalls, 0.5, 1-rev(smalls))){
# 	plotFats2(i, distrib="normal")
# }


# what were the minimum roots of the ARMA time series for which the Xi calculation failed? All were nonstationary.
# failI <- is.na(simXi[,"Xi"])
# dev.new()
# hist(simXi[failI,"minRoot"])

summary(lm(Xi~xiResid*Lambda, data=simXi))
summary(lm(Xi~xiResid*Lambda, data=simXi[simXi[,"minRoot"]>1,]))

summary(lm(Xi~xiResid+Distribution, data=simXi))

plot(simXi[,"xiResid"], simXi[, "Xi"])
plot(simXi[simXi[,"minRoot"]>1, "xiResid"], simXi[simXi[,"minRoot"]>1, "Xi"]); abline(a=0, b=1)



