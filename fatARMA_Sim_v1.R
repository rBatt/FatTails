# Simulate ARMA time series of varying order and noise, determine if their Xi values differ.
#_v0 (21-Jan-2014)

rm(list=ls())

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("fatARMA_v1.RData") #this is the data file containing the completed ARMA analysis. Note that _v2 is the same as _v1, because tonyARMA_short _v4 and _v5 compute the ARMA the same, but differ in the way they compute the sigma metrics. fatARMA_vX.RData only contains the ARMA fit, not the sigma metrics.
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
load("finalFrame_v2.RData")

source("FatTails_Functions_v7.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("TonySuggestions/tonyARMA_short_v5.R") #also loads GenSA and DEoptim packages
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")

test <- gev.fit(rnorm(100))$mle["sh_0"]

simPars <- expand.grid(P=0:3, Q=0:3, SD=seq(0.5,50, by=1))

# mySim <- function(model, n, sd){
# 	P <- length(model$ar)
# 	y0 <- rnorm(n=P)
# 	
# 	for(i in (P+1):n){
# 		
# 	}
# 	
# }

myFatSim <- function(x){
	# AR Coefficients
	if(x[,"P"]==0){
		tries0 <- 1
		arC <- NULL
		Lambda <- NA
	}else{
		tries0 <- 0
		repeat({
			# arC0 <- runif(min=-0.99, max=0.99, n=(x[,"P"]-1))
			# rSum <- runif(1, -0.99, 0.99)
			# arC00 <- rSum - sum(arC0)
			# arC <- c(arC0, arC00)
			tries0 <- tries0 + 1
			arC <- runif(-2, 2, n=(x[,"P"]))
			
			if(min(Mod(polyroot(c(1, -c(arC)))))>1){
				# tries <- tries0
				break
			}	
		})
		Lambda <- Eig(arC)
	}
	# if(sum(arC)>=0.99){
	# 	print(arC0)
	# 	print(rSum)
	# 	print(arC00)
	# 	print(arC)
	# 	print(sum(arC))
	# }
	
	# MA Coefficients
	if(x[,"Q"]<=1){
		if(x[,"Q"]==1){
			maC <- runif(min=-0.99, max=0.99, n=1)
		}else{
			maC <- NULL
		}
	}else{
		maC0 <- runif(min=-0.99, max=0.99, n=(x[,"Q"]-1))
		rSum <- runif(1, -0.99, 0.99)
		maC00 <- rSum - sum(maC0)
		maC <- c(maC0, maC00)
	}
	
	simOrder <- c(x[,"P"], 0, x[,"Q"])
	simTS <- c(arima.sim(model=list(order=simOrder, ar=arC, ma=maC), n=200, sd=x[,"SD"]))
	
	Xi <- gev.fit(simTS)$mle["sh_0"]
	Order <- x[,"P"] + x[,"Q"]
	
	data.frame(x, "Xi"=Xi, "tries"=tries0, "Lambda"=Lambda, "Order"=Order)
}


simXi <- ddply(.data=simPars, .variables=c("P","Q","SD"), .fun=myFatSim)

dev.new(height=5, width=5)
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), ps=9, tcl=-0.45, mgp=c(1.5, 0.5, 0))
hist(simXi[,"Xi"], main="")
plot(simXi[,"SD"], simXi[,"Xi"])
plot(simXi[,"Lambda"], simXi[,"Xi"])
plot(simXi[,"Order"], simXi[,"Xi"])

# minroots <- min(Mod(polyroot(c(1, -model$ar))))



# min(Mod(polyroot(c(1, -c(coeffs)))))



