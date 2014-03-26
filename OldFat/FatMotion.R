rm(list=ls())
graphics.off()


#Brownian Motion
Value <- c(0, rep(NA,99999))
Steps <- rnorm(100000, mean=0, sd=1)
for(i in 2:100000){
	Value[i] <- Value[i-1] + Steps[i]
	}

Plump <- c(0, rep(NA,99999))
AutoPlump <- c(0, rep(NA,99999))
Rotund <- c(0, rep(NA,99999))
F1 <- round(runif(100000,0,1),0)
F2 <- round(runif(100000,0,1),0)
F3 <- round(runif(100000,0,1),0)
for(i in 2:100000){
	Plump[i] <-  F1[i]*F2[i]*F3[i]*Steps[i] + F1[i]*F2[i]*Steps[i] + F1[i]*F3[i]*Steps[i] + F2[i]*F3[i]*Steps[i] + F1[i]*Steps[i] + F2[i]*Steps[i] + F3[i]*Steps[i]
	AutoPlump[i] <- Plump[i] + 0.85*Plump[i-1]
	Rotund[i] <- F1[i]*Steps[i] + F2[i]*Steps[i] + F3[i]*Steps[i]
	}

#dev.new(height=8, width=8)
#par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
#plot(1:100000, Value, type="l")

#dev.new(height=8, width=8)
#par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
#hist(Value)

#White Noise
#dev.new(height=8, width=8)
#par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
#hist(Steps)

Fatties <- rcauchy(100000, location=0, scale=1)
#dev.new(height=8, width=8)
#par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
#hist(Fatties)

#dev.new(height=8, width=8)
#par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
#hist(Plump)
load("/Users/Battrd/Documents/School&Work/LabGroup/FatTails/KFM_DOTempData.RData")


ValueExcess <- c()
StepsExcess <- c()
FattyExcess <- c()
PlumpExcess <- c()
AutoPlumpExcess <- c()
RotundExcess <- c()
DOExcess <- c()

mu <- seq(0,10,length.out=50)
for(i in 1:length(mu)){
	ValueExcess[i] <- mean(Value[which(Value>mu[i])]-mu[i])
	StepsExcess[i] <- mean(Steps[which(Steps>mu[i])]-mu[i])
	FattyExcess[i] <- mean(Fatties[which(Fatties>mu[i])]-mu[i])
	PlumpExcess[i] <- mean(Plump[which(Plump>mu[i])]-mu[i])
	AutoPlumpExcess[i] <- mean(AutoPlump[which(AutoPlump>mu[i])]-mu[i])
	RotundExcess[i] <- mean(Rotund[which(Rotund>mu[i])]-mu[i])
}

muDO <- seq(200, 250, length.out=100)
for(i in 1:length(muDO)){
	DOExcess[i] <- mean(Sonde1DOX[which(Sonde1DOX>muDO[i])]-muDO[i])
	}



dev.new(height=8, width=8)
par(family="Times", las=1) #, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, ValueExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, StepsExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, FattyExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, PlumpExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, AutoPlumpExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(mu, RotundExcess)

dev.new(height=8, width=8)
par(family="Times", las=1)#, cex=1.5, mar=c(3,3,2,1), oma=c(1,1,1,1))
plot(muDO, DOExcess)

#Spectrum function from Steve
# Compute spectra of X  from autoregression function
# X is an ordered vector of observations in space or time
# Output is a list which is unpacked as follows:
# LogFreq <- OutList[[1]]

 # End function

LogARspec <- function(Xin, n.per.day){
	X <- Xin
	SpecOut <- spectrum(X, plot=FALSE, method="ar")
	freq.days <- SpecOut$freq*n.per.day
	OutList <- list(log10(freq.days), log10(SpecOut$spec))
	return(OutList)
	
	}
	
	
#load("~/Desktop/KFM_DOTempData.RData")
DOSpec <- LogARspec(Xin=Sonde1DOX, n.per.day=360)
plot(DOSpec[[1]], DOSpec[[2]], type="l")
summary(lm(DOSpec[[2]][-1]~DOSpec[[1]][-1]))

BrownSpec <- LogARspec(Xin=Value, n.per.day=1)
plot(BrownSpec[[1]], BrownSpec[[2]], type="l")
summary(lm(BrownSpec[[2]][-1]~ BrownSpec[[1]][-1]))


DOSpecDiff <- LogARspec(Xin=diff(Sonde1DOX), n.per.day=360)
plot(DOSpecDiff[[1]], DOSpecDiff[[2]], type="l")
summary(lm(DOSpecDiff[[2]][-1]~ DOSpecDiff[[1]][-1]))


BrownSpecDiff <- LogARspec(Xin=Steps, n.per.day=1)
plot(BrownSpecDiff[[1]], BrownSpecDiff[[2]], type="l")
summary(lm(BrownSpecDiff[[2]][-1]~ BrownSpecDiff[[1]][-1]))

