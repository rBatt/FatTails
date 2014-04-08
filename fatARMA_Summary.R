#RDB
#_v0 (01-Dec-2013): Summarize the ARMA analysis of fatFrame by 1) selecting the best model by AICc; 2) getting the AICc-weighted average of the eigenvalue. Also, this might be a good place to begin comparing to the Fat Tails analysis.

# _v2 (13-Jan-2014): Making some corrections to the computation of the sigmas (changes in tonyARMA_short)

# _v3 (28-Jan-2013): Get the residuals from the ARMA fit. Delete old code that I used to make figures for CFL Seminar.

# _v4 (?): I think this version just tried to do the same thing as _v3, but with data that hadn't been made stationary.

#_v5 (23-Feb-2014): fatARMA_Summary will now non-graphically summarize the fatARMA results.I will create a new script to make the fatARMA figures. I am omitting the analysis on the "nonstationary" time series. See FatFrame_v3.R for that analysis, and fatARMA_noStat_v3.RData for the results.

#_v6 (12-Mar-2014): Calculate return times for stable distribution

rm(list=ls())
graphics.off()

library("plyr")
library("rpart")
library("party")
library("RColorBrewer")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("/DatafatARMA.RData") #this is the data file containing the completed ARMA analysis. Note that _v2 is the same as _v1, because tonyARMA_short _v4 and _v5 compute the ARMA the same, but differ in the way they compute the sigma metrics. fatARMA_vX.RData only contains the ARMA fit, not the sigma metrics.
# load("fatARMA_noStat_v3.RData") # note that there isn't a 'stationary' version of _v3 ... the change in _v3 was to leave in the linear trend. I am being redundant in renaming the objects and .RData files in addition to the version number to avoid confusion (running this analysis takes so long, it would suck to overwrite something etc).
# load("/Data/All_Params_TurnExtreme_Fat_Data.RData")
load("/Data/finalFrame.RData")

source("FatTails_Functions.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("ARMAFunctions.R") #also loads GenSA and DEoptim packages
# source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")

eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))


sigARMA0 <- dlply(fatARMA, .variables=c("variable", "location", "P", "Q"), .fun=getSE, data=finalFrame, .progress="time")
grabFatARMA <- function(x)x$sigs
# grabResidXi <- function(x) tryCatch({gev.fit(x$Resids)$mle["sh_0"]}, error=function(cond)NA)
grabResidGEV <- function(x) tryCatch({gev.fit(x$Resids)$mle[c("mu_0","sig_0","sh_0")]}, error=function(cond)rep(NA,3))
grabResidLvl2 <- function(x) lvl2(x$Resids, thresh=setThresh)

sigARMA00 <- ldply(sigARMA0, grabFatARMA)
# residXi <- ldply(sigARMA0, grabResidXi, .progress="time")
residGEV <- ldply(sigARMA0, grabResidGEV, .progress="time")
# names(residXi)[5] <- c("residual_sh_0")
names(residGEV)[5:7] <- c("residual_mu_0","residual_sig_0","residual_sh_0")

residLvl2 <- ldply(sigARMA0, grabResidLvl2, .progress="time")
names(residLvl2)[5] <- "Level2_residual"

# sigARMA <- merge(residXi, sigARMA00, all=TRUE)
sigARMA <- merge(residGEV, sigARMA00, all=TRUE)
sigARMA <- merge(sigARMA, residLvl2, all=TRUE)
fatARMA1 <- merge(fatARMA, sigARMA, all=TRUE)

# metVars <- c("")
# physVars <- c("DaysOpen", "extcoef", "LakeLevel", "max_air_temp", "min_air_temp", "o2", "o2sat", "precip_mm", "range_air_temp", "Secchi")
# chemVars <- c("alk", "brsiuf", "ca", "cl", "cond", "dic", "doc", "drsif", "fe", "k", "mg", "mn", "na", "nh4", "no3no2", "ph",)
# bioVars <- c("avg_length", "avg_zoop_mass", "chlor", "cpue1_Sum", "cpue3_WeiEff", "density")

fatARMA2 <- Inf2NA(data.frame(fatARMA1, "Order"=fatARMA1[,"P"]+fatARMA1[,"Q"], "PQ"=paste("(", fatARMA1[,"P"], ",", fatARMA1[,"Q"], ")",sep="")))
fatARMA3 <- ddply(fatARMA2, .variables=c("variable", "location"), .fun=fWeighted)


bestARMA <- ddply(.data=fatARMA3, .variables=c("variable", "location"), .fun=minaicc)
names(bestARMA)[1:2] <- c("Variable", "fitBy")

# plot(bestARMA[,"Lambda"], bestARMA[,"wLambda"])
sapLog <- sAP[,"Type"]!="Cosmic" & sAP[,"fitBy"]!="Formula"
final0 <- merge(bestARMA,sAP[sapLog,], all=TRUE)
final0[,"Order"] <- as.factor(final0[,"Order"])
final0[,"P"] <- as.factor(final0[,"P"])
final0[,"Q"] <- as.factor(final0[,"Q"])
final0[,"Variable"] <- as.factor(final0[,"Variable"])
final0[,"fitBy"] <- as.factor(final0[,"fitBy"])
final0[,"Type"] <- factor(final0[,"Type"], levels=c("Bio", "Chem", "Phys", "Met"))

# ===========================================================================================
# = Need to look into these calculations to make sure they are correct in _v2 (13-Jan-2014) =
# ===========================================================================================
InfE <- final0[,"sigInf"]/final0[,"sigE"]
Einf <- (final0[,"sigE"])^2/(final0[,"sigInf"])^2

final0[,"InfE"] <- InfE
final0[,"Einf"] <- Einf

# ============================================================
# = Calculated return time for Level 2 for Residuals of ARMA =
# ============================================================
final0[,"Level2_time_residual"] <- apply(final0, 1, lvl_return_res, level=2)

#interesting note, the reason TB alkalinity has NA for its mean and other "lvl" values is b/c it has some negative values, which means the logMean etc couldn't be computed, and I used complete.cases to remove "lvl" info from variables that have any NA's in the "lvl" categories. Actually, I'm not sure thsi could have been the case, because convNeg should take care of this.Ah, but it couldn't have more than 20% of the data set (after NA's and Inf removed) be negative.
# final_logical <- complete.cases(final0[,c("Type","Variable","fitBy","P","Q","Lambda","InfE","Einf")])
# final <- final0[final_logical,]
# sum(!is.finite(fatARMA2[,"Period"]))


# =============================================================
# = Compute return time for Level 2 using Stable Distribution =
# =============================================================
if(runStable){
	source("stableTime.R")
}else{
	load("/Data/stableTime.RData")
}


final <- final0[!is.na(final0[,"N"])&final0[,"N"]>=15,]
save(final, file="/Data/final_fatARMA_Summary.RData")

