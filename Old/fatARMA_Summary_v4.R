#RDB
#_v0 (01-Dec-2013): Summarize the ARMA analysis of fatFrame by 1) selecting the best model by AICc; 2) getting the AICc-weighted average of the eigenvalue. Also, this might be a good place to begin comparing to the Fat Tails analysis.

# _v2 (13-Jan-2014): Making some corrections to the computation of the sigmas (changes in tonyARMA_short)

# _v3 (28-Jan-2013): Get the residuals from the ARMA fit. Delete old code that I used to make figures for CFL Seminar.

rm(list=ls())
graphics.off()

library("plyr")
library("rpart")
library("party")
library("RColorBrewer")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("fatARMA_v1.RData") #this is the data file containing the completed ARMA analysis. Note that _v2 is the same as _v1, because tonyARMA_short _v4 and _v5 compute the ARMA the same, but differ in the way they compute the sigma metrics. fatARMA_vX.RData only contains the ARMA fit, not the sigma metrics.
load("fatARMA_noStat_v3.RData") # note that there isn't a 'stationary' version of _v3 ... the change in _v3 was to leave in the linear trend. I am being redundant in renaming the objects and .RData files in addition to the version number to avoid confusion (running this analysis takes so long, it would suck to overwrite something etc).
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
load("finalFrame_v3.RData")

source("FatTails_Functions_v7.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("TonySuggestions/tonyARMA_short_v7.R") #also loads GenSA and DEoptim packages
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")

eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))

grabFatARMA <- function(x)x$sigs
grabResidGEV <- function(x) tryCatch({gev.fit(x$Resids)$mle[c("mu_0","sig_0","sh_0")]}, error=function(cond)rep(NA,3))
grabResidLvl2 <- function(x) lvl2(x$Resids)

# ===================================================================
# = # =============================================================
# = summarize ARMA and combine with fat: Linear Trend removed =
# ============================================================= =
# ===================================================================
sigARMA0 <- dlply(fatARMA, .variables=c("variable", "location", "P", "Q"), .fun=getSE, data=finalFrame, .progress="time")
sigARMA00 <- ldply(sigARMA0, grabFatARMA)

residGEV <- ldply(sigARMA0, grabResidGEV, .progress="time")
names(residGEV)[5:7] <- c("residual_mu_0","residual_sig_0","residual_sh_0")

residLvl2 <- ldply(sigARMA0, grabResidLvl2, .progress="time")
names(residLvl2)[5] <- "Level2_residual"

sigARMA <- merge(residGEV, sigARMA00, all=TRUE)
sigARMA <- merge(sigARMA, residLvl2, all=TRUE)
fatARMA1 <- merge(fatARMA, sigARMA, all=TRUE)
fatARMA2 <- Inf2NA(data.frame(fatARMA1, "Order"=fatARMA1[,"P"]+fatARMA1[,"Q"], "PQ"=paste("(", fatARMA1[,"P"], ",", fatARMA1[,"Q"], ")",sep="")))
fatARMA3 <- ddply(fatARMA2, .variables=c("variable", "location"), .fun=fWeighted)

bestARMA <- ddply(.data=fatARMA3, .variables=c("variable", "location"), .fun=minaicc)
names(bestARMA)[1:2] <- c("Variable", "fitBy")

final0 <- merge(bestARMA,sAP, all=TRUE)
final0[,"Order"] <- as.factor(final0[,"Order"])
final0[,"P"] <- as.factor(final0[,"P"])
final0[,"Q"] <- as.factor(final0[,"Q"])
final0[,"Variable"] <- as.factor(final0[,"Variable"])
final0[,"fitBy"] <- as.factor(final0[,"fitBy"])
final0[,"Type"] <- factor(final0[,"Type"], levels=c("Bio", "Chem", "Phys", "Met", "Cosmic"))

# = Need to look into these calculations to make sure they are correct in _v2 (13-Jan-2014) =
InfE <- final0[,"sigInf"]/final0[,"sigE"]
Einf <- (final0[,"sigE"])^2/(final0[,"sigInf"])^2

final0[,"InfE"] <- InfE
final0[,"Einf"] <- Einf

# = Calculated return time for Level 2 for Residuals of ARMA =
final0[,"Level2_time_residual"] <- apply(final0, 1, lvl_return_res, level=2)
final <- final0[final0[,"N"]>=15,]


# ======================================================================
# = # ================================================================
# = summarize ARMA and combine with fat: NO linear trend removed =
# ================================================================ =
# ======================================================================
sigARMA0_nS <- dlply(fatARMA_noStat, .variables=c("variable", "location", "P", "Q"), .fun=getSE, data=finalFrame, .progress="time")
sigARMA00_nS <- ldply(sigARMA0_nS, grabFatARMA)
residGEV_nS <- ldply(sigARMA0_nS, grabResidGEV, .progress="time")
names(residGEV_nS)[5:7] <- c("residual_mu_0","residual_sig_0","residual_sh_0")
residLvl2_nS <- ldply(sigARMA0_nS, grabResidLvl2, .progress="time")
names(residLvl2_nS)[5] <- "Level2_residual"
sigARMA_nS <- merge(residGEV_nS, sigARMA00_nS, all=TRUE)
sigARMA_nS <- merge(sigARMA_nS, residLvl2_nS, all=TRUE)
fatARMA1_nS <- merge(fatARMA_noStat, sigARMA_nS, all=TRUE)

fatARMA2_nS <- Inf2NA(data.frame(fatARMA1_nS, "Order"=fatARMA1_nS[,"P"]+fatARMA1_nS[,"Q"], "PQ"=paste("(", fatARMA1_nS[,"P"], ",", fatARMA1_nS[,"Q"], ")",sep="")))
fatARMA3_nS <- ddply(fatARMA2_nS, .variables=c("variable", "location"), .fun=fWeighted)
bestARMA_nS <- ddply(.data=fatARMA3_nS, .variables=c("variable", "location"), .fun=minaicc)
names(bestARMA_nS)[1:2] <- c("Variable", "fitBy")

final0_nS <- merge(bestARMA_nS,sAP, all=TRUE)
final0_nS[,"Order"] <- as.factor(final0_nS[,"Order"])
final0_nS[,"P"] <- as.factor(final0_nS[,"P"])
final0_nS[,"Q"] <- as.factor(final0_nS[,"Q"])
final0_nS[,"Variable"] <- as.factor(final0_nS[,"Variable"])
final0_nS[,"fitBy"] <- as.factor(final0_nS[,"fitBy"])
final0_nS[,"Type"] <- factor(final0_nS[,"Type"], levels=c("Bio", "Chem", "Phys", "Met", "Cosmic"))

# Need to look into these calculations to make sure they are correct in _v2 (13-Jan-2014)
InfE_nS <- final0_nS[,"sigInf"]/final0_nS[,"sigE"]
Einf_nS <- (final0_nS[,"sigE"])^2/(final0_nS[,"sigInf"])^2

final0_nS[,"InfE"] <- InfE_nS
final0_nS[,"Einf"] <- Einf_nS

# Calculated return time for Level 2 for Residuals of ARMA
final0_nS[,"Level2_time_residual"] <- apply(final0_nS, 1, lvl_return_res, level=2)
final_nS <- final0_nS[final0_nS[,"N"]>=15,]

# ===============================
# = END BOTH FINAL and FINAL_NS =
# ===============================

# ========================================
# = 5) ARMA Results w/o Detrending First =
# ========================================
plot(final_nS[,"residual_sh_0"], final[,"residual_sh_0"], xlab="Xi of ARMA Residuals, no detrending", ylab="Xi of ARMA Residuals, detrended")
plot(final_nS[,"Level2_residual"], final[,"Level2_residual"], xlab="log10 Return time of 20% over record ARMA residual, no detrending", ylab="log10 Return time of 20% over record ARMA residual, detrended")
plot(final_nS[,"Period"], final[,"Period"], xlab="Period, no detrending", ylab="Period, detrended"); abline(a=0, b=1)
plot(final_nS[,"sigE"], final[,"sigE"], xlab="sigmaE, no detrending", ylab="sigmaE, detrended"); abline(a=0, b=1)
plot(final_nS[,"AICc"], final[,"AICc"], xlab="AICc, no detrending", ylab="AICc, detrended"); abline(a=0, b=1)
plot(final_nS[,"InfE"], final[,"InfE"], xlab="InfE, no detrending", ylab="InfE, detrended"); abline(a=0, b=1)
plot(final_nS[,"b1"], final[,"b1"], xlab="AR(1) Coeff, no detrending", ylab="AR(1) Coeff, detrended"); abline(a=0, b=1)
#Answer: It's about the same. Below I'll provide a summary of how it affects other results

# 1) sigE in the different groups shows a pretty similar pattern
dev.new(width=3.4, height=7); par(mfrow=c(2,1), mar=c(3,3,0.5, 0.5), ps=10, mgp=c(2, 0.5, 0), tcl=-0.4)
boxplot(final[,"sigE"]~final[,"Type"], ylab="detrended version")
boxplot(final_nS[,"sigE"]~final_nS[,"Type"], ylab="not detrending")

#2) The regressions have higher R^2 when detrending is performed. 
# W/o detrending, the Xi of Biology is more similar to the other Types, 
summary(lm(sh_0~residual_sh_0+Type, data=final))
summary(lm(sh_0~residual_sh_0+Type, data=final_nS))
# and sigInf is no longer correlated with the Xi of the time series
summary(lm(sh_0~residual_sh_0+sigInf, data=final))
summary(lm(sh_0~residual_sh_0+sigInf, data=final_nS))

#6) #How often are you surprised if you don't know where you are? (Xi from time series)
#How often are you suprised if you do know where you are? (Xi from ARMA residuals)
# This would be a lot of code to duplicate, so I'll just summarize the differences between detrending and not
# First, when comparing across categories of variables, the return times for breakign residual records is pretty similar between detrended/not detrended
# However, the correlation of return time for breaking the time series record and breaking the residual record becames slightly less correlated for the detrended arma
# Also, the time you would wait to break a residual record has a distribution that is more similar to that of the time series record waiting time -- i.e., it no longer takes a longer amount of time to break the residual record compared to the time series record. This was an interesting result w/ detrending, b/c in that case the longer waiting times to break residual records suggest that the arma models reduce frequency that we're surprised relative to just looking at the distribution of the time series. If the time series isn't detrended before applying the arma model, the waiting time to break a residual record decreases a bit (mode of distribution shifts left). This is probably because the arma model does not perform as well.

# 3, 4, & 7) Same results (no relation) for both detrended and not detrended


# ===================================
# = 1) Explore Reasons for sigma[E] =
# ===================================
# Answer:
dev.new(height=3.4, width=3.4)
par(mar=c(3,3,0.5,0.5), mgp=c(2,0.5,0), tcl=-0.4, ps=10, family="Times")

boxplot(final[,"sigE"]~final[,"Type"]) # biological time series were exposed to higher "environmental" variability (higher error variance from ARMA)
mtext(bquote(sigma[E]), side=2, line=2)
# boxplot(final[,"residual_sh_0"]~final[,"Type"]) # no pattern
# plot(final[,"residual_sh_0"], final[,"sigE"]) # no pattern



# ============================
# = 2) Xi of ARMA Residuals? =
# ============================
dev.new(width=7, height=3.4)
par(mfrow=c(1,2), mar=c(2.5,3,0.5,0.5), mgp=c(1.5,0.5,0), tcl=-0.4, ps=10, family="Times")
boxplot(sh_0~Type, data=final)
mtext(bquote(xi[Time~Series]), side=2, line=1.5)
boxplot(residual_sh_0~Type, data=final)
mtext(bquote(xi[ARMA~Residuals]), side=2, line=1.5)

# ==========================================================
# = 6) Xi of ARMA residuals compared to Xi of time series? =
# ==========================================================
gg_sh_logic <- !is.element(final[,"Type"],c("Met")) & !is.element(final[,"fitBy"],c("Formula", "KE", "LR", "WA")) & !is.na(final[,"N"])
png("FatFigures_v4/xiTs_xiResid.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(final[gg_sh_logic,"residual_sh_0"], final[gg_sh_logic,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~residual_sh_0, data=final[gg_sh_logic,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(xi[ARMA~Residuals]), side=1, line=1.5)
dev.off()


summary(lm(sh_0~residual_sh_0+Type, data=final)) # the Xi of the residuals explains a lot of the Xi in the full time series
summary(lm(sh_0~residual_sh_0+Type + sigInf, data=final)) # the variance of the stationary distribution also explains some of the Xi of the time series when the Xi of the residuals is considered
summary(lm(sh_0~residual_sh_0+sigInf, data=final)) # in fact, after the Xi of residuals and the variance of the stationary distribution are considered, the type of Variable does not explain much of the variability in the Xi of the time series
summary(lm(sigInf~Type, data=final)) # but the type of variable doesn't predict the stationary varince very well ... (_v4 ... yes it does? kinda?)
#boxplot(sigInf~Type, data=final) #there is 1 particularly large stationary varince value
smallLogic <- final[,"sigInf"]<1E4 #so let's see what happens when we remove it
summary(lm(sh_0~residual_sh_0+sigInf, data=final[smallLogic,]))
# dev.new(); par(mfrow=c(2,2)); plot(lm(sh_0~residual_sh_0+log10(sigInf), data=final[smallLogic,]))
# dev.new(); par(mfrow=c(2,2)); plot(lm(sh_0~residual_sh_0+Type, data=final))
#I thnk the diagnostic plots look better when Type is a predictor than when sigInf is a predictor
# Overall I would say that the Xi of the residuals are correlated with the Xi of the time series, with the additional explanation that biological and/or variable time series are also positively correlated with time series Xi.

#How often are you surprised if you don't know where you are? (Xi from time series) Answer:
dev.new(width=7, height=5)
par(mfrow=c(2,2), mar=c(2.5,3,0.5, 0.5), ps=9, mgp=c(1.5,0.5,0), tcl=-0.4, family="Times", cex=1)
boxplot(log10(Level2_time)~Type, data=final, outline=FALSE, ylab=("log10(Years) Surpise Obs"))

#How often are you suprised if you do know where you are? (Xi from ARMA residuals) Answer:
boxplot(log10(Level2_time_residual)~Type, data=final, outline=FALSE, ylab=("log10(Years) Surprise Resid"))

#How does frequency of surprise w/o knowledge of where you compare to that frequency w/ knowledge of where you are? Answer:
xfrm_lvl2TS <- 1/log10(final[,"Level2_time"])
xfrm_lvl2Res <- 1/log10(final[,"Level2_time_residual"])
plot(xfrm_lvl2Res, xfrm_lvl2TS, xlab=("1/log10(Return Time Residual)"), ylab=("1/log10(Return Time Obs)"))
abline(a=0, b=1, col="gray", lty="dashed")
retTime_model <- lm(xfrm_lvl2TS~xfrm_lvl2Res)
abline(retTime_model)
r2_ret_times <- round(summary(retTime_model)$r.squared,2)
legend("topleft", title=paste("R2 =", r2_ret_times), lty=c("solid","dashed"), col=c("black", "gray"), legend=c("regression", "1:1"))
summary(retTime_model)

tsDens <- density(log10(final[,"Level2_time"]), na.rm=TRUE, from=0, to=7)
residDens <- density(log10(final[,"Level2_time_residual"]), na.rm=TRUE, from=0, to=7)
plot(tsDens, main="", xlab="log10(Years until Record + 20%)", ylim=range(c(tsDens$y, residDens$y)), lwd=2, zero.line=FALSE)
lines(residDens, col="red", main="", xlab="", ylab="", lwd=2)
legend("topright", lty="solid", col=c("black", "red"), legend=c("Time Series Obs", "ARMA Residuals"), lwd=2)
segments(x0<-x1<-tsDens$x[which.max(tsDens$y)], y0=0, y1=max(tsDens$y), lty="dashed")
segments(x0<-x1<-residDens$x[which.max(residDens$y)], y0=0, y1=max(residDens$y), col="red", lty="dashed")

# =======================
# = 3) Xi vs. sigE/mean =
# =======================
summary(lm(sh_0~I(sigE/mean), data=final))#significant, but does not explain a lot of variance
summary(lm(sh_0~residual_sh_0+I(sigE/mean), data=final)) # not significant once the Xi of the residuals is considered

# =================================
# = 4.1) Xi vs Time series length =
# =================================
summary(lm(sh_0~N, data=final)) #nope
summary(lm(sh_0~residual_sh_0+N, data=final)) #nope
summary(lm(sh_0~residual_sh_0+N+Type, data=final)) # Nope

# =================================================
# = 4.2) Xi vs Time series length for simulations =
# =================================================
# ????

# =============================================================
# = 7) Xi[residuals]/Xi[time series] vs length of time series =
# =============================================================
relaXi <- final[,"residual_sh_0"]/final[,"sh_0"]
dev.new()
plot(final[,"N"], relaXi, ylim=c(-50, 50), xlim=c(15, 35))
summary(lm(relaXi[relaXi<50]~final[relaXi<50,"N"]))
# Xi of the residuals relative to Xi of the time series is not correlated with the length of the time series (same for Duration instead of N)

myWhite <- rgb(t(col2rgb(col="white")),alpha=100, maxColorValue=256)
myWhite2 <- rgb(t(col2rgb(col="white")),alpha=200, maxColorValue=256)


# =====================
# = Period Histrogram =
# =====================
dev.new(width=3.4, height=3.4)
# png("FatFigures_v3/periodHistogram.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
perHist_breaks <- hist(final[final[,"Period"]<150,"Period"], plot=FALSE)$breaks
hist(final[final[,"Period"]<150,"Period"], xlab="Period", ylab="Frequency", main="", col="white")
hist(final[final[,"Type"]=="Chem" & final[,"Period"]<150 ,"Period"], add=TRUE, col="gray", breaks=perHist_breaks)
hist(final[final[,"Type"]=="Bio" & final[,"Period"]<150,"Period"], add=TRUE, col="black", breaks=perHist_breaks)
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)
# dev.off()

# =================
# = Order Barplot =
# =================
tally <- prop.table(table(final[,c("Type", "Order")]), margin=1)[c("Bio", "Chem","Phys","Met"),]
# dev.new(width=10, height=5)
# par(mfrow=c(2,4), mar=c(3,2,1,0.5), cex=1, ps=8, mgp=c(1,0.5,0))
# Names <- c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
# for(i in 1:length(Names)){
# 	barplot(tally, beside=TRUE, col=brewer.pal(n=4, name=Names[i]), axisnames=TRUE, main=Names[i], legend.text=c("Bio", "Chem","Phys","Met"), args.legend=list(bty="n"))
# }
dev.new(width=3.4, height=3.4)
# png("FatFigures_v3/orderBar.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(mar=c(2.5,2.5,0.5,0.5), cex=1, ps=9, mgp=c(1.5,0.5,0), family="Times")
barplot(tally, beside=TRUE, col=brewer.pal(n=4, name="Set1"), axisnames=TRUE, legend.text=c("Bio", "Chem","Phys","Met"), args.legend=list(bty="n"), xlab="Order", ylab="Proportion")
# dev.off()


# ============================================================================================
# = Below this point there were 900+ lines of code in versions <4 for graphs and regressions =
# ============================================================================================


