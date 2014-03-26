rm(list=ls())
graphics.off()

library("plyr")

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails")
load("All_Params_TurnExtreme_Fat_Data_v7.RData")
load("FatTails_rawData/OrganizedFatData_Read_Fat_Data_v1.RData")
source("FatTails_Functions_v7.R")


# ========================================
# = For the Scenarios class presentation =
# ========================================
# dev.new(width=4, height=4)
# png("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Distributions.png", width=5, height=5, units="in", res=300)
pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Distributions.pdf", width=3.5, height=3.5, bg="white")
par(mar=c(4,4,1,1), ps=10)
plot(dnorm(seq(-6,6,by=0.05)), type="l", ylab="Probability", xlab="X")
lines(dt(seq(-6,6,by=0.05), df=1), type="l", lwd=4)
dev.off()
# boxplot(sh_0~Type, data=All_Params)
# 
# plot(1,1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
# text(x=1, y=1, labels="49 Variables \n from 'Cosmic' to 'Biological' \n 30+ years data & 13 lakes \n 432 values for 'shape' ", cex=1.5)
# 
# SignifShape(x=subset(All_Params, Type=="Phys")[order(subset(All_Params, Type=="Phys")[,"sh_0"]),])
# text(x=50 ,y=-1, "physical", cex=1.5)
# SignifShape(x=subset(All_Params, Type=="Chem")[order(subset(All_Params, Type=="Chem")[,"sh_0"]),])
# text(x=200 ,y=-0.5, "chemical", cex=1.5)
# SignifShape(x=subset(All_Params, Type=="Bio")[order(subset(All_Params, Type=="Bio")[,"sh_0"]),])
# text(x=80 ,y=-0.5, "biological", cex=1.5)





# =====================================================
# = Creat Return-Level plots for common distributions =
# =====================================================
StdShape <- matrix(data=NA, ncol=4, nrow=7, dimnames=list(c("Uniform", "Normal", "T-20", "Exp", "Log-Normal", "Cauchy", "Weibull"), c("N", "mu_0","sig_0", "sh_0")))

RanUnif <- runif(n=1E5)
TailUnif <- RanUnif[which(RanUnif> quantile(RanUnif, 0.90))]
StdShape["Uniform",] <- c("N"=length(TailUnif), gev.fit(TailUnif)$mle)

RanNorm <- rnorm(n=1E5)
TailNorm <- RanNorm[which(RanNorm> quantile(RanNorm, 0.90))]
StdShape["Normal",] <- c("N"=length(TailNorm), gev.fit(TailNorm)$mle)

RanT <- rt(n=1E5, df=15)
TailT <- RanT[which(RanT> quantile(RanT, 0.90))]
StdShape["T-20",] <- c("N"=length(TailT), gev.fit(TailT)$mle)

RanExp <- rexp(n=1E5, rate=1)
TailExp <- RanExp[which(RanExp> quantile(RanExp, 0.90))]
StdShape["Exp",] <- c("N"=length(TailExp), gev.fit(TailExp)$mle)

RanLNorm <- rlnorm(n=1E5)
# TailLNorm <- RanLNorm[which(RanLNorm> quantile(RanLNorm, 0.55))]
StdShape["Log-Normal",] <- c("N"=length(RanLNorm), gev.fit(RanLNorm)$mle)

RanCauchy <- rcauchy(n=1E5)
TailCauchy <- RanCauchy[which(RanCauchy> quantile(RanCauchy, 0.90))]
StdShape["Cauchy",] <- c("N"=length(TailCauchy), gev.fit(TailCauchy)$mle)

RanWeibull <- rweibull(n=1E5, shape=0.17)
TailWeibull <- RanWeibull[which(RanWeibull> quantile(RanWeibull, 0.90))]
StdShape["Weibull",] <- c("N"=length(TailWeibull), gev.fit(TailWeibull)$mle)


RetTimes <- seq(5, 200, by=5)
StdReturns <- matrix(NA, ncol=nrow(StdShape)+1, nrow=length(RetTimes), dimnames=list(NULL, c("RetTime", rownames(StdShape))))
StdReturns[,"RetTime"] <- RetTimes
for(j in 1:nrow(StdShape)){
	for(i in 1:length(RetTimes)){
		RetTime <- RetTimes[i]
		nExts <- StdShape[j,"N"]
		TS_Duration <- StdShape[j,"N"]
		mu <- StdShape[j, "mu_0"]
		sc <- StdShape[j, "sig_0"]
		xi <- mu <- StdShape[j, "sh_0"]
		StdReturns[i,j+1] <- xYr_Lvl(RetTime, nExts, TS_Duration, mu, sc, xi)
	}
}





# ===========================================
# = Create Return-Level Plots for LTER Data =
# ===========================================
LterNames1 <- do.call(paste, Lake_Params3[,c("fitBy", "Variable")])
LterNames <- gsub(" ", "_", LterNames1)
RetTimes <- seq(5, 200, by=5)
LterReturns <- matrix(NA, ncol=nrow(Lake_Params3)+1, nrow=length(RetTimes), dimnames=list(NULL, c("RetTime", LterNames)))
LterReturns[,"RetTime"] <- RetTimes
for(j in 1:nrow(Lake_Params3)){
	for(i in 1:length(RetTimes)){
		RetTime <- RetTimes[i]
		nExts <- Lake_Params3[j,"N"]
		TS_Duration <- Lake_Params3[j,"Duration"]
		mu <- Lake_Params3[j, "mu_0"]
		sc <- Lake_Params3[j, "sig_0"]
		xi <- Lake_Params3[j, "sh_0"]
		LterReturns[i,j+1] <- xYr_Lvl(RetTime, nExts, TS_Duration, mu, sc, xi)
	}
}
MyGray <- rgb(t(col2rgb("gray")), alpha=100, maxColorValue=255)
MyBlue <- rgb(t(col2rgb("blue")), alpha=80, maxColorValue=255)
MyRed <- rgb(t(col2rgb("red")), alpha=80, maxColorValue=255)
Bios <- which(Lake_Params3[,"Type"]=="Bio")
Physs <- which(Lake_Params3[,"Type"]=="Phys")
Chems <- which(Lake_Params3[,"Type"]=="Chem")
LterColors <- c()
LterColors[Bios] <- MyRed
LterColors[Physs] <- MyBlue
LterColors[Chems] <- MyGray


# ===========================================
# = Create a plot combining common and LTER =
# ===========================================
# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatFigures/")
# dev.new(width=3.5, height=6)
# png("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Returns.png", width=3, height=5, units="in", res=300, pointsize=9)
pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Returns.pdf", width=3, height=5, bg="white")
par(mfrow=c(2,1), mar=c(1,2.5,0.5,0.5), oma=c(1.5,0,0,0), ps=8, cex=1)
plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:6){
	par(new=TRUE)
	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
}
legend("topleft", c("Uniform", "Normal", "T (df=20)"), bty="n", inset=c(-0.11, -0.06))
legend("bottomright", c("Exponential", "Log-Normal", "Cauchy", "Weibull"), text.font=c(2,rep(1,3)), bty="n")

plot(LterReturns[,1], LterReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n", col=LterColors[1])
for(i in 1:(ncol(LterReturns)-2)){
	par(new=TRUE)
	plot(LterReturns[,1], LterReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=1, col=LterColors[i+2], bty="n")	
}
mtext("Return Time", side=1, line=-0.5, outer=TRUE)
mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)
legend("bottomright", c("Physical", "Chemical", "Biological"), text.col=c("blue", "gray", "red"))
dev.off()




# =====================
# = Time series plots =
# =====================
#3 panel graph, each panel is a different category
#pick the fattest variable from each category, and plot its time series
# then plot the time series of that variable from the lake for which the variable is the thinnest
# then draw abline() at the level of the N-year return time, where N is the number of years that the time series was observed
MaxChem_Ind <- which(Lake_Params3[,"Type"] == "Chem" & Lake_Params3[,"N"] > 20 )[which.max(Lake_Params3[which(Lake_Params3[,"Type"] == "Chem" & Lake_Params3[,"N"] > 20 ), "sh_0"])]
MaxChem <- Lake_Params3[MaxChem_Ind,]
Min_MaxChem_Ind <- which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"] & Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxChem <- Lake_Params3[Min_MaxChem_Ind,]

MaxPhys_Ind <- which(Lake_Params3[,"Type"] == "Phys"& Lake_Params3[,"N"] > 20)[which.max(Lake_Params3[which(Lake_Params3[,"Type"] == "Phys"& Lake_Params3[,"N"] > 20), "sh_0"])]
MaxPhys <- Lake_Params3[MaxPhys_Ind,]
Min_MaxPhys_Ind <- which(Lake_Params3[,"Variable"]==MaxPhys[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxPhys[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxPhys <- Lake_Params3[Min_MaxPhys_Ind,]

MaxBio_Ind <- which(Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20)[which.max(Lake_Params3[which(Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20), "sh_0"])]
MaxBio <- Lake_Params3[MaxBio_Ind,]
Min_MaxBio_Ind <- which(Lake_Params3[,"Variable"]==MaxBio[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxBio[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxBio <- Lake_Params3[Min_MaxBio_Ind,]

MaxBio_Ind2_Logic <- Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20 & Lake_Params3[,"Variable"] != "CPUE"
MaxBio_Ind2 <- which(MaxBio_Ind2_Logic)[which.max(Lake_Params3[MaxBio_Ind2_Logic, "sh_0"])]
MaxBio2 <- Lake_Params3[MaxBio_Ind2,]
Min_MaxBio_Ind2_Logic <- Lake_Params3[,"Variable"] == MaxBio2[,"Variable"]& Lake_Params3[,"N"] > 20
Min_MaxBio_Ind2 <- which(Min_MaxBio_Ind2_Logic)[which.min(Lake_Params3[which(Min_MaxBio_Ind2_Logic),"sh_0"])]
Min_MaxBio2 <- Lake_Params3[Min_MaxBio_Ind2,]

dateMax <- function(x, name="extcoef"){
	Xr <- x[which.max(x[,name]),"sampledate"]
}

pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/TimeSeries.pdf", width=5, height=7, bg="white")
# dev.new(height=7, width=5)
# pdf(file="EgTimeSeries_RetLvls_09May2013.pdf", width=5, height=7, pointsize=10)
par(mfrow=c(4,2), mar=c(2,2,0.5, 0.5), oma=c(1.25, 1.25, 0, 0), ps=9, cex=1)
# ================================
# = Light Extinction time series =
# ================================
LiExt_MaxTS <- subset(Data_X$LiExt, lakeid==as.character(MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MaxTS <- LiExt_MaxTS[complete.cases(LiExt_MaxTS),]
LiExt_MaxTS_years <- subset(Data_X$LiExt, lakeid==as.character(MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef", "year4"))
LiExt_MaxTS_years <- LiExt_MaxTS_years[complete.cases(LiExt_MaxTS_years),]
LiExt_MaxDates <- ddply(LiExt_MaxTS_years, .variables=c("year4"), .fun=dateMax)[,2]

Low_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"]-MaxPhys[,"se_mu_0"], MaxPhys[,"sig_0"]- MaxPhys[,"se_sig_0"], MaxPhys[,"sh_0"] - MaxPhys[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"]+MaxPhys[,"se_mu_0"], MaxPhys[,"sig_0"]+ MaxPhys[,"se_sig_0"], MaxPhys[,"sh_0"] + MaxPhys[,"se_sh_0"])
# Mean_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"], MaxPhys[,"sig_0"], MaxPhys[,"sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, LiExt_MaxTS[,2]))
plot(LiExt_MaxTS, type="l", col="gray", ylim=Ylim)
points(LiExt_MaxTS[is.element(LiExt_MaxTS[,"sampledate"], LiExt_MaxDates),c("sampledate", "extcoef")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(MaxPhys[,"fitBy"]), "ExtCoef", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxPhys[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(Light~Extinction), side=2, line=2)

LiExt_MinMaxTS <- subset(Data_X$LiExt, lakeid==as.character(Min_MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MinMaxTS <- LiExt_MinMaxTS[complete.cases(LiExt_MinMaxTS),]
LiExt_MinMaxTS_years <- subset(Data_X$LiExt, lakeid==as.character(Min_MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef", "year4"))
LiExt_MinMaxTS_years <- LiExt_MinMaxTS_years[complete.cases(LiExt_MinMaxTS_years),]
LiExt_MinMaxDates <- ddply(LiExt_MinMaxTS_years, .variables=c("year4"), .fun=dateMax)[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"]-Min_MaxPhys[,"se_mu_0"], Min_MaxPhys[,"sig_0"]- Min_MaxPhys[,"se_sig_0"], Min_MaxPhys[,"sh_0"] - Min_MaxPhys[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"]+Min_MaxPhys[,"se_mu_0"], Min_MaxPhys[,"sig_0"]+ Min_MaxPhys[,"se_sig_0"], Min_MaxPhys[,"sh_0"] + Min_MaxPhys[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, LiExt_MinMaxTS[,2]))
plot(LiExt_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(LiExt_MinMaxTS[is.element(LiExt_MinMaxTS[,"sampledate"], LiExt_MinMaxDates),c("sampledate", "extcoef")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(Min_MaxPhys[,"fitBy"]), "ExtCoef", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"], Min_MaxPhys[,"sig_0"], Min_MaxPhys[,"sh_0"])
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxPhys[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))


# ===============================
# = Nitrate/Nitrite time series =
# ===============================
NO3NO2_MaxTS <- subset(Data_X$Chem, lakeid==as.character(MaxChem[,"fitBy"]), select=c("sampledate", "no3no2"))
NO3NO2_MaxTS <- NO3NO2_MaxTS[complete.cases(NO3NO2_MaxTS),]
NO3NO2_MaxTS_years <- subset(Data_X$Chem, lakeid==as.character(MaxChem[,"fitBy"]), select=c("sampledate", "no3no2", "year4"))
NO3NO2_MaxTS_years <- NO3NO2_MaxTS_years[complete.cases(NO3NO2_MaxTS_years),]
NO3NO2_MaxDates <- ddply(NO3NO2_MaxTS_years, .variables=c("year4"), .fun=dateMax, name="no3no2")[,2]

Low_RetLvl <- xYr_Lvl(30, MaxChem[,"N"], MaxChem[,"Duration"], MaxChem[,"mu_0"]-MaxChem[,"se_mu_0"], MaxChem[,"sig_0"]- MaxChem[,"se_sig_0"], MaxChem[,"sh_0"] - MaxChem[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxChem[,"N"], MaxChem[,"Duration"], MaxChem[,"mu_0"]+MaxChem[,"se_mu_0"], MaxChem[,"sig_0"]+ MaxChem[,"se_sig_0"], MaxChem[,"sh_0"] + MaxChem[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, NO3NO2_MaxTS[,2]))
plot(NO3NO2_MaxTS, type="l", col="gray", ylim=Ylim)
points(NO3NO2_MaxTS[is.element(NO3NO2_MaxTS[,"sampledate"], NO3NO2_MaxDates),c("sampledate", "no3no2")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(MaxChem[,"fitBy"]), "NO3NO2", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxChem[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(NO[3] ~ "&" ~NO[2] ~(mu*g~L^-1)), side=2, line=2)


NO3NO2_MinMaxTS <- subset(Data_X$Chem, lakeid==as.character(Min_MaxChem[,"fitBy"]), select=c("sampledate", "no3no2"))
NO3NO2_MinMaxTS <- NO3NO2_MinMaxTS[complete.cases(NO3NO2_MinMaxTS),]
NO3NO2_MinMaxTS_years <- subset(Data_X$Chem, lakeid==as.character(Min_MaxChem[,"fitBy"]), select=c("sampledate", "no3no2", "year4"))
NO3NO2_MinMaxTS_years <- NO3NO2_MinMaxTS_years[complete.cases(NO3NO2_MinMaxTS_years),]
NO3NO2_MinMaxDates <- ddply(NO3NO2_MinMaxTS_years, .variables=c("year4"), .fun=dateMax, name="no3no2")[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxChem[,"N"], Min_MaxChem[,"Duration"], Min_MaxChem[,"mu_0"]-Min_MaxChem[,"se_mu_0"], Min_MaxChem[,"sig_0"]- Min_MaxChem[,"se_sig_0"], Min_MaxChem[,"sh_0"] - Min_MaxChem[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxChem[,"N"], Min_MaxChem[,"Duration"], Min_MaxChem[,"mu_0"]+Min_MaxChem[,"se_mu_0"], Min_MaxChem[,"sig_0"]+ Min_MaxChem[,"se_sig_0"], Min_MaxChem[,"sh_0"] + Min_MaxChem[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, NO3NO2_MinMaxTS[,2]))
plot(NO3NO2_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(NO3NO2_MinMaxTS[is.element(NO3NO2_MinMaxTS[,"sampledate"], NO3NO2_MinMaxDates),c("sampledate", "no3no2")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(Min_MaxChem[,"fitBy"]), "NO3NO2", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxChem[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))


# ============================
# = Chlorophyll Density time series =
# ============================
chlor_MaxTS <- subset(Data_X$Chl, lakeid==as.character(MaxBio2[,"fitBy"]), select=c("sampledate", "chlor"))
chlor_MaxTS <- chlor_MaxTS[complete.cases(chlor_MaxTS),]
chlor_MaxTS_years <- subset(Data_X$Chl, lakeid==as.character(MaxBio2[,"fitBy"]), select=c("sampledate", "chlor", "year4"))
chlor_MaxTS_years <- chlor_MaxTS_years[complete.cases(chlor_MaxTS_years),]
chlor_MaxDates <- ddply(chlor_MaxTS_years, .variables=c("year4"), .fun=dateMax, name="chlor")[,2]

Low_RetLvl <- xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]-MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]- MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] - MaxBio2[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+ MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + MaxBio2[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, chlor_MaxTS[,2]))
plot(chlor_MaxTS, type="l", col="gray", ylim=Ylim)
points(chlor_MaxTS[is.element(chlor_MaxTS[,"sampledate"], chlor_MaxDates),c("sampledate", "chlor")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+1.96*MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+ 1.96*MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + 1.96*MaxBio2[,"se_sh_0"])
xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + MaxBio2[,"se_sh_0"])
abline(h=LterReturns[6, paste(as.character(MaxBio2[,"fitBy"]), "Chla", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxBio2[,"sh_0"],2))), bty="n", inset=c(-.1,0.01))
mtext(expression(Chlorophyll~(mu*g~L^-1)), side=2, line=2)

chlor_MinMaxTS <- subset(Data_X$Chl, lakeid==as.character(Min_MaxBio2[,"fitBy"]), select=c("sampledate", "chlor"))
chlor_MinMaxTS <- chlor_MinMaxTS[complete.cases(chlor_MinMaxTS),]
chlor_MinMaxTS_years <- subset(Data_X$Chl, lakeid==as.character(Min_MaxBio2[,"fitBy"]), select=c("sampledate", "chlor", "year4"))
chlor_MinMaxTS_years <- chlor_MinMaxTS_years[complete.cases(chlor_MinMaxTS_years),]
chlor_MinMaxDates <- ddply(chlor_MinMaxTS_years, .variables=c("year4"), .fun=dateMax, name="chlor")[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]-Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]- Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] - Min_MaxBio2[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]+Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]+ Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] + Min_MaxBio2[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, chlor_MinMaxTS[,2]))
plot(chlor_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(chlor_MinMaxTS[is.element(chlor_MinMaxTS[,"sampledate"], chlor_MinMaxDates),c("sampledate", "chlor")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
# xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]+1.96*Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]+ 1.96*Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] + 1.96*Min_MaxBio2[,"se_sh_0"])
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio2[,"fitBy"]), "Chla", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi == .(round(Min_MaxBio2[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))



# =========================
# = Fish CPUE Time series =
# =========================
CPUE_MaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
CPUE_MaxTS <- CPUE_MaxTS[complete.cases(CPUE_MaxTS),]
Low_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]-MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]- MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] - MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]+MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]+ MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] + MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MaxTS[,2]))
plot(CPUE_MaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MaxTS, col="red", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(Fish~CPUE~(individuals)), side=2, line=2)

CPUE_MinMaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(Min_MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
Low_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]-Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]- Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] - Min_MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]+Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]+ Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] + Min_MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MinMaxTS[,2]))
plot(CPUE_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MinMaxTS, col="blue", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.07))
mtext(expression(Date), side=1, line=0, outer=TRUE)
dev.off()


# ====================================
# = Boxplots and return time scatter =
# ====================================

# ==================================================================
# = Return times predicted from 1) normal, 2) log-normal, & 3) GEV =
# ==================================================================
Inf2NA <- function(x) {x[x==-Inf | x==Inf] <- NA; x}

# dev.new(height=5, width=7)
pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/ScatterReturns.pdf", width=7, height=5, bg="white")
par(mfrow=c(2,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=10, cex=1, tcl=-0.25, mgp=c(3,0.5,0))

#FIRST LEVEL
fl_Ylim <- log10(range(Inf2NA(sAP[,c("Level1_normTime", "Level1_logNormTime", "Level1_time")]), na.rm=TRUE))
# First level, normal
plot(sAP[,"sh_0"], log10(sAP[,"Level1_normTime"]), xlab="", ylab="", ylim=fl_Ylim)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
mtext(bquote(.("[")~mu+2*sigma~.("]")), side=2, line=1.5)
mtext("Predicted by\nNormal", side=3, font=2)
# First Level, log-normal
plot(sAP[,"sh_0"], log10(sAP[,"Level1_logNormTime"]), xlab="", ylab="", ylim=fl_Ylim)
mtext("Predicted by\nLog-Normal", side=3, font=2)
# First level, GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level1_time"]), xlab="", ylab="", ylim=fl_Ylim)
mtext("Predicted by\nGEV", side=3, font=2)

#SECOND LEVEL
sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
# second level, normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
lT <- (setThresh*100)-100
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
dev.off()

# ========================
# = BLANK Scatter Return =
# ========================
pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/ScatterReturns_Blank.pdf", width=7, height=5, bg="white")
par(mfrow=c(2,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=10, cex=1, tcl=-0.25, mgp=c(3,0.5,0))

#FIRST LEVEL
fl_Ylim <- log10(range(Inf2NA(sAP[,c("Level1_normTime", "Level1_logNormTime", "Level1_time")]), na.rm=TRUE))
# First level, normal
plot(sAP[,"sh_0"], log10(sAP[,"Level1_normTime"]), xlab="", ylab="", ylim=fl_Ylim, pch=NA)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
mtext(bquote(.("[")~mu+2*sigma~.("]")), side=2, line=1.5)
mtext("Predicted by\nNormal", side=3, font=2)
# First Level, log-normal
plot(sAP[,"sh_0"], log10(sAP[,"Level1_logNormTime"]), xlab="", ylab="", ylim=fl_Ylim, pch=NA)
mtext("Predicted by\nLog-Normal", side=3, font=2)
# First level, GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level1_time"]), xlab="", ylab="", ylim=fl_Ylim, pch=NA)
mtext("Predicted by\nGEV", side=3, font=2)

#SECOND LEVEL
sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
# second level, normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
lT <- (setThresh*100)-100
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
dev.off()



# ==============================================================
# = Boxplot of shape & return time for Categories of variables =
# ==============================================================
# dev.new(width=3.5, height=6)
pdf("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Box_ReturnShape.pdf", width=3.5, height=6, bg="white")
par(mfrow=c(2,1), mar=c(3,4,0.5,0.5), ps=10, cex=1)
boxplot(log10(Level2_time)~Type, data=sAP, ylab="")

mtext(paste("Return Time (log10(yrs)) \n [record + ", ((setThresh*100)-100), "%]", sep=""), side=2, line=2)
boxplot((sh_0)~Type, data=sAP, ylab="")
mtext(expression(Value~~of~~xi~~from~~GEV), side=2, line=2.5)
dev.off()
# 
# # ===================================================================
# # = # =============================================================
# # = Boxplots & distribution comparison, but with no infinites =
# # ============================================================= =
# # ===================================================================
# sAP2 <- Inf2NA(sAP)
# sAP2 <- sAP2[complete.cases(sAP2),]
# 
# Inf2NA <- function(x) {x[x==-Inf | x==Inf] <- NA; x}
# 
# dev.new(height=5, width=7)
# par(mfrow=c(2,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=10, cex=1, tcl=-0.25, mgp=c(3,0.5,0))
# 
# #FIRST LEVEL
# fl_Ylim <- log10(range(Inf2NA(sAP2[,c("Level1_normTime", "Level1_logNormTime", "Level1_time")]), na.rm=TRUE))
# # First level, normal
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level1_normTime"]), xlab="", ylab="", ylim=fl_Ylim)
# mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
# mtext(bquote(.("[")~mu+2*sigma~.("]")), side=2, line=1.5)
# mtext("Predicted by\nNormal", side=3, font=2)
# # First Level, log-normal
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level1_logNormTime"]), xlab="", ylab="", ylim=fl_Ylim)
# mtext("Predicted by\nLog-Normal", side=3, font=2)
# # First level, GEV
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level1_time"]), xlab="", ylab="", ylim=fl_Ylim)
# mtext("Predicted by\nGEV", side=3, font=2)
# 
# #SECOND LEVEL
# sl_Ylim <- log10(range(Inf2NA(sAP2[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
# # second level, normal
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim)
# # mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
# mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)
# lT <- (setThresh*100)-100
# mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim)
# plot(sAP2[,"sh_0"], log10(sAP2[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim)
# mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
# 
# 
# 
# # ==============================================================
# # = Boxplot of shape & return time for Categories of variables =
# # ==============================================================
# 
# dev.new(width=3.5, height=6)
# par(mfrow=c(2,1), mar=c(3,4,0.5,0.5), ps=10, cex=1)
# boxplot(log10(Level2_time)~Type, data=sAP2, ylab="")
# 
# mtext(paste("Return Time (log10(yrs)) \n [record + ", ((setThresh*100)-100), "%]", sep=""), side=2, line=2)
# boxplot((sh_0)~Type, data=sAP2, ylab="")
# mtext(expression(Value~~of~~xi~~from~~GEV), side=2, line=2.5)
