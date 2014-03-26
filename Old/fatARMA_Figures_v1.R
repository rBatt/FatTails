
rm(list=ls())
graphics.off()

library("plyr")
library("rpart")
library("party")
library("RColorBrewer")
library("gridExtra")
library("fExtremes")
library("ggplot2")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
load("fina_fatARMA_Summary_v5.RData")
source("FatTails_Functions_v8.R")

myWhite <- rgb(t(col2rgb(col="white")),alpha=100, maxColorValue=256)
myWhite2 <- rgb(t(col2rgb(col="white")),alpha=200, maxColorValue=256)
myWhite3 <- rgb(t(col2rgb(col="white")),alpha=240, maxColorValue=256)
bgDef <- c(myWhite, myWhite2, "white")[3]

eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))

figFold <- "fatARMA_Figures_v1"

# Density Plot Color Scheme
# Densities will be estimated and plottted for each column of val (rows are multiple observations of same variable)
pDen <- function(vals=NULL, mu=NULL, sig=NULL){
	if(is.null(vals) & (is.null(mu)|is.null(sig))) stop("Must provide vals or (mu and sig), but not both.")
	if(any(sig<=0)) stop("Need positive sig")
	
	N <- length(mu)
	
	if(is.null(vals)){
		vals0 <- rnorm(n=100*N, mean=mu, sd=sig)
		vals <- matrix(vals0, ncol=N, byrow=TRUE)	
	}
	limX <- range(vals)
	
	dens <- apply(vals, 2, function(x)density(x, from=limX[1], to=limX[2])[c("x","y")])
	
	dX <- lapply(dens, function(z)z$x)
	dY <- lapply(dens, function(z)z$y)
	limY <- range(dY)
	
	cLine <- rainbow(n=N)
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=35, maxColorValue=255)
	
	# Need to be sure that the smallest and largest densities are 0 so that the bottom border of polygons are at 0 line
	xF <- 0.01 * diff(limX) # a "factor" by which to extend the range of X
	xA <- limX + c(-1,1)*xF # "add" this "adjustment" to the start and end of the dX

	# par(mar=c(1.5, 1.5, 0.5, 0.5), ps=10, cex=1, mgp=c(1.5, 0.1, 0), tcl=0.25, las=1)
	plot(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), type="l", col=cLine[1], xlab="", ylab="", xlim=limX, ylim=limY, lwd=2)
	polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[1]],0), col=cFill[1], border=NA)
	for(i in 2:N){
		polygon(c(xA[1],dX[[1]],xA[2]), c(0,dY[[i]],0), col=cFill[i], border=cLine[i], lwd=2)
	}
}





# =====================
# = Period Histrogram =
# =====================
# dev.new(width=3.4, height=3.4)
png(paste(figFold, "periodHistogram.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
perHist_breaks <- hist(final[final[,"Period"]<150,"Period"], plot=FALSE)$breaks
hist(final[final[,"Period"]<150,"Period"], xlab="Period", ylab="Frequency", main="", col="white")
hist(final[final[,"Type"]=="Chem" & final[,"Period"]<150 ,"Period"], add=TRUE, col="gray", breaks=perHist_breaks)
hist(final[final[,"Type"]=="Bio" & final[,"Period"]<150,"Period"], add=TRUE, col="black", breaks=perHist_breaks)
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)
dev.off()



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
png(paste(figFold, "Returns.png", sep="/"), width=3.4, height=6, units="in", res=300, bg=bgDef)
par(mfrow=c(2,1), mar=c(1.5,2.5,0.5,0.5), oma=c(0,0,0,0), ps=9, cex=1, family="Times")
plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:6){
	par(new=TRUE)
	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
}
legend("topleft", c("Uniform", "Normal", "T (df=20)"), bty="n", inset=c(-0.05, -0.03), bg="transparent")
legend("bottomright", c("Exponential", "Log-Normal", "Cauchy", "Weibull"), text.font=c(2,rep(1,3)), bty="n", bg="transparent")
# mtext("Return Time", side=1, line=0.5 , outer=F)

plot(LterReturns[,1], LterReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n", col=LterColors[1])
for(i in 1:(ncol(LterReturns)-2)){
	par(new=TRUE)
	plot(LterReturns[,1], LterReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=1, col=LterColors[i+2], bty="n")	
}
mtext("Return Time", side=1, line=0.25, outer=F)
mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)
legend("bottomright", c("Physical", "Chemical", "Biological"), text.col=c("blue", "gray", "red"), bg="transparent")
dev.off()

#Just the known distributions
# png(paste(figFold, "Returns_1.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
# par(mfrow=c(1,1), mar=c(1.5,2.5,0.5,0.5), oma=c(0,0,0,0), ps=9, cex=1, family="Times")
# plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
# for(i in 1:6){
# 	par(new=TRUE)
# 	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
# }
# legend("topleft", c("Uniform", "Normal", "T (df=20)"), bty="n", inset=c(-0.05, -0.03), bg="transparent")
# legend("bottomright", c("Exponential", "Log-Normal", "Cauchy", "Weibull"), text.font=c(2,rep(1,3)), bty="n", bg="transparent")
# mtext("Return Time", side=1, line=0.25 , outer=F)
# mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)
# 
# dev.off()



# =====================
# = Time series plots =
# =====================
#3 panel graph, each panel is a different category
#pick the fattest variable from each category, and plot its time series
# then plot the time series of that variable from the lake for which the variable is the thinnest
# then draw abline() at the level of the N-year return time, where N is the number of years that the time series was observed
chemLog_sig <- Lake_Params3[,"p_sh_0"]<0.1
chemLog0 <- Lake_Params3[,"Type"] == "Chem" & Lake_Params3[,"N"] > 20
chemLog <- chemLog0 & chemLog_sig
MaxChem_Ind <- which(chemLog)[which.max(Lake_Params3[chemLog, "sh_0"])]
MaxChem <- Lake_Params3[MaxChem_Ind,]
Min_MaxChem_Ind <- which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
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

png(paste(figFold, "fatTimeSeries.png", sep="/"), width=5, height=7, bg=bgDef, units="in", res=300)
# dev.new(height=7, width=5)
# pdf(file="EgTimeSeries_RetLvls_09May2013.pdf", width=5, height=7, pointsize=10)
par(mfrow=c(4,2), mar=c(1.5,1.5,0.5, 0.5), oma=c(0.75, 1, 0, 0), ps=9, cex=1, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.4)
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
mtext(expression(Light~Extinction), side=2, line=1.5)

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
mtext(expression(NO[3] ~ "&" ~NO[2] ~(mu*g~L^-1)), side=2, line=1.5)


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
mtext(expression(Chlorophyll~(mu*g~L^-1)), side=2, line=1.5)

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
mtext(expression(Fish~CPUE~(individuals)), side=2, line=1.5)

CPUE_MinMaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(Min_MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
Low_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]-Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]- Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] - Min_MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]+Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]+ Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] + Min_MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MinMaxTS[,2]))
plot(CPUE_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MinMaxTS, col="blue", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.07))
mtext(expression(Date), side=1, line=-0.25, outer=TRUE)
dev.off()


# ==========================
# = Stationary Time Series =
# ==========================
png(paste(figFold, "fatTimeSeries_stationary.png", sep="/"), width=5, height=7, bg=bgDef, units="in", res=300)
# dev.new(height=7, width=5)
# pdf(file="EgTimeSeries_RetLvls_09May2013.pdf", width=5, height=7, pointsize=10)
par(mfrow=c(4,2), mar=c(1.5,1.5,0.5, 0.5), oma=c(0.75, 1, 0, 0), ps=9, cex=1, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.4)
# ================================
# = Light Extinction time series =
# ================================
LiExt_MaxTS <- subset(Data_X$LiExt, lakeid==as.character(MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MaxTS <- LiExt_MaxTS[complete.cases(LiExt_MaxTS),]
LiExt_MaxTS <- Stationary3(LiExt_MaxTS)
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
mtext(expression(Light~Extinction), side=2, line=1.5)

LiExt_MinMaxTS <- subset(Data_X$LiExt, lakeid==as.character(Min_MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MinMaxTS <- LiExt_MinMaxTS[complete.cases(LiExt_MinMaxTS),]
LiExt_MinMaxTS <- Stationary3(LiExt_MinMaxTS)
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
NO3NO2_MaxTS <- Stationary3(NO3NO2_MaxTS)
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
mtext(expression(NO[3] ~ "&" ~NO[2] ~(mu*g~L^-1)), side=2, line=1.5)


NO3NO2_MinMaxTS <- subset(Data_X$Chem, lakeid==as.character(Min_MaxChem[,"fitBy"]), select=c("sampledate", "no3no2"))
NO3NO2_MinMaxTS <- NO3NO2_MinMaxTS[complete.cases(NO3NO2_MinMaxTS),]
NO3NO2_MinMaxTS <- Stationary3(NO3NO2_MinMaxTS)
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
chlor_MaxTS <- Stationary3(chlor_MaxTS)
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
mtext(expression(Chlorophyll~(mu*g~L^-1)), side=2, line=1.5)

chlor_MinMaxTS <- subset(Data_X$Chl, lakeid==as.character(Min_MaxBio2[,"fitBy"]), select=c("sampledate", "chlor"))
chlor_MinMaxTS <- chlor_MinMaxTS[complete.cases(chlor_MinMaxTS),]
chlor_MinMaxTS <- Stationary3(chlor_MinMaxTS)
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
CPUE_MaxTS <- Stationary3(CPUE_MaxTS)
Low_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]-MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]- MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] - MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]+MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]+ MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] + MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MaxTS[,2]))
plot(CPUE_MaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MaxTS, col="red", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(Fish~CPUE~(individuals)), side=2, line=1.5)

CPUE_MinMaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(Min_MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
CPUE_MinMaxTS <- Stationary3(CPUE_MinMaxTS)
Low_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]-Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]- Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] - Min_MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]+Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]+ Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] + Min_MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MinMaxTS[,2]))
plot(CPUE_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MinMaxTS, col="blue", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.07))
mtext(expression(Date), side=1, line=-0.25, outer=TRUE)
dev.off()

# ==================================================================
# = Return times predicted from 1) normal, 2) log-normal, & 3) GEV =
# ==================================================================
png(paste(figFold, "ScatterReturns.png", sep="/"), width=6, height=2, bg=bgDef, units="in", res=300)

par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=9, cex=1, tcl=-0.25, mgp=c(3,0.5,0), family="Times")

sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
lT <- (setThresh*100)-100

#GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
mtext("Predicted by\nGEV", side=3, font=2)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)

#log normal
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim)
mtext("Predicted by\nLog-Normal", side=3, font=2)

# normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim)
mtext("Predicted by\nNormal", side=3, font=2)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
dev.off()




# ==========================
# = Shape Boxplots by Type =
# ==========================
png(paste(figFold, "Box_ReturnShape.png", sep="/"), width=3.4, height=6, units="in", res=300, bg=bgDef)
par(mfrow=c(2,1), ps=9, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
boxplot((sh_0)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext(expression(Value~~of~~xi~~from~~GEV), side=2, line=2)

boxplot(log10(Level2_time)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext(paste("Return Time (log10(yrs)) \n [record + ", ((setThresh*100)-100), "%]", sep=""), side=2, line=1.5)
dev.off()

# =============================
# = Shape of TS vs. Residuals =
# =============================
png(paste(figFold, "Box_Surprise_TSvsResiduals.png", sep="/"), width=3.4, height=6, units="in", res=300, bg=bgDef)
par(mfrow=c(2,1), ps=9, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
boxplot(log10(Level2_time)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext("Time Series Surprise Rate", side=2, line=2)

boxplot(log10(Level2_time_residual)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext("Residual Surprise Rate", side=2, line=1.5)
dev.off()

# =====================================
# = TS-Residuals: Shape & Return Time =
# =====================================
png(paste(figFold, "Box_Diff_TS&Resid_Shape&Ret.png", sep="/"), width=3.4, height=6, units="in", res=300, bg=bgDef)
par(mfrow=c(2,1), ps=9, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
boxplot(log10(Level2_time)-log10(Level2_time_residual)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext("Time Series Surprise Rate", side=2, line=2)

boxplot((sh_0-residual_sh_0)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext("Residual Surprise Rate", side=2, line=1.5)
dev.off()



# MANY OF THE FIGURES BELOW WERE REVISED IN _V3 TO INCLUDE PHYS

# dev.new(width=3.4, height=3.4)

png(paste(figFold, "Xi_InfE.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
dscat(log10(final[,"InfE"]), final[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(InfE), data=final[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=1, line=1.5)
dev.off()

# dev.new(width=3.4, height=3.4)
# revised in _v3 to include physical

png(paste(figFold, "Xi_Lambda.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
dscat(sqrt(final[,"Lambda"]), final[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(Lambda), data=final[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(.(lNota2))), side=1, line=1.5)
dev.off()

# New _v2


# New _v3

png(paste(figFold, "xiTs_xiResid.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
dscat(final[,"residual_sh_0"], final[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~residual_sh_0, data=final[,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(xi[ARMA~Residuals]), side=1, line=1.5)
dev.off()

# New _v3 ... what happens if you include Physical variables for sigmaE
png(paste(figFold, "Xi_sigE.png", sep="/"), width=3.4, height=3.4, units="in", res=300, bg=bgDef)
dscat(sqrt(final[,"sigE"]), final[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(sigE), data=final[,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(sigma[E])), side=1, line=1.5)
dev.off()



# ===============================
# = Useful Numbers/ Regressions =
# ===============================
summary(lm(log10(Level2_time) ~Type, data=final))
summary(lm(log10(Level2_time) ~log10(Level2_time_residual)+Type, data=final))
summary(lm(sh_0~residual_sh_0+Type, data=final))
log10(colMeans(final[,c("Level2_time_residual","Level2_time")], na.rm=TRUE))
nrow(final)
table(final[,c("Type")])
table(final[,c("Type","fitBy")])
table(final[final[,"sh_0"]>0.5,"Type"])
table(final[final[,"sh_0"]>0.5,"Type"])/table(final[,c("Type")])
table(final[,c("Variable","Type")])
length(unique((final[,"Variable"])))




# ==========================================
# = Comparing ARMA, Norm, Log-Norm, to GEV =
# ==========================================
dTimes <- c("Level2_time","Level2_normTime","Level2_logNormTime","Level2_time_residual")
cutoff <- 0.25
boundLog <- final[,"sh_0"]<= -cutoff
bound <- final[boundLog,]
fatLog <- final[,"sh_0"] => cutoff
fat <- final[fatLog,]
thinLog <- final[,"sh_0"] > -cutoff & final[,"sh_0"] < cutoff #!boundLog & !fatLog
thin <- final[thinLog,]

# Why weren't return time estimated for some time series?
noTime2_log <- is.na(final[,"Level2_time"]) # the time series for which the time wasn't estimated
noLvl_log <- is.na(final[,"Level2"]) # one reason would be that the level 2 threshold was NA
boundedTail_log <- final[,"sh_0"]<0 # for the time series that had bounded upper tails ...
noProb_log <- final[,"Level2"]>(final[,"mu_0"]-final[,"sig_0"]/final[,"sh_0"]) # ... there is the restriction that the threshold may have been beyond that boudning
sum(boundedTail_log & noProb_log, na.rm=TRUE) # 82 cases where the time couldn't be estimated b/c distribution was bounded and threshold was above the bounding
sum(noLvl_log) # 18 cases where the threshold was NA
sum(noTime2_log) # 82 + 18 = the 102 cases where a return time wasn't estimated

# I'm pretty sure that the 18 cases where the threshold was NA resulted from there being fewer than 3 non-NA and non-negative observations. When calculating the mean and the sd, NA was return if fewer than 3 pos values (convNeg is called in lvlMean and lvlSd). Then, when the reutrn levels are merged with the rest of the GEV estimates, complete.cases is applied to the return level data frame.  Therefore, if mean or sd is NA, then the entire row goes NA. As a result, having fewer than 3 pos values makes the mean and sd NA, which makes the level 1 and level 2 NA, which results in the Level2_time also being NA.



dens4 <- function(x){
	df <- -1
	dt <- 11
	d1 <- matrix(unlist(density(log10(x[,"Level2_time"]), na.rm=TRUE, from=df, to=dt)[c("x","y")]), ncol=2)
	d2 <- density(log10(x[,"Level2_normTime"]), na.rm=TRUE, from=df, to=dt)$y
	d3 <- density(log10(x[,"Level2_logNormTime"]), na.rm=TRUE, from=df, to=dt)$y
	d4 <- density(log10(x[,"Level2_time_residual"]), na.rm=TRUE, from=df, to=dt)$y
	
	matrix(c(d1[,1], d1[,2]/max(d1[,2]), d2/max(d2), d3/max(d3), d4/max(d4)), ncol=5, dimnames=list(NULL, c("time","gev","norm","log","resid")))
}
# bCol <- brewer.pal(4,"Set1")
# tCols <- rgb(cbind(t(col2rgb(bCol))), maxColorValue=255, alpha=100)
# tCols2 <- rgb(cbind(t(col2rgb(bCol))), maxColorValue=255, alpha=255)
cLine <- rainbow(n=4, v=0.8, s=1)
cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)
# bCol <- tim.colors(n=4)

pDens <- function(x, go=TRUE, ...){
	myLwd <- 1.25
	tD <- dens4(x)
	tD <- rbind(matrix(c(min(tD[,1]-0.01), 0, 0, 0, 0),ncol=5), tD, matrix(c(max(tD[,1])+0.01, 0, 0, 0, 0), ncol=5))
	# plot(tD[,"time"], tD[,"gev"], pch=NA, ylim=range(tD[,c("gev","norm","log","resid")]), xlim=range(tD[,"time"]), ylab="", ...)
	plot(tD[,"time"], tD[,"gev"], pch=NA, ylim=c(-0.035,1), xlim=range(tD[,"time"]), ylab="", lwd=myLwd, ...)
	if(go){
		polygon(tD[,"time"], tD[,"gev"], col=cFill[1], border=cLine[1], lwd=myLwd)
		polygon(tD[,"time"], tD[,"norm"], col=cFill[2], border=cLine[2], lwd=myLwd)
		polygon(tD[,"time"], tD[,"log"], col=cFill[3], border=cLine[3], lwd=myLwd)
		polygon(tD[,"time"], tD[,"resid"], col=cFill[4], border=cLine[4], lwd=myLwd)

		modes <- tD[apply(tD[,-1], 2, which.max),1]

		if(max(tail(tD[,"norm"]))>0.5){
			nID <- colnames(tD)[-1]!="norm"
		}else{
			nID <- rep(TRUE,4)
		}
		points(modes[nID], rep(-0.0375,sum(nID)), pch=21, bg=cFill[nID], col=cLine[nID])
	}

}

# dev.new(width=3.5, height=6)
png(paste(figFold, "compareTimes.png", sep="/"), width=3.5, height=6, units="in", res=300, bg=bgDef)
par(mfrow=c(3,1), mar=c(1, 1.5, 0.25, 0.25), mgp=c(1, 0.4, 0), tcl=-0.35, cex=1, ps=10, family="Times", oma=c(1,1,0,0), xpd=T)
pDens(bound, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Bounded~(xi<=-0.25)))
text(8, c(0.4, 0.5, 0.6, 0.7), c("GEV", "Normal", "Log-Normal", "GEV ARMA Residuals"), col=cLine)
pDens(thin, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Thin~(-0.25<xi~phantom()<0.25)))
pDens(fat, xaxt="n", xlab="")
axis(1)
text(5, 0.95, bquote(Fat~(xi>=0.25)))
mtext("Waiting Time", side=1, outer=TRUE, line=0)
mtext("Relative Density", side=2, line=0, outer=TRUE)
dev.off()


# dev.new(width=3.5, height=6)
png(paste(figFold, "compareTimes1.png", sep="/"), width=10, height=3.5, units="in", res=300, bg="white")
par(mfrow=c(1,3), mar=c(1.1, 1.5, 0.25, 0.25), mgp=c(1, 0.4, 0), tcl=-0.35, cex=1, ps=12, family="Times", oma=c(1,1,0,0), xpd=T)
pDens(bound, xaxt="s", xlab="")
mtext("Relative Density", side=2, line=0, outer=TRUE)
text(7.85, 0.95, bquote(Bounded~(xi<=-0.25)))
text(7.85, c(0.4, 0.5, 0.6, 0.7), c("GEV", "Normal", "Log-Normal", "GEV ARMA Residuals"), col=cLine)
pDens(thin, xaxt="s", xlab="", go=FALSE)
text(7.85, 0.95, bquote(Thin~(-0.25<xi~phantom()<0.25)))
pDens(fat, xaxt="s", xlab="", go=FALSE)
text(5, 0.95, bquote(Fat~(xi>=0.25)))
mtext("Waiting Time", side=1, outer=TRUE, line=0)
dev.off()

# dev.new(width=3.5, height=6)
png(paste(figFold, "compareTimes2.png", sep="/"), width=10, height=3.5, units="in", res=300, bg="white")
par(mfrow=c(1,3), mar=c(1.1, 1.5, 0.25, 0.25), mgp=c(1, 0.4, 0), tcl=-0.35, cex=1, ps=12, family="Times", oma=c(1,1,0,0), xpd=T)
pDens(bound, xaxt="s", xlab="")
mtext("Relative Density", side=2, line=0, outer=TRUE)
text(7.85, 0.95, bquote(Bounded~(xi<=-0.25)))
text(7.85, c(0.4, 0.5, 0.6, 0.7), c("GEV", "Normal", "Log-Normal", "GEV ARMA Residuals"), col=cLine)
pDens(thin, xaxt="s", xlab="", go=TRUE)
text(7.85, 0.95, bquote(Thin~(-0.25<xi~phantom()<0.25)))
pDens(fat, xaxt="s", xlab="", go=FALSE)
text(5, 0.95, bquote(Fat~(xi>=0.25)))
mtext("Waiting Time", side=1, outer=TRUE, line=0)
dev.off()



# dev.new(width=3.5, height=6)
png(paste(figFold, "compareTimes3.png", sep="/"), width=10, height=3.5, units="in", res=300, bg="white")
par(mfrow=c(1,3), mar=c(1.1, 1.5, 0.25, 0.25), mgp=c(1, 0.4, 0), tcl=-0.35, cex=1, ps=12, family="Times", oma=c(1,1,0,0), xpd=T)
pDens(bound, xaxt="s", xlab="")
mtext("Relative Density", side=2, line=0, outer=TRUE)
text(7.85, 0.95, bquote(Bounded~(xi<=-0.25)))
text(7.85, c(0.4, 0.5, 0.6, 0.7), c("GEV", "Normal", "Log-Normal", "GEV ARMA Residuals"), col=cLine)
pDens(thin, xaxt="s", xlab="", go=TRUE)
text(7.85, 0.95, bquote(Thin~(-0.25<xi~phantom()<0.25)))
pDens(fat, xaxt="s", xlab="", go=TRUE)
text(5, 0.95, bquote(Fat~(xi>=0.25)))
mtext("Waiting Time", side=1, outer=TRUE, line=0)
dev.off()



egQuants2 <- seq(0, 8, length.out=250)
norm2Eg <- dnorm(egQuants2, mean=3, sd=1)
gumbelEg <- dgev(egQuants2, xi=0, mu=3, beta=1)
frechetEg <- dgev(egQuants2, xi=0.51, mu=3, beta=1)
weibullEg <- dgev(egQuants2, xi=-0.51, mu=3, beta=1)

cols3 <- c("normal"="black","gumbel"="#e41a1c", "weibull"="#377eb8", "frechet"="#4daf4a")
fatEg <- cbind(data.frame("quant"=egQuants2, "normal"=norm2Eg), data.frame("gumbel"=gumbelEg), data.frame("weibull"=weibullEg), data.frame("frechet"=frechetEg))
fatGG <- ggplot(fatEg, aes(x=quant))

fatTheme <- theme(title=element_blank(), text=element_blank(), panel.background=element_rect(fill="transparent", colour=NA), plot.background=element_rect(fill="transparent", colour=NA), panel.grid=element_blank(), axis.ticks=element_blank())

# fatGG + geom_line(aes(), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=density, group=type, fill=type), alpha=0.1, size=0, show_guide=FALSE) + fatTheme +annotate("text", x=1, y=0.35, label="xi", parse=TRUE, size=20)

fatColor <- scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3)
fatNote <- annotate("text", x=1, y=0.35, label="xi", parse=TRUE, size=20)

fatNorm <- fatGG + geom_line(aes(y=normal, colour="normal"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=normal, colour="normal", fill="normal"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme + annotate("text", x=1, y=0.35, label="normal", parse=TRUE, size=10)

fatGumbel <- fatGG + geom_line(aes(y=gumbel, colour="gumbel"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=gumbel, ymin=0, ymax=gumbel, colour="gumbel", fill="gumbel"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=5.7, y=0.15, label="xi==0", parse=TRUE, size=10)

fatWeibull <- fatGG + geom_line(aes(y=weibull, colour="weibull"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=weibull, ymin=0, ymax=weibull, colour="weibull", fill="weibull"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=5, y=0.38, label="xi<0", parse=TRUE, size=10)

fatFrechet <- fatGG + geom_line(aes(y=frechet, colour="frechet"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=frechet, ymin=0, ymax=frechet, colour="frechet", fill="frechet"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=7.4, y=0.07, label="xi>0", parse=TRUE, size=10)

	
png("/Users/Battrd/Documents/School&Work/Presentations/Rutgers_IMCS_07March2014/rFigs/eg_fats_grid.png", width=10, height=7, units="in", res=300, bg=myWhite3)
print(grid.arrange(fatNorm, fatGumbel, fatWeibull, fatFrechet, ncol=2))
dev.off()


# Sideways!
# ==========================
# = Shape Boxplots by Type =
# ==========================
png(paste(figFold, "Box_ReturnShape_sidways.png", sep="/"), width=7, height=3.5, units="in", res=300, bg=myWhite2)
par(mfrow=c(1,2), ps=12, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
boxplot((sh_0)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext(expression(Value~~of~~xi~~from~~GEV), side=2, line=2)

boxplot(log10(Level2_time)~Type, data=final, ylab="", outline=FALSE, col="white")
mtext(bquote(Waiting~~Time~~(log[10]*Years)), side=2, line=1.5, outer=FALSE)
dev.off()
	
# =============================================
# = Boxplots: Shape w/ residual shape removed =
# =============================================
png(paste(figFold, "Box_Shape_removeResidual.png", sep="/"), width=3.5, height=3.4, units="in", res=300, bg=myWhite2)
par(mfrow=c(1,1), ps=12, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
sh_rem <- as.numeric(lm(sh_0~residual_sh_0, data=final)$res)
boxplot(sh_rem~final[,"Type"], ylab="", outline=FALSE, col="white")
mtext(bquote(Residuals~(xi[ts]*'~'*xi[arma])), side=2, line=2)
dev.off()


catShape0 <- reshape(aggregate(final[,"sh_0"], by=list("lake"=final[,"fitBy"], "type"=final[,"Type"]), median, na.rm=TRUE), direction="wide", timevar="type", idvar=c("lake"))
catShape <- catShape0[1:11,1:4]
names(catShape) <- c("lake", "bio", "chem", "phys")
summary(lm(bio~chem+phys, data=catShape))

r_catShape0 <- reshape(aggregate(sh_rem, by=list("lake"=final[,"fitBy"], "type"=final[,"Type"]), median, na.rm=TRUE), direction="wide", timevar="type", idvar=c("lake"))
r_catShape <- r_catShape0[1:11,1:4]
names(r_catShape) <- c("lake", "bio", "chem", "phys")
summary(lm(bio~chem+phys, data=r_catShape))


summary(lm(catShape[,"bio"]~r_catShape[,"bio"]))

