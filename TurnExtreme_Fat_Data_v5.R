#Version 3 (02-Spet-2013) I am making these changes after the initial effort for Steve and Monica's class --- previous versions weren't documented because I was under a deadline!
#Version 5 (05-Sept-2013) I've gone through the time series and added two "levels" to the Params object.  The first level is the 2*sd+mean of the long-term time series.  The other is 2*mean.  These summary statistics are not calculated to weight each year equally. Nothing fancy.  I also cleaned up a lot of the old code, exported the Functions to a new script, and exported the figures to a new script.

rm(list=ls())
graphics.off()
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Fat_dGEV_v1.0.R")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails")
load("FatTails_rawData/OrganizedFatData_Read_Fat_Data_v1.RData")
library("wmtsa")
library("reshape")
source("FatTails_Functions_v4.R")

ShowPlots <- c(TRUE, FALSE)[2]


# ============
# = Sunspots =
# ============
Sunspot_Ext1 <- MIS(Data_X$SunSpots[,"SpotNum"], Thresh=150)
Sunspot_Ext2 <- PeakCycle(Data=Data_X$SunSpots[,c("SpotNum")], SearchFrac=0.026)
if(ShowPlots){
	dev.new(width=6, height=4)
	par(mar=c(4,4,0.5,0.5))
	plot(Data_X$SunSpots[,c("SpotNum")], type="l")
	points(Sunspot_Ext2, col="blue", pch=20)
}
Sunspot_gev1 <- gev.fit(Sunspot_Ext1[,2])
Sunspot_gev2 <- gev.fit(Sunspot_Ext2[,2])
Sunspot_gev <- Sunspot_gev2


Sunspot_Ext2 <- PeakCycle(Data=Data_X$SunSpots[,c("SpotNum")], SearchFrac=0.026)
SunSpot_Ext <- cbind("One"=1, "year4"=1, "Sunspots"=Sunspot_Ext2[,2])
SunSpot_fit <- calcGEV("SunSpot", datCols="Sunspots", fitBy="One")
SunSpot_Params <- SunSpot_fit[[1]]
SunSpot_gev <- SunSpot_fit[[2]]
SunSpot_Params <- SunSpot_Params[order(SunSpot_Params[,"sh_0"]),]

names(Data_X$SunSpots)[6] <- "Sunspots"
SunSpot_vars <- as.character(unique(SunSpot_Params$Variable))
SunSpot_lvl <- lvl(x=Data_X$SunSpot, variables=SunSpot_vars, fitby=NULL)
SunSpot <- merge(SunSpot_Params, SunSpot_lvl, all=TRUE)
SunSpot[,"fitBy"] <- "None"



# 1 / ((1-0.9) * (79/195)) # 79 is the number of extrema, 195 is the duration of the study (years), and 0.9 is the quantile
# 1 / ( (probability of observing this return level, given that extremes are being observed) * (probability of observing an extreme in the time series)
# 1 / ( (1 - quantile) * (# extrema / duration of time series) )
# So we get a return level of ~182 sunspots per ~25 years.
# Therefore, during a 200-year study, we should have ~8 sunspots that are greater than ~180
# length(which(test[,2]>170)) ---- I'd say that this calculation is approximately correct
# By that I mean that my calculation fo the return level is correct, not necessarily the threshold that I chose before I computed the gev.fit().





# ==================
# = Meteorological =
# ==================
Met_Ext1_01 <- aggregate(Data_X$Met[,c("min_air_temp","ave_air_temp", "max_air_temp", "range_air_temp")], Data_X$Met[,c("year4","location")], FUN=max, na.rm=TRUE)
Met_Ext1_02 <- aggregate(Data_X$Met[,c("precip_mm", "snow_cm")], Data_X$Met[,c("year4","location")], FUN=max, na.rm=TRUE)
Met_Ext1 <- merge(Met_Ext1_01, Met_Ext1_02, all=TRUE)
Met_Ext1 <- Inf2NA(Met_Ext1)
Met_Ext1_Mad <- subset(Met_Ext1, location=="Madison")

if(ShowPlots){
	dev.new(width=6, height=5)
	par(mfrow=c(2,2), mar=c(3.5,4,0.5,0.5))
	AirTrange <- c(min(Met_Ext1_Mad[,"min_air_temp"], na.rm=TRUE), max(Met_Ext1_Mad[,"max_air_temp"], na.rm=TRUE))
	plot(Met_Ext1_Mad[,c("year4","ave_air_temp")], xlab="", ylab="", type="l", ylim=AirTrange)
	lines(Met_Ext1_Mad[,c("year4","min_air_temp")], type="l", col="lightblue")
	lines(Met_Ext1_Mad[,c("year4","max_air_temp")], type="l", col="salmon")
	mtext("Madi. Air Temp Annual Max(daily min mu max)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Mad[,c("year4","range_air_temp")], xlab="", ylab="", type="l")
	mtext("Madi. Air Temp Annual Max(daily range)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Mad[,c("year4","precip_mm")], xlab="", ylab="", type="l")
	mtext("Madi. Annual Max Precip (mm)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Mad[,c("year4","snow_cm")], xlab="", ylab="", type="l")
	mtext("Madi. Annual Max Snow (cm)", side=2, line=2, cex=0.8)

	Met_Ext1_Min <- subset(Met_Ext1, location=="Minocqua")
	dev.new(width=6, height=5)
	par(mfrow=c(2,2), mar=c(3.5,4,0.5,0.5))
	AirTrange <- c(min(Met_Ext1_Min[,"min_air_temp"], na.rm=TRUE), max(Met_Ext1_Min[,"max_air_temp"], na.rm=TRUE))
	plot(Met_Ext1_Min[,c("year4","ave_air_temp")], xlab="", ylab="", type="l", ylim=AirTrange)
	lines(Met_Ext1_Min[,c("year4","min_air_temp")], type="l", col="lightblue")
	lines(Met_Ext1_Min[,c("year4","max_air_temp")], type="l", col="salmon")
	mtext("Minoc. Air Temp Annual Max(daily min mu max)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Min[,c("year4","range_air_temp")], xlab="", ylab="", type="l")
	mtext("Minoc. Air Temp Annual Max(daily range)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Min[,c("year4","precip_mm")], xlab="", ylab="", type="l")
	mtext("Minoc. Annual Max Precip (mm)", side=2, line=2, cex=0.8)
	plot(Met_Ext1_Min[,c("year4","snow_cm")], xlab="", ylab="", type="l")
	mtext("Minoc. Annual Max Snow (cm)", side=2, line=2, cex=0.8)
}
# Met_Ext1[which(Met_Ext1[,"snow_cm"]==0),"snow_cm"] <- NA #to be used if I had the sum of snow, but then a sum doesn't converge to GEV, so max it is!

Met_gev <- list()
MetFrame <- matrix(model.matrix(Met_Ext1[,"min_air_temp"]~Met_Ext1[,"location"])[,2], ncol=1, dimnames=list(NULL, "Location"))
Met_gev$min_air_temp <- gev.fit(Met_Ext1[,"min_air_temp"], ydat=MetFrame, mul="Location", sigl="Location")[c("conv","nllh","mle","se")]
Met_gev$max_air_temp <- gev.fit(Met_Ext1[,"max_air_temp"], ydat=MetFrame, mul="Location", sigl="Location")[c("conv","nllh","mle","se")]
AirRange_NoNA <- Met_Ext1[!is.na(Met_Ext1[,"range_air_temp"]),"range_air_temp"]
MetFrame2 <- matrix(model.matrix(AirRange_NoNA~Met_Ext1[!is.na(Met_Ext1[,"range_air_temp"]),"location"])[,2], ncol=1, dimnames=list(NULL, "Location"))
Met_gev$range_air_temp <- gev.fit(AirRange_NoNA, ydat=MetFrame2, mul="Location", sigl="Location")[c("conv","nllh","mle","se")]
Met_gev$precip_mm <- gev.fit(Met_Ext1[,"precip_mm"], ydat=MetFrame, mul="Location")[c("conv","nllh","mle","se")]
Snow_NoNA <- Met_Ext1[!is.na(Met_Ext1[,"snow_cm"]),"snow_cm"]
MetFrame3 <- matrix(model.matrix(Snow_NoNA~Met_Ext1[!is.na(Met_Ext1[,"snow_cm"]),"location"])[,2], ncol=1, dimnames=list(NULL, "Location"))
Met_gev$snow_cm <- gev.fit(Snow_NoNA, ydat=MetFrame3, mul="Location", shl="Location")[c("conv","nllh","mle","se")]

Met_Ext <- Met_Ext1
AllMet <- c("min_air_temp", "max_air_temp", "range_air_temp", "precip_mm", "snow_cm")
Met_fit <- calcGEV("Met", fitBy="location")
Met_Params <- Met_fit[[1]]
Met_gev <- Met_fit[[2]]
Met_Params <- Met_Params[order(Met_Params[,"sh_0"]),]

Met_vars <- as.character(unique(Met_Params$Variable))
Met_lvl <- lvl(x=Data_X$Met, variables=Met_vars, fitby="location")
Met <- merge(Met_Params, Met_lvl, all=TRUE)


# =======
# = Ice =
# =======
Ice_Ext <- Data_X$Ice
Ice_Ext[,"Region"] <- c("North","South")[(as.integer(is.element(Ice_Ext[,"lakeid"], c("ME","MO","WI")))+1)]
if(ShowPlots){
	dev.new(width=7, height=6)
	par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
	for(i in 1:length(unique(Ice_Ext[,"lakeid"]))){
		ThisIndex <- which(Ice_Ext[,"lakeid"]==unique(Ice_Ext[,"lakeid"])[i])
		plot(Ice_Ext[ThisIndex,"year4"], Ice_Ext[ThisIndex,"DaysOpen"], type="l")
		legend("topright", legend=unique(Ice_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
	mtext("Ice", side=3, line=0.5, outer=TRUE)

}
#I want to KEEP the commented code below; the new way of doing the calculations makes a separate fit for each lake, but the below code did it per region, which I kind of liked.

Ice_fit <- calcGEV(nameVarbl="Ice", datCols="DaysOpen", fitForm="DaysOpen~Region", MUl="Region")
Ice_Params <- Ice_fit[[1]]
Ice_gev <- Ice_fit[[2]]
Ice_Params <- Ice_Params[order(Ice_Params[,"sh_0"]),]

Ice_vars <- as.character(unique(Ice_Params$Variable))
Ice_lvl <- lvl(x=Ice_Ext, variables=Ice_vars, fitby="Region")
Ice <- merge(Ice_Params, Ice_lvl[1,3:4], all=TRUE)

# ==============
# = Lake Level =
# ==============
LakLev_Ext <- aggregate(Data_X$LakLev[,"LakeLevel"], by=Data_X$LakLev[,c("year4", "lakeid")], max, na.rm=TRUE)
names(LakLev_Ext) <- c("year4", "lakeid", "LakeLevel")
if(ShowPlots){
	dev.new(width=7, height=5)
	par(mfrow=c(3,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
	for(i in 1:length(unique(LakLev_Ext[,"lakeid"]))){
		ThisIndex <- which(LakLev_Ext[,"lakeid"]==unique(LakLev_Ext[,"lakeid"])[i])
		plot(LakLev_Ext[ThisIndex,"year4"], LakLev_Ext[ThisIndex,"LakeLevel"], type="l")
		legend("topright", legend=unique(LakLev_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
	mtext("LakLev", side=3, line=0.5, outer=TRUE)
}

LakLev_fit <- calcGEV("LakLev", datCols="LakeLevel", fitBy="lakeid")
LakLev_Params <- LakLev_fit[[1]]
LakLev_gev <- LakLev_fit[[2]]
LakLev_Params <- LakLev_Params[order(LakLev_Params[,"sh_0"]),]

LakLev_vars <- as.character(unique(LakLev_Params$Variable))
LakLev_lvl <- lvl(x=Data_X$LakLev, variables=LakLev_vars, fitby="lakeid")
LakLev <- merge(LakLev_Params, LakLev_lvl, all=TRUE)

# ========
# = Zmix =
# ========
Zmix_Ext <- aggregate(Data_X$Zmix[,"Zmix"], by=Data_X$Zmix[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
names(Zmix_Ext) <- c("year4", "lakeid", "Zmix")
Zmix_Ext[which(Zmix_Ext[,"Zmix"]==-Inf),"Zmix"] <- NA

if(ShowPlots){
	dev.new(width=7, height=6)
	par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
	for(i in 1:length(unique(Zmix_Ext[,"lakeid"]))){
		ThisIndex <- which(Zmix_Ext[,"lakeid"]==unique(Zmix_Ext[,"lakeid"])[i])
		plot(Zmix_Ext[ThisIndex,"year4"], Zmix_Ext[ThisIndex,"Zmix"], type="l")
		legend("topright", legend=unique(Zmix_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
	mtext("Zmix", side=3, line=0.5, outer=TRUE)
}
Zmix_fit <- calcGEV("Zmix", datCols="Zmix", fitBy="lakeid")
Zmix_Params <- Zmix_fit[[1]]
Zmix_gev <- Zmix_fit[[2]]
Zmix_Params <- Zmix_Params[order(Zmix_Params[,"sh_0"]),]

Zmix_vars <- as.character(unique(Zmix_Params$Variable))
Zmix_lvl <- lvl(x=Data_X$Zmix, variables=Zmix_vars, fitby="lakeid")
Zmix <- merge(Zmix_Params, Zmix_lvl, all=TRUE)

# =========
# = LiExt =
# =========
LiExt_Ext <- aggregate(Data_X$LiExt[,"extcoef"], by=Data_X$LiExt[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
names(LiExt_Ext) <- c("year4", "lakeid", "extcoef")
LiExt_Ext[which(LiExt_Ext[,"extcoef"]==-Inf),"extcoef"] <- NA

if(ShowPlots){
	dev.new(width=6, height=5)
	par(mfrow=c(3,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
	for(i in 1:length(unique(LiExt_Ext[,"lakeid"]))){
		ThisIndex <- which(LiExt_Ext[,"lakeid"]==unique(LiExt_Ext[,"lakeid"])[i])
		plot(LiExt_Ext[ThisIndex,"year4"], LiExt_Ext[ThisIndex,"extcoef"], type="l")
		legend("topright", legend=unique(LiExt_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
	mtext("Light Extinction", side=3, line=0.5, outer=TRUE)
}

LiExt_fit <- calcGEV("LiExt", datCols="extcoef", fitBy="lakeid")
LiExt_Params <- LiExt_fit[[1]]
LiExt_gev <- LiExt_fit[[2]]
LiExt_Params <- LiExt_Params[order(LiExt_Params[,"sh_0"]),]

LiExt_vars <- as.character(unique(LiExt_Params$Variable))
LiExt_lvl <- lvl(x=Data_X$LiExt, variables=LiExt_vars, fitby="lakeid")
LiExt <- merge(LiExt_Params, LiExt_lvl, all=TRUE)

# ==========
# = Secchi =
# ==========
Secchi_Ext <- aggregate(Data_X$Secchi[,"Secchi"], by=Data_X$Secchi[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
names(Secchi_Ext) <- c("year4", "lakeid", "Secchi")
Secchi_Ext[which(Secchi_Ext[,"Secchi"]==-Inf),"Secchi"] <- NA

if(ShowPlots){
	dev.new(width=7, height=6)
	par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
	for(i in 1:length(unique(Secchi_Ext[,"lakeid"]))){
		ThisIndex <- which(Secchi_Ext[,"lakeid"]==unique(Secchi_Ext[,"lakeid"])[i])
		plot(Secchi_Ext[ThisIndex,"year4"], Secchi_Ext[ThisIndex,"Secchi"], type="l")
		legend("topright", legend=unique(Secchi_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
	mtext("Secchi", side=3, line=0.5, outer=TRUE)
}
Secchi_fit <- calcGEV("Secchi", datCols="Secchi", fitBy="lakeid")
Secchi_Params <- Secchi_fit[[1]]
Secchi_gev <- Secchi_fit[[2]]
Secchi_Params <- Secchi_Params[order(Secchi_Params[,"sh_0"]),]

Secchi_vars <- as.character(unique(Secchi_Params$Variable))
Secchi_lvl <- lvl(x=Data_X$Secchi, variables=Secchi_vars, fitby="lakeid")
Secchi <- merge(Secchi_Params, Secchi_lvl, all=TRUE)



# ========
# = Phys =
# ========
AllPhys <- c("wtemp","o2","o2sat")
Phys_Ext <- aggregate(Data_X$Phys[,AllPhys], by=Data_X$Phys[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
if(ShowPlots){
	for(j in 1:length(AllPhys)){
		PhysName <- AllPhys[j]
		dev.new(width=7, height=5)
		par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
		for(i in 1:length(unique(Phys_Ext[,"lakeid"]))){
			ThisIndex <- which(Phys_Ext[,"lakeid"]==unique(Phys_Ext[,"lakeid"])[i])
			if(all(is.na(Phys_Ext[ThisIndex,PhysName]))){
				plot(1,1, type="l")
				legend("topright", legend=unique(Phys_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}else{
				plot(Phys_Ext[ThisIndex,"year4"], Phys_Ext[ThisIndex,PhysName], type="l")
				legend("topright", legend=unique(Phys_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}
		}
		mtext(PhysName, side=3, line=0.5, outer=TRUE)
	}
}
Phys_fit <- calcGEV("Phys", datCols=c("wtemp", "o2", "o2sat"), fitBy="lakeid")
Phys_Params <- Phys_fit[[1]]
Phys_gev <- Phys_fit[[2]]
Phys_Params <- Phys_Params[order(Phys_Params[,"sh_0"]),]

Phys_vars <- as.character(unique(Phys_Params$Variable))
Phys_lvl <- lvl(x=Data_X$Phys, variables=Phys_vars, fitby="lakeid")
Phys <- merge(Phys_Params, Phys_lvl, all=TRUE)

# ========
# = Ions =
# ========

AllIons <- c("cl","so4","ca", "mg", "na", "k", "fe", "mn", "cond")
Ions_Ext <- aggregate(Data_X$Ions[,AllIons], by=Data_X$Ions[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
Ions_Ext <- Inf2NA(Ions_Ext)
if(ShowPlots){
	for(j in 1:length(AllIons)){
		IonName <- AllIons[j]
		dev.new(width=7, height=5)
		par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
		for(i in 1:length(unique(Ions_Ext[,"lakeid"]))){
			ThisIndex <- which(Ions_Ext[,"lakeid"]==unique(Ions_Ext[,"lakeid"])[i])
		
			if(all(is.na(Ions_Ext[ThisIndex,IonName]))){
				plot(1,1, type="l")
				legend("topright", legend=unique(Ions_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}else{
				plot(Ions_Ext[ThisIndex,"year4"], Ions_Ext[ThisIndex,IonName], type="l")
				legend("topright", legend=unique(Ions_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}
		}
		mtext(IonName, side=3, line=0.5, outer=TRUE)
	}
}
#cond: AL, BM, SP, TB, and TR have time trend; FI, ME, MO, WI are NA.  CR and CB are stable.  CB maybe fat tail.
# mn: No apparent temporal trends. AL, BM, CR, FI, SP all seem to have very high starting values.
# fe: TB and maybe ME have downward trend.  BM, AL, and definitely SP have high starting value.
# k: ME, MO, SP(?) have shallow trends.
# na: AL, BM, FI, ME, MO, SP, TB, TR all ahve evident upward trends.
# mg: AL, BM, SP, and TR(?) have upward trends.
# ca: AL, BM, CB, SP, TR have upward trends.
# so4: CR, FI, TB have downward trends.  SP has a downward step change. BM might have a downward step change. MO has fairly high first value, then stable at low.
# cl: ME, MO, SP, TR have clear upward trends. FI higher at start, low in middle, then up again. BM has high first value.  AL has a peak.

Ions_fit <- calcGEV("Ions", fitBy="lakeid")
Ions_Params <- Ions_fit[[1]]
Ions_gev <- Ions_fit[[2]]
Ions_Params <- Ions_Params[order(Ions_Params[,"sh_0"]),]

Ions_vars <- as.character(unique(Ions_Params$Variable))
Ions_lvl <- lvl(x=Data_X$Ions, variables=Ions_vars, fitby="lakeid")
Ions <- merge(Ions_Params, Ions_lvl, all=TRUE)

# ========
# = Chem =
# ========
AllChem <- c("ph","alk","dic", "tic", "doc", "toc", "no3no2","nh4", "totnf", "totnuf", "totpf", "totpuf", "drsif", "brsif", "brsiuf", "tpm")
Chem_Ext <- aggregate(Data_X$Chem[,AllChem], by=Data_X$Chem[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
Chem_Ext <- Inf2NA(Chem_Ext)

if(ShowPlots){
	for(j in 1:length(AllChem)){
		ChemName <- AllChem[j]
		dev.new(width=7, height=5)
		par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
		for(i in 1:length(unique(Chem_Ext[,"lakeid"]))){
			ThisIndex <- which(Chem_Ext[,"lakeid"]==unique(Chem_Ext[,"lakeid"])[i])
			if(all(is.na(Chem_Ext[ThisIndex,ChemName]))){
				plot(1,1, type="l")
				legend("topright", legend=unique(Chem_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}else{
				plot(Chem_Ext[ThisIndex,"year4"], Chem_Ext[ThisIndex,ChemName], type="l")
				legend("topright", legend=unique(Chem_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}
		}
		mtext(ChemName, side=3, line=0.5, outer=TRUE)
	}
}
#tpm: fairly spikey, TB looks like it might have a positive linear trend
# brsiuf: very spikey CR, CB, AL look like maybe trend.
# brsif: spikes, for some reason the scaling of the axes seems to be off for SP, TB, and TR (y axis goes to zero, values don't drop below 3k, e.g.)
# drsif: patterns are very different among lakes. CR, AL(?), CB(?) linear trend. SP has huge spike in middle. SP, TB, WI possible quadratic trend.  only 2 data for MO.
# totpuf: huge spike at end for AL, CR and TR. The TR and AL spikes my be data errors. MO and ME almost no data. 
# totpf: huge spikes for AL and CR at end. Maybe data errors.  TR looks OK (no spike like for totpUF). TB has spike in middle.  ME MO little data.
# totnuf: FI and TB might have a slight positive trend.  FI ends in a huge spike. WI spike in middle. ME MO almost no data.
# totnf: CR has huge spike near end, might be data error. FI has spike at end. ME and MO almost no data. Trend in TB, WI?
# nh4: TB trend, maybe. WI spike at end. FI increases for last few points.
# no2: this should be dropped, almost no data.
# no3no2: AL, BM, CB, and CR all show noticable spikes, but they don't look fake. Slight trend in TR, maybe.
# toc: AL, BM, CB, FI, ME, TR, SP, WI all have spikes near end that are suspicious. trend in CR.
# doc: AL, BM, FI, ME, SP, MO, TR, WI have apikes at end that are suspecious. trend in CR.
# tic: Trend in BM, CR, ME, MO, TB, TR. AL, BM, FI, WI go up at end, not too suspicious.
# dic: Trend in BM, ME, maybe TR. a few spikes, nothing terribly crazy.
# alk: Trend in AL, BM, CB(?), CR, FI, SP, TB.  Huge spikes in middle of ME MO, huge spike at start of CB.
# ph: Trend in BM, CR, SP, TB (all increasing). Huge spike in second year for CR and TB (suspicious). More believable spike at start of BM and CB. Maybe downward trend in FI.

Chem_fit <- calcGEV("Chem", fitBy="lakeid")
Chem_Params <- Chem_fit[[1]]
Chem_gev <- Chem_fit[[2]]
Chem_Params <- Chem_Params[order(Chem_Params[,"sh_0"]),]

Chem_vars <- as.character(unique(Chem_Params$Variable))
Chem_lvl <- lvl(x=Data_X$Chem, variables=Chem_vars, fitby="lakeid")
Chem <- merge(Chem_Params, Chem_lvl, all=TRUE)



# =======
# = Chl =
# =======
AllChl <- c("chlor")
Chl_Ext <- aggregate(Data_X$Chl[,AllChl], by=Data_X$Chl[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
names(Chl_Ext) <- c("year4", "lakeid", "chlor")
Chl_Ext <- Inf2NA(Chl_Ext)

Chl_fit <- calcGEV("Chl", fitBy="lakeid")
Chl_Params <- Chl_fit[[1]]
Chl_gev <- Chl_fit[[2]]
Chl_Params <- Chl_Params[order(Chl_Params[,"sh_0"]),]

Chl_vars <- as.character(unique(Chl_Params$Variable))
Chl_lvl <- lvl(x=Data_X$Chl, variables=Chl_vars, fitby="lakeid")

Chl <- merge(Chl_Params, Chl_lvl, na.rm=TRUE)

# ========
# = Zoop =
# ========
Zoop_Ext <- Data_X$Zoop
AllZoop <- c("density", "avg_zoop_mass", "tot_zoop_mass", "avg_length")
Zoop_Ext <- Zoop_Ext[order(Zoop_Ext[,"lakeid"], Zoop_Ext[,"year4"], Zoop_Ext[,"daynum"]),]
Zoop_Ext <- subset(Zoop_Ext, lakeid!="FI" & lakeid!="WI")
Zoop_Ext <- aggregate(Zoop_Ext[,AllZoop], by=Zoop_Ext[,c("year4", "lakeid")], FUN=max, na.rm=TRUE)
Zoop_Ext[,"avg_length"] <- replace(Zoop_Ext[,"avg_length"], list=which(Zoop_Ext[,"avg_length"]==0), NA)
Zoop_Ext[,"avg_zoop_mass"] <- replace(Zoop_Ext[,"avg_zoop_mass"], list=which(is.na(Zoop_Ext[,"avg_length"])), NA) #using the average length here b/c that's how mass is calc'd
Zoop_Ext[,"tot_zoop_mass"] <- replace(Zoop_Ext[,"tot_zoop_mass"], list=which(is.na(Zoop_Ext[,"avg_length"])), NA)

if(ShowPlots){
	dev.new(width=7, height=5)
	par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), ps=9)
	for(i in 1:length(unique(Zoop_Ext[,"lakeid"]))){
		ThisIndex <- which(Zoop_Ext[,"lakeid"]==unique(Zoop_Ext[,"lakeid"])[i])
		plot(Zoop_Ext[ThisIndex,"year4"], Zoop_Ext[ThisIndex,"density"], type="l")
		legend("topright", legend=unique(Zoop_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}

	dev.new(width=7, height=5)
	par(mfrow=c(4,4), mar=c(3,3,0.5,0.5), ps=9)
	Xlim_LakLev <- range(Zoop_Ext[,"year4"], na.rm=TRUE)
	Ylim_LakLev <- range(Zoop_Ext[,"tot_zoop_mass"], na.rm=TRUE)
	for(i in 1:length(unique(Zoop_Ext[,"lakeid"]))){
		ThisIndex <- which(Zoop_Ext[,"lakeid"]==unique(Zoop_Ext[,"lakeid"])[i])
		plot(Zoop_Ext[ThisIndex,"year4"], Zoop_Ext[ThisIndex,"tot_zoop_mass"], type="l")
		legend("topright", legend=unique(Zoop_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}

	dev.new(width=7, height=5)
	par(mfrow=c(4,4), mar=c(3,3,0.5,0.5), ps=9)
	Xlim_LakLev <- range(Zoop_Ext[,"year4"], na.rm=TRUE)
	Ylim_LakLev <- range(Zoop_Ext[,"avg_zoop_mass"], na.rm=TRUE)
	for(i in 1:length(unique(Zoop_Ext[,"lakeid"]))){
		ThisIndex <- which(Zoop_Ext[,"lakeid"]==unique(Zoop_Ext[,"lakeid"])[i])
		plot(Zoop_Ext[ThisIndex,"year4"], Zoop_Ext[ThisIndex,"avg_zoop_mass"], type="l")
		legend("topright", legend=unique(Zoop_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}

	dev.new(width=7, height=5)
	par(mfrow=c(4,4), mar=c(3,3,0.5,0.5), ps=9)
	Xlim_LakLev <- range(Zoop_Ext[,"year4"], na.rm=TRUE)
	Ylim_LakLev <- range(Zoop_Ext[,"avg_length"], na.rm=TRUE)
	for(i in 1:length(unique(Zoop_Ext[,"lakeid"]))){
		ThisIndex <- which(Zoop_Ext[,"lakeid"]==unique(Zoop_Ext[,"lakeid"])[i])
		plot(Zoop_Ext[ThisIndex,"year4"], Zoop_Ext[ThisIndex,"avg_length"], type="l")
		legend("topright", legend=unique(Zoop_Ext[,"lakeid"])[i], cex=0.8, bty="n")
	}
}

Zoop_fit <- calcGEV("Zoop", fitBy="lakeid")
Zoop_Params <- Zoop_fit[[1]]
Zoop_gev <- Zoop_fit[[2]]

Zoop_Params <- Zoop_Params[order(Zoop_Params[,"sh_0"]),]
Zoop_SignifFacts <- SignifShape(x=Zoop_Params, main="Zooplankton Abundance")

Zoop_vars <- as.character(unique(Zoop_Params$Variable))
Zoop_lvl <- lvl(x=Data_X$Zoop, variables=Zoop_vars, fitby="lakeid")

Zoop <- merge(Zoop_Params, Zoop_lvl, all=TRUE)

# ==============
# = Fish Notes =
# ==============
#Meaning of "By Gear": Each year, what is the biggest, heaviest, longest, or most abundant fish you caught with gear type X?
#Meaning of "BySpec" (if left unaggregated): What was the biggest, heaviest... you caught regardless of gear type?  This might be useless, b/c it's most useful features are redundant with ByGear, and the thigns that don't make sense with it can only be fixed w/ BySpec.  An important exception is the potential for introducing species as a covariate.
#Meaning of "BySpec" (if analyzed for a specific species): What was the biggest bass caught?  How many bass were caught? What was the bass cpue?

#Variables that are meaningful for "ByGear": c("total_caught", "cpue1_Sum", "cpue3_WeiEff", "cpue4_LengEff", "max_Wei", "max_Leng") #the length/mass statistics might not work as well b/c it could just be large dependent, e.g., on whether or not you caught a musky that day.  This might be an OK thing to do, but it's probably a more relevant statistic for the 
#Variables that are meaningful for "BySpec"(un agg'd): 
#Variables that are meaningful for "BySpec"(truly by species): 

#For the individual-based metrics, the gearid and spname could matter a lot.  For metrics that are more community-oreiented

#I can calculate the gev by a) aggregating the species together but using the gearid as a covariate; or b) leaving gearid and spname separate, and using them both as covariates.
# "total_caught", "cpue1_Sum", "SumWei", "cpue3_WeiEff", "cpue4_LengEff", "max_Wei", "max_Leng" are probably the best variables to use.  I might have to drop "total_caught" and "SumWei" b/c they're just unnormalized versions of their CPUE counterparts.

# FishVars <- c("cpue1_Sum", "cpue3_WeiEff", "cpue4_LengEff", "max_Wei", "max_Leng")
FishVars <- c("cpue1_Sum", "cpue3_WeiEff")
# ===============
# = Fish_ByGear =
# ===============
AllFish_ByGear <- FishVars
# AllFish_ByGear <- c("total_caught","cpue1_Sum","Nfish", "SumWei", "cpue3_WeiEff", "cpue4_LengEff", "mean_Leng","mean_Wei", "max_Leng", "max_Wei", "min_Leng", "min_Wei")
Fish_ByGear_Ext <- aggregate(Data_X$Fish_ByGear[,AllFish_ByGear], by=Data_X$Fish_ByGear[,c("gearid","year4", "lakeid")], FUN=max, na.rm=TRUE)
Fish_ByGear_Ext <- Inf2NA(Fish_ByGear_Ext)
Fish_ByGear_Ext <- subset(Fish_ByGear_Ext, gearid=="ELFISH") #I should eventually remove this once I can introduce gearid as a covariate

if(ShowPlots){ #note that the peak in 2010 in BM max_Leng is the 1.3 m walleye lol
	for(j in 1:length(AllFish_ByGear)){
		Fish_ByGearName <- AllFish_ByGear[j]
		dev.new(width=7, height=5)
		par(mfrow=c(4,3), mar=c(3,3,0.5,0.5), oma=c(0,0,2,0), ps=9)
		for(i in 1:length(unique(Fish_ByGear_Ext[,"lakeid"]))){
			ThisIndex <- which(Fish_ByGear_Ext[,"lakeid"]==unique(Fish_ByGear_Ext[,"lakeid"])[i])
			if(all(is.na(Fish_ByGear_Ext[ThisIndex,Fish_ByGearName]))){
				plot(1,1, type="l")
				legend("topright", legend=unique(Fish_ByGear_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}else{
				plot(Fish_ByGear_Ext[ThisIndex,"year4"], Fish_ByGear_Ext[ThisIndex,Fish_ByGearName], type="l")
				legend("topright", legend=unique(Fish_ByGear_Ext[,"lakeid"])[i], cex=0.8, bty="n")
			}
		}
		mtext(Fish_ByGearName, side=3, line=0.5, outer=TRUE)
	}
}
 
Fish_ByGear_fit <- calcGEV("Fish_ByGear", fitBy="lakeid")
Fish_ByGear_Params <- cbind("gearid"="EleFish", Fish_ByGear_fit[[1]])
Fish_ByGear_gev <- Fish_ByGear_fit[[2]]
Fish_ByGear_Params <- Fish_ByGear_Params[order(Fish_ByGear_Params[,"sh_0"]),]

Fish_ByGear_vars <- as.character(unique(Fish_ByGear_Params$Variable))
gearInd <- Data_X$Fish_ByGear[,"gearid"]=="ELFISH"
Fish_ByGear_lvl <- lvl(x=Data_X$Fish_ByGear[gearInd,], variables=Fish_ByGear_vars, fitby="lakeid")

AllFish_ByGear <- FishVars
# AllFish_ByGear <- c("total_caught","cpue1_Sum","Nfish", "SumWei", "cpue3_WeiEff", "cpue4_LengEff", "mean_Leng","mean_Wei", "max_Leng", "max_Wei", "min_Leng", "min_Wei")
Fish_ByGear_Ext <- aggregate(Data_X$Fish_ByGear[,AllFish_ByGear], by=Data_X$Fish_ByGear[,c("gearid","year4", "lakeid")], FUN=max, na.rm=TRUE)
Fish_ByGear_Ext <- Inf2NA(Fish_ByGear_Ext)
Fish_ByGear_Ext <- subset(Fish_ByGear_Ext, gearid=="VGN") #I should eventually remove this once I can introduce gearid as a covariate

Fish_ByGear_fit2 <- calcGEV("Fish_ByGear", fitBy="lakeid")
Fish_ByGear_Params2 <- cbind("gearid"="GillNet", Fish_ByGear_fit2[[1]])
Fish_ByGear_gev2 <- Fish_ByGear_fit2[[2]]
Fish_ByGear_Params2 <- Fish_ByGear_Params2[order(Fish_ByGear_Params2[,"sh_0"]),]
Fish_ByGear2_vars <- as.character(unique(Fish_ByGear_Params2$Variable))
gearInd2 <- Data_X$Fish_ByGear[,"gearid"]=="ELFISH"
Fish_ByGear2_lvl <- lvl(x=Data_X$Fish_ByGear[gearInd2,], variables=Fish_ByGear2_vars, fitby="lakeid")

Fish1 <- merge(Fish_ByGear_Params, Fish_ByGear_lvl, all=TRUE)
Fish2 <- merge(Fish_ByGear_Params2, Fish_ByGear2_lvl, all=TRUE)

# Fish_ByGear_SignifFacts <- SignifShape(x=Fish_ByGear_Params, main="Fish", cex=1.5)
# Zoop_SignifFacts <- SignifShape(x=Zoop_Params, main="Zooplankton Abundance", cex=1.5)
# Chl_SignifFacts <- SignifShape(x=Chl_Params, main="Chlorophyll", cex=1.5)
# Chem_SignifFacts <- SignifShape(x=Chem_Params, main="Chemistry", cex=1.5)
# Ions_SignifFacts <- SignifShape(x=Ions_Params, main="Ions", cex=1.5)
# Phys_SignifFacts <- SignifShape(x=Phys_Params, main="Physical Limnology", cex=1.5)
# Secchi_SignifFacts <- SignifShape(x=Secchi_Params, main="Secchi", cex=1.5)
# LiExt_SignifFacts <- SignifShape(x=LiExt_Params, main="Light Extinction", cex=1.5)
# Zmix_SignifFacts <- SignifShape(x=Zmix_Params, main="Light Extinction", cex=1.5)
# LakLev_SignifFacts <- SignifShape(x=LakLev_Params, main="Light Extinction", cex=1.5)
# Ice_SignifFacts <- SignifShape(x=Ice_Params, main="Light Extinction", cex=1.5)
# Met_SignifFacts <- SignifShape(x=Met_Params, main="Meteorological", cex=1.5)
# SunSpots_SignifFacts <- SignifShape(x=SunSpot_Params, main="Sunspots", cex=1.5)

All_Params <- rbind(cbind("Type"="Bio", Fish1[,-3]), cbind("Type"="Bio",  Zoop), cbind("Type"="Bio", Chl), cbind("Type"="Chem", Chem), cbind("Type"="Chem", Ions), cbind("Type"="Phys", Phys), cbind("Type"="Phys", Secchi), cbind("Type"="Phys", LiExt), cbind("Type"="Phys", Zmix),cbind("Type"="Phys", LakLev), cbind("Type"="Phys",  Ice[,-c(6,10,14)]), cbind("Type"="Met", Met), cbind("Type"="Cosmic", SunSpot))

All_Params <- All_Params[order(All_Params[,"sh_0"]),]

Lake_Params <- rbind(cbind("Type"="Bio", Fish_ByGear_Params[,-1]), cbind("Type"="Bio", Zoop_Params), cbind("Type"="Bio", Chl_Params), cbind("Type"="Chem", Chem_Params), cbind("Type"="Chem", Ions_Params), cbind("Type"="Phys", Phys_Params), cbind("Type"="Phys", Secchi_Params), cbind("Type"="Phys", LiExt_Params), cbind("Type"="Phys", Zmix_Params),cbind("Type"="Phys", LakLev_Params))

Lake_Params <- Lake_Params[complete.cases(Lake_Params),]
DropUseless <- !is.element(as.character(Lake_Params[,"Variable"]), c("cpue4_LengEff", "max_Wei", "max_Leng", "avg_zoop_mass", "avg_length", "o2sat"))
Lake_Params3 <- Lake_Params[DropUseless,]
Lake_Params3 <- subset(Lake_Params, !is.element(Variable, c("cpue4_LengEff", "max_Wei", "max_Leng", "avg_zoop_mass", "avg_length", "o2sat")))[,]
Lake_Params3 <- droplevels(Lake_Params3)
levels(Lake_Params3[,"Variable"]) <- c("CPUE", "Weight_PUE", "Zoop_Density", "Zoop_Mass", "Chla", "pH", "Alk", "DIC", "TIC", "DOC", "TOC", "NO3NO2", "NH4", "TotNf", "TotNuf", "TotPf", "TotPuf", "drSif", "BrSif", "BrSiuf", "TPM", "Cl", "SO4", "Ca", "Mg", "Na", "K", "Fe", "Mn", "Cond", "wTemp", "O2", "Secchi", "ExtCoef", "Zmix", "LakeLevel")
Lake_Params3 <- Lake_Params3[order(Lake_Params3[,"sh_0"]),]

row.names(All_Params) <- NULL


save(All_Params, Lake_Params3, file="All_Params_TurnExtreme_Fat_Data_v5.RData")


# ==========================================================================
# = What is the return time for the Level1 (2*se+mean) & Level2 (2*mean) ? =
# ==========================================================================

lvl_return <- function(x, level=1){
	lvl0 <- as.numeric(x[paste("Level",level,sep="")])
	a0 <- as.numeric(x[c("mu_0","sig_0","sh_0")])
	nExts0 <- as.numeric(x["N"])
	TS_Duration0 <- as.numeric(x["Duration"])
	result <- lvlX_ReturnTime(lvl=lvl0, a=a0, nExts=nExts0, TS_Duration=TS_Duration0)
	# names(result) <- NULL
	# row.names(result) <- NULL
	result
}

FirstLevel <- apply(All_Params, 1, lvl_return)
SecondLevel <- apply(All_Params, 1, lvl_return, level=2)




All_Params[,"Level1_time"] <- FirstLevel
All_Params[,"Level2_time"] <- SecondLevel



lvlX_ReturnTime(lvl=lvl0, a=a0, nExts=nExts0, Ts_Duration=Ts_Duration0)

lvlX_ReturnTime(lvl=51, a=c(46.2, 1.55, 0.01297), nExts=15, TS_Duration=15)

dev.new()
boxplot(log10(FirstLevel[All_Params[,"sh_0"]>0 & All_Params[,"N"]>20])~All_Params[All_Params[,"sh_0"]>0 & All_Params[,"N"]>20,c("Type")])
boxplot(log10(FirstLevel[All_Params[,"sh_0"]!=0 & All_Params[,"N"]>20])~All_Params[All_Params[,"sh_0"]!=0 & All_Params[,"N"]>20,c("Type")])

hist(All_Params[All_Params[,"N"]<100,"N"]) #ah, it's just the ice one that's really high, and air temp is the one with N=144

which.max(All_Params[All_Params[,"N"]<400,"N"])
All_Params[All_Params[,"N"]<400,][164,]

mean(Data_X$Chem[Data_X$Chem[,"lakeid"]=="MO","tic"][Data_X$Chem[Data_X$Chem[,"lakeid"]=="MO","tic"]>0], na.rm=TRUE) #monona TIC was so high b/c there were bullshit negative values etc
