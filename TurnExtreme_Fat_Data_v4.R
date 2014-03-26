#Version 3 (02-Spet-2013) I am making these changes after the initial effort for Steve and Monica's class --- previous versions weren't documented because I was under a deadline!


rm(list=ls())
graphics.off()
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Fat_dGEV_v1.0.R")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_rawData")
load("OrganizedFatData_Read_Fat_Data_v1.RData")
library("wmtsa")
library("reshape")

ShowPlots <- c(TRUE, FALSE)[2]



# ====================
# = Extrema Functions =
# ====================
MIS <- function(x, Thresh=NULL, quant=0.9){
	if(is.null(Thresh)){
		Thresh <- quantile(x, quant)
	}
	Above <- as.integer(x>=Thresh)
	udF <- ifelse(Above[1]==1, 1, 0)
	udL <- ifelse(rev(Above)[1]==1, -1, 0)
	UpDown <- c(udF,diff(Above),udL)
	StartStop <- matrix(sort(c(which(UpDown==1), (which(UpDown==-1)-1))), ncol=2, byrow=TRUE)
	Storms <- rep(NA, nrow(StartStop))
	StormIndex <- Storms
	for(i in 1:nrow(StartStop)){
		ThisInd <- StartStop[i,1]:StartStop[i,2]
		StormIndex[i] <- which.max(x[ThisInd]) + ThisInd[1] - 1
		Storms[i] <- x[StormIndex[i]]

	}
	return(matrix(c(StormIndex, Storms), ncol=2, dimnames=list(NULL, c("Index", "Max"))))
}


PeakCycle <- function(Data, SearchFrac=0.028){
	# using package "wmtsa"
	#the SearchFrac parameter just controls how much to look to either side 
	#of wavCWTPeaks()'s estimated maxima for a bigger value
	#see dRange
	Wave <- wavCWT(Data)
	WaveTree <- wavCWTTree(Wave)
	WavePeaks <- wavCWTPeaks(WaveTree, snr.min=5)
	WavePeaks_Times <- attr(WavePeaks, which="peaks")[,"iendtime"]
	
	NewPeakTimes <- c()
	dRange <- round(SearchFrac*length(Data))
	for(i in 1:length(WavePeaks_Times)){
		NewRange <- max(c(WavePeaks_Times[i]-dRange, 1)):min(c(WavePeaks_Times[i]+dRange, length(Data)))
		NewPeakTimes[i] <- which.max(Data[NewRange])+NewRange[1]-1
	}
	
	return(matrix(c(NewPeakTimes, Data[NewPeakTimes]), ncol=2, dimnames=list(NULL, c("PeakIndices", "Peaks"))))
}

Inf2NA <- function(x) {x[which(x==-Inf | x==Inf, arr.ind=TRUE)] <- NA; x}

Stationary <- function(x, coluY, coluX){
	tx <- 0:(nrow(x)-1)
	x2 <- x[,coluX]
	Reg <- lm(x[,coluY]~tx*x2)
	xd <- residuals(Reg) #+ mean(x[,coluY])
	
	return(xd)
}

Stationary2 <- function(x, coluY, coluX){
	tx <- 0:(length(x)-1)
	Reg <- lm(x~tx)
	TrendPval <- summary(Reg)$coef["tx","Pr(>|t|)"]
	if(TrendPval>0.5){
		return(x)
	}else{
		xd <- residuals(Reg) #+ mean(x[,coluY])
		return(xd)
	}
}

# =================
# = GEV Functions =
# =================
ReturnLevel <- function(RetQ, mu, sc, xi){
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi)) 	
}
ReturnQuant <- function(RetTime, nExts, TS_Duration){
	-TS_Duration/(RetTime*nExts) + 1
}
ReturnTime <- function(RetQ, nExts, TS_Duration){
	1/((1-RetQ)*(nExts/TS_Duration))
}

xYr_Lvl <- function(RetTime, nExts, TS_Duration, mu, sc, xi){
	RetQ <- -TS_Duration/(RetTime*nExts) + 1
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi))
}



lvlX_ReturnTime <- function(lvl, mu, sc, xi, nExts, Ts_Duration){ #Added _v4
	tx <- (1 + xi * (xdat - mu)/sc)^(-1/xi) #the t(x) from Wikipedia GEV article
	#pdf <- (1/sc)*tx^(xi+1)*exp(-tx) #this is the pdf, if I ever want to use it
	cdf <- exp(-tx) #this is the CDF, which gives the percentile for a certain return level
	1/((1-cdf)*(nExts/TS_Duration))
}


# Data_Levels <- data.frame("Variable"=NA, "Lake"=NA, "Mean"=NA, "SD"=NA, "uMean"=NA, "uSD"=NA) #For each variable I need to calculate its long-term mean and SD. I will need to break it into per-year first to make sure that I am weight years equally.

lvl1 <- function(x){
	mean(x, na.rm=TRUE)+2*sd(x, na.rm=TRUE)
}

lvl2 <- function(x){
	2*mean(x, na.rm=TRUE)
}
reshape2 <- function(...){
	a <- reshape(...)
	row.names(a) <- NULL
	a <- a[,!names(a)=="id"]
	a
}

# Fish_ByGear2_vars <- as.character(unique(Fish_ByGear_Params2$Variable))
# gearInd <- Data_X$Fish_ByGear[,"gearid"]=="VGN"
# Fish_ByGear2_lvl1 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear2_vars], by=list("fitby"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl1)
# # names(Fish_ByGear_lvl1)[2] <- "Fish_ByGear"
# Fish_ByGear2_lvl2 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear2_vars], by=list("fitby"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl2)
# # names(Fish_ByGear_lvl2)[2] <- "Fish_ByGear"
# Fish_ByGear_lvl1_2 <- reshape2(Fish_ByGear_lvl1, varying=Fish_ByGear_vars, direction="long", v.names="Level1", times=Fish_ByGear_vars, timevar="Variable")
# Fish_ByGear_lvl2_2 <- reshape2(Fish_ByGear_lvl2, varying=Fish_ByGear_vars, direction="long", v.names="Level2", times=Fish_ByGear_vars, timevar="Variable")
# merge(Fish_ByGear_lvl1_2, Fish_ByGear_lvl2_2)

lvl <- function(x, fitby, variables, ...){
	# stopifnot(is.list(fitby))
	if(!is.null(fitby)){
		one0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvl1)
		two0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvl2)
	}else{
		one0 <- data.frame(lvl1(x[,variables])); names(one0) <- variables
		two0 <- data.frame(lvl2(x[,variables])); names(two0) <- variables
	}
	if(length(variables)==1 & !is.null(fitby)){
		names(one0)[2] <- variables
		names(two0)[2] <- variables
	}
	
	one <- reshape2(one0, varying=variables, direction="long", v.names="Level1", times=variables, timevar="Variable")
	two <- reshape2(two0, varying=variables, direction="long", v.names="Level2", times=variables, timevar="Variable")
	
	rlvl <- merge(one, two)
}

# ==========================
# = The GEV Wrapper for Me =
# ==========================

# calcGEV(nameVarbl="Ice", datCols="DaysOpen", fitForm="DaysOpen~Region", MUl="Region")
# calcGEV(nameVarbl="Ice", datCols="DaysOpen", fitBy="Region")

calcGEV <- function(nameVarbl, datCols=NULL, fitBy=NULL, fitForm=NULL, MUl=NULL, SIGl=NULL, SHl=NULL){
	if(is.null(fitBy)&is.null(fitForm) || !is.null(fitBy)&!is.null(fitForm)){stop("One & only one of fitBy and fitForm should be NULL")}
	if(is.null(datCols)){datCols <- get(paste("All",nameVarbl,sep=""))}
	Tempo_gev <- list()
	AllVarbl <- datCols
	for(j in 1:length(AllVarbl)){
		VarblName <- AllVarbl[j]
		Varbl_Ext <- get(paste(nameVarbl,"Ext",sep="_"))
		if(is.null(fitBy)){
			LoopThru <- "Formula" #need to change the formula arugment to just be responses, then paste it together with datCols[i] (i didn't write this thinking about multiple responses)
		}else{
			LoopThru <- unique(Varbl_Ext[,fitBy])
		}
		for(i in seq_along(LoopThru)){
			
			if(!is.null(fitBy)){
				fitByID <- as.character(LoopThru[i])
				ThisFitBy <- which(Varbl_Ext[,fitBy]==fitByID & !is.na(Varbl_Ext[,VarblName]))
				if(length(ThisFitBy)<5){next}
				ThisVarbl <- Varbl_Ext[ThisFitBy,VarblName]
				ThisVarbl <- as.vector(Stationary2(ThisVarbl)) + mean(ThisVarbl)
				Nobs <- length(ThisVarbl)
				dRange <- range(Varbl_Ext[ThisFitBy,"year4"])
				dDuration <- diff(dRange) + 1
				
				Tempo_gev <- gev.fit(ThisVarbl)[c("conv","nllh","mle","se")]
			}
			
			
			if(!is.null(fitForm)){
				ModForm <- with(Varbl_Ext, formula(fitForm))
				AllVars <- all.vars(ModForm)
				RespName <- head(AllVars, 1)
				PredNames <- tail(AllVars, -1)
				TempoFullMat <- model.matrix(ModForm)
				colnames(TempoFullMat) <- AllVars
				ThisVarbl <- as.vector(Varbl_Ext[,RespName])
				ThisVarbl_NoNA <- which(!is.na(ThisVarbl))
				ThisVarbl <- ThisVarbl[ThisVarbl_NoNA]
				TempoYdat <- matrix(TempoFullMat[, PredNames], ncol=length(PredNames), dimnames=list(NULL, PredNames)) #i think model.matrix automatically removes the NA rows
				Nobs <- length(ThisVarbl)
				dRange <- range(Varbl_Ext[,"year4"])
				dDuration <- diff(dRange) + 1
				
				Tempo_gev <- gev.fit(ThisVarbl, ydat=TempoYdat, mul=MUl, sigl=SIGl, shl=SHl)[c("conv","nllh","mle","se")]
			}
			
			tstat <- Tempo_gev$mle/Tempo_gev$se
			tstat2 <- (Tempo_gev$mle["sh_0"]-ifelse(tstat[3]>0, 0.5, -0.5))/(Tempo_gev$se["sh_0"]) #testing to see if the shape parameter is greater than 0.5
		
			pvalue <- c()
			for(k in 1:length(tstat)){
				pvalue[k] <- pt(tstat[k], df=Nobs-1, lower.tail=(tstat[k]<0))
			}
			pvalue2 <-  pt(tstat2, df=Nobs-1, lower.tail=(tstat2<0))#the p-value to see if shape is greater than 0.5
			Tempo_gev$Dates <- c("Nobs"=Nobs, "Dates"=dRange, "Duration"=dDuration)
			Tempo_gev$Pvalues <- c(pvalue, pvalue2)
			names(Tempo_gev$Pvalues) <- c(names(Tempo_gev$mle), "sh_0.5")
		
			TempoNames <- c(names(Tempo_gev$mle), paste("se_",names(Tempo_gev$se), sep=""), paste("p_",names(Tempo_gev$Pvalues), sep=""))
		
			if(all(c(j,i)==1)){
				tTempo_Params <- matrix(c(Tempo_gev$mle, Tempo_gev$se, Tempo_gev$Pvalues), nrow=1, dimnames=list(NULL, TempoNames)) 
				Tempo_Params <- data.frame("fitBy"=LoopThru[i], "Variable"=VarblName, "Duration"=dDuration, "N"=Nobs, tTempo_Params)
				# LoopThru, "Variable"=VarblName, 
			}else{
				tTempo_Params <- matrix(c(Tempo_gev$mle, Tempo_gev$se, Tempo_gev$Pvalues), nrow=1, dimnames=list(NULL, TempoNames))
				Tempo_Params <- rbind(Tempo_Params, data.frame("fitBy"=LoopThru[i], "Variable"=VarblName, "Duration"=dDuration, "N"=Nobs, tTempo_Params))
				# LoopThru, "Variable"=VarblName, 
			}
		}
	}
	return(list(Tempo_Params, Tempo_gev))
}


SignifShape <- function(x, ...){
	cneg <- as.integer(x[,"sh_0"]<0&x[,"p_sh_0"]<0.05)*-1
	c1 <- as.integer(x[,"sh_0"]>0&x[,"p_sh_0"]<0.05)
	c2 <- as.integer(x[,"sh_0"]>0.5&x[,"p_sh_0.5"]<0.05)
	cmat <- matrix(c(cneg, c1, c2), ncol=3)
	rs <- rowSums(cmat)
	plot(x[,"sh_0"], col=c("blue","black","red", "red")[rs+2], pch=ifelse((rs+2)>3, 10, 21), ylab="Shape", ...)
	invisible(rs)
}




# ============
# = Sunspots =
# ============
Sunspot_Ext1 <- MIS(Data_X$SunSpots[,"SpotNum"], Thresh=150)
Sunspot_Ext2 <- PeakCycle(Data=Data_X$SunSpots[,c("SpotNum")], SearchFrac=0.026)
if(ShowPlots){
	dev.new(width=6, height=4)
	par(mar=c(4,4,0.5,0.5))
	plot(Data_X$SunSpots[,c("SpotNum")], type="l")
	# points(Sunspot_Ext1, col="red", pch=20)
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

# Met_Ext1 <- aggregate(Data_X$Met[,c("min_air_temp","ave_air_temp", "max_air_temp", "range_air_temp", "precip_mm", "snow_cm")], Data_X$Met[,c("year4","location")], FUN=max, na.rm=TRUE)
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
# Met_mean <- aggregate(Data_X$Met[,Met_vars], by=list("location"=Data_X$Met[,c("location")]), FUN=mean, na.rm=TRUE)
# Met_sd <- aggregate(Data_X$Met[,Met_vars], by=list("location"=Data_X$Met[,c("location")]), FUN=sd, na.rm=TRUE)
# Met_lvl1 <- aggregate(Data_X$Met[,Met_vars], by=list("fitby"=Data_X$Met[,c("location")]), FUN=lvl1)
# Met_lvl2 <- aggregate(Data_X$Met[,Met_vars], by=list("fitby"=Data_X$Met[,c("location")]), FUN=lvl2)

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
# Ice_Ext <- Ice_Ext[which(complete.cases(Ice_Ext)),]
# IceNewYear <- Ice_Ext[,"year4"]-min(Ice_Ext[,"year4"], na.rm=TRUE) #1853 is the first year
# IceFrame <- matrix(model.matrix(Ice_Ext[,"DaysOpen"]~IceNewYear+Ice_Ext[,"Region"]+Ice_Ext[,"lakeid"])[,-1], ncol=12, dimnames=list(NULL, c("Year", "Region",as.character(unique(Ice_Ext[,"lakeid"]))[-1])))
# Ice_gev <- gev.fit(Ice_Ext[,"DaysOpen"], ydat=IceFrame, mul=c("Year", "Region"))[c("conv","nllh","mle","se")]

# Ice_fit <- calcGEV("Ice", datCols="DaysOpen")
Ice_fit <- calcGEV(nameVarbl="Ice", datCols="DaysOpen", fitForm="DaysOpen~Region", MUl="Region")
# calcGEV(nameVarbl="Ice", datCols="DaysOpen", fitBy="Region")
Ice_Params <- Ice_fit[[1]]
Ice_gev <- Ice_fit[[2]]
Ice_Params <- Ice_Params[order(Ice_Params[,"sh_0"]),]

Ice_vars <- as.character(unique(Ice_Params$Variable))
Ice_lvl <- lvl(x=Ice_Ext, variables=Ice_vars, fitby="Region")
# Ice_lvl1 <- aggregate(Data_X$Ice[,Ice_vars], by=list("Region"=Ice_Ext[,"Region"]), FUN=lvl1)
# Ice_lvl2 <- aggregate(Data_X$Ice[,Ice_vars], by=list("Region"=Ice_Ext[,"Region"]), FUN=lvl2)

Ice <- merge(Ice_Params, Ice_lvl[1,3:4], all=TRUE)

# ==============
# = Lake Level =
# ==============
LakLev_Ext <- aggregate(Data_X$LakLev[,"LakeLevel"], by=Data_X$LakLev[,c("year4", "lakeid")], max, na.rm=TRUE)
names(LakLev_Ext) <- c("year4", "lakeid", "LakeLevel")
# dLakLev <- Stationary(x=LakLev, coluY="LakeLevel", coluX="lakeid")
# LakLev_Ext[,"LakeLevel"] <- dLakLev
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
# LakLevForm <- with(LakLev_Ext, formula(LakeLevel~lakeid))
# LakLev_ModMat <- model.matrix(LakLevForm)[,-1]
# LakLev_gev <- gev.fit(LakLev_Ext[,"LakeLevel"], ydat=LakLev_ModMat, sigl=colnames(LakLev_ModMat))[c("conv","nllh","mle","se")]

LakLev_fit <- calcGEV("LakLev", datCols="LakeLevel", fitBy="lakeid")
LakLev_Params <- LakLev_fit[[1]]
LakLev_gev <- LakLev_fit[[2]]
LakLev_Params <- LakLev_Params[order(LakLev_Params[,"sh_0"]),]

LakLev_vars <- as.character(unique(LakLev_Params$Variable))
LakLev_lvl <- lvl(x=Data_X$LakLev, variables=LakLev_vars, fitby="lakeid")
# LakLev_lvl1 <- aggregate(Data_X$LakLev[,LakLev_vars], by=list("lakeid"=Data_X$LakLev[,"lakeid"]), FUN=lvl1)
# LakLev_lvl2 <- aggregate(Data_X$LakLev[,LakLev_vars], by=list("lakeid"=Data_X$LakLev[,"lakeid"]), FUN=lvl2)

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
# Zmix_lvl1 <- aggregate(data.frame("Zmix"=Data_X$Zmix[,Zmix_vars]), by=list("lakeid"=Data_X$Zmix[,"lakeid"]), FUN=lvl1)
# Zmix_lvl2 <- aggregate(data.frame("Zmix"=Data_X$Zmix[,Zmix_vars]), by=list("lakeid"=Data_X$Zmix[,"lakeid"]), FUN=lvl2)

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
# LiExt_lvl1 <- aggregate(Data_X$LiExt[,LiExt_vars], by=list("lakeid"=Data_X$LiExt[,"lakeid"]), FUN=lvl1)
# names(LiExt_lvl1)[2] <- "LiExt"
# LiExt_lvl2 <- aggregate(Data_X$LiExt[,LiExt_vars], by=list("lakeid"=Data_X$LiExt[,"lakeid"]), FUN=lvl2)
# names(LiExt_lvl2)[2] <- "LiExt"

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
# Secchi_lvl1 <- aggregate(Data_X$Secchi[,Secchi_vars], by=list("lakeid"=Data_X$Secchi[,"lakeid"]), FUN=lvl1)
# names(Secchi_lvl1)[2] <- "Secchi"
# Secchi_lvl2 <- aggregate(Data_X$Secchi[,Secchi_vars], by=list("lakeid"=Data_X$Secchi[,"lakeid"]), FUN=lvl2)
# names(Secchi_lvl2)[2] <- "Secchi"

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
# Phys_lvl1 <- aggregate(Data_X$Phys[,Phys_vars], by=list("lakeid"=Data_X$Phys[,"lakeid"]), FUN=lvl1)
# # names(Phys_lvl1)[2] <- "Phys"
# Phys_lvl2 <- aggregate(Data_X$Phys[,Phys_vars], by=list("lakeid"=Data_X$Phys[,"lakeid"]), FUN=lvl2)
# # names(Phys_lvl2)[2] <- "Phys"

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
# Ions_lvl1 <- aggregate(Data_X$Ions[,Ions_vars], by=list("lakeid"=Data_X$Ions[,"lakeid"]), FUN=lvl1)
# # names(Ions_lvl1)[2] <- "Ions"
# Ions_lvl2 <- aggregate(Data_X$Ions[,Ions_vars], by=list("lakeid"=Data_X$Ions[,"lakeid"]), FUN=lvl2)
# # names(Ions_lvl2)[2] <- "Ions"

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
# Chem_lvl1 <- aggregate(Data_X$Chem[,Chem_vars], by=list("lakeid"=Data_X$Chem[,"lakeid"]), FUN=lvl1)
# # names(Chem_lvl1)[2] <- "Chem"
# Chem_lvl2 <- aggregate(Data_X$Chem[,Chem_vars], by=list("lakeid"=Data_X$Chem[,"lakeid"]), FUN=lvl2)
# # names(Chem_lvl2)[2] <- "Chem"

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
# Chl_lvl1 <- aggregate(Data_X$Chl[,Chl_vars], by=list("lakeid"=Data_X$Chl[,"lakeid"]), FUN=lvl1)
# names(Chl_lvl1)[2] <- "Chl"
# Chl_lvl2 <- aggregate(Data_X$Chl[,Chl_vars], by=list("lakeid"=Data_X$Chl[,"lakeid"]), FUN=lvl2)
# names(Chl_lvl2)[2] <- "Chl"

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
# Zoop_lvl1 <- aggregate(Data_X$Zoop[,Zoop_vars], by=list("lakeid"=Data_X$Zoop[,"lakeid"]), FUN=lvl1)
# # names(Zoop_lvl1)[2] <- "Zoop"
# Zoop_lvl2 <- aggregate(Data_X$Zoop[,Zoop_vars], by=list("lakeid"=Data_X$Zoop[,"lakeid"]), FUN=lvl2)
# # names(Zoop_lvl2)[2] <- "Zoop"

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
# Fish_ByGear_lvl1 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear_vars], by=list("lakeid"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl1)
# # names(Fish_ByGear_lvl1)[2] <- "Fish_ByGear"
# Fish_ByGear_lvl2 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear_vars], by=list("lakeid"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl2)
# # names(Fish_ByGear_lvl2)[2] <- "Fish_ByGear"

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
# Fish_ByGear2_lvl1 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear2_vars], by=list("fitby"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl1)
# # names(Fish_ByGear_lvl1)[2] <- "Fish_ByGear"
# Fish_ByGear2_lvl2 <- aggregate(Data_X$Fish_ByGear[gearInd,Fish_ByGear2_vars], by=list("fitby"=Data_X$Fish_ByGear[gearInd,"lakeid"]), FUN=lvl2)
# # names(Fish_ByGear_lvl2)[2] <- "Fish_ByGear"

Fish1 <- merge(Fish_ByGear_Params, Fish_ByGear_lvl, all=TRUE)
Fish2 <- merge(Fish_ByGear_Params2, Fish_ByGear2_lvl, all=TRUE)

# ===============
# = Fish_BySpec =
# ===============


Fish_ByGear_SignifFacts <- SignifShape(x=Fish_ByGear_Params, main="Fish", cex=1.5)
Zoop_SignifFacts <- SignifShape(x=Zoop_Params, main="Zooplankton Abundance", cex=1.5)
Chl_SignifFacts <- SignifShape(x=Chl_Params, main="Chlorophyll", cex=1.5)
Chem_SignifFacts <- SignifShape(x=Chem_Params, main="Chemistry", cex=1.5)
Ions_SignifFacts <- SignifShape(x=Ions_Params, main="Ions", cex=1.5)
Phys_SignifFacts <- SignifShape(x=Phys_Params, main="Physical Limnology", cex=1.5)
Secchi_SignifFacts <- SignifShape(x=Secchi_Params, main="Secchi", cex=1.5)
LiExt_SignifFacts <- SignifShape(x=LiExt_Params, main="Light Extinction", cex=1.5)
Zmix_SignifFacts <- SignifShape(x=Zmix_Params, main="Light Extinction", cex=1.5)
LakLev_SignifFacts <- SignifShape(x=LakLev_Params, main="Light Extinction", cex=1.5)
Ice_SignifFacts <- SignifShape(x=Ice_Params, main="Light Extinction", cex=1.5)
Met_SignifFacts <- SignifShape(x=Met_Params, main="Meteorological", cex=1.5)
SunSpots_SignifFacts <- SignifShape(x=SunSpot_Params, main="Sunspots", cex=1.5)

# All_Params <- rbind(cbind("Type"="Bio", Fish_ByGear_Params[,-1]), cbind("Type"="Bio", Zoop_Params), cbind("Type"="Bio", Chl_Params), cbind("Type"="Chem", Chem_Params), cbind("Type"="Chem", Ions_Params), cbind("Type"="Phys", Phys_Params), cbind("Type"="Phys", Secchi_Params), cbind("Type"="Phys", LiExt_Params), cbind("Type"="Phys", Zmix_Params),cbind("Type"="Phys", LakLev_Params), cbind("Type"="Phys", Ice_Params[,-c(4,8,12)]), cbind("Type"="Met.", Met_Params), cbind("Type"="Cosmic", "fitBy"="Sun", SunSpot_Params[,-1]))
# 
# Ice_Params0 <- Ice_Params
# Ice_Params <- cbind("Type"="Phys", "gearid"=NA,Ice_Params0[,-c(6,10,14)])

# All_Params <- rbind(cbind("Type"="Bio", Fish_ByGear_Params), cbind("Type"="Bio", Fish_ByGear_Params2), cbind("Type"="Bio", "gearid"=NA, Zoop_Params), cbind("Type"="Bio", "gearid"=NA,Chl_Params), cbind("Type"="Chem", "gearid"=NA,Chem_Params), cbind("Type"="Chem", "gearid"=NA,Ions_Params), cbind("Type"="Phys", "gearid"=NA,Phys_Params), cbind("Type"="Phys", "gearid"=NA,Secchi_Params), cbind("Type"="Phys", "gearid"=NA,LiExt_Params), cbind("Type"="Phys", "gearid"=NA,Zmix_Params),cbind("Type"="Phys", "gearid"=NA,LakLev_Params), cbind(Ice_Params), cbind("Type"="Met.", "gearid"=NA,Met_Params), cbind("Type"="Cosmic", "fitBy"="Sun", "gearid"=NA,SunSpot_Params[,-1]))

# All_Params <- rbind(cbind("Type"="Bio", Fish_ByGear_Params), cbind("Type"="Bio", "gearid"=NA, Zoop_Params), cbind("Type"="Bio", "gearid"=NA,Chl_Params), cbind("Type"="Chem", "gearid"=NA,Chem_Params), cbind("Type"="Chem", "gearid"=NA,Ions_Params), cbind("Type"="Phys", "gearid"=NA,Phys_Params), cbind("Type"="Phys", "gearid"=NA,Secchi_Params), cbind("Type"="Phys", "gearid"=NA,LiExt_Params), cbind("Type"="Phys", "gearid"=NA,Zmix_Params),cbind("Type"="Phys", "gearid"=NA,LakLev_Params), cbind(Ice_Params), cbind("Type"="Met.", "gearid"=NA,Met_Params), cbind("Type"="Cosmic", "fitBy"="Sun", "gearid"=NA,SunSpot_Params[,-1]))

# All_Params <- rbind(cbind("Type"="Bio", Fish2), cbind("Type"="Bio", "gearid"=NA, Zoop), cbind("Type"="Bio", "gearid"=NA,Chl), cbind("Type"="Chem", "gearid"=NA,Chem), cbind("Type"="Chem", "gearid"=NA,Ions), cbind("Type"="Phys", "gearid"=NA,Phys), cbind("Type"="Phys", "gearid"=NA,Secchi), cbind("Type"="Phys", "gearid"=NA,LiExt), cbind("Type"="Phys", "gearid"=NA,Zmix),cbind("Type"="Phys", "gearid"=NA,LakLev), cbind("Type"="Phys", "gearid"=NA, Ice[,-c(6,10,14)]), cbind("Type"="Met", "gearid"=NA,Met), cbind("Type"="Cosmic", "gearid"=NA,SunSpot))


All_Params <- rbind(cbind("Type"="Bio", Fish2[,-3]), cbind("Type"="Bio",  Zoop), cbind("Type"="Bio", Chl), cbind("Type"="Chem", Chem), cbind("Type"="Chem", Ions), cbind("Type"="Phys", Phys), cbind("Type"="Phys", Secchi), cbind("Type"="Phys", LiExt), cbind("Type"="Phys", Zmix),cbind("Type"="Phys", LakLev), cbind("Type"="Phys",  Ice[,-c(6,10,14)]), cbind("Type"="Met", Met), cbind("Type"="Cosmic", SunSpot))

rbind(  cbind("Type"="Met", Met), cbind("Type"="Cosmic", "fitBy"="Sun", SunSpot))


All_Params <- All_Params[order(All_Params[,"sh_0"]),]
All_SignifFacts <- SignifShape(All_Params)



dev.new(width=11, height=7)
par(mfrow=c(2,3), mar=c(4,4,1,1), ps=11)
plot(dnorm(seq(-6,6,by=0.05)), type="l", ylab="Probability", xlab="X")
lines(dt(seq(-6,6,by=0.05), df=1), type="l", lwd=4)
boxplot(sh_0~Type, data=All_Params)

plot(1,1, pch="", xlab="", ylab="", xaxt="n", yaxt="n", bty="n")
text(x=1, y=1, labels="49 Variables \n from 'Cosmic' to 'Biological' \n 30+ years data & 13 lakes \n 432 values for 'shape' ", cex=1.5)

SignifShape(x=subset(All_Params, Type=="Phys")[order(subset(All_Params, Type=="Phys")[,"sh_0"]),])
text(x=50 ,y=-1, "physical", cex=1.5)
SignifShape(x=subset(All_Params, Type=="Chem")[order(subset(All_Params, Type=="Chem")[,"sh_0"]),])
text(x=200 ,y=-0.5, "chemical", cex=1.5)
SignifShape(x=subset(All_Params, Type=="Bio")[order(subset(All_Params, Type=="Bio")[,"sh_0"]),])
text(x=80 ,y=-0.5, "biological", cex=1.5)



# ==========
# = BoTraf =
# ==========

#eh, heck with it




save(All_Params, file="All_Params_TurnExtreme_Fat_Data_v4.RData")
# 
# 
# summary(lm(sh_0~Variable, data=Zoop_Params)) #mean zoop mass and mean zoop length are more thin-tailed than zoop density
# summary(lm(sh_0~Variable, data=Fish_ByGear_Params)) #maximum fish length is more thin-tailed than the others (cpue1_Sum as reference level)
# summary(lm(sh_0~Variable, data=Chem_Params))
# summary(lm(sh_0~Variable, data=Ions_Params))
# summary(lm(sh_0~Variable, data=Met_Params)) #the shapes of the meteorological params aren't different from each other
# 
# summary(lm(sh_0~fitBy+Variable, data=All_Params))
# 
# summary(lm(sh_0~fitBy+Variable, weights=(1/(se_sh_0)^2), data=All_Params))
# 
# 
# Lake_Params <- rbind(cbind("Type"="Bio", Fish_ByGear_Params[,-1]), cbind("Type"="Bio", Zoop_Params), cbind("Type"="Bio", Chl_Params), cbind("Type"="Chem", Chem_Params), cbind("Type"="Chem", Ions_Params), cbind("Type"="Phys", Phys_Params), cbind("Type"="Phys", Secchi_Params), cbind("Type"="Phys", LiExt_Params), cbind("Type"="Phys", Zmix_Params),cbind("Type"="Phys", LakLev_Params))
# Lake_Params <- Lake_Params[complete.cases(Lake_Params),]
# 
# # =====================================
# # = Try some simple linear regression =
# # =====================================
# summary(lm(sh_0~fitBy-1, data=Lake_Params))
# summary(lm(sh_0~Variable, data=Lake_Params))
# summary(lm(sh_0~Type, data=Lake_Params))
# summary(lm(sh_0~fitBy+Type -1 , data=Lake_Params))
# summary(lm(sh_0~fitBy+Variable, weights=(1/(se_sh_0)^2), data=Lake_Params))
# 
# 
# # =======================
# # = Looks at some plots =
# # =======================
# plot(Lake_Params[,"fitBy"], residuals(lm(sh_0~Variable, data=Lake_Params)))
# plot(Lake_Params[,"Variable"], residuals(lm(sh_0~fitBy, data=Lake_Params)))
# boxplot(sh_0~Variable, data=Lake_Params)
# 
# 
# 
# # =============
# # = Try a PCA =
# # =============
# Lake_Shapes <- Lake_Params[,c(1:3,6)]
# Wide_LS <- reshape(Lake_Shapes, v.names=c("sh_0"), direction="wide", timevar="Variable", idvar="fitBy", drop="Type")
# Wide_LS <- subset(Wide_LS, !is.element(fitBy, c("KE", "WA")))
# PCAdat <- as.matrix(x=Wide_LS[,-1])
# rownames(PCAdat) <- Wide_LS[,1]
# complete.cases(PCAdat)
# PCA <- princomp(~PCAdat, na.action=na.omit)
# biplot(PCA)
# 
# 
# # ======================
# # = Back to regression =
# # ======================
# summary(lm(sh_0~Type+fitBy-1, data=Lake_Params))
# Lake_Params2 <- Lake_Params
# Lake_Params2[,"sh_0"] <- scale(Lake_Params2[,"sh_0"], scale=FALSE)
# summary(lm(sh_0~fitBy*Type -1, data=Lake_Params2)) #the chemistry and the physics were more fat-tailed in crystal bog than in the other lakes (but CB Bio less FT'd); WI, CR, MO, TB AL, more fat-tailed.
# #My overall interpretation: Some lakes (AL, CR, MO, TB, WI) tended to have variables that were more heavy tailed (higher than average shape parameter).  No lake appeared to be particularly more thin-tailed than the others (i.e., significantly lower shape parameter than average). The biological variables were more fat-tailed than the chemistry or physics variables.  However, the physics and chemistry variables were more heavy-tailed in CB than in the other lakes, and the estimates on these parameters puts CB phys and chemistry close to CB biology.  I.e., In crystal bog the three categories have similar tailedness.  So 
# 
# # ===================================================
# # = I need to drop different version of same variable =
# # ===================================================
# summary(lm(sh_0~Variable-1, data=Lake_Params2))
# #drop fish max length
# 
# 
# # DropUseless <- which(!is.element(as.character(Lake_Params[,"Variable"]), c("cpue4_LengEff", "max_Wei", "max_Leng", "avg_zoop_mass", "avg_length", "o2sat")))
# # Lake_Params3 <- Lake_Params[DropUseless,]
# Lake_Params3 <- subset(Lake_Params, !is.element(Variable, c("cpue4_LengEff", "max_Wei", "max_Leng", "avg_zoop_mass", "avg_length", "o2sat")))[,]
# Lake_Params3 <- droplevels(Lake_Params3)
# 
# summary(lm(sh_0~Variable-1, data=Lake_Params3))
# 
# levels(Lake_Params3[,"Variable"]) <- c("CPUE", "Weight_PUE", "Zoop_Density", "Zoop_Mass", "Chla", "pH", "Alk", "DIC", "TIC", "DOC", "TOC", "NO3NO2", "NH4", "TotNf", "TotNuf", "TotPf", "TotPuf", "drSif", "BrSif", "BrSiuf", "TPM", "Cl", "SO4", "Ca", "Mg", "Na", "K", "Fe", "Mn", "Cond", "wTemp", "O2", "Secchi", "ExtCoef", "Zmix", "LakeLevel")
# 
# summary(lm(sh_0~fitBy*Type-1, data=Lake_Params3))
# test <- Lake_Params3
# test[,"sh_0"] <- scale(Lake_Params3[,"sh_0"], scale=FALSE)
# plot(lm(sh_0~fitBy*Type-1, data=test))
# summary(lm(sh_0~fitBy+Type-1, data=test))
# 
# Test2 <- subset(Lake_Params, !is.element(Variable, c("cpue4_LengEff", "max_Wei", "max_Leng", "avg_zoop_mass", "avg_length", "o2sat")) & !is.element(fitBy, c("KE", "WA")))[,]
# Test2 <- droplevels(Test2)
# summary(lm(sh_0~fitBy+Type, data=Test2))
# 
# save(All_Params, Lake_Params3, file="All_Params_TurnExtreme_Fat_Data_v2.RData")
# 
# Lake_Params3 <- Lake_Params3[order(Lake_Params3[,"sh_0"]),]
# 
# summary(lm(sh_0~Type+fitBy-1, data=Lake_Params3))



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

# dev.new(width=3.5, height=3.5)



make.rgb(t(col2rgb("gray"))[1])



# ===========================================
# = Create a plot combining common and LTER =
# ===========================================
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatFigures/")
# dev.new(width=3.5, height=6)
# setEPS()
# postscript(file="RetLevel_vs_RetTime_09May2013.eps", width=3.5, height=6, pointsize=10)
# tiff(filename="RetLevel_vs_RetTime_09May2013.tiff", width=3.5, height=6, units="in", compression="lzw", res=600, pointsize=10)
# png(filename="RetLevel_vs_RetTime_09May2013.png", width=3.5, height=6, units="in", res=200, pointsize=10)
# cairo_ps(filename="RetLevel_vs_RetTime_09May2013.eps", width=3.5, height=6, pointsize=10) #produces a bitmap with low resolution
# pdf(file="RetLevel_vs_RetTime_09May2013.pdf", width=3.5, height=6, pointsize=10)
dev.new(width=3.5, height=6)
par(mfrow=c(2,1), mar=c(1,2.5,0.5,0.5), oma=c(1.5,0,0,0), ps=10)
plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:6){
	par(new=TRUE)
	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
	
}


plot(LterReturns[,1], LterReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n", col=LterColors[1])
for(i in 1:(ncol(LterReturns)-2)){
	par(new=TRUE)
	plot(LterReturns[,1], LterReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=1, col=LterColors[i+2], bty="n")
	
}
mtext("Return Time", side=1, line=-0.5, outer=TRUE)
mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)
legend("bottomright", c("Physical", "Chemical", "Biological"), text.col=c("blue", "gray", "red"))
# dev.off()


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



dev.new(height=7, width=5)
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
# par(new=TRUE)
# plot(0,0, pch="", xaxt="n", yaxt="n", ylab="", xlab="", bty="n")
# text(x=-0.9, y=0.9, bquote(xi ==  .(round(MaxPhys[,"sh_0"],2))))
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

# dev.off()

# =====================================
# = A table summarizing the variables =
# =====================================


SummaryParams <- ddply(Lake_Params3, "Variable", colwise(range, .cols=c("N","Duration", "mu_0", "sig_0", "sh_0")))

