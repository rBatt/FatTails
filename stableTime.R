library("fBasics")
library("stabledist")
library("doParallel")

registerDoParallel()

setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/")
load("finalFrame_noXform_Stat.RData")

stableTime2 <- function(x){
	ts <- x[,"Data"]
	
	if(sum(!is.na(ts))<10){
		parsTime <- rep(NA,5)
		names(parsTime) <- c("alpha", "beta", "gamma", "delta", "Level2_stableTime")
		return(parsTime)
	}
	
	f0I <- final0[,"fitBy"]%in%x[,"location"] & final0[,"Variable"]%in%x[,"variable"]
	level2 <- final0[f0I,"Level2"]
	nExts <- final0[f0I, "N"]
	TS_Duration <- final0[f0I,"Duration"]
	
		
	sFit <- tryCatch({stable.fit(ts)}, error=function(cond){return(NA)})
	if(is.na(sFit)){
		parsTime <- rep(NA,5)
		names(parsTime) <- c("alpha", "beta", "gamma", "delta", "Level2_stableTime")
		return(parsTime)
	}
	#sFit <- c(0.75, 0.5, 5, 100)
	
	if(any(is.na(c(level2, nExts, TS_Duration)))){
		sTime <- NA
	}else{
		sTime <- stableTime(level2, a=sFit, nExts=nExts, TS_Duration=TS_Duration)	
	}
	
	
	
	parsTime <- c(sFit[1], sFit[2], sFit[3], sFit[4], sTime)
	names(parsTime) <- c("alpha", "beta", "gamma", "delta", "Level2_stableTime")
	
	return(parsTime)
}
# finalFrame_noXform[finalFrame_noXform[,"variable"]=="cpue1_Sum",] 
# finalFrame_noXform_Stat[finalFrame_noXform_Stat[,"variable"]=="cpue1_Sum",]

# test <- finalFrame_noXform[finalFrame_noXform[,"variable"]=="cpue1_Sum"& finalFrame_noXform[,"location"]=="AL",]
# plot(test[,"Data"], type="l")

# stableCalc <- ddply(finalFrame_noXform_Stat[finalFrame_noXform_Stat[,"variable"]=="cpue1_Sum",], c("variable","location"), stableTime2, .parallel=FALSE)

# stableCalc <- ddply(finalFrame_noXform_Stat, c("variable","location"), stableTime2, .parallel=TRUE, .progress="text")



