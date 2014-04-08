#RDB
#FatFrame: This file will also apply the ARMA analysis (from tonyARMA_short_v2.R) to the fatFrame

# _v0 (30-Nov-2013): Take the data from TurnExtreme_Fat_Data_v7.R and turn it into a data frame that has a column for Lake (or region), a column for the date (year), and a column for Data.  This would only be the extremes, not the full data set. I am going to split the data into the categories (cosmic, meteorological, physical, chemical, biological) first, just so that I can run only part of the data if I want to.

#_v1 (01-Dec-2013): Adding the "Type" column to the fatARMA data frame. Adding the "runARMA" option at the start, to safeguard against this long run and potentially overwriting the results.
# _v3 (30-Jan-2013): loading the newest version of ARMA_short, and saving a data frame that has fit the ARMA models to time series that have retained any linear trend

rm(list=ls())
graphics.off()
setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")
load("Data/All_Params_TurnExtreme_Fat_Data.RData")
source("ARMAFunctions.R")
library("plyr")


runARMA <- c(TRUE, FALSE)[2]

finalFrame0 <- data.max




phRows <- finalFrame0[,"variable"]=="ph"
finalFrame0[!phRows,"Data"] <- log(finalFrame0[!phRows,"Data"])

finalFrame_noXform_Stat <- ddply(finalFrame_noXform, c("variable", "location"), Stat)
save(file="finalFrame_noXform_Stat.RData", finalFrame_noXform_Stat)

finalFrame <- ddply(finalFrame0, c("variable", "location"), logStat)
finalFrame_noStat <- ddply(finalFrame0, c("variable", "location"), logStat, doStat=FALSE)
save(finalFrame, finalFrame_noStat, file="/Data/finalFrame.RData")

# ==================
# = Create fatARMA =
# ==================
if(runARMA){
	library("doParallel")
	registerDoParallel()
	#run ARMA analysis
	system.time(fatARMA0_noStat <- ddply(.data=finalFrame_noStat, .variables=c("variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, Method="Evolve"))
	
	#Add the "Type" column, using the All_Params data frame (from TurnExtreme_Fat_Data_v8.R) as a reference
	fatARMA_noStat <- data.frame("Type"=NA, fatARMA0_noStat)
	catKey <- unique(All_Params[,c("Type", "Variable")])
	uCats <- as.character(unique(catKey[,"Type"]))
	for(i in 1:length(uCats)){
		tCat <- uCats[i] #current category
		tVars <- catKey[catKey[,"Type"]==tCat,"Variable"]#variables in this category
		tI <- is.element(fatARMA_noStat[,"variable"], tVars) #rows of fatARMA that correspond to variables in the current category
		fatARMA_noStat[tI,"Type"] <- tCat
	}
	
	#save fatARMA
	save(file="/Data/fatARMA_noStat.RData", fatARMA_noStat)
}
