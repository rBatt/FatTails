

rm(list=ls())
graphics.off()
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R")
fnNames <- ls()
pkgNames <- c("DEoptim","GenSA")
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")

library("plyr")
library("doParallel")


runARMA <- c(TRUE, FALSE)[1]

finalFrame0 <- data.stat




phRows <- finalFrame0[,"variable"]=="ph"
big0 <- is.na(finalFrame0[,"Data"]) | finalFrame0[,"Data"]>0
finalFrame0[(!phRows)&big0,"Data"] <- log(finalFrame0[(!phRows)&big0,"Data"])

finalFrame <- finalFrame0 # used to be the step where I took the log and made it stationary – already done.

# buildFrame <- list()
# uVar <- unique(finalFrame[,"variable"])
# for(i in 1:length(uVar)){
# 	tFrame <- finalFrame[finalFrame[,"variable"]==uVar[i],]
# 	buildFrame[[i]] <- ddply(.data=tFrame, .variables=c("variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, .progress="none", Method="Evolve", .paropts=list(.export=fnNames, .packages=pkgNames))
# }
# ==================
# = Create fatARMA =
# ==================
if(runARMA){ 
	if(Sys.info()["sysname"]=="Windows"){
		nC <- floor(detectCores()*0.75)
		registerDoParallel(cores=nC)
	}else{
		registerDoParallel()
	}
	
	#run ARMA analysis
	system.time(fatARMA0 <- ddply(.data=finalFrame, .variables=c("variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, Method="Evolve", .paropts=list(.export=fnNames, .packages=pkgNames)))
	
	#Add the "Type" column, using the All_Params data frame (from TurnExtreme_Fat_Data_v8.R) as a reference
	fatARMA <- data.frame("Type"=NA, fatARMA0)
	catKey <- unique(data.fat[,c("Type", "variable")])
	uCats <- as.character(unique(catKey[,"Type"]))
	for(i in 1:length(uCats)){
		tCat <- uCats[i] #current category
		tVars <- catKey[catKey[,"Type"]==tCat,"variable"]#variables in this category
		tI <- is.element(fatARMA[,"variable"], tVars) #rows of fatARMA that correspond to variables in the current category
		fatARMA[tI,"Type"] <- tCat
	}
	
	#save fatARMA
	save(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA.RData", fatARMA)
}





