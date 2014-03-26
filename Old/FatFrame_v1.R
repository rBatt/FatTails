#RDB
#FatFrame: This file will also apply the ARMA analysis (from tonyARMA_short_v2.R) to the fatFrame

# _v0 (30-Nov-2013): Take the data from TurnExtreme_Fat_Data_v7.R and turn it into a data frame that has a column for Lake (or region), a column for the date (year), and a column for Data.  This would only be the extremes, not the full data set. I am going to split the data into the categories (cosmic, meteorological, physical, chemical, biological) first, just so that I can run only part of the data if I want to.

#_v1 (01-Dec-2013): Adding the "Type" column to the fatARMA data frame. Adding the "runARMA" option at the start, to safeguard against this long run and potentially overwriting the results.


rm(list=ls())
graphics.off()
setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
source("TonySuggestions/tonyARMA_short_v2.R")
library("plyr")

runARMA <- c(TRUE, FALSE)[2]

# ==========
# = Cosmic =
# ==========
#SunSpot_Ext[,3]
# spotFrame <- 


# ==================
# = Meteorological =
# ==================
metFrame <- reshape2(Met_Ext[,!is.element(names(Met_Ext), "ave_air_temp")], varying=AllMet, times=AllMet, v.names="Data", direction="long")
names(metFrame) <- c("year4", "location", "variable", "Data")


# ============
# = Physical =
# ============
iceFrame0 <- Ice_Ext[,-4]
names(iceFrame0) <- c("year4", "location", "Data")
iceFrame <- cbind(iceFrame0, "variable"="DaysOpen")

levelFrame0 <- LakLev_Ext
names(levelFrame0) <- c("year4", "location", "Data")
levelFrame <- cbind(levelFrame0, "variable"="LakeLevel")

zmixFrame0 <- Zmix_Ext
names(zmixFrame0) <- c("year4", "location", "Data")
zmixFrame <- cbind(zmixFrame0, "variable"="Zmix")

lightFrame0 <- LiExt_Ext
names(lightFrame0) <- c("year4", "location", "Data")
lightFrame <- cbind(lightFrame0, "variable"="extcoef")

secchiFrame0 <- Secchi_Ext
names(secchiFrame0) <- c("year4", "location", "Data")
secchiFrame <- cbind(secchiFrame0, "variable"="Secchi")

miscPhysFrame <- reshape2(Phys_Ext, varying=AllPhys, times=AllPhys, v.name="Data", direction="long")
names(miscPhysFrame) <- c("year4", "location", "variable", "Data")

physFrame <- rbind(iceFrame, levelFrame, zmixFrame, lightFrame, secchiFrame, miscPhysFrame)


# =============
# = Chemistry =
# =============
ionFrame <- reshape2(Ions_Ext, varying=AllIons, times=AllIons, v.name="Data", direction="long")
names(ionFrame) <- c("year4", "location", "variable", "Data")

miscChemFrame <- reshape2(Chem_Ext, varying=AllChem, times=AllChem, v.name="Data", direction="long")
names(miscChemFrame) <- c("year4", "location", "variable", "Data")

chemFrame <- rbind(ionFrame, miscChemFrame)


# ===========
# = Biology =
# ===========
chlFrame0 <- Chl_Ext
names(chlFrame0) <- c("year4", "location", "Data")
chlFrame <- cbind(chlFrame0, "variable"="chlor")

zoopFrame <- reshape2(Zoop_Ext, varying=AllZoop, times=AllZoop, v.name="Data", direction="long")
names(zoopFrame) <- c("year4", "location", "variable", "Data")

fishFrame <- reshape2(Fish_ByGear_Ext[,-1], varying=AllFish_ByGear, times=AllFish_ByGear, v.name="Data", direction="long")
names(fishFrame) <- c("year4", "location", "variable", "Data")

bioFrame <- rbind(chlFrame, zoopFrame, fishFrame)


# ============
# = Combined =
# ============
finalFrame0 <- rbind(metFrame, physFrame, chemFrame, bioFrame)
finalFrame_noXform <- finalFrame0
save(file="finalFrame_noXform.RData", finalFrame_noXform)

phRows <- finalFrame0[,"variable"]=="ph"
finalFrame0[phRows,"Data"] <- exp(finalFrame0[phRows,"Data"])


finalFrame <- ddply(finalFrame0, c("variable", "location"), logStat)
save(file="finalFrame_v1.RData", finalFrame)

# ==================
# = Create fatARMA =
# ==================
if(runARMA){
	
	#run ARMA analysis
	system.time(fatARMA0 <- ddply(.data=finalFrame, .variables=c("variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, Method="Evolve"))
	
	#Add the "Type" column, using the All_Params data frame (from TurnExtreme_Fat_Data_v8.R) as a reference
	fatARMA <- data.frame("Type"=NA, fatARMA0)
	catKey <- unique(All_Params[,c("Type", "Variable")])
	uCats <- as.character(unique(catKey[,"Type"]))
	for(i in 1:length(uCats)){
		tCat <- uCats[i] #current category
		tVars <- catKey[catKey[,"Type"]==tCat,"Variable"]#variables in this category
		tI <- is.element(fatARMA[,"variable"], tVars) #rows of fatARMA that correspond to variables in the current category
		fatARMA[tI,"Type"] <- tCat
	}
	
	#save fatARMA
	save(file="fatARMA_v1.RData", fatARMA)
}
