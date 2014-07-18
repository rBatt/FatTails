

rm(list=ls())
graphics.off()
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R")
fnNames <- ls()
pkgNames <- c("DEoptim","GenSA")
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")

library("plyr")
library("doParallel")

# ===========
# = Options =
# ===========
runARMA <- c(TRUE, FALSE)[2]


# =========================
# = Functions to be moved =
# =========================
grabFatARMA <- function(x)x$sigs
grabResidGEV <- function(x) tryCatch({gev.fit(x$Resids)$mle[c("mu_0","sig_0","sh_0")]}, error=function(cond)rep(NA,3))


# ======================
# = Set column classes =
# ======================
data.stat[,"Type"] <- as.character(data.stat[,"Type"])
data.stat[,"taxLvl"] <- as.character(data.stat[,"taxLvl"])
data.stat[,"taxID"] <- as.character(data.stat[,"taxID"])
data.stat[,"location"] <- as.character(data.stat[,"location"])
data.stat[,"variable"] <- as.character(data.stat[,"variable"])
data.stat[,"year4"] <- as.integer(data.stat[,"year4"])
data.stat[,"Data"] <- as.numeric(data.stat[,"Data"])
data.stat[,"N"] <- as.integer(data.stat[,"N"])

finalFrame0 <- data.stat

# ========================================
# = Take the log of variables, except pH =
# ========================================
phRows <- finalFrame0[,"variable"]=="ph"
big0 <- is.na(finalFrame0[,"Data"]) | finalFrame0[,"Data"]>0
finalFrame0[(!phRows)&big0,"Data"] <- log(finalFrame0[(!phRows)&big0,"Data"])

finalFrame <- finalFrame0 # used to be the step where I took the log and made it stationary – already done.
# testFrame <- finalFrame[finalFrame[,"variable"]%in%c("chlor","doc"),]
# test <- ddply(.data=testFrame, .variables=c("Type","taxLvl","taxID","variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, Method="Evolve", .paropts=list(.export=fnNames, .packages=pkgNames))
save(finalFrame, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/finalFrame.RData")

# ===================
# = Fit ARMA Models =
# ===================
if(runARMA){ 
	if(Sys.info()["sysname"]=="Windows"){
		nC <- floor(detectCores()*0.75)
		registerDoParallel(cores=nC)
	}else{
		registerDoParallel()
	}
	
	#run ARMA analysis
	system.time(fatARMA0 <- ddply(.data=finalFrame, .variables=c("Type","taxLvl","taxID","variable", "location"), .fun=ARMAfit, dName="Data", .parallel=TRUE, Method="Evolve", .paropts=list(.export=fnNames, .packages=pkgNames)))
	
	# #Add the "Type" column, using the All_Params data frame (from TurnExtreme_Fat_Data_v8.R) as a reference
	# fatARMA <- data.frame("Type"=NA, fatARMA0)
	# catKey <- unique(data.fat[,c("Type", "variable")])
	# uCats <- as.character(unique(catKey[,"Type"]))
	# for(i in 1:length(uCats)){
	# 	tCat <- uCats[i] #current category
	# 	tVars <- catKey[catKey[,"Type"]==tCat,"variable"]#variables in this category
	# 	tI <- is.element(fatARMA[,"variable"], tVars) #rows of fatARMA that correspond to variables in the current category
	# 	fatARMA[tI,"Type"] <- tCat
	# }
	
	#save fatARMA
	save(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA.RData", fatARMA0)
}else{
	load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA.RData")
}

# ==========================
# = Select best ARMA Model =
# ==========================
fatARMA <- ddply(.data=fatARMA0, .variables=c("Type","taxLvl","taxID","variable", "location"), .fun=minaicc)


# ==============================================
# = Calculate ARMA Variances and get Residuals =
# ==============================================
sigARMA0 <- dlply(fatARMA, .variables=c("Type","taxLvl","taxID","variable", "location"), .fun=getSE, data=finalFrame, .progress="time")
sigARMA00 <- ldply(sigARMA0, grabFatARMA)


# ==================================
# = Compute GEV for ARMA residuals =
# ==================================
residGEV <- ldply(sigARMA0, grabResidGEV, .progress="time")
names(residGEV)[6:8] <- c("residual_mu_0","residual_sig_0","residual_sh_0")

sigARMA <- merge(residGEV, sigARMA00, all=TRUE)
fatARMA1 <- merge(fatARMA, sigARMA, all=TRUE)

boxplot(residual_sh_0~Type, data=fatARMA1)

# ====================================================
# = Merge ARMA+GEV Analysis w/ Original GEV Analysis =
# ====================================================
data.2 <- merge(data.fat, fatARMA1, all.x=TRUE)

# ============================
# = Add ARMA Variance Ratios =
# ============================
InfE <- data.2[,"sigInf"]/data.2[,"sigE"]
Einf <- (data.2[,"sigE"])^2/(data.2[,"sigInf"])^2

data.2[,"InfE"] <- InfE
data.2[,"Einf"] <- Einf

# =============
# = Save data =
# =============
save(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData", data.2)

