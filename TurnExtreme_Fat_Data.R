#Version 3 (02-Spet-2013) I am making these changes after the initial effort for Steve and Monica's class --- previous versions weren't documented because I was under a deadline!
#Version 5 (05-Sept-2013) I've gone through the time series and added two "levels" to the Params object.  The first level is the 2*sd+mean of the long-term time series.  The other is 2*mean.  These summary statistics are not calculated to weight each year equally. Nothing fancy.  I also cleaned up a lot of the old code, exported the Functions to a new script, and exported the figures to a new script.
#_v6 (07-Sept-2013): I want to include the mean and sd in the params data frame.  I've also realized that there are outliers in the times series... but on the low end (e.g., tic in MO), which is why I didn't notice them while scanning for the extremes, b/c all of the extremes I've looked at have been on the high end of the distribution.  So I'm going to add in the mean and sd, and also the 1.2*max( )of the time series, and I might get to cleaning up the time series a bit more.
#_v7 (30-Sept-2013): When I last talked to Steve, he had 2 suggestions for the analysis. 1) If we pick a lower maximum-based threshold (e.g., 1.1*max instead of 1.2), do the disturbingly large return times go away? He didn't like return times of 10^12 years (even though this answer makes sense). 2) Perhaps we should use log-normal as the "straw man" rather than the normal. Normal is very thin-tailed, and people often use log instead, so we're knocking knocking down a straw man that everyone already knows is flimsy. In preparation for my committee meeting tomorrow, I'm going to try to implement these suggestions. I will create arguments to set threshold size, and attempt to use log-normal (the latter will be difficult in time series with many negative values). I might not get to the log-normal in this version.
#_v8 (30-Nov-2013) I think I was previously using Fish1 (electrofishing) for the parameter estimates, and Fish2 (gill net) for the threshold levels. NOT TRUE hmm, I was using electrofishing before, but it seems that gill nets are more fat-tailed. I'm going to stick with electroshocking for now.  However, I still had to make a new version because I had previously overwritten the FishByGear_Ext object (initially for Electrofishing) for the gill net run.
#_v9 (09-Dec-2013) Problems with false 0's (should be NA's) in the Zoop data set. I began investigating why certain variables/ data sets had NA listed for thigns like "N", even though I could compute the ARMA statistics for the extremes (indicating that there was in fact data to be had). It seems like the most common cause was 0's or negative values, which in some cases was combined with an unexpectedly aggressive complete.cases.

rm(list=ls())
graphics.off()
# source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Fat_dGEV.R")
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails")
# load("Data/OrganizedFatData_Read_Fat_Data.RData")
library("wmtsa")
library("reshape")
source("FatTails_Functions.R")

ShowPlots <- c(TRUE, FALSE)[2]
setThresh <- 1.1
threshObs <- 15


# 1 / ((1-0.9) * (79/195)) # 79 is the number of extrema, 195 is the duration of the study (years), and 0.9 is the quantile
# 1 / ( (probability of observing this return level, given that extremes are being observed) * (probability of observing an extreme in the time series)
# 1 / ( (1 - quantile) * (# extrema / duration of time series) )
# So we get a return level of ~182 sunspots per ~25 years.
# Therefore, during a 200-year study, we should have ~8 sunspots that are greater than ~180
# length(which(test[,2]>170)) ---- I'd say that this calculation is approximately correct
# By that I mean that my calculation fo the return level is correct, not necessarily the threshold that I chose before I computed the gev.fit().

# =======================
# = Read in time series =
# =======================
data.max0 <- read.table(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataMax.txt", sep="\t", header=TRUE)


# =================================
# = Subset the zoop and fish data =
# =================================
rmMass <- !data.max0[,"variable"]%in%c("tot_zoop_mass","cpue3_WeiEff")
rmHigher <- (!data.max0[,"variable"]%in%c("density","cpue1_Sum")) | (data.max0[,"variable"]%in%c("density","cpue1_Sum") & data.max0[,"taxLvl"]=="Genus")
data.max <- data.max0[rmMass&rmHigher,]


# ===================
# = make stationary =
# ===================
data.stat00 <- ddply(data.max, c("Type","taxLvl","taxID","location","variable"), function(x)Stationary(x[,c("year4","Data")]))


# ============================================
# = Remove time series with less than 15 obs =
# ============================================
data.stat0N <- ddply(data.stat00, c("Type","taxLvl","taxID","location","variable"), function(x)data.frame("N"=sum(is.finite(x[,"Data"]))))
data.stat0 <- merge(data.stat00, data.stat0N, all=TRUE)
data.stat <- data.stat0[data.stat0[,"N"]>=threshObs,]
row.names(data.stat) <- NULL
data.stat[,"taxID"] <- factor(as.character(data.stat[,"taxID"]))


# =======================================
# = Calculate Return level, duration, N =
# =======================================
data.stat2 <- Inf2NA(ddply(data.stat, c("Type","taxLvl","taxID","location","variable"), calc.level))
# data.stat2 <- data.stat2[data.stat2[,"N"]>=threshObs,]


# ===================================
# = Calculate mean sd (and for log) =
# ===================================
data.stat3 <- ddply(data.stat, c("Type","taxLvl","taxID","location","variable", "N"), musd)
data.stat1 <- merge(data.stat2, data.stat3, all=TRUE)


# ===========
# = Fit GEV =
# ===========
data.gev0 <- ddply(data.stat, c("Type","taxLvl","taxID","location","variable"), GEV)
data.gev <- data.gev0
data.gev[,"Type"] <- factor(data.gev[,"Type"], levels=c("Biological", "Chemical", "Physical", "Meteorological"))
data.gev[,"taxLvl"] <- factor(data.gev[,"taxLvl"], levels=c("Community","Phylum","Class","Order","Family","Genus","Species"))

# ==========================================
# = Calculate return times, combine w/ GEV =
# ==========================================
data.fat00 <- merge(data.gev, data.stat1, all.x=TRUE, all.y=FALSE)

# Return time calculated from GEV
data.fat01 <- ddply(data.fat00, c("Type","taxLvl","taxID","location","variable"), lvl_return, returnFull=TRUE)

# Calculate times for normal and log-normal
data.fat02a <- ddply(data.fat00, c("Type","taxLvl","taxID","location","variable"), normTime)
data.fat02b <- ddply(data.fat00, c("Type","taxLvl","taxID","location","variable"), lnormTime)
data.fat02 <- merge(data.fat02a, data.fat02b, all=TRUE)

# combine l/norm times w/ gev time
data.fat03 <- merge(data.fat01, data.fat02, all=TRUE)



data.fat <- data.fat03[!is.na(data.fat03[,"sh_0"]),]
row.names(data.fat) <- NULL
data.fat[,"shape.sig"] = fattestSig(mu=data.fat[,"sh_0"], se=data.fat[,"se.sh_0"])
data.fat <- data.fat[!is.na(data.fat[,"se.sh_0"]),]




#Split into Bio
bio.gev <- data.gev[data.gev[,"Type"]=="Biological",]

# =============================
# = Fix up Zoop for hierarchy =
# =============================
# Split into Zoop
zoop.gev0 <- bio.gev[bio.gev[,"variable"]%in%c("density","tot_zoop_mass"),]
zoop.gev0[,"taxLvl"] <- factor(zoop.gev0[,"taxLvl"], levels=c("Community","Phylum", "Class", "Order", "Family", "Genus", "Species"))
zoop.gev0[,"taxID"] <- factor(zoop.gev0[,"taxID"])

	# =============================
	# = Add taxa info to Zoop gev =
	# =============================
zTax.id <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatZoopTaxa.txt", sep="\t", header=TRUE, colClasses="character")
zTax.id[zTax.id==""] <- NA

subPhy <- !is.na(zTax.id[,"subphylum"])
zTax.id[subPhy,"phylum"] <- zTax.id[subPhy,"subphylum"]

subClass <- !is.na(zTax.id[,"subclass"])
zTax.id[subClass,"class"] <- zTax.id[subClass,"subclass"]

subSpec <- !is.na(zTax.id[,"species"])
zTax.id[subSpec,"species"] <- paste(substring(zTax.id[subSpec,"genus"], 1, 1), zTax.id[subSpec,"species"], sep=".")

zTax.id[,"Community"] <- "Zoop"

names(zTax.id) <- c("spname", "Phylum", "Subphylum", "Class", "Subclass", "Order", "Family", "Genus", "Species", "Community")
tlvls <- factor(c("Community","Phylum", "Class", "Order", "Family", "Genus", "Species"), levels=c("Community","Phylum", "Class", "Order", "Family", "Genus", "Species"))
	# ===============================================
	# = Combine Zooplankt classification w/ metrics =
	# ===============================================
zoop.gev <- cbind(zoop.gev0, "Community"=NA, "Phylum"=NA, "Class"=NA, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
uztax <- as.character(unique(zoop.gev[,"taxID"]))
uztax2 <- uztax #[!uztax%in%"Zoop"]
for(i in 1:length(uztax2)){
	gevI <- zoop.gev[,"taxID"]==uztax2[i]
	tlev <- unique(zoop.gev[gevI,"taxLvl"])
	
	getID <- as.character(tlvls[as.numeric(tlvls) <= as.numeric(tlev)])
	
	idI <- zTax.id[,as.character(tlev)] == as.character(uztax2[i]) & !is.na(zTax.id[,as.character(tlev)])
	
	if(length(getID)>1){
		zoop.gev[gevI, getID] <- zTax.id[idI,getID][1,]	
	}else{
		zoop.gev[gevI, getID] <- zTax.id[idI,getID][1]	
	}
	
}


# =============================
# = Fix up FISH for hierarchy =
# =============================
fish.gev0 <- bio.gev[bio.gev[,"variable"]%in%c("cpue1_Sum","cpue3_WeiEff"),]
fish.gev0[,"taxLvl"] <- factor(fish.gev0[,"taxLvl"], levels=c("Community","Order", "Family", "Genus", "Species"))
fish.gev0[,"taxID"] <- factor(fish.gev0[,"taxID"])

	# =============================
	# = Add taxa info to Fish gev =
	# =============================
fTax.id <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatFishTaxa.txt", sep="\t", header=TRUE, colClasses="character")
fTax.id[fTax.id==""] <- NA

fsubSpec <- !is.na(fTax.id[,"species"])
fTax.id[fsubSpec,"species"] <- paste(substring(fTax.id[fsubSpec,"genus"], 1, 1), fTax.id[fsubSpec,"species"], sep=".")

fTax.id[,"Community"] <- "Fish"

names(fTax.id) <- c("spname", "Order", "Family", "Genus", "Species", "Subspecies", "Community")
flvls <- factor(c("Community","Order", "Family", "Genus", "Species"), levels=c("Community","Order", "Family", "Genus", "Species"))
	# ===============================================
	# = Combine taxonomic classification w/ metrics =
	# ===============================================
fish.gev <- cbind(fish.gev0, "Community"=NA, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
uztax <- as.character(unique(fish.gev[,"taxID"]))
uztax2 <- uztax #[!uztax%in%"Fish"]
for(i in 1:length(uztax2)){
	gevI <- as.character(fish.gev[,"taxID"])==uztax2[i]
	tlev <- unique(fish.gev[gevI,"taxLvl"])
	
	getID <- as.character(flvls[as.numeric(flvls) <= as.numeric(tlev)])
	
	idI <- fTax.id[,as.character(tlev)] == as.character(uztax2[i]) & !is.na(fTax.id[,as.character(tlev)])
	
	if(length(getID)>1){
		fish.gev[gevI, getID] <- fTax.id[idI,getID][1,]	
	}else{
		fish.gev[gevI, getID] <- fTax.id[idI,getID][1]	
	}
	
}

save.image(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")