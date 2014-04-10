
rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails")
# library("wmtsa")
library("reshape")
source("FatTails_Functions.R")
source("Data_Functions.R")

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
data.max00 <- read.table(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataMax.txt", sep="\t", header=TRUE)


# =================================
# = Subset the zoop and fish data =
# =================================
rmMass <- !data.max00[,"variable"]%in%c("tot_zoop_mass","cpue3_WeiEff")
rmHigher <- (!data.max00[,"variable"]%in%c("density","cpue1_Sum")) | (data.max00[,"variable"]%in%c("density","cpue1_Sum") & data.max00[,"taxLvl"]%in%c("Species","Genus","Family","Order","Class","Phylum"))
data.max0 <- data.max00[rmMass&rmHigher,]

# =================================
# = Add in NA's for missing years =
# =================================
data.max.full <- ddply(data.max0, c("Type","taxLvl","taxID","location","variable"), fillMiss)
data.max <- sub.gen(data.max.full)

# ===================
# = make stationary =
# ===================
data.stat00 <- ddply(data.max.full, c("Type","taxLvl","taxID","location","variable"), function(x)Stationary(x[,c("year4","Data")]))

# ============================================
# = Remove time series with less than 15 obs =
# ============================================
data.stat0N <- ddply(data.stat00, c("Type","taxLvl","taxID","location","variable"), function(x)data.frame("N"=sum(is.finite(x[,"Data"]))))
data.stat0 <- merge(data.stat00, data.stat0N, all=TRUE)
data.stat.full <- data.stat0[data.stat0[,"N"]>=threshObs,]
row.names(data.stat.full) <- NULL
data.stat.full[,"taxID"] <- factor(as.character(data.stat.full[,"taxID"]))

data.stat <- sub.gen(data.stat.full)


# =======================================
# = Calculate Return level, duration, N =
# =======================================
data.stat2 <- Inf2NA(ddply(data.stat.full, c("Type","taxLvl","taxID","location","variable"), calc.level))
# data.stat2 <- data.stat2[data.stat2[,"N"]>=threshObs,]


# ===================================
# = Calculate mean sd (and for log) =
# ===================================
data.stat3 <- ddply(data.stat.full, c("Type","taxLvl","taxID","location","variable", "N"), musd)
data.stat1 <- merge(data.stat2, data.stat3, all=TRUE)


# ===========
# = Fit GEV =
# ===========
data.gev <- ddply(data.stat.full, c("Type","taxLvl","taxID","location","variable"), GEV)
# data.gev <- data.gev0
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



data.fat.full <- data.fat03[!is.na(data.fat03[,"sh_0"]),]
row.names(data.fat.full) <- NULL
data.fat.full[,"shape.sig"] = fattestSig(mu=data.fat.full[,"sh_0"], se=data.fat.full[,"se.sh_0"])
data.fat.full <- data.fat.full[!is.na(data.fat.full[,"se.sh_0"]),]
data.fat <- sub.gen(data.fat.full)

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
zoop.gev.full <- cbind(zoop.gev0, "Community"=NA, "Phylum"=NA, "Class"=NA, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
uztax <- as.character(unique(zoop.gev.full[,"taxID"]))
uztax2 <- uztax #[!uztax%in%"Zoop"]
for(i in 1:length(uztax2)){
	gevI <- zoop.gev.full[,"taxID"]==uztax2[i]
	tlev <- unique(zoop.gev.full[gevI,"taxLvl"])
	
	getID <- as.character(tlvls[as.numeric(tlvls) <= as.numeric(tlev)])
	
	idI <- zTax.id[,as.character(tlev)] == as.character(uztax2[i]) & !is.na(zTax.id[,as.character(tlev)])
	
	if(length(getID)>1){
		zoop.gev.full[gevI, getID] <- zTax.id[idI,getID][1,]	
	}else{
		zoop.gev.full[gevI, getID] <- zTax.id[idI,getID][1]	
	}
	
}
zoop.gev <- sub.gen(zoop.gev.full)


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
fish.gev.full <- cbind(fish.gev0, "Community"=NA, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
uztax <- as.character(unique(fish.gev.full[,"taxID"]))
uztax2 <- uztax #[!uztax%in%"Fish"]
for(i in 1:length(uztax2)){
	gevI <- as.character(fish.gev.full[,"taxID"])==uztax2[i]
	tlev <- unique(fish.gev.full[gevI,"taxLvl"])
	
	getID <- as.character(flvls[as.numeric(flvls) <= as.numeric(tlev)])
	
	idI <- fTax.id[,as.character(tlev)] == as.character(uztax2[i]) & !is.na(fTax.id[,as.character(tlev)])
	
	if(length(getID)>1){
		fish.gev.full[gevI, getID] <- fTax.id[idI,getID][1,]	
	}else{
		fish.gev.full[gevI, getID] <- fTax.id[idI,getID][1]	
	}
	
}
fish.gev <- sub.gen(fish.gev.full)

save.image(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")