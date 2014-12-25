
rm(list=ls())
graphics.off()

library("reshape")

source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/Data_Functions.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/FatTails_Functions.R")


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

# ===========================================
# = Only use meteorological data after 1980 =
# ===========================================
sm.logic <- (data.max.full[,"Type"]=="Meteorological" & data.max.full[,"year4"]>=1980) | data.max.full[,"Type"]!="Meteorological"
data.max.full <- data.max.full[sm.logic,]


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

	# Saving the time series that couldn't get se's for:
	bad.no.se <- data.fat.full[is.na(data.fat.full[,"se.sh_0"]),]
	bad.no.se.names <- paste(bad.no.se[,"taxID"], bad.no.se[,"location"], bad.no.se[,"variable"])
	save(bad.no.se.names, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/bad.no.se.names.RData")
	
data.fat.full <- data.fat.full[!is.na(data.fat.full[,"se.sh_0"]),] # this is where the number of time series drops from 597 to 595
data.fat.shortMet <- sub.gen(data.fat.full)


save(data.fat.shortMet, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.shortMet.RData")
save.image(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data_sm.RData")