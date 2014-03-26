#Fat tail for scenarios project/class.
#Organize the raw data into something useful.
#The goals of this script is to end up with data of which I can then take "block extrema".
#Version 0.0.0
#18-April-2013

#_v1 (02-Sept-2013): Transferring scripts from the scenario course to the FatTails folder.  I want to streamline or organize the progress I made in the class by deleting code that I don't want to use anymore, and by isolating old code into separate scripts and allow the analysis to progress in new scripts.
#_v2 (10-Sept-2013): Scanning the full time series for outliers or clearly bad values

#NOTE: "ERGASILUS" are known as "gill lice".  An interesting "fat tail" discussion point, potentially.  An outbreak of these suckers would be interesting.

rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_rawData")
library("plyr")
library("reshape")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Fat_dGEV_v1.0.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/CalcZoopBiomass_v2.R")

manClean <- function(x, varCols){
	x2 <- NULL
	dev.new(width=10, height=4)
	par(mar=c(3,3,3,0.5))
	for(i in 1:length(varCols)){
		tvar <- varCols[i]
		for(l in 1:length(unique(x[,"lakeid"]))){
			tlake <- unique(x[,"lakeid"])[l]
			tindex <- x[,"lakeid"]==tlake
			tdat <- x[tindex,tvar]
			if(any(!is.na(tdat))){
				epiInfo <- is.element(c("depth", "Zmix"), names(x))
				if(all(epiInfo)){ #if there is information about the epi depth, color the points red for epi, and NA otherwise
					MyRed <- rgb(t(col2rgb("red")), alpha=40, maxColorValue=255)
					inEpi <- x[tindex,"depth"]<=x[tindex,"Zmix"] | is.na(x[tindex,"Zmix"]) & x[tindex,"depth"] <=2
					plotCols <- c(NA, MyRed)[inEpi+1]
				}else{
					plotCols <- rgb(t(col2rgb("blue")), alpha=40, maxColorValue=255) #if the necessary columns for determining epi or not are missing, color blue
				}
				plot(tdat, main=paste(tvar, tlake), pch=21, bg=plotCols)
				bad <- identify(tdat)
				if(length(bad)>0){
					x2 <- c(x2, paste(tlake, tvar, x[x[,"lakeid"]==tlake,tvar][bad]))
					x[x[,"lakeid"]==tlake,tvar][bad] <- NA
				}
			}else{
				next
			}
		}
	}
	dev.off()
	return(list(x, x2))	
}

deemedBad <- list()

# Need to remove spuries "species" from Zoop data. (e.g., UNID)

#The Met data doesn't include wind speeds

# For the analysis of this data, I'm starting to think that I might want to use annual maxima as a solid first-cut.  This will simplify issues of irregular time series, making it easier to computer return-rate and return-levels.


# ======================
# = flag2NA, added _v2 =
# ======================
flag2NA <- function(x, summarizeBad=TRUE){
	flagNames <- names(x)[regexpr("^flag", names(x))==1] #which column names start with 'flag' ?
	flagVars <- gsub("flag", "", flagNames) #what are the column names that have a corresponding "flag___" column?
	hasFlag <- function(x){ #function to detect which elements contain an uppercase character (i.e., which have a flag value)
		grep("[:ALPHA:]", x)
	}
	bad <- apply(x[,flagNames], 2, hasFlag) #for each flag column, which rows are flagged?
	if(summarizeBad){
		print(summary(bad))	
	}
	for(i in 1:length(flagVars)){
		tbadRow <- do.call('$', list(bad, flagNames[i])) #which rows are bad in the ith flag column?
		if(length(tbadRow)>0){
			x[tbadRow, flagVars[i]] <- NA #replace those flagged elements in the parent column with NA
		}
	}
	x #return the data frame with flagged values replaced with NA's
}

# ======================
# = Physical Limnology =
# ======================
Phys000 <- read.csv("Phys_NrtnSrtn.csv")

zmix <- function(x){
		temp <- x[,"wtemp"]
		depth <- x[,"depth"]
        tempDiff <- temp[-length(temp)] - temp[-1]
        depthDiff <- depth[-1] - depth[-length(depth)]
        ratio <- tempDiff / depthDiff
		zhat <- depth[which(ratio >= 2)][1]
		# Zhat2 <- min(which(zhat>0))
		# Z <- Zhat2
		Z <- ifelse(!(zhat>0), NA, zhat)
       return(Z)
}
Phys00 <- subset(Phys000, sta==1 & rep==1)
Phys00 <- flag2NA(Phys00)
Phys00_clean <- manClean(Phys00, c("wtemp", "o2", "o2sat", "frlight"))
Phys00 <- Phys00_clean[[1]]
deemedBad$Phys <- Phys00_clean[[2]]


Zmix <- ddply(.data=Phys00, .variables=c("lakeid", "year4", "daynum"), .fun=zmix)
names(Zmix) <- c("lakeid", "year4", "daynum", "Zmix")
Zmix[,2] <- as.numeric(as.character(Zmix[,2]))
Zmix[,3] <- as.numeric(as.character(Zmix[,3]))
Zmix[,4] <- as.numeric(as.character(Zmix[,4]))

Phys0 <- merge(Phys00, Zmix, all=TRUE)
zd <- Phys0[,"Zmix"] - Phys0[,"depth"]
Phys0[,"InEpi"] <- as.numeric(zd>0 & !is.na(zd))
Phys0[,"Mixing"] <- as.numeric(is.na(zd))
Phys0[,"InMetaHypo"] <- as.numeric(Phys0[,"InEpi"]+Phys0[,"Mixing"] == 0)
Phys0 <- subset(Phys0, InMetaHypo==0, select=c("lakeid", "year4", "daynum", "sampledate", "depth", "wtemp", "o2", "o2sat", "deck", "light", "frlight", "Zmix", "InEpi", "Mixing"))

EpiMean <- function(x){
	xnames <- names(x)[!is.element(names(x), c("lakeid", "year4", "daynum", "sampledate"))]
	colMeans(x[,xnames], na.rm=TRUE)
}

Phys <- ddply(.data=Phys0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean)


# =======================
# = Chlorophyll Data =
# =======================
Chl000n <- read.csv("Chla_Nrtn.csv")
ChlorFlags <- strsplit("A AJ AL B BG BK BL D G GJ GK GL H I J JK JKL JL K L KL KO L O LB LG LK LO O OJ", split=" ")[[1]]
BadChlorN <- which(is.element(Chl000n[,"flagchlor"], ChlorFlags) | Chl000n[,"chlor"]<0)
PhaeoFlags <- strsplit("A AJ AK AL B BG BK BL D G GB GK GL H I J JK JKL JL K K L KL KLO KO L L K LB LG LK LO O OJ OK", split=" ")[[1]]
BadPhaeoN <- which(is.element(Chl000n[,"flagphaeo"], PhaeoFlags) | Chl000n[,"phaeo"]<0)
Chl000n[union(BadChlorN, BadPhaeoN),"chlor"] <- NA
Chl000n[union(BadChlorN, BadPhaeoN), "phaeo"] <- NA

# Chl000n <- manClean(Chl000n, c("chlor", "phaeo"))

#I'm not going to bother using the flag2NA function here, b/c I implemented a similar, more manual, method in the first version (before the flag2NA function was created)

Chl000s <- read.csv("Chla_Srtn.csv")
Chl000s[,"daynum"] <- as.numeric(as.character(format.Date(Chl000s[,"sampledate"], format="%j")))
nummean <- function(x){mean(as.numeric(x))}
Chl000s[,"depth"] <- unlist(lapply(strsplit(as.character(Chl000s[,"depth_range_m"]), split="-"), FUN=nummean))
BadChlorS <- which(is.element(Chl000s[,"flag_fluor"], c("C", "E", "CE", "EC")) | Chl000s[,"correct_chl_fluor"]<0)
Chl000s[BadChlorS,"correct_chl_fluor"] <- NA
Chl000s <- Chl000s[,c("lakeid", "year4", "sampledate", "daynum", "depth", "rep", "correct_chl_fluor")]
names(Chl000s) <- c("lakeid", "year4", "sampledate", "daynum", "depth", "rep", "chlor")

# Chl000s <- manClean(Chl000s, "chlor")

EpiMean2 <- function(x){
	stopifnot(all(is.element(c("depth","Zmix"), colnames(x))))
	xnames <- colnames(x)[!is.element(colnames(x), c("lakeid", "year4", "daynum", "sampledate"))]
	x2 <- x[,xnames]

	if(all(is.na(x2[,"Zmix"]))){ #if zmix is unknown, just return the shallowest values
		if(is.null(nrow(x2))){return(data.frame(x2))} #if zmix is unknown and there is only one depth of values, just return this depth
		
		xnames2 <- colnames(x2)[!is.element(colnames(x2), c("depth", "Zmix"))]
		if(all(apply(X=x2[,xnames2], MARGIN=2, FUN=is.na))){return(colMeans(x2, na.rm=TRUE))} #if every variable is NA for all depths, just return the mean of all columns (shortcut for getting the mean depth and retaining the format of the data frame)
		shallDepths <- c()
		# ShalVal <- data.frame("Zmix"=NA)
		ShalVal <- matrix(data=NA, ncol=length(xnames), dimnames=list(NULL, xnames))
		# ShalVal[,"Zmix"] <- NA #probably don't need this line
		for(i in 1:length(xnames2)){ #this elaborate loop is necessary b/c the shallowest non-NA value might not be the same for each column
			x_notNA_ind <- which(!is.na(x2[,xnames2[i]]))
			if(length(x_notNA_ind)==0){ #if all of the values in this column (i.e., this variable) are NA, record an NA for the depth and for the variable
				shallDepths[i] <- NA
				ShalVal[,xnames2[i]] <- NA
			}else{ #if this column contains non-NA values, record the shallowest depth at which the varialbe is non-NA, and record the variable's value at the aforementioned depth
				x_notNA <- x2[x_notNA_ind,]
				shallDepths[i] <- min(x_notNA[, "depth"]) # I'm intentionally leaving out na.rm=TRUE b/c none of these values should be NA
				ShalVal[,xnames2[i]] <- x_notNA[which.min(x_notNA[, "depth"]), xnames2[i]]
			}
		}
		UniqueShals <- unique(na.omit(shallDepths))
		if(length(UniqueShals)==1){
			DeepestShall <- UniqueShals
		}else{
			DeepestShall <- max(UniqueShals) #paste("<=", max(UniqueShals))
		}
		
		ShalVal[,"depth"] <- DeepestShall
		# ShalVal2 <- matrix(ShalVal, ncol=length(xnames), dimnames=list(NULL, xnames))
		return(ShalVal) #if zmix is unknown and values are at multiple depths, 
		# shall <- which.min(x2[,"depth"])
		# x3 <- colMeans(x2[shall,], na.rm=TRUE)
	}
	
	zd <- x2[,"Zmix"] - x2[,"depth"]
	InEpi <- as.numeric(zd>0 & !is.na(zd)) #if zmix is known, and values were from epi, return the mean of the epi
	if(any(InEpi==1)){
		return(colMeans(x2[which(InEpi==1),], na.rm=TRUE))
	}
	return(colMeans(x2, na.rm=TRUE)) #if zmix is known, but no values were in the epi, just return the mean of all values that were recorded.
}


Chl000 <- merge(Chl000n, Chl000s, all=TRUE)
Chl00 <- merge(Phys[,c("lakeid", "year4", "daynum", "sampledate", "Zmix")], Chl000, all=TRUE)
Chl00 <- Chl00[, c("lakeid", "year4", "daynum", "sampledate", "Zmix", "rep", "depth", "chlor", "sta", "phaeo")]
# Chl00 <- ddply(.data=Chl00, .variables=c("lakeid", "year4"), .fun=InterpZ) 
length(which(is.na(Chl00[,"Zmix"])))
# aggregate(Chl00[,"Zmix"], by=list(Chl00[,"daynum"], Chl00[,"year4"], Chl00[,"lakeid"]), FUN=mean, na.rm=TRUE)
# aggregate(Phys[,"Zmix"], by=list(Phys[,"daynum"], Phys[,"year4"], Phys[,"lakeid"]), FUN=mean, na.rm=TRUE)
Chl00 <- aggregate(Chl00[,c("Zmix", "chlor", "phaeo")], by=Chl00[,rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)
Chl0 <- ddply(.data=Chl00, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
Chl <- Chl0[,c("lakeid", "year4", "daynum", "sampledate", "depth", "chlor", "phaeo")]
names(Chl) <- c("lakeid", "year4", "daynum", "sampledate", "depth_chlor", "chlor", "phaeo")
# class(Chl[,"sampledate"]) <- "Date"


# ====================
# = Zooplankton Data =
# ====================
#not using flag2NA b/c these files don't have flag columns
Zoop000s1 <- read.csv("ZoopDens_Old_Srtn.csv")#density recorded as individuals/m^2
Zoop000s1_tow <- read.csv("ZoopDens_TowDepth_Old_Srtn.csv")#Waubesa = 9.75m; Kegonsa = 7.75m;
Zoop000s1 <- merge(Zoop000s1, Zoop000s1_tow, all.x=TRUE)
Zoop000s1[which(Zoop000s1[,"lakeid"]=="WA"),"towdepth"] <- 9.75
Zoop000s1[which(Zoop000s1[,"lakeid"]=="KE"),"towdepth"] <- 7.75
mean_td_ME <- mean(Zoop000s1[which(Zoop000s1[,"lakeid"]=="ME"),"towdepth"], na.rm=TRUE) #the average tow depth for mendota
mean_td_MO <- mean(Zoop000s1[which(Zoop000s1[,"lakeid"]=="MO"),"towdepth"], na.rm=TRUE) #the average tow depth for monona
na_td_ME <- which(Zoop000s1[,"lakeid"]=="ME" & is.na(Zoop000s1[,"towdepth"])) #which rows were missing a tow depth for mendota?
na_td_MO <- which(Zoop000s1[,"lakeid"]=="MO" & is.na(Zoop000s1[,"towdepth"])) #which rows were missing a tow depth for monona?
Zoop000s1[na_td_ME,"towdepth"] <- mean_td_ME #replace missing mendota tow depths with average mendota tow depth
Zoop000s1[na_td_MO,"towdepth"] <- mean_td_MO #replace missing monona tow depths with average monona tow depth
Zoop00s1 <- Zoop000s1
Zoop00s1[,"density"] <- (Zoop000s1[,"density"]/Zoop000s1[,"towdepth"])*0.001 #convert zooplankton density to individuals/L
names(Zoop00s1) <- c("lakeid", "year4", "sampledate", "taxon", "species_code", "density", "avg_length", "comments", "towdepth")
# Zoop00s1 <- manClean(Zoop00s1, c("density", "avg_length")) #no problems


Zoop000s2 <- read.csv("ZoopDens_New_Srtn.csv")#density recorded as individuals/m^2
Zoop00s2 <- Zoop000s2
Zoop00s2[,"density"] <- (Zoop000s2[,"density"]/Zoop000s2[,"towdepth"])*0.001 #convert zooplankton density to individuals/L
names(Zoop00s2) <- c("lakeid", "year4", "sampledate", "station", "towdepth", "species_code", "taxon", "density", "indiv_measured", "avg_length")
# Zoop00s2 <- manClean(Zoop00s2, c("density", "avg_length")) #no problems

Zoop00s <- merge(Zoop00s1, Zoop00s2, all=TRUE)
Zoop00s[,"daynum"] <- as.numeric(as.character(format.Date(Zoop00s[,"sampledate"], format="%j")))
for(i in 1:length(unique(Zoop00s[,"taxon"]))){
	TaxonRows <- which(is.element(Zoop00s[,"taxon"], unique(Zoop00s[,"taxon"])[i]))
	Zoop00s[TaxonRows, "avg_zoop_mass"] <- ZoopMass(Zoop00s[TaxonRows,])
}
Zoop00s[, "tot_zoop_mass"] <- Zoop00s[,"density"] * Zoop00s[,"avg_zoop_mass"]
Zoop0s <- Zoop00s[,c("lakeid", "year4", "daynum", "sampledate", "taxon", "density", "avg_length", "avg_zoop_mass", "tot_zoop_mass")]


Zoop000n <- read.csv("ZoopDens_Nrtn.csv")#density recorded as individuals/Liter
names(Zoop000n) <- c(c("lakeid", "year4", "sampledate", "station", "species_code", "taxon", "density", "indiv_measured", "avg_length"))
Zoop000n[,"daynum"] <- as.numeric(as.character(format.Date(Zoop000n[,"sampledate"], format="%j")))
for(i in 1:length(unique(Zoop000n[,"taxon"]))){
	TaxonRows <- which(is.element(Zoop000n[,"taxon"], unique(Zoop000n[,"taxon"])[i]))
	Zoop000n[TaxonRows, "avg_zoop_mass"] <- ZoopMass(Zoop000n[TaxonRows,])
}
Zoop000n[, "tot_zoop_mass"] <- Zoop000n[,"density"] * Zoop000n[,"avg_zoop_mass"]
Zoop000n[,"daynum"] <- as.numeric(as.character(format.Date(Zoop000n[,"sampledate"], format="%j")))
Zoop0n <- Zoop000n[,c("lakeid", "year4", "daynum", "sampledate", "taxon", "density", "avg_length", "avg_zoop_mass", "tot_zoop_mass")]
# Zoop0n <- manClean(Zoop0n, c("avg_length", "avg_zoop_mass", "tot_zoop_mass")) #no problems

Zoop0 <- merge(Zoop0n, Zoop0s, all=TRUE) #this contains zooplankton information broken down by taxon
Zoop <- aggregate(Zoop0[,c("density", "avg_length", "avg_zoop_mass", "tot_zoop_mass")], by=Zoop0[,c("year4", "daynum","sampledate", "lakeid")], FUN=sum, na.rm=TRUE)
# class(Zoop[,"sampledate"]) <- "Date"

# ========
# = Ions =
# ========
Ions000 <- flag2NA(read.csv("Ions_NrtnSrtn.csv")) #where are there so many rows with no data?
Ions00 <- merge(Ions000, Zmix, all.x=TRUE)
# Ions00 <- merge(x=Ions00, y=Phys[,c("lakeid", "year4", "daynum", "sampledate", "Zmix")], all.x=TRUE) #I'm not sure why I merged with zmix and with Phys[].  I'll leave it for now, but I might need to do one or the other.
IonTypes <- c("cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")
# IonFlags <- paste("flag", IonTypes, sep="")
# for(i in 1:length(IonFlags)){
# 	badions <- which(Ions00[,IonFlags[i]]!="")
# 	Ions00[badions,IonTypes[i]] <- NA
# }
# nonflag <- names(Ions00)[!is.element(names(Ions00), IonFlags)]
# Ions0 <- Ions00[,nonflag]
#I need to take the average for the rep's and stations first
allna <- function(x){all(is.na(x))}
Ions_Not_All_NA_Ind <- which(!apply(Ions00[,c("cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")], MARGIN=1, FUN=allna))
Ions00 <- Ions00[Ions_Not_All_NA_Ind,]
Ions0 <- aggregate(Ions00[,c("Zmix", "cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")], by=Ions00[,rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)

Ions0_clean <- manClean(Ions0, IonTypes)
Ions0 <- Ions0_clean[[1]]
deemedBad$Ions <- Ions0_clean[[2]]

Ions <- ddply(.data=Ions0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
names(Ions) <- c("lakeid", "year4", "daynum", "sampledate", "depth_ions", "zmix", "cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")
# class(Ions[,"sampledate"]) <- "Date"

# ==================
# = Fish Abundance =
# ==================
# c("VGN019", "VGN025", "VGN032", "VGN038", "VGN051", "VGN064", "VGN089", "VGN127", "VGN")
# gsub("[0123456789]", "", c("VGN019", "VGN025", "VGN032", "VGN038", "VGN051", "VGN064", "VGN089", "VGN127", "VGN"))
Fish000_abun <- read.csv("Fish_Abun_NrtnSrtn.csv")
Fish000_abun[,"gearid"] <- gsub("[0123456789]", "", Fish000_abun[,"gearid"]) #hell yeah, i learned a new function. (well, how to use reular expressions a bit better)
Fish000_abun[,"spname"] <- gsub(" ", "", Fish000_abun[,"spname"]) #remove white space from species names
BsSpecies <- c("GUPPY", "VIRILIS", "RUSTICUS", "PROPINQUUS", "UNIDCHUB", "UNIDDARTER", "UNIDENTIFIED", "UNIDMINNOW", "CRAYFISH", "DARTER", "LARVALFISH")
Fish00_abun <- subset(Fish000_abun, !is.element(spname, BsSpecies))

Fish00_abun[,"cpue1_Sum"] <- Fish00_abun[,"total_caught"]/Fish00_abun[,"effort"]
Fish0_abun <- aggregate(Fish00_abun[,c("effort", "total_caught", "cpue1_Sum")], by=list(Fish00_abun[,"spname"], Fish00_abun[,"lakeid"], Fish00_abun[,"year4"]), FUN=sum)
TotFish0_abun <- aggregate(Fish00_abun[,c("total_caught", "cpue1_Sum")], by=list(Fish00_abun[,"gearid"], Fish00_abun[,"lakeid"], Fish00_abun[,"year4"]), FUN=sum)
Fish0_abun[,"cpue2_AvgMeths"] <- Fish0_abun[,"total_caught"]/Fish0_abun[,"effort"]
names(Fish0_abun) <- c("spname", "lakeid", "year4", "effort", "total_caught", "cpue1_Sum", "cpue2_AvgMeths")
names(TotFish0_abun) <- c("gearid", "lakeid", "year4", "total_caught", "cpue_SumSpecies")

# =============
# = Fish Size =
# =============
Fish000_size <- read.csv("Fish_LenWei_NrtnSrtn.csv")
Fish000_size[,"spname"] <- gsub(" ", "", Fish000_size[,"spname"]) 
Fish000_size[,"gearid"] <- gsub("[0123456789]", "", Fish000_size[,"gearid"])
Fish00_size <- subset(Fish000_size, !is.element(spname, BsSpecies))
FixPike <- which(is.element(Fish00_size[,"spname"], "SILVERPIKE"))
Fish00_size[FixPike, "spname"] <- "NORTHERNPIKE"

SizSumry <- function(x){
	Nfish <- length(x[,"spname"])
	
	exLeng <- x[,"length"][!is.na(x[,"length"])]
	exWei <- x[,"weight"][!is.na(x[,"weight"])]
	
	Nleng <- length(exLeng)
	Nwei <- length(exWei)
	
	if(Nleng==0){
		maxLeng <- NA
		minLeng <- NA
	}else{
		maxLeng <- max(exLeng, na.rm=TRUE)
		minLeng <- min(exLeng, na.rm=TRUE)
	}
	
	if(Nwei==0){
		maxWei <- NA
		minWei <- NA
	}else{
		maxWei <- max(exWei, na.rm=TRUE)
		minWei <- min(exWei, na.rm=TRUE)
	}
	
	sumryLeng <- c("mean_"=mean(exLeng, na.rm=TRUE), "max_"=maxLeng, "min_"=minLeng) #summary(x[,"length"])[[c(1,4,6)]]
	names(sumryLeng) <- paste(names(sumryLeng), "Leng", sep="")
	sumryWei <- c("mean_"=mean(exWei, na.rm=TRUE), "max_"=maxWei, "min_"=minWei) #summary(x[,"weight"])[[c(1,4,6)]]
	names(sumryWei) <- paste(names(sumryWei), "Wei", sep="")
	
	x2names <- c(names(sumryLeng), names(sumryWei), "Nfish", "Nleng", "Nwei")
	# x2 <- matrix(data=c(sumryLeng, sumryWei, Nfish, Nleng, Nwei), nrow=1, dimnames=list(NULL, x2names))
	x2 <- c(sumryLeng, sumryWei, Nfish, Nleng, Nwei)
	names(x2) <- x2names
	return(x2)
}
Fish0_size <- ddply(Fish00_size, .variables=c("lakeid", "year4", "spname"), .fun=SizSumry)
# names(Fish0_size) <- c("lakeid", "year4", "sampledate", "spname", "sampletype", "depth", "rep", "indid", "weight", "sex", "fishpart", "spseq")

TotFish0_size <- ddply(Fish00_size, .variables=c("lakeid", "year4", "spname", "gearid"), .fun=SizSumry)
# names(TotFish0_size) <- c("lakeid", "year4", "sampledate", "gearid", "spname", "sampletype", "depth", "rep", "indid", "weight", "sex", "fishpart", "spseq")

# ================
# = Combine Fish =
# ================
Fish <- merge(Fish0_abun, Fish0_size, all=TRUE, by=c("lakeid", "year4", "spname"))
Fish[,"SumWei"] <- Fish[,"Nwei"]*Fish[,"mean_Wei"]
Fish[,"SumLeng"] <- Fish[,"Nleng"]*Fish[,"mean_Leng"]
Fish[,"cpue3_WeiEff"] <- Fish[,"SumWei"]/Fish[,"effort"]
Fish[,"cpue4_LengEff"] <- (Fish[,"Nleng"]*Fish[,"mean_Leng"])/Fish[,"effort"]

# TotFish <- merge(TotFish0_abun, TotFish0_size, all=TRUE, by=c("lakeid", "year4", "spname"))
TotFish <- merge(Fish00_abun, TotFish0_size, all=TRUE, by=c("lakeid", "year4", "spname", "gearid"))
TotFish[,"SumWei"] <- TotFish[,"Nwei"]*TotFish[,"mean_Wei"]
TotFish[,"SumLeng"] <- TotFish[,"Nleng"]*TotFish[,"mean_Leng"]
TotFish[,"cpue3_WeiEff"] <- TotFish[,"SumWei"]/TotFish[,"effort"]
TotFish[,"cpue4_LengEff"] <- (TotFish[,"Nleng"]*TotFish[,"mean_Leng"])/TotFish[,"effort"]

# FishCats <- c("total_caught", "cpue1_Sum", "mean_Leng", "max_Leng", "min_Leng", "mean_Wei", "max_Wei", "min_Wei", "Nfish", "Nleng", "Nwei", "SumWei", "SumLeng", "cpue3_WeiEff", "cpue4_LengEff")
FishCats1 <- c("total_caught", "cpue1_Sum", "Nfish", "Nleng", "Nwei", "SumWei", "SumLeng", "cpue3_WeiEff", "cpue4_LengEff")
FishCats2 <- c("mean_Leng", "mean_Wei")
FishCats3 <- c("max_Leng", "max_Wei")
FishCats4 <- c("min_Leng", "min_Wei")

# Inf2NA <- function(x) {x[which(x==-Inf | x==Inf, arr.ind=TRUE)] <- NA; x}
Inf2NA <- function(x) {x[x==-Inf | x==Inf] <- NA; x}

Fish_BySpec1 <- aggregate(TotFish[, FishCats1], by=TotFish[,c("spname", "year4", "lakeid")], sum, na.rm=TRUE)
Fish_BySpec2 <- aggregate(TotFish[, FishCats2], by=TotFish[,c("spname", "year4", "lakeid")], mean, na.rm=TRUE)
Fish_BySpec3 <- aggregate(TotFish[, FishCats3], by=TotFish[,c("spname", "year4", "lakeid")], max, na.rm=TRUE)
Fish_BySpec4 <- aggregate(TotFish[, FishCats4], by=TotFish[,c("spname", "year4", "lakeid")], min, na.rm=TRUE)
Fish_BySpec <- merge_recurse(list(Fish_BySpec1, Fish_BySpec2, Fish_BySpec3, Fish_BySpec4), all.x=TRUE, all.y=TRUE)
Fish_BySpec <- Inf2NA(Fish_BySpec)

Fish_ByGear1 <- aggregate(TotFish[, FishCats1], by=TotFish[,c("gearid", "year4", "lakeid")], sum, na.rm=TRUE)
Fish_ByGear2 <- aggregate(TotFish[, FishCats2], by=TotFish[,c("gearid", "year4", "lakeid")], mean, na.rm=TRUE)
Fish_ByGear3 <- aggregate(TotFish[, FishCats3], by=TotFish[,c("gearid", "year4", "lakeid")], max, na.rm=TRUE)
Fish_ByGear4 <- aggregate(TotFish[, FishCats4], by=TotFish[,c("gearid", "year4", "lakeid")], min, na.rm=TRUE)
Fish_ByGear <- merge_recurse(list(Fish_ByGear1, Fish_ByGear2, Fish_ByGear3, Fish_ByGear4), all.x=TRUE, all.y=TRUE)
Fish_ByGear <- Inf2NA(Fish_ByGear)

Fish_GearSpec1 <- aggregate(TotFish[, FishCats1], by=TotFish[,c("gearid","spname", "year4", "lakeid")], sum, na.rm=TRUE)
Fish_GearSpec2 <- aggregate(TotFish[, FishCats2], by=TotFish[,c("gearid","spname" ,"year4", "lakeid")], mean, na.rm=TRUE)
Fish_GearSpec3 <- aggregate(TotFish[, FishCats3], by=TotFish[,c("gearid","spname", "year4", "lakeid")], max, na.rm=TRUE)
Fish_GearSpec4 <- aggregate(TotFish[, FishCats4], by=TotFish[,c("gearid","spname" ,"year4", "lakeid")], min, na.rm=TRUE)
Fish_GearSpec <- merge_recurse(list(Fish_GearSpec1, Fish_GearSpec2, Fish_GearSpec3, Fish_GearSpec4), all.x=TRUE, all.y=TRUE)
Fish_GearSpec <- Inf2NA(Fish_GearSpec)


# ==================
# = Meteorological =
# ==================
# http://lter.limnology.wisc.edu/data/filter/11001
# above link for woodruff airport data from 1989 on, with wind and par etc.
Met000_Mad <- read.csv("MetDat_Madison.csv")
Met00_Mad <- Met000_Mad[,c("year4", "daynum", "sampledate", "max_air_temp_adjusted", "min_air_temp_adjusted", "ave_air_temp_adjusted", "range_air_temp_adjusted", "precip_raw_mm", "snow_raw_cm", "snow_depth_cm", "data_status")]
names(Met00_Mad) <- c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "ave_air_temp", "range_air_temp", "precip_mm", "snow_cm", "snow_depth_cm", "data_status")
Met00_Mad[,"location"] <- "Madison"

Met000_Min <- read.csv("MetDat_Minocqua.csv")
Met000_Min[,"daynum"] <- as.integer(as.character(format.Date(as.Date(Met000_Min[,"sampledate"]), format="%j")))
# Met000_Min[,"daynum"] <- as.numeric(as.character(format.Date(Met000_Min[,"sampledate"], format="%j"))) #interesting, the above line runs **a lot** faster than this one....
Met000_Min[,"range_air_temp"] <- Met000_Min[,"max_air_temp"] - Met000_Min[,"min_air_temp"]
Met00_Min <- Met000_Min[,c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "range_air_temp", "precip", "snow", "snow_depth", "data_status")]
names(Met00_Min) <- c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "range_air_temp", "precip_mm", "snow_cm", "snow_depth_cm", "data_status")
Met00_Min[,"location"] <- "Minocqua"
Met00_Min[,"snow_cm"] <- Met00_Min[,"snow_cm"]/10 #the minocqua snowfall was originally in mm

Met0 <- merge(Met00_Mad, Met00_Min, all=TRUE)
Met <- Met0[,c("location", "year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "ave_air_temp", "range_air_temp", "precip_mm", "snow_cm", "snow_depth_cm", "data_status")]
# class(Met[,4]) <- "Date"

# ==============
# = Lake Level =
# ==============
LakLev000 <- read.csv("LakeLevel_Nrtn.csv") # WHY ONLY NORTHERN?; level is in "meters above sea level"
LakLev <- LakLev000[,c("lakeid", "year4", "daynum", "sampledate", "llevel_elevation")]
names(LakLev) <- c("lakeid", "year4", "daynum", "sampledate", "LakeLevel")


# ======================
# = Chemical Limnology =
# ======================
Chem000 <- flag2NA(read.csv("Chem_NrtnSrtn.csv"))
# chemClasses <- c("factor", rep("integer", 2), "Date", "numeric", rep("integer", 2), rep("numeric", 24), rep("character",25))
# chemNums <- which(chemClasses=="numeric")
# for(i in chemNums){
# 	Chem000[,i] <- as.numeric(as.character(Chem000[,i]))
# 	Chem000[,i][Chem000[,i]<0] <- NA
# }
#remove flagged data (turn to NA)
# flags A-D, F, I, J, K, L, M should be removed
# A = Sample suspect; B = standard curve/ reduction suspect; C = No sample taken; D = sample lost; F = Duplicate analyses in error; I = outside of data entry constraints; J = nonstandard routine followed; K = data suspect; L = data point and blind value differ by more than 15%; M = "More than t" (I don't know what that means).
#E = average of duplicate analysis, G = analyzed late, H = outside of standard range.

ChemTypes <- c("depth", "ph", "phair", "alk", "dic", "tic", "doc", "toc", "no3no2", "no2", "nh4", "totnf", "totnuf", "totpf", "totpuf", "drsif", "brsif", "brsiuf", "tpm", "no3no2_sloh", "nh4_sloh", "kjdl_n_sloh", "totpuf_sloh", "drp_sloh", "drsif_sloh")
# ChemFlags <- paste("flag", ChemTypes, sep="")
# for(i in 1:length(ChemFlags)){
# 	badchems <- which(ChemFlags[i]!="") #wow, this really failed! hahaha, oops. Plus, just check length(which(Chem000<0))
# 	Chem000[badchems,ChemTypes[i]] <- NA
# }

#using flag2NA instead of my super fail method (which basically removed the first length(ChemTypes) rows...) reduced the number of negative values from 481 to 416... hmm.


#average reps, stations
Chem00 <- aggregate(Chem000[,c("ph", "phair", "alk", "dic", "tic", "doc", "toc", "no3no2", "no2", "nh4", "totnf", "totnuf", "totpf", "totpuf", "drsif", "brsif", "brsiuf", "tpm", "no3no2_sloh", "nh4_sloh", "kjdl_n_sloh", "totpuf_sloh", "drp_sloh", "drsif_sloh")], by=Chem000[, rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)
Chem0 <- merge(Chem00, Zmix, by=c("lakeid", "year4", "daynum"), all.x=TRUE)
# Chem0 <- merge(x=Chem00, y=Phys[,c("lakeid", "year4", "daynum", "sampledate", "Zmix")], by=c("lakeid", "year4", "daynum", "sampledate"), all.x=TRUE) #there is one more NA in Zmix when merging with Phys than when merging w/ Zmix.  Don't know why.
Chem_Not_All_NA_Ind <- which(!apply(Chem0[,c("ph", "phair", "alk", "dic", "tic", "doc", "toc", "no3no2", "no2", "nh4", "totnf", "totnuf", "totpf", "totpuf", "drsif", "brsif", "brsiuf", "tpm", "no3no2_sloh", "nh4_sloh", "kjdl_n_sloh", "totpuf_sloh", "drp_sloh", "drsif_sloh")], MARGIN=1, FUN=allna))
Chem0 <- Chem0[Chem_Not_All_NA_Ind,]

Chem0_clean <- manClean(Chem0, ChemTypes)
Chem0 <- Chem0_clean[[1]]
deemedBad$Chem <- Chem0_clean[[2]]

Chem0_clean2 <- manClean(Chem0, ChemTypes)
Chem0 <- Chem0_clean2[[1]]
deemedBad$Chem2 <- Chem0_clean2[[2]]

Chem0_clean3 <- manClean(Chem0, ChemTypes)
Chem0 <- Chem0_clean3[[1]]
deemedBad$Chem3 <- Chem0_clean3[[2]]


#average epilimnion
Chem <- ddply(Chem0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
names(Chem) <- c("lakeid", "year4", "daynum", "sampledate", names(Chem0)[!is.element(names(Chem0), c("lakeid", "year4", "daynum", "sampledate"))])
names(Chem)[5] <- "depth_chem"
# class(Chem[,4]) <- "Date"

# ==========
# = Secchi =
# ==========
Secchi000 <- read.csv("Secchi_Nrtn.csv", colClasses=c("factor", "integer", "integer", "factor", "integer", "numeric", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", "numeric")) #WHY ONLY NORTHERN?
Secchi00 <- aggregate(Secchi000[,c("secnview", "waveht", "cloud", "ice")], by=Secchi000[,rev(c("lakeid", "year4", "daynum", "sampledate"))], FUN=mean, na.rm=TRUE)
Sec_notNA <- which(!is.na(Secchi00[,"secnview"]))
Secchi0 <- Secchi00[Sec_notNA, c("lakeid", "year4", "daynum", "sampledate", "secnview", "waveht", "cloud", "ice")]
names(Secchi0) <- c("lakeid", "year4", "daynum", "sampledate", "Secchi", "waveht_secchi", "cloud_secchi", "ice_secchi")
Secchi0_clean <- manClean(Secchi0, c("Secchi"))
Secchi <- Secchi0_clean[[1]]
deemedBad$Secchi <- Secchi0_clean[[2]]







# ====================
# = Light Extinction =
# ====================
LiExt000 <- read.csv("LightExt_Nrtn.csv", colClasses=c("factor", "integer", "integer", "factor", "numeric", "character", "character")) #WHY ONLY NORTHERN?
LightNoFlag <- which(LiExt000[, "lightext_flag"]=="")
LiExt <- LiExt000[LightNoFlag, c("lakeid", "year4", "daynum", "sampledate", "extcoef")]





# =======
# = Ice =
# =======
Ice_Nrtn000 <- read.csv("Ice_Nrtn.csv")
Ice_Nrtn00 <- Ice_Nrtn000[which(!is.na(Ice_Nrtn000[,"firstopen"]) & !is.na(Ice_Nrtn000[,"lastopen"])),]
OpenDuration_Nrtn <- as.integer(difftime(as.Date(Ice_Nrtn00[,"datelastopen"]), as.Date(Ice_Nrtn00[,"datefirstopen"])))
Ice_Nrtn0 <- data.frame(Ice_Nrtn00[,c("lakeid", "sta", "year")], "DaysOpen"=OpenDuration_Nrtn)
Ice_Nrtn <- aggregate(Ice_Nrtn0[,"DaysOpen"], by=Ice_Nrtn0[,rev(c("lakeid", "year"))], mean, na.rm=TRUE)
names(Ice_Nrtn) <- c("year4", "lakeid", "DaysOpen")
# Ice_Nrtn00[order(format.Date(as.Date(Ice_Nrtn00[,"datefirstopen"]), format="%j")),]

Ice_Srtn000 <- read.csv("Ice_Srtn.csv")
FixIceDates <- function(x){
	# x <- x[complete.cases(x),]
	on_NA <- is.na(x[,"iceon"])
	x_on_NoNA <- x[!on_NA, c("lakeid", "season", "iceon", "ice_on", "ice_duration", "year4")]
	off_NA <- is.na(x[,"iceoff"])
	x_off_NoNA <- x[!off_NA, c("lakeid", "season", "iceoff", "ice_off", "ice_duration", "year4")]
	
	take4 <- function(x) paste(x[1:4], collapse="")
	IceOn_Year <- as.integer(unlist(lapply(strsplit(as.character(x_on_NoNA[,"iceon"]), split=""), take4)))
	IceOff_Year <- as.integer(unlist(lapply(strsplit(as.character(x_off_NoNA[,"iceoff"]), split=""), take4)))
	
	ensure2_month <- function(x){sprintf("%02s", x[1])}
	ensure2_day <- function(x){sprintf("%02s", x[2])}
	IceOn_Month <- unlist(lapply(strsplit(as.character(x_on_NoNA[,"ice_on"]), split="/"), ensure2_month))
	IceOn_Day <- unlist(lapply(strsplit(as.character(x_on_NoNA[,"ice_on"]), split="/"), ensure2_day))
	IceOff_Month <- unlist(lapply(strsplit(as.character(x_off_NoNA[,"ice_off"]), split="/"), ensure2_month))
	IceOff_Day <- unlist(lapply(strsplit(as.character(x_off_NoNA[,"ice_off"]), split="/"), ensure2_day))
	
	Ice_On <- paste(IceOn_Year, IceOn_Month, IceOn_Day, sep="-")
	Ice_Off <- paste(IceOff_Year, IceOff_Month, IceOff_Day, sep="-")
	
	x_on_NoNA[,"ice_on"] <- Ice_On
	x_off_NoNA[,"ice_off"] <- Ice_Off
	x2 <- merge(x_on_NoNA, x_off_NoNA, all=TRUE)
	return(x2[,c("lakeid", "year4", "ice_on", "ice_off", "ice_duration")])
}
Ice_Srtn00 <- FixIceDates(Ice_Srtn000)


CalcDaysOpen <- function(x){
	# OffYears <- format.Date(as.Date(x[,"ice_off"]), "%Y")
	OnYears <- x[,"year4"]
	for(i in 1:length(OnYears)){
		OffRow <- which(x[,"year4"]==(OnYears[i]-1))
		if(length(OffRow)==0 || is.na(x[OffRow,"ice_off"]) || is.na(x[i,"ice_on"])){ #should maybe set to length(OffRow)!=1 || ...x
			x[i,"DaysOpen"] <- NA
		}else{
			x[i,"DaysOpen"] <- as.integer(difftime(as.Date(x[i, "ice_on"]), as.Date(x[OffRow,"ice_off"])))
		}
	}
	return(x)
}
Ice_Srtn0 <- ddply(.data=Ice_Srtn00, .variables="lakeid", .fun=CalcDaysOpen)
Ice_Srtn <- Ice_Srtn0[,c("year4", "lakeid", "DaysOpen")]
Ice <- rbind(Ice_Nrtn, Ice_Srtn)

# Ice_Srtn0[order(format.Date(as.Date(Ice_Srtn0[,"ice_off"]), format="%j")),]

# IceModel <- lm(DaysOpen~I(year4-1850)*lakeid, data=Ice)
# summary(IceModel)
# plot(Ice[,"year4"], Ice[,"DaysOpen"], col=Ice[,"lakeid"])





# ================
# = Boat Traffic =
# ================
BoTraf000 <- read.csv("BoatTraffic.csv")
BoTraf000[,"sampledate"] <- as.Date(BoTraf000[,"sampledate"])
# ContDates <- data.frame("sampledate"=seq(as.Date("1976-04-11"), as.Date("2010-10-31"), by=1))
# Spaced_BoTraf <- merge(BoTraf000, ContDates, all=TRUE)
BoTraf000[,"Sun_Week"] <- format.Date(BoTraf000[,"sampledate"], format="%U")
BoTraf000[,"Mon_Week"] <- format.Date(BoTraf000[,"sampledate"], format="%W")
BoTraf000[,"daynum"] <- as.integer(format.Date(BoTraf000[,"sampledate"], format="%j"))
BoTraf000[,"lakeid"] <- "ME_MO"
BoTraf000 <- BoTraf000[,c("lakeid", "year4", "Mon_Week", "daynum", "sampledate", "total_boats")]
BoTraf00 <- aggregate(BoTraf000[,c("sampledate", "daynum")], by=BoTraf000[,c("Mon_Week", "year4", "lakeid")], FUN=max, na.rm=TRUE) #the date at the end of the weekend
BoTraf0 <- aggregate(BoTraf000[,"total_boats"], by=BoTraf000[,c("Mon_Week", "year4", "lakeid")], FUN=sum, na.rm=TRUE)
names(BoTraf0) <- c("Mon_Week", "year4", "lakeid", "total_boats")
BoTraf <- merge(BoTraf0, BoTraf00, all=TRUE)
BoTraf <- BoTraf[order(BoTraf[,"sampledate"]),]
BoTraf[,"Mon_Week"] <- as.integer(BoTraf[,"Mon_Week"])

plot(BoTraf[,"sampledate"], BoTraf[,"total_boats"], type="o", pch=20)





# ===============================================
# = Boat Traffic (introduce NA's, weekly means)
# ===============================================

BoTraf2000 <- read.csv("BoatTraffic.csv")
BoTraf2000[,"sampledate"] <- as.Date(BoTraf2000[,"sampledate"])
ContDates <- data.frame("sampledate"=seq(as.Date("1976-04-11"), as.Date("2010-10-31"), by=1))
BoTraf2000 <- merge(BoTraf2000, ContDates, all=TRUE)
BoTraf2000[,"Sun_Week"] <- format.Date(BoTraf2000[,"sampledate"], format="%U")
BoTraf2000[,"Mon_Week"] <- format.Date(BoTraf2000[,"sampledate"], format="%W")
BoTraf2000[,"daynum"] <- as.integer(format.Date(BoTraf2000[,"sampledate"], format="%j"))
BoTraf2000[,"year4"] <- as.integer(format.Date(BoTraf2000[,"sampledate"], format="%Y"))
BoTraf2000[,"lakeid"] <- "ME_MO"
BoTraf2000 <- BoTraf2000[,c("lakeid", "year4", "Mon_Week", "daynum", "sampledate", "total_boats")]
BoTraf200 <- aggregate(BoTraf2000[,c("sampledate", "daynum")], by=BoTraf2000[,c("Mon_Week", "year4", "lakeid")], FUN=max, na.rm=TRUE) #the date at the end of the weekend
BoTraf20 <- aggregate(BoTraf2000[,"total_boats"], by=BoTraf2000[,c("Mon_Week", "year4", "lakeid")], FUN=mean, na.rm=TRUE)
names(BoTraf20) <- c("Mon_Week", "year4", "lakeid", "total_boats")
BoTraf2 <- merge(BoTraf20, BoTraf200, all=TRUE)
BoTraf2 <- BoTraf2[order(BoTraf2[,"sampledate"]),]
BoTraf2[,"Mon_Week"] <- as.integer(BoTraf2[,"Mon_Week"])
BoTraf2 <- subset(BoTraf2, !is.element(year4, c(1976, 2010)))

# dev.new(width=7, height=5)
# par(mar=c(4,4,0.5,0.5))
# plot(BoTraf2[,"sampledate"], BoTraf2[,"total_boats"], type="o", pch=20, xlab="Date", ylab="Num. Boats (daily averages for each week)")
# acf(BoTraf2[,"total_boats"])

BoTraf3 <- aggregate(BoTraf2[,"total_boats"], by=list(BoTraf2[,"year4"]), max, na.rm=TRUE)
names(BoTraf3) <- c("year", "total_boats")
# dev.new(width=7, height=5)
# par(mar=c(4,4,0.5,0.5))
# plot(BoTraf3[,"year"], BoTraf3[,"total_boats"], type="o", pch=20, xlab="Year", ylab="Num. Boats (yearly peak of daily averages)")

# Traffic <- BoTraf2[,"total_boats"]
# plot(BoTraf2[,"sampledate"], scale(BoTraf2[,"total_boats"]), type="o", pch=20)
# test <- arima(scale(BoTraf2[,"total_boats"]), order=c(1,0,0), seasonal=list(order=c(1,0,0), period=52))

# =============
# = Sun Spots =
# =============
SunSpots000 <- matrix(scan("DailySunspots_1818_to_2013.txt"), byrow=TRUE, ncol=5, dimnames=list(NULL, c("year4","month","dayOmonth","DecYear","SpotNum")))
SpotDates <- paste(SunSpots000[,"year4"], sprintf("%02s", SunSpots000[,"month"]), sprintf("%02s", SunSpots000[,"dayOmonth"]), sep="-")
SunSpots000 <- as.data.frame(SunSpots000)
SunSpots000[,"sampledate"] <- as.Date(SpotDates)
SunSpots000[,"daynum"] <- as.integer(format.Date(SunSpots000[,"sampledate"], format="%j"))

SunSpots00 <- SunSpots000[-c(1:(which(!is.na(SunSpots000[,"SpotNum"]))[1]-1)),]
SunSpots0_Dates <- aggregate(SunSpots00[,c("daynum", "DecYear", "sampledate")], by=SunSpots00[,c("year4", "month")], mean, na.rm=TRUE)
SunSpots0 <- aggregate(SunSpots00[,"SpotNum"], by=SunSpots00[,c("year4", "month")], mean, na.rm=TRUE)
names(SunSpots0) <- c("year4", "month", "SpotNum")
SunSpots <- merge(SunSpots0_Dates, SunSpots0, all=TRUE)
SunSpots <- SunSpots[order(SunSpots[,"sampledate"]),]

SpotTS <- ts(data=SunSpots[,"SpotNum"], start=c(SunSpots[1,"year4"], SunSpots[1,"month"]), frequency=12)
SpotStruc <- StructTS(SpotTS)
SmoothSpotTS <- ts(rowSums(tsSmooth(SpotStruc)))
tsp(SmoothSpotTS) <- tsp(SpotTS)
SpotTS[74] <- SmoothSpotTS[74]

SpotSpec <- spectrum(SpotTS[1:length(SpotTS)])
PeakFreq <- (SpotSpec$freq[which.max(SpotSpec$spec)])
SolarCycle <- 1/(PeakFreq*12)

SunSpots[,"SpotNum"] <- as.integer(SpotTS)
# arima(SpotTS, order=c(1,0,0), seasonal=list(order=c(1,0,0), period=12*SolarCycle))






# =============================================
# = Export data for Organization into extrema =
# =============================================
Met[,"sampledate"] <- as.Date(Met[,"sampledate"])
LakLev[,"sampledate"] <- as.Date(LakLev[,"sampledate"])
Zmix[,"sampledate"] <- as.Date(format.Date(as.POSIXct(paste(Zmix[,"year4"], sprintf("%03d",Zmix[,"daynum"]), sep="-"), format="%Y-%j"), format="%Y-%m-%d"))
LiExt[,"sampledate"] <- as.Date(LiExt[,"sampledate"])
Secchi[,"sampledate"] <- as.Date(Secchi[,"sampledate"])
Phys[,"sampledate"] <- as.Date(Phys[,"sampledate"])
names(Phys)[5] <- "depth_phys"
Ions[,"sampledate"] <- as.Date(Ions[,"sampledate"])
Chem[,"sampledate"] <- as.Date(Chem[,"sampledate"])
Chl[,"sampledate"] <- as.Date(Chl[,"sampledate"])
Zoop[,"sampledate"] <- as.Date(Zoop[,"sampledate"])


Data_X <- list("SunSpots"=SunSpots, "Met"=Met, "Ice"=Ice, "LakLev"=LakLev, "Zmix"=Zmix, "LiExt"=LiExt, "Secchi"=Secchi, "Phys"=Phys, "Ions"=Ions, "Chem"=Chem, "Chl"=Chl, "Zoop"=Zoop, "Fish"=Fish_GearSpec, "BoTraf"=BoTraf, "Fish_ByGear"=Fish_ByGear, "Fish_BySpec"=Fish_BySpec)

save(Data_X, file="OrganizedFatData_Read_Fat_Data_v2.RData")

