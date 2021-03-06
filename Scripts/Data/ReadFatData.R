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
# library("reshape")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/Data_Functions.R")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Fat_dGEV.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/CalcZoopBiomass.R")


# deemedBad <- list()

# ======================================
# = Define subsets of variables to use =
# ======================================
AllPhys <- c( "wtemp", "o2", "o2sat", "frlight", "Zmix")
AllIons <- c("cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")
Zoop_Cats <- c("density", "tot_zoop_mass")
Met_Cats <- c("min_air_temp", "max_air_temp", "range_air_temp", "precip_mm", "snow_cm")
AllChem <- c("ph","alk","dic", "tic", "doc", "toc", "no3no2","nh4", "totnuf", "totpuf", "drsif", "brsiuf", "tpm")


# Need to remove spuries "species" from Zoop data. (e.g., UNID)
#The Met data doesn't include wind speeds
# For the analysis of this data, I'm starting to think that I might want to use annual maxima as a solid first-cut.  This will simplify issues of irregular time series, making it easier to computer return-rate and return-levels.


# ======================
# = Physical Limnology =
# ======================
# Read in raw data and deal with flags
Phys000 <- read.csv("Phys_NrtnSrtn.csv")
Phys00 <- subset(Phys000, sta==1 & rep==1)
Phys00 <- flag2NA(Phys00)

# Manually clean, save bad indices
# Phys00_clean <- manClean(Phys00, c("wtemp", "o2", "o2sat", "frlight"))
# Phys00 <- Phys00_clean[[1]]
# deemedBad$Phys <- Phys00_clean[[2]]
# deemedBad.phys <- which(is.na(Phys00_clean[[1]]))
# write.table(deemedBad.ions, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.phys.txt", sep="\t", row.names=FALSE)

# Read in Bad indices from manual clean
deemedBad.phys0 <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.phys.txt", sep="\t", header=TRUE)[,]
deemedBad.phys <- arrayInd(deemedBad.phys0, dim(Phys00), dimnames(Phys00))
Phys00[deemedBad.phys] <- NA

# Get Zmix
Zmix <- ddply(.data=Phys00, .variables=c("lakeid", "year4", "daynum"), .fun=zmix)
names(Zmix) <- c("lakeid", "year4", "daynum", "Zmix")
Zmix[,2] <- as.numeric(as.character(Zmix[,2]))
Zmix[,3] <- as.numeric(as.character(Zmix[,3]))
Zmix[,4] <- as.numeric(as.character(Zmix[,4]))

# Merge Zmix with other Physical data
Phys0 <- merge(Phys00, Zmix, all=TRUE)
zd <- Phys0[,"Zmix"] - Phys0[,"depth"] # the difference between the depth of the measurements and zmix (positive means the obs is above zmix)
Phys0[,"InEpi"] <- as.numeric(zd>0 & !is.na(zd)) # if the difference is positive, mark as in epi
Phys0[,"Mixing"] <- as.numeric(is.na(zd)) # If the difference between zmix and depth of obs is NA, assume mixing (no zmix)
Phys0[,"InMetaHypo"] <- as.numeric(Phys0[,"InEpi"]+Phys0[,"Mixing"] == 0) # if the observation isn't in epi and the lake isn't mixing, define as meta or hypo
Phys0 <- subset(Phys0, InMetaHypo==0, select=c("lakeid", "year4", "daynum", "sampledate", "depth", "wtemp", "o2", "o2sat", "deck", "light", "frlight", "Zmix", "InEpi", "Mixing"))

# Calculate epilimnetic means
Phys1 <- ddply(.data=Phys0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean)

# Reshape into long format w/ correct names
Phys1 <- Phys1[,c("lakeid", "year4", "daynum", AllPhys)]
Phys <- cbind("taxID"="Phys", "taxLvl"="Community", reshape2(Phys1, varying=AllPhys, times=AllPhys, v.names="Data", direction="long"))
names(Phys) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# =======================
# = Chlorophyll Data =
# =======================
# Read in the data from the northern lakes, deal with flags
Chl000n <- read.csv("Chla_Nrtn.csv")
ChlorFlags <- strsplit("A AJ AL B BG BK BL D G GJ GK GL H I J JK JKL JL K L KL KO L O LB LG LK LO O OJ", split=" ")[[1]]
BadChlorN <- which(is.element(Chl000n[,"flagchlor"], ChlorFlags) | Chl000n[,"chlor"]<0)
PhaeoFlags <- strsplit("A AJ AK AL B BG BK BL D G GB GK GL H I J JK JKL JL K K L KL KLO KO L L K LB LG LK LO O OJ OK", split=" ")[[1]]
BadPhaeoN <- which(is.element(Chl000n[,"flagphaeo"], PhaeoFlags) | Chl000n[,"phaeo"]<0)
Chl000n[union(BadChlorN, BadPhaeoN),"chlor"] <- NA
Chl000n[union(BadChlorN, BadPhaeoN), "phaeo"] <- NA

# Chl000n <- manClean(Chl000n, c("chlor", "phaeo"))
#I'm not going to bother using the flag2NA function here, b/c I implemented a similar, more manual, method in the first version (before the flag2NA function was created)

# Read in data from the southern lakes, deal with flags
Chl000s <- read.csv("Chla_Srtn.csv")
Chl000s[,"daynum"] <- as.numeric(as.character(format.Date(Chl000s[,"sampledate"], format="%j")))
nummean <- function(x){mean(as.numeric(x))}
Chl000s[,"depth"] <- unlist(lapply(strsplit(as.character(Chl000s[,"depth_range_m"]), split="-"), FUN=nummean))
useUncorr <- Chl000s[,"flag_fluor"]%in%c("C", "CE", "EC")
Chl000s[useUncorr,"correct_chl_fluor"] <- Chl000s[useUncorr,"uncorrect_chl_fluor"] 
BadChlorS <- which(is.element(Chl000s[,"flag_fluor"], c("E", "CE", "EC")) | Chl000s[,"correct_chl_fluor"]<0)
Chl000s[BadChlorS,"correct_chl_fluor"] <- NA
Chl000s <- Chl000s[,c("lakeid", "year4", "sampledate", "daynum", "depth", "rep", "correct_chl_fluor")]
names(Chl000s) <- c("lakeid", "year4", "sampledate", "daynum", "depth", "rep", "chlor")

# Chl000s <- manClean(Chl000s, "chlor")

# Merge chlorophyll from northern and southern lakes
Chl000 <- merge(Chl000n, Chl000s, all=TRUE)
Chl00 <- merge(Phys1[,c("lakeid", "year4", "daynum", "Zmix")], Chl000, all=TRUE)
Chl00 <- Chl00[, c("lakeid", "year4", "daynum", "sampledate", "Zmix", "rep", "depth", "chlor", "sta", "phaeo")]
Chl00 <- aggregate(Chl00[,c("Zmix", "chlor", "phaeo")], by=Chl00[,rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)
Chl0 <- ddply(.data=Chl00, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
Chl1 <- Chl0[,c("lakeid", "year4", "daynum", "chlor")]
names(Chl1) <- c("location", "year4", "daynum", "chlor")

# Reshape and rename data
Chl <- cbind("taxID"="Chl", "taxLvl"="Community", reshape2(Chl1, varying="chlor", times="chlor", v.name="Data", direction="long"))
names(Chl) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# ====================
# = Zooplankton Data =
# ====================
#not using flag2NA b/c these files don't have flag columns

	# =====================
	# = Southern Zoops #1 =
	# =====================
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

	# =====================
	# = Southern Zoops #2 =
	# =====================
Zoop000s2 <- read.csv("ZoopDens_New_Srtn.csv")#density recorded as individuals/m^2
Zoop00s2 <- Zoop000s2
Zoop00s2[,"density"] <- (Zoop000s2[,"density"]/Zoop000s2[,"towdepth"])*0.001 #convert zooplankton density to individuals/L
names(Zoop00s2) <- c("lakeid", "year4", "sampledate", "station", "towdepth", "species_code", "taxon", "density", "indiv_measured", "avg_length")

	# ========================
	# = Merge southern Zoops =
	# ========================
Zoop00s <- merge(Zoop00s1, Zoop00s2, all=TRUE)
Zoop00s[,"daynum"] <- as.numeric(as.character(format.Date(Zoop00s[,"sampledate"], format="%j")))
for(i in 1:length(unique(Zoop00s[,"taxon"]))){
	TaxonRows <- which(is.element(Zoop00s[,"taxon"], unique(Zoop00s[,"taxon"])[i]))
	Zoop00s[TaxonRows, "avg_zoop_mass"] <- ZoopMass(Zoop00s[TaxonRows,])
}
Zoop00s[, "tot_zoop_mass"] <- Zoop00s[,"density"] * Zoop00s[,"avg_zoop_mass"]
Zoop0s <- Zoop00s[,c("lakeid", "year4", "daynum", "sampledate", "taxon", "density", "avg_length", "avg_zoop_mass", "tot_zoop_mass")]

	# ==================
	# = Northern Zoops =
	# ==================
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

	# =====================================
	# = Merge northern and southern zoops =
	# =====================================
Zoop0 <- merge(Zoop0n, Zoop0s, all=TRUE) #this contains zooplankton information broken down by taxon
names(Zoop0)[names(Zoop0)=="taxon"] <- "spname"

	# ================================
	# = Save unique zooplankton taxa =
	# ================================
# zTax <- data.frame("spname"=as.character(unique(Zoop0[,"spname"])))
# write.table(zTax, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/uniqueZoopTaxa.txt", sep="\t", row.names=FALSE)

	# ==============================
	# = Sum zooplankton by 'taxon' =
	# ==============================
# Zoop_Cats <- c("density", "avg_length", "avg_zoop_mass", "tot_zoop_mass")
Zoop_Factrs <-  c("year4", "daynum","sampledate", "spname", "lakeid")
Zoop_ByTax0 <- aggregate(Zoop0[,Zoop_Cats], by=Zoop0[,Zoop_Factrs], FUN=Zoop_Sum)


	# =========================================
	# = Read in zoop taxonomic classification =
	# =========================================
zTax.id <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatZoopTaxa.txt", sep="\t", header=TRUE, colClasses="character")
zTax.id[zTax.id==""] <- NA

subPhy <- !is.na(zTax.id[,"subphylum"])
zTax.id[subPhy,"phylum"] <- zTax.id[subPhy,"subphylum"]

subClass <- !is.na(zTax.id[,"subclass"])
zTax.id[subClass,"class"] <- zTax.id[subClass,"subclass"]

subSpec <- !is.na(zTax.id[,"species"])
zTax.id[subSpec,"species"] <- paste(substring(zTax.id[subSpec,"genus"], 1, 1), zTax.id[subSpec,"species"], sep=".")


	# ===============================================
	# = Combine taxonomic classification w/ metrics =
	# ===============================================
Zoop_ByTax <- cbind(Zoop_ByTax0, "Phylum"=NA, "Class"=NA, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
for(i in 1:nrow(zTax.id)){
	ind <- Zoop_ByTax[,"spname"] == zTax.id[i,"spname"]
	Zoop_ByTax[ind, c("Phylum", "Class", "Order","Family","Genus","Species")] <- zTax.id[i,c("phylum", "class", "order", "family", "genus", "species")]
}

	# ====================================================
	# = Create Zoop time series for each taxonomic level =
	# ====================================================
zf <- Zoop_Factrs[Zoop_Factrs!="spname"]
zcf <- c(zf, Zoop_Cats)
	
# Zoop Species	
zoop.species0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Species"]), c(zcf, "Species")]
zoop.species0s <- aggregate(zoop.species0[,Zoop_Cats], by=zoop.species0[,c("Species", zf)], Zoop_Sum) # don't need this step for species level
names(zoop.species0s)[names(zoop.species0s)=="Species"] <- "taxID"
zoop.species <- cbind(zoop.species0s, "taxLvl"="Species")

# Zoop Genus	
zoop.genus0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Genus"]), c(zcf, "Genus")]
zoop.genus0s <- aggregate(zoop.genus0[,Zoop_Cats], by=zoop.genus0[,c("Genus", zf)], Zoop_Sum)
names(zoop.genus0s)[names(zoop.genus0s)=="Genus"] <- "taxID"
zoop.genus <- cbind(zoop.genus0s, "taxLvl"="Genus")

# Zoop Family	
zoop.family0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Family"]), c(zcf, "Family")]
zoop.family0s <- aggregate(zoop.family0[,Zoop_Cats], by=zoop.family0[,c("Family", zf)], Zoop_Sum)
names(zoop.family0s)[names(zoop.family0s)=="Family"] <- "taxID"
zoop.family <- cbind(zoop.family0s, "taxLvl"="Family")

# Zoop Order	
zoop.order0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Order"]), c(zcf, "Order")]
zoop.order0s <- aggregate(zoop.order0[,Zoop_Cats], by=zoop.order0[,c("Order", zf)], Zoop_Sum)
names(zoop.order0s)[names(zoop.order0s)=="Order"] <- "taxID"
zoop.order <- cbind(zoop.order0s, "taxLvl"="Order")

# Zoop Class	
zoop.class0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Class"]), c(zcf, "Class")]
zoop.class0s <- aggregate(zoop.class0[,Zoop_Cats], by=zoop.class0[,c("Class", zf)], Zoop_Sum)
names(zoop.class0s)[names(zoop.class0s)=="Class"] <- "taxID"
zoop.class <- cbind(zoop.class0s, "taxLvl"="Class")

# Zoop Phylum	
zoop.phylum0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"Phylum"]), c(zcf, "Phylum")]
zoop.phylum0s <- aggregate(zoop.phylum0[,Zoop_Cats], by=zoop.phylum0[,c("Phylum", zf)], Zoop_Sum)
names(zoop.phylum0s)[names(zoop.phylum0s)=="Phylum"] <- "taxID"
zoop.phylum <- cbind(zoop.phylum0s, "taxLvl"="Phylum")

# Zoop All	
zoop.all0 <- Zoop_ByTax[!is.na(Zoop_ByTax[,"spname"]), c(zcf, "spname")]
zoop.all0s <- aggregate(zoop.all0[,Zoop_Cats], by=zoop.all0[,zf], Zoop_Sum)
zoop.all <- cbind("taxID"="Zoop", zoop.all0s, "taxLvl"="Community")


	# ==================================================================
	# = Combine zooplankton time series (were separated by tax. level) =
	# ==================================================================
Zoop1 <- rbind(zoop.all, zoop.phylum, zoop.class, zoop.order, zoop.family, zoop.genus, zoop.species)
Zoop <- reshape2(Zoop1, varying=Zoop_Cats, times=Zoop_Cats, v.name="Data", direction="long")
names(Zoop) <- c("taxID", "year4", "daynum", "sampledate", "location", "taxLvl", "variable", "Data")
Zoop <- Zoop[,names(Zoop)[names(Zoop)!="sampledate"]]



# ========
# = Ions =
# ========
# Read in ion data (both N and S), reshape
Ions000 <- flag2NA(read.csv("Ions_NrtnSrtn.csv")) #why are there so many rows with no data?
Ions00 <- merge(Ions000, Zmix, all.x=TRUE)

#I need to take the average for the rep's and stations first
Ions_Not_All_NA_Ind <- which(!apply(Ions00[,c("cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")], MARGIN=1, FUN=function(x){all(is.na(x))}))
Ions00 <- Ions00[Ions_Not_All_NA_Ind,]
Ions0 <- aggregate(Ions00[,c("Zmix", "cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")], by=Ions00[,rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)

# Manually clean, save bad indices
# Ions0_clean <- manClean(Ions0, IonTypes)
# Ions0 <- Ions0_clean[[1]]
# deemedBad$Ions <- Ions0_clean[[2]]
# deemedBad.ions <- which(is.na(Ions0_clean[[1]]))
# write.table(deemedBad.ions, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.ions.txt", sep="\t", row.names=FALSE)

# Read in bad indices defined by manual clean
deemedBad.ions0 <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.ions.txt", sep="\t", header=TRUE)[,]
deemedBad.ions <- arrayInd(deemedBad.ions0, dim(Ions0), dimnames(Ions0))
Ions0[deemedBad.ions] <- NA

Ions1 <- ddply(.data=Ions0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
names(Ions1) <- c("location", "year4", "daynum", "sampledate", "depth_ions", "zmix", "cl", "so4", "ca", "mg", "na", "k", "fe", "mn", "cond")
Ions1 <- Ions1[,c("location","year4","daynum", AllIons)]

Ions <- cbind("taxID"="Ions", "taxLvl"="Community", reshape2(Ions1, varying=AllIons, times=AllIons, v.name="Data", direction="long"))
names(Ions) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# ========
# = Fish =
# ========
	# ==================
	# = Fish Abundance =
	# ==================
# Read in fish abundance data
Fish000_abun <- read.csv("Fish_Abun_NrtnSrtn.csv")
Fish000_abun[,"gearid"] <- gsub("[0123456789]", "", Fish000_abun[,"gearid"]) #hell yeah, i learned a new function. (well, how to use reular expressions a bit better)
Fish000_abun[,"gearid"] <- gsub("(?<=FYKNE)[LD]", "T", Fish000_abun[,"gearid"], perl=TRUE)
Fish000_abun[,"spname"] <- gsub(" ", "", Fish000_abun[,"spname"]) #remove white space from species names
BsSpecies <- c("GUPPY", "VIRILIS", "RUSTICUS", "PROPINQUUS", "UNIDENTIFIED", "CRAYFISH", "LARVALFISH")
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
# Read in Fish size data
Fish000_size <- read.csv("Fish_LenWei_NrtnSrtn.csv")
Fish000_size[,"spname"] <- gsub(" ", "", Fish000_size[,"spname"]) 
Fish000_size[,"gearid"] <- gsub("[0123456789]", "", Fish000_size[,"gearid"])
Fish000_size[,"gearid"] <- gsub("(?<=FYKNE)[LD]", "T", Fish000_size[,"gearid"], perl=TRUE)
Fish00_size <- subset(Fish000_size, !is.element(spname, BsSpecies))
FixPike <- which(is.element(Fish00_size[,"spname"], "SILVERPIKE"))
Fish00_size[FixPike, "spname"] <- "NORTHERNPIKE"

# Summarize fish size data
Fish0_size <- ddply(Fish00_size, .variables=c("lakeid", "year4", "spname"), .fun=SizSumry)
TotFish0_size <- ddply(Fish00_size, .variables=c("lakeid", "year4", "spname", "gearid"), .fun=SizSumry)

	# ================
	# = Combine Fish =
	# ================
# Fish <- merge(Fish0_abun, Fish0_size, all=TRUE, by=c("lakeid", "year4", "spname"))
# Fish[,"SumWei"] <- Fish[,"Nwei"]*Fish[,"mean_Wei"]
# Fish[,"cpue3_WeiEff"] <- Fish[,"SumWei"]/Fish[,"effort"]

TotFish0 <- merge(Fish00_abun, TotFish0_size, all=TRUE, by=c("lakeid", "year4", "spname", "gearid"))
TotFish0[,"SumWei"] <- TotFish0[,"Nwei"]*TotFish0[,"mean_Wei"]
TotFish0[,"cpue3_WeiEff"] <- TotFish0[,"SumWei"]/TotFish0[,"effort"]

	# ==================================
	# = Write table with unique spname =
	# ==================================
# fTax <- data.frame("spname"=as.character(unique(TotFish0[,"spname"])))


	# ===========================================================
	# = Look for gears that were consistently used in each lake =
	# ===========================================================
# I can combine data from the different gear types, but if a gear type wasn't used in a partiular year, 
# then the sum catch per effort might be skewed. Taking the average could have the same effect because
# the different methods could have different efficiencies.
# My solution is to figure out which methods were used in "most" years in a given lake,
# and only take the sum catch for those methods.
# When one of the methods was not used, that year cannot be analyzed. 
# HOWEVER, if a particular taxon is only ever seen in Method B, then 
# in a year when Method A is missing, the time series for this particular taxon is still unbroken.


gly <- TotFish0[,c("lakeid","year4", "gearid")] # gear lake year
uLake <- unique(TotFish0[,"lakeid"]) # unique lakes
gl.obs <- paste(TotFish0[,"gearid"], TotFish0[,"lakeid"], sep="") # observed combinations of gear-lake-year (non-unique)
gl.valid <- c() # to store valid gear-lake combinations (gear observed for > 96% of years in that lake)
for(i in 1:length(uLake)){
	gearYear <- table(gly[gly[,1]==uLake[i],2], gly[gly[,1]==uLake[i],3]) # count up the number of observations of each gear type in each year for this lake
	gy <- gearYear[!apply(gearYear, 1, function(x)all(x==0 | !is.finite(x))),] # for this lake, which years didn't have any sampling (no gears observed)
	propPres <- apply(gy, 2, function(x){sum(x>0)/length(x)}) # discounting years w/ no sampling, what proportion of years was each gear type used in this lake
	validNames <- names(propPres)[propPres>0.96] # if a gear type was used in > 96% of years in this lake, it is a valid gear for this lake
	gl.valid <- c(gl.valid, paste(validNames, as.character(uLake[i]), sep="")) # store & accumulate the valid gear-lake combinations
}
TotFish <- TotFish0[gl.obs%in%gl.valid,] # remove any rows that used a non-valid gear type for that lake



	# ===========================================
	# = Aggregate spgyl dups, subset to metrics =
	# ===========================================
# FishCats1 <- c("total_caught", "cpue1_Sum", "cpue3_WeiEff")
FishCats1 <- c("cpue1_Sum", "cpue3_WeiEff")
FishFactrs <- c("spname", "gearid", "year4", "lakeid")
Fish_ByGearSpec <- aggregate(TotFish[, FishCats1], by=TotFish[,FishFactrs], sum, na.rm=TRUE)


	# ========================================================================
	# = Determine lake-years when species weren't caught due to missing gear =
	# ========================================================================
# This part was tough
# m1 shows all possible year-lake-gear-spname combinations, given past obs of year-lake and species-gear-lake combintations
f1 <- Fish_ByGearSpec[!duplicated(paste(Fish_ByGearSpec[,"year4"], Fish_ByGearSpec[,"lakeid"])), c("year4","lakeid")]
f4 <- Fish_ByGearSpec[!duplicated(paste(Fish_ByGearSpec[,"spname"], Fish_ByGearSpec[,"gearid"], Fish_ByGearSpec[,"lakeid"])), c("spname","gearid", "lakeid")]
m1 <- merge(f1, f4, all=TRUE)
m1.lyg <- apply(m1[,c("lakeid", "year4", "gearid")], 1, paste, collapse=" ")

# f3 shows past combinations of year-gear-lake
f3 <- Fish_ByGearSpec[!duplicated(paste(Fish_ByGearSpec[,"year4"], Fish_ByGearSpec[,"gearid"], Fish_ByGearSpec[,"lakeid"])), c("year4","gearid", "lakeid")]
row.names(f3) <- NULL
f3.lyg <- apply(f3[,c("lakeid", "year4", "gearid")], 1, paste, collapse=" ")

# given the possible lake-year-gear combos (m1.lyg), 
# and the observed lake-year-gear combos (f3.lyg),
# determine which lake-year-species (or lake-year-gear-species) combinations could be missing
lygs.miss <- m1[!m1.lyg%in%f3.lyg,]


	# ===============================================================
	# = Add NA's for species who weren't caught b/c of missing gear =
	# ===============================================================
Fish_ByGearSpec <- merge(Fish_ByGearSpec, lygs.miss, all=TRUE)


	# ======================================
	# = Sum Fish metrics across gear types =
	# ======================================
Fish_BySpec0 <- aggregate(Fish_ByGearSpec[, FishCats1], by=Fish_ByGearSpec[,c("spname", "year4", "lakeid")], sum, na.rm=FALSE)


	# ==============================================
	# = Add taxonomic classification to data frame =
	# ==============================================
fTax.id <- read.table(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatFishTaxa.txt", sep="\t", header=TRUE, colClasses="character")
fTax.id[fTax.id==""] <- NA

fsubSpec <- !is.na(fTax.id[,"species"])
fTax.id[fsubSpec,"species"] <- paste(substring(fTax.id[fsubSpec,"genus"], 1, 1), fTax.id[fsubSpec,"species"], sep=".")

Fish_BySpec <- cbind(Fish_BySpec0, "Order"=NA, "Family"=NA, "Genus"=NA, "Species"=NA)
for(i in 1:nrow(fTax.id)){
	ind <- Fish_BySpec[,"spname"] == fTax.id[i,"spname"]
	Fish_BySpec[ind, c("Order","Family","Genus","Species")] <- fTax.id[i,c("order", "family", "genus", "species")]
}



	# =================================================
	# = Create time series for each tax level of fish =
	# =================================================
ff <- FishFactrs[!FishFactrs%in%c("spname","gearid")]
fcf <- c(FishCats1, ff) # fish categories and factors

# Fish Species	
fish.species0 <- Fish_BySpec[!is.na(Fish_BySpec[,"Species"]), c(fcf, "Species")]
fish.species0s <- aggregate(fish.species0[,FishCats1], by=fish.species0[,c("Species", ff)], sum) # don't need this step for species level
names(fish.species0s)[names(fish.species0s)=="Species"] <- "taxID"
fish.species <- cbind(fish.species0s, "taxLvl"="Species")

# Fish Genus	
fish.genus0 <- Fish_BySpec[!is.na(Fish_BySpec[,"Genus"]), c(fcf, "Genus")]
fish.genus0s <- aggregate(fish.genus0[,FishCats1], by=fish.genus0[,c("Genus", ff)], sum)
names(fish.genus0s)[names(fish.genus0s)=="Genus"] <- "taxID"
fish.genus <- cbind(fish.genus0s, "taxLvl"="Genus")

# Fish Family	
fish.family0 <- Fish_BySpec[!is.na(Fish_BySpec[,"Family"]), c(fcf, "Family")]
fish.family0s <- aggregate(fish.family0[,FishCats1], by=fish.family0[,c("Family", ff)], sum)
names(fish.family0s)[names(fish.family0s)=="Family"] <- "taxID"
fish.family <- cbind(fish.family0s, "taxLvl"="Family")

# Fish Order	
fish.order0 <- Fish_BySpec[!is.na(Fish_BySpec[,"Order"]), c(fcf, "Order")]
fish.order0s <- aggregate(fish.order0[,FishCats1], by=fish.order0[,c("Order", ff)], sum)
names(fish.order0s)[names(fish.order0s)=="Order"] <- "taxID"
fish.order <- cbind(fish.order0s, "taxLvl"="Order")


# Fish All	
fish.all0 <- Fish_BySpec[!is.na(Fish_BySpec[,"spname"]), c(fcf, "spname")]
fish.all0s <- aggregate(fish.all0[,FishCats1], by=fish.all0[,ff], sum)
fish.all <- cbind("taxID"="Fish", fish.all0s, "taxLvl"="Community")

	# =============================================
	# = Combine fish time series into long format =
	# =============================================
Fish0 <- rbind(fish.all, fish.order, fish.family, fish.genus, fish.species)
Fish <- cbind("daynum"=1, reshape2(Fish0, varying=FishCats1, times=FishCats1, v.name="Data", direction="long"))
names(Fish) <- c("daynum", "taxID", "year4", "location", "taxLvl", "variable", "Data")





# ==================
# = Meteorological =
# ==================
# http://lter.limnology.wisc.edu/data/filter/11001
# above link for woodruff airport data from 1989 on, with wind and par etc.
# Read in Madison met
Met000_Mad <- read.csv("MetDat_Madison.csv")
Met00_Mad <- Met000_Mad[,c("year4", "daynum", "sampledate", "max_air_temp_adjusted", "min_air_temp_adjusted", "ave_air_temp_adjusted", "range_air_temp_adjusted", "precip_raw_mm", "snow_raw_cm", "snow_depth_cm", "data_status")]
names(Met00_Mad) <- c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "ave_air_temp", "range_air_temp", "precip_mm", "snow_cm", "snow_depth_cm", "data_status")
Met00_Mad[,"location"] <- "Madison"

# Read in Minocqua met
Met000_Min <- read.csv("MetDat_Minocqua.csv") # from Minocqua Dam
Met000_Min[,"daynum"] <- as.integer(as.character(format.Date(as.Date(Met000_Min[,"sampledate"]), format="%j")))
Met000_Min[,"range_air_temp"] <- Met000_Min[,"max_air_temp"] - Met000_Min[,"min_air_temp"]
Met00_Min <- Met000_Min[,c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "range_air_temp", "precip", "snow", "snow_depth", "data_status")]
names(Met00_Min) <- c("year4", "daynum", "sampledate", "max_air_temp", "min_air_temp", "range_air_temp", "precip_mm", "snow_cm", "snow_depth_cm", "data_status")
Met00_Min[,"location"] <- "Minocqua"
Met00_Min[,"snow_cm"] <- Met00_Min[,"snow_cm"]/10 #the minocqua snowfall was originally in mm

# Merge Madison and Minocqua
Met_Factrs <- c("location", "year4", "daynum")
Met0 <- merge(Met00_Mad, Met00_Min, all=TRUE)
Met1 <- Met0[,c(Met_Factrs, Met_Cats)]

Met <- cbind("taxID"="Meteorological", "taxLvl"="Community", reshape2(Met1, varying=Met_Cats, times=Met_Cats, v.name="Data", direction="long"))
names(Met) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# ========
# = Wind =
# ========
	# ================
	# = Madison Wind =
	# ================
# Using the "corrected" values (sensor location changed, etc)
Wind00_Mad <- read.table("wind_Madison.txt", sep="\t", header=TRUE, colClasses=c("integer", "integer", "double"))
Wind0_Mad <- Wind00_Mad
names(Wind0_Mad) <- c("year4", "daynum", "Data")
Wind0_Mad["location"] <- "Madison"
Wind0_Mad["taxID"] <- "Wind"
Wind0_Mad["taxLvl"] <- "Community"
Wind0_Mad["variable"] <- "wind"

	# =================
	# = Minocqua Wind =
	# =================
# From Woodruff airport
Wind000_Min <- read.table("wind_Minocqua.txt", sep="\t", header=TRUE, colClasses=c("integer", "integer", "double", "double", "character", "character"))
Wind00_Min <- flag2NA(Wind000_Min)[,c(1,2,4)] # remove flagged values and subset to year, day, daily average wind speed (no peak speed for madison, so being consistent)
Wind0_Min <- Wind00_Min
names(Wind0_Min) <- c("year4", "daynum", "Data")
Wind0_Min["location"] <- "Minocqua"
Wind0_Min["taxID"] <- "Wind"
Wind0_Min["taxLvl"] <- "Community"
Wind0_Min["variable"] <- "wind"

# Combine minocqua and madison
Wind <- rbind(Wind0_Mad, Wind0_Min)




# ==============
# = Lake Level =
# ==============
LakLev000 <- read.csv("LakeLevel_Nrtn.csv") # WHY ONLY NORTHERN?; level is in "meters above sea level"
LakLev1 <- LakLev000[,c("lakeid", "year4", "daynum", "llevel_elevation")]
names(LakLev1) <- c("location", "year4", "daynum", "LakeLevel")
LakLev <- cbind("taxID"="LakeLevel","taxLvl"="Community", reshape2(LakLev1, varying="LakeLevel", times="LakeLevel", v.name="Data", direction="long"))
names(LakLev) <- c("taxID", "taxLvl","location","year4","daynum","variable","Data")


# ======================
# = Chemical Limnology =
# ======================
Chem000 <- flag2NA(read.csv("Chem_NrtnSrtn.csv"))
# flags A-D, F, I, J, K, L, M should be removed
# A = Sample suspect; B = standard curve/ reduction suspect; C = No sample taken; D = sample lost; F = Duplicate analyses in error; I = outside of data entry constraints; J = nonstandard routine followed; K = data suspect; L = data point and blind value differ by more than 15%; M = "More than t" (I don't know what that means).
#E = average of duplicate analysis, G = analyzed late, H = outside of standard range.

# ChemTypes <- c("depth", "ph", "phair", "alk", "dic", "tic", "doc", "toc", "no3no2", "no2", "nh4", "totnf", "totnuf", "totpf", "totpuf", "drsif", "brsif", "brsiuf", "tpm", "no3no2_sloh", "nh4_sloh", "kjdl_n_sloh", "totpuf_sloh", "drp_sloh", "drsif_sloh")
# 
# AllChem <- c("ph","alk","dic", "tic", "doc", "toc", "no3no2","nh4", "totnuf", "totpuf", "drsif", "brsiuf", "tpm")
#average reps, stations
Chem00 <- aggregate(Chem000[,AllChem], by=Chem000[, rev(c("lakeid", "year4", "daynum", "sampledate", "depth"))], FUN=mean, na.rm=TRUE)
Chem0 <- merge(Chem00, Zmix, by=c("lakeid", "year4", "daynum"), all.x=TRUE)
# Chem0 <- merge(x=Chem00, y=Phys[,c("lakeid", "year4", "daynum", "sampledate", "Zmix")], by=c("lakeid", "year4", "daynum", "sampledate"), all.x=TRUE) #there is one more NA in Zmix when merging with Phys than when merging w/ Zmix.  Don't know why.
Chem_Not_All_NA_Ind <- which(!apply(Chem0[,AllChem], MARGIN=1, FUN=function(x)all(is.na(x))))
Chem0 <- Chem0[Chem_Not_All_NA_Ind,]

# Chem0_clean <- manClean(Chem0, ChemTypes)
# Chem0 <- Chem0_clean[[1]]
# # which(is.na(Chem0_clean[[1]]), arr.ind=TRUE)
# deemedBad$Chem <- Chem0_clean[[2]]
# 
# Chem0_clean2 <- manClean(Chem0, ChemTypes)
# Chem0 <- Chem0_clean2[[1]]
# # which(is.na(Chem0_clean2[[1]]), arr.ind=TRUE)
# deemedBad$Chem2 <- Chem0_clean2[[2]]
# 
# Chem0_clean3 <- manClean(Chem0, ChemTypes)
# Chem0 <- Chem0_clean3[[1]]
# # which(is.na(Chem0_clean3[[1]]), arr.ind=TRUE)
# deemedBad$Chem3 <- Chem0_clean3[[2]]
# 
# deemedBad.chem <- which(is.na(Chem0_clean3[[1]]))
# write.table(deemedBad.chem, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.chem.txt", sep="\t", row.names=FALSE)
# deemedBad.chem <- which(is.na(Chem0_clean3[[1]]))
deemedBad.chem0 <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/deemedBad.chem.txt", sep="\t", header=TRUE)[,]
deemedBad.chem <- arrayInd(deemedBad.chem0, dim(Chem0), dimnames(Chem0))
Chem0[deemedBad.chem] <- NA


#average epilimnion
Chem1 <- ddply(Chem0, .variables=c("lakeid", "year4", "daynum", "sampledate"), .fun=EpiMean2)
names(Chem1) <- c("lakeid", "year4", "daynum", "sampledate", names(Chem0)[!is.element(names(Chem0), c("lakeid", "year4", "daynum", "sampledate"))])
names(Chem1)[5] <- "depth_chem"

Chem1 <- Chem1[,c("lakeid", "year4", "daynum", AllChem)]
Chem <- cbind("taxID"="Chem", "taxLvl"="Community", reshape2(Chem1, varying=AllChem, times=AllChem, v.names="Data", direction="long"))
names(Chem) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")



# ==========
# = Secchi =
# ==========
Secchi000 <- read.csv("Secchi_Nrtn.csv", colClasses=c("factor", "integer", "integer", "factor", "integer", "numeric", "character", "character", "numeric", "character", "numeric", "numeric", "numeric", "numeric")) #WHY ONLY NORTHERN?
Secchi00 <- aggregate(Secchi000[,c("secnview", "waveht", "cloud", "ice")], by=Secchi000[,rev(c("lakeid", "year4", "daynum", "sampledate"))], FUN=mean, na.rm=TRUE)
Sec_notNA <- which(!is.na(Secchi00[,"secnview"]))
Secchi0 <- Secchi00[Sec_notNA, c("lakeid", "year4", "daynum", "sampledate", "secnview", "waveht", "cloud", "ice")]
names(Secchi0) <- c("lakeid", "year4", "daynum", "sampledate", "Secchi", "waveht_secchi", "cloud_secchi", "ice_secchi")

# Secchi0_clean <- manClean(Secchi0, c("Secchi"))
# Secchi <- Secchi0_clean[[1]]
# deemedBad$Secchi <- Secchi0_clean[[2]]

Secchi0 <- Secchi0[,c("lakeid", "year4", "daynum", "Secchi")]
Secchi <- cbind("taxID"="Secchi", "taxLvl"="Community", reshape2(Secchi0, varying="Secchi", times="Secchi", v.names="Data", direction="long"))
names(Secchi) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# ====================
# = Light Extinction =
# ====================
LiExt000 <- read.csv("LightExt_Nrtn.csv", colClasses=c("factor", "integer", "integer", "factor", "numeric", "character", "character")) #WHY ONLY NORTHERN?
LightNoFlag <- which(LiExt000[, "lightext_flag"]=="")
LiExt1 <- LiExt000[LightNoFlag, c("lakeid", "year4", "daynum", "extcoef")]
LiExt <- cbind("taxID"="extcoef", "taxLvl"="Community", reshape2(LiExt1, varying="extcoef", times="extcoef", v.names="Data", direction="long"))
names(LiExt) <- c("taxID", "taxLvl", "location", "year4", "daynum", "variable", "Data")


# =======
# = Ice =
# =======
Ice_Nrtn000 <- read.csv("Ice_Nrtn.csv")
Ice_Nrtn00 <- Ice_Nrtn000[which(!is.na(Ice_Nrtn000[,"firstopen"]) & !is.na(Ice_Nrtn000[,"lastopen"])),]
OpenDuration_Nrtn <- as.integer(difftime(as.Date(Ice_Nrtn00[,"datelastopen"]), as.Date(Ice_Nrtn00[,"datefirstopen"])))
Ice_Nrtn0 <- data.frame(Ice_Nrtn00[,c("lakeid", "sta", "year")], "DaysOpen"=OpenDuration_Nrtn)
Ice_Nrtn <- aggregate(Ice_Nrtn0[,"DaysOpen"], by=Ice_Nrtn0[,rev(c("lakeid", "year"))], mean, na.rm=TRUE)
names(Ice_Nrtn) <- c("year4", "lakeid", "DaysOpen")

Ice_Srtn000 <- read.csv("Ice_Srtn.csv")

# Fix ice dates
Ice_Srtn00 <- FixIceDates(Ice_Srtn000)

# Calculate the days open
Ice_Srtn0 <- ddply(.data=Ice_Srtn00, .variables="lakeid", .fun=CalcDaysOpen)
Ice_Srtn <- Ice_Srtn0[,c("year4", "lakeid", "DaysOpen")]
Ice1 <- rbind(Ice_Nrtn, Ice_Srtn)
Ice <- cbind("taxID"="Ice", "taxLvl"="Community","daynum"=1, reshape2(Ice1, varying="DaysOpen", times="DaysOpen", v.names="Data", direction="long"))
names(Ice) <- c("taxID", "taxLvl", "daynum", "year4", "location", "variable", "Data")


# =====================================
# = Combine data sets and define Type =
# =====================================
DataFull <- rbind(
	cbind("Type"="Meteorological",Met),
	cbind("Type"="Meteorological", Wind),
	cbind("Type"="Physical",Phys),
	cbind("Type"="Physical", Ice),
	cbind("Type"="Physical", LakLev),
	cbind("Type"="Physical", LiExt),
	cbind("Type"="Physical", Secchi),
	cbind("Type"="Chemical", Ions), 
	cbind("Type"="Chemical", Chem), 
	cbind("Type"="Biological", Chl), 
	cbind("Type"="Biological", Zoop), 
	cbind("Type"="Biological", Fish)
)
DataFull <- DataFull[order(DataFull[,"Type"], DataFull[,"taxLvl"], DataFull[,"taxID"], DataFull[,"variable"], DataFull[,"location"], DataFull[,"year4"], DataFull[,"daynum"]),]
row.names(DataFull) <- NULL
DataFull <- DataFull[,c("Type","taxLvl","taxID","location","variable","year4","daynum","Data")]

write.table(DataFull, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataFull.txt", sep="\t", row.names=FALSE)

# ======================
# = Take annual maxima =
# ======================
DataMax <- aggregate(DataFull[,"Data"], by=DataFull[,c("Type","taxID","taxLvl","location","year4","variable")], FUN=max, na.rm=TRUE)
DataMax <- Inf2NA(DataMax)
names(DataMax)[7] <- "Data"
DataMax <- DataMax[order(DataMax[,"Type"], DataMax[,"taxLvl"], DataMax[,"taxID"], DataMax[,"variable"], DataMax[,"location"], DataMax[,"year4"]),]
row.names(DataMax) <- NULL
DataMax <- DataMax[,c("Type","taxLvl","taxID","location","variable","year4","Data")]

write.table(DataMax, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataMax.txt", sep="\t", row.names=FALSE)





