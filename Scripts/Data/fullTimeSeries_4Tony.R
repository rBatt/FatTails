
library(data.table)
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/FatTails_Functions.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/Data_Functions.R")


full0 <- fread("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/FatTailsDataFull.txt")


# =================================
# = Subset the zoop and fish data =
# =================================
rmMass <- !full0[,variable]%in%c("tot_zoop_mass","cpue3_WeiEff")
rmHigher <- (!full0[,variable]%in%c("density","cpue1_Sum")) | (full0[,variable]%in%c("density","cpue1_Sum") & full0[,taxLvl]%in%c("Species","Genus","Family","Order","Class","Phylum"))
full <- full0[rmMass&rmHigher,]



# ================================
# = Remove all but Genus for Bio =
# ================================
full <- full[!taxLvl%in%c("Species","Family","Order","Class","Phylum"),]




# ==================
# = Handy Function =
# ==================
lu <- function(x){
	length(unique(x))
}

full <- full[is.finite(Data),]
full[,n.yr:=lu(year4), by=c("Type","taxLvl","taxID","location","variable")]
setkey(full, Type, taxLvl, taxID, location, variable, year4, daynum)

# ============================================
# = Trim to time series w/ at least 15 years =
# ============================================
full <- full[n.yr>=15,]


# ==========================
# = Trim to useful columns =
# ==========================
full <- full[,list(Type, taxID, location, variable, year4, daynum, n.yr, Data)]



# =========================
# = Write to txt for Tony =
# =========================
write.table(full, file="~/Documents/School&Work/WiscResearch/FatTails/Data/fullTimeSeries_4Tony.txt", sep="\t", row.names=FALSE)




# ===================================
# = Pad to regular interval w/ NA's =
# ===================================
# 
# full[,minDiff:=min(dist(unique(daynum))), by=c("Type","taxID","location","variable")]
# 
# filler <- full[,CJ(as.numeric(seq(min(daynum),max(daynum),by=unique(minDiff))), as.numeric(year4)),by=c("location","variable")]











