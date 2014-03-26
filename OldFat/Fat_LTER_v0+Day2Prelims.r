#Fat analysis of some LTER database data
rm(list=ls())
graphics.off()
setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails")
Data_Chl_Phaeo <- read.table("Fat_LTER_Data/LTER_Chl_Phaeo.txt", sep=",", header=TRUE)
Data_GWaterLvl <- read.table("Fat_LTER_Data/LTER_Groundwater_Level.txt", sep=",", header=TRUE)
Data_Kd <- read.table("Fat_LTER_Data/LTER_Light_Extinction.txt", sep=",", header=TRUE)
Data_Secchi <- read.table("Fat_LTER_Data/LTER_Secchi_Base.txt", sep=",", header=TRUE)
Data_Mad_Meteo <- read.csv("Fat_LTER_Data/madison_daily_meteorological_data_.csv", sep=",", header=TRUE)
Data_Mad_Ice <- read.csv("Fat_LTER_Data/north_temperate_lakes_lter__ice_duration_-_madison_lakes_area.csv", sep=",", header=TRUE)
Data_Trout_Ice <- read.csv("Fat_LTER_Data/north_temperate_lakes_lter__ice_duration_-_trout_lake_area.csv", sep=",", header=TRUE)
Data_Mad_Zoop <- read.csv("Fat_LTER_Data/north_temperate_lakes_lter__zooplankton_-_madison_lakes_area.csv", sep=",", header=TRUE)

Data_Chl_Phaeo <- subset(Data_Chl_Phaeo, CHLOR >0 & PHAEO>-10000 & DEPTH<=10)# &LAKEID=="TB"
AvgChl <- aggregate(Data_Chl_Phaeo$CHLOR, by=list(Data_Chl_Phaeo$LAKEID, Data_Chl_Phaeo$YEAR4, Data_Chl_Phaeo$DAYNUM), FUN=mean)
plot.density(density(AvgChl$x))

AvgPhaeo <- aggregate(Data_Chl_Phaeo$PHAEO, by=list(Data_Chl_Phaeo$LAKEID, Data_Chl_Phaeo$YEAR4, Data_Chl_Phaeo$DAYNUM), FUN=mean)
plot.density(density(AvgPhaeo$x))

Data_GWaterLvl <- subset(Data_GWaterLvl,  FLAG!="N" & FLAG!="D" & FLAG!="F")
#AvgGWLvl <- aggregate(Data_GWaterLvl$CHLOR, by=list(Data_GWaterLvl$SAMPLEDATE, Data_GWaterLvl$YEAR4, Data_GWaterLvl$WELLID), FUN=mean)
plot.density(density(Data_GWaterLvl$WELL_LEVEL))

Data_Secchi <- subset(Data_Secchi, SECNVIEW >-10000)
AvgSecchi <- aggregate(Data_Secchi$SECNVIEW, by=list(Data_Secchi$LAKEID, Data_Secchi$YEAR4, Data_Secchi$DAYNUM), FUN=mean)
plot.density(density(AvgSecchi$x))

Data_Kd <- subset(Data_Kd, EXTCOEF >-10000)
AvgKd <- aggregate(Data_Kd$EXTCOEF, by=list(Data_Kd$LAKEID, Data_Kd$YEAR4, Data_Kd$DAYNUM), FUN=mean)
plot.density(density(AvgKd$x))


MeanExcess <- function(x){
	T  <- x# <- seq(min(x), max(x), length.out=200)
	MEs <- rep(NA, length(T))
	for(i in 1:length(T)){
		MEs[i] <- mean(x[which(x>T[i])]-T[i])
	}
	ME <- data.frame("T"=T, "MEs"=MEs)
	return(ME[-200,])
}


Data_Mad_Zoop <- subset(Data_Mad_Zoop, density>0 & avg_length>0)
Data_Mad_Zoop[,"dens_length"] <- Data_Mad_Zoop[,"avg_length"] * Data_Mad_Zoop[,"density"]
Total_DensLength <- aggregate(Data_Mad_Zoop$dens_length, by=list(Data_Mad_Zoop$lakeid, Data_Mad_Zoop$sample_date, Data_Mad_Zoop$year4), FUN=sum)
TotalDensity <- aggregate(Data_Mad_Zoop$density, by=list(Data_Mad_Zoop$lakeid, Data_Mad_Zoop$sample_date, Data_Mad_Zoop$year4), FUN=sum)
Mad_Zoop_WeiAvg_Length <- Total_DensLength$x / TotalDensity$x

ME_Mad_Zoop_WeiAvg_Length <- MeanExcess(Mad_Zoop_WeiAvg_Length)
dev.new()
plot(ME_Mad_Zoop_WeiAvg_Length[,1], ME_Mad_Zoop_WeiAvg_Length[,2])

ME_Mad_Zoop_dens_length <- MeanExcess(Data_Mad_Zoop[,"dens_length"])
dev.new()
plot(ME_Mad_Zoop_dens_length[,1], ME_Mad_Zoop_dens_length[,2])

# ME_Mad_Zoop <- MeanExcess(subset(Data_Mad_Zoop, species_name=="DAPHNIA PULICARIA" & density>0 & lakeid=="MO")[,"density"])
# dev.new()
# plot(ME_Mad_Zoop[,1], ME_Mad_Zoop[,2])

# ME_Mad_Precip <- MeanExcess(Data_Mad_Meteo[,"precip_raw_mm"])
# dev.new()
# plot(ME_Mad_Precip[,1], ME_Mad_Precip[,2])

# ME_Mad_SnowFall <- MeanExcess(subset(Data_Mad_Meteo, snow_raw_cm>-1)[,"snow_raw_cm"])
# dev.new()
# plot(ME_Mad_SnowFall[,1], ME_Mad_SnowFall[,2])

# ME_Mad_SnowDepth <- MeanExcess(subset(Data_Mad_Meteo, snow_depth_cm>-1)[,"snow_depth_cm"])
# dev.new()
# plot(ME_Mad_SnowDepth[,1], ME_Mad_SnowDepth[,2])

# ME_Mad_AirTempMax <- MeanExcess(Data_Mad_Meteo[,"max_air_temp_adjusted"])
# dev.new()
# plot(ME_Mad_AirTempMax[,1], ME_Mad_AirTempMax[,2])

# ME_Mad_AirTempMin <- MeanExcess(Data_Mad_Meteo[,"min_air_temp_adjusted"])
# dev.new()
# plot(ME_Mad_AirTempMin[,1], ME_Mad_AirTempMin[,2])

# ME_Mad_AirTempRange <- MeanExcess(Data_Mad_Meteo[,"range_air_temp_adjusted"])
# dev.new()
# plot(ME_Mad_AirTempRange[,1], ME_Mad_AirTempRange[,2])

# ME_Mad_IceDur <- MeanExcess(subset(Data_Mad_Ice, ice_duration>-1)[,"ice_duration"])
# dev.new()
# plot(ME_Mad_IceDur[,1], ME_Mad_IceDur[,2])

# ME_Mad_IceFree <- MeanExcess(365-subset(Data_Mad_Ice, ice_duration>-1)[,"ice_duration"])
# dev.new()
# plot(ME_Mad_IceDur[,1], ME_Mad_IceDur[,2])









ME_Chl <- MeanExcess(AvgChl$x)
dev.new()
par()
plot(ME_Chl[,1], ME_Chl[,2])

ME_Phaeo <- MeanExcess(AvgPhaeo$x)
dev.new()
plot(ME_Phaeo[,1], ME_Phaeo[,2])

ME_Kd <- MeanExcess(AvgKd$x)
dev.new()
plot(ME_Kd[,1], ME_Kd[,2])

ME_Secchi <- MeanExcess(AvgSecchi$x)
dev.new()
plot(ME_Secchi[,1], ME_Secchi[,2])

ME_GWLvl <- MeanExcess(Data_GWaterLvl$WELL_LEVEL)
dev.new()
plot(ME_GWLvl[,1], ME_GWLvl[,2])

ME_Nrml <- MeanExcess(rnorm(n=10000))
dev.new()
plot(ME_Nrml[,1], ME_Nrml[,2])


# UniqueLakes <- unique(subset(Data_Chl_Phaeo, select=LAKEID)[,1])
# for(i in 1:length(UniqueLakes)){
	# TempoChl <- subset(Data_Chl_Phaeo, LAKEID==UniqueLakes[i])
	# TempoDatesChl <- subset(TempoChl, select=YEAR4)[,1] + subset(TempoChl, select=DAYNUM)[,1]/1000
	# TempoChl[,"TempoDatesChl"] <- TempoDates
	# FunIndexChl <- 
	
	# for(j in 1:length(unique(TempoDatesChl))){
		
	# }
	
# }
# UniqueDates <- unique(subset(Data_Chl_Phaeo, select=YEAR4)[,1] + subset(Data_Chl_Phaeo, select=DAYNUM)[,1]/1000)
# Unique



dev.new()
par(mfrow=c(2,2), family="Times", las=1, mar=c(5,3,3,1), oma=c(2,2,0,0))

#ME_Secchi <- MeanExcess(AvgSecchi$x)
plot(ME_Secchi[,1], ME_Secchi[,2], xlab="Secchi Depth (m)", ylab="", bty="l")
text(x=15, y=3.5, labels="A", cex=1.5)

#ME_Mad_Zoop_WeiAvg_Length <- MeanExcess(Mad_Zoop_WeiAvg_Length)
plot(ME_Mad_Zoop_WeiAvg_Length[,1], ME_Mad_Zoop_WeiAvg_Length[,2], xlab="Mean Zooplankton Length (mm)", ylab="", bty="l")
text(x=2.2, y=0.62, labels="B", cex=1.5)

#ME_Chl <- MeanExcess(AvgChl$x)
plot(ME_Chl[,1], ME_Chl[,2], xlab="Chlorophyll (mg/L)", ylab="", bty="l")
text(x=187, y=46, labels="C", cex=1.5)

#ME_Mad_Zoop_dens_length <- MeanExcess(Data_Mad_Zoop[,"dens_length"])
plot(ME_Mad_Zoop_dens_length[,1]/100000, ME_Mad_Zoop_dens_length[,2]/100000, xlab="Total Zooplankton Mass (index)", ylab="", bty="l")
text(x=15, y=5.7, labels="D", cex=1.5)
mtext("Mean Excess", outer=TRUE, las=0, side=2, line=0, cex=1.5)
mtext("Threshold", outer=TRUE, las=0, side=1, line=0, cex=1.5)

# setwd("/Users/Battrd/Documents/School&Work/GradSchool/DissertationProposal")
# dev2bitmap("MeanExcess_4_DissProp.tif", type="tifflzw", height=6, width=6, pointsize=12, res=300, method="pdf")

#below added on day2 of prelims
UniqueLakes <- unique(Data_Mad_Zoop[,"lakeid"])
UniqueLakes_Chl <- unique(AvgChl[,"Group.1"])
for(i in 1:length(UniqueLakes)){
	assign(paste(UniqueLakes[i],"MeanExcess_ZoopBiomass",sep="_"), MeanExcess(subset(Data_Mad_Zoop, lakeid==UniqueLakes[i], select="dens_length")[,1]))
	
}

dev.new()
par(mfrow=c(2,2), las=0)
plot(ME_MeanExcess_ZoopBiomass)
plot(MO_MeanExcess_ZoopBiomass)
plot(WI_MeanExcess_ZoopBiomass)
plot(FI_MeanExcess_ZoopBiomass)

UniqueLakes_Chl <- unique(AvgChl[,"Group.1"])
for(i in 1:length(UniqueLakes_Chl)){
	assign(paste(UniqueLakes_Chl[i],"MeanExcess_Chl",sep="_"), MeanExcess(subset(AvgChl, Group.1==UniqueLakes_Chl[i], select="x")[,1]))
	
}

dev.new()
par(mfrow=c(3,3), las=0)
plot(AL_MeanExcess_Chl)
plot(BM_MeanExcess_Chl)
plot(CB_MeanExcess_Chl)
plot(CR_MeanExcess_Chl)
plot(SP_MeanExcess_Chl)
plot(TB_MeanExcess_Chl)
plot(TR_MeanExcess_Chl)



ME_Mad_Zoop_dens_length <- MeanExcess(Data_Mad_Zoop[,"dens_length"])
dev.new()
plot(ME_Mad_Zoop_dens_length[,1], ME_Mad_Zoop_dens_length[,2])