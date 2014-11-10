# Extra plots I'm making on the airplane for the JASM talk 
# Mostly just interested in creating time series for the fattest Fish, Zoop, and Chla time series
# From manuscript draft:
# The fish with the largest estimate of xi were Lepomis spp. in Lake Mendota (xi = 1.2; 90% CI: 0.7 < xi < 1.8), the fattest-tailed zooplankton were Kellicottia spp. in Crystal Bog (xi = 1.2; 90% CI: 0.7 < xi < 1.7), and Crystal Lake had the most fat-tailed chlorophyll time series (xi = 0.5; 90% CI: 0.2 < xi < 0.8).

load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatForest.RData") # also contains data.2
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/fatPlot_Functions.R")
library("beanplot")


# ====================
# = Fattest Examples =
# ====================
lepLog <- data.stat[,"location"]=="MO" & data.stat[,"variable"]=="cpue1_Sum" & data.stat[,"taxID"]=="Lepomis"
fLep <- data.stat[lepLog,]

kelLog <- data.stat[,"location"]=="CB" & data.stat[,"variable"]=="density" & data.stat[,"taxID"]=="Kellicottia"
fKel <- data.stat[kelLog,]

chlLog <- data.stat[,"location"]=="CR" & data.stat[,"variable"]=="chlor"
fChl <- data.stat[chlLog,]

# dev.new(width=3.5, height=7)
png(filename="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Other/egFat_lep_kel_chl.png", res=200, width=7, height=2, units="in")
par(mfrow=c(1,3), ps=10, mar=c(2,3,0.5,0.5), oma=c(0.5,0,0,0), cex=1, family="Times", mgp=c(1.1, 0.35, 0), tcl=-0.35)
plot(fLep[,"year4"], fLep[,"Data"], xlab="", ylab=bquote(Lake~~Mendota~~italic(Lepomis)), type="l", col="gray")
points(fLep[,"year4"], fLep[,"Data"], pch=20)
text(1994, max(fLep[,"Data"])*0.85, labels=bquote(atop(xi~'='~1.2,({0.6<xi}<1.8))), cex=0.85, pos=4)

plot(fKel[,"year4"], fKel[,"Data"], xlab="", ylab=bquote(Crystal~~Bog~~italic(Kellicottia)), type="l", col="gray")
points(fKel[,"year4"], fKel[,"Data"], pch=20)
text(1982, max(fKel[,"Data"])*0.85, labels=bquote(atop(xi~'='~1.2,({0.7<xi}<1.7))), cex=0.85, pos=4)

plot(fChl[,"year4"], fChl[,"Data"], xlab="", ylab=bquote(Crystal~~Lake~~Chlorophyll), type="l", col="gray")
points(fChl[,"year4"], fChl[,"Data"], pch=20)
text(1980, max(fChl[,"Data"])*0.85, labels=bquote(atop(xi~'='~0.5,({0.2<xi}<0.8))), cex=0.85, pos=4)
mtext("Year", side=1, line=-0.5, outer=TRUE)
dev.off()





# ==============================================
# = Compare Waiting Times using significant Xi =
# ==============================================
dTimes <- c("Level2_time","Level2_normTime","Level2_logNormTime")
boundLog <- data.fat[,"shape.sig"]< 0
bound <- data.fat[boundLog,]
fatLog <- data.fat[,"shape.sig"] > 0
fat <- data.fat[fatLog,]
thinLog <- data.fat[,"shape.sig"] == 0
thin <- data.fat[thinLog,]

cLine <- rainbow(n=3, v=0.8, s=1)
cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)

# dev.new(width=3.5, height=6)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Other/compareTimes_sigShape_JASM.png", width=6, height=2, units="in", res=300, bg="white")
par(mfrow=c(1,3), mar=c(0.5, 1.25, 0.25, 0.25), mgp=c(1, 0.3, 0), tcl=-0.3, cex=1, ps=9, family="Times", oma=c(1.5,1,0,0), xpd=T)
pDens(bound, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8.25, 0.95, bquote(Bounded~(xi<0)))
text(8.5, c(0.4, 0.5, 0.6), c("GEV", "Normal", "Log-Normal"), col=cLine)
pDens(thin, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8, 0.95, bquote(Thin~(xi==0)))
pDens(fat, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8, 0.95, bquote(Fat~(xi>0)))
mtext("Waiting Time (years)", side=1, outer=TRUE, line=0.5)
mtext("Relative Density", side=2, line=0.0, outer=TRUE)
dev.off()





