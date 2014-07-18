

library(plyr)
library(beanplot)
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/finalFrame.RData")
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")


fatLog00 <- ddply(finalFrame, c("Type","taxLvl","taxID","location","variable"), GEV)
fatLog00[,"Type"] <- factor(fatLog00[,"Type"], levels=c("Biological", "Chemical", "Physical", "Meteorological"))
fatLog00[,"taxLvl"] <- factor(fatLog00[,"taxLvl"], levels=c("Community","Phylum","Class","Order","Family","Genus","Species"))

data.stat.full.log <- data.stat.full
data.stat.full.log[,"Data"] <- log(data.stat.full[,"Data"])
data.statLog <- Inf2NA(ddply(data.stat.full.log, c("Type","taxLvl","taxID","location","variable"), calc.level))

test <- data.stat.full.log
names(test)[names(test)=="Data"] <- "dataLog"
test <- merge()


fatLog0 <- merge(fatLog00, data.statLog, all.x=TRUE)

# Return time calculated from GEV
fatLog <- ddply(fatLog0, c("Type","taxLvl","taxID","location","variable"), lvl_return, returnFull=TRUE)
# fatLog <- fatLog[complete.cases(fatLog),]

bLine <- rainbow(n=5, v=0.8, s=1)
bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
beanCol <- list(c(bFill[1]),
				c(bFill[2]),
				c(bFill[3]),
				c(bFill[4])
				)


par(mfrow=c(2,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(sh_0~Type, data=fatLog, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(xi~~from~~GEV), side=2, line=1.5)

beanplot(log10(Level2_time)~Type, data=fatLog[is.finite(fatLog[,"Level2_time"]),], log="", ylab="", xaxt="n", yaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
wtBase <- axTicks(2)
wtLab <- parse(text=paste(10,wtBase,sep="^"))
axis(side=2, at=wtBase, labels=wtLab)
axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)


