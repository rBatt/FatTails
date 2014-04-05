
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
library("beanplot")

# =======================
# = Some quick plotting =
# =======================
# V = "density"
# D = zoop.gev
# Stat = "sh_0"
# taxCol = "Phylum"
# dev.new()
# par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(1.5, 0.4, 0), tcl=-0.25, family="Times", cex=1, ps=10)
# plotTax(D=zoop.gev, V="density", Stat="sh_0", taxCol="Phylum", legendTitle="# Individuals", ylim=c(-0.5, 1.5))
# ztLab <- c("Spec","Genus","Fam","Ord","Class","Phy","Comm")
# text(1:length(ztLab), par("usr")[3]-0.09, labels=ztLab, pos=1, xpd=TRUE, srt=45)
# mtext(bquote(xi), side=2, line=1.5)
# 
# plotTax(D=zoop.gev, V="tot_zoop_mass", Stat="sh_0", taxCol="Phylum", ylim=c(-1.5, 2), legendTitle="Biomass")
# text(1:length(ztLab), par("usr")[3]-0.09, labels=ztLab, pos=1, xpd=TRUE, srt=45)
# mtext(bquote(xi), side=2, line=1.5)
# 
# plotTax(D=fish.gev, V="cpue1_Sum", Stat="sh_0", taxCol="Community", legendTitle="# Individuals")
# ftLab <- c("Spec","Genus","Fam","Ord","Comm")
# text(1:length(ftLab), par("usr")[3]-0.09, labels=ftLab, pos=1, xpd=TRUE, srt=45)
# mtext(bquote(xi), side=2, line=1.5)
# 
# plotTax(D=fish.gev, V="cpue3_WeiEff", Stat="sh_0", taxCol="Community", legendTitle="Biomass")
# text(1:length(ftLab), par("usr")[3]-0.09, labels=ftLab, pos=1, xpd=TRUE, srt=45)
# mtext(bquote(xi), side=2, line=1.5)
# 
# density.gev <- bio.gev[bio.gev[,"variable"]=="density",]
# density.gev[,"taxID"] <- as.character(density.gev[,"taxID"])
# bioComm.gev <- bio.gev[bio.gev[,"taxLvl"]=="Community",]
# bioComm.gev[,"taxID"] <- factor(bioComm.gev[,"taxID"])
# ===========================
# = Xi Waiting Time Boxplot =
# ===========================
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fatBoxXiWaiting.png", res=150, units="in", height=6, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
boxplot(sh_0~Type, data=data.fat, outline=FALSE, ylab="", yaxt="n", xaxt="n")
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(xi~~from~~GEV), side=2, line=1.5)

boxplot(log10(Level2_time)~Type, data=data.fat, outline=FALSE, ylab="", xaxt="n", yaxt="n")
wtBase <- pretty(log10(data.fat[,"Level2_time"]))
wtLab <- parse(text=paste(10,wtBase,sep="^"))
axis(side=2, at=wtBase, labels=wtLab)
axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)
dev.off()


# Instead of Xi, uses the 5th percentile from the norm w/ a mean of Xi (from the GEV fit) and a sd of the se from the GEV fit
# par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
# boxplot(shape.sig~Type, data=data.fat, outline=FALSE, ylab="", yaxt="n", xaxt="n")
# axis(side=2)
# axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(xi~~from~~GEV), side=2, line=1.5)
# 
# boxplot(log10(Level2_time)~Type, data=data.fat, outline=FALSE, ylab="", xaxt="n", yaxt="n")
# wtLab <- pretty(log10(data.fat[,"Level2_time"]))
# axis(side=2, at=wtLab, labels=10^wtLab)
# axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)



# ==============================================
# = Compare Waiting Times across Distributions =
# ==============================================
dTimes <- c("Level2_time","Level2_normTime","Level2_logNormTime")
cutoff <- 0.35
boundLog <- data.fat[,"sh_0"]<= -cutoff
bound <- data.fat[boundLog,]
fatLog <- data.fat[,"sh_0"] > cutoff
fat <- data.fat[fatLog,]
thinLog <- data.fat[,"sh_0"] > -cutoff & data.fat[,"sh_0"] < cutoff #!boundLog & !fatLog
thin <- data.fat[thinLog,]

cLine <- rainbow(n=3, v=0.8, s=1)
cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)

# dev.new(width=3.5, height=6)
png(paste("Figures","compareTimes.png",sep="/"), width=3.5, height=6, units="in", res=300, bg="white")
par(mfrow=c(3,1), mar=c(0.5, 1.25, 0.25, 0.25), mgp=c(1, 0.3, 0), tcl=-0.3, cex=1, ps=10, family="Times", oma=c(1.5,1,0,0), xpd=T)
pDens(bound, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Bounded~(xi<=-.(cutoff))))
text(8, c(0.4, 0.5, 0.6), c("GEV", "Normal", "Log-Normal"), col=cLine)
pDens(thin, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Thin~(-.(cutoff)<xi~phantom()<.(cutoff))))
pDens(fat, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8, 0.95, bquote(Fat~(xi>=.(cutoff))))
mtext("Waiting Time (years)", side=1, outer=TRUE, line=0.5)
mtext("Relative Density", side=2, line=0.0, outer=TRUE)
dev.off()


# ================
# = Try Beanplot =
# ================
bLine <- rainbow(n=5, v=0.8, s=1)
bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
beanCol <- list(c(bFill[1]),
				c(bFill[2]),
				c(bFill[3]),
				c(bFill[4])
				)

# dev.new(width=3.5, height=6)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fatBeanXiWaiting.png", res=150, units="in", height=6, width=3.5)
par(mfrow=c(2,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(sh_0~Type, data=data.fat, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(xi~~from~~GEV), side=2, line=1.5)

beanplot(log10(Level2_time)~Type, data=data.fat[is.finite(data.fat[,"Level2_time"]),], log="", ylab="", xaxt="n", yaxt="n", border=bLine, col=beanCol, ll=0.01)
wtBase <- axTicks(2)
wtLab <- parse(text=paste(10,wtBase,sep="^"))
axis(side=2, at=wtBase, labels=wtLab)
axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)
dev.off()




