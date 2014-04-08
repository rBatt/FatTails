
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

# =========================
# = ARIMA Simulation Plot =
# =========================
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatARMA_Sim.R")
	# ========================
	# = Prepare for plotting =
	# ========================
# Calculate densities for GEV fit (Panel E)
adat <- c(tfTS[tmTS], ffTS[fmTS]) # adat is 'all data' – both thin and fat maxima time series
ma <- min(adat) # find the smallest value of the maxima
dSeq <- seq(ma-sign(ma)*1.5*ma, max(adat)*1.1, by=0.1) # create a sequence of values over which to calculate the density
fatGEV <- dgev(dSeq, xi=simXiS[fattestI,"Xi"], simXiS[fattestI,"mu"], simXiS[fattestI,"sig"]) # probs for fat time series using GEV fit (package 'evir')
fatGEV[!is.finite(fatGEV)] <- 0 
thinGEV <- dgev(dSeq, xi=simXiS[thinnestI,"Xi"], simXiS[thinnestI,"mu"], simXiS[thinnestI,"sig"]) # probs for thin time series
thinGEV[!is.finite(thinGEV)] <- 0

# Prepare the Xi values to be plotted on figure (Panel E)
tXi <- round(simXiS[thinnestI,"Xi"], 2)
fXi <- round(simXiS[fattestI,"Xi"], 2)
tLabXi <- parse(text=paste("xi", tXi, sep=" = "))
fLabXi <- parse(text=paste("xi", fXi, sep=" = "))

	# ===============================
	# = Plot Example w/ reversed xy =
	# ===============================
# Set up figure space
# dev.new(width=3.5, height=5) # open graphical device
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fat_conceptFig.png", res=150, units="in", height=5, width=3.5)
cols1 <- rep(rep(1:3, each=3), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(4,5,3), each=3), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2), ncol=7) # create layout matrix
layout(lmat) # define graphical device layout
par(mar=c(2.5,0.0,0.5,0.1), oma=c(0, 2, 0, 0.25), ps=10, cex=1, mgp=c(1, 0.3, 0), tcl=-0.25, family="Times") # set graphical parameters

# Plot thin-tailed time series
ylim1 <- range(tfTS)*c(1, 1.15)
plot(tfTS, type="l", col="gray", xlab="", xaxt="n", ylab="", ylim=ylim1, bty="l")
text(0.025*length(tfTS), y=max(tfTS)*1.05, "A", font=2)
ats1 <- axTicks(1)
labs1 <- ats1/nPerYear
axis(side=1, labels=labs1, at=ats1)
points(tmTS, tfTS[tmTS], col="blue")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot fat-tailed time series
ylim2 <- range(ffTS)*c(1, 1.15)
plot(ffTS, type="l", col="gray", xlab="", ylab="", xaxt="n", ylim=ylim2, bty="l")
text(0.025*length(ffTS), y=max(ffTS)*1.05, "B", font=2)
ats2 <- axTicks(1)
labs2 <- ats2/nPerYear
axis(side=1, labels=labs2, at=ats2)
points(fmTS, ffTS[fmTS], col="red")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot GEV densities
colorPoly(quants=dSeq, dents=list(thinGEV,fatGEV), cols=c("blue", "red"), bty="l")
cm0 <- min(dSeq)
text(x=cm0+sign(cm0)*cm0*0.05, y=max(c(thinGEV,fatGEV))*0.95, "E", font=2)
mtext("density", side=2, line=1.25)
mtext("y", side=1, line=1.25)
cma0 <- max(dSeq)
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.75, bquote(xi~"="~.(tXi)), pos=4, col="blue")
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.58, bquote(xi~"="~.(fXi)), pos=4, col="red")

# Plot empirical densities for thin-tailed
colorDens(vals=list(tfTS, tfTS[tmTS]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1)
cm1 <- min(tfTS)
text(y=cm1+sign(cm1)*cm1*0.15, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
mtext("density", side=1, line=1.25)

# Plot empirical densities for fat-tailed
cm2 <- min(ffTS)
colorDens(vals=list(ffTS, ffTS[fmTS]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim2)
text(y=cm2+sign(cm2)*cm2*0.15, x=0.25*max(density(ffTS)$y, density(ffTS[fmTS])$y), "D", font=2)
mtext("density", side=1, line=1.25)
dev.off()







