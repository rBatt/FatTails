

load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
library("beanplot")

eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))


# ==================================
# = Plot Xi Across taxonomic level =
# ==================================
# dev.new(width=3.5, height=6)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/xi_taxClass.png", res=150, units="in", height=6, width=3.5)
par(mfrow=c(2,1), mar=c(3,3,0.5,0.5), mgp=c(1.5, 0.4, 0), tcl=-0.25, family="Times", cex=1, ps=10)
plotTax(D=zoop.gev.full, V="density", Stat="sh_0", taxCol="Phylum", ylim=c(-0.5, 1.5))
ztLab <- c("Spec","Genus","Fam","Ord","Class","Phy","Comm")
text(1:length(ztLab), par("usr")[3]-0.09, labels=ztLab, pos=1, xpd=TRUE, srt=45)
mtext(bquote(xi), side=2, line=1.5)

plotTax(D=fish.gev.full, V="cpue1_Sum", Stat="sh_0", taxCol="Order")
ftLab <- c("Spec","Genus","Fam","Ord","Comm")
text(1:length(ftLab), par("usr")[3]-0.09, labels=ftLab, pos=1, xpd=TRUE, srt=45)
mtext(bquote(xi), side=2, line=1.5)
dev.off()

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
png(paste("Figures","compareTimes_sigShape.png",sep="/"), width=3.5, height=6, units="in", res=300, bg="white")
par(mfrow=c(3,1), mar=c(0.5, 1.25, 0.25, 0.25), mgp=c(1, 0.3, 0), tcl=-0.3, cex=1, ps=10, family="Times", oma=c(1.5,1,0,0), xpd=T)
pDens(bound, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Bounded~(xi<0)))
text(8, c(0.4, 0.5, 0.6), c("GEV", "Normal", "Log-Normal"), col=cLine)
pDens(thin, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Thin~(xi==0)))
pDens(fat, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8, 0.95, bquote(Fat~(xi>0)))
mtext("Waiting Time (years)", side=1, outer=TRUE, line=0.5)
mtext("Relative Density", side=2, line=0.0, outer=TRUE)
dev.off()


# ==============================
# = Xi / Waiting Time Beanplot =
# ==============================
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

beanplot(log10(Level2_time)~Type, data=data.fat[is.finite(data.fat[,"Level2_time"]),], log="", ylab="", xaxt="n", yaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
wtBase <- axTicks(2)
wtLab <- parse(text=paste(10,wtBase,sep="^"))
axis(side=2, at=wtBase, labels=wtLab)
axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)
dev.off()


# =============================
# = ARMA Residual Xi Beanplot =
# =============================
# dev.new(width=3.5, height=3.5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/ARMAxi_beanplot.png", res=150, units="in", height=3.5, width=3.5)
par(mfrow=c(1,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(residual_sh_0~Type, data=data.2, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(xi~~of~~ARMA~Residuals), side=2, line=1.5)
dev.off()


# ======================
# = Signif Xi Beanplot =
# ======================
# dev.new(width=3.5, height=3.5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/signifXi_beanplot.png", res=150, units="in", height=3.5, width=3.5)
par(mfrow=c(1,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(shape.sig~Type, data=data.fat, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5, bw=0.05, what=c(1,1,1,1), maxstripline=0.05)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Critical~~Value~~'for'~~xi), side=2, line=1.5)
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






# ===============================
# = Dscat of Xi vs. log10(InfE) =
# ===============================
dev.new(height=3.5, width=3.5)
dscat(log(data.2[,"InfE"], base=10), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(InfE), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=1, line=1.5)

# ==========================
# = Dscat of Xi vs. Lambda =
# ==========================
dev.new(height=3.5, width=3.5)
dscat(sqrt(data.2[,"Lambda"]), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(Lambda), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(.(lNota2))), side=1, line=1.5)

# ===================================
# = Dscat of Xi vs. Xi of residuals =
# ===================================
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/dscat_xiTS_xiRES.png", res=150, units="in", height=3.5, width=3.5)
# dev.new(height=3.5, width=3.5)
dscat(data.2[,"residual_sh_0"], data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~residual_sh_0, data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(xi[ARMA~Residuals]), side=1, line=1.5)
dev.off()

# =======================
# = Dscat of Xi vs SigE =
# =======================
dev.new(height=3.5, width=3.5)
dscat(sqrt(data.2[,"sigE"]), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(sigE), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(sigma[E])), side=1, line=1.5)




# ==========================================================
# = Proportion bounded/thin/fat in at each taxonomic level =
# ==========================================================

png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/tailedness_taxLevel.png", res=150, units="in", height=4, width=8)
par(mfrow=c(1,2), mar=c(2,2.5,0.5,0.5), mgp=c(1.5,0.5,0), tcl=-0.3, ps=10, family="Times")

sXi.fish <- tapply(fish.gev.full[,"shape.sig"], fish.gev.full[,"taxLvl"], sign)
prop.fat.fish0 <- sapply(sXi.fish, function(x)c("bounded"=sum(x<0), "thin"=sum(x==0), "fat"=sum(x>0))/length(x))[,-1]
prop.fat.fish <- cbind("taxLvl"=names(prop.fat.fish0),t(prop.fat.fish0))
barplot(prop.fat.fish, beside=TRUE, legend=TRUE, args.legend=list(x="topleft",title="Fish"), ylab="Proportion")

sXi.zoop <- tapply(zoop.gev.full[,"shape.sig"], zoop.gev.full[,"taxLvl"], sign)
prop.fat.zoop0 <- sapply(sXi.zoop, function(x)c("bounded"=sum(x<0), "thin"=sum(x==0), "fat"=sum(x>0))/length(x))[,-1]
prop.fat.zoop <- cbind("taxLvl"=names(prop.fat.zoop0),t(prop.fat.zoop0))
barplot(prop.fat.zoop, beside=TRUE, legend=TRUE, args.legend=list(x="topleft",title="Zoops"))
dev.off()



# ==================================
# = Plot precipitation time series =
# ==================================
mad.precip <- data.max[data.max[,"variable"]=="precip_mm"&data.max[,"location"]=="Madison",]
plot(mad.precip[,"year4"], mad.precip[,"Data"], type="l")



fatQQ <- function(vari="cpue1_Sum", loca="ME", taxid="Lepomis", norm=FALSE){
	oQ <- data.max[data.max[,"variable"]==vari&data.max[,"location"]==loca&data.max[,"taxID"]==taxid,"Data"]
	oQ <- oQ[!is.na(oQ)]
	oFat <- data.fat[data.fat[,"variable"]==vari&data.fat[,"location"]==loca&data.fat[,"taxID"]==taxid,]
	exi <- oFat[,"sh_0"]
	escale <- oFat[,"sig_0"]
	eloc <- oFat[,"mu_0"]
	oN <- oFat[,"N"]
	tPs <- (1:oN)/(oN+1)
	tQ <- qgev(tPs, xi=exi, mu=eloc, sigma=escale)
	plot(1:oN, oQ, type="l", col="gray", lwd=2, xaxt="n", yaxt="n", xlab="", ylab="")
	par(new=TRUE)
	if(!norm){
		plot(tQ, sort(oQ))	
	}else{
		qqnorm(scale(oQ))
	}
	
	abline(a=0, b=1, lty="dotted")
}

fatQQ("cpue1_Sum", "ME", "Pimephales", norm=FALSE)








