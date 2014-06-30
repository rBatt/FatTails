

load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatForest.RData") # also contains data.2
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



# ====================================
# = Fat Forest: Zoop Marginal Effect =
# ====================================
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/zoop_Marginal.png", res=150, units="in", height=5, width=3.5)
# dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35)
plot(zoop.pp.sesh0, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Zooplankton~~xi))
plot(zoop.pp.N, type="l", xlab=bquote(Sample~~size),ylab=bquote(Zooplankton~~xi))
dev.off()


# ====================================
# = Fat Forest: Fish Marginal Effect =
# ====================================
# dev.new(width=3.5, height=5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fish_Marginal.png", res=150, units="in", height=5, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.5,0.3,0), tcl=-0.25)
plot(fish.pp.sesh0, , type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Fish~~xi))
plot(fish.pp.loc$y, xlab=bquote(Lake~~name), xaxt="n", pch=19, ylab=bquote(Fish~~xi))
axis(side=1, at=1:length(fish.pp.loc$x), labels=FALSE)
text(1:length(fish.pp.loc$x), par("usr")[3]-0.0025, labels=fish.pp.loc$x, pos=1, xpd=TRUE, srt=45)
dev.off()


# =========================================
# = Fat Forest: Full Data Marginal Effect =
# =========================================
# dev.new(width=3.5, height=5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fullData_Marginal.png", res=150, units="in", height=5, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35, cex=1)
plot(data.pp.sesh, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(xi))
plot(data.pp.type$y, xlab=bquote(Variable~~type), xaxt="n", pch=19, ylab=bquote(xi))
axis(side=1, at=1:length(data.pp.type$x), labels=FALSE)
lakeShort <- c("Physical"="Phys","Biological"="Bio", "Meteorological"="Met", "Chemical"="Chem")
text(1:length(data.pp.type$x), par("usr")[3]-0.0002, labels=lakeShort[as.character(data.pp.type$x)], pos=1, xpd=TRUE, srt=45)
dev.off()




