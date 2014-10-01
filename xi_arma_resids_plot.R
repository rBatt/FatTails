

library(beanplot)
# load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.RData")

bLine <- rainbow(n=5, v=0.8, s=1)
bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
beanCol <- list(c(bFill[1]),
				c(bFill[2]),
				c(bFill[3]),
				c(bFill[4])
				)


# =============================
# = ARMA Residual Xi Beanplot =
# =============================
# dev.new(width=3.5, height=3.5)
# png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/ARMAxi_beanplot.png", res=150, units="in", height=3.5, width=3.5)
# par(mfrow=c(1,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
# beanplot(residual_sh_0~Type, data=data.2, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
# axis(side=2)
# axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(xi~~of~~ARMA~Residuals), side=2, line=1.5)
# dev.off()


# =====================================
# = Updated ARMA Residual Xi Beanplot =
# =====================================
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/ARMAxi_beanplot.png", res=150, units="in", height=3.5, width=3.5)
par(mfrow=c(1,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(xi.resid~type, data=z, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(xi~~of~~ARMA~Residuals), side=2, line=1.5)
dev.off()



