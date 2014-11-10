

library(beanplot)
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
# source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

bLine <- rainbow(n=5, v=0.8, s=1)
bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
beanCol <- list(c(bFill[1]),
				c(bFill[2]),
				c(bFill[3]),
				c(bFill[4])
				)

# ======================
# = Signif Xi Beanplot =
# ======================
# dev.new(width=3.5, height=3.5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Other/signifXi_beanplot.png", res=150, units="in", height=3.5, width=3.5)
par(mfrow=c(1,1), mar=c(2,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(shape.sig~Type, data=data.fat, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5, bw=0.05, what=c(1,1,1,1), maxstripline=0.05)
axis(side=2)
axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
mtext(bquote(Critical~~Value~~'for'~~xi), side=2, line=1.5)
dev.off()
