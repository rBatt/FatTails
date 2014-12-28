

load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/fatPlot_Functions.R")


png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/ED_Fig6_xiTS_lambda.png", res=150, units="in", height=3.5, width=3.5)
# dev.new(height=3.5, width=3.5)
dscat(z[,"lambda"], z[,"xi2"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# abline(lm(xi~xi.resid, data=z[,]), lty="dashed", lwd=1)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(abs(abs({lambda}))), side=1, line=1.5)
dev.off()
