


load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

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