
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

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