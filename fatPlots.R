
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

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
wtLab <- pretty(log10(data.fat[,"Level2_time"]))
axis(side=2, at=wtLab, labels=10^wtLab)
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




dev.new()

boxplot(sh_0~taxLvl, data=density.gev)

dev.new()
boxplot(sh_0~taxID, data=density.gev[density.gev[,"taxLvl"]=="Order",])

dev.new()

boxplot(sh_0~taxID, data=bioComm.gev)

dev.new()
boxplot(sh_0~Type, data.fat, outline=FALSE, ylab=expression(xi))


dev.new()
boxplot(log10(Level2_time)~Type, data.fat, outline=FALSE, ylab=bquote(log[10]*Waiting~Time))



