

load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")

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





