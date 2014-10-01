
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/Diaptomus.eg.plot.RData")


png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Diaptomus.eg.arma.xi.png", res=150, units="in", height=6, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, cex=1, mgp=c(1.5, 0.4, 0), tcl=-0.3, family="Times")

plot(X, type="l",xlab="",ylab="X")
points(Y ~ ii)

plot(XX, type="l",xlab="Time",ylab="Residual")
points(YY ~ iii)

dev.off()