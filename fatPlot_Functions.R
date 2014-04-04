
plotTax <- function(D, V, Stat, taxCol, legendTitle=NULL, ...){
	
	D <- D[D[,"variable"]==V,]
	pLvls <- rev(levels(factor(D[,"taxLvl"])))
	# pLvls <- pLvls[!pLvls%in%"Community"]
	nLvl <- length(pLvls)

	nD <- names(D)
	lastLvl <- rev(pLvls)[1]
	llC <- which(nD==lastLvl)

	plot(c(1,length(pLvls)), range(D[,Stat], na.rm=TRUE), type="n", xlab="", ylab="", xaxt="n", ...)
	axis(side=1, labels=FALSE)

	Ds <- D[!is.na(D[,"Species"]),]
	# Dl <- !is.na(D[,names(D)%in%c(c("Phylum", "Class", "Order", "Family", "Genus", "Species"))])
	# startLvl <- rev(pLvls)[rowSums(Dl)]
	colFac <- as.numeric(as.factor(Ds[,taxCol]))
	Cols <- rainbow(n=max(colFac, na.rm=TRUE), alpha=0.25)
	Cols2 <- rainbow(n=max(colFac, na.rm=TRUE), alpha=1)
	holdAll <- matrix(NA, nrow=nrow(Ds), ncol=length(pLvls))
	for(l in 1:nrow(Ds)){
		shapes <- c() #Ds[l,"sh_0"]
		xLvl <- c() #rep(1, length(shapes))
		for(sh in 1:length(pLvls)){
			tlog <- as.character(D[,"taxID"])==as.character(Ds[l,pLvls[sh]])
			llog <- as.character(D[,"location"])==as.character(Ds[l,"location"])
			addShape <- D[tlog&llog,"sh_0"]
			# xLvl <- c(xLvl, rep(sh, length(addShape)))
			shapes <- c(shapes, addShape)
		}
		holdAll[l,] <- shapes
		tcol <- Cols[colFac[l]]
		lines(1:length(shapes), shapes, col=tcol)
		if(l==nrow(Ds)){
			points(1:length(shapes), colMeans(holdAll, na.rm=TRUE), pch=20, cex=3)
			legLabs <- unique(as.factor(Ds[,taxCol]))[unique(colFac)]
			nColumnLeg <- max(min(3, floor(length(legLabs)/3)),1)
			legend("topright", legend=legLabs, col=Cols2[unique(colFac)], lty="solid", title=legendTitle, ncol=nColumnLeg)
		}	
	}
}

# V = "density"
# D = zoop.gev
# Stat = "sh_0"
# taxCol = "Phylum"
# plotTax(D=zoop.gev, V="density", Stat="sh_0", taxCol="Phylum")


