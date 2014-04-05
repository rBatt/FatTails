
# ===============================================
# = Function for plotting taxa (no longer used) =
# ===============================================
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







# ==========================================
# = Comparing ARMA, Norm, Log-Norm, to GEV =
# ==========================================
# Calculate the densities for those 3 distributions
dens3 <- function(x){
	df <- -1
	dt <- 11
	d1 <- matrix(unlist(density(log10(x[,"Level2_time"]), na.rm=TRUE, from=df, to=dt)[c("x","y")]), ncol=2)
	d2 <- density(log10(x[,"Level2_normTime"]), na.rm=TRUE, from=df, to=dt)$y
	d3 <- density(log10(x[,"Level2_logNormTime"]), na.rm=TRUE, from=df, to=dt)$y
	matrix(c(d1[,1], d1[,2]/max(d1[,2]), d2/max(d2), d3/max(d3)), ncol=4, dimnames=list(NULL, c("time","gev","norm","log")))
}


# Make a graph showing the densities
pDens <- function(x, go=TRUE, ...){
	cLine <- rainbow(n=3, v=0.8, s=1)
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)
	myLwd <- 1.25
	tD <- dens3(x)
	tD <- rbind(matrix(c(min(tD[,1]-0.01), 0, 0, 0),ncol=4), tD, matrix(c(max(tD[,1])+0.01, 0, 0, 0), ncol=4))
	plot(tD[,"time"], tD[,"gev"], pch=NA, ylim=c(-0.035,1), xlim=range(tD[,"time"]), ylab="", lwd=myLwd,  ...)
	if(go){
		polygon(tD[,"time"], tD[,"gev"], col=cFill[1], border=cLine[1], lwd=myLwd)
		polygon(tD[,"time"], tD[,"norm"], col=cFill[2], border=cLine[2], lwd=myLwd)
		polygon(tD[,"time"], tD[,"log"], col=cFill[3], border=cLine[3], lwd=myLwd)

		centrals <- tD[apply(tD[,-1], 2, which.max),1]
		centrals <- tD[apply(tD[,-1], 2, which.max),1]

		if(max(tail(tD[,"norm"]))>0.5){
			nID <- colnames(tD)[-1]!="norm"
		}else{
			nID <- rep(TRUE,3)
		}
		points(centrals[nID], rep(-0.0375,sum(nID)), pch=21, bg=cFill[nID], col=cLine[nID])
	}

}



# =========================================================================
# = Generic function to create colored density plots from vectors of data =
# =========================================================================
colorDens <- function(vals=NULL, cols=NULL, revxy=FALSE, mu=NULL, sig=NULL, limX=NULL, ...){
	if(is.null(vals) & (is.null(mu)|is.null(sig))) stop("Must provide vals or (mu and sig), but not both.")
	if(any(sig<=0)) stop("Need positive sig")
	if(!is.null(mu) | !is.null(sig)){
		lsig <- length(sig)
		lmu <- length(mu)
		
		if( any(c(lsig,lmu)<2) | (lsig!=lmu)){
			stop("mu and sig must be of same length, and both must have length >= 2")
		}
	}
	
	if(is.null(vals)){
		N <- length(mu)
		vals0 <- rnorm(n=100*N, mean=mu, sd=sig)
		vals <- matrix(vals0, ncol=N, byrow=TRUE)	
	}else{
		if(!is.null(dim(vals))){
			N <- dim(vals)[2]
		}else{
			N <- length(vals)
		}
	}
	if(is.null(limX)){
		limX <- range(vals)
	}
	
	if(!is.null(dim(vals))){
		dens <- apply(vals, 2, function(x)density(x, from=limX[1], to=limX[2])[c("x","y")])
	}
	if(is.list(vals)){
		dens <- lapply(vals, function(x)density(x, from=limX[1], to=limX[2])[c("x","y")])
	}
	
	dX <- lapply(dens, function(z)z$x)
	dY <- lapply(dens, function(z)z$y)
	limY <- range(dY)
	
	if(is.null(cols)){
		cLine <- rainbow(n=N)
	}else{
		cLine <- rgb(t(col2rgb(cols)), maxColorValue=255) # color for the line
	}
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=35, maxColorValue=255) # color for the fill
	
	# Need to be sure that the smallest and largest densities are 0 so that the bottom border of polygons are at 0 line
	xF <- 0.0 * diff(limX) # a "factor" by which to extend the range of X
	xA <- limX + c(-1,1)*xF # "add" this "adjustment" to the start and end of the dX

	if(revxy){
		plot(c(0,dY[[1]],0), c(xA[1],dX[[1]],xA[2]), type="l", col=cLine[1], xlab="", ylab="", ylim=limX, xlim=limY, lwd=2, ...)
		polygon(c(0,dY[[1]],0), c(xA[1],dX[[1]],xA[2]), col=cFill[1], border=NA)
		for(i in 2:N){
			polygon(c(0,dY[[i]],0), c(xA[1],dX[[i]],xA[2]), col=cFill[i], border=cLine[i], lwd=2)
		}
	}else{
		plot(c(xA[1],dX,xA[2]), c(0,dY[[1]],0), type="l", col=cLine[1], xlab="", ylab="", xlim=limX, ylim=limY, lwd=2, ...)
		polygon(c(xA[1],dX,xA[2]), c(0,dY[[1]],0), col=cFill[1], border=NA)
		for(i in 2:N){
			polygon(c(xA[1],dX,xA[2]), c(0,dY[[i]],0), col=cFill[i], border=cLine[i], lwd=2)
		}
	}
}

# ======================================================================================
# = Created a colored polygon given a vector of quantiles and values (e.g., densities) =
# ======================================================================================
# This is a generic version of colorDens
colorPoly <- function(quants, dents, cols=NULL, revxy=FALSE, ...){

	if(!is.null(dim(dents))){
		N <- dim(dents)[2]
	}else{
		N <- length(dents)
	}
	limX <- range(quants)
	
	dX <- quants
	dY <- dents
	limY <- range(dY)
	
	if(is.null(cols)){
		cLine <- rainbow(n=N)
	}else{
		cLine <- rgb(t(col2rgb(cols)), maxColorValue=255) # color for the line
	}
	cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=35, maxColorValue=255) # color for the fill
	
	# Need to be sure that the smallest and largest densities are 0 so that the bottom border of polygons are at 0 line
	xF <- 0.01 * diff(limX) # a "factor" by which to extend the range of X
	xA <- limX + c(-1,1)*xF # "add" this "adjustment" to the start and end of the dX

	if(revxy){
		plot(c(0,dY[[1]],0), c(xA[1],dX,xA[2]), type="l", col=cLine[1], xlab="", ylab="", ylim=limX, xlim=limY, lwd=2, ...)
		polygon(c(0,dY[[1]],0), c(xA[1],dX,xA[2]), col=cFill[1], border=NA)
		for(i in 2:N){
			polygon(c(0,dY[[i]],0), c(xA[1],dX,xA[2]), col=cFill[i], border=cLine[i], lwd=2)
		}
	}else{
		plot(c(xA[1],dX,xA[2]), c(0,dY[[1]],0), type="l", col=cLine[1], xlab="", ylab="", xlim=limX, ylim=limY, lwd=2, ...)
		polygon(c(xA[1],dX,xA[2]), c(0,dY[[1]],0), col=cFill[1], border=NA)
		for(i in 2:N){
			polygon(c(xA[1],dX,xA[2]), c(0,dY[[i]],0), col=cFill[i], border=cLine[i], lwd=2)
		}
	}
}
