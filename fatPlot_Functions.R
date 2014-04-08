
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







# ================================================
# = Beanplot w/ configurable weights for density =
# ================================================
bp2 <- function (..., bw = "SJ-dpi", kernel = "gaussian", cut = 3, cutmin = -Inf, 
    cutmax = Inf, grownage = 10, what = c(TRUE, TRUE, TRUE, TRUE), 
    add = FALSE, col, axes = TRUE, log = "auto", handlelog = NA, 
    ll = 0.16, wd = NA, maxwidth = 0.8, maxstripline = 0.96, 
    method = "stack", names, overallline = "mean", beanlines = overallline, 
    horizontal = FALSE, side = "no", jitter = NULL, beanlinewd = 2, 
    frame.plot = axes, border = NULL, innerborder = NA, at = NULL, 
    boxwex = 1, ylim = NULL, xlim = NULL, show.names = NA, bweights=NULL) 
{
    mdensityxy <- function(x) {
        if (length(x) > 0) {
            from <- max(cutmin, (min(mlog(x)) - cut * bw))
            to <- min(cutmax, max(mlog(x)) + cut * bw)
            density(mlog(x), bw = bw, kernel = kernel, from = from, 
                to = to, weights=bweights)[c("x", "y")]
        }
        else list(x = numeric(), y = numeric())
    }
    args <- match.call()
    mcall <- as.list(args)
    method <- pmatch(method, c("overplot", "stack", "jitter"))
    if (is.na(method) || method == 0) 
        stop("invalid plotting method")
    beanlines <- pmatch(beanlines, c("mean", "median", "quantiles"))
    if (is.na(beanlines) || beanlines == 0) 
        stop("invalid beanlines")
    overallline <- pmatch(overallline, c("mean", "median"))
    if (is.na(overallline) || overallline == 0) 
        stop("invalid overallline")
    side <- pmatch(side, c("no", "first", "second", "both"))
    if (is.na(side) || side == 0) 
        stop("invalid side")
    groups <- getgroupsfromarguments(args)
    groups <- lapply(groups, na.omit)
    n <- length(groups)
    displayn <- if (side == 4) 
        ceiling(n/2)
    else n
    if (n == 0) 
        stop("no data found to beanplot")
    if (missing(names)) {
        if (is.null(base::names(groups))) 
            attr(groups, "names") = 1:displayn
        names <- base::names(groups)
    }
    else {
        attr(groups, "names") <- names
        if (is.na(show.names)) 
            show.names <- TRUE
    }
    if (is.null(at)) {
        at <- 1:displayn
    }
    if ((side == 4) && (length(names) > length(at))) {
        for (i in 1:length(at)) {
            names[i] <- makecombinedname(names[i * 2 - 1], names[i * 
                2])
        }
        length(names) <- length(at)
    }
    combinedpolygons <- ((side == 4) && (length(border) < 2) && 
        (n%%2 == 0))
    if (missing(col)) 
        col <- par("fg")
    if (!is.list(col)) 
        col <- list(col)
    else combinedpolygons <- FALSE
    col <- lapply(col, fixcolorvector)
    col <- rep(col, length.out = n)
    border <- rep(border, length.out = n)
    if (!add && log == "auto") {
        if (seemslog(groups)) {
            log <- "y"
            message("log=\"y\" selected")
        }
        else log <- ""
    }
    if (is.na(handlelog)) 
        if (add && ((horizontal & par()$xlog) || (!horizontal & 
            par()$ylog))) 
            handlelog <- TRUE
        else if (!add && (log != "")) 
            handlelog <- TRUE
        else handlelog <- FALSE
    if (handlelog) {
        mlog <- base::log
        mexp <- base::exp
    }
    else {
        mlog <- function(x) {
            x
        }
        mexp <- mlog
    }
    if (!is.numeric(bw)) {
        bw <- mean(sapply(groups, function(x) {
            ifelse(length(x) > 1, density(mlog(x), kernel = kernel, 
                bw = bw, weights=bweights)$bw, NA)
        }), na.rm = TRUE)
        if (is.nan(bw)) 
            bw <- 0.5
    }
    dens <- sapply(groups, mdensityxy)
    for (i in 1:n) dens[["y", i]] <- dens[["y", i]] * min(1, 
        length(groups[[i]])/grownage)
    if (is.na(wd)) 
        wd <- maxwidth/max(unlist(dens["y", ]))
    wd2 <- wd * boxwex/2
    axespars <- lapply(mcall[base::names(mcall) %in% c("xaxt", 
        "yaxt", "las", "cex.axis", "col.axis", "format", "tick", 
        "xaxp", "yaxp")], eval, parent.frame())
    if (!add) {
        if (!is.numeric(xlim)) {
            if (side == 2) 
                xlim <- c(0, displayn)
            else if (side == 3) 
                xlim <- c(1, displayn + 1)
            else xlim <- c(0.5, displayn + 0.5)
        }
        if (!is.numeric(ylim)) 
            ylim <- range(groups, mexp(unlist(dens["x", ])))
        plot.new()
        windowpars <- lapply(mcall[base::names(mcall) %in% c("yaxs", 
            "xaxs")], eval)
        if (horizontal) {
            names(windowpars)[names(windowpars) %in% c("xaxs", 
                "yaxs")] <- rev(names(windowpars)[names(windowpars) %in% 
                c("xaxs", "yaxs")])
            if (log == "y") 
                log <- "x"
            do.call("plot.window", c(list(xlim = ylim, ylim = xlim, 
                log = log), windowpars))
        }
        else {
            do.call("plot.window", c(list(xlim = xlim, ylim = ylim, 
                log = log), windowpars))
        }
        if (frame.plot) 
            box()
        if (axes) 
            do.call("axis", c(list(side = 2 - horizontal), axespars))
    }
    if (axes) {
        if (is.na(show.names)) 
            show.names <- (n > 1)
        if (show.names) 
            do.call("axis", c(list(1 + horizontal, at = at, labels = names), 
                axespars))
    }
    if (what[1]) {
        if (overallline == 2) 
            val <- mexp(median(mlog(unlist(groups))))
        else val <- mexp(mean(mlog(unlist(groups))))
        if (horizontal) 
            abline(v = val, lty = 3)
        else abline(h = val, lty = 3)
    }
    if (what[2]) {
        beanplotpolyshapes(side, dens, at, wd2, combinedpolygons, 
            displayn, n, col, border, horizontal, mlog, mexp)
    }
    if (what[3]) {
        beanplotbeanlines(groups, side, beanlines, beanlinewd, 
            at, boxwex, n, col, horizontal, mlog, mexp)
    }
    if (what[4]) {
        beanplotscatters(groups, side, method, jitter, dens, 
            at, wd2, boxwex, n, ll, maxstripline, col, horizontal, 
            mlog, mexp)
    }
    if (any(!is.na(innerborder))) {
        beanplotinnerborders(innerborder, at, dens, side, displayn, 
            n, horizontal, mexp)
    }
    titlepars <- lapply(mcall[base::names(mcall) %in% c("main", 
        "sub", "xlab", "ylab", "cex.main", "col.main", "cex.lab", 
        "col.lab", "cex.sub", "col.sub")], eval, parent.frame())
    do.call("title", titlepars)
    invisible(list(bw = bw, wd = wd))
}


getgroupsfromarguments2 <- function (args = match.call(sys.function(sys.parent()), sys.call(sys.parent())), 
    envir = parent.frame(2)) 
{
    nextargpos <- function(name, pos) {
        if (any(pos == length(args))) 
            return(pos)
        posnext <- match(name, base::names(args[max(pos + 1, 
            3):length(args)])) + max(pos + 1, 3) - 1
        if (!is.na(posnext)) 
            pos <- posnext
        pos
    }
    if (is.null(base::names(args))) 
        vars <- 1:length(args)
    else vars <- c(1, which(base::names(args)[2:length(args)] %in% 
        c("formula", "x", "data", "")) + 1)
    if (length(vars) < 2) 
        return(list())
    options <- which(base::names(args) %in% c("subset", "na.action", 
        "drop.unused.levels", "xlev"))
    args <- args[c(vars, options)]
    args[[1]] <- quote(model.frame)
    hashad <- rep(FALSE, length(vars))
    groups <- list()
    notnamed <- 0
    subsetno <- numeric()
    naactno <- numeric()
    dulno <- numeric()
    xlevno <- numeric()
    datano <- numeric()
    argsvals <- lapply(as.list(args[2:length(vars)]), eval, envir)
    islists <- lapply(argsvals, function(x) {
        is.list(x) || is.null(x)
    })
    for (i in 2:length(vars)) {
        if (hashad[i]) 
            next
        x <- argsvals[[i - 1]]
        if (inherits(x, "formula")) {
            datanonext <- match(TRUE, islists[max(datano, i):length(islists)]) + 
                max(datano, i)
            if (!is.na(datanonext)) {
                hashad[datanonext] <- TRUE
                datano <- datanonext
            }
            subsetno <- nextargpos("subset", subsetno)
            naactno <- nextargpos("na.action", naactno)
            dulno <- nextargpos("drop.unused.levels", dulno)
            xlevno <- nextargpos("xlev", xlevno)
            attr(args, "names")[i] <- "formula"
            m <- args[c(1, i, datano, subsetno, naactno, dulno, 
                xlevno)]
            mf <- eval(m, envir)
            response <- attr(attr(mf, "terms"), "response")
            groups <- c(groups, split(mf[[response]], mf[-response]))
        }
        else if (is.list(x)) {
            groups <- c(groups, x)
        }
        else {
            x <- list(x)
            notnamed <- notnamed + 1
            attr(x, "names") <- notnamed
            groups <- c(groups, x)
        }
    }
    groups
}
