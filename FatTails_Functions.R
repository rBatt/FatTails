#Ryan Batt
#24-April-2013
#Modified from Katz et al. 2005, Coles 2001
	#changed the way that initial values were guessed (they had used the equations for Xi=0, I used equations for Xi < 1 b/c our guess for Xi is 0.1)
	#I changed the likelihood penalty for having y<=0 or sc<=0; originally 10^6.  Changed to increase the penalty as y or sc became more negative.
	# Note to self: use model.matrix(), to get the ydat matrix --- very helpful for turning a column of factors (w/ n levels) into a matrix with n columns that contains 0 or 1
# _v0.2 (02-May-2013): Added names to the estimates in the output

#FatTail_Functions_v6.R (07-Sept-2013): Skipped version 5.  Added in mean and sd into the lvl function.
#_v7 (09-Dec-2013) I noticed some problems with returning NA from the lvl() function, because in the stationary2 function problems with log(0) are not accounted for. This was especially problematic for the zooplankton data, which for some reason contained a lot of 0's that seem like they should have been NA's.
	
gev.fit <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, show = FALSE, method = "Nelder-Mead", maxit = 10000, ...){

    z <- list() #these are going to be aspects of the output, more later
    npmu <- length(mul) + 1 #number of parameters to be fitted for mu.  if mul is NULL, then length(mul) == 0, and npmu is 1
    npsc <- length(sigl) + 1 #number of params for sigma
    npsh <- length(shl) + 1 # n param for shape
    z$trans <- FALSE #transformation?

	Xi_init <- 0.1
	g1 <- gamma(1+Xi_init)
	g2 <- gamma(1+2*Xi_init)
	sigma_init <- sqrt((var(xdat)/(g2-g1^2))*Xi_init^2)
	mu_init <- mean(xdat) + sigma_init/Xi_init - sigma_init/(Xi_init*g1)
	
    # in2 <- sqrt(6 * var(xdat))/pi #the initial guess for sigma to be used in the optimization
    # in1 <- mean(xdat) - 0.57722 * in2 #initial guess for the location

    if(is.null(mul)){
        mumat <- as.matrix(rep(1, length(xdat))) #if mul is NULL, then only 1 parameter, the intercept, is needed for mu, and the predictor variable is 1 
        muinit <- mu_init
    }else{
        z$trans <- TRUE #if we do have non-NULL mul, then the transformation (transfer function?) is TRUE
        mumat <- cbind(rep(1, length(xdat)), ydat[, mul]) #the predictor variable matrix will include a column of 1's (for the intercept), and then all of the columns in ydat that correspond to other desired predictor variables.  E.g., mul could be equal to "Year", in which case our predictor matrix, mumat, would have a column of 1's of length xdat, and a column of the Years, also of length xdat.  The model would be similar to Y[i] = 1*B0 + Year[i]*B1, or Y = mumat%*%params, where Y is actually xdat, xdat has N rows and 1 column, mumat has N rows and 2 columns, and params is a matrix with 2 rows and 1 column
        muinit <- c(mu_init, rep(0, length(mul))) #now we have to make more guesses for the initial values --- we calculate the intercept as the "default" estimate of the parameter, and 0 (or no effect) as the default estimate of the other parameters
    }
    if (is.null(sigl)){ #same idea as for mul, and subsequent estiamtes of initials and matrices
        sigmat <- as.matrix(rep(1, length(xdat)))
        siginit <- sigma_init
    }else{
        z$trans <- TRUE
        sigmat <- cbind(rep(1, length(xdat)), ydat[, sigl])
        siginit <- c(sigma_init, rep(0, length(sigl)))
    }

    if (is.null(shl)){ #same idea as for mu and sigma
        shmat <- as.matrix(rep(1, length(xdat)))
        shinit <- 0.1
    }else{
        z$trans <- TRUE
        shmat <- cbind(rep(1, length(xdat)), ydat[, shl])
        shinit <- c(0.1, rep(0, length(shl)))
    }

    z$model <- list(mul, sigl, shl) #this is a list of 3 elements, where each element is a character vector whose elements in turn specify a column (name) of ydat to be used as a predictor variables for the respective GEV parameter
    z$link <- deparse(substitute(c(mulink, siglink, shlink))) #turns the link functions used into character strings
    init <- c(muinit, siginit, shinit) #Consider log-transforming the siginit, and then inside gev.like exponentiating. This would avoid the problem of sigma being <= 0.
		#however, problem with transforming is that some of the sigma parameters can be negative, so long as none of the sc are negative. 
		#essentially, the sum of some parameters related to sigma can't be negative, but any one of them *may* be allowed negative
		#if we leave the 10^6 solution as-is, the optimization can get "lost" on a massively-elevated NLL plane
		#perhaps a solution would be to slope the plane such that the NLL gets worse and worse as sc gets further below 0. same idea for Y.  
		#an initial jump in NLL when either sc or y reach 0, and then a steady increase afterwards.  maybe I'll try this.

    gev.lik <- function(a){
        mu <- mulink(mumat %*% (a[1:npmu]))#unpack the mu parameters from the a vector
        sc <- siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)])) #unpack the sigma parameters from the a vector
        xi <- shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])) #unpack the shape parameters from the a vector
        y <- (xdat - mu)/sc #this line is part of t(x)
        y <- 1 + xi * y #this line is more of t(x), and is now only missing the ^(-1/Xi) part
        if (any(y <= 0) || any(sc <= 0)){
			if(any(y <= 0)){ypenal <- (1/exp(y[which(y<=0)]))}else{ypenal <- 0}
			if(any(sc <= 0)){scpenal <- (1/exp(sc[which(sc<=0)]))}else{scpenal <- 0}
			BadReturn <- 10^5 * (mean(ypenal) + mean(scpenal))
			return(BadReturn)
            # return(10^6) #if you get an impossible value for y or sc, just return a huge nll, which will discourage the optimizer from picking the current parameters
			#also, this is a shortcut for the support that X is an element of (mu-sigma/Xi to INF) when Xi >0, and X is an element of (-INF to mu-sigma/Xi) when Xi < 0
				#there is only support for those values of X because t_x must be defined... 
				#in the case where Xi!=0, (1+((x-mu)/sigma)*xi), which is defined as "y" above, ends up being in the denominator when Y^(-1/xi)
				# Y^(-1/xi) is t(x)
				#while Y can't be zero, b/c denom, it can't be negative either b/c we take the log of y in the likelihood equation
		}
        sum(log(sc)) + sum(y^(-1/xi)) + sum(log(y)*(1/xi + 1)) #negative log likelihood (from the pdf of the gev)
			#sum(log(sc)) is the -log() of 1/sigma
			#sum(y^(-1/xi)) is the -log() of exp(-t(x))
			#sum(log(y) * (1/xi + 1)) is the -log() of t(x)^(xi+1) 
				#because had (y^-1/xi)^*(xi+1)  
				#==> y^(-1 + -1/xi) 
				#==> log(y) * -(1 + 1/xi) 
				#==> -log(y)*(1/xi + 1) (and then take neg for nll)
    }

    x <- optim(init, gev.lik, hessian = TRUE, method = method, control = list(maxit = maxit, ...))
    z$conv <- x$convergence
    mu <- mulink(mumat %*% (x$par[1:npmu]))
    sc <- siglink(sigmat %*% (x$par[seq(npmu + 1, length = npsc)]))
    xi <- shlink(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)]))
    z$nllh <- x$value
    z$data <- xdat

    if (z$trans){
        z$data <- -log(as.vector((1 + (xi * (xdat - mu))/sc)^(-1/xi)))
    }
    z$mle <- x$par
	names(z$mle) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))
    z$cov <- solve(x$hessian)#I should use the hessian() function in package numDeriv to calculate this, and compare it to the hessian from the optimization algorithm.
    z$se <- sqrt(diag(z$cov))
	names(z$se) <- c(paste("mu",c(0,mul),sep="_"), paste("sig",c(0,sigl),sep="_"), paste("sh",c(0,shl),sep="_"))
    z$vals <- cbind(mu, sc, xi)
    if(show){
        if(z$trans){
            print(z[c(2, 3, 4)])
        }else{print(z[4])}
        if(!z$conv){
            print(z[c(5, 7, 9)])
		}
    }
    invisible(z)
}



# ====================
# = Extrema Functions =
# ====================
# ================================
# = Method of independent storms =
# ================================
MIS <- function(x, Thresh=NULL, quant=0.9){
	if(is.null(Thresh)){
		Thresh <- quantile(x, quant)
	}
	Above <- as.integer(x>=Thresh)
	udF <- ifelse(Above[1]==1, 1, 0)
	udL <- ifelse(rev(Above)[1]==1, -1, 0)
	UpDown <- c(udF,diff(Above),udL)
	StartStop <- matrix(sort(c(which(UpDown==1), (which(UpDown==-1)-1))), ncol=2, byrow=TRUE)
	Storms <- rep(NA, nrow(StartStop))
	StormIndex <- Storms
	for(i in 1:nrow(StartStop)){
		ThisInd <- StartStop[i,1]:StartStop[i,2]
		StormIndex[i] <- which.max(x[ThisInd]) + ThisInd[1] - 1
		Storms[i] <- x[StormIndex[i]]

	}
	return(matrix(c(StormIndex, Storms), ncol=2, dimnames=list(NULL, c("Index", "Max"))))
}

# ========================================
# = Find the peak of a cycle via wavelet =
# ========================================
PeakCycle <- function(Data, SearchFrac=0.028){
	# using package "wmtsa"
	#the SearchFrac parameter just controls how much to look to either side 
	#of wavCWTPeaks()'s estimated maxima for a bigger value
	#see dRange
	Wave <- wavCWT(Data)
	WaveTree <- wavCWTTree(Wave)
	WavePeaks <- wavCWTPeaks(WaveTree, snr.min=5)
	WavePeaks_Times <- attr(WavePeaks, which="peaks")[,"iendtime"]
	
	NewPeakTimes <- c()
	dRange <- round(SearchFrac*length(Data))
	for(i in 1:length(WavePeaks_Times)){
		NewRange <- max(c(WavePeaks_Times[i]-dRange, 1)):min(c(WavePeaks_Times[i]+dRange, length(Data)))
		NewPeakTimes[i] <- which.max(Data[NewRange])+NewRange[1]-1
	}
	
	return(matrix(c(NewPeakTimes, Data[NewPeakTimes]), ncol=2, dimnames=list(NULL, c("PeakIndices", "Peaks"))))
}

# ==============================
# = Replace non-finite w/ NA's =
# ==============================
Inf2NA <- function(x) {x[which(x==-Inf | x==Inf, arr.ind=TRUE)] <- NA; x}

Stationary <- function(x, coluY, coluX){
	tx <- 0:(nrow(x)-1)
	x2 <- x[,coluX]
	Reg <- lm(x[,coluY]~tx*x2)
	xd <- residuals(Reg) #+ mean(x[,coluY])
	
	return(xd)
	
}


Stationary2 <- function(x, coluY, coluX){
	tx <- 0:(length(x)-1)
	Reg <- lm(x~tx)
	TrendPval <- summary(Reg)$coef["tx","Pr(>|t|)"]
	xd <- residuals(Reg) #+ mean(x[,coluY])
	return(xd)
}

Stationary3 <- function(x){
	tx <- 0:(length(x[,2])-1)
	Reg <- lm(x[,2]~tx)
	# TrendPval <- summary(Reg)$coef["tx","Pr(>|t|)"]
	xd <- residuals(Reg) + mean(x[,2], na.rm=TRUE)
	xNew <- x
	xNew[,2] <- xd
	return(xNew)
}

# =================
# = GEV Functions =
# =================
ReturnLevel <- function(RetQ, mu, sc, xi){
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi)) 	
}
ReturnQuant <- function(RetTime, nExts, TS_Duration){
	-TS_Duration/(RetTime*nExts) + 1
}
ReturnTime <- function(RetQ, nExts, TS_Duration){
	1/((1-RetQ)*(nExts/TS_Duration))
}

xYr_Lvl <- function(RetTime, nExts, TS_Duration, mu, sc, xi){
	RetQ <- -TS_Duration/(RetTime*nExts) + 1
	mu-(sc/xi)*(1-(-log(RetQ))^(-xi))
}



lvlX_ReturnTime <- function(lvl, a, nExts, TS_Duration){ #Added _v4
	if(any(is.na(c(lvl,a)) | is.nan(c(lvl,a)))){return(NA)}
	
	mu <- a[1]
	sc <- a[2]
	xi <- a[3]
	
	if(xi < 0 & lvl > (mu-sc/xi)){return(NA)}
	
	tx <- (1 + xi * (lvl - mu)/sc)^(-1/xi) #the t(x) from Wikipedia GEV article
	#pdf <- (1/sc)*tx^(xi+1)*exp(-tx) #this is the pdf, if I ever want to use it
	cdf <- exp(-tx) #this is the CDF, which gives the percentile for a certain return level
	1/((1-cdf)*(nExts/TS_Duration))
}


lvl1 <- function(x){
	if(sum(!is.na(x))<3){return(NA)}
	x <- Stationary2(x) + mean(x, na.rm=TRUE)
	mean(x, na.rm=TRUE)+2*sd(x, na.rm=TRUE)
}

lvl2 <- function(x, thresh=1.2){
	if(sum(!is.na(x))<3){return(NA)}
	x <- Stationary2(x) + mean(x, na.rm=TRUE)
	thresh*max(x, na.rm=TRUE)
}
reshape2 <- function(...){
	a <- reshape(...)
	row.names(a) <- NULL
	a <- a[,!names(a)=="id"]
	a
}

convNeg <- function(x){ #changing in _v8 (09-Dec-2013) to also change 0's to NA
	a <- !is.na(x) & is.finite(x) #also adding the finite stipulation, although I think that is.finite() more than covers !is.na().
	n <- length(x[a])
	if(sum(x[a]<=0)<(0.2*n)){
		xr <- x
		xr[x<=0] <- NA
		return(xr)
	}else{
		return(NA)
		# m0 <- min(x, na.rm=TRUE)
		# xr <- (exp(1)-m0) + x
		# return(xr)
	}
}

lvlmean <- function(x, Log=FALSE){
	if(sum(!is.na(x))<3){return(NA)}
	x <- Stationary2(x) + mean(x, na.rm=TRUE)
	if(Log){
		x <- Stationary2(x) + mean(x, na.rm=TRUE)
		cx <- convNeg(x)
		lx <- log(cx)
		mx <- mean(lx, na.rm=TRUE) #note that this line was returning NA in cases where x contained 0's; this has now been fixed via update to convNeg (_v8)
		return(mx)
	}
	mean(x, na.rm=TRUE)
}

lvlsd <- function(x, Log=FALSE){
	if(sum(!is.na(x))<3){return(NA)}
	x <- Stationary2(x) + mean(x, na.rm=TRUE)
	if(Log){
		x <- Stationary2(x) + mean(x, na.rm=TRUE)
		cx <- convNeg(x)
		lx <- log(cx)
		sx <- sd(lx, na.rm=TRUE)
		return(sx)
	}
	sd(x, na.rm=TRUE)
}

lvlObs <- function(x){
	sum(!is.na(x))
}

lvl <- function(x, fitby, variables, lvl2Thresh=1.2){
	# stopifnot(is.list(fitby))
	if(!is.null(fitby)){
		one0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvl1)
		two0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvl2, thresh=lvl2Thresh)
		mean0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvlmean)
		sd0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvlsd)
		obs0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvlObs)
		lmean0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvlmean, Log=TRUE)
		lsd0 <- aggregate(x[,variables], by=list("fitBy"=x[,fitby]), FUN=lvlsd, Log=TRUE)
	}else{
		one0 <- data.frame(lvl1(x[,variables])); names(one0) <- variables
		two0 <- data.frame(lvl2(x[,variables])); names(two0) <- variables
		mean0 <- data.frame(lvlmean(x[,variables])); names(mean0) <- variables
		sd0 <- data.frame(lvlsd(x[,variables])); names(sd0) <- variables
		obs0 <- data.frame(lvlObs(x[,variables])); names(obs0) <- variables
		lmean0 <- data.frame(lvlmean(x[,variables], Log=TRUE)); names(lmean0) <- variables
		lsd0 <- data.frame(lvlsd(x[,variables], Log=TRUE)); names(lsd0) <- variables

	}
	if(length(variables)==1 & !is.null(fitby)){
		names(one0)[2] <- variables
		names(two0)[2] <- variables
		names(mean0)[2] <- variables
		names(sd0)[2] <- variables
		names(obs0)[2] <- variables
		names(lmean0)[2] <- variables
		names(lsd0)[2] <- variables
	}
	
	one <- reshape2(one0, varying=variables, direction="long", v.names="Level1", times=variables, timevar="Variable")
	two <- reshape2(two0, varying=variables, direction="long", v.names="Level2", times=variables, timevar="Variable")
	mean1 <- reshape2(mean0, varying=variables, direction="long", v.names="mean", times=variables, timevar="Variable")
	sd1 <- reshape2(sd0, varying=variables, direction="long", v.names="sd", times=variables, timevar="Variable")
	obs1 <- reshape2(obs0, varying=variables, direction="long", v.names="nObs", times=variables, timevar="Variable")
	lmean1 <- reshape2(lmean0, varying=variables, direction="long", v.names="logMean", times=variables, timevar="Variable")
	lsd1 <- reshape2(lsd0, varying=variables, direction="long", v.names="logSd", times=variables, timevar="Variable")
	
	# rlvl <- merge(merge(one, two), merge(mean1, sd1))
	m1 <- merge(one, two)
	m2 <- merge(mean1, sd1)
	m3 <- merge(lmean1, lsd1)
	
	ma <- merge(m1, m2)
	mb <- merge(m3, obs1)
	
	rlvl <- merge(ma, mb)
	# rlvl <- merge_all(one, two, mean1, sd1, obs1, by=c(""))
	
	return(rlvl)
}

normTime <- function(x, Level=1){
	lvlX <- paste("Level",Level,sep="")
	RetQ <- pnorm(x[,lvlX], x[,"mean"], x[,"sd"])
	1/((1-RetQ)*(x[,"nObs"]/x[,"Duration"]))
}

lnormTime <- function(x, Level=1){
	lvlX <- paste("Level",Level,sep="")
	RetQ <- plnorm(x[,lvlX], x[,"logMean"], x[,"logSd"])
	1/((1-RetQ)*(x[,"nObs"]/x[,"Duration"]))
}

# ==========================
# = The GEV Wrapper for Me =
# ==========================
calcGEV <- function(nameVarbl, datCols=NULL, fitBy=NULL, fitForm=NULL, MUl=NULL, SIGl=NULL, SHl=NULL){
	if(is.null(fitBy)&is.null(fitForm) || !is.null(fitBy)&!is.null(fitForm)){stop("One & only one of fitBy and fitForm should be NULL")}
	if(is.null(datCols)){datCols <- get(paste("All",nameVarbl,sep=""))}
	Tempo_gev <- list()
	AllVarbl <- datCols
	for(j in 1:length(AllVarbl)){
		VarblName <- AllVarbl[j]
		Varbl_Ext <- get(paste(nameVarbl,"Ext",sep="_"))
		
		# ================================================================
		# = Set number of for() iterations depending on fitBy or fitForm =
		# ================================================================
		if(is.null(fitBy)){
			LoopThru <- "Formula" #need to change the formula arugment to just be responses, then paste it together with datCols[i] (i didn't write this thinking about multiple responses)
		}else{
			LoopThru <- unique(Varbl_Ext[,fitBy])
		}
		
		
		for(i in seq_along(LoopThru)){ # 1 iteration for formula, or number of factor levels for fitBy
			# ==============================================
			# = Fit GEV for each level of a factor (fitBy) =
			# ==============================================
			if(!is.null(fitBy)){
				fitByID <- as.character(LoopThru[i])
				ThisFitBy <- which(Varbl_Ext[,fitBy]==fitByID & !is.na(Varbl_Ext[,VarblName]))
				if(length(ThisFitBy)<5){next}
				ThisVarbl <- Varbl_Ext[ThisFitBy,VarblName]
				ThisVarbl <- as.vector(Stationary2(ThisVarbl)) + mean(ThisVarbl)
				Nobs <- length(ThisVarbl)
				dRange <- range(Varbl_Ext[ThisFitBy,"year4"])
				dDuration <- diff(dRange) + 1
				
				Tempo_gev <- gev.fit(ThisVarbl)[c("conv","nllh","mle","se")] # Fit GEV to this factor level
			}
			
			
			# ============================================
			# = Fit GEV according to a formula (fitForm) =
			# ============================================
			if(!is.null(fitForm)){
				ModForm <- with(Varbl_Ext, formula(fitForm))
				AllVars <- all.vars(ModForm)
				RespName <- head(AllVars, 1)
				PredNames <- tail(AllVars, -1)
				TempoFullMat <- model.matrix(ModForm)
				colnames(TempoFullMat) <- AllVars
				ThisVarbl <- as.vector(Varbl_Ext[,RespName])
				ThisVarbl_NoNA <- which(!is.na(ThisVarbl))
				ThisVarbl <- ThisVarbl[ThisVarbl_NoNA]
				TempoYdat <- matrix(TempoFullMat[, PredNames], ncol=length(PredNames), dimnames=list(NULL, PredNames)) #i think model.matrix automatically removes the NA rows
				Nobs <- length(ThisVarbl)
				dRange <- range(Varbl_Ext[,"year4"])
				dDuration <- diff(dRange) + 1
				
				Tempo_gev <- gev.fit(ThisVarbl, ydat=TempoYdat, mul=MUl, sigl=SIGl, shl=SHl)[c("conv","nllh","mle","se")] # Fit GEV according to formula
			}
			
			# =======================================================
			# = Calculate t-statistic & p value for shape parameter =
			# =======================================================
			tstat <- Tempo_gev$mle/Tempo_gev$se
			tstat2 <- (Tempo_gev$mle["sh_0"]-ifelse(tstat[3]>0, 0.5, -0.5))/(Tempo_gev$se["sh_0"]) #testing to see if the shape parameter is greater than 0.5
		
			pvalue <- c()
			for(k in 1:length(tstat)){
				pvalue[k] <- pt(tstat[k], df=Nobs-1, lower.tail=(tstat[k]<0))
			}
			pvalue2 <-  pt(tstat2, df=Nobs-1, lower.tail=(tstat2<0))#the p-value to see if shape is greater than 0.5
			Tempo_gev$Dates <- c("Nobs"=Nobs, "Dates"=dRange, "Duration"=dDuration)
			Tempo_gev$Pvalues <- c(pvalue, pvalue2)
			names(Tempo_gev$Pvalues) <- c(names(Tempo_gev$mle), "sh_0.5")
		
			TempoNames <- c(names(Tempo_gev$mle), paste("se_",names(Tempo_gev$se), sep=""), paste("p_",names(Tempo_gev$Pvalues), sep=""))
		
			# ==========================================================
			# = Prepare function output (updated for each fitBy level) =
			# ==========================================================
			if(all(c(j,i)==1)){
				tTempo_Params <- matrix(c(Tempo_gev$mle, Tempo_gev$se, Tempo_gev$Pvalues), nrow=1, dimnames=list(NULL, TempoNames)) 
				Tempo_Params <- data.frame("fitBy"=LoopThru[i], "Variable"=VarblName, "Duration"=dDuration, "N"=Nobs, tTempo_Params)
				# LoopThru, "Variable"=VarblName, 
			}else{
				tTempo_Params <- matrix(c(Tempo_gev$mle, Tempo_gev$se, Tempo_gev$Pvalues), nrow=1, dimnames=list(NULL, TempoNames))
				Tempo_Params <- rbind(Tempo_Params, data.frame("fitBy"=LoopThru[i], "Variable"=VarblName, "Duration"=dDuration, "N"=Nobs, tTempo_Params))
				# LoopThru, "Variable"=VarblName, 
			}
		}
	}
	return(list(Tempo_Params, Tempo_gev))
}


SignifShape <- function(x, ...){
	cneg <- as.integer(x[,"sh_0"]<0&x[,"p_sh_0"]<0.05)*-1
	c1 <- as.integer(x[,"sh_0"]>0&x[,"p_sh_0"]<0.05)
	c2 <- as.integer(x[,"sh_0"]>0.5&x[,"p_sh_0.5"]<0.05)
	cmat <- matrix(c(cneg, c1, c2), ncol=3)
	rs <- rowSums(cmat)
	plot(x[,"sh_0"], col=c("blue","black","red", "red")[rs+2], pch=ifelse((rs+2)>3, 10, 21), ylab="Shape", ...)
	invisible(rs)
}

stable.fit <- function(x){	
	sFit <- stableFit(x, type="mle", doplot=FALSE)
	sFit@fit$estimate
}

stableTime <- function(lvl, a, nExts, TS_Duration){ #Added _v4
	if(any(is.na(c(lvl,a)) | is.nan(c(lvl,a)))){return(NA)}
	
	alpha <- a[1]
	beta <- a[2]
	cc <- a[3]
	mu <- a[4]
	
	# if(xi < 0 & lvl > (mu-sc/xi)){return(NA)}
	# tx <- (1 + xi * (lvl - mu)/sc)^(-1/xi) #the t(x) from Wikipedia GEV article
	#pdf <- (1/sc)*tx^(xi+1)*exp(-tx) #this is the pdf, if I ever want to use it
	cdf <- pstable(lvl, alpha=alpha, beta=beta, gamma=cc, delta=mu) #this is the CDF, which gives the percentile for a certain return level
	1/((1-cdf)*(nExts/TS_Duration))
}

# Example
# x= rnorm(100)
# stableTime(3, a=stable.fit(x), nExts=100, TS_Duration=100)



