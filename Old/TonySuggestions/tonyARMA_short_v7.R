library("DEoptim")
library("GenSA")

# _v5 (13-Jan-2014) Making some corrections to the computation of the sigmas after conversation w/ Tony. Apparently Vstationary[1,1] is equivalent to sig2Inf/sig2E, so I was effectively using sig2Inf/sig2E/sig2E in place of sig2Inf/sig2E.  Tony also had some suggestions about detrending and log-transforming data --- I might look into the impact of these transformations on a later date.
# _v6 (28-Jan-2013) have getSE return the residuals as well
# _v7 (30-Jan-2013) changed logStat to have the option to only take the log, but not make stationary

#Edited by Ryan Batt (22-Nov-2013) for syntax readability (semicolons, indentation, return() for functions, spacing between assignments etc.)
ARMApqREMLfunct <- function(pars, nX, pea, queue){ 
	badRes <- 10^10
  tryCatch(
	{
		##generates LL for ARMApq model.
		##translated from matlab function ARMApqREMLfunct.m by Nic Ziebarth and Tony Ives
		##10 March, 2009

		##pars should have values for both AR and MA components.
		##nX is the data, p = # of AR, q = # of MA

		##reoccuring issue is what should diagExtend output when size = 0. This comes up
		##when p=1 and q=1. I think any 1 x 1 matrix will work. 

		## this is standardized so that a(1)=1 and is not estimated
		p <- pea
		q <- queue
		nP <- length(pars)
		k <- max(p,q)
		nObs <- length(nX)
		d1n <- diag(1,nObs)
		b <- pars[1:p]

		##this slight change relative to Matlab case. If q=0 in old version, then a =[], which ##was no good
		if(q > 0){
			a <- c(1, pars[(p+1):nP])
		}else{
			a <- 1
		}

		B <- diagExtend(p,-1)
		B[1,1:p] <- b
		# eigs <- eigen(B)$values
		eigB <- max(abs(eigen(B)$values))
 
		if(eigB>=1 | max(abs(a))>10){
			# print("bad Eig")
			# flush.console()
			return(badRes)
		}
		## solve for stationary distribution
		##deal with case of no MA terms separsately
# ==============
# = TEST 1 START =
# ==============
		if (q > 0){ #LINE 45
			# A <- array(0,dim=c(k,q))
			A <- matrix(0, nrow=k, ncol=q)
			# A[1,1:q] <- a[2:(q+1)]
			A[1,] <- a[-1] # CHANGED
			if(k > 1){
				B <- diagExtend(k,-1)
				B[1,1:p] <- b
			}else{
				B <- matrix(b)
			}
			
			# CC <- matrix(c(B, A), ncol=q+k) # CHANGED
			# otherPart <- matrix(c(rep(0,k*q), diagExtend(q,1)), nrow=k, byrow=TRUE) # CHANGED *then* EDITED (before, nrow=q)
			# CC <- rbind(CC,otherPart)
			
			CC <- cbind(B,A)
			otherPart <- rbind(array(0,dim=c(k,q)), diagExtend(q,1))
			CC <- rbind(CC,t(otherPart))
		
			# Ve <- array(0,dim=c(k+q,k+q))
			Ve <- matrix(0, nrow=k+q, ncol=k+q)
			
			Ve[1,1] <- 1
			Ve[1,k+1] <- 1
			Ve[k+1,1] <- 1
			Ve[k+1,k+1] <- 1
						
			# sizeVe <- dim(Ve)
			# vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1))
			vecVe <- matrix(Ve, ncol=1)
			# vecVe <- c(Ve) # This should work ... always ... # CHANGED
			C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe
			# C <- array(C,dim=c(k+q,k+q))
			C <- matrix(C, nrow=k+q)
			Vstationary <- C[1:k,1:k]
		}else{
# ============
# = TEST 1 END =
# ============

# ================
# = TEST 2 START =
# ================
		   ##this is parsticular example of issue when p=1
# ================
# = TEST 3 START =
# ================
			if(p > 1){
				B <- diagExtend(p,-1)
				# sizeB <- ncol(B)
				# B[1,1:sizeB] <- b
				# B[1,1:p] <- b
				B[1,] <- b
			}else{
				B <- matrix(b)
			}
					
			Ve <- matrix(0, nrow=p, ncol=p)
			Ve[1,1] <- 1
			vecVe <- matrix(Ve, ncol=1)
# ==============
# = TEST 3 END =
# ==============
	#Test 4 = do not wrap kronecker() in matrix()
			# C <- solve(diag(1,(p+q)^2) - matrix(kronecker(B,B)))%*%vecVe
			C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe
	#End test 4
			# Vstationary <- array(C,dim=c(p+q,p+q))
			Vstationary <- matrix(C, ncol=p+q, nrow=p+q)
		}
# ==============
# = TEST 2 END =
# ==============
 
		## solve for V_T
    
		##MA component
    
		##This was added to separsately deal with case of no MA lags
		if(q > 0){ # if MA terms ...
			AA <- diag(1,q+1)
			for (i in 1:q){ # for each MA term
				AA <- AA + a[i+1]*diagExtend(q+1,-i)
			}
			AA <- t(AA)%*%AA
			# aa <- AA[1,1:ncol(AA)]
			aa <- AA[1,]
			# AA <- array(0,dim=c(nObs,nObs))
			AA <- matrix(0, ncol=nObs, nrow=nObs)
			for (i in 1:q){ # iterate through MA again
				AA <- AA + aa[i+1]*diagExtend(nObs,-i)
			}
			AA <- AA+t(AA)+aa[1]*d1n
		}else{ # if there are NO MA terms ...
		   AA <- d1n
		}

		BB <- diag(1,k)
		if(k > 1){ # if model order is greater than 1 ...
			for (i in 1:p){
				BB <- BB - b[i]*diagExtend(k,-i)
			}
		}
		A1 <- BB%*%Vstationary%*%t(BB)
		AA[1:k,1:k] <- A1

		## AR component
		W <- d1n
		for (i in 1:p){
		   W <- W - b[i]*diagExtend(nObs,-i)
		}
		invW <- solve(W)
		V <- invW%*%AA%*%t(invW)
    
		## compute LL function for extant data
		pick <- !is.na(nX)
		nObs <- sum(pick, na.rm=TRUE)
		X <- nX[pick]
		V <- V[pick,pick]
    
		invV <- solve(V)
    
		# U <- array(1,dim=c(nObs,1))
		U <- matrix(1, nrow=nObs, ncol=1)
		tU <- t(U)
		uVu <- tU%*%invV%*%U
		mu <- solve(uVu)%*%tU%*%invV%*%X
		H <- X - mu
 
		## condensed ML LL
		##s2 <- (t(H)%*%invV%*%H)/T;
		##LL <- .5*(T*log(2*pi)+T*log(s2)+log(det(V))+T);

		## condensed REML LL
		s2 <- (t(H)%*%invV%*%H)/(nObs-1)
		LL <- as.vector(.5*((nObs-1)*log(2*pi)+(nObs-1)*log(s2)+log(det(V))+log(det(uVu))+(nObs-1)))


		if(abs(Im(LL)) > 10^-6){
			# print('bad imag')
			# flush.console()
			return(badRes)
		}
		if(is.na(LL)){
			# print('missing ll')
			# flush.console()
			return(badRes)
		}
		
		return(LL)
		
	}, #end expr for tryCatch
	error=function(cond){return(badRes)} # handle for error
	) # end tryCatch function call
} ## end ARMApqREMLfunct



# diagExtend <- function(p,diagInd){
# 
# 	##This is a function with the functionality of diag in Matlab
# 	##In particular, this generates a matrix of size p x p with ones on
# 	##the diagInd diagonal. Restriction is abs(diagInd) <= p. Returns 0 matrix otw.
# 
# 	diagNew <- abs(diagInd)
# 	if(diagNew > p){
# 		return(0)
# 	}else{
# 		firstPart <- diag(1,p,p-diagNew)
# 	   secondPart <- array(0,dim=c(p,diagNew))
# 	   A <- cbind(secondPart,firstPart)
# 		if(diagInd >= 0){
# 			return(A)
# 		}else{
# 			return(t(A))
# 		}
# 	}
# }
## end diagExtend

#http://stackoverflow.com/questions/7745363/r-equivalent-to-diagx-k-in-matlab
#Changed by RDB
diagExtend <- function(size, offset=0){ 
	M <- matrix(0, size, size)
	M[row(M)+offset == col(M)] <- 1
	return(M)
}


Eig <- function(b){
	nb <- length(b)
	B <- diagExtend(nb,-1)
	B[1,1:nb] <- b
	eigs <- eigen(B)$values
	# eigB <- max(abs(eigs))
	return(eigs)
}

Period <- function(Eigs){
	(2*pi)/Im(Eigs)[1]
}

AICc <- function(nll, Kay, N){
	L <- exp(-nll)
	AIC <- 2*Kay - 2*log(L)
	AICc <- AIC + (2*Kay*(Kay+1))/(N-Kay-1)
	return(AICc)
}


# ====================
# = Down to business =
# ====================
OptAll <- function(pq, Xs, pMax=3, qMax=3, Method){
	Pea <- pq[1]
	Queue <- pq[2]
	Kay <- sum(pq)
	lo <- rep(-2, Kay)
	up <- -lo
	nPop <- 20*Kay
	if(Method=="Evolve"){
		Eck <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=nPop, itermax=50, strategy=2, p=0.3, trace=FALSE), nX=Xs, pea=Pea, queue=Queue)
		ba <- Eck$optim$bestmem
		ll <- Eck$optim$bestval
	}
	if(Method=="Anneal"){
		startpars <- c(rep(0, Kay))
		Eck <- GenSA(par=startpars, fn=ARMApqREMLfunct, upper=up, lower=lo, control=list(max.call=1E3, verbose=FALSE, smooth=FALSE, temperature=5230*1), nX=Xs, pea=Pea, queue=Queue)
		ba <- Eck$par
		ll <- Eck$value
	}

	
	bs0 <- ba[1:Pea]
	bs <- c(bs0, rep(NA,(pMax-Pea)))
	if(Queue==0){
		as <- rep(NA,qMax)
	}else{
		as0 <- ba[(Pea+1):Kay]
		as <- c(as0, rep(NA,(qMax-Queue)))
	}
	eigs <- Eig(bs0)
	Lambda <- max(abs(eigs))
	if(is.complex(eigs)){
		cPer <- Period(eigs)
	}else{
		cPer <- NA
	}
	nobs <- sum(!is.na(Xs))
	aicc <- AICc(nll=ll, Kay=Kay, N=nobs)
	
	Nvals <- pMax+qMax+6
	list(matrix(c(Pea, Queue, Lambda, cPer, bs, as, ll, aicc), ncol=Nvals))
	
}

ARMAfit <- function(X0, dName, pMax=3, qMax=3, Method){
	X <- X0[,dName]
	pqO <- as.matrix(expand.grid("P"=1:pMax, "Q"=0:qMax))
	colnames(pqO) <- NULL
	
	Yikes <- apply(pqO, MARGIN=1, OptAll, Xs=X, pMax=pMax, qMax=qMax, Method=Method)
	Yikes2 <- matrix(unlist(Yikes), ncol=(pMax+qMax+6), byrow=TRUE)
	bNames <- paste(rep("b",pMax), 1:pMax, sep="")
	aNames <- paste(rep("a",pMax), 1:qMax, sep="")
	names <- c(bNames, aNames)
	names <- c("P", "Q", "Lambda", "Period", names, "nll", "AICc")
	dimnames(Yikes2) <- list(NULL, names)
	
	return(Yikes2)
}

logStat <- function(x, coluY, coluX, doStat=TRUE){
	ntot <- nrow(x)
	# fill <- rep(NA, ntot)
	
	if(sum(!is.na(x[,"Data"]))<5){return(x)}

	x2 <- Inf2NA(log(x[,"Data"]))
	# miss <- is.na(x2)
	if(doStat){
		tx <- 0:(ntot-1)
		xd <- residuals(lm(x2~tx, na.action=na.exclude))	
	}else{
		xd <- x2
	}
	
	# x[!miss,"Data"] <- xd
	x[,"Data"] <- xd
	return(x)
}

Stat <- function(x){
	ntot <- nrow(x)
	# fill <- rep(NA, ntot)
	
	if(sum(!is.na(x[,"Data"]))<5){return(x)}

	x2 <- Inf2NA(x[,"Data"])
	# miss <- is.na(x2)
	
	tx <- 0:(ntot-1)
	xd <- residuals(lm(x2~tx, na.action=na.exclude)) + mean(x2, na.rm=TRUE)
	
	# x[!miss,"Data"] <- xd
	x[,"Data"] <- xd
	return(x)
}



# Chl_Ext2 <- Chl_Ext[Chl_Ext[,"lakeid"]=="ME",]
# sChl_Ext2 <- Chl_Ext2
# sChl_Ext2[,"chlor"] <- logStat(Chl_Ext2[,"chlor"])
# Rprof("~/file.out")
# Chl_pq <- ddply(.data=sChl_Ext2, .variables=("lakeid"), .fun=ARMAfit, dName="chlor", .parallel=FALSE, .progress="time", Method="Evolve")
# Rprof(NULL)
# summaryRprof("~/file.out")



# Chl_Ext_z <- Chl_Ext
# Chl_Ext_z[,"Variable"] <- "chlor"
# names(Chl_Ext_z) <- c("year4", "lakeid", "Data", "Variable")
# Chl_Ext_z <- Chl_Ext_z[,c(1,2,4,3)]
# sChl_Ext <- ddply(Chl_Ext_z, .variables=c("lakeid","Variable"), .fun=logStat)
# Chl_pq <- ddply(.data=sChl_Ext, .variables=("lakeid"), .fun=ARMAfit, dName="Data", .parallel=FALSE, .progress="time", Method="Evolve")

# test <- sChl_Ext[sChl_Ext[,"lakeid"]=="TB",]

# ddply(.data=Chl_Ext, .variables=("lakeid"), .fun=mean)
minaicc <- function(x){
	x[which.min(x[,"AICc"]),]
}
# ddply(Chl_pq, .variables="lakeid", minaicc)

# =============
# = Test Data =
# =============
testFit <- function(Order=c(2,0,0), ps=c(0.4, -0.5), qs=NULL, Bounds=c(-1,1)){
	
	#Simulate data
	X <- as.numeric(arima.sim(model=list(order=Order, ar=ps), sd=0.1, n=30, n.start=100))
	X <- matrix(X, ncol=1, dimnames=list(NULL, "Data"))
	
	#Search for order and parameter values
	mapFit <- ARMAfit(X, Method="Evolve")
	
	dev.new(width=4, height=5)
	par(mfrow=c(3,1), mar=c(3,2,0.5,1.5), cex=1, ps=7, mgp=c(3, 0.4, 0), tcl=-0.35)
	D2 <- matrix(c(1,2,1,0,0,1), ncol=2) # 2-dimensional orders
	for(i in 1:nrow(D2)){
		
		arO <- D2[i,1] # AR order
		maO <- D2[i,2] # MA order
		
		lo <- rep(-2, sum(D2[i,])) # lower bounds for search
		up <- -lo
		
		#Fit the model for this particular order and save populations
		mapPop0 <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=NA, itermax=50, strategy=2, p=0.3, trace=FALSE, storepopfrom=1), nX=X, pea=arO, queue=maO)
		mapPop <- do.call(rbind.data.frame, mapPop0$member$storepop)
		generationCol <- rgb(t(col2rgb(tim.colors(n=nrow(mapPop)))), alpha=255, maxColorValue=255)
		
		map1 <- seq(Bounds[1], Bounds[2], length.out=100)
		if(length(lo)<2){
			map2 <- matrix(map1, ncol=1)
		}else{
			map2 <- expand.grid(map1, map1)
		}
		mapLL <- apply(map2, 1, ARMApqREMLfunct, nX=X, pea=arO, queue=maO)
		mapLL2 <- matrix(mapLL, ncol=100)
		trueM <- max(mapLL[mapLL<10^10]) # real maximum of the Z dimension (excluding arbitrarily high 10^10)
		
		if(length(lo)<2){
			plot(map1, mapLL2, type="l", xlab="", ylab="", lwd=2, ylim=c(min(mapLL2),trueM))
			mtext("ar1", side=1, line=1.25)
			mtext("nll", side=2, line=1.25)
			abline(v=ps, lty="dashed")
			points(mapPop0$optim$bestmem, mapPop0$optim$bestval, pch="X", cex=1.2)
		}else{ #else, if the model order is not <2 ..
			if(arO==2){
				Xlab <- "ar1"
				Ylab <- "ar2"
			}else{
				Xlab <- "ar1"
				Ylab <- "ma1"
			}
			par(mar=c(3,3,0.5,0.5), cex=1, ps=7, mgp=c(3, 0.4, 0), tcl=-0.35)
			image.plot(mapLL2, x=map1, y=map1, zlim=c(min(mapLL2),trueM), xlab="", ylab="", legend.width=0.5, legend.mar=4, cex=1, smallplot=c(0.87, 0.91, 0.3, 0.93))
			points(mapPop, pch=21, bg=generationCol, col="white", cex=1)
			abline(v=ps[1], h=ps[2], lwd=1, lty="dashed") # plot coordinates of true values
			par(mar=c(3,3,0.5,0.5), cex=1, ps=7, mgp=c(3, 0.4, 0), tcl=-0.35)
			mtext(Xlab, side=1, line=1.1)
			mtext(Ylab, side=2, line=1.25)
			
			

		}

	}

}
# =================
# = End test data =
# =================

# ============================================
# = Function to get the sigma^2 and sigma[E] =
# ============================================
getSE <- function(fat, data){
	#"fat" is "fatARMA" (from FatFrame_v1.R and fatARMA_v1.RData)
	# "data" is "finalFrame" from "FatFrame_v1.R"
	
	pars0 <- fat[,c("b1", "b2", "b3", "a1", "a2", "a3")]
	pars <- pars0[!is.na(pars0)]
	
	nxID <- fat[,c("location","variable")]
	nX <- data[is.element(data[,"location"],nxID[,"location"]) & is.element(data[,"variable"],nxID[,"variable"]), "Data"]
	pea <- fat[,"P"]
	queue <- fat[,"Q"]

	badRes <- data.frame("sigEps"=as.vector(NA), "sigE"=as.vector(NA), "sigInf"=as.vector(NA))
	tryCatch(
		{
			##generates LL for ARMApq model.
			##translated from matlab function ARMApqREMLfunct.m by Nic Ziebarth and Tony Ives
			##10 March, 2009

			##pars should have values for both AR and MA components.
			##nX is the data, p = # of AR, q = # of MA

			##reoccuring issue is what should diagExtend output when size = 0. This comes up
			##when p=1 and q=1. I think any 1 x 1 matrix will work. 

			## this is standardized so that a(1)=1 and is not estimated
			p <- pea
			q <- queue
			nP <- length(pars)
			k <- max(p,q)
			nObs <- length(nX)
			d1n <- diag(1,nObs)
			b <- pars[1:p]

			##this slight change relative to Matlab case. If q=0 in old version, then a =[], which ##was no good
			if(q > 0){
				a <- c(1, pars[(p+1):nP])
			}else{
				a <- 1
			}

			B <- diagExtend(p,-1)
			B[1,1:p] <- b
			# eigs <- eigen(B)$values
			eigB <- max(abs(eigen(B)$values))

			# if(eigB>=1 | max(abs(a))>10){
			# 	# print("bad Eig")
			# 	# flush.console()
			# 	return(badRes)
			# }
			## solve for stationary distribution
			##deal with case of no MA terms separsately

			if (q > 0){ #LINE 45
				A <- matrix(0, nrow=k, ncol=q)
				A[1,] <- a[-1]
				if(k > 1){
					B <- diagExtend(k,-1)
					B[1,1:p] <- b
				}else{
					B <- matrix(b)
				}

				CC <- cbind(B,A)
				otherPart <- rbind(array(0,dim=c(k,q)), diagExtend(q,1))
				CC <- rbind(CC,t(otherPart))

				Ve <- matrix(0, nrow=k+q, ncol=k+q)

				Ve[1,1] <- 1
				Ve[1,k+1] <- 1
				Ve[k+1,1] <- 1
				Ve[k+1,k+1] <- 1

				vecVe <- matrix(Ve, ncol=1)
				C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe #Ives 2010 supplement, Eq. A.3
				C <- matrix(C, nrow=k+q)
				Vstationary <- C[1:k,1:k]
			}else{

			   ##this is parsticular example of issue when p=1
				if(p > 1){
					B <- diagExtend(p,-1)
					B[1,] <- b
				}else{
					B <- matrix(b)
				}

				Ve <- matrix(0, nrow=p, ncol=p)
				Ve[1,1] <- 1
				vecVe <- matrix(Ve, ncol=1)
				C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe #Ives 2010 supplement, Eq. A.3 (case with no MA)
				Vstationary <- matrix(C, ncol=p+q, nrow=p+q)
			}

			## solve for V_T

			##MA component

			##This was added to separsately deal with case of no MA lags
			if(q > 0){ # if MA terms ...
				AA <- diag(1,q+1)
				for (i in 1:q){ # for each MA term
					AA <- AA + a[i+1]*diagExtend(q+1,-i)
				}
				AA <- t(AA)%*%AA
				aa <- AA[1,]
				AA <- matrix(0, ncol=nObs, nrow=nObs)
				for (i in 1:q){ # iterate through MA again
					AA <- AA + aa[i+1]*diagExtend(nObs,-i)
				}
				AA <- AA+t(AA)+aa[1]*d1n
			}else{ # if there are NO MA terms ...
			   AA <- d1n
			}

			BB <- diag(1,k)
			if(k > 1){ # if model order is greater than 1 ...
				for (i in 1:p){
					BB <- BB - b[i]*diagExtend(k,-i)
				}
			}
			A1 <- BB%*%Vstationary%*%t(BB)
			AA[1:k,1:k] <- A1

			## AR component
			W <- d1n
			for (i in 1:p){
			   W <- W - b[i]*diagExtend(nObs,-i)
			}
			invW <- solve(W)
			V <- invW%*%AA%*%t(invW)

			## compute LL function for extant data
			pick <- !is.na(nX)
			nObs <- sum(pick, na.rm=TRUE)
			X <- nX[pick]
			V <- V[pick,pick]

			invV <- solve(V)

			# U <- array(1,dim=c(nObs,1))
			U <- matrix(1, nrow=nObs, ncol=1)
			tU <- t(U)
			uVu <- tU%*%invV%*%U
			mu <- solve(uVu)%*%tU%*%invV%*%X
			H <- X - mu

			## condensed ML LL
			##s2 <- (t(H)%*%invV%*%H)/T;
			##LL <- .5*(T*log(2*pi)+T*log(s2)+log(det(V))+T);

			## condensed REML LL
			s2 <- (t(H)%*%invV%*%H)/(nObs-1)
			# LL <- as.vector(.5*((nObs-1)*log(2*pi)+(nObs-1)*log(s2)+log(det(V))+log(det(uVu))+(nObs-1)))
			
			sigEps <- sqrt(s2)
			sigE <- sigEps*sqrt(sum(a^2)) #the first element of a, a0, is 1. Ives et al. 2010, supplement after eq. A.4: "Here, without loss of generality, we have assumed a0=1."
			# sigInf <- c(Vstationary)[1] #stationary (co)variance of Xt (Xt is 1,1, Xt-1 is 2,2, ...); from online supplement to Ziebarth et al. 2010, "The variance of the stationary distribution of the ARMA(3,2) process, sigma[inf]^2, is given by the top-left element of V", where V is the covariance matrix of the stationary distribution of the ARMA(3,2) process. According to Ives 2010 supplement, the top k-by-k submatrix of W[inf] (covariance matrix of stationary distribution Z[inf], which includes xt, xt-1, ... xt-p, epst, epst-1, ... epst-q), is V[inf], the stationary distribution for Xt. It is unclear (to me, due to ignorance) if the "V" in Ziebarth is equivalent to Ives V[inf] or Ives W[inf], but as long as Ziebarth V corresponds to one of those, but I think the V's from the 2 papers are the same (indexing in this code suggests so)
			# ======================================
			# = Making correction to sigmaInf for _v5 =
			# ======================================
			sigInf <- sqrt(c(Vstationary)[1])*sigE # Tony said that Vstationary[1,1] was sig2Inf/sig2E
			
			sigs <- data.frame("sigEps"=as.vector(sigEps), "sigE"=as.vector(sigE), "sigInf"=as.vector(sigInf))
			
			# ==============================
			# = Compute Cholesky Residuals =
			# ==============================
			# Matlab code from Tony's email:
	        # % Cholesky residuals
	        # iD=chol(V)';
	        # D=iD\eye(length(iD));
	        # cH=D*H;
	
			# Example input/output from Tony:
			# tV <- matrix(c(2.6162, 2.0563, 1.6162, 1.2703, 0.9985, 2.0563, 2.6162, 2.0563, 1.6162, 1.2703, 1.6162, 2.0563, 2.6162, 2.0563, 1.6162, 1.2703, 1.6162, 2.0563, 2.6162, 2.0563, 0.9985, 1.2703, 1.6162, 2.0563, 2.6162), ncol=5, byrow=TRUE)
			# tH <- c(-2.0133, -0.2254, 1.5622, 2.3660, 1.2209)
			# tcH <- c(-1.2447, 1.3571, 1.7393, 1.1381, -0.6388)
	
			# Step 1: iD=chol(V)';
			iD <- t(chol(V)) #iD is a 30x30 matrix, with the top-right triangle being 0
		
			# Step 2: iD\eye(length(iD)) 
			# would be like A\B = A^1*B = solve(A)%*%b. This is "backslash-division" or "left division"
			# Note that in matlab, length() takes the length of the longest dimension of a matrix, whereas R returns the total number of elements
			# so in this context, length() from matlab becomes max(dim())
			# D <- solve(iD)%*%diag(max(dim(iD)))
			D <- solve(iD) #%*%diag(max(dim(iD))) the part that's commented-out is simply the identity matrix for the inverse of iD.
			
			# Step 3: cH=D*H;
			# First thing to note is that matlab's "*" is "%*%" in R
			cH <- D%*%H

			return(list("sigs"=sigs, "Resids"= cH))

		}, #end expr for tryCatch
		error=function(cond){return(badRes)} # handle for error
		) # end tryCatch function call
}





	gev.fit2 <- function (xdat, ydat = NULL, mul = NULL, sigl = NULL, shl = NULL, mulink = identity, siglink = identity, shlink = identity, show = FALSE, method = "Nelder-Mead", maxit = 10000, ...){

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
	    z$cov <- tryCatch({solve(x$hessian)}, error=function(cond)NA)#I should use the hessian() function in package numDeriv to calculate this, and compare it to the hessian from the optimization algorithm.
		if(all(!is.na(z$cov))){
			z$se <- sqrt(diag(z$cov))
		}else{
			z$se <- rep(NA, length(init))
		}
	    # z$se <- sqrt(diag(z$cov))
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


# ============================================================
# = The following functions used to be in fatARMA_summary_v4 =
# ============================================================
fWeighted <- function(x){
	#http://www.ssc.wisc.edu/~bhansen/718/NonParametrics15.pdf
	# more accurate at: http://machinelearning102.pbworks.com/w/file/fetch/47699411/aic_reg.pdf

	finiteLogic <- all(!is.finite(x[,"AICc"])) #If all of the AICc's are infinite
	missLogic <- all(is.na(x[,"AICc"])) #If all of the AICc's are missing
	if(finiteLogic | missLogic){  #don't compute the AICc-weighted averages (one of these may be true depending on wether I first converted Inf to NA or not)
		return(data.frame(x, "wLambda"=NA, "wOrder"=NA))
	}else{
		minAIC <- min(x[,"AICc"], na.rm=TRUE)

		eaic <- exp(-0.5*(x[,"AICc"]-minAIC))

		saicL <- sum(eaic[!is.na(x[,"Lambda"])], na.rm=TRUE)
		saicO <- sum(eaic[!is.na(x[,"Order"])], na.rm=TRUE)
		saicEps <- sum(eaic[!is.na(x[,"sigEps"])], na.rm=TRUE)
		saicE <- sum(eaic[!is.na(x[,"sigE"])], na.rm=TRUE)
		saicInf <- sum(eaic[!is.na(x[,"sigInf"])], na.rm=TRUE)

		wsL <- eaic/saicL
		wdO <- eaic/saicO
		wEps <- eaic/saicEps
		wE <- eaic/saicE
		wInf <- eaic/saicInf

		wLambda <- sum(wsL*x[,"Lambda"], na.rm=TRUE)
		wOrder <- sum(wdO*x[,"Order"], na.rm=TRUE)
		wSigEps <- sum(wEps*x[,"sigEps"], na.rm=TRUE)
		wSigE <- sum(wE*x[,"sigE"], na.rm=TRUE)
		wSigInf <- sum(wInf*x[,"sigInf"], na.rm=TRUE)

		return(data.frame(x, "wLambda"=wLambda, "wOrder"=wOrder, "wSigEps"=wSigEps, "wSigE"=wSigE, "wSigInf"=wSigInf))
	}
}

lvl_return_res <- function(x, level=1){
	lvl0 <- as.numeric(x[paste("Level",level, "_residual", sep="")])
	a0 <- as.numeric(x[c("residual_mu_0","residual_sig_0","residual_sh_0")])
	nExts0 <- as.numeric(x["N"])
	TS_Duration0 <- as.numeric(x["Duration"])
	result <- lvlX_ReturnTime(lvl=lvl0, a=a0, nExts=nExts0, TS_Duration=TS_Duration0)
	# names(result) <- NULL
	# row.names(result) <- NULL
	result
}
# ====================================================
# = End functions that used to be in fatARMA_Summary =
# ====================================================
