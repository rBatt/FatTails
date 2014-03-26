library("DEoptim")
library("GenSA")


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

logStat <- function(x, coluY, coluX){
	ntot <- nrow(x)
	# fill <- rep(NA, ntot)
	
	if(sum(!is.na(x[,"Data"]))<5){return(x)}

	tx <- 0:(ntot-1)
	x2 <- Inf2NA(log(x[,"Data"]))
	miss <- is.na(x2)
	
	xd <- residuals(lm(x2~tx, na.action=na.exclude))
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


