
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
		if (q > 0){
			# A <- array(0,dim=c(k,q))
			A <- matrix(0, nrow=k, ncol=q)
			# A[1,1:q] <- a[2:(q+1)]
			A[1,] <- a[-1]
			if(k > 1){
				B <- diagExtend(k,-1)
				B[1,1:p] <- b
			}else{
				B <- b
			}
			# =======
			# = ONE =
			# =======
			# CC <- cbind(B,A)
			# otherPart <- rbind(matrix(0, nrow=k, ncol=q), diagExtend(q,1))
			# CC <- rbind(CC,t(otherPart))
			
			CC <- matrix(c(B, A), ncol=q+k)
			otherPart <- matrix(c(rep(0,k*q),diagExtend(q,1)), nrow=q, byrow=TRUE)
			CC <- rbind(CC,otherPart)
			


			# Ve <- array(0,dim=c(k+q,k+q))
			Ve <- matrix(0, nrow=k+q, ncol=k+q)
			
			Ve[1,1] <- 1
			Ve[1,k+1] <- 1
			Ve[k+1,1] <- 1
			Ve[k+1,k+1] <- 1
						
			# sizeVe <- dim(Ve)
			# vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1))
			# vecVe <- matrix(Ve, ncol=1)
			vecVe <- c(Ve) # This should work ... always ...
			C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe
			# C <- array(C,dim=c(k+q,k+q))
			C <- matrix(C, nrow=k+q)
			Vstationary <- C[1:k,1:k]
		}else{
		   ##this is parsticular example of issue when p=1
		   if(p > 1){
				B <- diagExtend(p,-1)
				# sizeB <- ncol(B)
				# B[1,1:sizeB] <- b
				B[1,1:p] <- b
			}else{
				B <- b
			}

			# Ve <- matrix(0, nrow=p, ncol=p)
			Ve <- matrix(0, nrow=p, ncol=p)
			Ve[1,1] <- 1
			#  			sizeVe <- dim(Ve)
			# vecVe <- matrix(Ve, nrow=(sizeVe[1]*sizeVe[2]), ncol=1)
			# diag(1,(p+q)^2) - matrix(kronecker(B,B))
			# vecVe <- matrix(Ve, ncol=1)
			vecVe <- c(Ve)
			C <- solve(diag(1,(p+q)^2) - matrix(kronecker(B,B)))%*%vecVe
			# Vstationary <- array(C,dim=c(p+q,p+q))
			Vstationary <- matrix(C, ncol=p+q)
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
			aa <- AA[1,1:ncol(AA)]
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
	if(Method=="Evolve"){
		Eck <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=NA, itermax=50, strategy=2, p=0.3, trace=FALSE), nX=Xs, pea=Pea, queue=Queue)
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
	fill <- rep(NA, ntot)

	tx <- 0:(ntot-1)
	x2 <- log(x[,"Data"])
	miss <- is.na(x2)
	
	xd <- residuals(lm(x2~tx))
	x[!miss,"Data"] <- xd
	return(x0)
}

library("DEoptim")
library("GenSA")

# Chl_Ext2 <- Chl_Ext[Chl_Ext[,"lakeid"]=="ME",]
# sChl_Ext2 <- Chl_Ext2
# sChl_Ext2[,"chlor"] <- logStat(Chl_Ext2[,"chlor"])
# Rprof("~/file.out")
# Chl_pq <- ddply(.data=sChl_Ext2, .variables=("lakeid"), .fun=ARMAfit, dName="chlor", .parallel=FALSE, .progress="time", Method="Evolve")
# Rprof(NULL)
# summaryRprof("~/file.out")

Chl_Ext_z <- Chl_Ext
Chl_Ext_z[,"Variable"] <- "chlor"
names(Chl_Ext_z) <- c("year4", "lakeid", "Data", "Variable")
Chl_Ext_z <- Chl_Ext_z[,c(1,2,4,3)]
# sChl_Ext <- ddply(Chl_Ext_z, .variables=c("lakeid","Variable"), .fun=logStat)
Chl_pq <- ddply(.data=Chl_Ext_z, .variables=("lakeid"), .fun=ARMAfit, dName="Variable", .parallel=FALSE, .progress="time", Method="Evolve")

# ddply(.data=Chl_Ext, .variables=("lakeid"), .fun=mean)
minaicc <- function(x){
	x[which.min(x[,"AICc"]),]
}
ddply(Chl_pq, .variables="lakeid", minaicc)




# #plot out likelihood landscape
map1 <- seq(-1, 1, length.out=100)
map2 <- expand.grid(map1, map1)
mapLL <- apply(map2, 1, ARMApqREMLfunct, nX=Chl_Ext2[,3], pea=2, queue=0)
mapLL2 <- matrix(mapLL, ncol=100)
trueM <- max(mapLL[mapLL<10^10])
image.plot(mapLL2, x=map1, y=map1, zlim=c(0,trueM))
