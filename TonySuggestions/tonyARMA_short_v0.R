
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
		eigs <- eigen(B)$values
		eigB <- max(abs(eigs))
 
		if(eigB>=1 | max(abs(a))>10){
			return(badRes)
			print("bad Eig")
			flush.console()
		}
		## solve for stationary distribution
		##deal with case of no MA terms separsately
# ================
# = TEST 1 START =
# ================
		if (q > 0){
			A <- array(0,dim=c(k,q))
			A[1,1:q] <- a[2:(q+1)]
			if(k > 1){
				B <- diagExtend(k,-1)
				B[1,1:p] <- b
			}else{
				B <- matrix(b)
			}
			CC <- cbind(B,A)
			otherPart <- rbind(array(0,dim=c(k,q)), diagExtend(q,1))
			CC <- rbind(CC,t(otherPart))

			Ve <- array(0,dim=c(k+q,k+q))
			Ve[1,1] <- 1
			Ve[1,k+1] <- 1
			Ve[k+1,1] <- 1
			Ve[k+1,k+1] <- 1

			sizeVe <- dim(Ve)
			vecVe <- array(Ve,dim=c(sizeVe[1]*sizeVe[2],1))
			C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe
			C <- array(C,dim=c(k+q,k+q))
			Vstationary<-C[1:k,1:k]
		}else{
# ==============
# = Test 1 END =
# ==============

# ================
# = TEST 2 START =
# ================
		   ##this is parsticular example of issue when p=1
# ================
# = TEST 3 START =
# ================
		   if(p > 1){
				B <- diagExtend(p,-1)
				sizeB <- ncol(B)
				B[1,1:sizeB] <- b
			}else{
				B <- matrix(b)
			}

			Ve <- matrix(0, nrow=p, ncol=p)
			Ve[1,1] <- 1
 			sizeVe <- dim(Ve)
			vecVe <- matrix(Ve, nrow=(sizeVe[1]*sizeVe[2]), ncol=1)
# ==============
# = TEST 3 END =
# ==============
			C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe
			Vstationary <- array(C,dim=c(p+q,p+q))
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
			sizeAA <- ncol(AA)
			aa <- AA[1,1:sizeAA]
			AA <- array(0,dim=c(nObs,nObs))
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
    
		U <- array(1,dim=c(nObs,1))
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
			print('bad imag')
			flush.console()
			return(badRes)
		}
		if(is.na(LL)){
			print('missing ll')
			flush.console()
			return(badRes)
		}
		
		return(LL)
		
	}, #end expr for tryCatch
	error=function(cond){return(badRes)} # handle for error
	) # end tryCatch function call
} ## end ARMApqREMLfunct



diagExtend <- function(p,diagInd){

	##This is a function with the functionality of diag in Matlab
	##In particular, this generates a matrix of size p x p with ones on
	##the diagInd diagonal. Restriction is abs(diagInd) <= p. Returns 0 matrix otw.

	diagNew <- abs(diagInd)
	if(diagNew > p){
		return(0)
	}else{
		firstPart <- diag(1,p,p-diagNew)
	   secondPart <- array(0,dim=c(p,diagNew))
	   A <- cbind(secondPart,firstPart)
		if(diagInd >= 0){
			return(A)
		}else{
			return(t(A))
		}
	}
}
## end diagExtend

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

#Generate fake data for optimization testing
X <- as.numeric(arima.sim(model=list(order=c(1,0,0), ar=c(0.75)), sd=0.5, n=50, n.start=100))
lo <- rep(-2, 1)
up <- rep(2, 1)
# optim(par=startpars, fn=ARMApqREMLfunct, method="Nelder-Mead", nX=X, p=2, q=1)
# optim(par=startpars, fn=ARMApqREMLfunct, method="SANN", nX=X, p=2, q=1, control=list(maxit=1E5))

# library("GenSA")
# system.time(
# 	blah <- GenSA(par=startpars, fn=ARMApqREMLfunct, upper=up, lower=lo, control=list(max.call=100, verbose=FALSE), nX=X, pea=1, queue=1)
# 	blah$value
# 	blah$par
# )

library("DEoptim")
blah <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=NA, itermax=50, strategy=2, p=0.3, trace=FALSE), nX=X, pea=1, queue=1)
# blah$optim$bestmem

pqO <- as.matrix(expand.grid(1:3, 1:3))
colnames(pqO) <- NULL

# test <- blah$member$storepop
# test2 <- do.call(rbind.data.frame, test)
# genCols <- rgb(t(col2rgb(tim.colors(n=nrow(test2)))), alpha=255, maxColorValue=255)
# dev.new()
# plot(test2[,1:2], col=genCols)
# abline(h=0.4, v=-0.8)


# #plot out likelihood landscape
# map1 <- seq(-1, 1, length.out=100)
# map2 <- expand.grid(map1, map1)
# mapLL <- apply(map2, 1, ARMApqREMLfunct, nX=X, pea=1, queue=1)
# mapLL2 <- matrix(mapLL, ncol=100)
# trueM <- max(mapLL[mapLL<10^10])
# image.plot(mapLL2, x=map1, y=map1, zlim=c(0,trueM))


# ====================
# = Down to business =
# ====================
#Step 1: Generate range of ARMA orders to be fitted
pqO <- as.matrix(expand.grid("P"=1:3, "Q"=0:3))
colnames(pqO) <- NULL
# pq2 <- as.matrix(expand.grid("P"=1:3, "Q"=0:3))
# pqO <- cbind(pq2, pqO0)


#Step 2: Optimize parameters for each set of ARMA models
OptAll <- function(pq, Xs){
	Pea <- pq[1]
	Queue <- pq[2]
	Kay <- sum(pq)
	lo <- rep(-2, Kay)
	up <- -lo
	Eck <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=NA, itermax=50, strategy=2, p=0.3, trace=FALSE), nX=X, pea=Pea, queue=Queue)
	ba <- Eck$optim$bestmem
	ll <- Eck$optim$bestval
	
	bs0 <- ba[1:Pea]
	bs <- c(bs0, rep(NA,(3-Pea)))
	if(Queue==0){
		as <- rep(NA,3)
	}else{
		as0 <- ba[(Pea+1):Kay]
		as <- c(as0, rep(NA,(3-Queue)))
	}
	# names <- paste(rep(c("a","b"),each=3), 1:3, sep="")
	# names <- c("P", "Q", names, "nll")
	eigs <- Eig(bs0)
	Lambda <- max(abs(eigs))
	if(is.complex(eigs)){
		cPer <- Period(eigs)
	}else{
		cPer <- NA
	}
	nobs <- sum(!is.na(X))
	aicc <- AICc(nll=ll, Kay=Kay, N=nobs)
	
	list(matrix(c(Pea, Queue, Lambda, cPer, bs, as, ll, aicc), ncol=12))
	
}

# names <- paste(rep(c("a","b"),each=3), 1:3, sep="")
# names <- c(names, "nll")
# Yikes <- matrix(NA, ncol=7, dimnames=list(NULL, c(names)))
# for(i in 1:nrow(pqO)){
# 	Yikes <- rbind(Yikes,OptAll(pq=pqO[i,], Xs=X))	
# }
# OptAll(pq=pqO[i,], Xs=X)

Yikes <- apply(pqO, MARGIN=1, OptAll, Xs=X)
Yikes2 <- matrix(unlist(Yikes), ncol=12, byrow=TRUE)
names <- paste(rep(c("b","a"),each=3), 1:3, sep="")
names <- c("P", "Q", "Lambda", "Period", names, "nll", "AICc")
dimnames(Yikes2) <- list(NULL, names)

#Step 3: Select the best-fitting ARMA model using AICc
#Step 4: Calculate the eigenvalue and period (if any) for the best-fitting model
#Step 5: Repeat for 400 variables from LTER ...
