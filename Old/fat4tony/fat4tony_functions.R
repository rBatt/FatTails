# Functions dervied from tonyARMA_short_v4
# ARMApqREMLfunct from Ives et al. 2010 Ecology, "Analysis of ecological time series with ARMA(p,q) models"
# Values reported by getSE based on Appendix to Ives et al. 2010, as well as Ziebarth et al. (2010) Ecology Letters "Weak population regulation in ecological time series", as well as the appendix of Ziebarth et al. (2010). 

# ====================
# = Function Summary =
# ====================
# Note: ARMAfit <== OptAll <== ARMApqREMLfunct. (where '<==' means the function on left calls the function on the right)

# 1) ARMApqREMLfunct: loss function for an ARMA model, based on Ives et al. 2010, modified by RDB. Called for fitting ARMA models.

# 2) diagExtend: Custom R function that mimics Matlab's diagx function. Reprogrammed by RDB to increase speed (based on stackoverflow post)

# 3) Eig: extract eigenvalues when supplied with the matrix of AR coefficients.

# 4) AICc: Corrected AICc based on NLL, # of parameters, and N

# 5) Period: Compute the period of a time series from the eigenvalues (which are from the AR coefficients)

# 6) OptAll: Given the ARMA order and data, find & return the best value for those ARMA coefficients, also return the NLL of the fit, the Period, and the Eigenvalue 

# 7) ARMAfit: top-level fitting, calls OptAll. OptAll only fits for a specific ARMA order. ARMAfit constructs a matrix of ARMA orders (1,0 to 3,3), and for each order calls OptAll to fit the parameter values. Meant for fitting several ARMA models (differing in order) to a single variable contained in a data frame that also contains many other variables (imagine a data frame with a column for time, a column for data, and a column for "variableName", such that the column "data" actually contains the observed values for several types of variables, with each type of variable named in the column "variableName".).  Combined with the function ddply(), ARMAfit makes it easy to fit a series of arma models to many time series whose observations are contained within a single data frame.

# 8) getSE: Sam function as ARMApqREMLfunct, but instead of returning the NLL, getSE returns the sigma[epsilon], sigma[infinity], and sigma[E]

# 9) minaicc: finds and returns the row of a data frame with the smallest AICc (returns actual contents, not just row index)


# ===========================
# = Begin loading functions =
# ===========================

#Edited by Ryan Batt (22-Nov-2013) for syntax readability (semicolons, indentation, return() for functions, spacing between assignments etc.). Other edits made later to improve speed, etc.
ARMApqREMLfunct <- function(pars, nX, pea, queue){ 
	badRes <- 10^10
  tryCatch( # I wrapped this in a tryCatch() so that the fitting process wouldn't crash if one time series threw an error. When there is an error, a likelihood of 10^10 is returned.
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
		eigB <- max(abs(eigen(B)$values))
 
		if(eigB>=1 | max(abs(a))>10){
			# print("bad Eig")
			# flush.console()
			return(badRes)
		}
		## solve for stationary distribution
		##deal with case of no MA terms separsately

		if (q > 0){
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
			C <- solve(diag(1,(k+q)^2) - kronecker(CC,CC))%*%vecVe
			C <- matrix(C, nrow=k+q)
			Vstationary <- C[1:k,1:k]
		}else{

		   ##this is parsticular example of issue when p=1

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

			C <- solve(diag(1,(p+q)^2) - kronecker(B,B))%*%vecVe
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



OptAll <- function(pq, Xs, pMax=3, qMax=3, Method){
	Pea <- pq[1] # the AR order
	Queue <- pq[2] # the MA order
	Kay <- sum(pq) # total order (not the k used in other functions to mean max(pq), so careful here)
	lo <- rep(-2, Kay) #the lower end of where the optimization routine should guess for the parameters
	up <- -lo #the uppend end of possible parameter values to be guessed by the optimization routine
	nPop <- 20*Kay # the "population size" for differential evolution; the number of starting points for this form of optimization
	if(Method=="Evolve"){ # implements differential evolution
		# "Eck" is a sound I make when programming things, especially when I know that the result of a command is a mess
		Eck <- DEoptim(fn=ARMApqREMLfunct, lower=lo, upper=up, control=list(NP=nPop, itermax=50, strategy=2, p=0.3, trace=FALSE), nX=Xs, pea=Pea, queue=Queue)
		ba <- Eck$optim$bestmem # get the b's and the a's from the optimization; i.e., the best estimates of the AR and MA coefficients
		ll <- Eck$optim$bestval # get the final (best) negative log likelihood for the fit
	}
	if(Method=="Anneal"){ #implements simulated annealing
		startpars <- c(rep(0, Kay))
		Eck <- GenSA(par=startpars, fn=ARMApqREMLfunct, upper=up, lower=lo, control=list(max.call=1E3, verbose=FALSE, smooth=FALSE, temperature=5230*1), nX=Xs, pea=Pea, queue=Queue)
		ba <- Eck$par # ARMA parameters
		ll <- Eck$value # the NLL for the best fit
	}

	
	bs0 <- ba[1:Pea]
	bs <- c(bs0, rep(NA,(pMax-Pea))) # AR coefficients
	if(Queue==0){
		as <- rep(NA,qMax)
	}else{
		as0 <- ba[(Pea+1):Kay]
		as <- c(as0, rep(NA,(qMax-Queue))) # MA coefficients
	}
	eigs <- Eig(bs0)
	Lambda <- max(abs(eigs)) # ||λ||
	if(is.complex(eigs)){
		cPer <- Period(eigs) #the period of cyclic time series
	}else{
		cPer <- NA
	}
	nobs <- sum(!is.na(Xs)) # the number of non-na observations
	aicc <- AICc(nll=ll, Kay=Kay, N=nobs) # corrected AIC
	
	Nvals <- pMax+qMax+6 #the number of columns in the returned data frame
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





minaicc <- function(x){
	x[which.min(x[,"AICc"]),]
}


# ============================================
# = Function to get the sigma^2 and sigma[E] =
# ============================================
getSE <- function(fat, data){
	#"fat" is the result of the arma fitting (e.g., from OptAll)
	# "data" is the data frame containing the time series used in the ARMA fitting
	
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
			sigInf <- c(Vstationary)[1] #stationary (co)variance of Xt (Xt is 1,1, Xt-1 is 2,2, ...); from online supplement to Ziebarth et al. 2010, "The variance of the stationary distribution of the ARMA(3,2) process, sigma[inf]^2, is given by the top-left element of V", where V is the covariance matrix of the stationary distribution of the ARMA(3,2) process. According to Ives 2010 supplement, the top k-by-k submatrix of W[inf] (covariance matrix of stationary distribution Z[inf], which includes xt, xt-1, ... xt-p, epst, epst-1, ... epst-q), is V[inf], the stationary distribution for Xt. It is unclear (to me, due to ignorance) if the "V" in Ziebarth is equivalent to Ives V[inf] or Ives W[inf], but as long as Ziebarth V corresponds to one of those, but I think the V's from the 2 papers are the same (indexing in this code suggests so)
			
			sigs <- data.frame("sigEps"=as.vector(sigEps), "sigE"=as.vector(sigE), "sigInf"=as.vector(sigInf))

			return(sigs)

		}, #end expr for tryCatch
		error=function(cond){return(badRes)} # handle for error
		) # end tryCatch function call
}


