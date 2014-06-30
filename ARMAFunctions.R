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
		eigB <- max(abs(eigen(B)$values))
 
		if(eigB>=1 | max(abs(a))>10){
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
		
			# Ve <- array(0,dim=c(k+q,k+q))
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

			if(p > 1){
				B <- diagExtend(p,-1)
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
			return(badRes)
		}
		if(is.na(LL)){
			return(badRes)
		}
		
		return(LL)
		
	}, #end expr for tryCatch
	error=function(cond){return(badRes)} # handle for error
	) # end tryCatch function call
} ## end ARMApqREMLfunct


# ===========================================================================
# = Place 1's on diagonal, or on superdiagonal 'offset' rows above diagonal =
# ===========================================================================
diagExtend <- function(size, offset=0){ 
	M <- matrix(0, size, size)
	M[row(M)+offset == col(M)] <- 1
	return(M)
}

# ========================
# = Calculate eigenvalue =
# ========================
Eig <- function(b){
	nb <- length(b)
	B <- diagExtend(nb,-1)
	B[1,1:nb] <- b
	eigs <- eigen(B)$values
	# eigB <- max(abs(eigs))
	return(eigs)
}

# ====================================
# = Calculate period from eigenvalue =
# ====================================
Period <- function(Eigs){
	(2*pi)/Im(Eigs)[1]
}

# ===========================
# = Calculate corrected AIC =
# ===========================
AICc <- function(nll, Kay, N){
	L <- exp(-nll)
	AIC <- 2*Kay - 2*log(L)
	AICc <- AIC + (2*Kay*(Kay+1))/(N-Kay-1)
	return(AICc)
}


# =============================================================================
# = Intermediate ARMA fitting function (above ARMApqREMLfunct, below ARMAfit) =
# =============================================================================
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

# =========================================
# = WRAPPER to fit all ARMA (call OptAll) =
# =========================================
ARMAfit <- function(X0, dName, pMax=3, qMax=3, Method){
	if(!library("GenSA", logical.return=TRUE)){stop("install package 'GenSA'")}
	if(!library("DEoptim", logical.return=TRUE)){stop("install package 'DEoptim'")}
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

# ============================================================================================
# = Take the log() of the series, return residuals of time trend (neg to NA, don't add mean) =
# ============================================================================================
logStat <- function(x, doStat=TRUE){
	ntot <- nrow(x)	
	if(sum(!is.na(x[,"Data"]))<5){return(x)}
	x2 <- Inf2NA(log(x[,"Data"]))
	if(doStat){
		tx <- 0:(ntot-1)
		xd <- residuals(lm(x2~tx, na.action=na.exclude))	
	}else{
		xd <- x2
	}	
	x[,"Data"] <- xd
	return(x)
}

# ===========================================
# = Remove linear time trend, add back mean =
# ===========================================
Stat <- function(x){
	ntot <- nrow(x)
	if(sum(!is.na(x[,"Data"]))<5){return(x)}
	x2 <- Inf2NA(x[,"Data"])	
	tx <- 0:(ntot-1)
	xd <- residuals(lm(x2~tx, na.action=na.exclude)) + mean(x2, na.rm=TRUE)
	x[,"Data"] <- xd
	return(x)
}

# =========================
# = Which is smallest AIC =
# =========================
minaicc <- function(x){
	x[which.min(x[,"AICc"]),]
}


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

	
	
	

# ==================================
# = Fit the GEV, w/ tryCatch added =
# ==================================
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
	
	# ==========================
	# = Different from GEV fit =
	# ==========================
	# The tryCatch and logicals below are the parts that are different between gev.fit and gev.fit2
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

# ==================================================================================================
# = Calculate return level for a data frame containing ARMA residuals (labeled as Level2_residual) =
# ==================================================================================================
lvl_return_res <- function(x, level=1){ # Called by fatARMA_Summary
	lvl0 <- as.numeric(x[paste("Level",level, "_residual", sep="")])
	a0 <- as.numeric(x[c("residual_mu_0","residual_sig_0","residual_sh_0")])
	nExts0 <- as.numeric(x["N"])
	TS_Duration0 <- as.numeric(x["Duration"])
	result <- lvlX_ReturnTime(lvl=lvl0, a=a0, nExts=nExts0, TS_Duration=TS_Duration0)
	result
}



# ===============================================================
# = Which element corresponds to a given quantile (probability) =
# ===============================================================
which.quantile <- function (x, probs, na.rm = FALSE){
  if (! na.rm & any (is.na (x)))
  return (rep (NA_integer_, length (probs)))

  o <- order (x)
  n <- sum (! is.na (x))
  o <- o [seq_len (n)]

  nppm <- n * probs - 0.5
  j <- floor(nppm)
  h <- ifelse((nppm == j) & ((j%%2L) == 0L), 0, 1)
  j <- j + h

  j [j == 0] <- 1
  o[j]
}

# =======================================================
# = Modify arima.sim to allow nonstationary simulations =
# =======================================================
arima.sim2 <- function (model, n, rand.gen = rnorm, innov = rand.gen(n, ...), 
    n.start = NA, start.innov = rand.gen(n.start, ...), ...) 
{
    if (!is.list(model)) 
        stop("'model' must be list")
    if (n <= 0L) 
        stop("'n' must be strictly positive")
    p <- length(model$ar)
    if (p) {
        minroots <- min(Mod(polyroot(c(1, -model$ar))))
		# CHANGED removed this stop to permit nonstationary simulation
        # if (minroots <= 1) 
        #     stop("'ar' part of model is not stationary")
    }
    q <- length(model$ma)
    if (is.na(n.start)) 
		#CHANGED added the max(), because the logic breaks if minroots is <= 1
        n.start <- p + q + ifelse(p > 0, ceiling(6/log(max(minroots, 1.1))), 
            0)
    if (n.start < p + q) 
        stop("burn-in 'n.start' must be as long as 'ar + ma'")
    d <- 0
    if (!is.null(ord <- model$order)) {
        if (length(ord) != 3L) 
            stop("'model$order' must be of length 3")
        if (p != ord[1L]) 
            stop("inconsistent specification of 'ar' order")
        if (q != ord[3L]) 
            stop("inconsistent specification of 'ma' order")
        d <- ord[2L]
        if (d != round(d) || d < 0) 
            stop("number of differences must be a positive integer")
    }
    if (!missing(start.innov) && length(start.innov) < n.start) 
        stop(sprintf(ngettext(n.start, "'start.innov' is too short: need %d point", 
            "'start.innov' is too short: need %d points"), n.start), 
            domain = NA)
    x <- ts(c(start.innov[seq_len(n.start)], innov[1L:n]), start = 1 - 
        n.start)
    if (length(model$ma)) {
        x <- filter(x, c(1, model$ma), sides = 1L)
        x[seq_along(model$ma)] <- 0
    }
    if (length(model$ar)) 
        x <- filter(x, model$ar, method = "recursive")
    if (n.start > 0) 
        x <- x[-(seq_len(n.start))]
    if (d > 0) 
        x <- diffinv(x, differences = d)
    as.ts(x)
}

# =============================
# = Random Numbers via Hurdle =
# =============================
r.jump.diff <- function(n, mu=0, sdev=1, lmu=0, lsdev=1, qFunc="rbinom", qArgs=list(n=n, size=1, prob=0.2)){
	rnorm(n=n, mean=mu, sd=sdev) + do.call(qFunc, qArgs)*rlnorm(n=n, meanlog=lmu, sdlog=lsdev)
}

# =================
# = Coin Toss RNG =
# =================
# set.seed(1); rTest <- rnorm(4)

rcoin <- function(N, nCoins=4, phi=0.5, acc=FALSE){
	# N is the number of random numbers to generate
	# nCoins is the number of coin flips ... it behaves somewhat like a variance parameter. The higher it is, the bigger the shocks (in either direction)
	# If acc is set to TRUE, the print-out will be a matrix where the first row of each column contains output that would be equivalent to acc=FALSE. The last row is the shock (startShock) that initializes the process, and the values above the bottom row show how the shocks accmulate after each coin flip. If N=100 and nCoins=4, you'll get a 4x100 matrix. Note that if you look from the last to first row in a column, successive rows that have the same value indicate that the coin flip was a 0.
	Funcall <- function(f, ...){
		f(...)	
	} 
	Iterate <- function(f, n = 1){
	    function(x){
			Reduce(Funcall, rep.int(list(f), n), x, right=TRUE, accumulate=acc)
		} 
	}
	startShock <- rnorm(N)
	if(acc){
		matrix(c(unlist(
			Iterate(f=function(x){coin <- rbinom(n=N, size=1, prob=0.5); shock <- rnorm(N, sd=1); Z <- rnorm(N, sd=0.01); x * (phi + Z + coin*shock)}, n=nCoins)(startShock)
			)), ncol=N, byrow=T)
	}else{
		Iterate(f=function(x){coin <- rbinom(n=N, size=1, prob=0.5); shock <- rnorm(N, sd=1); Z <- rnorm(N, sd=0.01); x * (phi + Z + coin*shock)}, n=nCoins)(startShock)	
	}
}

coin2 <- function(N, nCoins, nInter=nCoins){
	nShocks <- matrix(NA, nrow=N, ncol=min(nCoins, nInter))
	
	for(nsim in 1:N){
		coins <- rbinom(n=nCoins, size=1, prob=0.5)
		shocks <- rnorm(n=nCoins)
		for(i in 1:min(nCoins, nInter)){
			 nShocks[nsim, i] <- sum(combn(1:min(nCoins, nInter), i, function(x)sum(coins[x]*shocks[x])))
		}
	}
	rowSums(nShocks)
}

coin3 <- function(N, nCoins, nInter=nCoins){
	nShocks <- matrix(NA, nrow=N, ncol=min(nCoins, nInter))
	
	for(nsim in 1:N){
		coins <- rbinom(n=nCoins, size=1, prob=0.5)
		for(i in 1:min(nCoins, nInter)){
			 nShocks[nsim, i] <- sum(combn(1:min(nCoins, nInter), i, function(x)prod(coins[x])*(rnorm(1))))
		}
	}
	rowSums(nShocks)
}


# ======================================================================================
# = Simulate a variety of ARMA time series and fit GEV to maxima (1 maximum per nYear) =
# ======================================================================================
myFatSim <- function(x, nYear=35){
	N <- x[,"N"]
	# AR Coefficients
	if(x[,"P"]>=0){
		# arC <- runif(-1, 1, n=(x[,"P"]))
		arC <- x[,"P"] #seq(1/(x[,"P"]+1), x[,"P"]/(x[,"P"]+1), length.out=x[,"P"])
		Lambda <- max(abs(Eig(arC)))
		minRoot <- min(Mod(polyroot(c(1, -arC))))
	}else{
		arC <- NULL
		Lambda <- NA
		minRoot <- NA
	}

	# MA Coefficients
	if(x[,"Q"]>=1){
		maC <- runif(min=-1, max=1, n=(x[,"Q"]))
	}else{
		maC <- NULL
	}
	
	simOrder <- c(1,0,x[,"Q"]) #c(x[,"P"], 0, x[,"Q"])
	
	
	fullTS <- rep(NA, N*nYear)
	maxTS <- rep(NA, nYear)
	
	if(x[,"Distribution"]!="JDD"){
		simResid <- switch(as.character(x[,"Distribution"]),
			normal=rnorm(N*nYear, 0, 1),
			t=rt(N*nYear, 5),
			cauchy=rcauchy(n=N*nYear),
			lnorm=rlnorm(N*nYear, sdlog=0.65),
			r.jump.diff=r.jump.diff(n=N*nYear),
			rcoin=rcoin(N=N*nYear, nCoin=3),
			coin3=coin3(N=N*nYear, nCoin=4)
		)
		maxResid <- rep(NA, nYear)

		for(i in 1:nYear){
			simIndex <- (i*N-(N-1)):(i*N)
			tsimResid <- simResid[simIndex]
			maxResid[i] <- max(tsimResid)

			tsimTS <- c(arima.sim2(model=list(order=simOrder, ar=arC, ma=maC), n=N, innov=tsimResid))+0

			fullTS[simIndex] <- tsimTS
			maxTS[i] <- simIndex[which.max(tsimTS)]
		}

		gevRes0 <- gev.fit2(fullTS[maxTS])
		gevRes <- c(gevRes0$mle[c("sh_0", "mu_0", "sig_0")], "se"=gevRes0$se["sh_0"])
		Xi_resid <- gev.fit2(maxResid)$mle["sh_0"]
	
	}else{
	
		simResid <- NA
		maxResid <- rep(NA, nYear)

		# fullTS <- Stationary(log(jdd.sim(N=N*nYear)))[,2]
		# fullTS <- log(jdd.sim(N=N*nYear))
		# fullTS <- jdd.sim(N=N*nYear)
		for(i in 1:nYear){
			simIndex <- (i*N-(N-1)):(i*N)
			tsimResid <- NA
			maxResid[i] <- NA

			tsimTS <- jdd.sim(N=N, S0=2, lambda=0.01)
			# tsimTS <- fullTS[simIndex]
			fullTS[simIndex] <- tsimTS
			maxTS[i] <- simIndex[which.max(tsimTS)]
		}

		gevRes0 <- gev.fit2(fullTS[maxTS])
		gevRes <- c(gevRes0$mle[c("sh_0", "mu_0", "sig_0")], "se"=gevRes0$se["sh_0"])
		Xi_resid <- NA #gev.fit2(maxResid)$mle["sh_0"]
	}
	
	Order <- 1 + x[,"Q"] #x[,"P"] + x[,"Q"]
	
	arC2 <- c("AR1"=NA,"AR2"=NA,"AR3"=NA)
	arC2[0:x[,"P"]] <- arC
	maC2 <- c("MA1"=NA,"MA2"=NA,"MA3"=NA)
	maC2[0:x[,"Q"]] <- maC
	
	adf <- data.frame(x, "Order"=Order, "Lambda"=Lambda, "minRoot"=minRoot, t(arC2), t(maC2), "Xi"=gevRes["sh_0"], "mu"=gevRes["mu_0"], "sig"=gevRes["sig_0"], "residXi"=Xi_resid, "xi.se"=gevRes["se.sh_0"])
	list("summary"=adf, "fullTS"=fullTS, "maxTS"=maxTS)
}

