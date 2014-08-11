

# ================
# = Load scripts =
# ================
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") #also loads GenSA and DEoptim packages
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R")

set.seed(5)

# =================
# = Tony Function =
# =================
# =========================
# Series of yearly maxima =
# =========================
yearly.Max <- function(X=NULL, n.per.year=NULL){
  n <- length(X)
  n.year <- n/n.per.year
  
  Y <- array(0,c(n.year,1))
  for(i in 1:n.year){
    Y[i] <- max(X[(1+(i-1)*n.per.year):(i*n.per.year)])
  }
  return(Y)
}

yearly.Max2 <- function(X=NULL, n.per.year=NULL){
  n <- length(X)
  n.year <- n/n.per.year
  
  Y <- array(0,c(n.year,2))
  Y0 <- c()
  x.index <- 1:length(X)
  for(i in 1:n.year){
	t.window <- (1+(i-1)*n.per.year):(i*n.per.year)
	Y0[i] <- x.index[t.window][which.max(X[t.window])]
  }
  Y[,1] <- Y0
  Y[,2] <- X[Y0]
  return(Y)
}

# ===========================
# = First half of Tony code =
# ===========================
B <- function(b,phi) max(1,exp(b*(1+phi)))
S <- function(s,phi) (1+exp(-s*(1+phi)))^(-1)
FF <- function(a,X) (1+a*X)^(-1)

n.per.sample <- 25 # number of time steps per sample
samples.per.year <- 20 # number of samples per year
years <- 30 # number of years

tmax <- n.per.sample * samples.per.year * years
Tmax <- samples.per.year * years

variables <- array(0,c(tmax,4))
samples <- array(0,c(Tmax,3))
samples2 <- array(0,c(Tmax,5), dimnames=list(NULL, c("phi1","phi2","rB","rS","X")))

b <- 0.2
s <- 2
a <- 0.1

sd1 <- 0.25
sd2 <- 0.25

X <- 25
cR <- 1
Time <- 1

X.t <- X
compare <- c()

for(time in 1:tmax){
	
	phi1 <- rnorm(n=1,m=0,sd=sd1)
	phi2 <- rnorm(n=1,m=0,sd=sd2)
	
	#rB <- B(b,phi1 + phi2)
	#rS <- S(s,phi1 - phi2)
	rB <- B(b,phi1)
	rS <- S(s,phi2)
	
	variables[time,] <- c(phi1, phi2, rB, rS)
	
	cR <- cR*rB*rS

	# =================
	# = Added by Ryan =
	# =================
	X.t[time+1] <- rB*rS*X.t[time]*FF(a, X)^(1/n.per.sample) # Ryan
	# ============
	# = End Ryan =
	# ============
	
	
	if(time %% n.per.sample == 0){
		X <- X * cR * FF(a,X)
		samples[Time,] <- c(phi1, phi2, X)
		samples2[Time,] <- c(phi1, phi2, rB, rS, X)
		
		# =================
		# = Added by Ryan =
		# =================
		compare[Time] <- X.t[time+1] # Ryan
		# ============
		# = End Ryan =
		# ============
		
		cR <- 1
		Time <- Time + 1
	}
}



# # ##############################################
# # ############ Iterate to check    #############
# # ##############################################
simXis <- function(nreps=100){
	n.per.sample <- 25 # number of time steps per sample
	samples.per.year <- 20 # number of samples per year
	years <- 30 # number of years

	tmax <- n.per.sample * samples.per.year * years
	Tmax <- samples.per.year * years

	b <- 0.2
	s <- 2
	a <- 0.1

	sd1 <- 0.25
	sd2 <- 0.25

	xis <- array(0,c(nreps,3))
	for(i in 1:nreps){
		variables <- array(0,c(tmax,4))
		samples <- array(0,c(Tmax,3))

		X <- 25
		cR <- 1
		t.T <- 1
		for(t in 1:tmax){

			phi1 <- rnorm(n=1,m=0,sd=sd1)
			phi2 <- rnorm(n=1,m=0,sd=sd2)

			rB <- B(b,phi1)
			rS <- S(s,phi2)

			variables[t,] <- c(phi1, phi2, rB, rS)

			cR <- cR*rB*rS

			if(t %% n.per.sample == 0){
				X <- X * cR * FF(a,X)
				samples[t.T,] <- c(phi1, phi2, X)

				cR <- 1
				t.T <- t.T + 1
			}
		}

		Y <- yearly.Max(samples[,1],n.per.year=samples.per.year)
		gev.envir1 <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

		Y <- yearly.Max(samples[,2],n.per.year=samples.per.year)
		gev.envir2 <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

		Y <- yearly.Max(samples[,3],n.per.year=samples.per.year)
		gev.pop <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

		xis[i,] <- c(gev.envir1, gev.envir2, gev.pop)
	}
	return(xis)
}

xis <- simXis(500)


save.image(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatBirth_Sim.RData")


