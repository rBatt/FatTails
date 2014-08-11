

# ================
# = Load scripts =
# ================
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") #also loads GenSA and DEoptim packages
# source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R")

set.seed(35)

# =================
# = Tony Function =
# =================
# ======================================================================================
# Series of yearly maxima
# ======================================================================================
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
	# lastT <- which((1:time) %% n.per.sample == 0) # what is the last value of X.t that lined up with X?
	# last.XT <- ifelse(length(lastT)==0L, X, X.t[time]) # select the last value of X.t that lined up with X, or, if time < n.per.sample, then just use X
	# X.t[time+1] <- rB*rS*X.t[time]*FF(a, last.XT)^(1/n.per.sample) # Ryan
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

# =================
# = Added by Ryan =
# =================
# plot(head(samples[,3], 10), type="l"); lines(head(compare, 10), type="l", col="red")
# ============
# = End Ryan =
# ============


# ============================
# = Second half of Tony Code =
# ============================
# dev.new()
# par(mfrow=c(2,2))
# 
# # values of output for year 2
# annual.variables <- variables[(n.per.sample * samples.per.year+1):(2*(n.per.sample * samples.per.year)),]
# annual.samples <- samples[(samples.per.year+1):(2*samples.per.year),]
# 
# # A: within-year environmental variables (phi1 and phi2)
# plot((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,1], type='l', col='black', xlab='Time t within a year', ylab='Standardized value')
# lines((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,2], col='blue')
# points(1:samples.per.year,annual.samples[,1], col='black')
# points(1:samples.per.year,annual.samples[,2], col='blue')
# 
# # B: within-year population variables (B, S, and X)
# plot((1:(n.per.sample * samples.per.year))/n.per.sample, annual.variables[,3], type='l', col='green', xlab='Time t within a year', ylab='Standardized value', ylim=c(0.5,1.5), xlim=c(0, samples.per.year))
# lines((1:(n.per.sample * samples.per.year))/n.per.sample,annual.variables[,4], col='red')
# # points(1:samples.per.year, 3*annual.samples[,3]/max(annual.samples[,3]), col='black')
# par(new=TRUE);plot(1:samples.per.year, annual.samples[,3], col='black', type="p", yaxt="n", xlim=c(0, samples.per.year)); axis(side=4)
# 
# # C: among-year environmental variables (phi1 and phi2)
# plot(samples[,1],type='l')
# Y <- yearly.Max2(samples[,1],n.per.year=samples.per.year)
# points(Y)
# 
# gev.fit2(xdat=Y[,2], show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
# 
# # D: among-year environmental variables (phi1 and phi2)
# plot(samples[,3],type='l')
# Y <- yearly.Max2(samples[,3],n.per.year=samples.per.year)
# points(Y)
# 
# gev.fit2(xdat=Y[,2], show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
# # 
# # ##############################################
# # ############ Iterate to check    #############
# # ##############################################
# 
nreps <- 100

n.per.sample <- 25 # number of time steps per sample
samples.per.year <- 20 # number of samples per year
years <- 30 # number of years

tmax <- n.per.sample * samples.per.year * years
Tmax <- samples.per.year * years
# 
# b <- 0.2
# s <- 2
# a <- 0.1
# 
# sd1 <- 0.25
# sd2 <- 0.25

xis <- array(0,c(nreps,3))
for(i in 1:nreps){
	variables <- array(0,c(tmax,4))
	samples <- array(0,c(Tmax,3))
	
	X <- 25
	cR <- 1
	T <- 1
	for(t in 1:tmax){
		
		phi1 <- rnorm(n=1,m=0,sd=sd1)
		phi2 <- rnorm(n=1,m=0,sd=sd2)
		
		#rB <- B(b,phi1 + phi2)
		#rS <- S(s,phi1 - phi2)
		rB <- B(b,phi1)
		rS <- S(s,phi2)
		
		variables[t,] <- c(phi1, phi2, rB, rS)
		
		cR <- cR*rB*rS
		
		if(t %% n.per.sample == 0){
			X <- X * cR * FF(a,X)
			samples[T,] <- c(phi1, phi2, X)
			
			cR <- 1
			T <- T + 1
		}
	}
	
	Y <- yearly.Max(samples[,1],n.per.year=samples.per.year)
	gev.envir1 <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
	
	Y <- yearly.Max(samples[,2],n.per.year=samples.per.year)
	gev.envir2 <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
	
	Y <- yearly.Max(samples[,3],n.per.year=samples.per.year)
	gev.pop <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
	
	xis[i,] <- c(gev.envir1, gev.envir2, gev.pop)
	
	# show(xis[i,])
}

dev.new()
par(mfrow=c(4,1))
hist(xis[,1],main='Envir 1', xlim=c(-0.8,0.8), breaks=seq(-0.8,0.8,0.1))
hist(xis[,2],main='Envir 1', xlim=c(-0.8,0.8), breaks=seq(-0.8,0.8,0.1))
hist(xis[,3],main='Pop', xlim=c(-0.8,0.8), breaks=seq(-0.8,0.8,0.1))
hist(xis[,3] - (xis[,1]+xis[,2]), main='Pop-Envir', xlim=c(-1,1), breaks=seq(-1,1,0.1))

colMeans(xis)
mean(xis[,2] - xis[,1])
cor(xis)



