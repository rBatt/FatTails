# install.packages("ismev")
library(ismev)

# ============
# = Set Seed =
# ============
set.seed(1337)


# ===============================================
# = 1-1 Ricker Map (Simulation for Figure 1 (4) =
# ===============================================

sd <- 0

r <- 3.1
p <- 3

n <- 1000
step <- 1

x <- .5
X <- array(0, dim=n)
E <- array(0, dim=n)
for(t in 1:(step*n)){
	# fat-tailed residuals
	if(sd != 0) {
		e <- rlnorm(1,sdlog=sd)^expon
		E[t] <- e
		x <- x*exp(r*(1-x))*e
	}else{
		x <- x*exp(r*(1-x))
	}
	X[t] <- x
}

X <- X[step*(1:n)]
E <- E[step*(1:n)]
gX <- gev.fit(as.numeric(X), show=F)

# ================
# = Figure 1 (7) =
# ================
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS4_DeterRickerChaos.pdf", height=3.42, width=3.42)
# dev.new()
par(mfrow=c(2,2), mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)

plot(X, type="l", xlab="Time")
hist(X, main=bquote(GEV~xi==.(.001*round(1000*gX$mle[3]))),30)

w <- arima(X, order=c(p,0,p))
w
we <- matrix(0,nrow=p, ncol=p)
we[1,] <- w$coef[1:p]
if(p > 1) for(j in 2:p) we[j,j-1] <- 1
abs(eigen(we)$values)

Y <- w$residuals

gY <- gev.fit(as.numeric(Y), show=F)
plot(Y, type="l", ylab="residuals")
hist(Y, main=bquote(GEV~xi==.(.001*round(1000*gY$mle[3]))),30)

dev.off()


# ===========================
# = Figure 2 (5) Simulation =
# ===========================
# ---- Simulation ----
sd <- .1
expon <- 1

r <- 3.1
p <- 3

n <- 1000
step <- 1

x <- .5
X <- array(0, dim=n)
E <- array(0, dim=n)
for(t in 1:(step*n)){
	# fat-tailed residuals
	if(sd != 0) {
		e <- rlnorm(1,sdlog=sd)^expon
		E[t] <- e
		x <- x*exp(r*(1-x))*e
	}else{
		x <- x*exp(r*(1-x))
	}
	X[t] <- x
}

X <- X[step*(1:n)]
E <- E[step*(1:n)]


gX <- gev.fit(as.numeric(X), show=F)

# ---- Plot Figure 2 (8) ----
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS5_StochRickerChaos.pdf", height=5, width=3.42)
# dev.new()
par(mfrow=c(3,2), mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)
if(sd == 0) par(mfrow=c(2,2))

plot(X, type="l", xlab="Time")
hist(X, main=bquote(GEV~xi==.(.001*round(1000*gX$mle[3]))),30)

w <- arima(X, order=c(p,0,p))
w
we <- matrix(0,nrow=p, ncol=p)
we[1,] <- w$coef[1:p]
if(p > 1) for(j in 2:p) we[j,j-1] <- 1
abs(eigen(we)$values)

Y <- w$residuals

gY <- gev.fit(as.numeric(Y), show=F)
plot(Y, type="l", ylab="residuals")
hist(Y, main=bquote(GEV~xi==.(.001*round(1000*gY$mle[3]))),30)

if(sd != 0){
	
	gE <- gev.fit(as.numeric(E), show=F)	
	plot(E, type="l", ylab="environment", xlab="Time")
	hist(E, main=bquote(GEV~xi==.(.001*round(1000*gE$mle[3]))),30)
}
dev.off()



# =======================================
# = 'Repeats' Simulation / Figure 3 (6) =
# =======================================
nreps <- 10

n <- 1000
step <- 1

sd <- 0.1

expon <- 1

output <- data.frame(r=0, expon=0, gX = array(0, dim=6),gY = 0,gE = 0)
i <- 0
for(r in c(1, 2, 2.5, 2.6, 2.9, 3.1)){
	i <- i+1
	g <- data.frame(gX = array(0, dim=nreps),gY = 0,gE = 0)
	
	for(rep in 1:nreps){
	
		x <- 1
		X <- array(0, dim=n)
		E <- array(0, dim=n)
		for(t in 1:(step*n)){
			e <- 1
			# fat-tailed residuals
			if(sd != 0) {
				e <- rlnorm(1,sdlog=sd)^expon
				E[t] <- e
				x <- x*exp(r*(1-x))*e
			}else{
				x <- x*exp(r*(1-x))
			}
			X[t] <- x
		}
		X <- X[step*(1:n)]
		E <- E[step*(1:n)]
		
		gX <- gev.fit(as.numeric(X), show=F)
		g$gX[rep] = gX$mle[3]
		
		w <- arima(X, order=c(3,0,3))
		Y <- w$residuals
		
		gY <- gev.fit(as.numeric(Y), show=F)
		g$gY[rep] = gY$mle[3]
			
		gE <- gev.fit(as.numeric(E), show=F)
		g$gE[rep] = gE$mle[3]
	}
	output[i,] <- c(r, expon, colMeans(g))
}
output
xx <- as.factor(output$r)

# ---- Figure 3 (9) ----
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS6_xi_vs_r_3lines.pdf", height=3.42, width=3.42)
par(mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)
# dev.new()
plot(gX ~ xx, data=output, lty="blank", xlab="r", ylab=bquote(GEV~xi))
lines(gX ~ xx, data=output)
lines(gY ~ xx, data=output, col='blue')
lines(gE ~ xx, data=output, col='red')
dev.off()


# ================================================
# = 1-1 Ricker Map (Simulation for Figure 4(7)) =
# ================================================
# ---- Simulation ----
sd <- .15
expon <- 6

r <- 1
p <- 3

n <- 1000
step <- 1

x <- .5
X <- array(0, dim=n)
E <- array(0, dim=n)
for(t in 1:(step*n)){
	# fat-tailed residuals
	if(sd != 0) {
		e <- rlnorm(1,sdlog=sd)^expon
		E[t] <- e
		x <- x*exp(r*(1-x))*e
	}else{
		x <- x*exp(r*(1-x))
	}
	X[t] <- x
}

X <- X[step*(1:n)]
E <- E[step*(1:n)]


gX <- gev.fit(as.numeric(X), show=F)

# ---- Plot Figure 4 (10) ----
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS7_StochRickerStable.pdf", height=5, width=3.42)
# dev.new()
par(mfrow=c(3,2), mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)
if(sd == 0) par(mfrow=c(2,2))

plot(X, type="l", xlab="Time")
hist(X, main=bquote(GEV~xi==.(.001*round(1000*gX$mle[3]))),30)

w <- arima(X, order=c(p,0,p))
w
we <- matrix(0,nrow=p, ncol=p)
we[1,] <- w$coef[1:p]
if(p > 1) for(j in 2:p) we[j,j-1] <- 1
abs(eigen(we)$values)

Y <- w$residuals

gY <- gev.fit(as.numeric(Y), show=F)
plot(Y, type="l", ylab="residuals")
hist(Y, main=bquote(GEV~xi==.(.001*round(1000*gY$mle[3]))),30)

if(sd != 0){
	
	gE <- gev.fit(as.numeric(E), show=F)	
	plot(E, type="l", ylab="environment", xlab="Time")
	hist(E, main=bquote(GEV~xi==.(.001*round(1000*gE$mle[3]))),30)
}
dev.off()



# =====================================================
# = Repeats (?) on Exponent; Simulation Figure 5 (8) =
# =====================================================
nreps <- 10

n <- 1000
step <- 1

sd <- 0.15
r <- 1

output <- data.frame(r=0, expon=0, gX = array(0, dim=6),gY = 0,gE = 0)
i <- 0
for(expon in 1:6){
	i <- i+1
	g <- data.frame(gX = array(0, dim=nreps),gY = 0,gE = 0)
	
	for(rep in 1:nreps){
	
		x <- .5
		X <- array(0, dim=n)
		E <- array(0, dim=n)
		for(t in 1:(step*n)){
			e <- 1
			# fat-tailed residuals
			if(sd != 0) {
				e <- rlnorm(1,sdlog=sd)^expon
				E[t] <- e
				x <- x*exp(r*(1-x))*e
			}else{
				x <- x*exp(r*(1-x))
			}
			X[t] <- x
		}
		X <- X[step*(1:n)]
		E <- E[step*(1:n)]
		
		gX <- gev.fit(as.numeric(X), show=F)
		g$gX[rep] = gX$mle[3]
		
		w <- arima(X, order=c(3,0,3))
		Y <- w$residuals
		
		gY <- gev.fit(as.numeric(Y), show=F)
		g$gY[rep] = gY$mle[3]
			
		gE <- gev.fit(as.numeric(E), show=F)
		g$gE[rep] = gE$mle[3]
	}
	output[i,] <- c(r, expon, colMeans(g))
}
output

xx <- as.factor(output$expon)

# ---- Figure 5 (8) ----
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS8_xi_vs_E_3lines.pdf", height=3.42, width=3.42)
# dev.new()
par(mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)
plot(gX ~ xx, data=output, lty="blank", xlab="E exponential", ylab=bquote(GEV~xi))
lines(gX ~ xx, data=output)
lines(gY ~ xx, data=output, col='blue')
lines(gE ~ xx, data=output, col='red')
dev.off()



# ===========================================================
# = Ricker Map (slow approach to equilibrium) Figure 6 (9) =
# ===========================================================
# ---- Simulation ----
sd <- .15
expon <- 6

r <- 0.1
p <- 3

n <- 1000
step <- 1

x <- .5
X <- array(0, dim=n)
E <- array(0, dim=n)
for(t in 1:(step*n)){
	# fat-tailed residuals
	if(sd != 0) {
		e <- rlnorm(1,sdlog=sd)^expon
		E[t] <- e
		x <- x*exp(r*(1-x))*e
	}else{
		x <- x*exp(r*(1-x))
	}
	X[t] <- x
}

X <- X[step*(1:n)]
E <- E[step*(1:n)]


gX <- gev.fit(as.numeric(X), show=F)

# ---- Plot Figure 6 (9) ----
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS9_StochRickerLowStable.pdf", height=5, width=3.42)
# dev.new()
par(mfrow=c(3,2), mar=c(2.5,2.5,1,0), cex=1, ps=8, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.35)
if(sd == 0) par(mfrow=c(2,2))

plot(X, type="l", xlab="Time")
hist(X, main=bquote(GEV~xi==.(.001*round(1000*gX$mle[3]))),30)

w <- arima(X, order=c(p,0,p))
w
we <- matrix(0,nrow=p, ncol=p)
we[1,] <- w$coef[1:p]
if(p > 1) for(j in 2:p) we[j,j-1] <- 1
abs(eigen(we)$values)

Y <- w$residuals

gY <- gev.fit(as.numeric(Y), show=F)
plot(Y, type="l", ylab="residuals")
hist(Y, main=bquote(GEV~xi==.(.001*round(1000*gY$mle[3]))),30)

if(sd != 0){
	
	gE <- gev.fit(as.numeric(E), show=F)	
	plot(E, type="l", ylab="environment", xlab="Time")
	hist(E, main=bquote(GEV~xi==.(.001*round(1000*gE$mle[3]))),30)
}
dev.off()

#
# # ======================================
# # = 1-1 Hassell (not a figure in supp) =
# # ======================================
# sd <- 0.15
# expon <- 1
#
# r <- .5
# wnorm <- 0
#
# n <- 1000
#
# x <- .2
# X <- array(0, dim=n)
# E <- array(0, dim=n)
# for(t in 1:n){
# 	# fat-tailed residuals
# 	if(sd != 0) {
# 		e <- rlnorm(1,sdlog=sd)^expon
# 		E[t] <- e
# 		x <- x*exp(r)/(1+x)*e
# 	}else{
# 		x <- x*exp(r)/(1+x)
# 	}
# 	X[t] <- x
# }
#
# X <- X[10*(1:100)]
# E <- E[10*(1:100)]
#
# gX <- gev.fit(as.numeric(X), show=F)
#
#
# # ---- Figure (not in supp) ----
# dev.new()
# par(mfrow=c(3,2), oma=c(0.1,0.1,1,0.1))
# if(sd == 0) par(mfrow=c(2,2))
#
# plot(X, type="l", xlab="Time")
# hist(X, main=paste('GEV.sh = ',.001*round(1000*gX$mle[3])),30)
#
# w <- arima(X, order=c(1,0,1))
# Y <- w$residuals
#
# gY <- gev.fit(as.numeric(Y), show=F)
# plot(Y, type="l", ylab="residuals")
# hist(Y, main=paste('GEV.sh = ',.001*round(1000*gY$mle[3])),30)
#
# if(sd != 0){
#
# 	gE <- gev.fit(as.numeric(E), show=F)
# 	plot(E, type="l", ylab="environment")
# 	hist(E, main=paste('GEV.sh = ',.001*round(1000*gE$mle[3])),30)
# }
# mtext("1-1 Hassell", side=3, outer=TRUE, line=-0.5)
#
#
#
# # =====================================
# # = Cushing et al. 2001 (Not in supp) =
# # =====================================
# sd <- .2
# expon <- 6
#
# b <- 10.45
# muL <- .2
# cL <- .01731
# cA <- .0131
#
# n <- 1000
#
# muA <- .96
# cP <- .35
#
#
# s <- .8
# cP <- 0
#
# L <- 10
# P <- 10
# A <- 10
#
# X <- array(0, dim=n+100)
# E <- array(0, dim=n+100)
# for(t in 1:(n+100)){
# 	# fat-tailed residuals
# 	if(sd != 0) {
# 		e <- rlnorm(1,sdlog=sd)^expon
# 		E[t] <- e
#
# 		nextL <- b*A*exp(-cA*A-cL*L)
# 		nextP <- (1-muL)*L
# 		nextA <- e*s*P*exp(-cP*A)+(1-muA)*A
# 		L <- nextL
# 		P <- nextP
# 		A <- nextA
# 	}else{
# 		nextL <- b*A*exp(-cA*A-cL*L)
# 		nextP <- (1-muL)*L
# 		nextA <- s*P*exp(-cP*A)+(1-muA)*A
# 		L <- nextL
# 		P <- nextP
# 		A <- nextA
# 	}
# 	X[t] <- L
# }
# X <- X[101:(n+100)]
# E <- E[101:(n+100)]
#
# gX <- gev.fit(as.numeric(X), show=F)
#
# dev.new()
# par(mfrow=c(3,2), oma=c(0.1,0.1,1,0.1))
# if(sd == 0) par(mfrow=c(2,2))
#
# plot(X, type="l", xlab="Time")
# hist(X, main=paste('GEV.sh = ',.001*round(1000*gX$mle[3])),30)
#
# w <- arima(X, order=c(3,0,3))
# Y <- w$residuals
#
# gY <- gev.fit(as.numeric(Y), show=F)
# plot(Y, type="l", ylab="residuals")
# hist(Y, main=paste('GEV.sh = ',.001*round(1000*gY$mle[3])),30)
#
# if(sd != 0){
#
# 	gE <- gev.fit(as.numeric(E), show=F)
# 	plot(E, type="l", ylab="environment")
# 	hist(E, main=paste('GEV.sh = ',.001*round(1000*gE$mle[3])),30)
# }
# mtext("Cushing 2001", side=3, outer=TRUE, line=-0.5)
