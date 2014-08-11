
# ===============
# = Ryan Figure =
# ===============
library("evir")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")

pm.1 <- c(0.75, 0, 0.5, 0)
pm.blank <- c(0, 0, 1, 0)
pm.2 <- c(0.75, 0, 0.5, 0)
pm.3 <- c(0.75, 0, 0.5, 0)

dev.new(width=3.5, height=7)
# cols1 <- rep(rep(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7), each=4), 4) # first set of columns for layout matrix
# cols2 <- rep(rep(c(8, 8, 8, 9, 9, 9, 0, 0, 10, 10, 10, 10, 0, 0, 11, 11, 11, 11, 7, 7, 7, 7), each=4), 3) # second set of column for layout matrix
cols1 <- rep(rep(c(1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5, 6, 6, 7, 7, 7, 7, 8, 8, 8, 8), each=4), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(9, 9, 9, 10, 10, 10, 11, 11, 12, 12, 12, 13, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 16, 16), each=4), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2), ncol=7) # create layout matrix
layout(lmat) # define graphical device layout
# par(mfcol=c(6,2), mar=c(1,2,0.5,0.5), ps=8, cex=1, mgp=c(0.25, 0.0, 0), tcl=0.15)
par(mar=pm.1, oma=c(0.75, 1.75, 0.1, 0.1), ps=8, cex=1, mgp=c(0.25, 0.0, 0), tcl=0.15)

	# ==========
	# = part 0 =
	# ==========
# part 0: annual maxima of environmental variables, birth/ survival rates, and population size
samples2.max <- matrix(c(apply(samples2, 2, yearly.Max2, n.per.year=samples.per.year)), nrow=years)
colnames(samples2.max) <- paste(rep(c("phi1","phi2","rB","rS","X"), each=2), rep(c(".i",""), 5), sep="")

	# ============
	# = part 0.5 =
	# ============
# part 0.5: Calculate densities for GEV fit (Panel E)
adat <- c(apply(samples2.max[,c("phi1","phi2","X")], 2, scale)) # adat is 'all data' – both thin and fat maxima time series
ma <- min(adat) # find the smallest value of the maxima
dSeq <- seq(ma, max(adat)*1.5, by=0.1) # create a sequence of values over which to calculate the density
# dSeq <- seq(2, 7, by=0.1) # create a sequence of values over which to calculate the density

phi1.gev <- gev.fit2(scale(samples2.max[,"phi1"]))$mle
phi1GEV <- dgev(dSeq, xi=phi1.gev[3], phi1.gev[1], phi1.gev[2]) # probs for fat time series using GEV fit (package 'evir')
phi1GEV[!is.finite(phi1GEV)] <- 0 

phi2.gev <- gev.fit2(scale(samples2.max[,"phi2"]))$mle
phi2GEV <- dgev(dSeq, xi=phi2.gev[3], phi2.gev[1], phi2.gev[2]) # probs for thin time series
phi2GEV[!is.finite(phi2GEV)] <- 0

X.gev <- gev.fit2(scale(samples2.max[,"X"]))$mle
XGEV <- dgev(dSeq, xi=X.gev[3], X.gev[1], X.gev[2]) # probs for thin time series
XGEV[!is.finite(XGEV)] <- 0

	# ==========
	# = part 1 =
	# ==========
# part 1: time series of all samples of environmental variables, phi1 and phi2
ylim1.1 <- range(samples2[,"phi1"])*c(1, 1.15)
par(mar=pm.1)
plot(samples2[,"phi1"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim1.1) # plot full samples of variable phi1
points(samples2.max[,"phi1.i"], samples2.max[,"phi1"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1
mtext(bquote(phi1[1]), side=2, line=0.75)

par(mar=pm.1)
ylim1.2 <- range(samples2[,"phi2"])*c(1, 1.15)
plot(samples2[,"phi2"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim1.2) # plot full samples of variable phi2
points(samples2.max[,"phi2.i"], samples2.max[,"phi2"], col="red", xpd=TRUE) # add circles for annual maxima of phi2
mtext(bquote(phi1[2]), side=2, line=0.75)


par(mar=pm.blank)
plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region
# text(1, 1, labels=quote(atop({phantom() %down% phantom()}, {B(phi1[1])*","~S(phi1[2])})), cex=2)
text(1, 1, labels=quote(phantom() %dbldown% phantom()), cex=3, xpd=TRUE) # add double arrow to visually indicate that part 1 gives rise to part 2

	# ==========
	# = part 2 =
	# ==========
# part 2: time series of all samples of birth and survival rates, B(phi1) and S(phi2)
# pt2.ylim <- range(samples2[,c("rB", "rS")])
par(mar=pm.2)
ylim2.1 <- range(samples2[,"rB"])*c(1, 1.15)
plot(samples2[,"rB"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim2.1) # plot full samples of variable phi1
points(samples2.max[,"rB.i"], samples2.max[,"rB"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1
mtext(bquote(B(phi1[1])), side=2, line=0.75)

par(mar=pm.2)
ylim2.2 <- range(samples2[,"rS"])*c(1, 1.15)
plot(samples2[,"rS"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim2.2) # plot full samples of variable phi2
points(samples2.max[,"rS.i"], samples2.max[,"rS"], col="red", xpd=TRUE) # add circles for annual maxima of phi2
mtext(bquote(S(phi1[2])), side=2, line=0.75)

par(mar=pm.blank)
plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region
# text(1, 1, labels=quote(atop({phantom() %down% phantom()}, {B(phi1[1])*","~S(phi1[2])})), cex=2)
text(1, 1, labels=quote(phantom() %dbldown% phantom()), cex=3, xpd=TRUE) # add double arrow to visually indicate that part 1 gives rise to part 2

# 
# ylim2 <- range(samples2[,"phi2"])*c(1, 1.15)
# plot(samples2[,"rB"], type="l", col="gray", xlab="", ylab="", bty="l", ylim=ylim2, xpd=TRUE) # plot full samples of variable phi1
# points(samples2.max[,"rB.i"], samples2.max[,"rB"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1
# 
# lines(samples2[,"rS"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi2
# points(samples2.max[,"rS.i"], samples2.max[,"rS"], col="red", xpd=TRUE) # add circles for annual maxima of phi2
# 
# mtext(bquote(B(phi1[1])*","~S(phi1[2])), side=2, line=0.75)
# 
# plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region
# text(1, 1, labels=quote(phantom() %dbldown% phantom()), cex=3, xpd=TRUE) # add double arrow to visually indicate that part 2 gives rise to part 3

	# ==========
	# = part 3 =
	# ==========
# part 3: time series of all samples of the population size, X
par(mar=pm.3)
plot(samples2[,"X"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi1
points(samples2.max[,"X.i"], samples2.max[,"X"], col="black", xpd=TRUE) # add circles for annual maxima of phi1
mtext(bquote(X), side=2, line=0.75)


	# ==========
	# = part 4 =
	# ==========
# part 4: GEV Distributions
phi1.xi <- round(phi1.gev[3], 2)
phi2.xi <- round(phi2.gev[3], 2)
X.xi <- round(X.gev[3], 2)
phi1.LabXi <- parse(text=paste("xi", phi1.xi, sep=" = "))
phi2.LabXi <- parse(text=paste("xi", phi2.xi, sep=" = "))
X.LabXi <- parse(text=paste("xi", X.xi, sep=" = "))

# par(xaxt="n", yaxt="n")
colorPoly(quants=dSeq, dents=list(phi1GEV,phi2GEV, XGEV), cols=c("blue", "red", "black"), bty="l")
# axis(side=2, labels=FALSE)
cm0 <- min(dSeq)
text(x=cm0+sign(cm0)*cm0*0.05, y=max(c(phi1GEV,phi2GEV, XGEV))*0.95, "I", font=2)
mtext("density", side=2, line=0.75)
mtext(bquote(standardized~phi1[1]*","~phi1[2]*","~or~X), side=1, line=0.75)
cma0 <- max(dSeq)
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.75, bquote(xi~"="~.(phi1.xi)), pos=4, col="blue")
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.55, bquote(xi~"="~.(phi2.xi)), pos=4, col="red")
text(x=cma0*0.5, y=max(c(phi1GEV,phi2GEV, XGEV))*0.35, bquote(xi~"="~.(X.xi)), pos=4, col="black")


	# ==========
	# = part 5 =
	# ==========
# part 5: Parent and Maxima distributions for parts 1-3
# part 5.1
par(mar=pm.1)
colorDens(vals=list(samples2[,"phi1"], samples2.max[,"phi1"]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1.1)
cm1 <- min(samples2[,"phi1"])
# text(y=cm1+sign(cm1)*cm1*0.5, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
text(y=0.25*sum(range(density(samples2[,"phi1"])$x)), x=0.75*sum(range(density(samples2.max[,"phi1"])$y)), "E", font=2)
mtext("density", side=1, line=0.5)

par(mar=pm.1)
colorDens(vals=list(samples2[,"phi2"], samples2.max[,"phi2"]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1.2)
cm1 <- min(samples2[,"phi2"])
# text(y=cm1+sign(cm1)*cm1*0.5, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
text(y=0.25*sum(range(density(samples2[,"phi2"])$x)), x=0.75*sum(range(density(samples2.max[,"phi2"])$y)), "F", font=2)
mtext("density", side=1, line=0.5)

par(mar=pm.blank)
plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region

# part 5.2
par(mar=pm.2)
colorDens(vals=list(samples2[,"rB"], samples2.max[,"rB"]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim2.1)
cm1 <- min(samples2[,"rB"])
# text(y=cm1+sign(cm1)*cm1*0.5, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
rB.dens <- density(samples2[,"rB"])$x
rB.max.dens <- density(samples2.max[,"rB"])$y
text(y=0.45*sum(range(rB.dens)), x=0.75*sum(range(rB.max.dens)), "G", font=2)
mtext("density", side=1, line=0.5)

par(mar=pm.2)
colorDens(vals=list(samples2[,"rS"], samples2.max[,"rS"]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim2.2)
cm1 <- min(samples2[,"rS"])
# text(y=cm1+sign(cm1)*cm1*0.5, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
rS.dens <- density(samples2[,"rS"])$x
rS.max.dens <- density(samples2.max[,"rS"])$y
text(y=0.45*sum(range(rS.dens)), x=0.75*sum(range(rS.max.dens)), "H", font=2)
mtext("density", side=1, line=0.5)

par(mar=pm.blank)
plot(1, type="n", xaxt="n", yaxt="n", xlab="", ylab="", bty="n") # plot blank region

# part 5.3
par(mar=pm.3)
colorDens(vals=list(samples2[,"X"], samples2.max[,"X"]), cols=c("gray","black"), revxy=TRUE, yaxt="n", bty="n", limX=range(samples2[,"X"]))
cm1 <- min(samples2[,"rS"])
# text(y=cm1+sign(cm1)*cm1*0.5, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
X.dens <- density(samples2[,"X"])$x
X.max.dens <- density(samples2.max[,"X"])$y
text(y=0.225*sum(range(X.dens)), x=0.85*sum(range(X.max.dens)), "I", font=2)
mtext("density", side=1, line=0.5)


# part 6: difference in xi's between population and the environment
hist(xis[,3] - (xis[,1]+xis[,2]), main='Pop-Envir', xlim=c(-1,1), breaks=seq(-1,1,0.1))



