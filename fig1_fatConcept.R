
# =========================
# = ARIMA Simulation Plot =
# =========================
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatARMA_Sim.R")
	# ========================
	# = Prepare for plotting =
	# ========================
# Calculate densities for GEV fit (Panel E)
adat <- c(tfTS[tmTS], ffTS[fmTS]) # adat is 'all data' – both thin and fat maxima time series
ma <- min(adat) # find the smallest value of the maxima
dSeq <- seq(ma-sign(ma)*1.5*ma, max(adat)*1.1, by=0.1) # create a sequence of values over which to calculate the density
fatGEV <- dgev(dSeq, xi=simXiS[fattestI,"Xi"], simXiS[fattestI,"mu"], simXiS[fattestI,"sig"]) # probs for fat time series using GEV fit (package 'evir')
fatGEV[!is.finite(fatGEV)] <- 0 
thinGEV <- dgev(dSeq, xi=simXiS[thinnestI,"Xi"], simXiS[thinnestI,"mu"], simXiS[thinnestI,"sig"]) # probs for thin time series
thinGEV[!is.finite(thinGEV)] <- 0

# Prepare the Xi values to be plotted on figure (Panel E)
tXi <- round(simXiS[thinnestI,"Xi"], 2)
fXi <- round(simXiS[fattestI,"Xi"], 2)
tLabXi <- parse(text=paste("xi", tXi, sep=" = "))
fLabXi <- parse(text=paste("xi", fXi, sep=" = "))

	# ===============================
	# = Plot Example w/ reversed xy =
	# ===============================
# Set up figure space
dev.new(width=3.5, height=5) # open graphical device
# png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/fat_conceptFig.png", res=150, units="in", height=5, width=3.5)
cols1 <- rep(rep(1:3, each=3), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(4,5,3), each=3), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2), ncol=7) # create layout matrix
layout(lmat) # define graphical device layout
par(mar=c(2.5,0.0,0.5,0.1), oma=c(0, 2, 0, 0.25), ps=10, cex=1, mgp=c(1, 0.3, 0), tcl=-0.25, family="Times") # set graphical parameters

# Plot thin-tailed time series
ylim1 <- range(tfTS)*c(1, 1.15)
plot(tfTS, type="l", col="gray", xlab="", xaxt="n", ylab="", ylim=ylim1, bty="l")
text(0.025*length(tfTS), y=max(tfTS)*1.05, "A", font=2)
ats1 <- axTicks(1)
labs1 <- ats1/nPerYear
axis(side=1, labels=labs1, at=ats1)
points(tmTS, tfTS[tmTS], col="blue")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot fat-tailed time series
ylim2 <- range(ffTS)*c(1, 1.15)
plot(ffTS, type="l", col="gray", xlab="", ylab="", xaxt="n", ylim=ylim2, bty="l")
text(0.025*length(ffTS), y=max(ffTS)*1.05, "B", font=2)
ats2 <- axTicks(1)
labs2 <- ats2/nPerYear
axis(side=1, labels=labs2, at=ats2)
points(fmTS, ffTS[fmTS], col="red")
mtext("y", side=2, line=1.25)
mtext("time", side=1, line=1.25)

# Plot GEV densities
colorPoly(quants=dSeq, dents=list(thinGEV,fatGEV), cols=c("blue", "red"), bty="l")
cm0 <- min(dSeq)
text(x=cm0+sign(cm0)*cm0*0.05, y=max(c(thinGEV,fatGEV))*0.95, "E", font=2)
mtext("density", side=2, line=1.25)
mtext("y", side=1, line=1.25)
cma0 <- max(dSeq)
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.75, bquote(xi~"="~.(tXi)), pos=4, col="blue")
text(x=cma0*0.5, y=max(c(thinGEV,fatGEV))*0.58, bquote(xi~"="~.(fXi)), pos=4, col="red")

# Plot empirical densities for thin-tailed
colorDens(vals=list(tfTS, tfTS[tmTS]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1)
cm1 <- min(tfTS)
text(y=cm1+sign(cm1)*cm1*0.15, x=0.25*max(density(tfTS)$y, density(tfTS[tmTS])$y), "C", font=2)
mtext("density", side=1, line=1.25)

# Plot empirical densities for fat-tailed
cm2 <- min(ffTS)
colorDens(vals=list(ffTS, ffTS[fmTS]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim2)
text(y=cm2+sign(cm2)*cm2*0.15, x=0.25*max(density(ffTS)$y, density(ffTS[fmTS])$y), "D", font=2)
mtext("density", side=1, line=1.25)
# dev.off()
