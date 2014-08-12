


source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatBirth_Sim.RData")

# set Par Margins (pm)
pm.1 <- c(0.75, 0.75, 0.5, 0)
pm.2 <- c(1.5, 0.75, 0.5, 0)
pm.3 <- c(1, 0.75, 0.5, 0)

# dev.new(width=3.5, height=5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Fig1_fat_conceptFig.png", res=150, units="in", height=5, width=3.5)
cols1 <- rep(rep(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 7, 7, 7, 7), each=4), 4) # first set of columns for layout matrix
cols2 <- rep(rep(c(1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 7, 7, 7, 7), each=4), 1) # first set of columns for layout matrix
cols3 <- rep(rep(c(4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 7, 7, 7, 7), each=4), 3) # second set of column for layout matrix
lmat <- matrix(c(cols1, cols2, cols3), ncol=8) # create layout matrix
layout(lmat) # define graphical device layout
par(mar=pm.1, oma=c(0.75, 1, 0.1, 0.1), ps=8, cex=1, mgp=c(0.25, 0.0, 0), tcl=0.15, family="Times")

	# ==========
	# = part 0 =
	# ==========
# part 0: annual maxima of environmental variables, birth/ survival rates, and population size
samples2.max <- matrix(c(apply(samples2, 2, yearly.Max2, n.per.year=samples.per.year)), nrow=years)
colnames(samples2.max) <- paste(rep(c("phi1","phi2","rB","rS","X"), each=2), rep(c(".i",""), 5), sep="")

	# ==========
	# = part 1 =
	# ==========
# part 1: time series of all samples of environmental variables, phi1 and phi2
ylim1.1 <- range(samples2[,"phi1"])*c(1, 1.15)
par(mar=pm.1)
plot(samples2[,"phi1"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim1.1) # plot full samples of variable phi1
points(samples2.max[,"phi1.i"], samples2.max[,"phi1"], col="blue", xpd=TRUE) # add circles for annual maxima of phi1
mtext(bquote(phi), side=2, line=0.75)
text(x=25, y=0.85, "A", font=2)

par(mar=pm.1)
ylim1.2 <- range(samples2[,"phi2"])*c(1, 1.15)
plot(samples2[,"phi2"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE, ylim=ylim1.2) # plot full samples of variable phi2
points(samples2.max[,"phi2.i"], samples2.max[,"phi2"], col="red", xpd=TRUE) # add circles for annual maxima of phi2
mtext(bquote(gamma), side=2, line=0.75)
text(x=25, y=0.85, "B", font=2)

	# ==========
	# = part 2 =
	# ==========
# part 2: time series of all samples of the population size, X
par(mar=pm.2)
plot(samples2[,"X"], type="l", col="gray", xlab="", ylab="", bty="l", xpd=TRUE) # plot full samples of variable phi1
points(samples2.max[,"X.i"], samples2.max[,"X"], col="black", xpd=TRUE) # add circles for annual maxima of phi1
mtext(bquote(X), side=2, line=0.75)
mtext(bquote(Time), side=1, line=0.75)
text(x=25, y=115, "C", font=2)



	# ==========
	# = part 3 =
	# ==========
# part 3: Parent and Maxima distributions for parts 1-3
# part 3.1
par(mar=pm.1)
colorDens(vals=list(samples2[,"phi1"], samples2.max[,"phi1"]), cols=c("gray","blue"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1.1)
cm1 <- min(samples2[,"phi1"])
text(y=-0.65, x=0.85*sum(range(density(samples2.max[,"phi1"])$y)), "D", font=2)

par(mar=pm.1)
colorDens(vals=list(samples2[,"phi2"], samples2.max[,"phi2"]), cols=c("gray","red"), revxy=TRUE, yaxt="n", bty="n", limX=ylim1.2)
cm1 <- min(samples2[,"phi2"])
text(y=-0.65, x=0.85*sum(range(density(samples2.max[,"phi2"])$y)), "E", font=2)

# part 3.2
par(mar=pm.2)
colorDens(vals=list(samples2[,"X"], samples2.max[,"X"]), cols=c("gray","black"), revxy=TRUE, yaxt="n", bty="n", limX=range(samples2[,"X"]))
cm1 <- min(samples2[,"X"])
X.dens <- density(samples2[,"X"])$x
X.max.dens <- density(samples2.max[,"X"])$y
text(y=25, x=0.85*sum(range(X.max.dens)), "F", font=2)
mtext(bquote(Density), side=1, line=0.75)

	# ==========
	# = Part 4 =
	# ==========
# part 4: difference in xi's between population and the environment
par(mar=pm.3)
colorDens(vals=list(xis[,1], xis[,2], xis[,3], (xis[,3]-(xis[,1]+xis[,2]))), cols=c("blue","red","black","green"), revxy=FALSE, yaxt="s", bty="l")
mtext(bquote(xi), side=1, line=0.75)
mtext(bquote(Density), side=2, line=0.75)
text(y=2.5, x=-0.75, "G", font=2)

dev.off()





