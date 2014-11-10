


source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/JDD.R")
source("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/FatTails_Functions.R")




# ==============================
# = Simulations of both models =
# ==============================
n.maxima <- 50
n.obs.per.block <- seq(10, 1E3, length.out=199) # I did 199 originally

# Big simulation involving 50 simulations of each of 199 time series that individually range between 10 and 1E3 time steps
set.seed(439)
jdd.johan.sims <- mapply(function(x)apply(mapply(jdd.johan, rep(x, n.maxima)), 2, max), n.obs.per.block)
jdd.johan.xi <- apply(jdd.johan.sims, 2, function(x)gev.fit(x)$mle["sh_0"])
jdd.johan.xi.log <- apply(jdd.johan.sims, 2, function(x)gev.fit(log(x))$mle["sh_0"])

set.seed(439)
jdd.brigo.sims <- mapply(function(x)apply(mapply(jdd.brigo, rep(x, n.maxima)), 2, max), n.obs.per.block)
jdd.brigo.xi <- apply(jdd.brigo.sims, 2, function(x)gev.fit(x)$mle["sh_0"])
jdd.brigo.xi.log <- apply(jdd.brigo.sims, 2, function(x)gev.fit(log(x))$mle["sh_0"])

# Small simulation of a single time series (50 years of 500 observations per year)
set.seed(2)
jdd.johan.1ts <- jdd.johan(N=500)

set.seed(4)
jdd.brigo.1ts <- jdd.brigo(N=500)





# ===============================
# = Plotting Simulation Results =
# ===============================
# Create hot to cold color scheme
ramp.cols <- col2rgb(colorRampPalette(c("blue", "cyan", "green", "yellow", "orange", "red"))(length(n.obs.per.block)), alpha=FALSE) # create color ramp
mycols <- rgb(red=ramp.cols["red",], green=ramp.cols["green",], blue=ramp.cols["blue",], alpha=50, maxColorValue=256) # add transparency

# dev.new(width=3.9, height=6) # new graphical device
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Other/JDD_2mods.png", width=3.9, height=6, units="in", res=150) 
# par(mfcol=c(4, 2), mar=c(2.25, 3, 0.5, 1.5), oma=c(0.1, 0.1, 8, 0.1), ps=10, mgp=c(0.9, 0.15, 0), tcl=-0.20, family="Times", cex=1) # graphical params
par(mfcol=c(4, 2), mar=c(1.9, 1.75, 0.5, 1.6), oma=c(0.0, 0.01, 5.5, 0.01), ps=8, mgp=c(0.65, 0.1, 0), tcl=-0.20, family="Times", cex=1)
# Plotting Johannes 2004 Eq 1
plot(jdd.johan.1ts, type="l", xlab="time (within 1 year)", ylab="Change (r)")
abline(h=0)

# Add Johannes equation
mtext(bquote( # begin mtext of equations
	atop({
				{r[t+Delta]-r[t]} == {mu*(r[t])*Delta+sigma*(r[t])*epsilon[t+Delta]*sqrt(Delta)+J[t+Delta]*Z[t+Delta]}
			},
			{
				{J[t+Delta]*Z[t+Delta]} == {sum({N*(0*','*sigma[z])[i]}, i==1, {Pois*(lambda*Delta)})}
			}
	))
, line=2.5)
mtext("Model for Change", line=5.25, font=2) # mtext for model label

# Plot time series of max returns
jdd.johan.ylim <- log(range(jdd.johan.sims)) # calculate ylims for plot of maxima time series
plot(log(jdd.johan.sims[,length(n.obs.per.block)]), type="l", ylim=jdd.johan.ylim, col=mycols[1], xlab="time (across 50 years)", ylab="log Max Change (r)") # plot first maxima
for(i in (length(n.obs.per.block)-1):1){ # loop through other maxima, plotting with lines()
	lines(log(jdd.johan.sims[,i]), type="l", col=mycols[i])
}
legend("topleft", legend=c("N=10", "N=1000"), text.col=c("blue", "red"), cex=1, bty="n", inset=c(-0.2, -0.2), xjust=0, horiz=TRUE) # add color scheme legend


# scatter of xi vs. number of observations used per each maximum (x axis is analagous to # obs per year when calculating xi from annual maxima)
plot(n.obs.per.block, jdd.johan.xi, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~r~~Maxima))
abline(lm(jdd.johan.xi~n.obs.per.block), lwd=3) # add solid black regression line
abline(lm(jdd.johan.xi~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.johan.xi), lty="dashed", col="blue", lwd=4) # add dotted blue median line
abline(h=median(jdd.johan.xi), lty="dashed", col="lightblue", lwd=2) # add dotted blue median line

# same as above, but for log(maxima)
plot(n.obs.per.block, jdd.johan.xi.log, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~log[e]*(r)~~Maxima))
abline(lm(jdd.johan.xi.log~n.obs.per.block), lwd=3)
abline(lm(jdd.johan.xi.log~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.johan.xi.log), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.johan.xi.log), lty="dashed", col="lightblue", lwd=2)


# Plotting Brigo 2007 Code 6
# Most of this code is same as above – only commenting on differences
# Plotting Johannes 2004 Eq 1
par(mar=c(2, 1.85, 0.5, 1.75))
plot(jdd.brigo.1ts, type="l", xlab="time (within 1 year)", ylab="State (S)", ylim=c(-0.1, max(jdd.brigo.1ts)))
abline(h=0)

# Add Brigo Equation
mtext(bquote(
	atop(
		{
			{r[t+Delta]-r[t]} == {mu*(r[t])*Delta+sigma*(r[t])*epsilon[t+Delta]*sqrt(Delta)+J[t+Delta]*Z[t+Delta]}
		},
		{
			{J[t+Delta]*Z[t+Delta]} == {sum({N*(0*','*sigma[z]*sqrt(n*(t+Delta)))[i]}, i==1, {n*(t+Delta)})} # note that this equation is slightly different from previous (multiply sigma[z] by sqrt(n(t+Delta)))
		}
		)
	)
, line=2.5)
mtext(bquote( # only difference in plotting is that I need 2 more equations, which are inserted here
	atop(
		{
			{n*(t+Delta)}=={Pois*(lambda*Delta)}
		},
		{
			{S*(T)} == {prod({exp*(r[t]*Delta)}, t==1, T)}
		}
		)
	)
, line=0)
mtext("Model for State", line=5.25, font=2)

# Plot max prices
jdd.brigo.ylim <- log(range(jdd.brigo.sims))
plot(log(jdd.brigo.sims[,length(n.obs.per.block)]), type="l", ylim=jdd.brigo.ylim, col=mycols[1], xlab="time (across 50 years)", ylab="log Max State (S)") # I think that this is the Price, and others are Returns
for(i in (length(n.obs.per.block)-1):1){
	lines(log(jdd.brigo.sims[,i]), type="l", col=mycols[i])
}
legend("topleft", legend=c("N=10", "N=1000"), text.col=c("blue", "red"), cex=1, bty="n", inset=c(-0.2, -0.2), xjust=0, horiz=TRUE)

# Scatter plots
plot(n.obs.per.block, jdd.brigo.xi, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~S~~Maxima))
abline(lm(jdd.brigo.xi~n.obs.per.block), lwd=3)
abline(lm(jdd.brigo.xi~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.brigo.xi), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.brigo.xi), lty="dashed", col="lightblue", lwd=2)

plot(n.obs.per.block, jdd.brigo.xi.log, xlab="N Obs per Maximum", ylab=bquote(xi~~of~~log[e]*(S)~~Maxima)) 
abline(lm(jdd.brigo.xi.log~n.obs.per.block), lwd=3)
abline(lm(jdd.brigo.xi.log~n.obs.per.block), lwd=1, col="white")
abline(h=median(jdd.brigo.xi.log), lty="dashed", col="blue", lwd=4)
abline(h=median(jdd.brigo.xi.log), lty="dashed", col="lightblue", lwd=2)
dev.off() # use this if you want to save the png

