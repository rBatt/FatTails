rm(list=ls())
graphics.off()


Plump <- c(0, rep(NA,999))
Slim <- c(0, rep(NA,999))

F1 <- round(runif(1000,0,1),0)
F2 <- round(runif(1000,0,1),0)
F3 <- round(runif(1000,0,1),0)
Steps <- rnorm(1000, mean=0, sd=1)

for(i in 2:10000){
	Plump[i] <-  F1[i]*F2[i]*F3[i]*Steps[i] + F1[i]*F2[i]*Steps[i] + F1[i]*F3[i]*Steps[i] + F2[i]*F3[i]*Steps[i] + F1[i]*Steps[i] + F2[i]*Steps[i] + F3[i]*Steps[i]
	Slim[i] <- F1[i]*Steps[i] + F2[i]*Steps[i] + F3[i]*Steps[i]
	}

Plump <- F1*F2*F3*Steps + F1*F2*Steps + F1*F3*Steps + F2*F3*Steps + F1*Steps + F2*Steps +F3*Steps
#Slim <- F1*Steps + F2*Steps + F3*Steps

Slim <- apply(cbind(F1,F2,F3), 2, max)*Steps

dev.new(height=8, width=8)
par(family="Times", las=1, cex=1.5, mar=c(4,4,1,1), oma=c(0,0,0,0))
plot.density(density(Plump), main="", ylab="Density", xlab="Value of 'Plump'", bty="l")

dev.new(height=8, width=8)
par(family="Times", las=1, cex=1.5, mar=c(4,4,1,1), oma=c(0,0,0,0))
plot.density(density(Slim), main="", ylab="Density", xlab="Value of 'Slim'", bty="l")


MeanExcess <- function(x){
	T  <- x #seq(min(x), max(x), length.out=200)
	MEs <- rep(NA, length(T))
	for(i in 1:length(T)){
		MEs[i] <- mean(x[which(x>T[i])]-T[i])
	}
	ME <- data.frame("T"=T, "MEs"=MEs)
	return(ME)
}

StepsME <- MeanExcess(Steps)
PlumpME <- MeanExcess(Plump)
SlimME <- MeanExcess(Slim)

Fatty <- rweibull(n=500, shape=1, scale=1) #rcauchy(1000, location=-2, scale=1)

meplot(Fatty, omit=0)
points(MeanExcess(Fatty), col="red")


