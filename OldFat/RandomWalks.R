rm(list=ls())
graphics.off()

#Druankard's Walk
StepOver <- ifelse(round(runif(8, 0, 1), 0)==1,1,-1)
StepUp <- ifelse(round(runif(8, 0, 1), 0)==1,1,-1)
Start <- 0
Loc <- data.frame("X"=0,"Y"=0)

for(i in 2:9){
	Loc[i,1] <- Loc[i-1,1] + StepOver[i-1]
	Loc[i,2] <- Loc[i-1,2] + StepUp[i-1] 
	}
	
dev.new(height=8, width=8)
par(family="Times", las=1, cex=2, mar=c(1,1,1,1), oma=c(1,1,1,1))
plot(Loc[,1],Loc[,2], pch=as.character(1:9), type="b", lty=3, ylim=c(-6,6), xlim=c(-6,6), ann=F, labels=F)
mtext("0", cex=2, side=1, line=1)
mtext("0", cex=2, side=2, line=1)

#Brownian Motion
Value <- c(0, rep(NA,999))
Steps <- rnorm(1000, mean=0, sd=10)
for(i in 2:1000){
	Value[i] <- Value[i-1] + Steps[i]
	}

dev.new(height=8, width=8)
par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
plot(1:1000, Value, type="l")
mean(Steps)
dev.new(height=8, width=8)
par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
hist(Steps)

#White Noise
dev.new(height=8, width=8)
par(family="Times", las=1, cex=1.5, mar=c(2,2,1,1), oma=c(1,1,1,1))
plot(1:1000, Steps, type="l")