

load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/fatPlot_Functions.R")



# ==============================================
# = Compare Waiting Times using significant Xi =
# ==============================================
dTimes <- c("Level2_time","Level2_normTime","Level2_logNormTime")
boundLog <- data.fat[,"shape.sig"]< 0
bound <- data.fat[boundLog,]
fatLog <- data.fat[,"shape.sig"] > 0
fat <- data.fat[fatLog,]
thinLog <- data.fat[,"shape.sig"] == 0
thin <- data.fat[thinLog,]

cLine <- rainbow(n=3, v=0.8, s=1)
cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)

# dev.new(width=3.5, height=6)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/Fig3_compareTimes.png", width=3.5, height=6, units="in", res=300, bg="white")
par(mfrow=c(3,1), mar=c(0.5, 1.25, 0.25, 0.25), mgp=c(1, 0.3, 0), tcl=-0.3, cex=1, ps=10, family="Times", oma=c(1.5,1,0,0), xpd=T)
pDens(bound, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Bounded~(xi<0)))
text(8, c(0.4, 0.5, 0.6), c("GEV", "Normal", "Log-Normal"), col=cLine)
pDens(thin, xaxt="n", xlab="")
axis(1, labels=FALSE)
text(8, 0.95, bquote(Thin~(xi==0)))
pDens(fat, xaxt="n", xlab="")
xat <- axTicks(1)
xlab <- parse(text=paste(10,xat,sep="^"))
axis(side=1, at=xat, labels=xlab)
text(8, 0.95, bquote(Fat~(xi>0)))
mtext("Waiting Time (years)", side=1, outer=TRUE, line=0.5)
mtext("Relative Density", side=2, line=0.0, outer=TRUE)
dev.off()



# ==============================================
# = Compare Waiting Times across Distributions =
# ==============================================
# dTimes <- c("Level2_time","Level2_normTime","Level2_logNormTime")
# cutoff <- 0.35
# boundLog <- data.fat[,"sh_0"]<= -cutoff
# bound <- data.fat[boundLog,]
# fatLog <- data.fat[,"sh_0"] > cutoff
# fat <- data.fat[fatLog,]
# thinLog <- data.fat[,"sh_0"] > -cutoff & data.fat[,"sh_0"] < cutoff #!boundLog & !fatLog
# thin <- data.fat[thinLog,]
# 
# cLine <- rainbow(n=3, v=0.8, s=1)
# cFill <- rgb(t(col2rgb(cLine, alpha=TRUE)), alpha=40, maxColorValue=255)
# 
# # dev.new(width=3.5, height=6)
# png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/Fig3_compareTimes.png", width=3.5, height=6, units="in", res=300, bg="white")
# par(mfrow=c(3,1), mar=c(0.5, 1.25, 0.25, 0.25), mgp=c(1, 0.3, 0), tcl=-0.3, cex=1, ps=10, family="Times", oma=c(1.5,1,0,0), xpd=T)
# pDens(bound, xaxt="n", xlab="")
# axis(1, labels=FALSE)
# text(8, 0.95, bquote(Bounded~(xi<=-.(cutoff))))
# text(8, c(0.4, 0.5, 0.6), c("GEV", "Normal", "Log-Normal"), col=cLine)
# pDens(thin, xaxt="n", xlab="")
# axis(1, labels=FALSE)
# text(8, 0.95, bquote(Thin~(-.(cutoff)<xi~phantom()<.(cutoff))))
# pDens(fat, xaxt="n", xlab="")
# xat <- axTicks(1)
# xlab <- parse(text=paste(10,xat,sep="^"))
# axis(side=1, at=xat, labels=xlab)
# text(8, 0.95, bquote(Fat~(xi>=.(cutoff))))
# mtext("Waiting Time (years)", side=1, outer=TRUE, line=0.5)
# mtext("Relative Density", side=2, line=0.0, outer=TRUE)
# dev.off()



