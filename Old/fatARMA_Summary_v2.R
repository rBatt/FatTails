#RDB
#_v0 (01-Dec-2013): Summarize the ARMA analysis of fatFrame by 1) selecting the best model by AICc; 2) getting the AICc-weighted average of the eigenvalue. Also, this might be a good place to begin comparing to the Fat Tails analysis.

# _v2 (13-Jan-2014): Making some corrections to the computation of the sigmas (changes in tonyARMA_short)

rm(list=ls())
graphics.off()

library("plyr")
library("rpart")
library("party")
library("RColorBrewer")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("fatARMA_v1.RData") #this is the data file containing the completed ARMA analysis. Note that _v2 is the same as _v1, because tonyARMA_short _v4 and _v5 compute the ARMA the same, but differ in the way they compute the sigma metrics. fatARMA_vX.RData only contains the ARMA fit, not the sigma metrics.
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
load("finalFrame_v2.RData")

source("FatTails_Functions_v7.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("TonySuggestions/tonyARMA_short_v5.R") #also loads GenSA and DEoptim packages
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")

eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))

fWeighted <- function(x){
	#http://www.ssc.wisc.edu/~bhansen/718/NonParametrics15.pdf
	# more accurate at: http://machinelearning102.pbworks.com/w/file/fetch/47699411/aic_reg.pdf
	
	finiteLogic <- all(!is.finite(x[,"AICc"])) #If all of the AICc's are infinite
	missLogic <- all(is.na(x[,"AICc"])) #If all of the AICc's are missing
	if(finiteLogic | missLogic){  #don't compute the AICc-weighted averages (one of these may be true depending on wether I first converted Inf to NA or not)
		return(data.frame(x, "wLambda"=NA, "wOrder"=NA))
	}else{
		minAIC <- min(x[,"AICc"], na.rm=TRUE)

		eaic <- exp(-0.5*(x[,"AICc"]-minAIC))

		saicL <- sum(eaic[!is.na(x[,"Lambda"])], na.rm=TRUE)
		saicO <- sum(eaic[!is.na(x[,"Order"])], na.rm=TRUE)
		saicEps <- sum(eaic[!is.na(x[,"sigEps"])], na.rm=TRUE)
		saicE <- sum(eaic[!is.na(x[,"sigE"])], na.rm=TRUE)
		saicInf <- sum(eaic[!is.na(x[,"sigInf"])], na.rm=TRUE)

		wsL <- eaic/saicL
		wdO <- eaic/saicO
		wEps <- eaic/saicEps
		wE <- eaic/saicE
		wInf <- eaic/saicInf

		wLambda <- sum(wsL*x[,"Lambda"], na.rm=TRUE)
		wOrder <- sum(wdO*x[,"Order"], na.rm=TRUE)
		wSigEps <- sum(wEps*x[,"sigEps"], na.rm=TRUE)
		wSigE <- sum(wE*x[,"sigE"], na.rm=TRUE)
		wSigInf <- sum(wInf*x[,"sigInf"], na.rm=TRUE)

		return(data.frame(x, "wLambda"=wLambda, "wOrder"=wOrder, "wSigEps"=wSigEps, "wSigE"=wSigE, "wSigInf"=wSigInf))
	}
}

sigARMA <- ddply(fatARMA, .variables=c("variable", "location", "P", "Q"), .fun=getSE, data=finalFrame, .progress="time")

fatARMA1 <- merge(fatARMA, sigARMA, all=TRUE)

# metVars <- c("")
# physVars <- c("DaysOpen", "extcoef", "LakeLevel", "max_air_temp", "min_air_temp", "o2", "o2sat", "precip_mm", "range_air_temp", "Secchi")
# chemVars <- c("alk", "brsiuf", "ca", "cl", "cond", "dic", "doc", "drsif", "fe", "k", "mg", "mn", "na", "nh4", "no3no2", "ph",)
# bioVars <- c("avg_length", "avg_zoop_mass", "chlor", "cpue1_Sum", "cpue3_WeiEff", "density")

fatARMA2 <- Inf2NA(data.frame(fatARMA1, "Order"=fatARMA1[,"P"]+fatARMA1[,"Q"], "PQ"=paste("(", fatARMA1[,"P"], ",", fatARMA1[,"Q"], ")",sep="")))
fatARMA3 <- ddply(fatARMA2, .variables=c("variable", "location"), .fun=fWeighted)


bestARMA <- ddply(.data=fatARMA3, .variables=c("variable", "location"), .fun=minaicc)
names(bestARMA)[1:2] <- c("Variable", "fitBy")

# plot(bestARMA[,"Lambda"], bestARMA[,"wLambda"])

final0 <- merge(bestARMA,sAP, all=TRUE)
final0[,"Order"] <- as.factor(final0[,"Order"])
final0[,"P"] <- as.factor(final0[,"P"])
final0[,"Q"] <- as.factor(final0[,"Q"])
final0[,"Variable"] <- as.factor(final0[,"Variable"])
final0[,"fitBy"] <- as.factor(final0[,"fitBy"])
final0[,"Type"] <- as.factor(final0[,"Type"])

# ===========================================================================================
# = Need to look into these calculations to make sure they are correct in _v2 (13-Jan-2014) =
# ===========================================================================================
InfE <- final0[,"sigInf"]/final0[,"sigE"]
Einf <- (final0[,"sigE"])^2/(final0[,"sigInf"])^2

final0[,"InfE"] <- InfE
final0[,"Einf"] <- Einf


#interesting note, the reason TB alkalinity has NA for its mean and other "lvl" values is b/c it has some negative values, which means the logMean etc couldn't be computed, and I used complete.cases to remove "lvl" info from variables that have any NA's in the "lvl" categories. Actually, I'm not sure thsi could have been the case, because convNeg should take care of this.Ah, but it couldn't have more than 20% of the data set (after NA's and Inf removed) be negative.
final_logical <- complete.cases(final0[,c("Type","Variable","fitBy","P","Q","Lambda","InfE","Einf")])
final <- final0[final_logical,]
sum(!is.finite(fatARMA2[,"Period"]))

myWhite <- rgb(t(col2rgb(col="white")),alpha=100, maxColorValue=256)
myWhite2 <- rgb(t(col2rgb(col="white")),alpha=200, maxColorValue=256)


# =====================
# = Period Histrogram =
# =====================
# dev.new(width=3.4, height=3.4)
png("FatFigures_v2/periodHistogram.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
perHist_breaks <- hist(final[final[,"Period"]<150,"Period"], plot=FALSE)$breaks
hist(final[final[,"Period"]<150,"Period"], xlab="Period", ylab="Frequency", main="", col="white")
hist(final[final[,"Type"]=="Chem" & final[,"Period"]<150 ,"Period"], add=TRUE, col="gray", breaks=perHist_breaks)
hist(final[final[,"Type"]=="Bio" & final[,"Period"]<150,"Period"], add=TRUE, col="black", breaks=perHist_breaks)
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)
dev.off()

# =================
# = Order Barplot =
# =================
tally <- prop.table(table(final[,c("Type", "Order")]), margin=1)[c("Bio", "Chem","Phys","Met"),]
# dev.new(width=10, height=5)
# par(mfrow=c(2,4), mar=c(3,2,1,0.5), cex=1, ps=8, mgp=c(1,0.5,0))
# Names <- c("Accent", "Dark2", "Paired", "Pastel1", "Pastel2", "Set1", "Set2", "Set3")
# for(i in 1:length(Names)){
# 	barplot(tally, beside=TRUE, col=brewer.pal(n=4, name=Names[i]), axisnames=TRUE, main=Names[i], legend.text=c("Bio", "Chem","Phys","Met"), args.legend=list(bty="n"))
# }
# dev.new(width=3.4, height=3.4)
png("FatFigures_v2/orderBar.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(mar=c(2.5,2.5,0.5,0.5), cex=1, ps=9, mgp=c(1.5,0.5,0), family="Times")
barplot(tally, beside=TRUE, col=brewer.pal(n=4, name="Set1"), axisnames=TRUE, legend.text=c("Bio", "Chem","Phys","Met"), args.legend=list(bty="n"), xlab="Order", ylab="Proportion")
dev.off()




# =======================
# = Variance Parameters =
# =======================
# InfE <- final[,"sigInf"]/final[,"sigE"]
# Einf <- (final[,"sigE"])^2/(final[,"sigInf"])^2
# 
# wInfE <- final[,"wSigInf"]/final[,"wSigE"]
# wEinf <- final[,"wSigE"]/final[,"wSigInf"]


#a lot here is significant ...
summary(lm(log10(Einf)~sh_0+Lambda, data=final))

summary(lm(log10(Einf)~sh_0+Type+Lambda, data=final))

summary(lm(log10(Einf)~sh_0+Type+Lambda+log10(Level1_time), data=final))


# let's see what significance we can get with Xi as the response
summary(lm(sh_0~log10(Einf)+Lambda, data=final))

summary(lm(sh_0~log10(Einf)+Lambda+Type, data=final))#not so good

summary(lm(sh_0~log10(Einf)+Lambda+Type+log10(Level1_time), data=final))# real bad




library("lme4")
library("MuMIn")
meFrame0 <- data.frame("Type"=final[,"Type"],"Variable"=final[,"Variable"], "Order"=as.factor(final[,"Order"]), "P"=as.factor(final[,"P"]), "fitBy"=final[,"fitBy"], "Einf"=final[,"Einf"], "sh_0"=final[,"sh_0"], "Lambda"=final[,"Lambda"], "l10L1time"=log10(final[,"Level1_time"]), "InfE"=final[,"InfE"], l10Einf=log10(final[,"Einf"]), "l10InfE" = log10(final[,"InfE"]), "N"=final[,"N"])
meFrame <- meFrame0[complete.cases(meFrame0),]
meFrameMeans <- colMeans(meFrame[,6:13])
# meFrame[,4:9] <- apply(meFrame[,4:9], 2, scale, center=TRUE, scale=FALSE)
meFrame2 <- meFrame[meFrame[,"l10Einf"]>-5,]

# meMod10 <- lmer(log10(Einf)~sh_0+Lambda+(1|fitBy), data=meFrame)
# summary(meMod10)
# r.squaredGLMM(meMod10)

# meMod1 <- lmer(log10(Einf)~sh_0+(1|fitBy)+(sh_0|Variable)+Lambda, data=meFrame)
# summary(meMod1)
# r.squaredGLMM(meMod1)

summary(lm(l10Einf~sh_0+Lambda+Type, data=meFrame))
summary(lm(l10Einf~sh_0*Lambda+Type, data=meFrame))
meMod1 <- lmer(l10Einf~sh_0+(1|Type)+Lambda, data=meFrame)
summary(meMod1)
r.squaredGLMM(meMod1)

# meMod2 <- lmer(log10(InfE)~sh_0+(1|fitBy)+(sh_0|Variable)+Lambda, data=meFrame)
# summary(meMod2)
# r.squaredGLMM(meMod2)

summary(lm(sh_0~l10InfE, data=meFrame))
meMod2 <- lmer((sh_0)~l10InfE+Lambda+(l10InfE|Type), data=meFrame)
summary(meMod2)
r.squaredGLMM(meMod2)


# 
# #replicate estimates from the lmer model
# newV <- 0.3672*final[,"sh_0"] + -0.7369*final[,"Lambda"]
# dev.new(width=3.4, height=3.4)
# dscat((newV), log10(final[,"Einf"]), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# # dscat((newV)[log10(Einf)>-5], log10(Einf)[log10(Einf)>-5], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# mtext("Envi explained variance", side=2, line=1.5)
# mtext("newV", side=1, line=1.5)
# 
# #replicate estimates from the lm model that only had Lambda and sh_0
# newV <- 1.29104*final[,"sh_0"] + -1.30992*final[,"Lambda"]
# dev.new(width=3.4, height=3.4)
# # dscat((newV), log10(Einf), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# dscat((newV)[log10(final[,"Einf"])>-5], log10(final[,"Einf"])[log10(final[,"Einf"])>-5], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# mtext("Envi explained variance", side=2, line=1.5)
# # mtext("newV", side=1, line=1.5)
# mtext(bquote((1.3%*%xi)~~~-~~~(1.3%*%lambda)), side=1, line=1.5)

# 
# g0 <- ggplot(meFrame, aes(y=sh_0,x=l10Einf, group=Order))
# g1 <- g0 + facet_grid(.~Type, margin=TRUE, scales="free") + geom_point(aes(size=Order))
# g2 <- g1 + stat_smooth(method="lm")
# 
# g1 + continuous_scale( aes(is.element(Order, c(1,2))), scale_name="size") 
# 
# z0 <- ggplot(meFrame, aes(x=N, y=Lambda), group=Type)
# z1 <- z0 + geom_point()




# ==============================
# = Reproduce Ziebarth figures =
# ==============================


#Figure 5
sN <- final[,"N"]<50

c3 <- function(x, y, funct="mean", N=20, ...){
	require("zoo")
	n2 <- floor(sum(!is.na(x)&!is.na(y))/N)
	
	funct <- get(funct)
	ox <- order(x)
	x1 <- x[ox]
	y1 <- as.numeric(y[ox])
	fx <- rollapplyr(x1, width=n2, by=n2, FUN=funct, ...) 
	fy <- rollapplyr(y1, width=n2, by=n2, FUN=funct, ...)
	# fc <- rollapplyr(y1, width=n2, by=n2, FUN=function(x)sum(!is.na(x)))
	fc <- rollapplyr(y1, width=n2, by=n2, FUN=sd, na.rm=TRUE)
	z <- cbind(fx, fy, fc)
	return(z)
}
bNL <- c3(final[sN,"N"], final[sN,"Lambda"], N=19, na.rm=TRUE)
bNO <- c3(final[sN,"N"], final[sN,"Order"], N=19, na.rm=TRUE)



# dev.new(width=3.4, height=5)
png("FatFigures_v2/lambda_order_N.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(mfrow=c(2,1), mar=c(2,2,0.5, 0.5), mgp=c(1,0.1,0), tcl=0.25, ps=9, pch=19, bty="l")
plot(bNL[,1], bNL[,2], xlab="", ylab=lNota)
plot(bNO[,1], bNO[,2], xlab="Time series length", ylab=bquote(ave*.~p+q))
dev.off()

#Figure 4

# dev.new(width=6, height=3.4)
png("FatFigures_v2/lambda_variance_hist.png", width=6, height=3.4, units="in", res=300, bg=myWhite)

par(family="Times", mfrow=c(1,2), ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
hist(final[,"Lambda"], xlab=lNota2, ylab="Frequency", main="", col="white")
hist(final[final[,"Type"]=="Chem","Lambda"], add=TRUE, col="gray")
hist(final[final[,"Type"]=="Bio","Lambda"], add=TRUE, col="black")
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)

ieBreaks <- hist(log10(InfE), plot=FALSE)$breaks
hist(log10(InfE), xlab=bquote(log[10]~(.(ieNota))), ylab="Frequency", main="", col="white")
hist(log10(InfE)[final[,"Type"]=="Chem"], add=TRUE, col="gray")
hist(log10(InfE)[final[,"Type"]=="Bio"], add=TRUE, col="black", breaks=ieBreaks)

dev.off()

# g0 <- ggplot((data.frame(bNL)), aes(fx, fy))
# g0+geom_point(aes(size=fc))+scale_size(range=c(5,2)) + stat_smooth(aes(fx,fy),method="lm") + theme_classic()



#Figure 2
# dev.new(width=6, height=3.4)
png("FatFigures_v2/variance_lambda.png", width=6, height=3.4, units="in", res=300, bg=myWhite)
library("gridExtra")
library("ggplot2")

f2_xlb <- xlab(lNota2)
f2_theme <- theme_classic() + theme(legend.position="none", text=element_text(family="Times", size=10)) + theme(panel.background=element_rect(fill="transparent", colour=NA), plot.background=element_rect(fill="transparent", colour=NA))
f2_col <- scale_colour_brewer(palette="Dark2")

f20 <- ggplot(final, aes(x=Lambda, group=P)) + f2_theme

f21_fit <- stat_smooth(aes(y=log10(InfE), x=(Lambda), colour=P), formula=y~exp(x), fill=NA, na.value=NA, method="lm", size=1.1)
# f21_theme <- theme(legend.position=c(0.2, 0.8), plot.margin=unit(c(0.2,0.2,-3,2), "points")) 
f21_theme <- theme(legend.position=c(0.2, 0.8)) 
f21 <- f20 + geom_point(aes(y=log10(InfE), colour=P), alpha=0.5) + f21_fit + f2_col + f2_xlb + ylab(bquote(log[10](.(ieNota)))) + f21_theme

f22_fit <- stat_smooth(aes(y=(Einf), x=(Lambda), colour=P), formula=(y)~(x), fill=NA, na.value=NA, method="lm", size=1.1)
# f22_theme <- theme(plot.margin=unit(c(0.1,0.1,2,2), "points"))
f22 <- f20 + geom_point(aes(y=(Einf), colour=P), alpha=0.5) + f2_col + f2_xlb + ylab(eiNota) + f22_fit#  + scale_x_continuous(expand=c(0.001,0)) + scale_y_continuous()#+ f22_theme 

print(grid.arrange(f21, f22, ncol=2))
dev.off()

final <- transform(final, nCat=N>=25, Ext_mu_Rat=final[,"mu_0"]/final[,"mean"])

# ==========================
# = Explore patterns again =
# ==========================
# dev.new(width=8, height=3.4)
gg_sh_logic <- !is.element(final[,"Type"],c("Cosmic", "Met")) & !is.element(final[,"fitBy"],c("Formula", "KE", "LR", "WA")) & !is.na(final[,"N"])
# sh0 <- ggplot(final[gg_sh_logic,], aes(y=sh_0)) + f2_theme + theme(legend.position="right")
# sh1 <- sh0 + geom_point(aes(x=log10(InfE), color=Lambda), na.rm=TRUE, alpha=0.5) + scale_size(range=c(2, 6))
# # sh_guide <- guide_legend("N")
# sh_stat <- stat_smooth(aes(y=sh_0,x=log10(InfE)), fill=NA, na.value=NA, method="lm", formula=y~x, na.rm=TRUE, color="black")
# sh1 + facet_grid(.~Type, margins=TRUE, scales="free") + sh_stat #+ guides(size=sh_guide, colour=sh_guide)


# dev.new(width=8, height=3.4)
# ret0 <- ggplot(final[gg_sh_logic,], aes(y=log10(Level2_time))) + f2_theme + theme(legend.position="right")
# ret1 <- ret0 + geom_point(aes(x=log10(InfE), color=N, size=N), na.rm=TRUE, alpha=0.5) + scale_size(range=c(2, 6))
# ret_stat <- stat_smooth(aes(y=log10(Level2_time),x=log10(InfE)), fill=NA, na.value=NA, method="lm", formula=y~x, na.rm=TRUE, color="black")
# ret1 + facet_grid(.~Type, margins=TRUE, scales="free") + ret_stat

# dev.new(width=8, height=3.4)
# ret0 <- ggplot(final[gg_sh_logic,], aes(y=sh_0)) + f2_theme + theme(legend.position="right")
# ret1 <- ret0 + geom_point(aes(x=(Ext_mu_Rat), color=N, size=N), na.rm=TRUE, alpha=0.5) + scale_size(range=c(2, 6))
# ret_stat <- stat_smooth(aes(y=sh_0,x=(Ext_mu_Rat)), fill=NA, na.value=NA, method="lm", formula=y~x, na.rm=TRUE, color="black")
# ret1 + facet_grid(.~Type, margins=TRUE, scales="free") + ret_stat

# dev.new(width=3.4, height=3.4)
png("FatFigures_v2/Xi_InfE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","InfE"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(InfE), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=1, line=1.5)
dev.off()

# dev.new(width=3.4, height=3.4)
png("FatFigures_v2/Xi_Lambda.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(sqrt(final[gg_sh_logic&final[,"Type"]!="Phys","Lambda"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(Lambda), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.5, cex=1)
mtext(bquote(sqrt(.(lNota2))), side=1, line=1.5)
dev.off()

# New _v2
png("FatFigures_v2/Xi_sigE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","sigE"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(sigE), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[E])), side=1, line=1.5)
dev.off()


# ================================
# = Dive into cyclic time series =
# ================================
cl1 <- is.finite(final[,"Period"])&is.finite(final[,"sh_0"]) & final[,"Type"]=="Bio" & final[,"N"]>=30
cl2 <- is.element(finalFrame[,"location"],unique(final[cl1,"fitBy"])) & is.element(finalFrame[,"variable"],unique(final[cl1,"Variable"]))
cycl <- finalFrame[cl2,]

c0 <- ggplot(cycl, aes(x=year4, y=exp(Data), group=location)) + f2_theme
dev.new()
c0 + geom_line() + facet_grid(variable~location,scales="free")



# 
# quickAvg <- function(x, ...){two <- apply(x[,-c(1,2)], 2, as.numeric); colMeans(two, ...)}
# qaLogic <- !is.element(final[,"Type"],c("Met","Cosmic")) & 
# 
# avgFinal <- ddply(final[,c("fitBy","Type","P","Q","Lambda","Order","N","sh_0","Level2_time","InfE","Einf")], .variables=c("fitBy","Type"), .fun=quickAvg, na.rm=TRUE)
# avgFinal2 <- avgFinal[avgFinal[,"InfE"]<40,]
# 
# summary(lm())
# plot(avgFinal2[,"Order"], avgFinal2[,"N"])


# 
# durs <- exp(seq(1, 15, by=0.5))
# exFlood <- matrix(rep(NA, 4*length(durs)), ncol=4, dimnames=list(NULL, c("meme", "mema", "mama", "mame")))
# for(d in seq_along(durs)){
# 	td_me <- rep(NA,40)
# 	td_ma <- rep(NA, 40)
# 	for(i in 1:40){
# 		tt <- rt(n=durs[d], df=3)
# 		td_me[i] <- mean(tt)
# 		td_ma[i] <- max(tt)
# 	}
# 	exFlood[d,"meme"] <- mean(td_me)
# 	exFlood[d,"mema"] <- max(td_me)
# 	exFlood[d,"mama"] <- max(td_ma)
# 	exFlood[d,"mame"] <- mean(td_ma)
# }
# dev.new(); par(mfrow=c(2,2))
# for(i in 1:4){
# 	plot(log(durs), exFlood[,i], ylab=colnames(exFlood)[i])
# }
# 
# 
# durs <- exp(seq(1, 15, by=0.5))
# exRain <- matrix(rep(NA, 4*length(durs)), ncol=4, dimnames=list(NULL, c("meme", "mema", "mama", "mame")))
# for(d in seq_along(durs)){
# 	td_me <- rep(NA,40)
# 	td_ma <- rep(NA, 40)
# 	for(i in 1:40){
# 		tt <- rt(n=durs[d], df=30)
# 		td_me[i] <- mean(tt)
# 		td_ma[i] <- max(tt)
# 	}
# 	exRain[d,"meme"] <- mean(td_me)
# 	exRain[d,"mema"] <- max(td_me)
# 	exRain[d,"mama"] <- max(td_ma)
# 	exRain[d,"mame"] <- mean(td_ma)
# }
# dev.new(); par(mfrow=c(2,2))
# for(i in 1:4){
# 	plot(log(durs), exRain[,i], ylab=colnames(exRain)[i])
# }


# 
# nLakes0 <- readOGR(dsn="/Users/Battrd/Desktop/lakeShape/nhld_study_lakes/", layer="nhld_study_lakes")
# nLakes0@data$id <- nLakes0@data[,"LAKEID"]
# nLakes_points <- fortify(nLakes0, region="id")
# nLakes <- join(nLakes_points, nLakes0@data, by="id")
# 
# sLakes0 <- readOGR(dsn="/Users/Battrd/Desktop/lakeShape/yld_study_lakes/", layer="yld_study_lakes")
# sLakes0@data$id <- sLakes0@data[,"LAKEID"]
# sLakes_points <- fortify(sLakes0, region="id")
# sLakes <- join(sLakes_points, sLakes0@data, by="id")
# 
# Lakes <- rbind(nLakes, sLakes)
# 
# lakeTheme <- theme_classic() + theme(line=element_blank(), text=element_blank(), plot.background=element_rect(fill="transparent",colour=NA), panel.background=element_rect(fill="transparent",colour=NA))
# 
# ggplot(nLakes) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.0,-10)) + scale_y_continuous(expand=c(0.0,-10)) + lakeTheme
# 
# ggplot(nLakes) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.0,-10)) + scale_y_continuous(expand=c(0.0,-10)) + lakeTheme
# 
# 
# 	unl <- unique(Lakes[,"id"])
# 		for(i in 1:length(unl)){
# 			lpng <- ggplot(Lakes[Lakes[,"id"]==unl[i],]) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.0,-10)) + scale_y_continuous(expand=c(0.0,-10)) + lakeTheme
# 			png(file=paste(unl[i],"outline.png",sep=""), width=1, height=1, units="in", bg="transparent", res=300)
# 			print(lpng)
# 			dev.off()
# 		}
# 		
# 		
# 		
# 		
# unl <- unique(Lakes[,"id"])
# 	for(i in 1:length(unl)){
# 		if(is.element(unl[i], c("CR", "CB", "FI", "TB"))){
# 			lpng <- ggplot(Lakes[Lakes[,"id"]==unl[i],]) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.0,0)) + scale_y_continuous(expand=c(0.0,0)) + lakeTheme
# 			png(file=paste(unl[i],"_outline_2.png",sep=""), width=10, height=10, units="in", bg="transparent", res=72)
# 			print(lpng)
# 			dev.off()
# 		}
# 		
# 		if(is.element(unl[i], c("SP"))){
# 			lpng <- ggplot(Lakes[Lakes[,"id"]==unl[i],]) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.5,0)) + scale_y_continuous(expand=c(0.0,0)) + lakeTheme
# 			png(file=paste(unl[i],"_outline_2.png",sep=""), width=10, height=10, units="in", bg="transparent", res=72)
# 			print(lpng)
# 			dev.off()
# 		}
# 		
# 		if(is.element(unl[i], c("TR"))){
# 			lpng <- ggplot(Lakes[Lakes[,"id"]==unl[i],]) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.2,0)) + scale_y_continuous(expand=c(0.0,0)) + lakeTheme
# 			png(file=paste(unl[i],"_outline_2.png",sep=""), width=10, height=10, units="in", bg="transparent", res=72)
# 			print(lpng)
# 			dev.off()
# 		}
# 		
# 		if(!is.element(unl[i], c("CR", "CB", "FI", "TR", "SP", "TB"))){
# 			lpng <- ggplot(Lakes[Lakes[,"id"]==unl[i],]) + aes(long, lat, group=group) + geom_polygon(fill="#386cb0") + scale_x_continuous(expand=c(0.0,-10)) + scale_y_continuous(expand=c(0.0,-10)) + lakeTheme
# 			png(file=paste(unl[i],"_outline_2.png",sep=""), width=10, height=10, units="in", bg="transparent", res=72)
# 			print(lpng)
# 			dev.off()
# 		}
# 				
# 				
# }
# 
# 
# 
myT <- pt(1:100, df=5, lower.tail=FALSE)
nT <- 100
mBe <- -2
mEn <- 50
N=100
gT <- function(N, DF, be, en){(dcauchy(seq(be, en, length.out=N), scale=DF))}
myTs <- data.frame("Quant"=c(seq(mBe, mEn, length.out=N), seq(mBe+0, mEn, length.out=N)), "Prob"=c(gT(N,25, be=mBe, en=mEn),gT(N, 0.25, be=mBe+4.75, en=mEn)), "Type"=rep(c("fat","thin"), each=N))

# png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_tails.png", width=8, height=3, units="in", res=300, bg=myWhite)
# ggplot(myTs, aes(x=Quant, y=Prob, color=Type)) + geom_line(size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=Prob, group=Type, fill=Type), alpha=0.5, size=0, show_guide=FALSE) + theme(title=element_blank(), text=element_blank(), panel.background=element_rect(fill="transparent", colour=NA), plot.background=element_rect(fill="transparent", colour=NA), panel.grid=element_blank(), axis.ticks=element_blank())
# dev.off()



egQuants <- seq(-10, 10, length.out=500)
normEg <- dnorm(egQuants)
tEg <- dt(egQuants, df=1)

# png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_norm.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
# par(family="Times", mfrow=c(1,1), ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
# plot(egQuants, normEg, type="l", xlab="X", ylab="Density", bty="l", lwd=2)
# dev.off()

# png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_norm_t.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
# par(family="Times", mfrow=c(1,1), ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
# plot(egQuants, normEg, type="l", xlab="X", ylab="Density", bty="l", lwd=2)
# lines(egQuants, tEg, lwd=5)
# dev.off()


# png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_norm_t_Tail.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
# par(family="Times", mfrow=c(1,1), ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
# plot(egQuants[301:500], tEg[301:500], type="l", xlab="X", ylab="Density", bty="l", lwd=5, ylim=c(-0.0001, max(tEg[301:500])))
# abline(h=0, lty="dashed")
# lines(egQuants[301:500], normEg[301:500], lwd=2)
# dev.off()

# 

# 
# library("fExtremes")
# egQuants2 <- seq(0, 8, length.out=250)
# norm2Eg <- dnorm(egQuants2, mean=3, sd=1)
# gumbelEg <- dgev(egQuants2, xi=0, mu=3, beta=1)
# frechetEg <- dgev(egQuants2, xi=0.51, mu=3, beta=1)
# weibullEg <- dgev(egQuants2, xi=-0.51, mu=3, beta=1)
# 
# cols3 <- c("normal"="black","gumbel"="#e41a1c", "weibull"="#377eb8", "frechet"="#4daf4a")
# fatEg <- cbind(data.frame("quant"=egQuants2, "normal"=norm2Eg), data.frame("gumbel"=gumbelEg), data.frame("weibull"=weibullEg), data.frame("frechet"=frechetEg))
# fatGG <- ggplot(fatEg, aes(x=quant))
# 
# fatTheme <- theme(title=element_blank(), text=element_blank(), panel.background=element_rect(fill="transparent", colour=NA), plot.background=element_rect(fill="transparent", colour=NA), panel.grid=element_blank(), axis.ticks=element_blank())
# 
# # fatGG + geom_line(aes(), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=density, group=type, fill=type), alpha=0.1, size=0, show_guide=FALSE) + fatTheme +annotate("text", x=1, y=0.35, label="xi", parse=TRUE, size=20)
# 
# fatColor <- scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3)
# fatNote <- annotate("text", x=1, y=0.35, label="xi", parse=TRUE, size=20)
# 
# fatNorm <- fatGG + geom_line(aes(y=normal, colour="normal"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=normal, colour="normal", fill="normal"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme + annotate("text", x=1, y=0.35, label="normal", parse=TRUE, size=10)
# 
# fatGumbel <- fatGG + geom_line(aes(y=gumbel, colour="gumbel"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=gumbel, ymin=0, ymax=gumbel, colour="gumbel", fill="gumbel"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=5.7, y=0.15, label="xi==0", parse=TRUE, size=10)
# 
# fatWeibull <- fatGG + geom_line(aes(y=weibull, colour="weibull"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=weibull, ymin=0, ymax=weibull, colour="weibull", fill="weibull"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=5, y=0.38, label="xi<0", parse=TRUE, size=10)
# 
# fatFrechet <- fatGG + geom_line(aes(y=frechet, colour="frechet"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=frechet, ymin=0, ymax=frechet, colour="frechet", fill="frechet"), alpha=0.1, size=0, show_guide=FALSE) + scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme+ annotate("text", x=7.4, y=0.07, label="xi>0", parse=TRUE, size=10)
# 
# #Normal and gumbel
# # fatGG + 
# # 	geom_line(aes(y=normal, colour="normal"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=normal, colour="normal", fill="normal"), alpha=0.1, size=0, show_guide=FALSE) +
# # 	geom_line(aes(y=gumbel, colour="gumbel"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=gumbel, colour="gumbel", fill="gumbel"), alpha=0.1, size=0, show_guide=FALSE) +
# # 	scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatNote + fatTheme
# 	
# 	# dev.new(width=10, height=7)
# 	png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_fats_grid.png", width=10, height=7, units="in", res=150, bg=myWhite)
# 	print(grid.arrange(fatNorm, fatGumbel, fatWeibull, fatFrechet, ncol=2))
# 	dev.off()
# 	
# 	png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_fats.png", width=10, height=10, units="in", res=150, bg=myWhite)
# fatGG + 
# 	geom_line(aes(y=normal, colour="normal"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=normal, colour="normal", fill="normal"), alpha=0.1, size=0, show_guide=FALSE) +annotate("text", x=1, y=0.2, label="normal", parse=TRUE, size=10)+
# 	geom_line(aes(y=gumbel, colour="gumbel"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(ymin=0, ymax=gumbel, colour="gumbel", fill="gumbel"), alpha=0.1, size=0, show_guide=FALSE) +annotate("text", x=5.35, y=0.14, label="xi==0", parse=TRUE, size=10)+
# 	geom_line(aes(y=weibull, colour="weibull"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=weibull, ymin=0, ymax=weibull, colour="weibull", fill="weibull"), alpha=0.1, size=0, show_guide=FALSE) +annotate("text", x=4.45, y=0.38, label="xi<0", parse=TRUE, size=10)+
# 	geom_line(aes(y=frechet, colour="frechet"), size=1.5,show_guide=FALSE) + geom_ribbon(aes(y=frechet, ymin=0, ymax=frechet, colour="frechet", fill="frechet"), alpha=0.1, size=0, show_guide=FALSE) +annotate("text", x=7.5, y=0.05, label="xi>0", parse=TRUE, size=10)+
# 	scale_colour_manual(values=cols3) + scale_fill_manual(values=cols3) + fatTheme
# dev.off()
# #
# 	
# 
# 
# 
# 
# png("/Users/Battrd/Documents/School&Work/Presentations/CFL_Seminar/Batt_CFL_Seminar_FatTails_11Dec2013/eg_fat.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
# par(family="Times", mfrow=c(1,1), ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
# plot(egQuants2, frechetEg, type="l", xlab="X", ylab="Density", bty="l", lwd=5, ylim=c(-0.0001, max(frechetEg)))
# abline(h=0, lty="dashed")
# lines(egQuants2, frechetEg, lwd=2)
# lines(egQuants2, gumbelEg, lwd=1)
# dev.off()





# ==========================
# = Shape Boxplots by Type =
# ==========================
png("FatFigures_v2/Box_ReturnShape.png", width=7, height=3.4, units="in", res=300, bg=myWhite)
par(mfrow=c(1,2), ps=9, mar=c(2.5,3.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0),family="Times")
boxplot((sh_0)~Type, data=sAP, ylab="", outline=FALSE, col="white")
mtext(expression(Value~~of~~xi~~from~~GEV), side=2, line=2)

boxplot(log10(Level2_time)~Type, data=sAP, ylab="", outline=FALSE, col="white")
mtext(paste("Return Time (log10(yrs)) \n [record + ", ((setThresh*100)-100), "%]", sep=""), side=2, line=1.5)
dev.off()






# ==================================================================
# = Return times predicted from 1) normal, 2) log-normal, & 3) GEV =
# ==================================================================
png("FatFigures_v2/ScatterReturns.png", width=6, height=2, bg=myWhite, units="in", res=300)

par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=9, cex=1, tcl=-0.25, mgp=c(3,0.5,0), family="Times")

sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
lT <- (setThresh*100)-100

#GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
mtext("Predicted by\nGEV", side=3, font=2)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)

#log normal
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim)
mtext("Predicted by\nLog-Normal", side=3, font=2)

# normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim)
mtext("Predicted by\nNormal", side=3, font=2)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
dev.off()

# ========================
# = BLANK nrml+l-nrml Scatter Return =
# ========================
png("FatFigures_v2/ScatterReturns_2blank.png", width=6, height=2, bg=myWhite, units="in", res=300)

par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=9, cex=1, tcl=-0.25, mgp=c(3,0.5,0), family="Times")

sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
lT <- (setThresh*100)-100

#GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
mtext("Predicted by\nGEV", side=3, font=2)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)

#log normal
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext("Predicted by\nLog-Normal", side=3, font=2)

# normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext("Predicted by\nNormal", side=3, font=2)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
dev.off()





# ========================
# = BLANK Scatter Return =
# ========================
png("FatFigures_v2/ScatterReturns_allBlank.png", width=6, height=2, bg=myWhite, units="in", res=300)
# dev.new(width=7, height=4)
par(mfrow=c(1,3), mar=c(2,2,0.5,0.5), oma=c(1,2,1.5,0), ps=9, cex=1, tcl=-0.25, mgp=c(3,0.5,0), family="Times")

sl_Ylim <- log10(range(Inf2NA(sAP[,c("Level2_normTime", "Level2_logNormTime", "Level2_time")]), na.rm=TRUE))
lT <- (setThresh*100)-100

#GEV
plot(sAP[,"sh_0"], log10(sAP[,"Level2_time"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext(expression(Value~~of~~xi~~from~~GEV), outer=TRUE, line=-0.25, side=1)
mtext("Predicted by\nGEV", side=3, font=2)
mtext(bquote(Return~~Time~~.("[")*log[10]*(yrs)*.("]")), side=2, line=2.5)

#log normal
mtext(bquote(.("[")~.(lT)*"%"~~over~~record~.("]")), side=2, line=1.5)
plot(sAP[,"sh_0"], log10(sAP[,"Level2_logNormTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext("Predicted by\nLog-Normal", side=3, font=2)

# normal
plot(sAP[,"sh_0"], log10(sAP[,"Level2_normTime"]), xlab="", ylab="", ylim=sl_Ylim, pch=NA)
mtext("Predicted by\nNormal", side=3, font=2)
# mtext(paste("Return Time (log10(yrs)) \n [record+", ((setThresh*100)-100), "%]", sep=""), side=2, line=2.25)
dev.off()






# =====================================================
# = Creat Return-Level plots for common distributions =
# =====================================================
StdShape <- matrix(data=NA, ncol=4, nrow=7, dimnames=list(c("Uniform", "Normal", "T-20", "Exp", "Log-Normal", "Cauchy", "Weibull"), c("N", "mu_0","sig_0", "sh_0")))

RanUnif <- runif(n=1E5)
TailUnif <- RanUnif[which(RanUnif> quantile(RanUnif, 0.90))]
StdShape["Uniform",] <- c("N"=length(TailUnif), gev.fit(TailUnif)$mle)

RanNorm <- rnorm(n=1E5)
TailNorm <- RanNorm[which(RanNorm> quantile(RanNorm, 0.90))]
StdShape["Normal",] <- c("N"=length(TailNorm), gev.fit(TailNorm)$mle)

RanT <- rt(n=1E5, df=15)
TailT <- RanT[which(RanT> quantile(RanT, 0.90))]
StdShape["T-20",] <- c("N"=length(TailT), gev.fit(TailT)$mle)

RanExp <- rexp(n=1E5, rate=1)
TailExp <- RanExp[which(RanExp> quantile(RanExp, 0.90))]
StdShape["Exp",] <- c("N"=length(TailExp), gev.fit(TailExp)$mle)

RanLNorm <- rlnorm(n=1E5)
# TailLNorm <- RanLNorm[which(RanLNorm> quantile(RanLNorm, 0.55))]
StdShape["Log-Normal",] <- c("N"=length(RanLNorm), gev.fit(RanLNorm)$mle)

RanCauchy <- rcauchy(n=1E5)
TailCauchy <- RanCauchy[which(RanCauchy> quantile(RanCauchy, 0.90))]
StdShape["Cauchy",] <- c("N"=length(TailCauchy), gev.fit(TailCauchy)$mle)

RanWeibull <- rweibull(n=1E5, shape=0.17)
TailWeibull <- RanWeibull[which(RanWeibull> quantile(RanWeibull, 0.90))]
StdShape["Weibull",] <- c("N"=length(TailWeibull), gev.fit(TailWeibull)$mle)


RetTimes <- seq(5, 200, by=5)
StdReturns <- matrix(NA, ncol=nrow(StdShape)+1, nrow=length(RetTimes), dimnames=list(NULL, c("RetTime", rownames(StdShape))))
StdReturns[,"RetTime"] <- RetTimes
for(j in 1:nrow(StdShape)){
	for(i in 1:length(RetTimes)){
		RetTime <- RetTimes[i]
		nExts <- StdShape[j,"N"]
		TS_Duration <- StdShape[j,"N"]
		mu <- StdShape[j, "mu_0"]
		sc <- StdShape[j, "sig_0"]
		xi <- mu <- StdShape[j, "sh_0"]
		StdReturns[i,j+1] <- xYr_Lvl(RetTime, nExts, TS_Duration, mu, sc, xi)
	}
}





# ===========================================
# = Create Return-Level Plots for LTER Data =
# ===========================================
LterNames1 <- do.call(paste, Lake_Params3[,c("fitBy", "Variable")])
LterNames <- gsub(" ", "_", LterNames1)
RetTimes <- seq(5, 200, by=5)
LterReturns <- matrix(NA, ncol=nrow(Lake_Params3)+1, nrow=length(RetTimes), dimnames=list(NULL, c("RetTime", LterNames)))
LterReturns[,"RetTime"] <- RetTimes
for(j in 1:nrow(Lake_Params3)){
	for(i in 1:length(RetTimes)){
		RetTime <- RetTimes[i]
		nExts <- Lake_Params3[j,"N"]
		TS_Duration <- Lake_Params3[j,"Duration"]
		mu <- Lake_Params3[j, "mu_0"]
		sc <- Lake_Params3[j, "sig_0"]
		xi <- Lake_Params3[j, "sh_0"]
		LterReturns[i,j+1] <- xYr_Lvl(RetTime, nExts, TS_Duration, mu, sc, xi)
	}
}
MyGray <- rgb(t(col2rgb("gray")), alpha=100, maxColorValue=255)
MyBlue <- rgb(t(col2rgb("blue")), alpha=80, maxColorValue=255)
MyRed <- rgb(t(col2rgb("red")), alpha=80, maxColorValue=255)
Bios <- which(Lake_Params3[,"Type"]=="Bio")
Physs <- which(Lake_Params3[,"Type"]=="Phys")
Chems <- which(Lake_Params3[,"Type"]=="Chem")
LterColors <- c()
LterColors[Bios] <- MyRed
LterColors[Physs] <- MyBlue
LterColors[Chems] <- MyGray


# ===========================================
# = Create a plot combining common and LTER =
# ===========================================
# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatFigures/")
# dev.new(width=3.5, height=6)
# png("/Users/Battrd/Documents/School&Work/Presentations/comMeeting/Returns.png", width=3, height=5, units="in", res=300, pointsize=9)
png("FatFigures_v2/Returns.png", width=7, height=3.4, units="in", res=300, bg=myWhite)
par(mfrow=c(1,2), mar=c(1.5,2.5,0.5,0.5), oma=c(0,0,0,0), ps=9, cex=1, family="Times")
plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:6){
	par(new=TRUE)
	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
}
legend("topleft", c("Uniform", "Normal", "T (df=20)"), bty="n", inset=c(-0.05, -0.03), bg="transparent")
legend("bottomright", c("Exponential", "Log-Normal", "Cauchy", "Weibull"), text.font=c(2,rep(1,3)), bty="n", bg="transparent")
mtext("Return Time", side=1, line=0.5 , outer=F)

plot(LterReturns[,1], LterReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n", col=LterColors[1])
for(i in 1:(ncol(LterReturns)-2)){
	par(new=TRUE)
	plot(LterReturns[,1], LterReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=1, col=LterColors[i+2], bty="n")	
}
mtext("Return Time", side=1, line=0.25, outer=F)
mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)
legend("bottomright", c("Physical", "Chemical", "Biological"), text.col=c("blue", "gray", "red"), bg="transparent")
dev.off()

#Just the known distributions
png("FatFigures_v2/Returns_1.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
par(mfrow=c(1,1), mar=c(1.5,2.5,0.5,0.5), oma=c(0,0,0,0), ps=9, cex=1, family="Times")
plot(StdReturns[,1], StdReturns[,2], type="l", xlab="", ylab="", xaxt="n", yaxt="n")
for(i in 1:6){
	par(new=TRUE)
	plot(StdReturns[,1], StdReturns[,i+2], type="l", xlab="", ylab="" ,xaxt="n", yaxt="n", lwd=ifelse(i==3, 3, 1), bty="n")
}
legend("topleft", c("Uniform", "Normal", "T (df=20)"), bty="n", inset=c(-0.05, -0.03), bg="transparent")
legend("bottomright", c("Exponential", "Log-Normal", "Cauchy", "Weibull"), text.font=c(2,rep(1,3)), bty="n", bg="transparent")
mtext("Return Time", side=1, line=0.25 , outer=F)
mtext("Return Level (relative)", side=2, line=-1.25, outer=TRUE)

dev.off()














# =====================
# = Time series plots =
# =====================
#3 panel graph, each panel is a different category
#pick the fattest variable from each category, and plot its time series
# then plot the time series of that variable from the lake for which the variable is the thinnest
# then draw abline() at the level of the N-year return time, where N is the number of years that the time series was observed
chemLog_sig <- Lake_Params3[,"p_sh_0"]<0.1
chemLog0 <- Lake_Params3[,"Type"] == "Chem" & Lake_Params3[,"N"] > 20
chemLog <- chemLog0 & chemLog_sig
MaxChem_Ind <- which(chemLog)[which.max(Lake_Params3[chemLog, "sh_0"])]
MaxChem <- Lake_Params3[MaxChem_Ind,]
Min_MaxChem_Ind <- which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxChem[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxChem <- Lake_Params3[Min_MaxChem_Ind,]

MaxPhys_Ind <- which(Lake_Params3[,"Type"] == "Phys"& Lake_Params3[,"N"] > 20)[which.max(Lake_Params3[which(Lake_Params3[,"Type"] == "Phys"& Lake_Params3[,"N"] > 20), "sh_0"])]
MaxPhys <- Lake_Params3[MaxPhys_Ind,]
Min_MaxPhys_Ind <- which(Lake_Params3[,"Variable"]==MaxPhys[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxPhys[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxPhys <- Lake_Params3[Min_MaxPhys_Ind,]

MaxBio_Ind <- which(Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20)[which.max(Lake_Params3[which(Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20), "sh_0"])]
MaxBio <- Lake_Params3[MaxBio_Ind,]
Min_MaxBio_Ind <- which(Lake_Params3[,"Variable"]==MaxBio[,"Variable"]& Lake_Params3[,"N"] > 20)[which.min(Lake_Params3[which(Lake_Params3[,"Variable"]==MaxBio[,"Variable"]& Lake_Params3[,"N"] > 20),"sh_0"])]
Min_MaxBio <- Lake_Params3[Min_MaxBio_Ind,]

MaxBio_Ind2_Logic <- Lake_Params3[,"Type"] == "Bio"& Lake_Params3[,"N"] > 20 & Lake_Params3[,"Variable"] != "CPUE"
MaxBio_Ind2 <- which(MaxBio_Ind2_Logic)[which.max(Lake_Params3[MaxBio_Ind2_Logic, "sh_0"])]
MaxBio2 <- Lake_Params3[MaxBio_Ind2,]
Min_MaxBio_Ind2_Logic <- Lake_Params3[,"Variable"] == MaxBio2[,"Variable"]& Lake_Params3[,"N"] > 20
Min_MaxBio_Ind2 <- which(Min_MaxBio_Ind2_Logic)[which.min(Lake_Params3[which(Min_MaxBio_Ind2_Logic),"sh_0"])]
Min_MaxBio2 <- Lake_Params3[Min_MaxBio_Ind2,]

dateMax <- function(x, name="extcoef"){
	Xr <- x[which.max(x[,name]),"sampledate"]
}

png("FatFigures_v2/fatTimeSeries.png", width=5, height=7, bg=myWhite, units="in", res=300)
# dev.new(height=7, width=5)
# pdf(file="EgTimeSeries_RetLvls_09May2013.pdf", width=5, height=7, pointsize=10)
par(mfrow=c(4,2), mar=c(1.5,1.5,0.5, 0.5), oma=c(0.75, 1, 0, 0), ps=9, cex=1, family="Times", mgp=c(1.5, 0.5, 0), tcl=-0.4)
# ================================
# = Light Extinction time series =
# ================================
LiExt_MaxTS <- subset(Data_X$LiExt, lakeid==as.character(MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MaxTS <- LiExt_MaxTS[complete.cases(LiExt_MaxTS),]
LiExt_MaxTS_years <- subset(Data_X$LiExt, lakeid==as.character(MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef", "year4"))
LiExt_MaxTS_years <- LiExt_MaxTS_years[complete.cases(LiExt_MaxTS_years),]
LiExt_MaxDates <- ddply(LiExt_MaxTS_years, .variables=c("year4"), .fun=dateMax)[,2]

Low_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"]-MaxPhys[,"se_mu_0"], MaxPhys[,"sig_0"]- MaxPhys[,"se_sig_0"], MaxPhys[,"sh_0"] - MaxPhys[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"]+MaxPhys[,"se_mu_0"], MaxPhys[,"sig_0"]+ MaxPhys[,"se_sig_0"], MaxPhys[,"sh_0"] + MaxPhys[,"se_sh_0"])
# Mean_RetLvl <- xYr_Lvl(30, MaxPhys[,"N"], MaxPhys[,"Duration"], MaxPhys[,"mu_0"], MaxPhys[,"sig_0"], MaxPhys[,"sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, LiExt_MaxTS[,2]))
plot(LiExt_MaxTS, type="l", col="gray", ylim=Ylim)
points(LiExt_MaxTS[is.element(LiExt_MaxTS[,"sampledate"], LiExt_MaxDates),c("sampledate", "extcoef")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(MaxPhys[,"fitBy"]), "ExtCoef", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxPhys[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(Light~Extinction), side=2, line=1.5)

LiExt_MinMaxTS <- subset(Data_X$LiExt, lakeid==as.character(Min_MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef"))
LiExt_MinMaxTS <- LiExt_MinMaxTS[complete.cases(LiExt_MinMaxTS),]
LiExt_MinMaxTS_years <- subset(Data_X$LiExt, lakeid==as.character(Min_MaxPhys[,"fitBy"]), select=c("sampledate", "extcoef", "year4"))
LiExt_MinMaxTS_years <- LiExt_MinMaxTS_years[complete.cases(LiExt_MinMaxTS_years),]
LiExt_MinMaxDates <- ddply(LiExt_MinMaxTS_years, .variables=c("year4"), .fun=dateMax)[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"]-Min_MaxPhys[,"se_mu_0"], Min_MaxPhys[,"sig_0"]- Min_MaxPhys[,"se_sig_0"], Min_MaxPhys[,"sh_0"] - Min_MaxPhys[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"]+Min_MaxPhys[,"se_mu_0"], Min_MaxPhys[,"sig_0"]+ Min_MaxPhys[,"se_sig_0"], Min_MaxPhys[,"sh_0"] + Min_MaxPhys[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, LiExt_MinMaxTS[,2]))
plot(LiExt_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(LiExt_MinMaxTS[is.element(LiExt_MinMaxTS[,"sampledate"], LiExt_MinMaxDates),c("sampledate", "extcoef")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(Min_MaxPhys[,"fitBy"]), "ExtCoef", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
xYr_Lvl(30, Min_MaxPhys[,"N"], Min_MaxPhys[,"Duration"], Min_MaxPhys[,"mu_0"], Min_MaxPhys[,"sig_0"], Min_MaxPhys[,"sh_0"])
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxPhys[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))


# ===============================
# = Nitrate/Nitrite time series =
# ===============================
NO3NO2_MaxTS <- subset(Data_X$Chem, lakeid==as.character(MaxChem[,"fitBy"]), select=c("sampledate", "no3no2"))
NO3NO2_MaxTS <- NO3NO2_MaxTS[complete.cases(NO3NO2_MaxTS),]
NO3NO2_MaxTS_years <- subset(Data_X$Chem, lakeid==as.character(MaxChem[,"fitBy"]), select=c("sampledate", "no3no2", "year4"))
NO3NO2_MaxTS_years <- NO3NO2_MaxTS_years[complete.cases(NO3NO2_MaxTS_years),]
NO3NO2_MaxDates <- ddply(NO3NO2_MaxTS_years, .variables=c("year4"), .fun=dateMax, name="no3no2")[,2]

Low_RetLvl <- xYr_Lvl(30, MaxChem[,"N"], MaxChem[,"Duration"], MaxChem[,"mu_0"]-MaxChem[,"se_mu_0"], MaxChem[,"sig_0"]- MaxChem[,"se_sig_0"], MaxChem[,"sh_0"] - MaxChem[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxChem[,"N"], MaxChem[,"Duration"], MaxChem[,"mu_0"]+MaxChem[,"se_mu_0"], MaxChem[,"sig_0"]+ MaxChem[,"se_sig_0"], MaxChem[,"sh_0"] + MaxChem[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, NO3NO2_MaxTS[,2]))
plot(NO3NO2_MaxTS, type="l", col="gray", ylim=Ylim)
points(NO3NO2_MaxTS[is.element(NO3NO2_MaxTS[,"sampledate"], NO3NO2_MaxDates),c("sampledate", "no3no2")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(MaxChem[,"fitBy"]), "NO3NO2", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxChem[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(NO[3] ~ "&" ~NO[2] ~(mu*g~L^-1)), side=2, line=1.5)


NO3NO2_MinMaxTS <- subset(Data_X$Chem, lakeid==as.character(Min_MaxChem[,"fitBy"]), select=c("sampledate", "no3no2"))
NO3NO2_MinMaxTS <- NO3NO2_MinMaxTS[complete.cases(NO3NO2_MinMaxTS),]
NO3NO2_MinMaxTS_years <- subset(Data_X$Chem, lakeid==as.character(Min_MaxChem[,"fitBy"]), select=c("sampledate", "no3no2", "year4"))
NO3NO2_MinMaxTS_years <- NO3NO2_MinMaxTS_years[complete.cases(NO3NO2_MinMaxTS_years),]
NO3NO2_MinMaxDates <- ddply(NO3NO2_MinMaxTS_years, .variables=c("year4"), .fun=dateMax, name="no3no2")[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxChem[,"N"], Min_MaxChem[,"Duration"], Min_MaxChem[,"mu_0"]-Min_MaxChem[,"se_mu_0"], Min_MaxChem[,"sig_0"]- Min_MaxChem[,"se_sig_0"], Min_MaxChem[,"sh_0"] - Min_MaxChem[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxChem[,"N"], Min_MaxChem[,"Duration"], Min_MaxChem[,"mu_0"]+Min_MaxChem[,"se_mu_0"], Min_MaxChem[,"sig_0"]+ Min_MaxChem[,"se_sig_0"], Min_MaxChem[,"sh_0"] + Min_MaxChem[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, NO3NO2_MinMaxTS[,2]))
plot(NO3NO2_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(NO3NO2_MinMaxTS[is.element(NO3NO2_MinMaxTS[,"sampledate"], NO3NO2_MinMaxDates),c("sampledate", "no3no2")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
abline(h=LterReturns[6, paste(as.character(Min_MaxChem[,"fitBy"]), "NO3NO2", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxChem[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))


# ============================
# = Chlorophyll Density time series =
# ============================
chlor_MaxTS <- subset(Data_X$Chl, lakeid==as.character(MaxBio2[,"fitBy"]), select=c("sampledate", "chlor"))
chlor_MaxTS <- chlor_MaxTS[complete.cases(chlor_MaxTS),]
chlor_MaxTS_years <- subset(Data_X$Chl, lakeid==as.character(MaxBio2[,"fitBy"]), select=c("sampledate", "chlor", "year4"))
chlor_MaxTS_years <- chlor_MaxTS_years[complete.cases(chlor_MaxTS_years),]
chlor_MaxDates <- ddply(chlor_MaxTS_years, .variables=c("year4"), .fun=dateMax, name="chlor")[,2]

Low_RetLvl <- xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]-MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]- MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] - MaxBio2[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+ MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + MaxBio2[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, chlor_MaxTS[,2]))
plot(chlor_MaxTS, type="l", col="gray", ylim=Ylim)
points(chlor_MaxTS[is.element(chlor_MaxTS[,"sampledate"], chlor_MaxDates),c("sampledate", "chlor")],  col="red", pch=20) #it took me forever to do this... like at least 60 minutes
xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+1.96*MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+ 1.96*MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + 1.96*MaxBio2[,"se_sh_0"])
xYr_Lvl(30, MaxBio2[,"N"], MaxBio2[,"Duration"], MaxBio2[,"mu_0"]+MaxBio2[,"se_mu_0"], MaxBio2[,"sig_0"]+MaxBio2[,"se_sig_0"], MaxBio2[,"sh_0"] + MaxBio2[,"se_sh_0"])
abline(h=LterReturns[6, paste(as.character(MaxBio2[,"fitBy"]), "Chla", sep="_")], lty="dotted", lwd=2)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
legend("topleft", legend=bquote(xi ==  .(round(MaxBio2[,"sh_0"],2))), bty="n", inset=c(-.1,0.01))
mtext(expression(Chlorophyll~(mu*g~L^-1)), side=2, line=1.5)

chlor_MinMaxTS <- subset(Data_X$Chl, lakeid==as.character(Min_MaxBio2[,"fitBy"]), select=c("sampledate", "chlor"))
chlor_MinMaxTS <- chlor_MinMaxTS[complete.cases(chlor_MinMaxTS),]
chlor_MinMaxTS_years <- subset(Data_X$Chl, lakeid==as.character(Min_MaxBio2[,"fitBy"]), select=c("sampledate", "chlor", "year4"))
chlor_MinMaxTS_years <- chlor_MinMaxTS_years[complete.cases(chlor_MinMaxTS_years),]
chlor_MinMaxDates <- ddply(chlor_MinMaxTS_years, .variables=c("year4"), .fun=dateMax, name="chlor")[,2]

Low_RetLvl <- xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]-Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]- Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] - Min_MaxBio2[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]+Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]+ Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] + Min_MaxBio2[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, chlor_MinMaxTS[,2]))
plot(chlor_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(chlor_MinMaxTS[is.element(chlor_MinMaxTS[,"sampledate"], chlor_MinMaxDates),c("sampledate", "chlor")],  col="blue", pch=20) #it took me forever to do this... like at least 60 minutes
# xYr_Lvl(30, Min_MaxBio2[,"N"], Min_MaxBio2[,"Duration"], Min_MaxBio2[,"mu_0"]+1.96*Min_MaxBio2[,"se_mu_0"], Min_MaxBio2[,"sig_0"]+ 1.96*Min_MaxBio2[,"se_sig_0"], Min_MaxBio2[,"sh_0"] + 1.96*Min_MaxBio2[,"se_sh_0"])
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio2[,"fitBy"]), "Chla", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi == .(round(Min_MaxBio2[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))



# =========================
# = Fish CPUE Time series =
# =========================
CPUE_MaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
CPUE_MaxTS <- CPUE_MaxTS[complete.cases(CPUE_MaxTS),]
Low_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]-MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]- MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] - MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, MaxBio[,"N"], MaxBio[,"Duration"], MaxBio[,"mu_0"]+MaxBio[,"se_mu_0"], MaxBio[,"sig_0"]+ MaxBio[,"se_sig_0"], MaxBio[,"sh_0"] + MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MaxTS[,2]))
plot(CPUE_MaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MaxTS, col="red", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.03))
mtext(expression(Fish~CPUE~(individuals)), side=2, line=1.5)

CPUE_MinMaxTS <- subset(Data_X$Fish_ByGear, lakeid==as.character(Min_MaxBio[,"fitBy"]) & gearid=="ELFISH", select=c("year4", "cpue1_Sum"))
Low_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]-Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]- Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] - Min_MaxBio[,"se_sh_0"])
High_RetLvl <- xYr_Lvl(30, Min_MaxBio[,"N"], Min_MaxBio[,"Duration"], Min_MaxBio[,"mu_0"]+Min_MaxBio[,"se_mu_0"], Min_MaxBio[,"sig_0"]+ Min_MaxBio[,"se_sig_0"], Min_MaxBio[,"sh_0"] + Min_MaxBio[,"se_sh_0"])
Ylim <- range(c(Low_RetLvl, High_RetLvl, CPUE_MinMaxTS[,2]))
plot(CPUE_MinMaxTS, type="l", col="gray", ylim=Ylim)
points(CPUE_MinMaxTS, col="blue", pch=20)
abline(h=c(Low_RetLvl, High_RetLvl), lty="dotted", lwd=0.75)
abline(h=LterReturns[6, paste(as.character(Min_MaxBio[,"fitBy"]), "CPUE", sep="_")], lty="dotted", lwd=2)
legend("topleft", legend=bquote(xi ==  .(round(Min_MaxBio[,"sh_0"],2))), bty="n", inset=c(-.1,-0.07))
mtext(expression(Date), side=1, line=-0.25, outer=TRUE)
dev.off()







