#RDB
#_v0 (01-Dec-2013): Summarize the ARMA analysis of fatFrame by 1) selecting the best model by AICc; 2) getting the AICc-weighted average of the eigenvalue. Also, this might be a good place to begin comparing to the Fat Tails analysis.

rm(list=ls())
graphics.off()

library("plyr")
library("rpart")
library("party")
library("RColorBrewer")

setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")

load("fatARMA_v1.RData") #this is the data file containing the completed ARMA analysis
load("All_Params_TurnExtreme_Fat_Data_v8.RData")
load("finalFrame_v1.RData")

source("FatTails_Functions_v7.R") #the logStat function in tonyARMA_short needs the Inf2NA function
source("TonySuggestions/tonyARMA_short_v3.R") #also loads GenSA and DEoptim packages
source("/Users/Battrd/Documents/School&Work/WiscResearch/dscat_v0.R")


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
final0[,"Variable"] <- as.factor(final0[,"Variable"])
final0[,"fitBy"] <- as.factor(final0[,"fitBy"])
final0[,"Type"] <- as.factor(final0[,"Type"])

final <- final0
sum(!is.finite(fatARMA2[,"Period"]))


# blah <- rpart(sh_0~Type+Variable+fitBy, data=final)
# plot(blah)
# text(blah, use.n=FALSE, all=FALSE, cex=0.8)


# blah2 <- ctree(Lambda~Order, data=final[!is.na(final[,"Lambda"]),], control=ctree_control(mincriterion=0.85, testtype="Univariate", mtry=0))
# plot(blah2, pt.cex=0.8)
# 
# 
# blah2 <- ctree(wOrder~Type+Lambda, data=final[!is.na(final[,"Order"]),], control=ctree_control(mincriterion=0.85, testtype="Univariate", mtry=0))
# plot(blah2, pt.cex=0.8)
# 
# 
# 
# plot(final[,"wOrder"], log10(final[,"Level1_time"]))
# 
# plot(final[,"Lambda"], log10(final[,"Level1_logNormTime"]))
# 
# plot(final[,"N"], final[,"Lambda"], xlim=c(10, 40))



# =================================
# = Lambda and wLambda Histograms =
# =================================
dev.new(width=3.4, height=3.4)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
hist(final[,"Lambda"], xlab=bquote(abs(abs(~lambda~phantom()))), ylab="Frequency", main="", xlim=c(0,1))
hist(final[final[,"Type"]=="Chem","Lambda"], add=TRUE, col="gray")
hist(final[final[,"Type"]=="Bio","Lambda"], add=TRUE, col="black")
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)


dev.new(width=3.4, height=3.4)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
hist(final[,"wLambda"], xlab=bquote(bar(abs(abs(~lambda~phantom())))), ylab="Frequency", main="", xlim=c(0,1))
hist(final[final[,"Type"]=="Chem","Lambda"], add=TRUE, col="gray")
hist(final[final[,"Type"]=="Bio","Lambda"], add=TRUE, col="black")
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)


# =====================
# = Period Histrogram =
# =====================
dev.new(width=3.4, height=3.4)
par(family="Times", ps=9, mar=c(2.5,2.5,0.5,0.5), tcl=c(-0.4), mgp=c(1.5, 0.5, 0))
hist(final[final[,"Period"]<150,"Period"], xlab="Period", ylab="Frequency", main="")
hist(final[final[,"Type"]=="Chem" & final[,"Period"]<150 ,"Period"], add=TRUE, col="gray")
hist(final[final[,"Type"]=="Bio" & final[,"Period"]<150,"Period"], add=TRUE, col="black")
legend("topright", legend=c("All Types", "Chem", "Bio"), pt.bg=c("white","gray", "black"), pch=22)


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
dev.new(width=3.4, height=3.4)
par(mar=c(2.5,2.5,0.5,0.5), cex=1, ps=9, mgp=c(1.5,0.5,0), family="Times")
barplot(tally, beside=TRUE, col=brewer.pal(n=4, name="Set1"), axisnames=TRUE, legend.text=c("Bio", "Chem","Phys","Met"), args.legend=list(bty="n"), xlab="Order", ylab="Poportion")




# =======================
# = Variance Parameters =
# =======================
InfE <- final[,"sigInf"]/final[,"sigE"]
Einf <- (final[,"sigE"])^2/(final[,"sigInf"])^2

wInfE <- final[,"wSigInf"]/final[,"wSigE"]
wEinf <- final[,"wSigE"]/final[,"wSigInf"]


dev.new(width=3.4, height=3.4)
dscat(final[,"Lambda"][InfE < 100], InfE[InfE < 100], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext("Lambda", side=1, line=1.5)
mtext(bquote(sigma[infinity]/sigma[E]), side=2, line=1.25)

dev.new(width=3.4, height=3.4)
dscat(final[,"sh_0"][InfE < 100], InfE[InfE < 100], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext(bquote(hat(xi)), side=1, line=1.5)
mtext(bquote(sigma[infinity]/sigma[E]), side=2, line=1.25)

dev.new(width=3.4, height=3.4)
dscat(final[,"sh_0"][InfE < 100], log10(InfE[InfE < 100]), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext(bquote(hat(xi)), side=1, line=1.5)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=2, line=1.25)
x1 <- final[,"sh_0"][InfE < 100]
y1 <- log10(InfE[InfE < 100])
crossing <- lm(y1~x1)
# summary(crossing)
abline(crossing, lty="dashed")


dev.new(width=3.4, height=3.4)
dscat(final[,"sh_0"][log10(Einf) >- 4], log10(Einf)[log10(Einf) > -4], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext(bquote(hat(xi)), side=1, line=1.5)
mtext("Envi explained variance", side=2, line=1.5)

dev.new(width=3.4, height=3.4)
dscat(log10(final[,"Level2_time"][InfE < 100]), (InfE[InfE < 100]), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext(bquote(Return~time~(log[10](years)*';'~10*'%'~over~current~record)), side=1, line=1.5)
mtext(bquote(sigma[infinity]/sigma[E]), side=2, line=1.25)


dev.new(width=3.4, height=3.4)
dscat((Einf), (log10(final[,"Level2_time"])), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext("Envi explained variance", side=1, line=1.5)
mtext(bquote(Return~time~(log[10](years)*';'~10*'%'~over~current~record)), side=2, line=1.5)




dev.new(width=3.4, height=3.4)
dscat(final[,"Lambda"], (Einf), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext("Envi explained variance", side=2, line=1.5)
mtext(bquote(abs(~abs(~lambda~"")~"")), side=1, line=1.5)


dev.new(width=3.4, height=3.4)
dscat(final[,"Lambda"], final[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext(bquote(xi), side=2, line=1.5)
mtext(bquote(abs(~abs(~lambda~"")~"")), side=1, line=1.5)

summary(lm(final[,"sh_0"]~log10(Einf)+final[,"Lambda"]))


#a lot here is significant ...
summary(lm(log10(Einf)~final[,"sh_0"]+final[,"Lambda"]))

summary(lm(log10(Einf)~final[,"sh_0"]+final[,"Type"]+final[,"Lambda"]))

summary(lm(log10(Einf)~final[,"sh_0"]+final[,"Type"]+final[,"Lambda"]+log10(final[,"Level1_time"])))


library("lme4")
library("MuMIn")
meFrame0 <- data.frame("Type"=final[,"Type"], "Einf"=Einf, "sh_0"=final[,"sh_0"], "Lambda"=final[,"Lambda"], "l10L1time"=log10(final[,"Level1_time"]), "fitBy"=final[,"fitBy"], "Variable"=final[,"Variable"], "InfE"=InfE)
meFrame <- meFrame0[complete.cases(meFrame0),]

# meMod10 <- lmer(log10(Einf)~sh_0+Lambda+(1|fitBy), data=meFrame)
# summary(meMod10)
# r.squaredGLMM(meMod10)

meMod1 <- lmer(log10(Einf)~sh_0+(1|fitBy)+(sh_0|Variable)+Lambda, data=meFrame)
summary(meMod1)
r.squaredGLMM(meMod1)

meMod1 <- lmer(log10(Einf)~sh_0+(1|Type)+Lambda, data=meFrame)
summary(meMod1)
r.squaredGLMM(meMod1)

meMod2 <- lmer(log10(InfE)~sh_0+(1|fitBy)+(sh_0|Variable)+Lambda, data=meFrame)
summary(meMod2)
r.squaredGLMM(meMod2)

#replicate estimates from the lmer model
newV <- 0.3672*final[,"sh_0"] + -0.7369*final[,"Lambda"]
dev.new(width=3.4, height=3.4)
dscat((newV), log10(Einf), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
# dscat((newV)[log10(Einf)>-5], log10(Einf)[log10(Einf)>-5], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext("Envi explained variance", side=2, line=1.5)
mtext("newV", side=1, line=1.5)

#replicate estimates from the lm model that only had Lambda and sh_0
newV <- 1.29104*final[,"sh_0"] + -1.30992*final[,"Lambda"]
dev.new(width=3.4, height=3.4)
# dscat((newV), log10(Einf), mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
dscat((newV)[log10(Einf)>-5], log10(Einf)[log10(Einf)>-5], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
mtext("Envi explained variance", side=2, line=1.5)
# mtext("newV", side=1, line=1.5)
mtext(bquote((1.3%*%xi)~~~-~~~(1.3%*%lambda)), side=1, line=1.5)


