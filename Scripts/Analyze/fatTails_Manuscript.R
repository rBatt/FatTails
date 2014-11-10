
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Simulate/fatARMA_Sim.R")

library(lme4)
library(multcomp)

# =========================
# = Data start/ end dates =
# =========================
dataDates <- ddply(data.stat, c("location","variable","taxID","Type"), function(x){data.frame("year.start"=min(x[,"year4"], na.rm=TRUE), "year.end"=max(x[,"year4"], na.rm=TRUE))})


# ======================================
# = Relationship between Xi and its SE =
# ======================================
summary(lm(sh_0~se.sh_0, data=data.fat))

# ================================================
# = Relationship between Xi se and variable Type =
# ================================================
summary(aov(se.sh_0~Type, data=data.fat))
1-(3.1860/(3.1860+ 0.1482)) # R-squared = 1-(RSS/TSS)

# ======================
# = Regressions for Xi =
# ======================

	# =====================
	# = Linear regression =
	# =====================
mod1 <- lm(sh_0~Type, data=data.fat)
summary(mod1)
summary(glht(mod1, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0"))))
TukeyHSD(aov(mod1))
	# ================================
	# = Linear regression w/ weights =
	# ================================
xiWeights0 <- 1/(data.fat[,"se.sh_0"]^2)
xiWeights <- xiWeights0/sum(xiWeights0)
mod2 <- lm(sh_0~Type, data=data.fat, weights=xiWeights)
summary(mod2)
TukeyHSD(aov(mod2))
# summary(mod2) # all significantly less than intercept, except biological, which serves as the intercept and is greater than 0
# TukeyHSD(aov(mod2))

 # "Biological variables were more fat-tailed than the other types of variables, which had significantly thinner tails"
summary(glht(mod2, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0"))))

summary(glht(mod2, linfct=mcp(Type="Dunnett")))

 # "Furthermore, tail thickness (xi) tended to decrease across the gradient of biological to meteorological variables"
summary(glht(mod2, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0", "Chemical-Physical=0", "Chemical-Meteorological=0","Physical-Meteorological=0")))) # WRONG - AFTER WEIGHT, THE REVERSE IS ACTUALLY TRUE! all significant except difference between physical and meteorological

	# ===================================
	# = Tail thickness across locations =
	# ===================================
notMet <- !data.fat[,"Type"]=="Meteorological"
mod3 <- lm(sh_0~Type+location, data=data.fat[notMet,], weights=xiWeights[notMet])
summary(mod3)
TukeyHSD(aov(mod3))
summary(aov(mod3))
summary(aov(lm(sh_0~Type+location, data=data.fat, weights=xiWeights)))


fGen <- fish.gev[fish.gev[,"taxLvl"]=="Genus",]
summary(lmer(sh_0~ Family + (1|Orde/Genus), data=fGen))


# =========================
# = Screenshot regression =
# =========================
summary(lm(sh_0~Type*residual_sh_0, data=data.2))


# ======================================================================
# = Proportion of time series from each variable type in each category =
# ======================================================================
xiCat <- rep(NA, nrow(data.fat))
xiCat[data.fat[,"shape.sig"]>0] <- "Fat"
xiCat[data.fat[,"shape.sig"]<0] <- "Bounded"
xiCat[data.fat[,"shape.sig"]==0] <- "Thin"
xiCat <- factor(xiCat)

table(data.fat[,"Type"], xiCat)
table(data.fat[,"Type"], xiCat)/rowSums(table(data.fat[,"Type"], xiCat))

# # ==============================================================
# # = Proportion fish time series in each category of tailedness =
# # ==============================================================
# xiCat.fish <- rep(NA, nrow(data.fat[data.fat[,"variable"]=="cpue1_Sum",]))
# xiCat.fish[data.fat[data.fat[,"variable"]=="cpue1_Sum","shape.sig"]>0] <- "Fat"
# xiCat.fish[data.fat[data.fat[,"variable"]=="cpue1_Sum","shape.sig"]<0] <- "Bounded"
# xiCat.fish[data.fat[data.fat[,"variable"]=="cpue1_Sum","shape.sig"]==0] <- "Thin"
# xiCat.fish <- factor(xiCat.fish)
# 
# xiCat.fish.table <- table(factor(as.character(data.fat[data.fat[,"variable"]=="cpue1_Sum","taxID"])), xiCat.fish)
# xiCat.fish.table <- xiCat.fish.table[order(xiCat.fish.table[,"Fat"]),]
# order(xiCat.fish.table/rowSums(xiCat.fish.table)[,"Fat"])
# (xiCat.fish.table/rowSums(xiCat.fish.table))[order((xiCat.fish.table/rowSums(xiCat.fish.table))[,"Fat"]),]



# ===================
# = Rank of fish Xi =
# ===================
xi.fish <- data.fat[data.fat[,"variable"]=="cpue1_Sum",c("taxID","location","sh_0","se.sh_0","shape.sig")]
xi.fish[,"lowCI"] <- xi.fish[,"sh_0"] - qnorm(0.05, lower.tail=FALSE)*xi.fish[,"se.sh_0"]
xi.fish[,"highCI"] <- xi.fish[,"sh_0"] + qnorm(0.05, lower.tail=FALSE)*xi.fish[,"se.sh_0"]

head(xi.fish[order(xi.fish[,"sh_0"], decreasing=TRUE),], 10)
head(xi.fish[order(xi.fish[,"shape.sig"], decreasing=TRUE),], 10)


# ===================
# = Rank of zoop Xi =
# ===================
xi.zoop <- data.fat[data.fat[,"variable"]=="density",c("taxID","location","sh_0","se.sh_0","shape.sig")]
xi.zoop[,"lowCI"] <- xi.zoop[,"sh_0"] - qnorm(0.05, lower.tail=FALSE)*xi.zoop[,"se.sh_0"]
xi.zoop[,"highCI"] <- xi.zoop[,"sh_0"] + qnorm(0.05, lower.tail=FALSE)*xi.zoop[,"se.sh_0"]

head(xi.zoop[order(xi.zoop[,"sh_0"], decreasing=TRUE),], 10)
head(xi.zoop[order(xi.zoop[,"shape.sig"], decreasing=TRUE),], 10)


# ===================
# = Rank of chlor Xi =
# ===================
xi.chlor <- data.fat[data.fat[,"variable"]=="chlor",c("taxID","location","sh_0","se.sh_0","shape.sig")]
xi.chlor[,"lowCI"] <- xi.chlor[,"sh_0"] - qnorm(0.05, lower.tail=FALSE)*xi.chlor[,"se.sh_0"]
xi.chlor[,"highCI"] <- xi.chlor[,"sh_0"] + qnorm(0.05, lower.tail=FALSE)*xi.chlor[,"se.sh_0"]

head(xi.chlor[order(xi.chlor[,"sh_0"], decreasing=TRUE),], 10)
head(xi.chlor[order(xi.chlor[,"shape.sig"], decreasing=TRUE),], 10)



# ============================
# = Compare distribution AIC =
# ============================

fat.v.norm <- function(x){
	vari <- as.character(x[,"variable"])
	loca <- as.character(x[,"location"])
	taxid <- as.character(x[,"taxID"])
	
	oQ <- data.max[data.max[,"variable"]==vari&data.max[,"location"]==loca&data.max[,"taxID"]==taxid,"Data"]
	
	
	exi <- x[,"sh_0"]
	escale <- x[,"sig_0"]
	eloc <- x[,"mu_0"]
	
	emean <- x[,"mean"]
	esd <- x[,"sd"]
	
	elmean <- x[,"logMean"]
	elsd <- x[,"logSd"]
	
	aicNorm <- 2*2 - 2*sum(dnorm(oQ, mean=emean, sd=esd, log=TRUE),na.rm=TRUE)
	aicLnorm <- 2*2 - 2*sum(dlnorm(oQ, meanlog=elmean, sdlog=elsd, log=TRUE),na.rm=TRUE)
	aicGEV <- 2*3 - 2*sum(log(dgev(oQ, xi=exi, mu=eloc, sigma=escale)), na.rm=TRUE)
	
	aics <- c("aic.norm"=aicNorm, "aic.lnorm"=aicLnorm, "aic.gev"=aicGEV)
	bestDist <- gsub("aic\\.", "", names(aics)[which.min(aics)])
	return(data.frame(x, t(aics), "bestDist"=bestDist))
}
bestDist <- ddply(data.fat, c("Type","taxLvl","taxID","location","variable"), fat.v.norm)
table(bestDist[,"bestDist"])

# "...selecting distributions by AIC reveals that most time series are best ch"
table(bestDist[,"bestDist"])/sum(table(bestDist[,"bestDist"]))

table(bestDist[,"bestDist"], bestDist[,"Type"])
bestDist[bestDist[,"Type"]=="Biological"&bestDist[,"bestDist"]=="gev",]


gevWin <- bestDist[bestDist[,"bestDist"]=="gev",]
gevWin <- gevWin[complete.cases(gevWin),]

summary((gevWin[, "Level2_logNormTime"] - gevWin[, "Level2_time"]))

summary(abs(gevWin[gevWin[,"sh_0"]>0, "Level2_logNormTime"] - gevWin[gevWin[,"sh_0"]>0, "Level2_time"]))


colorDens(vals=list(gevWin[, "Level2_logNormTime"], gevWin[, "Level2_time"]))


# ===========================================
# = Waiting Time, log-norm and GEV, example =
# ===========================================
fatLogic <- data.fat[,"shape.sig"]>0
bioLogic <- data.fat[,"Type"]=="Biological"

long.gev.bio <- which.max(data.fat[fatLogic&bioLogic,"Level2_time"])
long.gev <- which.max(data.fat[fatLogic,"Level2_time"])

long.lnorm.bio <- which.max(data.fat[fatLogic&bioLogic,"Level2_logNormTime"])
long.lnorm <- which.max(data.fat[fatLogic,"Level2_logNormTime"])


data.fat[fatLogic&bioLogic,][long.gev.bio,]

data.fat[fatLogic&bioLogic,][long.lnorm.bio,]















