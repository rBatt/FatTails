
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatARMA_Sim.R")

library(lme4)
library(multcomp)

# ======================
# = Regressions for Xi =
# ======================
	# =====================
	# = Linear regression =
	# =====================
mod1 <- lm(sh_0~Type, data=data.fat)
summary(mod1)

	# ================================
	# = Linear regression w/ weights =
	# ================================
xiWeights0 <- 1/(data.fat[,"se.sh_0"]^2)
xiWeights <- xiWeights0/sum(xiWeights0)
mod2 <- lm(sh_0~Type, data=data.fat, weights=xiWeights)
# summary(mod2) # all significantly less than intercept, except biological, which serves as the intercept and is greater than 0
# TukeyHSD(aov(mod2))

 # "Biological variables were more fat-tailed than the other types of variables, which had significantly thinner tails"
summary(glht(mod2, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0"))))

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



