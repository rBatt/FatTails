
# ================
# = Script Notes =
# ================
# "location" column is lake or region
# "sh_0" is xi (shape parameter from GEV)
# "se.sh_0" is the estimated standard error of xi
# "N" is the number of annual maxima
# "Duration" is the duration of the time series – not the same as N if there was a missing year (1986, 1988 would be duration = 3, N=2)
# "Type" is the variable category – Biological, Chemical, Physical, Meteorological


# ===================
# = Load Data Files =
# ===================
# NOTE: Tony, you can just double click the .RData files instead of running the ntext two lines
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.shortMet.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.RData")


# ==================
# = Load Libraries =
# ==================
library(lme4)
library(multcomp)
library(car)

# Weights
xiWeights <- 1/(data.fat[,"se.sh_0"]^2)
xiWeights.shortMet <- 1/(data.fat.shortMet[,"se.sh_0"]^2) 

hist(xiWeights^.5)
summary(data.fat)
cbind(data.fat[,c(1,4,5,9,10)],xiWeights)

# xiWeights.trun takes the one weight > 2000 and truncates it as 2000
xiWeights.trun <- xiWeights
xiWeights.trun[xiWeights.trun>2*10^3] <- 2*10^3
hist(xiWeights.trun^.5, breaks=40)

# Best analysis
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.fat))

# simple analysis
summary(lm(sh_0 ~ Type, data=data.fat))

# There is no N : Type interaction
Anova(lmer(sh_0 ~ Type * N + (1 | location), data=data.fat))

# Analyses with truncated weights
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.fat, weights=xiWeights.trun))
summary(lmer(sh_0 ~ Type * N + (1 | location), data=data.fat, weights=xiWeights.trun))
Anova(lmer(sh_0 ~ Type * N + (1 | location), data=data.fat, weights=xiWeights.trun))
summary(lmer(sh_0 ~ Type + N + (1 | location) + (0 + N | Type), data=data.fat, weights=xiWeights.trun))

# pairwise comparisons
data.BC <- data.fat[is.element(data.fat$Type,c("Biological","Chemical")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BC))
summary(lm(sh_0 ~ Type, data=data.BC))

data.BP <- data.fat[is.element(data.fat$Type,c("Biological","Physical")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BP))
summary(lm(sh_0 ~ Type, data=data.BP))

data.BM <- data.fat[is.element(data.fat$Type,c("Biological","Meteorological")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BM))
summary(lm(sh_0 ~ Type, data=data.BM))
