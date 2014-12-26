


load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.shortMet.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.RData")


library(lme4)
library(multcomp)

# ========================
# = Definte some objects =
# ========================
data.fat.noMet <- data.fat[data.fat[,"Type"]!="Meteorological",] # no meteorological time series

# Weights
xiWeights.noMet <- 1/(data.fat.noMet[,"se.sh_0"]^2) # weights for the data set with no met
xiWeights.shortMet <- 1/(data.fat.shortMet[,"se.sh_0"]^2)

# contrasts for comparison
contr <- c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0")
contr.noMet <- c("Chemical-Biological=0","Physical-Biological=0")

# ============================================================
# = Drop Meteorological, Compare Type|Location, with weights =
# ============================================================
# These results give the result that the non-met variables are less fat-tailed than biology
# Compare w/o weights
fm4 <- (lmer(sh_0~Type+N+(Type|location), data=data.fat.noMet))
summary(glht(fm4, linfct=mcp(Type=contr.noMet)))

# Do it with weights
fm4.w <- (lmer(sh_0~Type+N+(Type|location), weights=xiWeights.noMet, data=data.fat.noMet))
summary(glht(fm4.w, linfct=mcp(Type=contr.noMet)))



# ===================================
# = With Met, after Consulting Tony =
# ===================================
# tried again after consulting with Tony
re1 <- lmer(sh_0 ~ Type + N + (1|location), data=data.fat.shortMet) # random effect model 1
summary(glht(re1, linfct=mcp(Type=contr)))

re1.w <- lmer(sh_0 ~ Type + N + (1|location), data=data.fat.shortMet, weights=xiWeights.shortMet) # random effect w/ weights
summary(glht(re1.w, linfct=mcp(Type=contr)))

re2 <- lmer(sh_0 ~ Type + N + (N|location), data=data.fat.shortMet) # random effect model 2
summary(glht(re2, linfct=mcp(Type=contr)))

re2.w <- lmer(sh_0 ~ Type + N + (N|location), data=data.fat.shortMet, weights=xiWeights.shortMet) # re2 with weights
summary(glht(re2.w, linfct=mcp(Type=contr)))






