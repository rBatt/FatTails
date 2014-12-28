
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
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.RData")


# ==================
# = Load Libraries =
# ==================
library(lme4)
library(car)


# Best analysis (Table S2)
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.fat))

# simple analysis (Table S1)
summary(lm(sh_0 ~ Type, data=data.fat))

# There is no N : Type interaction
Anova(lmer(sh_0 ~ Type * N + (1 | location), data=data.fat))

# pairwise comparisons (Comparisons in Main Text)
data.BC <- data.fat[is.element(data.fat$Type,c("Biological","Chemical")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BC))
summary(lm(sh_0 ~ Type, data=data.BC)) # Chem part of first sentence of main result

data.BP <- data.fat[is.element(data.fat$Type,c("Biological","Physical")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BP))
summary(lm(sh_0 ~ Type, data=data.BP)) # Phys part of first sentence of main result

data.BM <- data.fat[is.element(data.fat$Type,c("Biological","Meteorological")),]
Anova(lmer(sh_0 ~ Type + N + (1 | location), data=data.BM))
summary(lm(sh_0 ~ Type, data=data.BM)) # Met part of first sentence of main result


# Second part of main result statement in main text
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.BC))
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.BP))
summary(lmer(sh_0 ~ Type + N + (1 | location), data=data.BM))


