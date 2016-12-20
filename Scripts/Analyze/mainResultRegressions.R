
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
regS2 <- lmer(sh_0 ~ Type + N + (1 | location), data=data.fat)
summary(regS2)

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

# =======
# = R2R =
# =======
library(multcomp)
regS2 <- lmer(sh_0 ~ Type + N + (1 | location), data=data.fat)
contr <- c("Chemical-Biological=0", "Physical-Biological=0" ,"Meteorological-Biological=0")
summary(multcomp::glht(regS2, linfct=mcp(Type=contr))) # , test=adjusted("bonferroni") # ~same
#
# 	 Simultaneous Tests for General Linear Hypotheses
#
# Multiple Comparisons of Means: User-defined Contrasts
#
#
# Fit: lmer(formula = sh_0 ~ Type + N + (1 | location), data = data.fat)
#
# Linear Hypotheses:
#                                  Estimate Std. Error z value Pr(>|z|)
# Chemical - Biological == 0       -0.25328    0.02536  -9.987   <1e-04 ***
# Physical - Biological == 0       -0.40655    0.03566 -11.401   <1e-04 ***
# Meteorological - Biological == 0 -0.27256    0.10267  -2.655   0.0234 *
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# (Adjusted p values reported -- single-step method)


