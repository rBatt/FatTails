





library(lme4)
library(multcomp)





contr <- c("Chemical-Biological=0","Physical-Biological=0","Meteorological-Biological=0")

# ============================================================
# = Drop Meteorological, Compare Type|Location, with weights =
# ============================================================
data.fat.noMet <- data.fat[data.fat[,"Type"]!="Meteorological",]

# Compare a few w/o weights
fm4 <- ((lmer(sh_0~Type+N+(Type|location), data=data.fat.noMet)))
summary(glht(fm4, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0"))))

# Do it with weights
xiWeights.noMet <- 1/(data.fat.noMet[,"se.sh_0"]^2)
fm4.w <- ((lmer(sh_0~Type+N+(Type|location), weights=xiWeights.noMet, data=data.fat.noMet)))
summary(glht(fm4.w, linfct=mcp(Type=c("Chemical-Biological=0","Physical-Biological=0"))))



# ===================================
# = With Met, after Consulting Tony =
# ===================================
# tried again after consulting with Tony
re1 <- lmer(sh_0 ~ Type + N + (1|location), data=data.fat) # random effect model 1
summary(glht(re1, linfct=mcp(Type=contr)))

re1.w <- lmer(sh_0 ~ Type + N + (1|location), data=data.fat, weights=1/(data.fat[,"se.sh_0"]^2)) # random effect w/ weights
summary(glht(re1.w, linfct=mcp(Type=contr)))

re2 <- lmer(sh_0 ~ Type + N + (N|location), data=data.fat) # random effect model 2
summary(glht(re2, linfct=mcp(Type=contr)))

re2.w <- lmer(sh_0 ~ Type + N + (N|location), data=data.fat, weights=1/(data.fat[,"se.sh_0"]^2)) # re2 with weights
summary(glht(re2.w, linfct=mcp(Type=contr)))





