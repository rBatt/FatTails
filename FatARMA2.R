# Fat ARMA from Tony Ives on 22-Sapt-2014
# Uses the arima() function in R, not Tony's
# Also, fits to full time series, not annual maxima like previous ARMA analysis

library(plyr)

source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data_Functions.R") # for tony.yearly.max, fill.Full
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R") # for fillMiss
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") # for gev.fit2, fitARMA2


load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/bad.no.se.names.RData")

d <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fullTimeSeries_4Tony.txt", header=TRUE)

bad.no.se.names.index <- paste(d[,"taxID"], d[,"location"], d[,"variable"]) %in% bad.no.se.names

bad.key <- paste(as.numeric(d[,"taxID"]), as.numeric(d[,"location"]), as.numeric(d[,"variable"]))[bad.no.se.names.index]

dd <- d



dd3 <- ddply(dd, c("Type","taxID","location","variable"), fill.Full, .progress="text")



z0 <- ddply(dd3, c("Type","taxID","location","variable"), fitARMA2)


# remove time series that don't have se's
# doing this after the ARMA fitting because Tony had to manually select the right orders for some of the time series
# so that the ARMA fits were stationary
# and in order to preserve his counting of those time series, I'll just remove the no-se time series afterwards
# z <- z[!(paste(z[,"taxID"], z[,"location"], z[,"variable"]) %in% bad.key),]
z <- z0[!(paste(z0[,"taxID"], z0[,"location"], z0[,"variable"]) %in% bad.no.se.names),]
z.full <- z
z[which.max(z[,"xi.resid"]),]
z <- z[z[,"xi.resid"]<3,]

z[,"Type"] <- factor(z[,"Type"], levels=c("Biological", "Chemical", "Physical", "Meteorological"))


save(z0, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z0.RData")
save(z, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.RData")
save(z.full, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.full.RData")



summary(lm(xi.resid ~ xi2*type, data=z))
# Coefficients:
             # Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.007243   0.017026  -0.425 0.670704    
# xi           0.808744   0.032770  24.679  < 2e-16 ***
# type2        0.077891   0.021664   3.595 0.000351 ***
# type3        0.013259   0.065263   0.203 0.839078    
# type4        0.022811   0.028669   0.796 0.426541    
# xi:type2     0.019840   0.048314   0.411 0.681481    
# xi:type3     0.004898   0.295698   0.017 0.986790    
# xi:type4    -0.311038   0.071018  -4.380 1.41e-05 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Residual standard error: 0.1916 on 589 degrees of freedom
# Multiple R-squared:  0.7168,	Adjusted R-squared:  0.7134 
# F-statistic: 212.9 on 7 and 589 DF,  p-value: < 2.2e-16

library(car)
Anova(lm(xi.resid ~ xi2*type, data=z), type=3)
             # Sum Sq  Df  F value    Pr(>F)    
# (Intercept)  0.0066   1   0.1810  0.670704    
# xi          22.3521   1 609.0768 < 2.2e-16 ***
# type         0.5193   3   4.7171  0.002911 ** 
# xi:type      0.8336   3   7.5716 5.632e-05 ***
# Residuals   21.6153 589  

Anova(lm(xi2 ~ xi2*type, data=z), type=3) # new 30-Sept-2014 Ryan
                     

Anova(lm(xi.resid ~ xi2*type+I(p+q)+lambda, data=z), type=3)
            # Sum Sq  Df  F value    Pr(>F)    
# (Intercept)  0.0001   1   0.0032  0.954571    
# xi          22.3132   1 607.0846 < 2.2e-16 ***
# type         0.5212   3   4.7265  0.002874 ** 
# p            0.0036   1   0.0975  0.754920    
# xi:type      0.8290   3   7.5179 6.069e-05 ***
# Residuals   21.6118 588                       

Anova(lm(xi2 ~ I(p+q)+lambda+Type, data=z), type=3)
            # Sum Sq  Df F value    Pr(>F)    
# (Intercept)  8.409   1 68.6572 7.931e-16 ***
# type         1.835   3  4.9940  0.001989 ** 
# p            0.064   1  0.5189  0.471612    
# type:p       0.344   3  0.9368  0.422497    
# Residuals   72.137 589                      

Anova(lm(xi.resid ~ type+lambda*I(p+q), data=z), type=3)
            # Sum Sq  Df F value    Pr(>F)    
# (Intercept)  6.623   1 59.2109 5.995e-14 ***
# type         1.392   3  4.1485  0.006345 ** 
# p            0.259   1  2.3155  0.128623    
# type:p       0.254   3  0.7577  0.518180

# ============================
# = New regression summaries =
# ============================
 summary(lm(xi.resid ~ xi*Type, data=z))
#                       Estimate Std. Error t value Pr(>|t|)    
# (Intercept)           -0.02201    0.01716  -1.283   0.2001    
# xi                     0.85712    0.03390  25.284  < 2e-16 ***
# TypeChemical           0.10033    0.02157   4.650  4.1e-06 ***
# TypePhysical           0.04834    0.02862   1.689   0.0917 .  
# TypeMeteorological     0.02807    0.06412   0.438   0.6618    
# xi:TypeChemical       -0.04753    0.04783  -0.994   0.3207    
# xi:TypePhysical       -0.34987    0.07646  -4.576  5.8e-06 ***
# xi:TypeMeteorological -0.04369    0.29015  -0.151   0.8803

summary(lm(xi2~Type+I(p+q)+lambda, data=z))
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.32664    0.02414  13.531  < 2e-16 ***
# TypeChemical       -0.20132    0.02545  -7.911 1.27e-14 ***
# TypePhysical       -0.34689    0.03582  -9.683  < 2e-16 ***
# TypeMeteorological -0.24604    0.08143  -3.021  0.00263 ** 
# I(p + q)            0.02696    0.01398   1.928  0.05434 .  
# lambda             -0.30993    0.06746  -4.594 5.32e-06 ***



summary(lm(xi2~xi.resid+Type+lambda, data=z))
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)         0.14080    0.02135   6.596 9.40e-11 ***
# xi.resid            0.49980    0.02692  18.566  < 2e-16 ***
# TypeChemical       -0.14680    0.02048  -7.169 2.28e-12 ***
# TypeMeteorological -0.10567    0.06531  -1.618 0.106191    
# TypePhysical       -0.20022    0.02957  -6.771 3.11e-11 ***
# lambda             -0.10386    0.02830  -3.671 0.000264 ***


summary(lm(xi.resid~xi2*Type+lambda, data=z))
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.15248    0.02838   5.372 1.12e-07 ***
# xi2                     0.70637    0.05824  12.129  < 2e-16 ***
# TypeChemical            0.03110    0.02760   1.127   0.2602    
# TypeMeteorological     -0.09731    0.09187  -1.059   0.2899    
# TypePhysical           -0.09648    0.04076  -2.367   0.0183 *  
# lambda                 -0.03803    0.03457  -1.100   0.2718    
# xi2:TypeChemical        0.15756    0.08226   1.915   0.0559 .  
# xi2:TypeMeteorological  0.19108    0.44618   0.428   0.6686    
# xi2:TypePhysical       -0.28981    0.12663  -2.289   0.0225 *


summary(lm(xi.resid~xi2*Type+I(p+q)+lambda, data=z))
#                        Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.14724    0.02887   5.100  4.6e-07 ***
# xi2                     0.70284    0.05835  12.045  < 2e-16 ***
# TypeChemical            0.03023    0.02761   1.095   0.2740    
# TypeMeteorological     -0.09440    0.09192  -1.027   0.3048    
# TypePhysical           -0.09806    0.04079  -2.404   0.0165 *  
# I(p + q)                0.01337    0.01349   0.991   0.3219    
# lambda                 -0.09394    0.06615  -1.420   0.1561    
# xi2:TypeChemical        0.15724    0.08227   1.911   0.0564 .  
# xi2:TypeMeteorological  0.20958    0.44658   0.469   0.6390    
# xi2:TypePhysical       -0.28644    0.12667  -2.261   0.0241 *


##############################################
## Single example
##############################################

i1 <- "Biological" #1
i2 <- "Diaptomus" #24
i3 <- "ME" #8
i4 <- "density"# 9
c(levels(d$Type)==i1,levels(d$taxID)==i2,levels(d$location)==i3,levels(d$variable)==i4)

index <- (dd$Type == i1 & dd$taxID == i2 & dd$location == i3 & dd$variable == i4)
x <- dd$Data[index]
year <- dd$year4[index]
u <- 1:length(x)
y <- lm(x ~ u)$residuals
y <- y/sd(y)

X <- y
Y <- tony.yearly.Max(X,year)[,1]
ii <- tony.yearly.Max(X,year)[,2]
xi <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

a <- list()
a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
a[[3]] <- arima(y,order=c(3,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
a[[4]] <- arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
a[[5]] <- arima(y,order=c(3,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
a[[6]] <- arima(y,order=c(3,0,2), method="CSS-ML", optim.control=list(maxit=10^3))
pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic,a[[4]]$aic,a[[5]]$aic,a[[6]]$aic))[1]
				
XX <- a[[pick]]$residuals
YY <- tony.yearly.Max(XX,year)[,1]
iii <- tony.yearly.Max(XX,year)[,2]

xi.resid <- gev.fit2(xdat=YY, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]

c(i1,i2,i3,i4, length(x), a[[pick]]$ar[1:2], a[[pick]]$loglik, a[[pick]]$aic, xi, xi.resid)
#   1.0000000   24.0000000    8.0000000    9.0000000  241.0000000    1.0000000    0.0000000 -272.6432545  551.2865090    0.4438000    0.3313518 


save(X, Y, ii, XX, YY, iii, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/Diaptomus.eg.plot.RData")




# =======================================
# = Double check with original analysis =
# =======================================
# load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
# plot(test[,c("xi2","sh_0")])
# plot(test[,c("xi2.se","se.sh_0")])



# =========================
# = To go into manuscript =
# =========================
summary(lm(xi2~Type+I(p+q)+lambda, data=z, weights=1/xi2.se^2))
dev.new(); par(mfrow=c(2,2), mar=c(2,2,2,0.5), ps=10, cex=1, mgp=c(0.75,0.15,0), tcl=-0.15, family="Times")
plot(lm(xi2~Type+I(p+q)+lambda, data=z, weights=1/xi2.se)) # these are surprisingly amazing diagnostics

# summary(lm(xi2~Type+I(p+q)+lambda+xi.resid, data=z, weights=1/xi2.se)) # shows that order isn't significant
summary(lm(xi2~Type+I(p+q)+lambda+xi.resid, data=z, weights=1/xi2.se^2))
dev.new(); par(mfrow=c(2,2), mar=c(2,2,2,0.5), ps=10, cex=1, mgp=c(0.75,0.15,0), tcl=-0.15, family="Times")
plot(lm(xi2~Type+lambda+xi.resid, data=z, weights=1/xi2.se^2)) # these are surprisingly amazing diagnostics


summary(lm(xi2~Type+I(p+q)+lambda+xi.resid, data=z, weights=1/xi2.se^2))$coef


# =====================================================
# = Why I excluded 1 ts due to crazy high residual xi =
# =====================================================
summary(lm(xi2~Type+lambda+xi.resid, data=z, weights=1/xi2.se^2))
dev.new(); par(mfrow=c(2,2), mar=c(2,2,2,0.5), ps=10, cex=1, mgp=c(0.75,0.15,0), tcl=-0.15, family="Times")
plot(lm(xi2~Type+I(p+q)+lambda+xi.resid, data=z0, weights=1/xi2.se^2)) # above regression, but w/o the outlier removed


