# Fat ARMA from Tony Ives on 22-Sapt-2014
# Uses the arima() function in R, not Tony's
# Also, fits to full time series, not annual maxima like previous ARMA analysis


source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data_Functions.R") # for tony.yearly.max
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/FatTails_Functions.R") # for fillMiss
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/ARMAFunctions.R") # for gev.fit2


load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/bad.no.se.names.RData")

d <- read.table("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fullTimeSeries_4Tony.txt", header=TRUE)

bad.no.se.names.index <- paste(d[,"taxID"], d[,"location"], d[,"variable"]) %in% bad.no.se.names

bad.key <- paste(as.numeric(d[,"taxID"]), as.numeric(d[,"location"]), as.numeric(d[,"variable"]))[bad.no.se.names.index]

dd <- d

# x <- dd[dd[,"variable"]=="chlor"&dd[,"location"]=="AL",]
test <- dd[dd[,"variable"]=="chlor"&dd[,"location"]=="TB",]
fill.Full <- function(x){
	require(zoo)
	require(plyr)
	
	x <- x[order(x[,"year4"], x[,"daynum"]),]
	
	xyr <- x[,"year4"]
	tyr <- table(xyr)
	max.obs <- max(tyr, na.rm=TRUE)
	max.yrs <- names(tyr)[tyr==max.obs]
	
	# =======================
	# = Handle simple cases =
	# =======================
	if(all(tyr==max.obs)){ # if all are the same, don't bother with complexities
		x2 <- x
		x2[,"rank"] <- 1:max.obs
		
		possYears <- do.call(":", as.list(range(xyr, na.rm=TRUE)))
		possRanks <- 1:max.obs
		refFrame0 <- expand.grid(possRanks, possYears)
		refFrame <- data.frame("Type"=unique(x2[,"Type"]), "taxID"=unique(x2[,"taxID"]), "location"=unique(x2[,"location"]), "variable"=unique(x2[,"variable"]), "year4"=refFrame0[,2], "rank"=refFrame0[,1])
	
		fillFrame0 <- merge(x2, refFrame, all=TRUE)
		fillFrame <- fillFrame0[,c("Type", "taxID","location","variable", "year4", "daynum", "rank", "n.yr", "Data")]
	
		return(fillFrame)
	}
	
	if(max.obs>60){ # if there's more than 1 a week, just fill in to daily resolution
		x2 <- x
		x2[,"rank"] <- x2[,"daynum"]
		
		possYears <- do.call(":", as.list(range(xyr, na.rm=TRUE)))
		possRanks <- 1:365
		refFrame0 <- expand.grid(possRanks, possYears)
		refFrame <- data.frame("Type"=unique(x2[,"Type"]), "taxID"=unique(x2[,"taxID"]), "location"=unique(x2[,"location"]), "variable"=unique(x2[,"variable"]), "year4"=refFrame0[,2], "rank"=refFrame0[,1])
	
		fillFrame0 <- merge(x2, refFrame, all=TRUE)
		fillFrame <- fillFrame0[,c("Type", "taxID","location","variable", "year4", "daynum", "rank", "n.yr", "Data")]
	
		return(fillFrame)
		
	}
	
	if(max.obs==1){ # if there's only 1 observation per year, just fill in years
		x2 <- x
		x2[,"rank"] <- 1
		
		possYears <- do.call(":", as.list(range(xyr, na.rm=TRUE)))
		possRanks <- 1
		refFrame0 <- expand.grid(possRanks, possYears)
		refFrame <- data.frame("Type"=unique(x2[,"Type"]), "taxID"=unique(x2[,"taxID"]), "location"=unique(x2[,"location"]), "variable"=unique(x2[,"variable"]), "year4"=refFrame0[,2], "rank"=refFrame0[,1])
	
		fillFrame0 <- merge(x2, refFrame, all=TRUE)
		fillFrame <- fillFrame0[,c("Type", "taxID","location","variable", "year4", "daynum", "rank", "n.yr", "Data")]
	
		return(fillFrame)
	}
	
	# ================================
	# = Handle more complicated case =
	# ================================
	max.dat <- x[x[,"year4"]%in%max.yrs,]	
	mu.day <- rollapply(sort(max.dat[,"daynum"]), mean, width=length(max.yrs), by=length(max.yrs))
	
	x.ranked00 <- cbind(x,rank00=findInterval(x[,"daynum"], c(mu.day), all.inside=TRUE))
	
	# plot(blah[,"rank"], blah[,"daynum"])
	# abline(h=rollapply(sort(max.dat[,"daynum"]), mean, width=3, by=3))
	
	x.ranked0 <- ddply(x.ranked00, "year4", function(x){x[,"rank0"] <- x[,"rank00"]+cumsum(duplicated(x[,"rank00"])); x})
	
	
	demote <- function(x){
		xr2 <- x[,"rank0"]
		excess <- max(xr2) - max.obs


		if(excess>0){
			rankBump <- integer(length(xr2))
			# rankDemote.index <- c()
			# rankDemote <- integer(length(xr2))
			for(i in 1:excess){
				rankDiffs <- c(NA, diff(xr2))
				dayDiffs <- c(NA, diff(x[,"daynum"]))
				rank.per.day <- dayDiffs/rankDiffs
				
				# updating xr2 iteratively just in case a certain index represented a >2 jump in rank, and even after being demoted once, deserves a 2nd demotion because its rank.per.day was still the smallest after 1st (etc) demote
				rankDemote.index <- which.min(rank.per.day) #which(rank.per.day <= sort(rank.per.day)[i])[i]
				rankDemote <- cumsum(seq_along(xr2) == rankDemote.index)
				xr2 <- xr2 - rankDemote
			}
		}
		
		x[,"rank"] <- xr2
		x
	}
	
	x2 <- ddply(x.ranked0, "year4", demote)
	
	
	possYears <- do.call(":", as.list(range(xyr, na.rm=TRUE)))
	possRanks <- 1:max.obs
	refFrame0 <- expand.grid(possRanks, possYears)
	refFrame <- data.frame("Type"=unique(x2[,"Type"]), "taxID"=unique(x2[,"taxID"]), "location"=unique(x2[,"location"]), "variable"=unique(x2[,"variable"]), "year4"=refFrame0[,2], "rank"=refFrame0[,1])
	
	fillFrame0 <- merge(x2, refFrame, all=TRUE)
	fillFrame <- fillFrame0[,c("Type", "taxID","location","variable", "year4", "daynum", "rank", "n.yr", "Data")]
	
	return(fillFrame)
	
}

# dd2 <- ddply(dd, c("Type","taxID","location","variable"), fill.Full, .progress="text")
dd3 <- ddply(dd, c("Type","taxID","location","variable"), fill.Full, .progress="text")



# "Type"	"taxID"	"location"	"variable"	"year4"	"n.yr"	"Data"
# dd$Type <- as.numeric(dd$Type)
# dd$taxID <- as.numeric(dd$taxID)
# dd$location <- as.numeric(dd$location)
# dd$variable <- as.numeric(dd$variable)

# dd$ID <- 1000000*dd$Type + 10000*dd$taxID + 100*dd$location + dd$variable

# dd$ID <- 1000000*as.numeric(dd$Type) + 10000*as.numeric(dd$taxID) + 100*as.numeric(dd$location) + as.numeric(dd$variable)



# c(length(unique(dd$ID[dd$Type==1])),length(unique(dd$ID[dd$Type==2])),length(unique(dd$ID[dd$Type==3])),length(unique(dd$ID[dd$Type==4])))
# colSums(table(dd$ID, dd$Type)>1)
# [1] 284 220  12  81

# length(unique(dd$ID))
# [1] 597 # after Batt removed 2 with no se's, down to 595

# w <- cbind(dd$Type,dd$ID,dd$year4,dd$Data)

# write.table(w,file="Batt_for_matlab_12Sep14",sep=",",row.names = F,col.names = F)

### using arima
# counter <- 0
# best.list <- matrix(0,nrow=597,ncol=12)
# best.list <- matrix(0,nrow=597,ncol=13)
# best.list <- data.frame("id"=NA, "type"=NA, "taxID"=NA, "location"=NA, "variable"=NA, "length"=NA, "p"=NA, "q"=NA, "loglik"=NA, "aic"=NA, "xi"=NA, "xi2"=NA, "xi.resid"=NA)

fitARMA2 <- function(DF){
	x <- DF[,"Data"]
	year <- DF[,"year4"]
	
	Y.gev000 <- fillMiss(DF)
	Y.gev00 <- tony.yearly.Max(Y.gev000[,"Data"], Y.gev000[,"year4"])
	Y.gev0 <- lm(Y.gev00[,1] ~ I(1:length(unique(Y.gev00[,2]))))$residuals
	Y.gev <- Y.gev0/sd(Y.gev0, na.rm=TRUE)
	
	u <- 1:length(x) # year-min(year)
	y <- lm(x ~ u)$residuals
	y <- y/sd(y)

	X <- y
	Y <- tony.yearly.Max(X,year[!is.na(x)])[,1]
	xi <- as.numeric(gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3])
	
	gev0 <- gev.fit2(xdat=Y.gev, show = FALSE, method = "Nelder-Mead", maxit = 10000)
	xi2 <- as.numeric(gev0$mle[3])
	xi2.se <- as.numeric(gev0$se[3])

	a <- list()
	if(length(y)>30){
		a[[1]] <- tryCatch(arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){NA})
		a[[2]] <- tryCatch(arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		a[[3]] <- tryCatch(arima(y,order=c(3,0,0), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		a[[4]] <- tryCatch(arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		a[[5]] <- tryCatch(arima(y,order=c(3,0,1), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		a[[6]] <- tryCatch(arima(y,order=c(3,0,2), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic,a[[4]]$aic,a[[5]]$aic,a[[6]]$aic))[1]
	}
	if(length(y)<=30){
		a[[1]] <- tryCatch(arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){NA})
		a[[2]] <- tryCatch(arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		a[[3]] <- tryCatch(arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3)), error=function(cond){list(aic=NA)})
		pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic))[1]
	}
	
	X <- a[[pick]]$residuals
	Y <- tony.yearly.Max(X,year[!is.na(x)])[,1]
	
	gev.resid0 <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)
	
	xi.resid <- as.numeric(gev.resid0$mle[3])
	xi.resid.se <- as.numeric(gev.resid0$se[3])
	
	
	
	best.pq <- a[[pick]]$ar[1:2]
	best.eig <- max(abs(Eig(a[[pick]]$coef[1:best.pq[1]])))
	data.frame("length"=length(x), "p"=best.pq[1], "q"=best.pq[2], "lambda"=best.eig, "loglik"=a[[pick]]$loglik, "aic"=a[[pick]]$aic, "xi"=xi, "xi2"=xi2, "xi.resid"=xi.resid, "xi2.se"=xi2.se, "xi.resid.se"=xi.resid.se)
	
	# best.list[counter,c("id","type","taxID","location","variable")] <- c(counter,i1,i2,i3,i4)
	# best.list[counter,c("length","p","q","lambda","loglik","aic","xi","xi2","xi.resid")] <- c(length(x), a[[pick]]$ar[1:2], Eig(a[[pick]]$) a[[pick]]$loglik, a[[pick]]$aic, xi, xi2, xi.resid)
	# 
	# return(best.list)
	
}
# fitARMA2(test)
# 
# for(i1 in unique(dd$Type)){
# 	w1 <- dd[dd$Type == i1,]
# 	for(i2 in unique(w1$taxID)){
# 		w2 <- w1[w1$taxID == i2,]
# 		for(i3 in unique(w2$location)){
# 			w3 <- w2[w2$location == i3,]
# 			for(i4 in unique(w3$variable)){
# 				counter <- counter+1
# 				
# 				if(counter<=600){
# 				x <- w3$Data[w3$variable == i4]
# 				year <- w3$year4[w3$variable == i4]
# 				
# 				
# 				# wtest <- w3[w3$variable==i4 & w3[w3$variable==i4,"year4"]!=2008L,]
# 				# 			Y.gev000 <- fillMiss(wtest) # works
# 				Y.gev000 <- fillMiss(w3[w3$variable==i4,])
# 				Y.gev00 <- tony.yearly.Max(Y.gev000[,"Data"], Y.gev000[,"year4"])
# 				Y.gev0 <- lm(Y.gev00[,1] ~ I(1:length(unique(Y.gev00[,2]))))$residuals
# 				Y.gev <- Y.gev0/sd(Y.gev0, na.rm=TRUE)
# 				
# 				u <- 1:length(x) # year-min(year)
# 				y <- lm(x ~ u)$residuals
# 				y <- y/sd(y)
# 
# 				X <- y
# 				Y <- tony.yearly.Max(X,year)[,1]
# 				xi <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
# 				xi2 <- gev.fit2(xdat=Y.gev, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
# 
# 				a <- list()
# 				if(!is.element(counter,c(53,54,152,163,221,254,255,330,345,353,369,382, 383,392,410,415,426,432,435,437,441,450,451,452,455,464,470,475,493,494))){
# 					if(length(y)>30){
# 						a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[3]] <- arima(y,order=c(3,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[4]] <- arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[5]] <- arima(y,order=c(3,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[6]] <- arima(y,order=c(3,0,2), method="CSS-ML", optim.control=list(maxit=10^3))
# 						pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic,a[[4]]$aic,a[[5]]$aic,a[[6]]$aic))[1]
# 					}
# 					if(length(y)<=30){
# 						a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 						a[[3]] <- arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 						pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic))[1]
# 					}
# 				}
# 				if(is.element(counter,c(53,54,163,441,451,475))){
# 					a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					pick <- order(c(a[[1]]$aic, a[[2]]$aic))[1]
# 				}
# 				if(is.element(counter,c(221))){
# 					a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					pick <- order(c(a[[1]]$aic))[1]
# 				}
# 				if(is.element(counter,c(152,254,255,330,345,353,369,382,383,392,410,432,450,452,494))){
# 					a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[3]] <- arima(y,order=c(3,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[4]] <- arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[5]] <- arima(y,order=c(3,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 					pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic,a[[4]]$aic,a[[5]]$aic))[1]
# 				}
# 				if(is.element(counter,c(415,426,435,437,455,464,470,493))){
# 					a[[1]] <- arima(y,order=c(1,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[2]] <- arima(y,order=c(2,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[3]] <- arima(y,order=c(3,0,0), method="CSS-ML", optim.control=list(maxit=10^3))
# 					a[[4]] <- arima(y,order=c(2,0,1), method="CSS-ML", optim.control=list(maxit=10^3))
# 					pick <- order(c(a[[1]]$aic, a[[2]]$aic,a[[3]]$aic,a[[4]]$aic))[1]
# 				}
# 				
# 				X <- a[[pick]]$residuals
# 				Y <- tony.yearly.Max(X,year)[,1]
# 				xi.resid <- gev.fit2(xdat=Y, show = FALSE, method = "Nelder-Mead", maxit = 10000)$mle[3]
# 
# 
# 				# out <- (c(counter,i1,i2,i3,i4, length(x), a[[pick]]$ar[1:2], a[[pick]]$loglik, a[[pick]]$aic, xi, xi2, xi.resid))
# 				# out <- data.frame(counter,i1,i2,i3,i4, length(x), a[[pick]]$ar[1:2], a[[pick]]$loglik, a[[pick]]$aic, as.numeric(xi), as.numeric(xi2), as.numeric(xi.resid))
# 				# best.list[counter,] <- out
# 				
# 				best.list[counter,c("id","type","taxID","location","variable")] <- c(counter,i1,i2,i3,i4)
# 				best.list[counter,c("length","p","q","loglik","aic","xi","xi2","xi.resid")] <- c(length(x), a[[pick]]$ar[1:2], a[[pick]]$loglik, a[[pick]]$aic, xi, xi2, xi.resid)
# 								
# 				# show(out)
# 				}
# 			}
# 		}
# 	}
# }
# 

# z <- best.list 
# z <- data.frame(best.list)
# names(z) <- c("id", "type", "taxID", "location", "variable", "length", "p", "q", "loglik", "aic", "xi", "xi2", "xi.resid")
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

save(z, file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.RData")


# par(mfrow=c(2,2))
# 
# # plot(xi ~ type, data=z)
# # plot(xi.resid ~ type, data=z)
# # 
# # beanplot(xi ~ type, data=z, what=c(1,1,1,0))
# # beanplot(xi.resid ~ type, data=z, what=c(1,1,1,0))
# 
# plot(xi2 ~ type, data=z)
# plot(xi.resid ~ type, data=z)
# 
# beanplot(xi2 ~ type, data=z, what=c(1,1,1,0))
# beanplot(xi.resid ~ type, data=z, what=c(1,1,1,0))

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
                   Estimate Std. Error t value Pr(>|t|)    
(Intercept)         0.14080    0.02135   6.596 9.40e-11 ***
xi.resid            0.49980    0.02692  18.566  < 2e-16 ***
TypeChemical       -0.14680    0.02048  -7.169 2.28e-12 ***
TypeMeteorological -0.10567    0.06531  -1.618 0.106191    
TypePhysical       -0.20022    0.02957  -6.771 3.11e-11 ***
lambda             -0.10386    0.02830  -3.671 0.000264 ***


# summary(lm(xi.resid~xi2*Type+lambda, data=z))
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


