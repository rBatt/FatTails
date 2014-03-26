# Run examples of fitting ARMA models
# RDB 17-Dec-2013


# =================
# = Notes to Tony =
# =================
# After setting working directory to the location of fat4tony.RData and fat4tony_function.R, this could will run as-is
# Note that this code will automatically download, install, and load 3 R packages (plyr, DEoptim, GenSA)

#-----
# sim1Final, chlor1Final are objects of interest
#-----

# if you want an example of the ARMA analysis applied for an array of model orders, set runCombs to TRUE.
# if you want to use lake time series that were log-transformed and then made stationary, set logStat to TRUE
# if you want to NOT see a graph of the time series, set vG to FALSE


rm(list=ls()) #remove objects from memory


# =================================================
# = Notes for RDB on how example data was created =
# =================================================
# setwd("/Users/battrd/Documents/School&Work/WiscResearch/FatTails")
# load("finalFrame_v1.RData")
# load("finalFrame_noXform.RData")
# egData_logStat <- rbind(finalFrame[finalFrame[,"variable"]=="chlor",])
# egData <- rbind(finalFrame_noXform[finalFrame_noXform[,"variable"]=="chlor",])
# setwd("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fat4tony/")

# ==============================
# = Install packages if needed =
# ==============================
if(!require("plyr")){install.packages("plyr")}
if(!require("DEoptim")){install.packages("DEoptim")}
if(!require("GenSA")){install.packages("GenSA")}

# ====================================
# = Load data, scripts, and packages =
# ====================================
load("fat4tony.RData") # this loads some of the data sets

source("fat4tony_functions.R") # this sources the R script that contains functions used in this script

library("plyr") # plyr package is used for conveniently applying a function to a subsets of data, then recombining results.
library("DEoptim") # Differential evolution form of optimization; method that I prefer. Set using Method="Evolve" in the OptAll and ARMAfit functions
library("GenSA") # Simulated annealing method of optimization; slightly faster than DE, but I haven't bothered to tune all the settings. Set using Method="Anneal"

# ============
# = OPTIONS! =
# ============
# OPTION 1
# Change the index to 1 to try multipe ARMA orders for chlorophyll 
runCombs <- c(TRUE, FALSE)[2] # option for exploring all ARMA orders. Takes 2 minutes to run on my machine. This is in the last part of the code

# OPTION 2
# Change the index to 1 to use log-xformed & stationary chlorophyll
# I'm not sure what stationary means anymore (thanks for that, Tony!), but here it means I removed a linear time trend
logStat <- c(TRUE, FALSE)[2]
if(logStat){
	egChlor <- egData_logStat
}else{
	egChlor <- egData
}

# OPTION 3
# Change the index to 2 to avoid viewing a graph of the time series
vG <- c(TRUE, FALSE)[1]

# ===========================
# = Some quick naming stuff =
# ===========================
bNames <- paste(rep("b",3), 1:3, sep="")
aNames <- paste(rep("a",3), 1:3, sep="")
names0 <- c(bNames, aNames)
namesOA <- c("P", "Q", "Lambda", "Period", names0, "nll", "AICc")

# =================================
# = Example of 1 fit w/ Fake Data =
# =================================
#Simulate data
set.seed(1) #set seed for random number generator
simDat <- as.numeric(arima.sim(model=list(ar=c(0.2, -0.9), ma=c(0.5)), n=100)) #simulate time series
simDat_x <- data.frame("location"="computer", "variable"="fake", "Data"=simDat) # put in a data.frame so compatible with functions that were designed for real data

#Fit coefficients of the specified model, and return Lambda, Period, nll, AICc in addition to coefficients
simFit1 <- OptAll(pq=c(2,1), Xs=simDat, Method="Evolve")[[1]] # run the "OptAll" function, which fits the specificed ARMA model to the data
simFit1 <- matrix(simFit1, nrow=1, dimnames=list(NULL, namesOA)) # label/name the elements returned by OptAll, and put the result in a matrix format

# Compute the sigmas (infinity, epsilon, E)
simFit1_x <- data.frame("location"="computer", "variable"="fake",simFit1) # put in a data.frame so compatible with functions that were designed for real d
simFit1_SE <- getSE(fat=simFit1_x, data=simDat_x)

# Combine the OptAll results and the getSE results
sim1Final <- merge(simFit1_x, simFit1_SE, all=TRUE)


# ==================================================
# = Example of 1 fit w/ Chlorophyll from Trout Bog =
# ==================================================
#The steps below are analagous to what was done for the simulated data above
chlorDat <- egChlor[egChlor[,"location"] == "TB",]

chlorFit1 <- OptAll(pq=c(2,0), Xs=chlorDat[,"Data"], Method="Evolve")[[1]]
chlorFit1 <- matrix(chlorFit1, nrow=1, dimnames=list(NULL, namesOA))

chlorFit1_x <- data.frame("location"="TB", "variable"="chlor",chlorFit1)
chlorFit1_SE <- getSE(fat=chlorFit1_x, data=chlorDat) # note that the SigE is way higher than the sigInf ... this can't be right, can it?

#combine the sigma results with the arma fit
chlor1Final <- merge(chlorFit1_x, chlorFit1_SE, all=TRUE)


if(runCombs | vG){
	egChlor_short <- egChlor[egChlor[,"location"] %in% c("CB", "TB"),] #get data from 2 lakes (more available)
}

if(runCombs){ # I am wrapping this in a logical statement so you don't accidentally run it ... it takes ~2 minutes on my machine
	# ============================================================
	# = Example of using model selection for multiple time seies =
	# ============================================================
	# apply the ARMAfit() function to all lake-variable combinations in the data set egChlor_short.  
	# the "dd" in the ddply means "input *d*ata frame, get *d*ata frame back"
	# In this case there is only one .variable ("chlor"), but putting this in there means that the function will return this column (convenient labeling)
	chlorFit <- ddply(egChlor_short, .variables=c("location", "variable"), .fun=ARMAfit, Method="Evolve", dName="Data", .progress="time")

	# apply the getSE function to the above arma results
	chlorFit_SE <- ddply(chlorFit, .variables=c("variable", "location", "P", "Q"), .fun=getSE, data=egChlor_short, .progress="time")

	# Combine the OptAll results and the getSE results
	chlorFinal <- merge(chlorFit, chlorFit_SE, all=TRUE)

	# for each lake-variable combination, select the row with the smallest AICc
	chlorFinal_best <- ddply(chlorFinal, .variables=c("location", "variable"), .fun=minaicc)
}

# =====================
# = Graph time series =
# =====================
if(vG){
	dev.new(width=3.4, height=6)
	par(mfrow=c(3,1), mar=c(2,3,0.5,0.5), oma=c(1,0,0,0), ps=9, cex=1, tcl=-0.25, mgp=c(1.55,0.5,0), family="Times", xpd=NA)
	plot(simDat_x[,3], type="l", xlab="", ylab="Simulated ARMA")
	plot(egChlor_short[egChlor_short[,"location"]=="CB",c("year4","Data")], xlab="", ylab="Crystal Bog Chlorophyll", type="l")
	plot(egChlor_short[egChlor_short[,"location"]=="TB",c("year4","Data")], xlab="Year/ Time Index", ylab="Trout Bog Chlorophyll", type="l")
}



