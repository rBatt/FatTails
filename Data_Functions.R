
# ===========================
# = My version of reshape() =
# ===========================
reshape2 <- function(...){
	a <- reshape(...)
	row.names(a) <- NULL
	a <- a[,!names(a)=="id"]
	a
}


# ================================================
# = Sum w/ removing NA's, except when all are NA =
# ================================================
Zoop_Sum <- function(x){
	if(all(is.na(x))){
		NA
	}else{
		sum(x, na.rm=TRUE)
	}
}


# ====================
# = Epilimnetic Mean =
# ====================
# Used for the physical data, which already has zmix calculated for it
# Because zmix was calculated for this data, can rely on consistency between zmix and variable
EpiMean <- function(x){
	xnames <- names(x)[!is.element(names(x), c("lakeid", "year4", "daynum", "sampledate"))]
	colMeans(x[,xnames], na.rm=TRUE)
}

# ====================
# = Epilimnetic Mean =
# ====================
# For other variables, need a function for Epi Mean that
# Computes a "surface mean" even when no zmix is available (value of shallowest non-NA depth)
# 
EpiMean2 <- function(x){
	stopifnot(all(is.element(c("depth","Zmix"), colnames(x))))
	xnames <- colnames(x)[!is.element(colnames(x), c("lakeid", "year4", "daynum", "sampledate"))]
	x2 <- x[,xnames]

	if(all(is.na(x2[,"Zmix"]))){ #if zmix is unknown, just return the shallowest values
		if(is.null(nrow(x2))){return(data.frame(x2))} #if zmix is unknown and there is only one depth of values, just return this depth
		
		xnames2 <- colnames(x2)[!is.element(colnames(x2), c("depth", "Zmix"))]
		if(all(apply(X=x2[,xnames2], MARGIN=2, FUN=is.na))){return(colMeans(x2, na.rm=TRUE))} #if every variable is NA for all depths, just return the mean of all columns (shortcut for getting the mean depth and retaining the format of the data frame)
		shallDepths <- c()
		ShalVal <- matrix(data=NA, ncol=length(xnames), dimnames=list(NULL, xnames))
		for(i in 1:length(xnames2)){ #this elaborate loop is necessary b/c the shallowest non-NA value might not be the same for each column
			x_notNA_ind <- which(!is.na(x2[,xnames2[i]]))
			if(length(x_notNA_ind)==0){ #if all of the values in this column (i.e., this variable) are NA, record an NA for the depth and for the variable
				shallDepths[i] <- NA
				ShalVal[,xnames2[i]] <- NA
			}else{ #if this column contains non-NA values, record the shallowest depth at which the varialbe is non-NA, and record the variable's value at the aforementioned depth
				x_notNA <- x2[x_notNA_ind,]
				shallDepths[i] <- min(x_notNA[, "depth"]) # I'm intentionally leaving out na.rm=TRUE b/c none of these values should be NA
				ShalVal[,xnames2[i]] <- x_notNA[which.min(x_notNA[, "depth"]), xnames2[i]]
			}
		}
		UniqueShals <- unique(na.omit(shallDepths))
		if(length(UniqueShals)==1){
			DeepestShall <- UniqueShals
		}else{
			DeepestShall <- max(UniqueShals) #paste("<=", max(UniqueShals))
		}
		
		ShalVal[,"depth"] <- DeepestShall
		# ShalVal2 <- matrix(ShalVal, ncol=length(xnames), dimnames=list(NULL, xnames))
		return(ShalVal) #if zmix is unknown and values are at multiple depths, 
		# shall <- which.min(x2[,"depth"])
		# x3 <- colMeans(x2[shall,], na.rm=TRUE)
	}
	
	zd <- x2[,"Zmix"] - x2[,"depth"]
	InEpi <- as.numeric(zd>0 & !is.na(zd)) #if zmix is known, and values were from epi, return the mean of the epi
	if(any(InEpi==1)){
		return(colMeans(x2[which(InEpi==1),], na.rm=TRUE))
	}
	return(colMeans(x2, na.rm=TRUE)) #if zmix is known, but no values were in the epi, just return the mean of all values that were recorded.
}




# ========
# = Zmix =
# ========
zmix <- function(x){
		temp <- x[,"wtemp"]
		depth <- x[,"depth"]
        tempDiff <- temp[-length(temp)] - temp[-1]
        depthDiff <- depth[-1] - depth[-length(depth)]
        ratio <- tempDiff / depthDiff
		zhat <- depth[which(ratio >= 2)][1]
		# Zhat2 <- min(which(zhat>0))
		# Z <- Zhat2
		Z <- ifelse(!(zhat>0), NA, zhat)
       return(Z)
}



# ==================================================
# = If a variable is flagged, turn its value to NA =
# ==================================================
flag2NA <- function(x, summarizeBad=TRUE){
	flagNames <- names(x)[regexpr("^flag", names(x))==1] #which column names start with 'flag' ?
	flagVars <- gsub("flag", "", flagNames) #what are the column names that have a corresponding "flag___" column?
	hasFlag <- function(x){ #function to detect which elements contain an uppercase character (i.e., which have a flag value)
		grep("[:ALPHA:]", x)
	}
	bad <- apply(x[,flagNames], 2, hasFlag) #for each flag column, which rows are flagged?
	if(summarizeBad){
		print(summary(bad))	
	}
	for(i in 1:length(flagVars)){
		tbadRow <- do.call('$', list(bad, flagNames[i])) #which rows are bad in the ith flag column?
		if(length(tbadRow)>0){
			x[tbadRow, flagVars[i]] <- NA #replace those flagged elements in the parent column with NA
		}
	}
	x #return the data frame with flagged values replaced with NA's
}




# ==============================
# = Changing Inf or -Inf to NA =
# ==============================
Inf2NA <- function(x) {x[x==-Inf | x==Inf] <- NA; x}



# =========================================================
# = Graphically inspect time series and remove bad values =
# =========================================================
manClean <- function(x, varCols){
	x2 <- NULL
	dev.new(width=10, height=4)
	par(mar=c(3,3,3,0.5))
	for(i in 1:length(varCols)){
		tvar <- varCols[i]
		for(l in 1:length(unique(x[,"lakeid"]))){
			tlake <- unique(x[,"lakeid"])[l]
			tindex <- x[,"lakeid"]==tlake
			tdat <- x[tindex,tvar]
			if(any(!is.na(tdat))){
				epiInfo <- is.element(c("depth", "Zmix"), names(x))
				if(all(epiInfo)){ #if there is information about the epi depth, color the points red for epi, and NA otherwise
					MyRed <- rgb(t(col2rgb("red")), alpha=40, maxColorValue=255)
					inEpi <- x[tindex,"depth"]<=x[tindex,"Zmix"] | is.na(x[tindex,"Zmix"]) & x[tindex,"depth"] <=2
					plotCols <- c(NA, MyRed)[inEpi+1]
				}else{
					plotCols <- rgb(t(col2rgb("blue")), alpha=40, maxColorValue=255) #if the necessary columns for determining epi or not are missing, color blue
				}
				plot(tdat, main=paste(tvar, tlake), pch=21, bg=plotCols)
				abline(h=0)
				bad <- identify(tdat)
				if(length(bad)>0){
					x2 <- c(x2, paste(tlake, tvar, x[x[,"lakeid"]==tlake,tvar][bad]))
					x[x[,"lakeid"]==tlake,tvar][bad] <- NA
				}
			}else{
				next
			}
		}
	}
	dev.off()
	return(list(x, x2))	
}


# =======================
# = Summarize fish size =
# =======================
SizSumry <- function(x){
	Nfish <- length(x[,"spname"])
	
	exLeng <- x[,"length"][!is.na(x[,"length"])]
	exWei <- x[,"weight"][!is.na(x[,"weight"])]
	
	Nleng <- length(exLeng)
	Nwei <- length(exWei)
	
	if(Nleng==0){
		maxLeng <- NA
		minLeng <- NA
	}else{
		maxLeng <- max(exLeng, na.rm=TRUE)
		minLeng <- min(exLeng, na.rm=TRUE)
	}
	
	if(Nwei==0){
		maxWei <- NA
		minWei <- NA
	}else{
		maxWei <- max(exWei, na.rm=TRUE)
		minWei <- min(exWei, na.rm=TRUE)
	}
	
	sumryLeng <- c("mean_"=mean(exLeng, na.rm=TRUE), "max_"=maxLeng, "min_"=minLeng) #summary(x[,"length"])[[c(1,4,6)]]
	names(sumryLeng) <- paste(names(sumryLeng), "Leng", sep="")
	sumryWei <- c("mean_"=mean(exWei, na.rm=TRUE), "max_"=maxWei, "min_"=minWei) #summary(x[,"weight"])[[c(1,4,6)]]
	names(sumryWei) <- paste(names(sumryWei), "Wei", sep="")
	
	x2names <- c(names(sumryLeng), names(sumryWei), "Nfish", "Nleng", "Nwei")
	# x2 <- matrix(data=c(sumryLeng, sumryWei, Nfish, Nleng, Nwei), nrow=1, dimnames=list(NULL, x2names))
	x2 <- c(sumryLeng, sumryWei, Nfish, Nleng, Nwei)
	names(x2) <- x2names
	return(x2)
}


# =================
# = Fix ice dates =
# =================
FixIceDates <- function(x){
	on_NA <- is.na(x[,"iceon"])
	x_on_NoNA <- x[!on_NA, c("lakeid", "season", "iceon", "ice_on", "ice_duration", "year4")]
	off_NA <- is.na(x[,"iceoff"])
	x_off_NoNA <- x[!off_NA, c("lakeid", "season", "iceoff", "ice_off", "ice_duration", "year4")]
	
	take4 <- function(x) paste(x[1:4], collapse="")
	IceOn_Year <- as.integer(unlist(lapply(strsplit(as.character(x_on_NoNA[,"iceon"]), split=""), take4)))
	IceOff_Year <- as.integer(unlist(lapply(strsplit(as.character(x_off_NoNA[,"iceoff"]), split=""), take4)))
	
	ensure2_month <- function(x){sprintf("%02s", x[1])}
	ensure2_day <- function(x){sprintf("%02s", x[2])}
	IceOn_Month <- unlist(lapply(strsplit(as.character(x_on_NoNA[,"ice_on"]), split="/"), ensure2_month))
	IceOn_Day <- unlist(lapply(strsplit(as.character(x_on_NoNA[,"ice_on"]), split="/"), ensure2_day))
	IceOff_Month <- unlist(lapply(strsplit(as.character(x_off_NoNA[,"ice_off"]), split="/"), ensure2_month))
	IceOff_Day <- unlist(lapply(strsplit(as.character(x_off_NoNA[,"ice_off"]), split="/"), ensure2_day))
	
	Ice_On <- paste(IceOn_Year, IceOn_Month, IceOn_Day, sep="-")
	Ice_Off <- paste(IceOff_Year, IceOff_Month, IceOff_Day, sep="-")
	
	x_on_NoNA[,"ice_on"] <- Ice_On
	x_off_NoNA[,"ice_off"] <- Ice_Off
	x2 <- merge(x_on_NoNA, x_off_NoNA, all=TRUE)
	return(x2[,c("lakeid", "year4", "ice_on", "ice_off", "ice_duration")])
}






# =======================
# = Calculate days open =
# =======================
# "Open" meaning ice-free
CalcDaysOpen <- function(x){
	OnYears <- x[,"year4"]
	for(i in 1:length(OnYears)){
		OffRow <- which(x[,"year4"]==(OnYears[i]-1))
		if(length(OffRow)==0 || is.na(x[OffRow,"ice_off"]) || is.na(x[i,"ice_on"])){ #should maybe set to length(OffRow)!=1 || ...x
			x[i,"DaysOpen"] <- NA
		}else{
			x[i,"DaysOpen"] <- as.integer(difftime(as.Date(x[i, "ice_on"]), as.Date(x[OffRow,"ice_off"])))
		}
	}
	return(x)
}





# =========================
# = Subset to Genus level =
# =========================
# Remove species, family, order, class, phylum
sub.gen <- function(x){
	if(!any(colnames(x)=="taxLvl")){
		stop("no taxLvl column")
	}
	x[!x[,"taxLvl"]%in%c("Species","Family","Order","Class","Phylum"),]
}





# ======================================================================================
# Series of yearly maxima
# ======================================================================================
tony.yearly.Max <- function(X=NULL, t=NULL){
    years <- unique(t)
	n.year <- length(years)
  
    Y <- array(0,c(n.year,1))
    index <- Y
    i <- 0
    counter <- 0
    for(y in years){
	    i <- i+1
    	# Y[i] <- ifelse(all(!is.finite(X[t==t])), NA, max(X[t==y], na.rm=TRUE))
		if(all(!is.finite(X[t==y]))){
			Y[i] <- NA
		}else{
			Y[i] <- max(X[t==y], na.rm=TRUE) 	
		}
    	index[i] <- counter + order(X[t==y],decreasing = T)[1]
    	counter <- counter + length(X[t==y])
    }
  	return(cbind(Y,index))
}




