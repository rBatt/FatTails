
library("plyr")
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")

fatCat <- factor(sign(data.fat[,"shape.sig"]), labels=c("bounded","thin","fat"))
dfat2 <- cbind(data.fat[,c("Type","taxLvl","taxID","location","variable")],"fatCat"=fatCat)
dstat.fat <- merge(data.stat, dfat2)

d2 <- merge(dstat.fat, data.fat)
test <- randomForest(sh_0~Type+location+N+year4+Data+fatCat)

maxstat <- ddply(dstat.fat, c("Type","taxLvl","taxID","location","variable"), function(x)x[which.max(x[,"Data"]),])
rsyvl <- rowSums(table(maxstat[,c("year4","variable","location")]))
rsyvl_years <- as.numeric(names(rsyvl))

# ==============================================================================
# = Plot proportion maxima in each year, aggregated across lakes and variables =
# ==============================================================================
ym.tot <- matrix(c(as.integer(names(rsyvl[rsyvl_years>1980])), rsyvl[rsyvl_years>1980]), byrow=FALSE, ncol=2, dimnames=list(NULL,c("year4","nMax")))
nObs.tot <- table(data.stat[!is.na(data.stat[,"Data"]),"year4"])
yo.tot <- matrix(c(as.integer(names(nObs.tot)), nObs.tot), byrow=FALSE, ncol=2, dimnames=list(NULL, c("year4","nObs")))
om.tot <- merge(ym.tot, yo.tot)
plot(om.tot[,1], om.tot[,2]/om.tot[,3], type="l")



# =====================================================================================
# = Plot proportion maxima in each year, aggregated across variables (lakes separate) =
# =====================================================================================
yearsMax <- adply(table(maxstat[,c("year4","Type","location")]), 1:3)
names(yearsMax) <- c("year4","Type","location","nMax")
yearsMax[,"year4"] <- as.integer(as.character(yearsMax[,"year4"]))
yearsMax <- yearsMax[yearsMax[,"year4"]>1980,]
yearsMax <- yearsMax[yearsMax[,"Type"]%in%c("Physical","Chemical","Biological"),]
# yearsMax2 <- ddply(yearsMax, c("year4","Type","location"), function(x)sum(x[,"nMax"]))

nObs <- adply(table(data.stat[!is.na(data.stat[,"Data"]),c("year4","location","Type")]), 1:3)
nObs[,"year4"] <- as.integer(as.character(nObs[,"year4"]))
nObs <- nObs[nObs[,"year4"]>1980 & nObs[,"V1"]>=2,]
names(nObs) <- c("year4", "location", "Type", "nObs")

yom <- merge(nObs, yearsMax)
yom[,"propMax"] <- yom[,"nMax"]/yom[,"nObs"]
yom[yom[,"nMax"]==0&yom[,"nObs"]==0,"propMax"] <- 0

uT.ym <- c("Physical","Chemical","Biological")#unique(yearsMax[,"Type"]) # unique type (yearsMax)
uL.ym <- as.character(unique(yom[,"location"])) # unique location (yearsMax)

dev.new(width=3.5, height=7)
par(mfrow=c(3,1), mar=c(2,2,0.5,0.5))
lcols <- rainbow(n=length(uL.ym))
yLim <- c(0,max(yom[,"propMax"], na.rm=TRUE))
for(ut in 1:length(uT.ym)){
	tut <- as.character(yom[,"Type"])==uT.ym[ut]
	for(ul in 1:length(uL.ym)){
		tul <- as.character(yom[,"location"])==uL.ym[ul]
		tdat <- yom[tut&tul,]
		if(ul==1){
			plot(tdat[,"year4"], tdat[,"propMax"], type="l", col=lcols[ul], ylim=yLim)
		}else{
			lines(tdat[,"year4"], tdat[,"propMax"], col=lcols[ul])
		}
	if(ut==2){
		legend("top", legend=uL.ym, text.col=lcols, ncol=3)
	}	
	}
	tdat.tot <- aggregate(yom[tut,"propMax"], list(yom[tut,"year4"]), mean)
	lines(tdat.tot[,"Group.1"], tdat.tot[,"x"], col="black", lwd=3)
	
}





dev.new(width=8, height=6)
par(mfrow=c(4,3), mar=c(2,2,0.5,0.5))
# lcols <- rainbow(n=length(uL.ym))
tcols <- rainbow(n=3)

for(ul in 1:length(uL.ym)){
	tul <- as.character(yom[,"location"])==as.character(uL.ym[ul])
	yLim2 <- c(0,max(yom[tul,"propMax"], na.rm=TRUE))
	for(ut in 1:length(uT.ym)){
		tut <- as.character(yom[,"Type"])==uT.ym[ut]
		tdat <- yom[tut&tul,]

		if(ut==1){
			plot(tdat[,"year4"], tdat[,"propMax"], type="l", col=tcols[ut], ylim=yLim2)
			legend("top", legend=uL.ym[ul], ncol=1)
		}else{
			lines(tdat[,"year4"], tdat[,"propMax"], col=tcols[ut])
		}
	}
	tdat.tot <- aggregate(yom[tul,"propMax"], list(yom[tul,"year4"]), mean)
	lines(tdat.tot[,"Group.1"], tdat.tot[,"x"], col="black", lwd=3)
	
}
plot(1,1, type="n", bty="n", xaxt="n", yaxt="n", ylab="", xlab="")
text(rep(1,3),seq(1.3, 0.7, length.out=3), strsplit(paste(uT.ym, collapse=" "), " ")[[1]], col=tcols, font=2, cex=1.1)


