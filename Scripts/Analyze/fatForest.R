


load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")
load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatARMA2.z.RData")

z[,"pq"] <- z[,"p"] + z[,"q"]
data.3 <- merge(data.2, z[,c("Type","taxID","location", "variable", "pq", "lambda", "xi.resid")], all=TRUE)

library(party)
library(rpart)
library(randomForest)

not4tree0 <- c("a1", "a2", "a3", "AICc", "b1", "b2", "b3", "Community", "Duration", "level", "Level2_logNormTime", "Level2_normTime", "Level2_time", "logMean", "logSd", "mean", "mu_0", "nll", "Period", "sd", "shape.sig", "sig_0", "Species", "taxID", "taxLvl", "variable", "residual_mu_0", "residual_sig_0", "residual_sh_0", "sigEps", "sigE", "sigInf","InfE","Einf","Lambda", "P", "Q", "se.sh_0")


# ============================================
# = Prepare a zoop data frame for the forest =
# ============================================
# Grab zoop GEV results and taxonomic info
bigZoop00 <- merge(zoop.gev[,c("taxID","location","variable","Phylum","Class","Order","Family")], data.3, all.x=TRUE)
bigZoop00 <- bigZoop00[,!names(bigZoop00)%in%not4tree.fz]
bigZoop00 <- bigZoop00[complete.cases(bigZoop00),]
row.names(bigZoop00) <- NULL

# ldply(bigZoop00, function(x)class(x))
bigZoop00[,"Phylum"] <- factor(bigZoop00[,"Phylum"])
bigZoop00[,"Class"] <- factor(bigZoop00[,"Class"])
bigZoop00[,"Order"] <- factor(bigZoop00[,"Order"])
bigZoop00[,"Family"] <- factor(bigZoop00[,"Family"])
# bigZoop00[,"Genus"] <- factor(bigZoop00[,"Genus"]) # 33 levels, doesn't work with randomForest
# bigZoop00[,"P"] <- factor(bigZoop00[,"P"])
# bigZoop00[,"Q"] <- factor(bigZoop00[,"Q"])
bigZoop00[,"pq"] <- factor(bigZoop00[,"pq"])
# bigZoop00[,"variable"] <- factor(bigZoop00[,"variable"])


# ===============
# = Zoop Forest =
# ===============

# Construct Random Forests
zoop.prox0 <- randomForest(sh_0~., bigZoop00, importance=TRUE, ntree=2E3)
zoop.party <- cforest(sh_0~., bigZoop00, controls=cforest_unbiased(ntree=2E3, mtry=ceiling(sqrt(dim(bigZoop00)[2]-1))))

# Get importances
zoop.prox.imp <- importance(zoop.prox0, type=1)
zoop.prox.imp[order(zoop.prox.imp, decreasing=TRUE),] # 1st = xi.resid, 2nd = N; 3rd = location; 4th = Family
z.party.imp <- sort(varimp(zoop.party), decreasing=TRUE) # 1st = xi.resid, 2nd = N; 3rd = lambda; 4th = Class

# Get partial plots
zoop.pp.xires <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="xi.resid", plot=FALSE)
zoop.pp.N <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="N", plot=FALSE)
zoop.pp.lambda <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="lambda", plot=FALSE)
zoop.pp.class <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="Class", plot=FALSE)

# zoop.pp.loc0 <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="location", plot=FALSE)
# zoop.pp.loc.or <- order(zoop.pp.loc0$y)
# zoop.pp.loc <- zoop.pp.loc0
# zoop.pp.loc$y <- zoop.pp.loc0$y[zoop.pp.loc.or]
# zoop.pp.loc$x <- factor(zoop.pp.loc0$x[zoop.pp.loc.or], levels=zoop.pp.loc0$x[zoop.pp.loc.or])




# Summarize top 4 zoop forests
zoopForest.range <- as.data.frame(cbind(
	"predictor"=c("xi.res","N","location","class"), 
	"min"=c(
		min(zoop.pp.xires$y), 
		min(zoop.pp.N$y), 
		min(zoop.pp.lambda$y), 
		min(zoop.pp.class$y)
	),
	"max"=c(
		max(zoop.pp.xires$y), 
		max(zoop.pp.N$y), 
		max(zoop.pp.lambda$y), 
		max(zoop.pp.class$y)
	),
	"range"=c(
		diff(range(zoop.pp.xires$y)), 
		diff(range(zoop.pp.N$y)), 
		diff(range(zoop.pp.lambda$y)), 
		diff(range(zoop.pp.class$y))
	)
))




# ============================================
# = Prepare a fish data frame for the forest =
# ============================================
# Grab fish GEV results and taxonomic info
# not4tree <- c("Type","taxLvl","Community","Species","a1","a2","a3","b1","b2","b3","nll","AICc","Period","taxID", "mu_0","sig_0", "Level2_time","Level2_normTime","logsd","Level2_logNormTime","Duration","logSd","level", "shape.sig", "mean", "sd", "logMean")
not4tree.fz <- c(not4tree0, "Type") 
bigFish0 <- merge(fish.gev[,c("taxID","location","Order","Family","Genus")], data.3, all.x=TRUE)
bigFish0 <- bigFish0[,!names(bigFish0)%in%not4tree.fz]
bigFish0 <- bigFish0[complete.cases(bigFish0),]
row.names(bigFish0) <- NULL
bigFish0[,"Order"] <- factor(bigFish0[,"Order"])
bigFish0[,"Family"] <- factor(bigFish0[,"Family"])
bigFish0[,"Genus"] <- factor(bigFish0[,"Genus"])
# bigFish0[,"P"] <- factor(bigFish0[,"P"])
# bigFish0[,"Q"] <- factor(bigFish0[,"Q"])
bigFish0[,"pq"] <- factor(bigFish0[,"pq"])
# bigFish0[,"variable"] <- factor(bigFish0[,"variable"])
bigFish0[,"location"] <- factor(bigFish0[,"location"])


# ===============
# = Fish Forest =
# ===============

# Construct Random Forest
fish.prox0 <- randomForest(sh_0~., bigFish0, importance=TRUE, ntree=2E3)
fish.party <- cforest(sh_0~., bigFish0, controls=cforest_unbiased(ntree=2E3, mtry=ceiling(sqrt(dim(bigFish0)[2]-1))))



# Get importances
fish.prox.imp <- importance(fish.prox0, type=1)
fish.prox.imp[order(fish.prox.imp, decreasing=TRUE),] # 1st=xi.resid; 2nd = location; 3rd = lambda; 4th = pq
fish.party.imp <- sort(varimp(fish.party), decreasing=TRUE) # 1st = xi.resid; 2nd = lambda; 3rd = location; 4th = pq

# Partial Plots
fish.pp.xires <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="xi.resid", plot=FALSE) # 1
fish.pp.lambda <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="lambda", plot=FALSE) # 2

fish.pp.loc0 <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="location", plot=FALSE) # 3
fish.pp.loc.or <- order(fish.pp.loc0$y)
fish.pp.loc <- fish.pp.loc0
fish.pp.loc$y <- fish.pp.loc0$y[fish.pp.loc.or]
fish.pp.loc$x <- factor(fish.pp.loc0$x[fish.pp.loc.or], levels=fish.pp.loc0$x[fish.pp.loc.or])

fish.pp.pq <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="pq", plot=FALSE) # 4


# Summarize top 4 fish forests
fishForest.range <- as.data.frame(cbind(
	"predictor"=c("xi.res","lambda","location","pq"), 
	"min"=c(
		min(fish.pp.xires$y),
		min(fish.pp.lambda$y), 
		min(fish.pp.loc$y), 
		min(fish.pp.pq$y)
	),
	"max"=c(
		max(fish.pp.xires$y),
		max(fish.pp.lambda$y),
		max(fish.pp.loc$y),
		max(fish.pp.pq$y)
	),
	"range"=c(
		diff(range(fish.pp.xires$y)),
		diff(range(fish.pp.lambda$y)),
		diff(range(fish.pp.loc$y)),
		diff(range(fish.pp.pq$y))
	)
))








# ====================
# = Try all fat data =
# ====================
# Organize "all" data
not4tree.full <- c(not4tree0)
data.tree00 <- data.3[,names(data.3)[!names(data.3)%in%c(not4tree.full)]]
data.tree00 <- data.tree00[complete.cases(data.tree00),]

# Create random forests
dt.prox <- randomForest(sh_0~., data.tree00, importance=TRUE, ntree=2E3)
dt.party <- cforest(sh_0~., data=data.tree00, controls=cforest_unbiased(ntree=2E3, mtry=ceiling(sqrt(dim(data.tree00)[2]-1))))

# Determine importances
dt.imp <- importance(dt.prox, type=1)
dt.imp[order(dt.imp, decreasing=TRUE),] # 1st = xi.resid; 2nd = Type; 3rd = Lambda; 4th = N
dt.party.imp <- sort(varimp(dt.party), decreasing=TRUE) # 1st = xi.resid; 2nd = Type; 3r = lambda; 4th = location

# Create partial plots
data.pp.xires <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("xi.resid"), plot=FALSE) # 1

data.pp.type0 <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("Type"), plot=FALSE) # 2
data.pp.type.or <- order(data.pp.type0$y)
data.pp.type <- data.pp.type0
data.pp.type$y <- data.pp.type0$y[data.pp.type.or]
data.pp.type$x <- factor(data.pp.type0$x[data.pp.type.or], levels=data.pp.type0$x[data.pp.type.or])

data.pp.lambda <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("lambda"), plot=FALSE) # 3

data.pp.loc0 <- partialPlot(dt.prox, pred.data=data.tree00, x.var="location", plot=FALSE) # 4
data.pp.loc.or <- order(data.pp.loc0$y)
data.pp.loc <- data.pp.loc0
data.pp.loc$y <- data.pp.loc0$y[data.pp.loc.or]
data.pp.loc$x <- factor(data.pp.loc0$x[data.pp.loc.or], levels=data.pp.loc0$x[data.pp.loc.or])

# Summarize top 4 data forests
dataForest.range <- as.data.frame(cbind(
	"predictor"=c("xi.res","Type","lambda","location"), 
	"min"=c(
		min(data.pp.xires$y), 
		min(data.pp.type$y), 
		min(data.pp.lambda$y), 
		min(data.pp.loc$y)
	),
	"max"=c(
		max(data.pp.xires$y), 
		max(data.pp.type$y), 
		max(data.pp.lambda$y), 
		max(data.pp.loc$y)
	),
	"range"=c(
		diff(range(data.pp.xires$y)), 
		diff(range(data.pp.type$y)), 
		diff(range(data.pp.lambda$y)), 
		diff(range(data.pp.loc$y))
	)
))



save.image("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatForest.RData")


# dev.new(width=6, height=6)
# par(mfrow=c(2,2), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35, cex=1)
# plot(data.pp.sesh, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(xi))
# plot(data.pp.type$y, xlab=bquote(Variable~~type), xaxt="n", pch=19, ylab=bquote(xi))
# axis(side=1, at=1:length(data.pp.type$x), labels=FALSE)
# lakeShort <- c("Physical"="Phys","Biological"="Bio", "Meteorological"="Met", "Chemical"="Chem")
# text(1:length(data.pp.type$x), par("usr")[3]-0.0002, labels=lakeShort[as.character(data.pp.type$x)], pos=1, xpd=TRUE, srt=45)
# plot(data.pp.ressh, type="l", xlab=bquote(ARMA~~residual~~xi),ylab=bquote(xi))
# plot(data.pp.ressig, type="l", xlab=bquote(ARMA~~residual~~sigma),ylab=bquote(xi))






