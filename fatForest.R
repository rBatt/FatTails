


load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")


library(party)
library(rpart)

library(randomForest)

# ============================================
# = Prepare a fish data frame for the forest =
# ============================================
# Grab fish GEV results and taxonomic info
not4tree <- c("Type","taxLvl","Community","Species","a1","a2","a3","b1","b2","b3","nll","AICc","Period","taxID", "mu_0","sig_0", "Level2_time","Level2_normTime","logsd","Level2_logNormTime","Duration","logSd","level", "shape.sig", "mean", "sd", "logMean")
bigFish0 <- merge(fish.gev[,c("taxID","location","Order","Family","Genus")], data.2, all.x=TRUE)
bigFish0 <- bigFish0[,!names(bigFish0)%in%not4tree]
bigFish0 <- bigFish0[complete.cases(bigFish0),]
row.names(bigFish0) <- NULL
bigFish0[,"Order"] <- factor(bigFish0[,"Order"])
bigFish0[,"Family"] <- factor(bigFish0[,"Family"])
bigFish0[,"Genus"] <- factor(bigFish0[,"Genus"])
bigFish0[,"P"] <- factor(bigFish0[,"P"])
bigFish0[,"Q"] <- factor(bigFish0[,"Q"])
bigFish0[,"variable"] <- factor(bigFish0[,"variable"])
bigFish0[,"location"] <- factor(bigFish0[,"location"])

# Prepare a fish data frame w/ other types of variables (really big)
# Grab the GEV results from all of the other variable types, organize into a per-lake wide format
# data.2.wide0 <- data.2[!data.2[,"Type"]%in%c("Meteorological","Biological"),]
# data.2.wide <- reshape(data.2.wide0[,c("location","variable", "mu_0","sig_0","sh_0")], timevar=c("variable"), idvar=c("location"), direction="wide")

# bigFish <- merge(bigFish0, data.2.wide, all.x=TRUE)
# colsNA <- ldply(bigFish, function(x)any(is.na(x)))[,2]
# bigFish <- bigFish[,!colsNA]

# ===============
# = Fish Forest =
# ===============
fish.prox0 <- randomForest(sh_0~., bigFish0, importance=TRUE)
varImpPlot(fish.prox0)
f.blah <- importance(fish.prox0, type=1)
f.blah <- f.blah[order(f.blah[,1], decreasing=TRUE),]
f.blah[1]/f.blah[2]



fish.party <- cforest(sh_0~., bigFish0)
fish.party.imp <- sort(varimp(fish.party), decreasing=TRUE)

fish.pp.N <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="N", plot=FALSE)

fish.pp.sesh0 <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="se.sh_0", plot=FALSE)
fish.pp.loc0 <- partialPlot(fish.prox0, pred.data=bigFish0, x.var="location", plot=FALSE)
fish.pp.loc.or <- order(fish.pp.loc0$y)
fish.pp.loc <- fish.pp.loc0
fish.pp.loc$y <- fish.pp.loc0$y[fish.pp.loc.or]
fish.pp.loc$x <- factor(fish.pp.loc0$x[fish.pp.loc.or], levels=fish.pp.loc0$x[fish.pp.loc.or])

dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.5,0.3,0), tcl=-0.25)
plot(fish.pp.sesh0, , type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Fish~~xi))
plot(fish.pp.loc$y, xlab=bquote(Lake~~name), xaxt="n", pch=19, ylab=bquote(Fish~~xi))
axis(side=1, at=1:length(fish.pp.loc$x), labels=FALSE)
text(1:length(fish.pp.loc$x), par("usr")[3]-0.0025, labels=fish.pp.loc$x, pos=1, xpd=TRUE, srt=45)


# fish.party.imp <- sort(varimp(fish.party, conditional=TRUE), decreasing=TRUE)
# fish.imp.party <- varimp(fish.party, conditional=TRUE)

# fish.party2 <- ctree(sh_0~., bigFish0)

# 
# 
# pt <- party:::prettytree(fish.party@ensemble[[5]], names(fish.party@data@get("input"))) 
# pt 
# nt <- new("BinaryTree") 
# nt@tree <- pt 
# nt@data <- fish.party@data 
# nt@responses <- fish.party@responses 
# nt 
# plot(nt, type="simple")

# fish.prox <- randomForest(sh_0~., bigFish, importance=TRUE, ntree=1E4)
# varImpPlot(fish.prox)




# ============================================
# = Prepare a zoop data frame for the forest =
# ============================================
# Grab zoop GEV results and taxonomic info
bigZoop00 <- merge(zoop.gev[,c("taxID","location","variable","Phylum","Class","Order","Family")], data.2, all.x=TRUE)
bigZoop00 <- bigZoop00[,!names(bigZoop00)%in%not4tree]
bigZoop00 <- bigZoop00[complete.cases(bigZoop00),]
row.names(bigZoop00) <- NULL

# ldply(bigZoop00, function(x)class(x))
bigZoop00[,"Phylum"] <- factor(bigZoop00[,"Phylum"])
bigZoop00[,"Class"] <- factor(bigZoop00[,"Class"])
bigZoop00[,"Order"] <- factor(bigZoop00[,"Order"])
bigZoop00[,"Family"] <- factor(bigZoop00[,"Family"])
# bigZoop00[,"Genus"] <- factor(bigZoop00[,"Genus"]) # 33 levels, doesn't work with randomForest
bigZoop00[,"P"] <- factor(bigZoop00[,"P"])
bigZoop00[,"Q"] <- factor(bigZoop00[,"Q"])
bigZoop00[,"variable"] <- factor(bigZoop00[,"variable"])


# expZoopGen <- model.matrix(~Genus, data=bigZoop00)[,-1]
# bigZoop001 <- bigZoop00[,names(bigZoop00)[!names(bigZoop00)%in%"Genus"]]
# bigZoop0 <- cbind(bigZoop001, expZoopGen)

# Prepare a zoop data frame w/ other types of variables (really big)
# Grab the GEV results from all of the other variable types, organize into a per-lake wide format
# zoop.wide0 <- data.2[!data.2[,"Type"]%in%c("Meteorological","Biological"),]
# zoop.wide <- reshape(zoop.wide0[,c("location","variable", "mu_0","sig_0","sh_0")], timevar=c("variable"), idvar=c("location"), direction="wide")

# bigZoop <- merge(bigZoop0, zoop.wide, all.x=TRUE)
# colsNA <- ldply(bigZoop, function(x)any(is.na(x)))[,2]
# bigZoop <- bigZoop[,!colsNA]


# ===============
# = Zoop Forest =
# ===============
zoop.prox0 <- randomForest(sh_0~., bigZoop00, importance=TRUE, ntree=1E3)
varImpPlot(zoop.prox0)
z.blah <- importance(zoop.prox0, type=1)
z.blah[order(z.blah[,1]),]

partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="N")

zoop.party <- cforest(sh_0~., bigZoop00)
z.party.imp <- sort(varimp(zoop.party), decreasing=TRUE)
# zoop.prox <- randomForest(sh_0~., bigZoop, importance=TRUE, ntree=1E3)
# varImpPlot(zoop.prox)

zoop.pp.sesh0 <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="se.sh_0")
zoop.pp.N <- partialPlot(zoop.prox0, pred.data=bigZoop00, x.var="N")

dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35)
plot(zoop.pp.sesh0, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Zooplankton~~xi))
plot(zoop.pp.N, type="l", xlab=bquote(Sample~~size),ylab=bquote(Zooplankton~~xi))

# 
# # ===========================
# # = quick combine fish zoop =
# # ===========================
# fz00 <- rbind(bigFish0, bigZoop00[,names(bigZoop00)[!names(bigZoop00)%in%c("Phylum","Class")]])
# zfTax <- model.matrix(~Genus+Order+Family, data=fz00)[,-1]
# fz0 <- fz00[,names(fz00)[!names(fz00)%in%c("Genus","Family","Order")]]
# fz <- cbind(fz0, zfTax)
# 
# 
# 
# fz.prox <- randomForest(sh_0~.,fz, importance=TRUE)
# varImpPlot(fz.prox)
# partialPlot(fz.prox, pred.data=fz, x.var="N")
# partialPlot(fz.prox, pred.data=fz, x.var="se.sh_0")
# partialPlot(fz.prox, pred.data=fz, x.var="se.sh_0")
# 
# fz.party <- cforest(sh_0~., data=fz00, controls=cforest_unbiased(ntree=1E3))
# sort(varimp(fz.party))

# ====================
# = Try all fat data =
# ====================
not4tree2 <- c("taxLvl","Species","a1","a2","a3","b1","b2","b3","nll","AICc","Period", "taxID", "mu_0","sig_0", "Level2_time","Level2_normTime","logsd","Level2_logNormTime","Duration","logSd","level", "shape.sig","mean", "logMean","sd", "sigEps", "residual_mu_0", "variable", "Genus")
data.tree00 <- data.2[,names(data.2)[!names(data.2)%in%c(not4tree2)]]
data.tree00 <- data.tree00[complete.cases(data.tree00),]
# dtVar <- model.matrix(~variable, data=data.tree00)
# data.tree0 <- cbind(data.tree00[,names(data.tree00)[!names(data.tree00)%in%c("variable")]], dtVar[,-1])

dt.prox <- randomForest(sh_0~., data.tree00, importance=TRUE)
varImpPlot(dt.prox)
dt.imp <- importance(dt.prox, type=1)
dt.imp[order(dt.imp, decreasing=TRUE),]

dt.party <- cforest(sh_0~., data=data.tree00)
dt.party.imp <- varimp(dt.party)
sort(dt.party.imp, decreasing=TRUE)



data.pp.sesh <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("se.sh_0"), plot=FALSE)
data.pp.type0 <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("Type"), plot=FALSE)
data.pp.ressh <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("residual_sh_0"), plot=FALSE)
data.pp.ressig <- partialPlot(dt.prox, pred.data=data.tree00, x.var=c("residual_sig_0"), plot=FALSE)

fish.pp.loc.or <- order(data.pp.type0$y)
data.pp.type <- data.pp.type0
data.pp.type$y <- data.pp.type0$y[fish.pp.loc.or]
data.pp.type$x <- factor(data.pp.type0$x[fish.pp.loc.or], levels=data.pp.type0$x[fish.pp.loc.or])

dev.new(width=6, height=6)
par(mfrow=c(2,2), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35)
plot(data.pp.sesh, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(xi))
plot(data.pp.type$y, xlab=bquote(Variable~~type), xaxt="n", pch=19, ylab=bquote(xi))
axis(side=1, at=1:length(data.pp.type$x), labels=FALSE)
text(1:length(data.pp.type$x), par("usr")[3]-0.0025, labels=substr(data.pp.type$x, 1, 5), pos=1, xpd=TRUE, srt=45)
plot(data.pp.ressh, type="l", xlab=bquote(ARMA~~residual~~xi),ylab=bquote(xi))
plot(data.pp.ressig, type="l", xlab=bquote(ARMA~~residual~~sigma),ylab=bquote(xi))


dt.party.imp.cond <- varimp(dt.party, conditional=TRUE)




