
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")


library(party)
library(rpart)

library(randomForest)

# ============================================
# = Prepare a fish data frame for the forest =
# ============================================
# Grab fish GEV results and taxonomic info
not4tree <- c("Type","taxLvl","Community","Species","a1","a2","a3","b1","b2","b3","nll","AICc","Period", "variable","taxID", "mu_0","sig_0", "Level2_time","Level2_normTime","logsd","Level2_logNormTime","Duration","logSd","level", "shape.sig")
bigFish0 <- merge(fish.gev, data.2, all.x=TRUE)
bigFish0 <- bigFish0[,!names(bigFish0)%in%not4tree]
bigFish0 <- bigFish0[complete.cases(bigFish0),]
row.names(bigFish0) <- NULL
bigFish0[,"Order"] <- factor(bigFish0[,"Order"])
bigFish0[,"Family"] <- factor(bigFish0[,"Family"])
bigFish0[,"Genus"] <- factor(bigFish0[,"Genus"])
bigFish0[,"P"] <- factor(bigFish0[,"P"])
bigFish0[,"Q"] <- factor(bigFish0[,"Q"])

# Prepare a fish data frame w/ other types of variables (really big)
# Grab the GEV results from all of the other variable types, organize into a per-lake wide format
data.2.wide0 <- data.2[!data.2[,"Type"]%in%c("Meteorological","Biological"),]
data.2.wide <- reshape(data.2.wide0[,c("location","variable", "mu_0","sig_0","sh_0")], timevar=c("variable"), idvar=c("location"), direction="wide")

bigFish <- merge(bigFish0, data.2.wide, all.x=TRUE)
colsNA <- ldply(bigFish, function(x)any(is.na(x)))[,2]
bigFish <- bigFish[,!colsNA]

# ===============
# = Fish Forest =
# ===============
fish.prox0 <- randomForest(sh_0~., bigFish0, importance=TRUE, ntree=1E3)
varImpPlot(fish.prox0)

fish.prox <- randomForest(sh_0~., bigFish, importance=TRUE, ntree=1E3)
varImpPlot(fish.prox)




# ============================================
# = Prepare a zoop data frame for the forest =
# ============================================
# Grab zoop GEV results and taxonomic info
bigZoop00 <- merge(zoop.gev, data.2, all.x=TRUE)
bigZoop00 <- bigZoop00[,!names(bigZoop00)%in%not4tree]
bigZoop00 <- bigZoop00[complete.cases(bigZoop00),]
row.names(bigZoop00) <- NULL

# ldply(bigZoop00, function(x)class(x))
bigZoop00[,"Phylum"] <- factor(bigZoop00[,"Phylum"])
bigZoop00[,"Class"] <- factor(bigZoop00[,"Class"])
bigZoop00[,"Order"] <- factor(bigZoop00[,"Order"])
bigZoop00[,"Family"] <- factor(bigZoop00[,"Family"])
bigZoop00[,"Genus"] <- factor(bigZoop00[,"Genus"])
bigZoop00[,"P"] <- factor(bigZoop00[,"P"])
bigZoop00[,"Q"] <- factor(bigZoop00[,"Q"])

expZoopGen <- model.matrix(~Genus, data=bigZoop00)[,-1]
bigZoop001 <- bigZoop00[,names(bigZoop00)[!names(bigZoop00)%in%"Genus"]]
bigZoop0 <- cbind(bigZoop001, expZoopGen)

# Prepare a zoop data frame w/ other types of variables (really big)
# Grab the GEV results from all of the other variable types, organize into a per-lake wide format
zoop.wide0 <- data.2[!data.2[,"Type"]%in%c("Meteorological","Biological"),]
zoop.wide <- reshape(zoop.wide0[,c("location","variable", "mu_0","sig_0","sh_0")], timevar=c("variable"), idvar=c("location"), direction="wide")

bigZoop <- merge(bigZoop0, zoop.wide, all.x=TRUE)
colsNA <- ldply(bigZoop, function(x)any(is.na(x)))[,2]
bigZoop <- bigZoop[,!colsNA]

# ===============
# = Zoop Forest =
# ===============
zoop.prox0 <- randomForest(sh_0~., bigZoop0, importance=TRUE, ntree=1E3)
varImpPlot(zoop.prox0)

zoop.prox <- randomForest(sh_0~., bigZoop, importance=TRUE, ntree=1E3)
varImpPlot(zoop.prox)



# ===========================
# = quick combine fish zoop =
# ===========================
fz00 <- rbind(bigFish0, bigZoop00[,names(bigZoop00)[!names(bigZoop00)%in%c("Phylum","Class")]])
zfTax <- model.matrix(~Genus+Order+Family, data=fz00)[,-1]
fz0 <- fz00[,names(fz00)[!names(fz00)%in%c("Genus","Family","Order")]]
fz <- cbind(fz0, zfTax)



fz.prox <- randomForest(sh_0~.,fz, importance=TRUE)
varImpPlot(fz.prox)


# ====================
# = Try all fat data =
# ====================
not4tree2 <- c("taxLvl","Species","a1","a2","a3","b1","b2","b3","nll","AICc","Period", "taxID", "mu_0","sig_0", "Level2_time","Level2_normTime","logsd","Level2_logNormTime","Duration","logSd","level", "shape.sig","mean", "logMean","sd", "sigEps")
data.tree00 <- data.2[,names(data.2)[!names(data.2)%in%c(not4tree2)]]
data.tree00 <- data.tree00[complete.cases(data.tree00),]
dtVar <- model.matrix(~variable, data=data.tree00)
data.tree0 <- cbind(data.tree00[,names(data.tree00)[!names(data.tree00)%in%c("variable")]], dtVar[,-1])

dt.prox <- randomForest(sh_0~., data.tree0, importance=TRUE)
varImpPlot(dt.prox)


