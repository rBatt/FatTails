
# ========================================================
# = Simulated effect of sampling frequency on tailedness =
# ========================================================
library(ismev)
set.seed(42)
nYear <- 30
Freqs <- 1:365
samps <- lapply(Freqs, function(f)replicate(nYear, max(rnorm(n=f))))
shapes <- sapply(samps, function(x)gev.fit(x)$mle[3])
correl <- cor(Freqs, shapes)
png("~/Desktop/xi_sampsPerYear.png")
plot(Freqs, shapes, xlab="Obervations per year (for 30 years)", ylab=bquote(xi))
mtext(paste("r",round(correl,2),sep=" = "), line=-1, adj=0.1, font=2)
abline(lm(shapes~Freqs))
dev.off()
summary(lm(shapes~Freqs))


# ======================================================
# = Effect of sampling frequency on tailedness in data =
# ======================================================
library(data.table)
d <- fread("~/Documents/School&Work/WiscResearch/FatTails/Data/fullTimeSeries_4Tony.txt")
setkey(d, taxID, location, variable, year4, daynum)
dsum <- d[,list(Freq=as.numeric(median(table(year4)))),by=c("location","taxID","variable")]

load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data.fat.RData")
data.fat <- as.data.table(data.fat)

df_freq <- merge(data.fat[,list(Type, taxLvl, taxID, location, variable, sh_0, N)], dsum, by=c("taxID","location","variable"), all=FALSE) # drops 2 w/ non-finite shape standard error
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/FigS1_fatFreqTest.pdf", width=3.5, height=3.5)
par(mfrow=c(2,2), mar=c(2.5,2.5,1,0.75), cex=1, ps=8, tcl=-0.15, mgp=c(1,0.25,0))
df_freq[,j={plot(Freq, sh_0, main=Type, xlab="Median samples per year", ylab=bquote(xi))},by="Type"]
# abline(v=c(1, 4, 26, 52, 365))
dev.off()
df_freq[,summary(lm(sh_0~Freq))]