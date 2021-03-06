
library(beanplot)
load("/Users/battrd/Documents/School&Work/WiscResearch/FatTails/Data/TurnExtreme_Fat_Data.RData")
source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Scripts/Functions/fatPlot_Functions.R")


# ==============================
# = Xi / Waiting Time Beanplot =
# ==============================
bLine <- rainbow(n=5, v=0.8, s=1)
bFill <- rgb(t(col2rgb(bLine, alpha=TRUE)), alpha=40, maxColorValue=255)
beanCol <- list(c(bFill[1]),
				c(bFill[2]),
				c(bFill[3]),
				c(bFill[4])
				)

# dev.new(width=3.5, height=6)
# png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/Fig2_fatBeanXiWaiting.png", res=150, units="in", height=5.862857, width=3.42)
pdf("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/Fig1_fatBeanXiWaiting.pdf", width=3.42, height=5.862857)
par(mfrow=c(2,1), mar=c(2,2.5,0.5,0.5), ps=8, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
beanplot(sh_0~Type, data=data.fat, ylab="", yaxt="n", xaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
axis(side=2)
tbl_nms <- c("Biological","Chemical","Physical","Meteorological") # just to make sure oder is correct in label
a1_labs <- paste0(c("Bio (","Chem (","Phys (","Met ("), table(data.fat[,"Type"])[tbl_nms], ")")
axis(side=1, at=1:4, labels=a1_labs)
mtext(bquote(Tailedness~(xi)), side=2, line=1.5)
text(0.5,1.5, "A", font=2)

fin_l2Time <- is.finite(data.fat[,"Level2_time"])
beanplot(log10(Level2_time)~Type, data=data.fat[fin_l2Time,], log="", ylab="", xaxt="n", yaxt="n", border=bLine, col=beanCol, ll=0.01, beanlinewd=1.5)
wtBase <- axTicks(2)
wtLab <- parse(text=paste(10,wtBase,sep="^"))
axis(side=2, at=wtBase, labels=wtLab)
tbl_nms <- c("Biological","Chemical","Physical","Meteorological") # just to make sure oder is correct in label
a2_labs <- paste0(c("Bio (","Chem (","Phys (","Met ("), table(data.fat[fin_l2Time,"Type"])[tbl_nms], ")")
axis(side=1, at=1:4, labels=a2_labs)
mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)
text(0.5, 8.75, "B", font=2)
dev.off()





# ===========================
# = Xi Waiting Time Boxplot =
# ===========================
# png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/fatBoxXiWaiting.png", res=150, units="in", height=6, width=3.5)
# par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
# boxplot(sh_0~Type, data=data.fat, outline=FALSE, ylab="", yaxt="n", xaxt="n")
# axis(side=2)
# axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(xi~~from~~GEV), side=2, line=1.5)
# 
# boxplot(log10(Level2_time)~Type, data=data.fat, outline=FALSE, ylab="", xaxt="n", yaxt="n")
# wtBase <- pretty(log10(data.fat[,"Level2_time"]))
# wtLab <- parse(text=paste(10,wtBase,sep="^"))
# axis(side=2, at=wtBase, labels=wtLab)
# axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)
# dev.off()


# Instead of Xi, uses the 5th percentile from the norm w/ a mean of Xi (from the GEV fit) and a sd of the se from the GEV fit
# par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, cex=1, mgp=c(2, 0.4, 0), tcl=-0.3, family="Times")
# boxplot(shape.sig~Type, data=data.fat, outline=FALSE, ylab="", yaxt="n", xaxt="n")
# axis(side=2)
# axis(side=1, at=1:4, labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(xi~~from~~GEV), side=2, line=1.5)
# 
# boxplot(log10(Level2_time)~Type, data=data.fat, outline=FALSE, ylab="", xaxt="n", yaxt="n")
# wtLab <- pretty(log10(data.fat[,"Level2_time"]))
# axis(side=2, at=wtLab, labels=10^wtLab)
# axis(side=1, at=1:4,labels=c("Bio","Chem","Phys","Met"))
# mtext(bquote(Waiting~Time~(years)), side=2, line=1.5)

# embedFonts("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/MainText/Fig2_fatBeanXiWaiting.pdf") # doing this screwed up the lines