

load("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/fatForest.RData") # also contains data.2


text.corner <- function(words, corner=2, nchar.buffer.width=2, ...){
	p.usr <- par("usr")
	word.width <- nchar(words)*par("cxy")[1] #nchar(words)*par("cin")[1]
	word.height <- par("cxy")[2] #par("cin")[2]
	width.buffer <- nchar.buffer.width*par("cxy")[1] # distance in # of characters between the word and the left or right border of the plot
		
	if(par("xaxs")=="r"){
		xaxs.adj <- 0.04*diff(range(p.usr[1:2]))
	}else{
		xaxs.adj <- 0
	}
	
	if(par("yaxs")=="r"){
		yaxs.adj <- 0.04*diff(range(p.usr[3:4]))
	}else{
		yaxs.adj <- 0
	}

	
	left <- p.usr[1] + (word.width/2) - xaxs.adj
	top <- p.usr[4] - (word.height/2 + word.height/2) + yaxs.adj
	right <- p.usr[2] - (word.width/2) - xaxs.adj
	bottom <- p.usr[3] + (word.height/2 + word.height/2) + yaxs.adj
	
	print(left)
	print(top)
	
	switch(corner,
		text(x=left, y=bottom, words, ...),
		text(x=left, y=top, words, ...),
		text(x=right, y=top, words, ...),
		text(x=right, y=top, words, ...)
		)
	
}


# ====================================
# = Fat Forest: Zoop Marginal Effect =
# ====================================
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/zoop_Marginal.png", res=150, units="in", height=5, width=3.5)
# dev.new(width=3.5, height=5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35)
plot(zoop.pp.sesh0, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Zooplankton~~xi))
text.corner("A")
# plot(zoop.pp.N, type="l", xlab=bquote(Sample~~size),ylab=bquote(Zooplankton~~xi))
plot(zoop.pp.xires, type="l", xlab=bquote(xi[Residual]),ylab=bquote(Zooplankton~~xi))
text.corner("B")
dev.off()


# ====================================
# = Fat Forest: Fish Marginal Effect =
# ====================================
# dev.new(width=3.5, height=5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/fish_Marginal.png", res=150, units="in", height=5, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.5,0.3,0), tcl=-0.25)
plot(fish.pp.xires, type="l", xlab=bquote(xi[Residual]),ylab=bquote(Fish~~xi))
plot(fish.pp.sesh0, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(Fish~~xi))
# plot(fish.pp.loc$y, xlab=bquote(Lake~~name), xaxt="n", pch=19, ylab=bquote(Fish~~xi))
# axis(side=1, at=1:length(fish.pp.loc$x), labels=FALSE)
# text(1:length(fish.pp.loc$x), par("usr")[3]-0.0025, labels=fish.pp.loc$x, pos=1, xpd=TRUE, srt=45)
dev.off()


# =========================================
# = Fat Forest: Full Data Marginal Effect =
# =========================================
# dev.new(width=3.5, height=5)
png("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Figures/Supplement/fullData_Marginal.png", res=150, units="in", height=5, width=3.5)
par(mfrow=c(2,1), mar=c(2.5,2.5,0.5,0.5), ps=10, family="Times", mgp=c(1.25,0.4,0), tcl=-0.35, cex=1)
plot(data.pp.xires, type="l", xlab=bquote(xi[Residual]),ylab=bquote(xi))
plot(data.pp.sesh, type="l", xlab=bquote(Standard~~error~~of~~xi),ylab=bquote(xi))
# plot(data.pp.type$y, xlab=bquote(Variable~~type), xaxt="n", pch=19, ylab=bquote(xi))
# axis(side=1, at=1:length(data.pp.type$x), labels=FALSE)
lakeShort <- c("Physical"="Phys","Biological"="Bio", "Meteorological"="Met", "Chemical"="Chem")
# text(1:length(data.pp.type$x), par("usr")[3]-0.0002, labels=lakeShort[as.character(data.pp.type$x)], pos=1, xpd=TRUE, srt=45)
dev.off()



