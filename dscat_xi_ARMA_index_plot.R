


source("/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/fatPlot_Functions.R")
load(file="/Users/Battrd/Documents/School&Work/WiscResearch/FatTails/Data/data2.RData")


eiNota <- bquote(sigma[E]^2~'/'~sigma[infinity]^2)
ieNota <- bquote(sigma[infinity]~'/'~sigma[E])
lNota <- bquote(ave*.~abs(~~abs(~~lambda~~phantom())~~phantom()))
lNota2 <- bquote(abs(~~abs(~~lambda~~phantom())~~phantom()))



# ===============================
# = Dscat of Xi vs. log10(InfE) =
# ===============================
dev.new(height=3.5, width=3.5)
dscat(log(data.2[,"InfE"], base=10), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(InfE), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=1, line=1.5)


# ==========================
# = Dscat of Xi vs. Lambda =
# ==========================
dev.new(height=3.5, width=3.5)
dscat(sqrt(data.2[,"Lambda"]), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(Lambda), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(.(lNota2))), side=1, line=1.5)



# =======================
# = Dscat of Xi vs SigE =
# =======================
dev.new(height=3.5, width=3.5)
dscat(sqrt(data.2[,"sigE"]), data.2[,"sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sqrt(sigE), data=data.2[,]), lty="dashed", lwd=2)
mtext(bquote(xi[Time~Series]), side=2, line=1.35, cex=1)
mtext(bquote(sqrt(sigma[E])), side=1, line=1.5)





