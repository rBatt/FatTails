png("Figures_v2/Xi_InfE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(final[gg_sh_logic&final[,"Type"]!="Phys","InfE"], final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~InfE, data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sigma[infinity]/sigma[E]), side=1, line=1.5)
dev.off()



png("Figures_v2/Xi_InfE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","InfE"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(InfE), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity]/sigma[E])), side=1, line=1.5)
dev.off()


png("Figures_v2/Xi_sigE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(final[gg_sh_logic&final[,"Type"]!="Phys","sigE"], final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~sigE, data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sigma[E]), side=1, line=1.5)
dev.off()

png("Figures_v2/Xi_sigE.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","sigE"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(sigE), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[E])), side=1, line=1.5)
dev.off()

png("Figures_v2/Xi_sigInf.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","sigInf"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(sigInf), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[infinity])), side=1, line=1.5)
dev.off()


png("Figures_v2/Xi_eInf.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(log10(final[gg_sh_logic&final[,"Type"]!="Phys","Einf"]), final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~log10(Einf), data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(log[10](sigma[E]/sigma[infinity])), side=1, line=1.5)
dev.off()

png("Figures_v2/Xi_eInf.png", width=3.4, height=3.4, units="in", res=300, bg=myWhite)
dscat(final[gg_sh_logic&final[,"Type"]!="Phys","Einf"], final[gg_sh_logic&final[,"Type"]!="Phys","sh_0"], mar=c(2.5,2.5,0,0), cex=1, ps=9, family="Times", mgp=c(3, 0.5, 0), tcl=-0.35)
abline(lm(sh_0~Einf, data=final[gg_sh_logic&final[,"Type"]!="Phys",]), lty="dashed", lwd=2)
mtext(bquote(xi), side=2, line=1.35, cex=1)
mtext(bquote(sigma[E]/sigma[infinity]), side=1, line=1.5)
dev.off()

# 1/(1 - lambda^2)
hypAR1_InfE <- sqrt(1/(1-(final[final[,"Order"]==1,"Lambda"])^2))
AR1_InfE <- final[final[,"Order"]==1, "InfE"]
png("Figures_v2/Corrected_sigma_InfE_justAR1.png", width=3.4, height=6, units="in", res=300, bg=myWhite)
par(mfrow=c(2,1), mar=c(3,2.5,0.1,0.1), cex=1, ps=9, mgp=c(1.5,0.4,0), family="Times", tcl=-0.35)
hypLab <- bquote(sqrt((1-lambda^2)^-1))
plot(hypAR1_InfE, AR1_InfE, xlab=hypLab, ylab=ieNota)

plot(final[final[,"Order"]==1, "Lambda"], final[final[,"Order"]==1, "InfE"], xlab=lNota2, ylab=ieNota)
dev.off()



