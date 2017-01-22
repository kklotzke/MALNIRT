rep.i <- 50
XG <- 1000
burnin <- 0.20
n1 <- 100
k <- 5
zeta1 <- 0
delta1 <- seq(from = -0.1, to=0.8, by=0.05)

FBF <- matrix(NA, nrow=rep.i, ncol=length(delta1))

set.seed(1)
for (ii in 1:rep.i)
{
  for(dd in 1:length(delta1))
  {
    cat(ii, ": ", dd, "\n")
    sim1 <- simdata(100, 5, c(delta1[dd], 0), zeta.offset=zeta1)
    out <- MALNRT(RT, Group = NULL, data = sim1, XG = XG, burnin = burnin, est.person = FALSE, silent = TRUE)
    FBF[ii, dd] <- out$FBF
  }
}


# Save results
setwd("~/Dropbox/Master UU Year 2/Master thesis/R/MALNIRT/sim")
save(FBF, file = "sim2.RData")

m.FBF <- matrix(NA, nrow=1, ncol=length(delta1))
for (dd in 1:length(delta1))
{
  tmp <- FBF[,dd]
  tmp <- tmp[!is.nan(tmp)]
  tmp <- tmp[!is.infinite(tmp)]
  m.FBF[1, dd] <- mean(tmp)
}
png(filename="figure2.png", width=2000, height=2000)
par(mgp = c(3, 3, 0), tcl=-1.0, mar=c(10,10,9,1) + 0.1)
par(oma=c(4,2,0,0))
par(bg=NA)
plot(1:17, m.FBF[1:17], col="orange", pch=20, cex.lab=5, cex.main = 6, cex=5, xaxt="n", bty="n", xlab="", ylab="", main="")
axis(1, at=c(1,3,5,7,9,11,13,15,17), labels=c(-0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), lwd=4, cex.axis = 4.5, cex=6, lwd.ticks=3.5, line=0)
axis(2, lwd=4, cex.axis = 4.5, lwd.ticks=3.5)
lines(1:17, m.FBF[1:17], lwd=3)
mtext(side=1, text=expression(paste("Covariance ", delta)), line=8, cex=5)
mtext(side=2, text="Log-BF", line=8, cex=5)
mtext(side=3, text=expression(paste("H0: ", delta, " = 0 vs. H1: ", delta," != 0")), line=0, cex=6)
dev.off()



