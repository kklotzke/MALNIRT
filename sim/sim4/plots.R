setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim3")

load(file = "simulation4_050417_1.RData")
#load(file = "simulation3_050417_2.RData")
#load(file = "simulation1_010417_3.RData")
#load(file = "simulation1_010417_4.RData")

out.BIC12.nu <- sim.BIC12.nu # cbind(sim.BIC12.nu1, sim.BIC12.nu2)


x <- 1:7
y <- y1 <- colMeans(out.BIC12.nu)
sd <- apply(out.BIC12.theta, FUN = sd, MARGIN = 2)

plot(x, y, col="blue", ylim=c(10, -40))
#points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")


x <- 1:7
y <- y2 <- colMeans(out.BIC12.zeta)
sd <- apply(out.BIC12.zeta, FUN = sd, MARGIN = 2)

plot(x, y, col="blue", ylim=c(10, -40))
#points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")

plot(x, y1, col="green", type="l", xaxt = "n", ylab=expression(paste(Delta,"BIC")), xlab=expression(paste(Delta,nu))
     , main=expression(paste("Unrestricted vs. ",nu[1],"= ... =",nu[p])))

#lines(x, y2, col="red")
points(x, y1, col="green")
#points(x, y2, col="red")
abline(a = 0, b = 0, lty=3)
axis(1, at=1:7, labels=c(0, 0.04, 0.08, 0.12, 0.16, 0.20, 0.24))
#legend(x = "topright", legend = c("Theta", "Zeta"), lty = c(1,1), col=c("blue", "red"), bty = "n", cex = 1,
#       x.intersp=0.5, y.intersp = 0.2, xjust=0.5, yjust=0, inset = c(0, -0.1))
