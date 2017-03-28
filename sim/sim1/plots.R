setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")

load(file = "simulation1_250328_1.RData")
load(file = "simulation1_250328_2.RData")
load(file = "simulation1_250328_3.RData")
load(file = "simulation1_250328_4.RData")

out.nu1 <- rbind(sim.nu1_1, sim.nu2_1, sim.nu3_1, sim.nu4_1)
out.nu2 <- rbind(sim.nu1_2, sim.nu2_2, sim.nu3_2, sim.nu4_2)


x <- 1:10
y <- colMeans(out.nu1)
sd <- apply(out.nu1, FUN = sd, MARGIN = 2)

plot(x, nu1_1, col="blue", ylim=c(-0.6, 0.4))
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")


x <- 1:10
y <- colMeans(out.nu2)
sd <- apply(out.nu2, FUN = sd, MARGIN = 2)

plot(x, nu1_2, col="blue", ylim=c(-0.6, 0.4))
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")
