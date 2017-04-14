setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim3")

load(file = "simulation3_050417_1.RData")
load(file = "simulation3_050417_2.RData")
#load(file = "simulation1_020417_3.RData")
#load(file = "simulation1_020417_4.RData")

out.BIC12.theta <- cbind(sim.BIC12.theta1, sim.BIC12.theta2)
out.BIC12.zeta <- cbind(sim.BIC12.zeta1, sim.BIC12.zeta2)

rownames(out.BIC12.theta) <- 1:50
colnames(out.BIC12.theta) <- c("0", ".04", ".08", ".12", ".16", ".20", ".24")
rownames(out.BIC12.zeta) <- 1:50
colnames(out.BIC12.zeta) <- c("0", ".04", ".08", ".12", ".16", ".20", ".24")

cat("BIC full model - restricted model \n\n")
cat("theta \n")
print(round(colMeans(out.BIC12.theta), digits=2))
cat("\n zeta \n")
print(round(colMeans(out.BIC12.zeta), digits=2))

cat("\n SD theta \n")
print(summary(apply(out.BIC12.theta, FUN = sd, MARGIN = 2)))

cat("\n SD zeta \n")
print(summary(apply(out.BIC12.zeta, FUN = sd, MARGIN = 2)))
