setwd("~/Desktop/github_kklotzke/MALNIRT/sim")

out.list1 <- load(file = "simMH210322_1.RData")
out.list2 <- load(file = "simMH210322_2.RData")

out.tau <- c(sim.tau1, sim.tau2)
out.delta <- c(sim.delta1, sim.delta2)
out.nu <- rbind(sim.nu1, sim.nu2)
out.cor.beta <- c(sim.cor.beta1, sim.cor.beta2)
out.cor.lambda <- c(sim.cor.lambda1, sim.cor.lambda2)

cat("Simulated, Estimated (SD) \n")
cat("Tau: 0.2, ", round(mean(out.tau), digits = 3), " (", round(sd(out.tau), digits = 3), ")\n", sep = "")
cat("Delta: 0.15, ", round(mean(out.delta), digits = 3), " (", round(sd(out.delta), digits = 3), ")\n", sep = "")
cat("Cor. beta: ", round(mean(out.cor.beta), digits = 3), " (", round(sd(out.cor.beta), digits = 3), ")\n", sep = "")
cat("Cor. lambda: ", round(mean(out.cor.lambda), digits = 3), " (", round(sd(out.cor.lambda), digits = 3), ")\n\n", sep = "")

cat("Estimated nu - Simulated nu (Means)\n")
print(summary(colMeans(out.nu)- (seq(from = -0.35, to = 0.2, length.out = 20))))
cat("\n Estimated nu (SD) \n")
cat(round(apply(out.nu, FUN = sd, MARGIN = 2), digits = 3), "\n")

#hist(out.tau, breaks=5)
#hist(out.delta, breaks=5)
