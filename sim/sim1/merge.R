setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")

load(file = "simulation1_020417_1.RData")
load(file = "simulation1_020417_2.RData")
load(file = "simulation1_020417_3.RData")
load(file = "simulation1_020417_4.RData")

out.tau1 <- na.omit(c(sim.tau1_1, sim.tau2_1, sim.tau3_1, sim.tau4_1))
out.delta1 <- na.omit(c(sim.delta1_1, sim.delta2_1, sim.delta3_1, sim.delta4_1))
out.nu1 <- na.omit(rbind(sim.nu1_1, sim.nu2_1, sim.nu3_1, sim.nu4_1))

out.tau2 <- na.omit(c(sim.tau1_2, sim.tau2_2, sim.tau3_2, sim.tau4_2))
out.delta2 <- na.omit(c(sim.delta1_2, sim.delta2_2, sim.delta3_2, sim.delta4_2))
out.nu2 <- na.omit(rbind(sim.nu1_2, sim.nu2_2, sim.nu3_2, sim.nu4_2))

out.theta2 <- na.omit(c(sim.theta1_2, sim.theta2_2, sim.theta3_2, sim.theta4_2))
out.zeta2 <- na.omit(c(sim.zeta1_2, sim.zeta2_2, sim.zeta3_2, sim.zeta4_2))

out.cor.beta <- na.omit(c(sim.cor.beta1, sim.cor.beta2, sim.cor.beta3, sim.cor.beta4))
out.cor.lambda <- na.omit(c(sim.cor.lambda1, sim.cor.lambda2, sim.cor.lambda3, sim.cor.lambda4))
out.cor.sig2k1 <- na.omit(c(sim.cor.sig2k1_1, sim.cor.sig2k2_1, sim.cor.sig2k3_1, sim.cor.sig2k4_1))
out.cor.sig2k2 <- na.omit(c(sim.cor.sig2k1_2, sim.cor.sig2k2_2, sim.cor.sig2k3_2, sim.cor.sig2k4_2))

out.mh.accept <- (tmp.mh1 + tmp.mh2 + tmp.mh3 + tmp.mh4)/4
rownames(out.mh.accept) <- c("Group 1", "Group 2")
colnames(out.mh.accept) <- c("tau", "delta", "nu")


cat("Simulated, Estimated (SD) \n")
cat("Tau_1: ", tau1_1[1], ", ", round(mean(out.tau1), digits = 3), " (", round(sd(out.tau1), digits = 3), ")\n", sep = "")
cat("Tau_2: ", tau1_2[1], ", ", round(mean(out.tau2), digits = 3), " (", round(sd(out.tau2), digits = 3), ")\n", sep = "")
cat("Delta_1: ", delta1_1[1], ", ", round(mean(out.delta1), digits = 3), " (", round(sd(out.delta1), digits = 3), ")\n", sep = "")
cat("Delta_2: ", delta1_2[1], ", ", round(mean(out.delta2), digits = 3), " (", round(sd(out.delta2), digits = 3), ")\n", sep = "")
cat("theta_2: ", theta1_2, ", ", round(mean(out.theta2), digits = 3), " (", round(sd(out.theta2), digits = 3), ")\n", sep = "")
cat("zeta_2: ", zeta1_2, ", ", round(mean(out.zeta2), digits = 3), " (", round(sd(out.zeta2), digits = 3), ")\n", sep = "")
cat("\n")
cat("Cor. beta: ", round(mean(out.cor.beta), digits = 6), " (", round(sd(out.cor.beta), digits = 6), ")\n", sep = "")
cat("Cor. lambda: ", round(mean(out.cor.lambda), digits = 6), " (", round(sd(out.cor.lambda), digits = 6), ")\n", sep = "")
cat("Cor. sig2k_1: ", round(mean(out.cor.sig2k1), digits = 6), " (", round(sd(out.cor.sig2k1), digits = 6), ")\n", sep = "")
cat("Cor. sig2k_2: ", round(mean(out.cor.sig2k2), digits = 6), " (", round(sd(out.cor.sig2k2), digits = 6), ")\n\n", sep = "")


cat("Estimated nu_1 - Simulated nu_1\n")
print(summary(colMeans(out.nu1)- nu1_1, digits = 3))
cat("\n Estimated nu_1 (SD) \n")
print(summary(apply(out.nu1, FUN = sd, MARGIN = 2)))
cat("\n")
cat("Estimated nu_2 - Simulated nu_2\n")
print(summary(colMeans(out.nu2)- nu1_2))
cat("\n Estimated nu_2 (SD) \n")
print(summary(apply(out.nu2, FUN = sd, MARGIN = 2)))
cat("\n")
cat("MH acceptance rates\n")
print(round(out.mh.accept, digits=2))
cat("\n")


#hist(out.tau, breaks=5)
#hist(out.delta, breaks=5)
