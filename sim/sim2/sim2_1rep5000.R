library(MALNIRT)
setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim2")

nu1 <- c(seq(0.2, -0.05, length.out = 3), seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3))
nu2 <- c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2))

sim.N <- 5000
sim.K <- 10
sim.XG <- 5000

delta1 <- c(0.2, 0)
delta2 <- c(0.3, 0)
tau1 <- c(0.2, 0)
tau2 <- c(0.3, 0)
theta2 <- -0.5
zeta2 <- 0.5

group <- c(rep(1, sim.N), rep(2, sim.N))
set.seed(50)
dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1, tau = tau1, nu = nu1)#, beta = beta, lambda = lambda)
dat2 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta2, tau = tau2, nu = nu2, beta = dat1$beta, lambda = dat1$lambda,
                     theta.offset = theta2, zeta.offset = zeta2)

y.all <- rbind(dat1$Y, dat2$Y)
rt.all <- rbind(dat1$RT, dat2$RT)

out <- MALNIRT(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE)
save(out, delta1, delta2, tau1, tau2, nu1, nu2, theta2, zeta2, dat1, dat2, file = "simulation2_280317_5000.Rdata")

load("simulation2_280317_5000.Rdata")
out.delta1 <- out$post.means[[1]]$delta
out.tau1 <- out$post.means[[1]]$tau
out.nu1 <- out$post.means[[1]]$nu
out.delta2 <- out$post.means[[2]]$delta
out.tau2 <- out$post.means[[2]]$tau
out.nu2 <- out$post.means[[2]]$nu
out.theta2 <- out$post.means[[2]]$theta
out.zeta2 <- out$post.means[[2]]$zeta
out.cor.beta <- cor(dat1$beta, (out$post.means[[1]]$beta + out$post.means[[2]]$beta)/2)
out.cor.lambda <- cor(dat1$lambda, (out$post.means[[1]]$lambda + out$post.means[[2]]$lambda)/2)
out.cor.sig2k1 <- cor(dat1$sig2k, out$post.means[[1]]$sig2k)
out.cor.sig2k2 <- cor(dat2$sig2k, out$post.means[[2]]$sig2k)


sink("simulation2_280317_5000.txt")
cat("Simulated, Estimated (SD) \n")
cat("Tau_1: ", tau1[1], ", ", round(mean(out.tau1), digits = 3), " (", round(out$post.means[[1]]$tau.sd, digits = 3), ")\n", sep = "")
cat("Tau_2: ", tau2[1], ", ", round(mean(out.tau2), digits = 3), " (", round(out$post.means[[2]]$tau.sd, digits = 3), ")\n", sep = "")
cat("Delta_1: ", delta1[1], ", ", round(mean(out.delta1), digits = 3), " (", round(out$post.means[[1]]$delta.sd, digits = 3), ")\n", sep = "")
cat("Delta_2: ", delta2[1], ", ", round(mean(out.delta2), digits = 3), " (", round(out$post.means[[2]]$delta.sd, digits = 3), ")\n", sep = "")
cat("theta_2: ", theta2, ", ", round(mean(out.theta2), digits = 3), " (", round(out$post.means[[2]]$theta.sd, digits = 3), ")\n", sep = "")
cat("zeta_2: ", zeta2, ", ", round(mean(out.zeta2), digits = 3), " (", round(out$post.means[[2]]$zeta.sd, digits = 3), ")\n", sep = "")
cat("\n")
cat("Cor. beta: ", round(mean(out.cor.beta), digits = 7), "\n", sep = "")
cat("Cor. lambda: ", round(mean(out.cor.lambda), digits = 7), "\n", sep = "")
cat("Cor. sig2k_1: ", round(mean(out.cor.sig2k1), digits = 7), "\n", sep = "")
cat("Cor. sig2k_2: ", round(mean(out.cor.sig2k2), digits = 7), "\n\n", sep = "")


cat("Estimated nu_1 - Simulated nu_1 (Means)\n")
print(summary(out.nu1- nu1))
cat("\n Estimated nu_1 (SD) \n")
cat(round(out$post.means[[1]]$nu.sd, digits = 3), "\n")
cat("\n")
cat("Estimated nu_2 - Simulated nu_2 (Means)\n")
print(summary(out.nu2 - nu2))
cat("\n Estimated nu_2 (SD) \n")
cat(round(out$post.means[[2]]$nu.sd, digits = 3), "\n")
sink()


#### Plots ####

## Density

par(mfrow=c(2,2))
breaks <- 50
hist(out$samples[[1]]$tau[c(a1:b1,a2:b2)], breaks=breaks, main = "tau", xlab = "")
hist(out$samples[[1]]$delta[c(a1:b1,a2:b2)], breaks=breaks, main = "delta", xlab = "")
hist(out$samples[[1]]$nu.1[c(a1:b1,a2:b2)], breaks=breaks, main = "nu (item 1)", xlab = "")
#hist(out$samples[[1]]$nu.2[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[1]]$sig2k.1[c(a1:b1,a2:b2)], breaks=breaks, main = "sig2 (item1)", xlab = "")
hist(out$samples[[1]]$sig2k.2[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[1]]$beta.1[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[1]]$lambda.1[c(a1:b1,a2:b2)], breaks=breaks)

hist(out$samples[[2]]$tau[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$delta[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$nu.1[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$nu.2[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$sig2k.1[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$sig2k.2[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$beta.1[c(a1:b1,a2:b2)], breaks=breaks)
hist(out$samples[[2]]$lambda.1[c(a1:b1,a2:b2)], breaks=breaks)

#par(mfrow=c(2,2))

## Trace plots
a1 <- 500
b1 <- 5000
a2 <- 5500
b2 <- 10000

par(mfrow=c(2,2))
plot(a1:b1, out$samples[[1]]$tau[a1:b1], type = "l", col = "red", main = "tau",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(a1:b1, out$samples[[1]]$tau[a2:b2], col = "blue")

plot(a1:b1, out$samples[[1]]$delta[a1:b1], type = "l", col = "red", main = "delta",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(a1:b1, out$samples[[1]]$delta[a2:b2], col = "blue")

plot(a1:b1, out$samples[[1]]$nu.1[a1:b1], type = "l", col = "red", main = "nu (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(a1:b1, out$samples[[1]]$nu.1[a2:b2], col = "blue")

plot(a1:b1, out$samples[[1]]$sig2k.1[a1:b1], type = "l", col = "red", main = "sig2 (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(a1:b1, out$samples[[1]]$sig2k.1[a2:b2], col = "blue")


par(mfrow=c(2,1))

## Nu
x <- 1:10
y <- out.nu1
sd <- 1.96*out$post.means[[1]]$nu.sd

plot(x, nu1, col="blue", ylim=c(-0.6, 0.4), main = "Nu (group 1) - Posterior means and 95% Credibility Intervals", xlab="", ylab="")
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")

y <- out.nu2
sd <- 1.96*out$post.means[[2]]$nu.sd

plot(x, nu2, col="blue", ylim=c(-0.6, 0.4), main = "Nu (group 2) - Posterior means and 95% Credibility Intervals", xlab="", ylab="")
points(y, col="black", pch=5)
segments(x, y-sd,x, y+sd, col="black")
epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd, col="black")
segments(x-epsilon,y+sd,x+epsilon,y+sd, col="black")




