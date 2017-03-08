sim.N <- 500
sim.K <- 20
sim.XG <- 2000
sim.rep <- 20
sim.delta <- sim.tau <- numeric(sim.rep)
sim.nu <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta <- sim.lambda <- matrix(NA, nrow = sim.rep, ncol = sim.K)
#dat.tmp <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.25,sim.K))
#be <- dat.tmp$beta
#la <- dat.tmp$lambda

for (ii in 1:sim.rep)
{
  print(ii)
  dat <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.25,sim.K))
  out <- MALNIRT.1StepZT2(Y = Y, RT = RT, data = dat, XG = sim.XG, est.person = FALSE)

  sim.delta[ii] <- out$delta[1]
  sim.tau[ii] <- out$tau[1]
  sim.nu[ii, ] <- out$nu
  sim.beta[ii, ] <- out$beta
  sim.lambda[ii, ] <- out$lambda
}

nu <- rep(-0.25,sim.K)
setwd("~/Desktop/github_kklotzke/MALNIRT/sim")
save(sim.tau, sim.delta, sim.nu, sim.beta, sim.lambda, nu = nu, file = "simMH080317.RData")

plot(sim.tau)
plot(sim.delta)
hist(out$data.chain1$tau.1[200:2000], breaks = 50)
hist(out$data.chain1$delta.1[200:2000], breaks = 50)
hist(out$data.chain1$nu[200:2000], breaks = 50)
mean(sim.tau)
mean(sim.delta)
sd(sim.tau)
sd(sim.delta)

plot(colMeans(sim.nu) - (-0.25))
plot(rowMeans(sim.nu) - (-0.25))
summary(colMeans(sim.nu)- (-0.25))
apply(sim.nu, FUN = sd, MARGIN = 2)

plot(colMeans(sim.beta[-bad, ]) - be)
summary(colMeans(sim.beta[-bad, ]) - be)
apply(sim.beta[-bad, ], FUN = sd, MARGIN = 2)

plot(colMeans(sim.lambda) - la)
summary(colMeans(sim.lambda) - la)
apply(sim.lambda, FUN = sd, MARGIN = 2)


plot(1:2000, out$data.chain1$nu, type = "l", col = "red", main = "nu_1",
     xlab = "", ylab = "", xaxt="n", frame.plot=F, cex.axis=1.1)
lines(1:2000, out$data.chain2$nu, col = "blue")

plot(1:2000, out$data.chain1$tau.1, type = "l", col = "red", main = "tau",
     xlab = "", ylab = "", xaxt="n", frame.plot=F, cex.axis=1.1)
lines(1:2000, out$data.chain2$tau.1, col = "blue")

plot(1:2000, out$data.chain1$delta.1, type = "l", col = "red", main = "delta",
     xlab = "", ylab = "", xaxt="n", frame.plot=F, cex.axis=1.1)
lines(1:2000, out$data.chain2$delta.1, col = "blue")

