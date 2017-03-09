sim.N <- 250
sim.K <- 20
sim.XG <- 1000
sim.rep <- 10
sim.delta <- sim.tau <- numeric(sim.rep)
sim.nu <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta <- sim.lambda <- matrix(NA, nrow = sim.rep, ncol = sim.K)
#dat.tmp <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.25,sim.K))
#be <- dat.tmp$beta
#la <- dat.tmp$lambda

system.time({
for (ii in 1:sim.rep)
{
  print(ii)
  dat <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.25,sim.K))#, beta = be, lambda = la)
  out <- MALNIRT.1StepZT2(Y = Y, RT = RT, data = dat, XG = sim.XG, est.person = FALSE)

  sim.delta[ii] <- out$delta[1]
  sim.tau[ii] <- out$tau[1]
  sim.nu[ii, ] <- out$nu
  sim.beta[ii, ] <- out$beta
  sim.lambda[ii, ] <- out$lambda
}
})

nu <- rep(-0.25,sim.K)
setwd("~/Desktop/github_kklotzke/MALNIRT/sim")
save(sim.tau, sim.delta, sim.nu, sim.beta, sim.lambda, nu = nu, be = be, la = la, file = "simMH090317_2.RData")

hist(sim.tau, breaks=10)
plot(sim.delta)
hist(c(out$data.chain1$tau.1[200:2000], out$data.chain2$tau.1[200:2000]), breaks = 50, main = "Tau = .15", xlab = "")
hist(c(out$data.chain1$delta.1[200:2000], out$data.chain2$delta.1[200:2000]), breaks = 50, main = "Delta = .1", xlab = "")
hist(c(out$data.chain1$nu[200:2000], out$data.chain2$nu[200:2000]), breaks = 50, main = "Nu (Item 1) = -.25", xlab ="")
mean(sim.tau)
mean(sim.delta)
sd(sim.tau)
sd(sim.delta)

plot(colMeans(sim.nu) - (-0.25))
plot(rowMeans(sim.nu) - (-0.25))
summary(colMeans(sim.nu)- (-0.25))
apply(sim.nu, FUN = sd, MARGIN = 2)

plot(colMeans(sim.beta[-bad, ]) - be)
summary(colMeans(sim.beta) - be)
apply(sim.beta, FUN = sd, MARGIN = 2)

plot(colMeans(sim.lambda) - la)
summary(colMeans(sim.lambda) - la)
apply(sim.lambda, FUN = sd, MARGIN = 2)


plot(200:1000, out$data.chain1$nu.1[200:1000], type = "l", col = "red", main = "nu (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:1000, out$data.chain2$nu.1[200:1000], col = "blue")

plot(200:1000, out$data.chain1$tau.1[200:1000], type = "l", col = "red", main = "tau",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:1000, out$data.chain2$tau.1[200:1000], col = "blue")

plot(200:1000, out$data.chain1$delta.1[200:1000], type = "l", col = "red", main = "delta",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:1000, out$data.chain2$delta.1[200:1000], col = "blue")

plot(200:1000, out$data.chain1$sig2k.1[200:1000], type = "l", col = "red", main = "sig2k (item 1)",
     xlab = "", ylab = "", frame.plot=F, cex.axis=1.1)
lines(200:1000, out$data.chain2$sig2k.1[200:1000], col = "blue")

