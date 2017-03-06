sim.N <- 200
sim.K <- 20
sim.XG <- 2000
sim.rep <- 20
sim.delta <- sim.tau <- numeric(sim.rep)
sim.nu <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta <- sim.lambda <- matrix(NA, nrow = sim.rep, ncol = sim.K)
dat.tmp <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.2,sim.K))
be <- dat.tmp$beta
la <- dat.tmp$lambda

for (ii in 1:sim.rep)
{
  print(ii)
  dat <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.1,0), tau = c(0.15,0), nu = rep(-0.2,sim.K), beta = be, lambda = la)
  out <- MALNIRT(Y = Y, RT = RT, data = dat, XG = sim.XG, est.person = FALSE)

  sim.delta[ii] <- out$delta[1]
  sim.tau[ii] <- out$tau[1]
  sim.nu[ii, ] <- out$nu
  sim.beta[ii, ] <- out$beta
  sim.lambda[ii, ] <- out$lambda
}

nu <- rep(-0.2,sim.K)
setwd("~/Desktop/github_kklotzke/MALNIRT/sim")
save(sim.tau, sim.delta, sim.nu, sim.beta, sim.lambda, be = beta, la = la, nu = nu, file = "simMH060317.RData")

plot(sim.tau)
plot(sim.delta)
mean(sim.tau)
mean(sim.delta)
sd(sim.tau)
sd(sim.delta)

bad <- c(5,13,17,20)
plot(colMeans(sim.nu[-bad, ]) - (-0.2))
plot(rowMeans(sim.nu[-bad, ]) - (-0.2))
summary(colMeans(sim.nu[-bad, ]) - (-0.2))
apply(sim.nu[-bad, ], FUN = sd, MARGIN = 2)

plot(colMeans(sim.beta[-bad, ]) - be)
summary(colMeans(sim.beta[-bad, ]) - be)
apply(sim.beta[-bad, ], FUN = sd, MARGIN = 2)

plot(colMeans(sim.lambda) - la)
summary(colMeans(sim.lambda) - la)
apply(sim.lambda, FUN = sd, MARGIN = 2)


