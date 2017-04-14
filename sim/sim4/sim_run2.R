library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim3")
#load("beta_lambda.Rdata")

set.seed(2)
sim.N <- 250
sim.K <- 10
sim.XG <- 600
sim.rep <- 1
nu2_1 <- seq(-0.2, 0.2, length.out = sim.K)
nu2_2 <- seq(-0.3, 0.3, length.out = sim.K)
delta2_1 <- c(0.3, 0)
delta2_2 <- c(0.4, 0)
tau2_1 <- c(0.3, 0)
tau2_2 <- c(0.4, 0)
theta2_2 <- c(0.16, 0.20, 0.24)
zeta2_2 <- c(0.16, 0.20, 0.24)
sim.BIC12.theta2 <- matrix(NA, nrow = sim.rep, ncol = length(theta2_2))
sim.BIC12.zeta2 <- matrix(NA, nrow = sim.rep, ncol = length(zeta2_2))

group <- c(rep(1, 250), rep(2, 250))
system.time({
  for(r in 1:length(theta2_2))
  {
    ii <- 1
    while (ii <= sim.rep)
    {
      print(ii)
      dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta2_1, tau = tau2_1, nu = nu2_1)
      dat2 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta2_2, tau = tau2_2, nu = nu2_2, beta = dat1$beta, lambda = dat1$lambda,
                           theta.offset = theta2_2[r], zeta.offset = zeta2_2[r])

      y.all <- rbind(dat1$Y, dat2$Y)
      rt.all <- rbind(dat1$RT, dat2$RT)

      out <- MALNIRT3Steps(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE, doBIC.zeta = TRUE, doBIC.theta = TRUE)

      if(!is.null(out)) {
        BIC.zeta <- out$BIC.zeta
        sim.BIC12.zeta2[ii, r] <- BIC.zeta$BIC1 - BIC.zeta$BIC2
        BIC.theta <- out$BIC.theta
        sim.BIC12.theta2[ii, r] <- BIC.theta$BIC1 - BIC.theta$BIC2

        save(sim.BIC12.theta2, sim.BIC12.zeta2, tau2_1, delta2_1, nu2_1, tau2_2, delta2_2, nu2_2, theta2_2, zeta2_2, file = "simulation3_050417_2.RData")
        ii <- ii + 1
      }
    }
  }
})

