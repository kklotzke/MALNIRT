library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim3")
#load("beta_lambda.Rdata")

set.seed(1)
sim.N <- 250
sim.K <- 10
sim.XG <- 600
sim.rep <- 1
nu1_1 <- seq(-0.2, 0.2, length.out = sim.K)
nu1_2 <- seq(-0.3, 0.3, length.out = sim.K)
delta1_1 <- c(0.3, 0)
delta1_2 <- c(0.4, 0)
tau1_1 <- c(0.3, 0)
tau1_2 <- c(0.4, 0)
theta1_2 <- c(0, 0.04, 0.08, 0.12)
zeta1_2 <- c(0, 0.04, 0.08, 0.12)
sim.BIC12.theta1 <- matrix(NA, nrow = sim.rep, ncol = length(theta1_2))
sim.BIC12.zeta1 <- matrix(NA, nrow = sim.rep, ncol = length(zeta1_2))

group <- c(rep(1, 250), rep(2, 250))
system.time({
  for(r in 1:length(theta1_2))
  {
    ii <- 1
    while (ii <= sim.rep)
    {
      print(ii)
      dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1_1, tau = tau1_1, nu = nu1_1)
      dat2 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1_2, tau = tau1_2, nu = nu1_2, beta = dat1$beta, lambda = dat1$lambda,
                           theta.offset = theta1_2[r], zeta.offset = zeta1_2[r])

      y.all <- rbind(dat1$Y, dat2$Y)
      rt.all <- rbind(dat1$RT, dat2$RT)

      out <- MALNIRT3Steps(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE, doBIC.zeta = TRUE, doBIC.theta = TRUE)

      if(!is.null(out)) {
        BIC.zeta <- out$BIC.zeta
        sim.BIC12.zeta1[ii, r] <- BIC.zeta$BIC1 - BIC.zeta$BIC2
        BIC.theta <- out$BIC.theta
        sim.BIC12.theta1[ii, r] <- BIC.theta$BIC1 - BIC.theta$BIC2

        save(sim.BIC12.theta1, sim.BIC12.zeta1, tau1_1, delta1_1, nu1_1, tau1_2, delta1_2, nu1_2, theta1_2, zeta1_2, file = "simulation3_050417_1.RData")
        ii <- ii + 1
      }
    }
  }
})

