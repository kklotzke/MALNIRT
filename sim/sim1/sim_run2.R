library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")
#load("beta_lambda.Rdata")

set.seed(2)
sim.N <- 500
sim.K <- 10
sim.XG <- 1000
sim.rep <- 5
sim.delta2_1 <- sim.delta2_2 <- sim.tau2_1 <- sim.tau2_2 <- numeric(sim.rep)
sim.theta2_2 <- sim.zeta2_2 <- numeric(sim.rep)
sim.nu2_1 <- sim.nu2_2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta2 <- sim.lambda2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.cor.beta2 <- sim.cor.lambda2 <- sim.cor.sig2k2_1 <- sim.cor.sig2k2_2 <- numeric(sim.rep)
sim.mh.accept2 <- matrix(0, nrow = 2, ncol = 3)
#nu <- seq(from = -0.4, to = 0.2, length.out = sim.K)# #rep(-0.15,sim.K)
nu2_1 <- c(seq(0.2, -0.05, length.out = 3), seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3))
nu2_2 <- c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2))
#nu2_1 <- nu2_2 <- rep(-0.1, sim.K)
delta2_1 <- c(0.5, 0)
delta2_2 <- c(0.6, 0)
tau2_1 <- c(0.5, 0)
tau2_2 <- c(0.6, 0)
theta2_2 <- 1
zeta2_2 <- 0.5

group <- c(rep(1, 500), rep(2, 500))
out.list2 <- list()
system.time({
  ii <- 1
  while (ii <= sim.rep)
  {
    print(ii)
    dat1 <- simdataLNIRTP(N = sim.N, K = sim.K, delta = delta2_1[1], tau = tau2_1[1], nu = nu2_1)#, beta = beta, lambda = lambda)
    dat2 <- simdataLNIRTP(N = sim.N, K = sim.K, delta = delta2_2[1], tau = tau2_2[1], nu = nu2_2, beta = dat1$beta, lambda = dat1$lambda,
                         theta.offset = theta2_2, zeta.offset = zeta2_2)

    y.all <- rbind(dat1$Y, dat2$Y)
    rt.all <- rbind(dat1$RT, dat2$RT)

    out <- MALNIRT3Steps(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE)

    if(!is.null(out)) {
      sim.delta2_1[ii] <- out$post.means[[1]]$delta
      sim.tau2_1[ii] <- out$post.means[[1]]$tau
      sim.nu2_1[ii, ] <- out$post.means[[1]]$nu

      sim.delta2_2[ii] <- out$post.means[[2]]$delta
      sim.tau2_2[ii] <- out$post.means[[2]]$tau
      sim.nu2_2[ii, ] <- out$post.means[[2]]$nu
      sim.theta2_2[ii] <- out$post.means[[2]]$theta
      sim.zeta2_2[ii] <- out$post.means[[2]]$zeta


      sim.beta2[ii, ] <- (out$post.means[[1]]$beta + out$post.means[[2]]$beta)/2
      sim.lambda2[ii, ] <- (out$post.means[[1]]$lambda + out$post.means[[2]]$lambda)/2
      sim.cor.beta2[ii] <- cor(dat1$beta, sim.beta2[ii, ])
      sim.cor.lambda2[ii] <- cor(dat1$lambda, sim.lambda2[ii, ])
      sim.cor.sig2k2_1[ii] <- cor(dat1$sig2k,  out$post.means[[1]]$sig2k)
      sim.cor.sig2k2_2[ii] <- cor(dat2$sig2k,  out$post.means[[2]]$sig2k)

      sim.mh.accept2 <- sim.mh.accept2 + out$mh.accept

      out.list2[[ii]] <- out
      tmp.mh2 <-sim.mh.accept2/ii
      save(out.list2, sim.tau2_1, sim.delta2_1, sim.nu2_1, sim.tau2_2, sim.delta2_2, sim.nu2_2, sim.theta2_2, sim.zeta2_2,
           sim.cor.beta2, sim.cor.lambda2, sim.cor.sig2k2_1, sim.cor.sig2k2_2, tmp.mh2,
           tau2_1, delta2_1, nu2_1, tau2_2, delta2_2, nu2_2, theta2_2, zeta2_2, file = "simulation1_140417_2.RData")
      ii <- ii + 1
    }
  }
})

# print(mean(sim.tau))
# print(mean(sim.delta))
# print(sd(sim.tau))
# print(sd(sim.delta))
# print(summary(colMeans(sim.nu)- (-0.15)))
# print(apply(sim.nu, FUN = sd, MARGIN = 2))

#summary(colMeans(sim.beta) - be)
#apply(sim.beta, FUN = sd, MARGIN = 2)

#summary(colMeans(sim.lambda) - la)
#apply(sim.lambda, FUN = sd, MARGIN = 2)

