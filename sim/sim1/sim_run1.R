library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")
#load("beta_lambda.Rdata")

sim.N <- 500
sim.K <- 10
sim.XG <- 1000
sim.rep <- 3
sim.delta1_1 <- sim.delta1_2 <- sim.tau1_1 <- sim.tau1_2 <- numeric(sim.rep)
sim.theta1_2 <- sim.zeta1_2 <- numeric(sim.rep)
sim.nu1_1 <- sim.nu1_2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta1 <- sim.lambda1 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.cor.beta1 <- sim.cor.lambda1 <- sim.cor.sig2k1_1 <- sim.cor.sig2k1_2 <- numeric(sim.rep)
#nu <- seq(from = -0.4, to = 0.2, length.out = sim.K)# #rep(-0.15,sim.K)
nu1_1 <- c(seq(0.2, -0.05, length.out = 3), seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3))
nu1_2 <- c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2))
delta1_1 <- c(0.2, 0)
delta1_2 <- c(0.3, 0)
tau1_1 <- c(0.2, 0)
tau1_2 <- c(0.3, 0)
theta1_2 <- 0.4
zeta1_2 <- 0.2

group <- c(rep(1, 500), rep(2, 500))
out.list1 <- list()
system.time({
  ii <- 1
  while (ii <= sim.rep)
  {
    print(ii)
    dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1_1, tau = tau1_1, nu = nu1_1)#, beta = beta, lambda = lambda)
    dat2 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1_2, tau = tau1_2, nu = nu1_2, beta = dat1$beta, lambda = dat1$lambda,
                         theta.offset = theta1_2, zeta.offset = zeta1_2)

    y.all <- rbind(dat1$Y, dat2$Y)
    rt.all <- rbind(dat1$RT, dat2$RT)

    out <- MALNIRT(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE)

    if(!is.null(out)) {
      sim.delta1_1[ii] <- out$post.means[[1]]$delta
      sim.tau1_1[ii] <- out$post.means[[1]]$tau
      sim.nu1_1[ii, ] <- out$post.means[[1]]$nu

      sim.delta1_2[ii] <- out$post.means[[2]]$delta
      sim.tau1_2[ii] <- out$post.means[[2]]$tau
      sim.nu1_2[ii, ] <- out$post.means[[2]]$nu
      sim.theta1_2[ii] <- out$post.means[[2]]$theta
      sim.zeta1_2[ii] <- out$post.means[[2]]$zeta


      sim.beta1[ii, ] <- out$post.means[[1]]$beta
      sim.lambda1[ii, ] <- out$post.means[[1]]$lambda
      sim.cor.beta1[ii] <- cor(dat1$beta, sim.beta1[ii, ])
      sim.cor.lambda1[ii] <- cor(dat1$lambda, sim.lambda1[ii, ])
      sim.cor.sig2k1_1[ii] <- cor(dat1$sig2k,  out$post.means[[1]]$sig2k)
      sim.cor.sig2k1_2[ii] <- cor(dat2$sig2k,  out$post.means[[2]]$sig2k)


      out.list1[[ii]] <- out
      save(out.list1, sim.tau1_1, sim.delta1_1, sim.nu1_1, sim.tau1_2, sim.delta1_2, sim.nu1_2, sim.theta1_2, sim.zeta1_2,
           sim.cor.beta1, sim.cor.lambda1, sim.cor.sig2k1_1, sim.cor.sig2k1_2,
           tau1_1, delta1_1, nu1_1, tau1_2, delta1_2, nu1_2, theta1_2, zeta1_2, file = "simulation1_250317_1.RData")
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

