library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim1")
#load("beta_lambda.Rdata")

sim.N <- 500
sim.K <- 10
sim.XG <- 600
sim.rep <- 15
sim.delta4_1 <- sim.delta4_2 <- sim.tau4_1 <- sim.tau4_2 <- numeric(sim.rep)
sim.theta4_2 <- sim.zeta4_2 <- numeric(sim.rep)
sim.nu4_1 <- sim.nu4_2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta4 <- sim.lambda4 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.cor.beta4 <- sim.cor.lambda4 <- sim.cor.sig2k4_1 <- sim.cor.sig2k4_2 <- numeric(sim.rep)
#nu <- seq(from = -0.4, to = 0.2, length.out = sim.K)# #rep(-0.15,sim.K)
nu4_1 <- c(seq(0.2, -0.05, length.out = 3), seq(-0.1, -0.2, length.out = 4), seq(-0.3, -0.4, length.out = 3))
nu4_2 <- c(seq(-0.05, 0.2, length.out = 3), c(0, -0.2, -0.4, 0.05, -0.1), seq(-0.15, 0.2, length.out = 2))
#nu4_1 <- nu4_2 <- rep(-0.1, sim.K)
delta4_1 <- c(0.2, 0)
delta4_2 <- c(0.3, 0)
tau4_1 <- c(0.2, 0)
tau4_2 <- c(0.3, 0)
theta4_2 <- -1
zeta4_2 <- 0.1

group <- c(rep(1, 500), rep(2, 500))
out.list1 <- list()
system.time({
  ii <- 1
  while (ii <= sim.rep)
  {
    print(ii)
    dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta4_1, tau = tau4_1, nu = nu4_1)#, beta = beta, lambda = lambda)
    dat2 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta4_2, tau = tau4_2, nu = nu4_2, beta = dat1$beta, lambda = dat1$lambda,
                         theta.offset = theta4_2, zeta.offset = zeta4_2)

    y.all <- rbind(dat1$Y, dat2$Y)
    rt.all <- rbind(dat1$RT, dat2$RT)

    out <- MALNIRT3Steps(Y = y.all, RT = rt.all, group = group, XG = sim.XG, est.person = FALSE)

    if(!is.null(out)) {
      sim.delta4_1[ii] <- out$post.means[[1]]$delta
      sim.tau4_1[ii] <- out$post.means[[1]]$tau
      sim.nu4_1[ii, ] <- out$post.means[[1]]$nu

      sim.delta4_2[ii] <- out$post.means[[2]]$delta
      sim.tau4_2[ii] <- out$post.means[[2]]$tau
      sim.nu4_2[ii, ] <- out$post.means[[2]]$nu
      sim.theta4_2[ii] <- out$post.means[[2]]$theta
      sim.zeta4_2[ii] <- out$post.means[[2]]$zeta


      sim.beta4[ii, ] <- (out$post.means[[1]]$beta + out$post.means[[2]]$beta)/2
      sim.lambda4[ii, ] <- (out$post.means[[1]]$lambda + out$post.means[[2]]$lambda)/2
      sim.cor.beta4[ii] <- cor(dat1$beta, sim.beta4[ii, ])
      sim.cor.lambda4[ii] <- cor(dat1$lambda, sim.lambda4[ii, ])
      sim.cor.sig2k4_1[ii] <- cor(dat1$sig2k,  out$post.means[[1]]$sig2k)
      sim.cor.sig2k4_2[ii] <- cor(dat2$sig2k,  out$post.means[[2]]$sig2k)


      out.list1[[ii]] <- out
      save(out.list1, sim.tau4_1, sim.delta4_1, sim.nu4_1, sim.tau4_2, sim.delta4_2, sim.nu4_2, sim.theta4_2, sim.zeta4_2,
           sim.cor.beta4, sim.cor.lambda4, sim.cor.sig2k4_1, sim.cor.sig2k4_2,
           tau4_1, delta4_1, nu4_1, tau4_2, delta4_2, nu4_2, theta4_2, zeta4_2, file = "simulation1_010417_4.RData")
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

