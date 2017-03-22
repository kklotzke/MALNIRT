library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim")
#load("beta_lambda.Rdata")

sim.N <- 200
sim.K <- 20
sim.XG <- 1000
sim.rep <- 15
sim.delta2 <- sim.tau2 <- numeric(sim.rep)
sim.nu2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta2 <- sim.lambda2 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.cor.beta2 <- sim.cor.lambda2 <- numeric(sim.rep)
nu <- seq(from = -0.35, to = 0.2, length.out = sim.K)# #rep(-0.15,sim.K)
out.list2 <- list()
system.time({
  ii <- 1
  while (ii <= sim.rep)
  {
    print(ii)
    dat <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.15,0), tau = c(0.20,0), nu = nu)#, beta = beta, lambda = lambda)
    out <- MALNIRT(Y = Y, RT = RT, data = dat, XG = sim.XG, est.person = FALSE)

    if(!is.null(out)) {
      sim.delta2[ii] <- out$post.means[[1]]$delta
      sim.tau2[ii] <- out$post.means[[1]]$tau
      sim.nu2[ii, ] <- out$post.means[[1]]$nu
      sim.beta2[ii, ] <- out$post.means[[1]]$beta
      sim.lambda2[ii, ] <- out$post.means[[1]]$lambda
      sim.cor.beta2[ii] <- cor(dat$beta, sim.beta2[ii, ])
      sim.cor.lambda2[ii] <- cor(dat$lambda, sim.lambda2[ii, ])

      out.list2[[ii]] <- out
      save(out.list2, sim.tau2, sim.delta2, sim.nu2, sim.beta2, sim.lambda2, sim.cor.beta2, sim.cor.lambda2, file = "simMH210322_2.RData")
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

