library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim")
#load("beta_lambda.Rdata")

sim.N <- 500
sim.K <- 20
sim.XG <- 1000
sim.rep <- 15
sim.delta1 <- sim.tau1 <- numeric(sim.rep)
sim.nu1 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.beta1 <- sim.lambda1 <- matrix(NA, nrow = sim.rep, ncol = sim.K)
sim.cor.beta1 <- sim.cor.lambda1 <- numeric(sim.rep)
nu <- seq(from = -0.35, to = 0.2, length.out = sim.K)# #rep(-0.15,sim.K)
out.list1 <- list()
system.time({
  ii <- 1
  while (ii <= sim.rep)
  {
    print(ii)
    dat <- simdataLNIRT(N = sim.N, K = sim.K, delta = c(0.15,0), tau = c(0.20,0), nu = nu)#, beta = beta, lambda = lambda)
    out <- MALNIRT(Y = Y, RT = RT, data = dat, XG = sim.XG, est.person = FALSE)

    if(!is.null(out)) {
      sim.delta1[ii] <- out$post.means[[1]]$delta
      sim.tau1[ii] <- out$post.means[[1]]$tau
      sim.nu1[ii, ] <- out$post.means[[1]]$nu
      sim.beta1[ii, ] <- out$post.means[[1]]$beta
      sim.lambda1[ii, ] <- out$post.means[[1]]$lambda
      sim.cor.beta1[ii] <- cor(dat$beta, sim.beta1[ii, ])
      sim.cor.lambda1[ii] <- cor(dat$lambda, sim.lambda1[ii, ])

      out.list1[[ii]] <- out
      save(out.list1, sim.tau1, sim.delta1, sim.nu1, sim.beta1, sim.lambda1, sim.cor.beta1, sim.cor.lambda1, file = "simMH210322_1.RData")
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

