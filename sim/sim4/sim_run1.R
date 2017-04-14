library(MALNIRT)

setwd("~/Desktop/github_kklotzke/MALNIRT/sim/sim4")
#load("beta_lambda.Rdata")

set.seed(2)
sim.N <- 800
sim.K <- 50
sim.XG <- 600
sim.rep <- 1
# nu1 <- matrix(c(0,0,
#                -0.02, 0.02,
#                -0.04, 0.04,
#                -0.06, 0.06,
#                -0.08, 0.08,
#                -0.10, 0.10,
#                -0.12, 0.12), nrow = 7, ncol = 2, byrow = TRUE)

nu1 <- matrix(c(-0.04, 0.04, -0.12, 0.12), nrow = 2, ncol = 2, byrow = TRUE)

delta1 <- c(0.3, 0)
tau1 <- c(0.5, 0)
sim.BIC12.nu1 <- matrix(NA, nrow = sim.rep, ncol = nrow(nu1))

system.time({
  for(r in 1:nrow(nu1))
  {
    print(r)
    nu <- seq(nu1[r, 1], nu1[r, 2], length.out = sim.K)
    #nu <- c(rep(nu1[r, 1], sim.K/2), rep(nu1[r, 2], sim.K/2))

    ii <- 1
    while (ii <= sim.rep)
    {
      print(ii)
      dat1 <- simdataLNIRT(N = sim.N, K = sim.K, delta = delta1, tau = tau1, nu = nu)

      out_unres <- MALNIRT3Steps(Y = Y, RT = RT, data = dat1, XG = sim.XG, est.person = FALSE, doBIC.nu = FALSE, doBIC.full = TRUE)
      out_res <- MALNIRT3Steps(Y = Y, RT = RT, data = dat1, XG = sim.XG, est.person = FALSE, doBIC.nu = TRUE, doBIC.full = FALSE)

      if(!is.null(out_unres) && !is.null(out_res)) {
        BIC.nu.unres <- out_unres$BIC.full[[1]]
        BIC.nu.res <- out_res$BIC.nu[[1]]

        sim.BIC12.nu1[ii, r] <- BIC.nu.unres$BIC - BIC.nu.res$BIC

        save(sim.BIC12.nu1, tau1, delta1, nu1, file = "simulation4_060417_1.RData")
        ii <- ii + 1
      }
    }
  }
})

print(sim.BIC12.nu1)
