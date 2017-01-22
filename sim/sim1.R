rep.i <- 100
XG <- 1000
burnin <- 0.20
n1 <- n2 <- n3 <- 200
k <- 20
zeta1 <- 0
zeta2 <- -2
zeta3 <- 3
delta1 <- 1.0
delta2 <- 0.5
delta3 <- 0.35

# Group membership
Group <- c(rep(1,n1), rep(2,n2), rep(3,n3))

zeta <- matrix(NA, nrow = rep.i, ncol = 3)
delta <- matrix(NA, nrow = rep.i, ncol = 3)
sig2 <- matrix(NA, nrow = rep.i, ncol = 3)
sig2.sim <- matrix(NA, nrow = rep.i, ncol = 3)
cor.lambda <- matrix(NA, nrow = rep.i, ncol = 1)
cor.zetai <- matrix(NA, nrow = rep.i, ncol = 3)

set.seed(1000)
for (ii in 1:rep.i)
{
  print(ii)
  sim1 <- simdata(n1, k, c(0,delta1), zeta.offset=zeta1)
  sim2 <- simdata(n2, k, c(0,delta2), lambda=sim1$lambda, zeta.offset=zeta2)
  sim3 <- simdata(n3, k, c(0,delta3), lambda=sim1$lambda, zeta.offset=zeta3)
  df <- rbind(sim1$RT, sim2$RT, sim3$RT)
  out <- MALNRT(df, Group = Group, data = NULL, XG = XG, burnin = burnin)

  zeta[ii,] <- out$zeta
  delta[ii,] <- out$delta[2:4]
  sig2[ii,] <- out$sig2[2:4]
  sig2.sim[ii, 1] <- sim1$sig2 / k
  sig2.sim[ii, 2] <- sim2$sig2 / k
  sig2.sim[ii, 3] <- sim3$sig2 / k
  cor.lambda[ii, 1] <- cor(out$lambda, sim1$lambda)
  cor.zetai[ii, 1] <- cor(out$zeta_i[,1], sim1$zeta_i)
  cor.zetai[ii, 2] <- cor(out$zeta_i[,2], sim2$zeta_i)
  cor.zetai[ii, 3] <- cor(out$zeta_i[,3], sim3$zeta_i)
}

# Save results
setwd("~/Dropbox/Master UU Year 2/Master thesis/R/MALNIRT/sim")
save(zeta, delta, sig2, sig2.sim, cor.lambda, cor.zetai, file = "sim1.RData")

# Zeta
colMeans(zeta)
sd(zeta[,1]); sd(zeta[,2]); sd(zeta[,3])

# Delta
colMeans(delta)
sd(delta[,1]); sd(delta[,2]); sd(delta[,3])

# Sig2
colMeans(sig2)
sd(sig2[,1]); sd(sig2[,2]); sd(sig2[,3])
colMeans(sig2.sim)
sd(sig2.sim[,1]); sd(sig2.sim[,2]); sd(sig2.sim[,3])

# Correlation lambda
mean(cor.lambda)
sd(cor.lambda)

# Correlation zeta_i
colMeans(cor.zetai)
sd(cor.zetai[,1]); sd(cor.zetai[,2]); sd(cor.zetai[,3])
