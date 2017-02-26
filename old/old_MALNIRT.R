#Z <- Z.group[[1]] <- data.lnirt$RTZ[,1:K]
#Z.RT <- Z.RT.group[[1]] <- data.lnirt$ZRT


# nu <- 0.6
# I <- diag(K)
# J <- matrix(1, K, K)
# A11 <- I + tau.s0[1]
# A22 <- diag(sig2k.s0) + delta.s0[1]
# A11.inv <- I - J/(1/tau.s0[1] + K)  # solve(A11)
# a <- diag(1/sig2k.s0)
# A22.inv <- a - (a %*% J %*% a) / (1/delta.s0[1] + sum(1/sig2k.s0)) # solve(A22)
#
# mu.z <- theta.s0[1] - mean(beta.s0)
# mu.RT <- mean(lambda.s0) - zeta.s0[1]
# Sigma.zRT <- A11 - nu^2 * A22.inv
#
# for (i in 1:N) {
#   mu.zRT <- as.numeric(mu.z + nu * A22.inv %*% (RT[i, ] - mu.RT))
#   tmp <- mvtnorm::rmvnorm(1, mean = mu.zRT, sigma = Sigma.zRT)
#   lower <- upper <- numeric(K)
#   for (kk in 1:K) {
#     if (Y[i, kk] == 0) {
#       lower[kk] <- -Inf
#       upper[kk] <- 0
#     }
#     else if (Y[i, kk] == 1) {
#       lower[kk] <- 0
#       upper[kk] <- Inf
#     }
#   }
#   #browser()
#
#   Z[i, ] <- tmvtnorm::rtmvnorm(1, mean = mu.zRT, sigma = Sigma.zRT, lower = lower, upper = upper, algorithm = "gibbs")
#
#   #  Z[Y[i, kk]==0, kk] <- qnorm(runif(1, 0, pnorm(0, tmp, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y[, kk] == 0]
#   #  Z[Y[i, kk]==1, kk] <- qnorm(runif(1, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y[, kk] == 1]
# }


# See Fox et al. (2016) Equation 14 and 15
Sjc <- matrix(tau.s0[1] / (1 + (K - 1) * tau.s0[1]), ncol = 1, nrow = K - 1)
var.Z <- (1 + K * tau.s0[1])/(1 + (K - 1) * tau.s0[1])
theta1 <- matrix(mean(theta.s0), ncol = K - 1, nrow = N)

for (kk in 1:K){
  beta1 <- beta.s0[-kk]
  Z1 <- Z[, -kk] # Latent responses to all but the current item
  mu.Z <- (mean(theta.s0) - beta.s0[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
  Z[Y[, kk]==0, kk] <- qnorm(runif(N, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y[, kk] == 0]
  Z[Y[, kk]==1, kk] <- qnorm(runif(N, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y[, kk] == 1]
}

for (gg in 1:G) {
  N.g <- nrow(Y.group[[gg]])

  Sjc <- matrix(tau.s0[gg+1] / (1 + (K - 1) * tau.s0[gg+1]), ncol = 1, nrow = K - 1)
  var.Z <- (1 + K * tau.s0[gg+1])/(1 + (K - 1) * tau.s0[gg+1])
  theta1 <- matrix(theta.s0[gg], ncol = K - 1, nrow = N.g)

  for (kk in 1:K){
    beta1 <- beta.s0[-kk]
    Z1 <- Z.group[[gg]][, -kk] # Latent responses to all but the current item
    mu.Z <- (theta.s0[gg] - beta.s0[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
    Z.group[[gg]][Y.group[[gg]][, kk]==0, kk] <- qnorm(runif(N.g, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 0]
    Z.group[[gg]][Y.group[[gg]][, kk]==1, kk] <- qnorm(runif(N.g, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 1]
  }
}

# Compute Z|RT

I <- diag(K)
J <- matrix(1, K, K)
A11 <- I + tau.s0[1]
A22 <- diag(sig2k.s0) + delta.s0[1]
A11.inv <- I - J/(1/tau.s0[1] + K)  # solve(A11)
a <- diag(1/sig2k.s0)
A22.inv <- a - (a %*% J %*% a) / (1/delta.s0[1] + sum(1/sig2k.s0)) # solve(A22)

mu.z <- mean(Z)
mu.RT <- mean(RT) #mean(lambda.s0) - zeta.s0[1]

Z.RT <- matrix(NA, nrow = N, ncol = K)
for (i in 1:N)
{
  Z.RT[i, ] <- mu.z + nu.s0 * A22.inv %*% (RT[i, ] - mu.RT)
  #Z.RT[i, ] <- mu.z - 0.25 * A22.inv %*% (RT[i, ] - mu.RT)
}
Z.RT.group[[1]] <- Z.RT



### Sample covariance responses-response times ###

# Hyper parameters
# flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
var0.b <- 10^10
mu0.b <- 0

sig2k <- this.sig2k
# print(sig2k)

I <- diag(K)
J <- matrix(1, K, K)
A11 <- I + tau # data.lnirt$Sigma.irt #
A22 <- diag(sig2k) + delta #data.lnirt$Sigma.lnrt # diag(sig2k) + delta[1] # diag(data.lnirt$sig2k.lnrt) + delta
A11.inv <- I - J/(1/tau + K)  # solve(A11)
a <- diag(1/sig2k)
A22.inv <- a - (a %*% J %*% a) / (1/delta + sum(1/sig2k)) # solve(A22)
#A11.inv <- solve(A11)
#A22.inv <- solve(A22)

mu.z <- mean(Z) #theta[1] - mean(beta)#
mu.RT <- mean(RT) #mean(lambda) - zeta[1] #
#X <- A22.inv %*% (RT - mu.RT) # N x K

#nu.s0 <- 0.6
Sigma.zRT <- A11 - nu.s0^2 * A22.inv
#Sigma.zRT <- A11 - (nu.s0*I) %*% A22.inv %*% (nu.s0*I)
Sigma.zRT.inv <- solve(Sigma.zRT)

tmp1 <- tmp2 <- 0
for (i in 1:N) {
  xi <- cbind(rep(1,K), A22.inv %*% (RT[i, ] - mu.RT)) #mu.RT))
  tmp1 <- tmp1 + t(xi) %*%  Sigma.zRT.inv %*% xi
  tmp2 <- tmp2 + t(xi) %*%  Sigma.zRT.inv %*% Z.RT[i, ]
}

Sigma.b <- solve(1/var0.b + tmp1)
mu.b <- Sigma.b %*% (mu0.b/var0.b + tmp2)

# Draw nu
b <- mvtnorm::rmvnorm(1, mean=mu.b, sigma=Sigma.b)
#nu <- rnorm(1, mu.nu, sqrt(Sigma.nu))
#nu <- mean(nu)
chain[[17]][ii, 1] <- b[1, 2]
#print(b)
#cat(mu.z, "|", b[1], "\n")




# Sigma_ZT2.inv <- 1/var_ZT2
# tmp1 <- tmp2 <- 0
# b1 <- b2 <- var.b2 <- mu.b2 <- numeric(K)
# for (k in 1:K) {
#   for (i in 1:N) {
#     xi <- cbind(1, x[i, k])
#
#     tmp1 <- tmp1 + t(xi) %*%  Sigma_ZT2.inv[k] %*% xi
#     tmp2 <- tmp2 + t(xi) %*%  Sigma_ZT2.inv[k] %*% (r[i, k] + mu_Z[k])
#   }
#
#   Sigma.b2 <- solve(1/var0.beta + tmp1)
#   mu.b2 <- Sigma.b2 %*% (mu0.beta/var0.beta + tmp2)
#   if(k == 1)
#     print(mu.b2)
#
#   b.out <- mvtnorm::rmvnorm(1, mean=mu.b2, sigma=Sigma.b2)
#   b1[k] <- b.out[1]
#   b2[k] <- b.out[2]
# }
# #
# # Sigma.b <- solve(1/var0.b + tmp1)
# # mu.b <- Sigma.b %*% (mu0.b/var0.b + tmp2)
# #
# # # Draw nu
# # b <- mvtnorm::rmvnorm(1, mean=mu.b, sigma=Sigma.b)
#
# chain[[1]][ii, ] <- beta <- -b1
# chain[[17]][ii, ] <- nu <- b2


# Sigma_ZT.inv <- solve(Sigma_ZT)
#
# tmp1 <- tmp2 <- 0
# b2 <- var.b2 <- mu.b2 <- numeric(K)
#
# xi <- matrix(0, nrow = N, ncol = K)
# for (i in 1:N) {
#   xi[i, ] <- (nu*I) %*% Sigma_T.inv %*% (RT[i, ] - mu_T)
# }
#
# for (k in 1:K) {
#    for (i in 1:N) {
#      tmp1 <- tmp1 + t(rep(1, K)) %*%  Sigma_ZT.inv %*% rep(1, K)
#      tmp2 <- tmp2 + t(rep(1, K)) %*%  Sigma_ZT.inv %*% (ZT[i, ] - xi[i, ])
#   }
#
#    Sigma.b2 <- solve(1/var0.beta + tmp1)
#    mu.b2 <- Sigma.b2 %*% (mu0.beta/var0.beta + tmp2)
# #   if(k == 1)
# #     print(mu.b2)
# #
#  b2[k] <- rnorm(1, mean=mu.b2, sd=sqrt(Sigma.b2))
# #   b1[k] <- b.out[1]
# #   b2[k] <- b.out[2]
#     }
# #
# # Sigma.b <- solve(1/var0.b + tmp1)
# # mu.b <- Sigma.b %*% (mu0.b/var0.b + tmp2)
# #
# # # Draw nu
# # b <- mvtnorm::rmvnorm(1, mean=mu.b, sigma=Sigma.b)
#
#      chain[[1]][ii, ] <- beta <- -b2
# chain[[17]][ii, ] <- nu <- b2
#print(beta[1])


#print(nu)

#  var.beta <- 1/(N*(var.gen.Z[1]) + 1/var0.beta)
#  mu.beta <- var.beta*((N*(mean(theta) - colMeans(ZT2)))*(var.gen.Z[1]) + mu0.beta/var0.beta)
# # #browser()
# #var.beta <- 1/(N*(1/var_ZT2) + 1/var0.beta)
# #mu.beta <- var.beta*((N*(mean(theta) - colMeans(ZT2)))*(1/var_ZT2) + mu0.beta/var0.beta)
#
#
# # Draw K item difficulty parameters
# #chain[[1]][ii, ] <- beta <- data$beta #rnorm(K, mu.beta, sqrt(var.beta))
# chain[[1]][ii, ] <- beta <- rnorm(K, mu.beta, sqrt(var.beta))
# #browser()
#print(beta[1])

