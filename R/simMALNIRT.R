#library(MASS)
#library(mvtnorm)

#' @importFrom mvtnorm rmvnorm
#' @importFrom MASS mvrnorm
#' @export
simdataLNRT <- function(N, K, delta, lambda, zeta.offset = 0){

  sig2k <- runif(K, 0.5, 1.5)
  Sigma  <- diag(sig2k) + delta[1]

  if(missing(lambda)) {
    lambda <- rnorm(K, mean=5, sd=1)
    #lambda <- lambda - mean(lambda)
  }
  zeta <- rnorm(N, mean=0, sd=sqrt(delta[2]))
  zeta <- zeta - mean(zeta) + zeta.offset

  RT <- matrix(0, ncol=K, nrow=N)
  #for(ii in 1:N) {
    #RT[ii,]	<- mvtnorm::rmvnorm(1, mean=lambda - zeta[ii], sigma=Sigma) #nr replications
  #  RT[ii,]	<- mvrnorm(1,mu=as.vector(lambda - zeta[ii]),Sigma=Sigma, empirical=TRUE) #nr replications
  #}
  RT <- MASS::mvrnorm(N, mu = lambda + zeta.offset, Sigma = Sigma, empirical = TRUE)

  mul <- mean(lambda)
  sigmal <- var(lambda)

  return(list(RT = RT, lambda = lambda, zeta = mean(zeta), zeta_i = zeta, Sigma = Sigma, sig2k = sig2k, sig2 = sum(sig2k), delta = delta))
}

#' @export
simdataIRT <- function(N, K, tau, beta, theta.offset = 0){

  sig2k <- runif(K, 1, 1)
  Sigma  <- diag(sig2k) + tau[1]

  if(missing(beta)) {
    beta <- rnorm(K, mean=0, sd=1)
    #lambda <- lambda - mean(lambda)
  }
  theta <- rnorm(N, mean=0, sd=sqrt(tau[2]))
  theta <- theta - mean(theta) + theta.offset

  Z <- matrix(0, ncol=K, nrow=N)
  for(ii in 1:N) {
    Z[ii,]	<- mvtnorm::rmvnorm(1, mean=theta[ii] - beta, sigma=Sigma) #nr replications
    #RT[ii,]	<- mvrnorm(1,mu=as.vector(lambda - zeta[ii]),Sigma=Sigma, empirical=FALSE) #nr replications
  }

  Y <- matrix(0, ncol=K, nrow=N)
  for(kk in 1:K){
    Y[Z[,kk]< 0,kk] <- 0
    Y[Z[,kk]> 0,kk] <- 1
  }

  mub <- mean(beta)
  sigmab <- var(beta)

  return(list(Y = Y, Z = Z, beta = beta, theta = mean(theta), theta_i = theta, Sigma = Sigma, sig2k = sig2k, sig2 = sum(sig2k), tau = tau))
}

#' @export
simdataLNIRT <- function(N, K, delta, tau, nu, lambda, beta, zeta.offset = 0, theta.offset = 0, rep = 1){
#data.lnirt <- simdataLNIRT(N = 10, K = 5, delta = c(0,0.1), tau = c(0,0.2), nu = rep(-0.25, 5))

  # Covariance matrix RTs
  sig2k.lnrt <- runif(K, 0.5, 1.5)
  Sigma.lnrt  <- diag(sig2k.lnrt) + delta[1]

  # Covariance matrix Responses
  sig2k.irt <- runif(K, 1, 1)
  Sigma.irt  <- diag(sig2k.irt) + tau[1]

  # Covariance matrix SAT
  Sigma.nu <- diag(nu)

  # Covariance matrix RTs + Respones + SAT
  #upper <- cbind(Sigma.lnrt, Sigma.nu)
  #lower <- cbind(Sigma.nu, Sigma.irt)
  upper <- cbind(Sigma.irt, Sigma.nu)
  lower <- cbind(Sigma.nu, Sigma.lnrt)
  Sigma <-rbind(upper, lower)


  # Item parameters
  if(missing(beta) && missing(lambda)) {
    #item.param <- mvtnorm::rmvnorm(K, mean=c(0, 5), sigma=matrix(c(1,0,0,1), nrow = 2, ncol = 2))
    item.param <- MASS::mvrnorm(K, mu=c(0, 5), Sigma=matrix(c(1,0,0,1), nrow = 2, ncol = 2), empirical = TRUE)
    beta <- item.param[,1]
    lambda <- item.param[, 2]

    #beta <- rnorm(K, mean = 0, sd = 1)
    #lambda <- rnorm(K, mean=5, sd=1)
  }

  # LNRT
  if (delta[2] == 0)
    sd = 0
  else
    sd = sqrt(delta[2])
  zeta <- rnorm(N, mean = 0, sd = sd)
  zeta <- zeta - mean(zeta) + zeta.offset

  # IRT
  if (tau[2] == 0)
    sd = 0
  else
    sd = sqrt(tau[2])
  theta <- rnorm(N, mean = 0, sd = sd)
  theta <- theta - mean(theta) + theta.offset

  # Sample data
  RTZ <- matrix(0, ncol=2*K, nrow=N) # first K columns: RTs, next K columns: Zs
  for(ii in 1:N) {
    #mu <- c(theta[ii], zeta[ii])
    #theta_k <- zeta_k <- numeric(K)
    #for(kk in 1:K) {
    #  sig <- matrix(c(0, nu[kk], nu[kk], 0), nrow = 2)
    #  tz <- mvtnorm::rmvnorm(1, mean=mu, sigma=sig)
    #  theta_k <- theta[ii]#tz[1]
    #  zeta_k <- zeta[ii]#tz[2]
    #}

    #means.irt <- theta_k - beta
    #means.lnrt <- lambda - zeta_k
    means.irt <- theta[ii] - beta
    means.lnrt <- lambda - zeta[ii]
    #means <- c(means.lnrt, means.irt)
    means <- c(means.irt, means.lnrt)

    #Z_i <- mvtnorm::rmvnorm(1, mean=means.irt, sigma=Sigma.irt)
    #RT_i <- mvtnorm::rmvnorm(1, mean=means.lnrt, sigma=Sigma.lnrt)
    #RTZ[ii, ] <- cbind(Z_i, RT_i)
    #RTZ[ii,]	<- colMeans(mvtnorm::rmvnorm(1000, mean=means, sigma=Sigma)) #nr replications
    #RTZ[ii,]	<- mvtnorm::rmvnorm(1, mean=means, sigma=Sigma) #nr replications
    #RTZ[ii,]	<- MASS::mvrnorm(1, mu=as.numeric(means), Sigma=Sigma, empirical = TRUE) #nr replications
    #RT[ii,]	<- mvrnorm(1,mu=as.vector(lambda - zeta[ii]),Sigma=Sigma, empirical=FALSE) #nr replications
  }
  RTZ	<- MASS::mvrnorm(N, mu=c(theta.offset - beta, lambda - zeta.offset), Sigma=Sigma, empirical = TRUE)


  # Truncate for binary outcome
  RTY <- RTZ
  for(kk in 1:K){
    RTY[RTZ[,kk]< 0,kk] <- 0
    RTY[RTZ[,kk]> 0,kk] <- 1
  }

  mul <- mean(lambda)
  sigmal <- var(lambda)
  mub <- mean(beta)
  sigmab <- var(beta)


  # Z|T
  I <- diag(K)
  A22.inv <- solve(Sigma.lnrt)
  muz <- mean(RTZ[,1:K])
  muT <- mean(RTZ[,(K+1):(2*K)])
  Sigma.ZRT <- Sigma.irt - diag(nu) %*% A22.inv %*% diag(nu)
  Sigma.ZRT.inv  <- solve(Sigma.ZRT)

  ZRT <- matrix(NA, nrow = N, ncol = K)
  for (i in 1:N) {
    ZRT[i, ] <- muz + (nu * I) %*% A22.inv %*% (RTZ[i, (K+1):(2*K)] - muT)
  }


  return(list(Y = RTY[, 1:K], RT = RTY[, (K+1):(2*K)], Z = RTZ[, 1:K], ZRT = ZRT[, 1:K], lambda = lambda, beta = beta, zeta = mean(zeta), zeta_i = zeta, theta = mean(theta), theta_i = theta, means = means,
              Sigma = Sigma, Sigma.nu = Sigma.nu, Sigma.lnrt = Sigma.lnrt, Sigma.irt = Sigma.irt, sig2k = sig2k.lnrt, delta = delta[1], tau = tau[1], nu = nu))
}


#' @export
simdataLNIRT2 <- function(N, K, delta, tau, nu, lambda, beta, zeta.offset = 0, theta.offset = 0){
  #data.lnirt <- simdataLNIRT(N = 10, K = 5, delta = c(0,0.1), tau = c(0,0.2), nu = rep(-0.25, 5))

  I <- diag(K)
  J <- matrix(1, nrow = K, ncol = K)

  # Measurement error variances RTs
  sig2k <- runif(K, 0.5, 1.5)

  # Sigma_11, Sigma_22 and Sigma_12 in joint-model covariance matrix
  Sigma_Z <- I + tau*J # Latent responses
  Sigma_T <- sig2k*I + delta*J # RTs
  Sigma_Z_T <- nu*I # Correlation

  # LNRT
  if(missing(lambda)) {
    lambda <- rnorm(K, mean=5, sd=1)
  }
  zeta <- rnorm(N, mean = 0, sd = sqrt(delta))
  zeta <- zeta - mean(zeta) + zeta.offset

  # IRT
  if(missing(beta)) {
    beta <- rnorm(K, mean = 0, sd = 1)
  }
  theta <- rnorm(N, mean = 0, sd = sqrt(tau))
  theta <- theta - mean(theta) + theta.offset

  # Means latent responses, response times
  mu_theta <- mean(theta)
  mu_zeta <- mean(zeta)
  mu_Z <- mu_theta - beta
  mu_T <- lambda - mu_zeta

  # Sigma_11.inv, Sigma_22.inv
  Sigma_Z.inv <- I - J / (1/tau + K) #solve(Sigma_Z)
  Sigma_T.inv <- (1/sig2k)*I - ((1/sig2k)%*%t(1/sig2k)) / (1/delta + sum(1/sig2k)) #solve(Sigma_T)

  # Sigma Z|T
  Sigma_ZT <- Sigma_Z - Sigma_Z_T %*% Sigma_T.inv %*% Sigma_Z_T

  # Temporary terms
  w <- nu/sig2k
  a <- 1/delta + sum(1/sig2k)
  b <- sig2k / (sig2k - nu^2)
  c <- 1/tau + sum(b)
  A.inv <- b*I - (b %*% t(b)) / c
  d <- a + t(w) %*% A.inv %*% w
  g <- sum(nu / (sig2k - nu^2)) / c
  A.inv_w <- (nu / (sig2k - nu^2) - g*b)

  # Inverse of Sigma Z|T #
  Sigma_ZT.inv <- A.inv - (A.inv_w %*% t(A.inv_w)) / d[1,1]

  # Sample response times
  RT <- mvtnorm::rmvnorm(N, mean = mu_T, sigma = Sigma_T)

  # Conditional mean of Z|T for each item
  mu_ZT <- matrix(0, ncol = K, nrow = N)
  for (k in 1:K) {
    mu_ZT[, k] <- mu_Z[k] + nu[k] * (1 / sig2k[k]) * ((RT[, k] - mu_T[k]) - (RT[, k] - mu_T[k]) / (a * sig2k[k]))
  }

  # Sample Z|T
  ZT <- matrix(0, ncol = K, nrow = N)
  for(ii in 1:N) {
    ZT[ii, ] <- mvtnorm::rmvnorm(1, mean = mu_ZT[ii, ], sigma = Sigma_ZT)
  }

  # Conditional mean of Zk|Z.mink, T1..p for each item
  I.min1 <- diag(K-1)
  J.min1 <- matrix(1, nrow = K-1, ncol = K-1)
  ones.min1 <- rep(1, K-1)

  Y <- matrix(0, ncol = K, nrow = N)
  ZT2 <- matrix(0, ncol = K, nrow = N) # Zk|Z.mink, T1..p
  mu_ZT2 <- matrix(NA, ncol = K, nrow = N)
  var_ZT2 <- numeric(K)
  for (k in 1:K)
  {
    w.min1 <- nu[-k]/sig2k[-k]
    b.min1 <- sig2k[-k] / (sig2k[-k] - nu[-k]^2)
    c.min1 <- 1/tau + sum(b.min1)
    A.min1.inv <- b.min1*I.min1 - (b.min1 %*% t(b.min1)) / c.min1
    d.min1 <- a + t(w.min1) %*% A.min1.inv %*% w.min1
    g.min1 <- sum(nu[-k] / (sig2k[-k] - nu[-k]^2)) / c.min1
    A.min1.inv_w <- (nu[-k] / (sig2k[-k] - nu[-k]^2) - g.min1*b.min1)

    B11 <- 1 + tau - nu[k]^2 * (1/sig2k[k] - ((1/sig2k[k])^2) / a)
    B22 <- (1 - (nu[-k]^2) / sig2k[-k])*I.min1 + tau*J.min1 + w.min1 %*% t(w.min1) / a
    B12 <- tau*ones.min1 + ((nu[k] / sig2k[k]) %*% t(w.min1)) / a
    B21 <- t(B12)
    B22.inv <- A.min1.inv - (A.min1.inv_w %*% t(A.min1.inv_w)) / d.min1[1,1]

    mu_ZT2[, k] <-  mu_ZT[, k] + B12 %*% B22.inv %*% t(ZT[, -k] - mu_ZT[, -k])
    var_ZT2[k] <- B11 - B12 %*% B22.inv %*% B21

    # Sample Zk|Z.mink, T1..p
    ZT2[, k] <- rnorm(N, mean = mu_ZT2[, k], sd = sqrt(var_ZT2[k]))

    # Truncate at mean of latent item response
    #Y[ZT2[, k] > mean(mu_ZT2[, k]), k] <- 1
    Y[ZT2[, k] > 0, k] <- 1
  }

  return(list(Y = Y, RT = RT, Z = ZT2, Z2 = ZT, beta = beta, lambda = lambda,
              mu_theta = mu_theta, mu_zeta = mu_zeta, theta_i = theta, zeta_i = zeta,
              Sigma_Z = Sigma_Z, Sigma_T = Sigma_T, Sigma_Z_T = Sigma_Z_T, Sigma_Z.inv = Sigma_Z.inv,
              Sigma_T.inv = Sigma_T.inv, Sigma_ZT = Sigma_ZT, Sigma_ZT.inv = Sigma_ZT.inv,
              sig2k = sig2k, delta = delta, tau = tau, nu = nu))
}
