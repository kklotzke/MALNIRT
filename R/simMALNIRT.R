#library(MASS)
#library(mvtnorm)

#' @importFrom mvtnorm rmvnorm
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
  for(ii in 1:N) {
    RT[ii,]	<- mvtnorm::rmvnorm(1, mean=lambda - zeta[ii], sigma=Sigma) #nr replications
    #RT[ii,]	<- mvrnorm(1,mu=as.vector(lambda - zeta[ii]),Sigma=Sigma, empirical=FALSE) #nr replications
  }
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
simdataLNIRT <- function(N, K, delta, tau, nu, lambda, beta, zeta.offset = 0, theta.offset = 0){
#data.lnirt <- simdataLNIRT(N = 10, K = 5, delta = c(0.1,0), tau = c(0.2,0), nu = rep(-0.25, 5))

  # Covariance matrix RTs
  sig2k.lnrt <- runif(K, 0.5, 1.5)
  Sigma.lnrt  <- diag(sig2k.lnrt) + delta[1]

  # Covariance matrix Responses
  sig2k.irt <- runif(K, 1, 1)
  Sigma.irt  <- diag(sig2k.irt) + tau[1]

  # Covariance matrix SAT
  Sigma.nu <- diag(nu)

  # Covariance matrix RTs + Respones + SAT
  upper <- cbind(Sigma.lnrt, Sigma.nu)
  lower <- cbind(Sigma.nu, Sigma.irt)
  Sigma <-rbind(upper, lower)

  # LNRT
  if(missing(lambda)) {
    lambda <- rnorm(K, mean=5, sd=1)
  }
  if (delta[2] == 0)
    sd = 0
  else
    sd = sqrt(delta[2])
  zeta <- rnorm(N, mean = 0, sd = sd)
  zeta <- zeta - mean(zeta) + zeta.offset

  # IRT
  if(missing(beta)) {
    beta <- rnorm(K, mean = 0, sd = 1)
  }
  if (tau[2] == 0)
    sd = 0
  else
    sd = sqrt(tau[2])
  theta <- rnorm(N, mean = 0, sd = sd)
  theta <- theta - mean(theta) + theta.offset

  # Sample data
  RTZ <- matrix(0, ncol=2*K, nrow=N) # first K columns: RTs, next K columns: Zs
  for(ii in 1:N) {
    means.lnrt <- lambda - zeta[ii]
    means.irt <- theta[ii] - beta
    means <- c(means.lnrt, means.irt)

    RTZ[ii,]	<- mvtnorm::rmvnorm(1, mean=means, sigma=Sigma) #nr replications
    #RT[ii,]	<- mvrnorm(1,mu=as.vector(lambda - zeta[ii]),Sigma=Sigma, empirical=FALSE) #nr replications
  }

  # Truncate for binary outcome
  RTY <- RTZ
  for(kk in (K+1):(2*K)){
    RTY[RTZ[,kk]< 0,kk] <- 0
    RTY[RTZ[,kk]> 0,kk] <- 1
  }

  mul <- mean(lambda)
  sigmal <- var(lambda)
  mub <- mean(beta)
  sigmab <- var(beta)

  return(list(RTY = RTY, RTZ = RTZ, lambda = lambda, beta = beta, zeta = mean(zeta), zeta_i = zeta, theta = mean(theta), theta_i = theta, means = means,
              Sigma = Sigma, Sigma.nu = Sigma.nu, Sigma.lnrt = Sigma.lnrt, Sigma.irt = Sigma.irt, sig2k.lnrt = sig2k.lnrt, sig2k.irt = sig2k.irt, delta = delta, tau = tau, nu = nu))
}
