#' @export
computeBIC.full <- function(RT, lambda, zeta, Z, muZ, sig2k, tau, delta, nu)
{
  #Z <- dat4$Z
  K <- ncol(RT)
  N <- nrow(RT)
  J <- matrix(1, nrow = K, ncol = K)

  lik <- 0
  sigma <- diag(sig2k - nu^2) + J*delta + (nu %*% t(nu)) / (1 / tau + K)

  # Full model
  mu <- t((lambda - zeta) + diag(nu) %*% (diag(K) - J/(1/tau + K)) %*% t(Z - matrix(muZ, nrow = N, ncol = K, byrow = TRUE)))

  for(i in 1:N) {
    lik <- lik + sum(mvtnorm::dmvnorm(RT[i, ], mean = mu[i, ], sigma = sigma, log = TRUE))
  }

  # Parameters full model
  # 1 Group speed means + 2*K item parameters + 1 Group ability means + K measurement error variances + 2 covariance parameters + K cross-covariance parameters
  d1 <- 1 + 2*K + 1 + K + 2 + K

  # Log-likelihood full model
  lik1 <- sum(lik)

  # BIC full model
  BIC1 <- -2*lik1 + d1*log(sum(N))

  return(list(BIC = BIC1, d = d1, lik = lik1))
}

#' @export
computeBIC.nu <- function(RT, lambda, zeta, Z, muZ, sig2k, tau, delta, nu)
{
  #Z <- dat4$Z
  K <- ncol(RT)
  N <- nrow(RT)
  J <- matrix(1, nrow = K, ncol = K)

  lik <- 0
  sigma <- diag(sig2k - nu^2) + J*delta + (nu %*% t(nu)) / (1 / tau + K)

  # Full model
  mu <- t((lambda - zeta) + diag(nu) %*% (diag(K) - J/(1/tau + K)) %*% t(Z - matrix(muZ, nrow = N, ncol = K, byrow = TRUE)))

  for(i in 1:N) {
    lik <- lik + sum(mvtnorm::dmvnorm(RT[i, ], mean = mu[i, ], sigma = sigma, log = TRUE))
  }

  # Parameters model restricted nu
  # 1 Group speed means + 2*K item parameters + 1 Group ability means + K measurement error variances + 2 covariance parameters + 1 cross-covariance parameter
  d1 <- 1 + 2*K + 1 + K + 2 + 1

  # Log-likelihood full model
  lik1 <- sum(lik)

  # BIC full model
  BIC1 <- -2*lik1 + d1*log(sum(N))

  return(list(BIC = BIC1, d = d1, lik = lik1))
}


#' @export
computeBIC.zeta <- function(RTg, lambda, zetag, Zg, muZg, sig2gk, taug, deltag, nug)
{
  G <- length(RTg)
  K <- ncol(RTg[[1]])
  J <- matrix(1, nrow = K, ncol = K)

  if(G <= 1) {
    print("BIC for group means requires at least 2 groups.")
    return(NULL)
  }

  zeta.mean <- 0#mean(zetag)
  lik <- matrix(0, nrow = G, ncol = 2) # Groups x Models
  Ng <- numeric(G)

  for(g in 1:G) {
    Ng[g] <- nrow(RTg[[g]])

    sigma <- diag(sig2gk[[g]] - nug[[g]]^2) + J*deltag[g] + (nug[[g]] %*% t(nug[[g]])) / (1 / taug[g] + K)

    # Model 1: G speed group parameters
    mu1 <- t((lambda - zetag[g]) + diag(nug[[g]]) %*% (diag(K) - J/(1/taug[g] + K)) %*% t(Zg[[g]] - matrix(muZg[[g]], nrow = Ng[g], ncol = K, byrow = TRUE)))

    # Model 2: one speed group parameter
    mu2 <- t((lambda - 0) + diag(nug[[g]]) %*% (diag(K) - J/(1/taug[g] + K)) %*% t(Zg[[g]] - matrix(muZg[[g]], nrow = Ng[g], ncol = K, byrow = TRUE)))

    for(i in 1:Ng[g]) {
      lik[g, 1] <- lik[g, 1] + sum(mvtnorm::dmvnorm(RTg[[g]][i, ], mean = mu1[i, ], sigma = sigma, log = TRUE))
      lik[g, 2] <- lik[g, 2] + sum(mvtnorm::dmvnorm(RTg[[g]][i, ], mean = mu2[i, ], sigma = sigma, log = TRUE))
    }
  }

  # Parameters model 1
  # G Group speed means + 2*K item parameters + G Group ability means + G*K measurement error variances + 2*G covariance parameters + G*K cross-covariance parameters
  d1 <- G + 2*K + G + G*K + 2*G + G*K

  # Parameters model 2
  # 1 Group speed mean + 2*K item parameters + G Group ability means + G*K measurement error variances + 2*G covariance parameters + G*K cross-covariance parameters
  d2 <- 1 + 2*K + G + G*K + 2*G + G*K

  # Log-likelihood model 1
  lik1 <- sum(lik[,1])

  # Log-likelihood model 2
  lik2 <- sum(lik[,2])

  # BIC model 1
  BIC1 <- -2*lik1 + d1*log(sum(Ng))

  # BIC model 2
  BIC2 <- -2*lik2 + d2*log(sum(Ng))

  return(list(BIC1 = BIC1, BIC2 = BIC2, d1 = d1, d2 = d2, lik1 = lik1, lik2 = lik2))
}

#' @export
computeBIC.theta <- function(Zg, Zg0, beta, thetag, RTg, muTg, sig2gk, taug, deltag, nug)
{
  G <- length(RTg)
  K <- ncol(RTg[[1]])
  J <- matrix(1, nrow = K, ncol = K)

  if(G <= 1) {
    print("BIC for group means requires at least 2 groups.")
    return(NULL)
  }

  theta.mean <- 0#mean(thetag)
  lik <- matrix(0, nrow = G, ncol = 2) # Groups x Models
  Ng <- numeric(G)

  for(g in 1:G) {
    Ng[g] <- nrow(RTg[[g]])

    sigma <- ( diag(K) + J*taug[g] ) - diag(nug[[g]]) %*% solve( diag(sig2gk[[g]]) + J*deltag[g] ) %*% t(diag(nug[[g]]))
    a <- 1/deltag[g] + sum(1/sig2gk[[g]])

    mu1 <- mu2 <- matrix(NA, ncol = K, nrow = Ng[g])
    for (k in 1:K)
    {
      x <- t( (1/sig2gk[[g]][k]) * ( t(RTg[[g]][, k] - muTg[[g]][k]) - (t(1/sig2gk[[g]]) %*% (t(RTg[[g]]) - matrix(muTg[[g]], nrow = K, ncol = Ng[g]))) / a ) )

      # Model 1: G ability group parameters
      mu1[, k] <- (thetag[g] - beta[k]) + nug[[g]][k] * x

      # Model 2: one ability group parameter
      mu2[, k] <- (0 - beta[k]) + nug[[g]][k] * x
    }


    for(i in 1:Ng[g]) {
      lik[g, 1] <- lik[g, 1] + sum(mvtnorm::dmvnorm(Zg[[g]][i, ], mean = mu1[i, ], sigma = sigma, log = TRUE))
      lik[g, 2] <- lik[g, 2] + sum(mvtnorm::dmvnorm(Zg0[[g]][i, ], mean = mu2[i, ], sigma = sigma, log = TRUE))
    }
  }

  # Parameters model 1
  # G Group ability means + 2*K item parameters + G Group speed means + G*K measurement error variances + 2*G covariance parameters + G*K cross-covariance parameters
  d1 <- G + 2*K + G + G*K + 2*G + G*K

  # Parameters model 2
  # 1 Group ability mean + 2*K item parameters + G Group speed means + G*K measurement error variances + 2*G covariance parameters + G*K cross-covariance parameters
  d2 <- 1 + 2*K + G + G*K + 2*G + G*K

  # Log-likelihood model 1
  lik1 <- sum(lik[,1])

  # Log-likelihood model 2
  lik2 <- sum(lik[,2])

  # BIC model 1
  BIC1 <- -2*lik1 + d1*log(sum(Ng))

  # BIC model 2
  BIC2 <- -2*lik2 + d2*log(sum(Ng))

  return(list(BIC1 = BIC1, BIC2 = BIC2, d1 = d1, d2 = d2, lik1 = lik1, lik2 = lik2))
}
