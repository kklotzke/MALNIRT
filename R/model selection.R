#' @export
computeBIC.zeta <- function(RTg, lambda, zetag, sig2gk, deltag)
{
  G <- length(RTg)
  K <- ncol(RTg[[1]])
  I <- matrix(1, nrow = K, ncol = K)

  if(G <= 1) {
    print("BIC for group means requires at least 2 groups.")
    return(NULL)
  }

  zeta.mean <- mean(zetag)
  lik <- matrix(NA, nrow = G, ncol = 2) # Groups x Models
  Ng <- numeric(G)

  for(g in 1:G) {
    Ng[g] <- nrow(RTg[[g]])

    # Model 1: G speed group parameters
    mu1 <- lambda - zetag[[g]]
    sigma <- I * deltag[g] + diag(sig2gk[[g]])
    lik[g, 1] <- sum(mvtnorm::dmvnorm(RTg[[g]], mean = mu1, sigma = sigma, log = TRUE))

    # Model 2: one speed group parameter
    mu2 <- lambda - zeta.mean
    lik[g, 2] <- sum(mvtnorm::dmvnorm(RTg[[g]], mean = mu2, sigma = sigma, log = TRUE))
  }

  # Parameters model 1
  # G Group means + K item difficulties + G*K measurement error variances + G covariance parameters
  d1 <- G + K + G*K + G

  # Parameters model 2
  # 1 Group mean + K item difficulties + G*K measurement error variances + G covariance parameters
  d2 <- 1 + K + G*K + G

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
