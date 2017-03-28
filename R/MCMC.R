#' @export
helmert <- function(p) {
  helmert <- matrix(0,ncol=p,nrow=p)

  diag(helmert)[2:p] <- -(1:(p-1))
  for(ii in 1:p){

    if(ii == 1){
      helmert[ii,1:p] <- rep(1,p)/sqrt(p)
    }
    else{
      helmert[ii,ii] <- helmert[ii,ii]/sqrt((ii-1)*ii)
      helmert[ii,1:(ii-1)] <- 1/sqrt((ii-1)*ii)
    }
  }
  return(helmert)
}

#' Sample person ability parameter
#' @export
sampleTheta_i <- function(Z.g, beta, tau.g, theta_i.min1)
{
  a.theta_i <- 1
  b.theta_i <- 1
  n0.theta_i <- 1

  N.g <- nrow(Z.g)
  K <- ncol(Z.g)
  sig2k.g <- rep(1, K)
  var.g <- 1/K + tau.g
  theta.g <- mean(theta_i.min1)

  # Hyper parameters
  SS <- b.theta_i + sum((theta_i.min1 - theta.g)^2) #+ (N.g*n0.theta_i*theta.g)/(2*(N.g + n0.theta_i))
  var0.theta <- 1 / rgamma(1,(N.g + a.theta_i)/2, SS/2)
  mu0.theta <- rnorm(1, (N.g*var0.theta*theta.g)/(var0.theta*(N.g + n0.theta_i)), sqrt(1/(var0.theta*(N.g + n0.theta_i))))

  var.theta_i <- 1/(N.g*(var.g) + 1/var0.theta)
  mu.theta_i <- rep(NA, N.g)
  for(zi in 1:N.g) {
    mu.theta_i[zi] <- var.theta_i*((N.g*(mean(beta) + mean(Z.g[zi,])))*(var.g) + mu0.theta/var0.theta)
  }

  # Draw N person ability parameters
  theta_i <- rnorm(N.g, mu.theta_i, sqrt(var.theta_i))
  theta_i <- theta_i - mean(theta_i)

  return(theta_i)
}

#' @export
sampleZ <- function(Y, Z, ZT, mu.Z, mu.ZT, partMatrix, u, theta, likelihood = FALSE)
{
  K <- ncol(Z)
  N <- nrow(Z)

  var.Z <- numeric(K)
  lik_k <- rep(NA, K) # Probability density for Z_k
  tmp <- vector("list", K)
  for (k in 1:K) {
    #tmp[[k]] <- partMatrix[[k]]$B12 %*% partMatrix[[k]]$B22.inv %*% t(Z[, -k] - mu.Z[, -k])
    tmp[[k]] <- partMatrix[[k]]$B12 %*% partMatrix[[k]]$B22.inv %*% t(ZT[, -k] - mu.ZT[, -k])
    mu.Z[, k] <-  mu.ZT[, k] +  tmp[[k]]
    var.Z[k] <- partMatrix[[k]]$B11 - partMatrix[[k]]$B12 %*% partMatrix[[k]]$B22.inv %*% partMatrix[[k]]$B21

    if(is.nan(var.Z[k]) || (var.Z[k] < 0))
       return(NULL)

    # Sample Zk|Z.mink, T1..p
    #Z[Y[,k]==0, k] <- truncnorm::rtruncnorm(n = length(mu.Z[Y[, k]==0, k]), mean = mu.Z[Y[, k]==0, k], sd = sqrt(var.Z[k]), a = -Inf, b = 0)
    #Z[Y[,k]==1, k] <- truncnorm::rtruncnorm(n = length(mu.Z[Y[, k]==1, k]), mean = mu.Z[Y[, k]==1, k], sd = sqrt(var.Z[k]), a = 0, b = Inf)
    Z[Y[,k]==0, k] <- extraDistr::rtnorm(n = length(mu.Z[Y[, k]==0, k]), mean = mu.Z[Y[, k]==0, k], sd = sqrt(var.Z[k]), a = -Inf, b = 0)
    Z[Y[,k]==1, k] <- extraDistr::rtnorm(n = length(mu.Z[Y[, k]==1, k]), mean = mu.Z[Y[, k]==1, k], sd = sqrt(var.Z[k]), a = 0, b = Inf)

    #if(any(Z[Y[,k]==0, k] > 0) || any(Z[Y[,k]==1, k] < 0))
    #  browser()

    # Compute probability density
    if(likelihood) {
      #lik0 <- sum(log(truncnorm::dtruncnorm(x = Z[Y[,k]==0, k], mean = mu.Z[Y[, k]==0, k], sd = sqrt(var.Z[k]), a = -Inf, b = 0)))
      #lik1 <- sum(log(truncnorm::dtruncnorm(x = Z[Y[,k]==1, k], mean = mu.Z[Y[, k]==1, k], sd = sqrt(var.Z[k]), a = 0, b = Inf)))
      lik0 <- sum(extraDistr::dtnorm(x = Z[Y[,k]==0, k], mean = mu.Z[Y[, k]==0, k], sd = sqrt(var.Z[k]), a = -Inf, b = 0, log = TRUE))
      lik1 <- sum(extraDistr::dtnorm(x = Z[Y[,k]==1, k], mean = mu.Z[Y[, k]==1, k], sd = sqrt(var.Z[k]), a = 0, b = Inf, log = TRUE))
      lik_k[k] <- lik0 + lik1
    }

  }
  return(list(Z = Z, mu.Z = mu.Z, var.Z = var.Z, tmp = tmp, lik_k = lik_k))
}

#' @export
sampleBeta <- function(Z, beta, theta = 0, tau, a.beta = 0.5, b.beta = 0.5, n0.beta = 1)
{
  K <- ncol(Z)
  N <- nrow(Z)

  #ones <- rep(1, K)
  #varinv <- diag(1/(rep(1, K) + tau))
  #var.gen <- (t(ones) %*% varinv) %*% ones
  var.gen <- 1/(1 + tau) #rep(1/(1 + tau), K)

  ### Sample item difficulty paramaters ###

  # Hyper parameters
  SS <- b.beta + sum((beta - mean(beta))^2) + (K*n0.beta*mean(beta))/(2*(K + n0.beta))
  var0.beta <- 1 / rgamma(1, (K + a.beta)/2, SS/2)
  mu0.beta <- rnorm(1, (-K*var0.beta*mean(beta))/(var0.beta*(K + n0.beta)), sqrt(1/(var0.beta*(K + n0.beta))))

  #print(mu0.beta)
  var0.beta <- 10^10
  mu0.beta <- 0

  var.beta <- 1/(N*(var.gen) + 1/var0.beta)
  mu.beta <- var.beta*((N*(mean(theta) - colMeans(Z)))*(var.gen) + mu0.beta/var0.beta)

  # Draw K item difficulty parameters
  beta <- rnorm(K, mu.beta, sqrt(var.beta))

  return(beta = beta)
}

#' @export
sampleLambda <- function(RT, lambda, zeta = 0, sig2k, delta, a.lambda = 0.5, b.lambda = 0.5, n0.lambda = 1)
{
  K <- ncol(RT)
  N <- nrow(RT)

  #ones <- rep(1, K)
  #varinv <- diag(1/(sig2k[1:K] + delta))
  #var.gen <- (t(ones) %*% varinv) %*% ones
  var.gen <- 1/(sig2k[1:K] + delta)

  ### Sample item time intensity paramaters ###

  # Hyper parameters
  SS <- b.lambda + sum((lambda - mean(lambda))^2) + (K*n0.lambda*mean(lambda))/(2*(K + n0.lambda))
  var0.lambda <- 1 / rgamma(1, (K + a.lambda)/2, SS/2)
  mu0.lambda <- rnorm(1, (K*var0.lambda*mean(lambda))/(var0.lambda*(K + n0.lambda)), sqrt(1/(var0.lambda*(K + n0.lambda))))

  var.lambda <- 1/(N*(var.gen) + 1/var0.lambda)
  mu.lambda <- var.lambda*((N*(colMeans(RT) + mean(zeta)))*(var.gen) + mu0.lambda/var0.lambda)

  # Draw K time intensity parameters
  lambda <- rnorm(K, mu.lambda, sqrt(var.lambda))

  return(lambda = lambda)
}

#' @export
sampleTheta <- function(Z, beta, tau)
{
  K <- ncol(Z)
  N <- nrow(Z)

  #ones <- rep(1, K)
  #varinv <- diag(1/(rep(1, K) + tau))
  #var.gen <- (t(ones) %*% varinv) %*% ones
  var.gen <- 1/((1 + K * tau) / K) #sum(rep(1/(1 + tau), K))


  ### Sample group ability parameter ###

  # Hyper parameters
  var0.theta <- 10^10
  mu0.theta <- 0

  var.theta <- 1/(N*(var.gen) + 1/var0.theta)
  mu.theta <- var.theta*((N*(mean(beta) + mean(Z)))*(var.gen) + mu0.theta/var0.theta)
  #mu.theta <- var.theta*((N*(mean(Z)))*(var.gen) + mu0.theta/var0.theta)


  #theta <- mean(Z)

  # Draw group ability mean
  theta <- rnorm(1, mu.theta, sqrt(var.theta))

  #print(theta)

  return(theta = theta)
}

#' @export
sampleZeta <- function(RT, lambda, sig2k, delta)
{
  K <- ncol(RT)
  N <- nrow(RT)

  #ones <- rep(1, K)
  #varinv <- diag(1/(sig2k[1:K] + delta))
  #var.gen <- (t(ones) %*% varinv) %*% ones
  #var.gen <- sum(1/(sig2k[1:K]/K^2 + delta/K))
  var.gen <- 1/(sum(sig2k[1:K]/K^2 + delta/K))
  #print(var.gen)

  ### Sample group speed parameter ###

  # Hyper parameters
  var0.zeta <- 10^10
  mu0.zeta <- 0

  var.zeta <- 1/(N*(var.gen) + 1/var0.zeta)
  mu.zeta <- var.zeta*((N*(mean(lambda) - mean(RT)))*(var.gen) + mu0.zeta/var0.zeta)

  # Draw group speed mean
  zeta <- rnorm(1, mu.zeta, sqrt(var.zeta))

  return(zeta = zeta)
}

#' @export
sampleMarAbilityModel <- function(Y, Z.mar, beta.mar, theta.mar, tau.mar, firstGroup = TRUE, init = TRUE, a.beta = 0.5, b.beta = 0.5, n0.beta = 1)
{
  K <- ncol(Z.mar)
  N <- nrow(Z.mar)

  #ones <- rep(1, K)
  #varinv <- diag(1/(rep(1, K) + tau.mar))
  #var.gen <- (t(ones) %*% varinv) %*% ones
  var.gen <- rep(1/(1 + tau.mar), K)

  # if(init) {

    ### Sample item difficulty parameters ###

    # Hyper parameters
    SS <- b.beta + sum((beta.mar - mean(beta.mar))^2) + (K*n0.beta*mean(beta.mar))/(2*(K + n0.beta))
    var0.beta <- 1 / rgamma(1, (K + a.beta)/2, SS/2)
    mu0.beta <- rnorm(1, (K*var0.beta*mean(beta.mar))/(var0.beta*(K + n0.beta)), sqrt(1/(var0.beta*(K + n0.beta))))

    var.beta <- 1/(N*(var.gen) + 1/var0.beta)
    mu.beta <- var.beta*((N*(mean(theta.mar) - colMeans(Z.mar)))*(var.gen) + mu0.beta/var0.beta)

    # Draw K time intensity parameters
    beta.mar <- rnorm(K, mu.beta, sqrt(var.beta))
  # }


  if(init) {

    ### Sample ability parameters group means ###

    if(firstGroup) {
      theta.mar <- 0 # First group
    }
    else {
      var.gen <- 1/((1 + K * tau.mar) / K) #sum(rep(1/(1 + tau), K))

      # Hyper parameters
      #flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
      var0.theta <- 10^10
      mu0.theta <- 0

      var.theta <- 1/(N/(var.gen) + 1/var0.theta)
      mu.theta <- var.theta*((N*(mean(beta.mar) + mean(Z.mar)))/(var.gen) + mu0.theta/var0.theta)

      # Draw ability parameter group mean
      theta.mar <- rnorm(1, mu.theta, sqrt(var.theta))
    }
  }


  #if(init) {

    ### Sample latent responses ###

    Sjc <- matrix(tau.mar / (1 + (K - 1) * tau.mar), ncol = 1, nrow = K - 1)
    var.Z.mar <- (1 + K * tau.mar)/(1 + (K - 1) * tau.mar)
    theta1 <- matrix(theta.mar, ncol = K - 1, nrow = N)
    for (kk in 1:K){
      beta1 <- beta.mar[-kk]
      Z1 <- Z.mar[, -kk] # Latent responses to all but the current item
      mu.Z.mar <- (- beta.mar[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
      Z.mar[Y[,kk]==0, kk] <- extraDistr::rtnorm(n = length(mu.Z.mar[Y[, kk]==0]), mean = mu.Z.mar[Y[, kk]==0], sd = sqrt(var.Z.mar), a = -Inf, b = 0)
      Z.mar[Y[,kk]==1, kk] <- extraDistr::rtnorm(n = length(mu.Z.mar[Y[, kk]==1]), mean = mu.Z.mar[Y[, kk]==1], sd = sqrt(var.Z.mar), a = 0, b = Inf)
#browser()
      #Z.mar[Y[, kk]==0, kk] <- qnorm(runif(N, 0, pnorm(0, mu.Z.mar, sqrt(var.Z.mar))), mu.Z.mar, sqrt(var.Z.mar))[Y[, kk] == 0]
      #Z.mar[Y[, kk]==1, kk] <- qnorm(runif(N, pnorm(0, mu.Z.mar, sqrt(var.Z.mar)),1), mu.Z.mar, sqrt(var.Z.mar))[Y[, kk] == 1]
    }
  #}

  ### Sample covariance tau ###

  mean.person <- apply(Z.mar, 1, mean)
  SSb.tau <- sum((mean.person - mean(Z.mar))^2)

  # Draw covariance parameter
  #tau.mar <- 1 / rgamma(1, N/2, SSb.tau/2) - 1/K
  tau.mar <- extraDistr::rinvgamma(1, N/2, SSb.tau/2) - 1/K

  return(list(Z.mar = Z.mar, beta.mar = beta.mar, theta.mar = theta.mar, tau.mar = tau.mar,
              q.tau = list(dist = "invgamma + 1/K", shape = N/2, rate = SSb.tau/2)))
}

#' @export
sampleMarSpeedModel <- function(RT, lambda, zeta = 0, delta.mar, sig2k.mar, a.sig2 = 0.001, b.sig2 = 0.001)
{
  K <- ncol(RT)
  N <- nrow(RT)

  errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) + mean(zeta)
  mean.person <- apply(errors, 1, mean)

  ### Sample measurement error variance parameters ###

  # Within sum of squares
  SSw.sig2 <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),1,sum))
  SSwk.sig2k <- apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),2,sum)

  # Draw mean measurement error variance
  #sig2.mar <- ( 1 / rgamma(1, a.sig2 + (N*(K-1))/2, b.sig2 + SSw.sig2/2) )
  sig2.mar <- extraDistr::rinvgamma(1, a.sig2 + (N*(K-1))/2, b.sig2 + SSw.sig2/2)

  # Draw K individual measurement error variances
  #sig2k.mar <- (1 / rgamma(K, a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk.sig2k/2))
  sig2k.mar <- extraDistr::rinvgamma(K, a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk.sig2k/2)

  ### Sample covariance parameter ###

  # Between sum of squares
  mean.person <- apply(errors,1,mean)
  SSb.delta <- sum((mean.person - mean(errors))^2)

  #print(mean(1 / rgamma(10000, N/2, SSb/2) - sig2/K))
  #delta.mar <- 1 / rgamma(1, N/2, SSb.delta/2) - sig2.mar/K
  delta.mar <- extraDistr::rinvgamma(1, N/2, SSb.delta/2) - sig2.mar/K

  return(list(sig2.mar = sig2.mar, sig2k.mar = sig2k.mar, delta.mar = delta.mar,
              q.delta = list(dist = "invgamma + sig2/K", shape = N/2, rate = SSb.delta/2),
              q.sig2k = list(dist = "invgamma", shape = a.sig2 + (N*(K-1)/K)/2, rate = b.sig2 + SSwk.sig2k/2)))
}

#' @export
sampleCorrelationNu <- function(r, x, var.Z)
{
  K <- ncol(r)
  N <- nrow(r)

  ### Sample covariance responses-response times ###
  var0.b <- 10^10
  mu0.b <- 0

  tmp1 <- tmp2 <- rep(0, K)
  nu <- var.b <- mu.b <- numeric(K)
  for (k in 1:K) {
    for (i in 1:N) {
      tmp1[k] <- tmp1[k] + x[i, k] * 1/var.Z[k] * x[i, k]
      tmp2[k] <- tmp2[k] + x[i, k] * 1/var.Z[k] * r[i, k]
    }

    var.b[k] <- (1/(1/var0.b + tmp1[k]))
    mu.b[k] <- (var.b[k] * (mu0.b/var0.b + tmp2[k]))
    nu[k] <- rnorm(1, mean = mu.b[k], sd = sqrt(var.b[k]))
  }
#browser()
#print(mu.b[1:3])
#print(var.b[1:3])
#print(mean(diag(cov(r,x))))
  #print(var.Z)

  #print(mean(nu))

  return(list(nu = nu,  q.nu = list(dist = "rnorm", mu = mu.b, sd = sqrt(var.b))))
}








