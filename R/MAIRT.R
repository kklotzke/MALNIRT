#' @export
MAIRT <- function(Y, Group = NULL, data, XG = 1000, burnin = 0.10, inits.1 = NULL, inits.2 = NULL, hyper.priors = NULL, est.person = TRUE, silent = FALSE) {

  ###### Initialize ######

  if (!missing(data)) {
    # Try to find Y and Group in the data set first
    tryCatch(Y <- eval(substitute(Y), data), error=function(e) NULL)
    if (!is.null(Group)) {
      tryCatch(Group <- eval(substitute(Group), data), error=function(e) NULL)
    }
  } else {
    data = NULL
  }

  N <- nrow(Y) # Number of persons
  K <- ncol(Y) # Number of items

  # Multiple groups
  if (!is.null(Group)){
    groups <- unique(Group)
    G <- length(groups)
    Y.group <- vector("list", G)
    for (i in 1:G)
    {
      Y.group[[i]] <- subset(Y, Group == groups[i])
      Y.group[[i]] <- Y.group[[i]]
    }
  }
  else {
    G <- 1
    Y.group <- vector("list", 1)
    Y.group[[1]] <- Y
  }

  # Read hyperpriors
  if (is.null(hyper.priors)) {
    hyper.priors <- c(0.5, 0.5, 1, 1, 1, 1, 0.001, 0.001)
  }
  a.beta <- hyper.priors[1]
  b.beta <- hyper.priors[2]
  n0.beta <- hyper.priors[3]
  a.theta <- hyper.priors[4]
  b.theta <- hyper.priors[5]
  n0.theta <- hyper.priors[6]
  a.sig2 <- hyper.priors[7]
  b.sig2 <- hyper.priors[8]

  # Initialize chains
  chains.list <- vector("list", 2) # Two chains
  chains.list[[1]] <- vector("list", 6)
  chains.list[[1]][[1]] <- matrix(NA, nrow = XG, ncol = K) # Item difficulty
  chains.list[[1]][[2]] <- matrix(NA, nrow = XG, ncol = G) # Ability group mean
  chains.list[[1]][[3]] <- array(1, dim = c(XG, K, G+1)) # Measurement variance per item, first block across all groups - fixed to 1
  chains.list[[1]][[4]] <- array(1, dim = c(XG, 1, G+1)) # Average measurement variance, first block across all groups - fixed to 1
  chains.list[[1]][[5]] <- array(NA, dim = c(XG, 1, G+1)) # Response covariance per group, first block across all groups
  chains.list[[1]][[6]] <- array(NA, dim = c(XG, N, G)) # Person ability parameter
  chains.list[[1]][[7]] <- matrix(0, ncol = K, nrow = N) # Z matrix
  chains.list[[1]][[8]] <- vector("list", G) # List to contain Z matrices per group
  chains.list[[2]] <- chains.list[[1]]
  if (is.null(inits.1)) {
    inits.1 <- vector("list", 6)
    inits.1[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.1[[2]] <- rnorm(G, 0, 1) # Abiliy group mean
    inits.1[[3]] <- runif(K, 0.5, 1.5) # Measurement variance per item
    inits.1[[4]] <- mean(inits.1[[3]]) # Average measurement variance
    inits.1[[5]] <- runif(1, 0, 1) # Response covariance
    inits.1[[6]] <- rnorm(N, 0, 1) # Person ability parameter
  }
  if (is.null(inits.2)) {
    inits.2 <- vector("list", 6)
    inits.2[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.2[[2]] <- rnorm(G, 0, 1) # Abiliy group mean
    inits.2[[3]] <- runif(K, 0.5, 1.5) # Measurement variance per item
    inits.2[[4]] <- mean(inits.1[[3]]) # Average measurement variance
    inits.2[[5]] <- runif(1, 0, 1) # Response covariance
    inits.2[[6]] <- rnorm(N, 0, 1) # Person ability parameter
  }
  # First row are the initial values
  for (i in 1:2) {
    chains.list[[1]][[i]][1, ] <- inits.1[[i]]
    chains.list[[2]][[i]][1, ] <- inits.2[[i]]

    for (gg in 1:G) {
      chains.list[[i]][[8]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
    }
  }

  for (gg in 1:(G)) {
    chains.list[[1]][[5]][1,,gg] <- inits.1[[5]]
    chains.list[[2]][[5]][1,,gg] <- inits.2[[5]]
    chains.list[[1]][[5]][1,,gg+1] <- inits.1[[5]]
    chains.list[[2]][[5]][1,,gg+1] <- inits.2[[5]]
    chains.list[[1]][[6]][1,,gg] <- inits.1[[6]]
    chains.list[[2]][[6]][1,,gg] <- inits.2[[6]]
  }

  # Create helmert matrix
  hmat <- helmert(K)

  # Fractional Bayes factor
  m0 <- m1 <- FBF <- matrix(NA, nrow = XG-1, ncol = 2)


  ###### Run MCMC ######
  # For each chain
  for (cc in 1:2) {
    chain <- chains.list[[cc]] # Read the cc-th chain

    Z <- chain[[7]]
    Z.group <- chain[[8]]

    # For each iteration
    for (ii in 2:XG) {

      # Read values from former chain (t-1)
      beta.min1 <- chain[[1]][ii-1, ]
      theta.min1 <- chain[[2]][ii-1, ]
      sig2k.min1 <- chain[[3]][ii-1,,1]
      sig2.min1 <- chain[[4]][ii-1,,1]
      tau.min1 <- chain[[5]][ii-1,,]

      # Generalize sigma^2 and tau
      ones <- rep(1, K)
      varinv <- diag(1/(sig2k.min1[1:K] + tau.min1[1]))
      var.gen <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          sig2k.min1.g <- chain[[3]][ii-1,,gg+1]
          varinv <- diag(1/(sig2k.min1.g[1:K]))
          tmp <- (t(ones) %*% varinv) %*% ones
          tmp <- tmp + tau.min1[gg+1]
          var.gen <- c(var.gen, tmp)
        }
      }

      ### Sample latent response variable ###

      # See Fox et al. (2016) Equation 14 and 15
      Sjc <- matrix(tau.min1[1] / (1 + (K - 1) * tau.min1[1]), ncol = 1, nrow = K - 1)
      var.Z <- (1 + K * tau.min1[1])/(1 + (K - 1) * tau.min1[1])
      theta1 <- matrix(mean(theta.min1), ncol = K - 1, nrow = N)

      for (kk in 1:K){
        beta1 <- beta.min1[-kk]
        Z1 <- Z[, -kk] # Latent responses to all but the current item
        mu.Z <- (mean(theta.min1) - beta.min1[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
        Z[Y[, kk]==0, kk] <- qnorm(runif(N, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y[, kk] == 0]
        Z[Y[, kk]==1, kk] <- qnorm(runif(N, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y[, kk] == 1]
        #boundary <- 0
        #Z[Y[, kk]==0, kk] <- truncnorm::rtruncnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), a = -Inf, b = boundary)[Y[, kk]==0]
        #Z[Y[, kk]==1, kk] <- truncnorm::rtruncnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), a = boundary, b = Inf)[Y[, kk]==1]
        #Z[Y[, kk]==0, kk] <- msm::rtnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), lower = -Inf, upper = boundary)[Y[, kk]==0]
        #Z[Y[, kk]==1, kk] <- msm::rtnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), lower = boundary, upper = Inf)[Y[, kk]==1]
      }

      #print(mean(Z[, 1]))

      for (gg in 1:G) {
        N.g <- nrow(Y.group[[gg]])

        Sjc <- matrix(tau.min1[gg+1] / (1 + (K - 1) * tau.min1[gg+1]), ncol = 1, nrow = K - 1)
        var.Z <- (1 + K * tau.min1[gg+1])/(1 + (K - 1) * tau.min1[gg+1])
        theta1 <- matrix(theta.min1[gg], ncol = K - 1, nrow = N.g)

        for (kk in 1:K){
          beta1 <- beta.min1[-kk]
          Z1 <- Z.group[[gg]][, -kk] # Latent responses to all but the current item
          mu.Z <- (theta.min1[gg] - beta.min1[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
          boundary <- 0#mu.Z
          #Z.group[[gg]][Y.group[[gg]][, kk]==0, kk] <- qnorm(runif(N.g, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 0]
          #Z.group[[gg]][Y.group[[gg]][, kk]==1, kk] <- qnorm(runif(N.g, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 1]
          Z.group[[gg]][Y.group[[gg]][, kk]==0, kk] <- truncnorm::rtruncnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), a = -Inf, b = boundary)[Y.group[[gg]][, kk]==0]
          Z.group[[gg]][Y.group[[gg]][, kk]==1, kk] <- truncnorm::rtruncnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), a = boundary, b = Inf)[Y.group[[gg]][, kk]==1]
          #Z.group[[gg]][Y.group[[gg]][, kk]==0, kk] <- msm::rtnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), lower = -Inf, upper = boundary)[Y.group[[gg]][, kk]==0]
          #Z.group[[gg]][Y.group[[gg]][, kk]==1, kk] <- msm::rtnorm(n = N, mean = mu.Z, sd = sqrt(var.Z), lower = boundary, upper = Inf)[Y.group[[gg]][, kk]==1]
        }
      }

      ### Sample ability parameters group means ###

      # Hyper parameters
      #flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
      var0.theta <- 10^10
      mu0.theta <- 0

      var.theta <- rep(NA, G)
      mu.theta <- rep(NA, G)
      var.theta[1] <- 0
      mu.theta[1] <- 0
      if(G > 1) {
        for(gg in 2:G) {
          N.g <- nrow(Z.group[[gg]])
          var.theta[gg] <- 1/(N.g/(var.gen[gg+1]) + 1/var0.theta)
          mu.theta[gg] <- var.theta[gg]*((N.g*(mean(beta.min1) + mean(Z.group[[gg]])))/(var.gen[gg+1]) + mu0.theta/var0.theta)
        }
      }

      # Draw G ability parameter group means
      theta <- rnorm(G, mu.theta, sqrt(var.theta))
      chain[[2]][ii, ] <- theta


      ### Sample item difficulty paramaters ###

      # Hyper parameters
      SS <- b.beta + sum((beta.min1 - mean(beta.min1))^2) + (K*n0.beta*mean(beta.min1))/(2*(K + n0.beta))
      var0.beta <- 1 / rgamma(1, (K + a.beta)/2, SS/2)
      mu0.beta <- rnorm(1, (K*var0.beta*mean(beta.min1))/(var0.beta*(K + n0.beta)), sqrt(1/(var0.beta*(K + n0.beta))))

      var.beta <- 1/(N*(var.gen[1]) + 1/var0.beta)
      mu.beta <- var.beta*((N*(mean(theta) - colMeans(Z)))*(var.gen[1]) + mu0.beta/var0.beta)

      # Draw K time intensity parameters
      chain[[1]][ii, ] <- beta <- rnorm(K, mu.beta, sqrt(var.beta))


      ## Across all groups

      ### Sample measurement error variances

      # Fixed to 1

      sig2 <- 1
      sig2k <- rep(1, K)

      ### Sample covariance parameter

      # Helmert transformation
      errors <- Z + matrix(beta, nrow = N, ncol = K, byrow = TRUE)
      tmat <- errors %*% t(hmat)

      # Between sum of squares
      mean.person <- apply(errors,1,mean)
      SSb <- sum((mean.person - mean(errors))^2)

      # Draw covariance parameter
      chain[[5]][ii,,1] <- tau <- 1 / rgamma(1, N/2, SSb/2) - sig2/K


      ## For each group
      for(gg in 1:G) {
        N.g <- nrow(Z.group[[gg]])

        # Helmert transformation
        errors <- Z.group[[gg]] + matrix(beta, nrow = N.g, ncol = K, byrow = TRUE) + theta[gg]
        tmat <- errors %*% t(hmat)

        # Between sum of squares
        SSb <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K

        # Draw total measurement error variance
        chain[[4]][ii,,gg+1] <- sig2 <- 1 #1 / rgamma(1, (N.g+a.sig2)/2, (SSw+b.sig2)/2)

        # Draw K individual measurement error variances
        chain[[3]][ii,,gg+1] <- sig2k <- rep(1, K) # 1 / rgamma(K, (N.g+a.sig2)/2, (SSwk+b.sig2)/2)

        # Draw covariance parameter
        chain[[5]][ii,,gg+1] <- tau <- 1 / rgamma(1, N.g/2, SSb/2) - sig2/K
      }


      if ((!silent) && (ii%%100 == 0))
        cat("Iteration ", ii, " ", "\n")
      flush.console()
    }

    chain[[7]] <- Z
    chain[[8]] <- Z.group
    chains.list[[cc]] <- chain
  }

  # Number of burnin iterations
  XG.burnin <- floor(XG * burnin)

  # Merge chains
  chain.1 <- chains.list[[1]]
  chain.2 <- chains.list[[2]]
  post.beta <- colMeans((chain.1[[1]][XG.burnin:XG, ] + chain.2[[1]][XG.burnin:XG, ]) / 2)
  if(G > 1) {
    post.theta <- colMeans((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
  }
  else {
    post.theta <- mean((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
  }
  post.sig2k <- colMeans((chain.1[[3]][XG.burnin:XG,,] + chain.2[[3]][XG.burnin:XG,,]) / 2)
  post.sig2 <- colMeans((chain.1[[4]][XG.burnin:XG,,] + chain.2[[4]][XG.burnin:XG,,]) / 2)
  post.tau <- colMeans((chain.1[[5]][XG.burnin:XG,,] + chain.2[[5]][XG.burnin:XG,,]) / 2)
  post.Z <- (chain.1[[7]] + chain.2[[7]]) / 2
  post.Z.group <- vector("list", G)
  for (gg in 1:G) {
    post.Z.group[[gg]] <- (chain.1[[8]][[gg]] + chain.2[[8]][[gg]]) / 2
  }

  ### Sample person ability parameter

  if(est.person)
  {
    a.theta_i <- 1
    b.theta_i <- 1
    n0.theta_i <- 1

    # For each chain
    for (cc in 1:2) {
      chain <- chains.list[[cc]] # Read the cc-th chain

      # For each iteration
      for (ii in 2:XG) {
        for (gg in 1:G)
        {
          Z.group <- post.Z.group
          N.g <- nrow(Z.group[[gg]])
          sig2k.g <- post.sig2k[,gg+1] #chain[[3]][ii,,gg+1]
          tau.g <- post.tau[gg+1]
          theta.g <- post.theta[gg]
          beta <- post.beta
          theta_i.min1 <- chain[[6]][ii-1,,gg]
          ones <- rep(1, K)
          varinv <- diag(1/(sig2k.g[1:K]))
          var.g <- (t(ones) %*% varinv) %*% ones

          # Hyper parameters
          SS <- b.theta_i + sum((theta_i.min1 - theta.g)^2) + (N.g*n0.theta_i*theta.g)/(2*(N.g + n0.theta_i))
          var0.theta <- 1 / rgamma(1,(N.g + a.theta_i)/2, SS/2)
          mu0.theta <- rnorm(1, (N.g*var0.theta*theta.g)/(var0.theta*(N.g + n0.theta_i)), sqrt(1/(var0.theta*(N.g + n0.theta_i))))

          var.theta_i <- 1/(N.g*(var.g) + 1/var0.theta)
          mu.theta_i <- rep(NA, N.g)
          for(zi in 1:N.g) {
            mu.theta_i[zi] <- var.theta_i*((N.g*(mean(beta) + mean(Z.group[[gg]][zi,])))*(var.g) + mu0.theta/var0.theta)
          }
          # Draw N person speed parameters
          theta_i <- rnorm(N.g, mu.theta_i, sqrt(var.theta_i))
          chain[[6]][ii,,gg] <- theta_i
        }

        if (ii%%100 == 0)
          cat("Post-Hoc: Iteration ", ii, " ", "\n")
        flush.console()
      }

      chains.list[[cc]] <- chain
    }
  }


  chain.1 <- chains.list[[1]]
  chain.2 <- chains.list[[2]]
  post.theta_i <- colMeans((chain.1[[6]][XG.burnin:XG,,] + chain.2[[6]][XG.burnin:XG,,]) / 2)

  return(list(beta = post.beta, theta = post.theta, sig2k = post.sig2k, sig2 = post.sig2, tau = post.tau, theta_i = post.theta_i))
}
