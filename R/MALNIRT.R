#' @export
MALNIRT <- function(Y, RT, Group = NULL, data, XG = 1000, burnin = 0.10, inits.1 = NULL, inits.2 = NULL, hyper.priors = NULL, est.person = FALSE, silent = FALSE) {

  ###### Initialize ######

  if (!missing(data)) {
    # Try to find Y, RT and Group in the data set first
    tryCatch(Y <- eval(substitute(Y), data), error=function(e) NULL)
    tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
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
    RT.group <- vector("list", G)
    for (i in 1:G)
    {
      Y.group[[i]] <- subset(Y, Group == groups[i])
      Y.group[[i]] <- Y.group[[i]]
      RT.group[[i]] <- subset(RT, Group == groups[i])
      RT.group[[i]] <- RT.group[[i]]
    }
  }
  else {
    G <- 1
    Y.group <- vector("list", 1)
    Y.group[[1]] <- Y
    RT.group <- vector("list", 1)
    RT.group[[1]] <- RT
  }

  # Read hyperpriors
  if (is.null(hyper.priors)) {
    hyper.priors <- c(0.5, 0.5, 1, 1, 1, 1, # Responses
                      0.5, 0.5, 1, 0.1, 0.1, 1, 0.001, 0.001) # RTs
  }
  a.beta <- hyper.priors[1] # Responses
  b.beta <- hyper.priors[2]
  n0.beta <- hyper.priors[3]
  a.theta <- hyper.priors[4]
  b.theta <- hyper.priors[5]
  n0.theta <- hyper.priors[6]
  a.lambda <- hyper.priors[7] # RTs
  b.lambda <- hyper.priors[8]
  n0.lambda <- hyper.priors[9]
  a.zeta <- hyper.priors[10]
  b.zeta <- hyper.priors[11]
  n0.zeta <- hyper.priors[12]
  a.sig2 <- hyper.priors[13]
  b.sig2 <- hyper.priors[14]

  # Initialize chains
  chains.list <- vector("list", 2) # Two chains
  chains.list[[1]] <- vector("list", 12)

  # Item parameters
  chains.list[[1]][[1]] <- matrix(NA, nrow = XG, ncol = K) # Item difficulty
  chains.list[[1]][[2]] <- matrix(NA, nrow = XG, ncol = K) # Time intensity

  # Group parameters
  chains.list[[1]][[3]] <- matrix(NA, nrow = XG, ncol = G) # Ability group means
  chains.list[[1]][[4]] <- matrix(NA, nrow = XG, ncol = G) # Speed group means

  # Person parameters
  chains.list[[1]][[5]] <- array(NA, dim = c(XG, N, G)) # Person ability parameter
  chains.list[[1]][[6]] <- array(NA, dim = c(XG, N, G)) # Person speed parameter

  # Measurement variance parameters
  chains.list[[1]][[7]] <- array(1, dim = c(XG, K, G+1)) # Measurement variance for RT per item, first block across all groups
  chains.list[[1]][[8]] <- array(1, dim = c(XG, 1, G+1)) # Average measurement variance for RT, first block across all groups

  # Covariance parameters
  chains.list[[1]][[9]] <- array(NA, dim = c(XG, 1, G+1)) # Response covariance per group, first block across all groups
  chains.list[[1]][[10]] <- array(NA, dim = c(XG, 1, G+1)) # RT covariance per group, first block across all groups

  # Latent response parameters
  chains.list[[1]][[11]] <- matrix(0, ncol = K, nrow = N) # Z matrix
  chains.list[[1]][[12]] <- vector("list", G) # List to contain Z matrices per group

  # SAT
  chains.list[[1]][[13]] <- matrix(NA, nrow = XG, ncol = 1) # nu

  chains.list[[2]] <- chains.list[[1]]
  if (is.null(inits.1)) {
    inits.1 <- vector("list", 10)
    inits.1[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.1[[2]] <- rnorm(K, 10, 5) # Time intensity
    inits.1[[3]] <- rnorm(G, 0, 1) # Abiliy group means
    inits.1[[4]] <- rnorm(G, 10, 5) # Speed group means
    inits.1[[5]] <- rnorm(N, 0, 1) # Person ability parameter
    inits.1[[6]] <- rnorm(N, 10, 5) # Person speed parameter
    inits.1[[7]] <- runif(K, 0.5, 1.5) # Measurement variance for RT per item
    inits.1[[8]] <- mean(inits.1[[3]]) # Average measurement variance for RT
    inits.1[[9]] <- runif(1, 0, 1) # Response covariance
    inits.1[[10]] <- runif(1, 0, 0.2) # RT covariance
  }
  if (is.null(inits.2)) {
    inits.2 <- vector("list", 10)
    inits.2[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.2[[2]] <- rnorm(K, 10, 5) # Time intensity
    inits.2[[3]] <- rnorm(G, 0, 1) # Abiliy group means
    inits.2[[4]] <- rnorm(G, 10, 5) # Speed group means
    inits.2[[5]] <- rnorm(N, 0, 1) # Person ability parameter
    inits.2[[6]] <- rnorm(N, 10, 5) # Person speed parameter
    inits.2[[7]] <- runif(K, 0.5, 1.5) # Measurement variance for RT per item
    inits.2[[8]] <- mean(inits.1[[3]]) # Average measurement variance for RT
    inits.2[[9]] <- runif(1, 0, 1) # Response covariance
    inits.2[[10]] <- runif(1, 0, 0.2) # RT covariance
  }
  # First row are the initial values
  for (i in 1:4) {
    chains.list[[1]][[i]][1, ] <- inits.1[[i]]
    chains.list[[2]][[i]][1, ] <- inits.2[[i]]
  }
  for (i in 7:10) {
    for (gg in 1:(G+1)) {
      chains.list[[1]][[i]][1,,gg] <- inits.1[[i]]
      chains.list[[2]][[i]][1,,gg] <- inits.2[[i]]
    }
  }

  for (gg in 1:(G)) {
    chains.list[[1]][[5]][1,,gg] <- inits.1[[5]]
    chains.list[[2]][[5]][1,,gg] <- inits.2[[5]]
    chains.list[[1]][[6]][1,,gg] <- inits.1[[6]]
    chains.list[[2]][[6]][1,,gg] <- inits.2[[6]]
    chains.list[[1]][[12]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
    chains.list[[2]][[12]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
  }

  chains.list[[1]][[13]][1,1] <- 0.5
  chains.list[[2]][[13]][1,1] <- -0.5

  # Create helmert matrix
  hmat <- helmert(K)
  tmp1 <- hmat[K,1]
  tmp2 <- hmat[K,K]
  hmat2 <- matrix(tmp1, nrow = K, ncol = K)
  diag(hmat2) <- rep(tmp2, K)

  # Fractional Bayes factor
  m0 <- m1 <- FBF <- matrix(NA, nrow = XG-1, ncol = 2)

  # SAT
  draw.nu2 <- numeric(0)

  ###### Run MCMC ######
  # For each chain
  for (cc in 1:2) {
    chain <- chains.list[[cc]] # Read the cc-th chain

    Z <- chain[[11]]
    Z.group <- chain[[12]]

    # For each iteration
    for (ii in 2:XG) {

      # Read values from former chain (t-1)
      beta.min1 <- chain[[1]][ii-1, ]
      lambda.min1 <- chain[[2]][ii-1, ]
      theta.min1 <- chain[[3]][ii-1, ]
      zeta.min1 <- chain[[4]][ii-1, ]
      sig2k.min1 <- chain[[7]][ii-1,,1]
      sig2.min1 <- chain[[8]][ii-1,,1]
      tau.min1 <- chain[[9]][ii-1,,]
      delta.min1 <- chain[[10]][ii-1,,]
      nu.min1 <- chain[[13]][ii-1, 1]

      # Generalize sigma^2 (fixed to 1) and tau
      ones <- rep(1, K)
      varinv <- diag(1/(rep(1,K) + tau.min1[1]))
      var.gen.Z <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          varinv <- diag(1/(rep(1,K) + tau.min1[gg+1]))
          tmp <- (t(ones) %*% varinv) %*% ones
          var.gen.Z <- c(var.gen.Z, tmp)
        }
      }

      # Generalize sigma^2 and delta
      ones <- rep(1, K)
      varinv <- diag(1/(sig2k.min1[1:K] + delta.min1[1]))
      var.gen.RT <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          sig2k.min1.g <- chain[[7]][ii-1,,gg+1]
          varinv <- diag(1/(sig2k.min1.g[1:K] + delta.min1[gg+1]))
          tmp <- (t(ones) %*% varinv) %*% ones
          var.gen.RT <- c(var.gen.RT, tmp)
        }
      }



      ##### Responses #####

      ### Sample latent response variable ###
Z <- Z.group[[1]] <- data.lnirt$ZRT

      # nu <- 0.6
      # I <- diag(K)
      # J <- matrix(1, K, K)
      # A11 <- I + tau.min1[1]
      # A22 <- diag(sig2k.min1) + delta.min1[1]
      # A11.inv <- I - J/(1/tau.min1[1] + K)  # solve(A11)
      # a <- diag(1/sig2k.min1)
      # A22.inv <- a - (a %*% J %*% a) / (1/delta.min1[1] + sum(1/sig2k.min1)) # solve(A22)
      #
      # mu.z <- theta.min1[1] - mean(beta.min1)
      # mu.RT <- mean(lambda.min1) - zeta.min1[1]
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


      # # See Fox et al. (2016) Equation 14 and 15
      # Sjc <- matrix(tau.min1[1] / (1 + (K - 1) * tau.min1[1]), ncol = 1, nrow = K - 1)
      # var.Z <- (1 + K * tau.min1[1])/(1 + (K - 1) * tau.min1[1])
      # theta1 <- matrix(mean(theta.min1), ncol = K - 1, nrow = N)
      #
      # for (kk in 1:K){
      #   beta1 <- beta.min1[-kk]
      #   Z1 <- Z[, -kk] # Latent responses to all but the current item
      #   mu.Z <- (mean(theta.min1) - beta.min1[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
      #   Z[Y[, kk]==0, kk] <- qnorm(runif(N, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y[, kk] == 0]
      #   Z[Y[, kk]==1, kk] <- qnorm(runif(N, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y[, kk] == 1]
      # }
      #
      # for (gg in 1:G) {
      #   N.g <- nrow(Y.group[[gg]])
      #
      #   Sjc <- matrix(tau.min1[gg+1] / (1 + (K - 1) * tau.min1[gg+1]), ncol = 1, nrow = K - 1)
      #   var.Z <- (1 + K * tau.min1[gg+1])/(1 + (K - 1) * tau.min1[gg+1])
      #   theta1 <- matrix(theta.min1[gg], ncol = K - 1, nrow = N.g)
      #
      #   for (kk in 1:K){
      #     beta1 <- beta.min1[-kk]
      #     Z1 <- Z.group[[gg]][, -kk] # Latent responses to all but the current item
      #     mu.Z <- (theta.min1[gg] - beta.min1[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
      #     Z.group[[gg]][Y.group[[gg]][, kk]==0, kk] <- qnorm(runif(N.g, 0, pnorm(0, mu.Z, sqrt(var.Z))), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 0]
      #     Z.group[[gg]][Y.group[[gg]][, kk]==1, kk] <- qnorm(runif(N.g, pnorm(0, mu.Z, sqrt(var.Z)),1), mu.Z, sqrt(var.Z))[Y.group[[gg]][, kk] == 1]
      #   }
      # }

      ### Sample ability parameters group means ###

      # Hyper parameters
      # flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
      var0.theta <- 10^10
      mu0.theta <- 0

      var.theta <- rep(NA, G)
      mu.theta <- rep(NA, G)
      var.theta[1] <- 0
      mu.theta[1] <- 0
      if(G > 1) {
        for(gg in 2:G) {
          N.g <- nrow(Z.group[[gg]])
          var.theta[gg] <- 1/(N.g/(var.gen.Z[gg+1]) + 1/var0.theta)
          mu.theta[gg] <- var.theta[gg]*((N.g*(mean(beta.min1) + mean(Z.group[[gg]])))/(var.gen.Z[gg+1]) + mu0.theta/var0.theta)
        }
      }

      # Draw G speed parameter group means
      theta <- rnorm(G, mu.theta, sqrt(var.theta))
      chain[[3]][ii, ] <- theta


      ### Sample item difficulty paramaters ###

      Z1 <- data.lnirt$RTZ[,1:K]
      # Hyper parameters
      SS <- b.beta + sum((beta.min1 - mean(beta.min1))^2) + (K*n0.beta*mean(beta.min1))/(2*(K + n0.beta))
      var0.beta <- 1 / rgamma(1, (K + a.beta)/2, SS/2)
      mu0.beta <- rnorm(1, (K*var0.beta*mean(beta.min1))/(var0.beta*(K + n0.beta)), sqrt(1/(var0.beta*(K + n0.beta))))

      var.beta <- 1/(N*(var.gen.Z[1]) + 1/var0.beta)
      mu.beta <- var.beta*((N*(mean(theta) - colMeans(Z1)))*(var.gen.Z[1]) + mu0.beta/var0.beta)

      # Draw K time intensity parameters
      chain[[1]][ii, ] <- beta <- rnorm(K, mu.beta, sqrt(var.beta))


      ### Sample covariance parameter

      # Helmert transformation
      errors <- Z1 + matrix(beta, nrow = N, ncol = K, byrow = TRUE)
      tmat <- errors %*% t(hmat)
      tmp.z2 <- tmat[, 2]

      # Between sum of squares
      mean.person <- apply(errors,1,mean)
      SSb <- sum((mean.person - mean(errors))^2)

      # Draw covariance parameter
      chain[[9]][ii,,1] <- tau <- 1 / rgamma(1, N/2, SSb/2) - 1/K


      ## For each group
      for(gg in 1:G) {
        N.g <- nrow(Z.group[[gg]])

        # Helmert transformation
        errors <- Z1 + matrix(beta, nrow = N.g, ncol = K, byrow = TRUE) + theta[gg] #Z.group[[gg]]
        tmat <- errors %*% t(hmat)

        # Between sum of squares
        SSb <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K

        # Draw covariance parameter
        chain[[9]][ii,,gg+1] <- tau <- 1 / rgamma(1, N.g/2, SSb/2) - 1/K
      }




      ##### Response times #####

      ### Sample speed parameters group means ###

      # Hyper parameters
      # flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
      var0.zeta <- 10^10
      mu0.zeta <- 0

      var.zeta <- rep(NA, G)
      mu.zeta <- rep(NA, G)
      var.zeta[1] <- 0
      mu.zeta[1] <- 0
      if(G > 1) {
        for(gg in 2:G) {
          N.g <- nrow(RT.group[[gg]])
          var.zeta[gg] <- 1/(N.g/(var.gen.RT[gg+1]) + 1/var0.zeta)
          mu.zeta[gg] <- var.zeta[gg]*((N.g*(mean(lambda.min1) - mean(RT.group[[gg]])))/(var.gen.RT[gg+1]) + mu0.zeta/var0.zeta)
        }
      }

      # Draw G speed parameter group means
      zeta <- rnorm(G, mu.zeta, sqrt(var.zeta))
      chain[[4]][ii, ] <- zeta


      ### Sample item time intensity paramaters ###

      # Hyper parameters
      SS <- b.lambda + sum((lambda.min1 - mean(lambda.min1))^2) + (K*n0.lambda*mean(lambda.min1))/(2*(K + n0.lambda))
      var0.lambda <- 1 / rgamma(1, (K + a.lambda)/2, SS/2)
      mu0.lambda <- rnorm(1, (K*var0.lambda*mean(lambda.min1))/(var0.lambda*(K + n0.lambda)), sqrt(1/(var0.lambda*(K + n0.lambda))))

      var.lambda <- 1/(N*(var.gen.RT[1]) + 1/var0.lambda)
      mu.lambda <- var.lambda*((N*(colMeans(RT) + mean(zeta)))*(var.gen.RT[1]) + mu0.lambda/var0.lambda)

      # Draw K time intensity parameters
      chain[[2]][ii, ] <- lambda <- rnorm(K, mu.lambda, sqrt(var.lambda))


      ### Sample measurement error variances

      # Within sum of squares
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta)
      tmat <- errors %*% t(hmat2)

      mean.person <- apply(errors,1,mean)
      mean.item <- apply(errors,2,mean)
      #SSw <- sum((errors - matrix(mean.item, ncol = K, nrow = N, byrow = TRUE))^2)
      #SSwk <- colSums((errors - matrix(mean.item, ncol = K, nrow = N, byrow = TRUE))^2)
      SSwk <- numeric(K)
      for (kk in 1:K) {
        for (i in 1:N) {
          a <- mean(tmat[i,])
          SSwk[kk] <- SSwk[kk] + (tmat[i, kk] - a)^2
        }
      }

      SSw <- sum(SSwk)

      # Draw total measurement error variance
      chain[[8]][ii,,1] <- sig2 <- (( 1 / rgamma(1, a.sig2 + N/2, b.sig2 + SSw/2) ) / K)

      # Draw K individual measurement error variances
      chain[[7]][ii,,1] <- sig2k <- (1 / rgamma(K, a.sig2 + N/2, b.sig2 + SSwk/2))
      #chain[[7]][ii,,1] <- sig2k <- 1 / rgamma(K, (N)/2, (SSwk)/2)
      this.sig2k <- sig2k #<- data.lnirt$sig2k.lnrt


      ### Sample covariance parameter

      # Helmert transformation
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta) #+ matrix(zeta_i, nrow = N, ncol = K)
      tmat <- errors %*% t(hmat)
      tmp.t2 <- tmat[, 2]
      draw.nu2 <- c(draw.nu2, cov(tmp.t2, tmp.z2))

      # Between sum of squares
      mean.person <- apply(errors,1,mean)
      #SSb1 <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K
      SSb <- sum((mean.person - mean(errors))^2)
      #cat(SSb, ",", SSb1, "\n")
      #SSb <- sum((mean.person - (mean(lambda) - mean(zeta_i)))^2)

      # Draw covariance parameter
      chain[[10]][ii,,1] <- delta <- 1 / rgamma(1, N/2, SSb/2) - sig2/K
      #browser()
      #chain[[10]][ii,,1] <- delta <- 1 / rgamma(1, N/2, SSb/2) - sig2/K
      #tmp <- (1 / rgamma(1, N/2, SSb/2))
      #chain[[8]][ii,,1] <- sig2 <- (tmp - delta) * K
      #print(sig2)


      ## For each group
      for(gg in 1:G) {
        N.g <- nrow(RT.group[[gg]])

        # Helmert transformation
        errors <- RT.group[[gg]] - matrix(lambda, nrow = N.g, ncol = K, byrow = TRUE) + zeta[gg]
        tmat <- errors %*% t(hmat)

        # Between sum of squares
        SSb <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K

        # Within sum of squares
        mean.person <- apply(errors,1,mean)
        mean.item <- apply(errors,2,mean)
        SSw <- sum((errors - matrix(mean.person, ncol = K, nrow = N, byrow = TRUE))^2)
        SSwk <- colSums((errors - matrix(mean.person, ncol = K, nrow = N, byrow = TRUE))^2)

        # Draw total measurement error variance
        chain[[8]][ii,,gg+1] <- sig2 <- (( 1 / rgamma(1, a.sig2 + N/2, b.sig2 + SSw/2) ) / K)

        # Draw K individual measurement error variances
        chain[[7]][ii,,gg+1] <- sig2k <- (1 / rgamma(K, a.sig2 + N/2, b.sig2 + SSwk/2))

        # Draw covariance parameter
        chain[[10]][ii,,gg+1] <- delta <- 1 / rgamma(1, N.g/2, SSb/2) - sig2/K
      }

      ##### Joint-model #####

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

      #nu.min1 <- 0.6
      Sigma.zRT <- A11 - nu.min1^2 * A22.inv
      #Sigma.zRT <- A11 - (nu.min1*I) %*% A22.inv %*% (nu.min1*I)
      #Sigma.zRT <- data.lnirt$Sigma.irt - (data.lnirt$Sigma.nu) %*% solve(data.lnirt$Sigma.lnrt) %*% (data.lnirt$Sigma.nu)

      Sigma.zRT.inv <- solve(Sigma.zRT)
#browser()

      tmp1 <- tmp2 <- 0
      for (i in 1:N) {
        xi <- cbind(rep(1,K), A22.inv %*% (RT[i, ] - mu.RT)) #mu.RT))
        #xi <- cbind(rep(1,K), solve(data.lnirt$Sigma.lnrt) %*% (RT[i, ] - mu.RT)) #mu.RT))
        tmp1 <- tmp1 + t(xi) %*%  Sigma.zRT.inv %*% xi
        tmp2 <- tmp2 + t(xi) %*%  Sigma.zRT.inv %*% Z[i, ]
      }

      Sigma.b <- solve(1/var0.b + tmp1)
      mu.b <- Sigma.b %*% (mu0.b/var0.b + tmp2)
#print(mu.b)
      #Sigma.b[1,1] <- 0.01
      #mu.b[1] <- mu.z #mean(Z)

      #Sigma.nu <- solve(N * Sigma.zRT.inv + 1 / var0.nu)

      #browser()
      #print(eigen(Sigma.zRT)$values)
      #print(eigen(Sigma.zRT.inv)$values)
      #print(eigen(Sigma.nu)$values)
      #mu.nu <- Sigma.nu %*% (N * Sigma.zRT.inv %*% ((colMeans(Z) - mu.z) / (Sigma.zRT.inv %*% (colMeans(RT) - mu.RT))) + mu0.nu/var0.nu)

#browser()
      # Draw nu
      b <- mvtnorm::rmvnorm(1, mean=mu.b, sigma=Sigma.b)
      #nu <- rnorm(1, mu.nu, sqrt(Sigma.nu))
      #nu <- mean(nu)
      chain[[13]][ii, 1] <- b[1, 2]
      #print(b)
      #cat(mu.z, "|", b[1], "\n")


      if ((!silent) && (ii%%100 == 0))
        cat("Iteration ", ii, " ", "\n")
      flush.console()
    }

    chain[[11]] <- Z
    chain[[12]] <- Z.group
    chains.list[[cc]] <- chain
  }

  # Number of burnin iterations
  XG.burnin <- floor(XG * burnin)

  # Merge chains
  chain.1 <- chains.list[[1]]
  chain.2 <- chains.list[[2]]
  post.beta <- colMeans((chain.1[[1]][XG.burnin:XG, ] + chain.2[[1]][XG.burnin:XG, ]) / 2)
  post.lambda <- colMeans((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
  if(G > 1) {
    post.theta <- colMeans((chain.1[[3]][XG.burnin:XG, ] + chain.2[[3]][XG.burnin:XG, ]) / 2)
    post.zeta <- colMeans((chain.1[[4]][XG.burnin:XG, ] + chain.2[[4]][XG.burnin:XG, ]) / 2)
  }
  else {
    post.theta <- mean((chain.1[[3]][XG.burnin:XG, ] + chain.2[[3]][XG.burnin:XG, ]) / 2)
    post.zeta <- mean((chain.1[[4]][XG.burnin:XG, ] + chain.2[[4]][XG.burnin:XG, ]) / 2)
  }
  post.sig2k <- colMeans((chain.1[[7]][XG.burnin:XG,,] + chain.2[[7]][XG.burnin:XG,,]) / 2)
  post.sig2 <- colMeans((chain.1[[8]][XG.burnin:XG,,] + chain.2[[8]][XG.burnin:XG,,]) / 2)
  post.tau <- colMeans((chain.1[[9]][XG.burnin:XG,,] + chain.2[[9]][XG.burnin:XG,,]) / 2)
  post.delta <- colMeans((chain.1[[10]][XG.burnin:XG,,] + chain.2[[10]][XG.burnin:XG,,]) / 2)
  post.Z <- (chain.1[[11]] + chain.2[[11]]) / 2
  post.Z.group <- vector("list", G)
  for (gg in 1:G) {
    post.Z.group[[gg]] <- (chain.1[[12]][[gg]] + chain.2[[12]][[gg]]) / 2
  }
  post.nu <- (chain.1[[13]][XG.burnin:XG, 1] + chain.2[[13]][XG.burnin:XG, 1]) / 2


  if(est.person)
  {

    a.theta_i <- 1
    b.theta_i <- 1
    n0.theta_i <- 1
    a.zeta_i <- 1
    b.zeta_i <- 1
    n0.zeta_i <- 1

    # For each chain
    for (cc in 1:2) {
      chain <- chains.list[[cc]] # Read the cc-th chain

      # For each iteration
      for (ii in 2:XG) {
        for (gg in 1:G)
        {
          ### Sample person ability parameter

          Z.group <- post.Z.group
          N.g <- nrow(Z.group[[gg]])
          tau.g <- post.tau[gg+1]
          theta.g <- post.theta[gg]
          beta <- post.beta
          theta_i.min1 <- chain[[5]][ii-1,,gg]
          ones <- rep(1, K)
          varinv <- diag(1/(rep(1, K)))
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
          chain[[5]][ii,,gg] <- theta_i


          ### Sample person speed parameter

          sig2k.g <- post.sig2k[,gg+1]
          delta.g <- post.delta[gg+1]
          zeta.g <- post.zeta[gg]
          lambda <- post.lambda
          zeta_i.min1 <- chain[[6]][ii-1,,gg]
          ones <- rep(1, K)
          varinv <- diag(1/(sig2k.g[1:K]))
          var.g <- (t(ones) %*% varinv) %*% ones

          # Hyper parameters
          SS <- b.zeta_i + sum((zeta_i.min1 - zeta.g)^2) + (N.g*n0.zeta_i*zeta.g)/(2*(N.g + n0.zeta_i))
          var0.zeta <- 1 / rgamma(1,(N.g + a.zeta_i)/2, SS/2)
          mu0.zeta <- rnorm(1, (N.g*var0.zeta*zeta.g)/(var0.zeta*(N.g + n0.zeta_i)), sqrt(1/(var0.zeta*(N.g + n0.zeta_i))))

          var.zeta_i <- 1/(N.g*(var.g) + 1/var0.zeta)
          mu.zeta_i <- rep(NA, N.g)
          for(zi in 1:N.g) {
            mu.zeta_i[zi] <- var.zeta_i*((N.g*(mean(lambda) - mean(RT.group[[gg]][zi,])))*(var.g) + mu0.zeta/var0.zeta)
          }
          # Draw N person speed parameters
          zeta_i <- rnorm(N.g, mu.zeta_i, sqrt(var.zeta_i))
          chain[[6]][ii,,gg] <- zeta_i
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
  post.theta_i <- colMeans((chain.1[[5]][XG.burnin:XG,,] + chain.2[[5]][XG.burnin:XG,,]) / 2)
  post.zeta_i <- colMeans((chain.1[[6]][XG.burnin:XG,,] + chain.2[[6]][XG.burnin:XG,,]) / 2)

  ### SAT ###
  nu.obs <- numeric(N)
  for (ii in 1:N) {
    select.RT <- RT[-ii, 2]
    select.Z <- Z[-ii, 2]
    nu.obs[ii] <- cov(select.RT, select.Z)
  }

  nu.prop <- numeric(N)
  mu.RT <- mean(post.lambda) - mean(post.zeta)
  mu.Z <- mean(post.theta) - mean(post.beta)
  for (dd in 1:N) {
    mu.t1 <- sqrt(K) * mu.RT
    var.t1 <- post.sig2 + K*post.delta
    draw.t1 <- rnorm(N, mu.t1, sqrt(var.t1))

    mu.z1 <- sqrt(K) * mu.Z
    var.z1 <- 1 + K*post.tau
    draw.z1 <- rnorm(N, mu.z1, sqrt(var.z1))

    nu.prop[dd] <- cov(draw.t1, draw.z1)
  }

  return(list(beta = post.beta, lambda = post.lambda, theta = post.theta, zeta = post.zeta, sig2k = post.sig2k, sig2 = post.sig2,
              tau = post.tau, delta = post.delta, theta_i = post.theta_i, zeta_i = post.zeta_i, nu = post.nu, nu.obs = nu.obs, nu.prop = draw.nu2, Z = Z))
}
