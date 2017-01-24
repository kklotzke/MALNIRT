#' @importFrom gmp asNumeric
#' @export
MALNRT <- function(RT, Group = NULL, data, XG = 1000, burnin = 0.10, inits.1 = NULL, inits.2 = NULL, hyper.priors = NULL, est.person = TRUE, silent = FALSE) {

  ###### Initialize ######

  if (!missing(data)) {
    # Try to find RT and Group in the data set first
    tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
    if (!is.null(Group)) {
      tryCatch(Group <- eval(substitute(Group), data), error=function(e) NULL)
    }
  } else {
    data = NULL
  }

  N <- nrow(RT) # Number of persons
  K <- ncol(RT) # Number of items

  # Multiple groups
  if (!is.null(Group)){
    groups <- unique(Group)
    G <- length(groups)
    RT.group <- vector("list", G)
    for (i in 1:G)
    {
      RT.group[[i]] <- subset(RT, Group == groups[i])
      RT.group[[i]] <- RT.group[[i]]
    }
  }
  else {
    G <- 1
    RT.group <- vector("list", 1)
    RT.group[[1]] <- RT
  }

  # Read hyperpriors
  if (is.null(hyper.priors)) {
    hyper.priors <- c(0.5, 0.5, 1, 0.1, 0.1, 1, 0.001, 0.001)
  }
  a.lambda <- hyper.priors[1]
  b.lambda <- hyper.priors[2]
  n0.lambda <- hyper.priors[3]
  a.zeta <- hyper.priors[4]
  b.zeta <- hyper.priors[5]
  n0.zeta <- hyper.priors[6]
  a.sig2 <- hyper.priors[7]
  b.sig2 <- hyper.priors[8]

  # Initialize chains
  chains.list <- vector("list", 2) # Two chains
  chains.list[[1]] <- vector("list", 6)
  chains.list[[1]][[1]] <- matrix(NA, nrow = XG, ncol = K) # Time intensity
  chains.list[[1]][[2]] <- matrix(NA, nrow = XG, ncol = G) # Speed group mean
  chains.list[[1]][[3]] <- array(NA, dim = c(XG, K, G+1)) # Measurement variance per item, first block across all groups
  chains.list[[1]][[4]] <- array(NA, dim = c(XG, 1, G+1)) # Average measurement variance, first block across all groups
  chains.list[[1]][[5]] <- array(NA, dim = c(XG, 1, G+1)) # RT covariance per group, first block across all groups
  chains.list[[1]][[6]] <- array(NA, dim = c(XG, N, G)) # Person speed parameter
  chains.list[[2]] <- chains.list[[1]]
  if (is.null(inits.1)) {
    inits.1 <- vector("list", 6)
    inits.1[[1]] <- rnorm(K, 10, 5) # Time intensity
    inits.1[[2]] <- rnorm(G, 10, 5) # Speed group mean
    inits.1[[3]] <- runif(K, 0.5, 1.5) # Measurement variance per item
    inits.1[[4]] <- mean(inits.1[[3]]) # Average measurement variance
    inits.1[[5]] <- runif(1, 0, 1.5) # RT covariance
    inits.1[[6]] <- rnorm(N, 10, 5) # Person speed parameter
  }
  if (is.null(inits.2)) {
    inits.2 <- vector("list", 6)
    inits.2[[1]] <- rnorm(K, 10, 5) # Time intensity
    inits.2[[2]] <- rnorm(G, 10, 5) # Speed group mean
    inits.2[[3]] <- runif(K, 0.5, 1.5) # Measurement variance per item
    inits.2[[4]] <- mean(inits.1[[3]]) # Average measurement variance
    inits.2[[5]] <- runif(1, 0, 1.5) # RT covariance
    inits.2[[6]] <- rnorm(N, 10, 5) # Person speed parameter
  }
  # First row are the initial values
  for (i in 1:2) {
    chains.list[[1]][[i]][1, ] <- inits.1[[i]]
    chains.list[[2]][[i]][1, ] <- inits.2[[i]]
  }
  for (i in 3:5) {
    for (g in 1:(G+1)) {
      chains.list[[1]][[i]][1,,g] <- inits.1[[i]]
      chains.list[[2]][[i]][1,,g] <- inits.2[[i]]
    }
  }
  for (g in 1:(G)) {
    chains.list[[1]][[6]][1,,g] <- inits.1[[6]]
    chains.list[[2]][[6]][1,,g] <- inits.2[[6]]
  }


  # Create helmert matrix
  hmat <- helmert(K)

  # Fractional Bayes factor
  m0 <- m1 <- FBF <- matrix(NA, nrow = XG-1, ncol = 2)


  ###### Run MCMC ######
  # For each chain
  for (cc in 1:2) {
    chain <- chains.list[[cc]] # Read the cc-th chain

    # For each iteration
    for (ii in 2:XG) {

      # Read values from former chain (t-1)
      lambda.min1 <- chain[[1]][ii-1, ]
      zeta.min1 <- chain[[2]][ii-1, ]
      sig2k.min1 <- chain[[3]][ii-1,,1]
      sig2.min1 <- chain[[4]][ii-1,,1]
      delta.min1 <- chain[[5]][ii-1,,1]
      mean.zeta <- mean(chain[[2]][1:(ii-1),])

      # Generalize sigma^2 and delta
      ones <- rep(1, K)
      #varinv <- matrix(1/delta.min1, nrow = K, ncol = K)
      #diag(varinv) <- 1/(sig2k.min1[1:K] + delta.min1)
      varinv <- diag(1/(sig2k.min1[1:K] + delta.min1))
      var.gen <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          sig2k.min1.g <- chain[[3]][ii-1,,gg+1]
          delta.min1.g <- chain[[5]][ii-1,,gg+1]
          varinv <- diag(1/(sig2k.min1.g[1:K]))# + delta.min1.g))
          tmp <- (t(ones) %*% varinv) %*% ones
          tmp <- tmp + delta.min1.g
          var.gen <- c(var.gen, tmp)
        }
      }
      #var.gen <- 1/(sig2.min1 + K*delta.min1)


      ### Sample speed parameters group means ###

      # Hyper parameters
      SS <- b.zeta + sum((zeta.min1 - mean(zeta.min1))^2) + (G*n0.zeta*mean(zeta.min1))/(2*(G + n0.zeta))
      var0.zeta <- 1 / rgamma(1, (G + a.zeta)/2, SS/2)
      mu0.zeta <- rnorm(1, (G*var0.zeta*mean(zeta.min1))/(var0.zeta*(G + n0.zeta)), sqrt(1/(var0.zeta*(G + n0.zeta))))
      #flat prior can be approximated by vague normal density prior with mean = 0 and variance = 10^10
      var0.zeta <- 10^10
      mu0.zeta <- 0

      var.zeta <- rep(NA, G)
      mu.zeta <- rep(NA, G)
      var.zeta[1] <- 0
      mu.zeta[1] <- 0
      if(G > 1) {
        for(gg in 2:G) {
          N.g <- nrow(RT.group[[gg]])
          var.zeta[gg] <- 1/(N.g/(var.gen[gg+1]) + 1/var0.zeta)
          mu.zeta[gg] <- var.zeta[gg]*((N.g*(mean(lambda.min1) - mean(RT.group[[gg]])))/(var.gen[gg+1]) + mu0.zeta/var0.zeta)
        }
      }

      # Draw G speed parameter group means
      zeta <- rnorm(G, mu.zeta, sqrt(var.zeta))
      chain[[2]][ii, ] <- zeta


      ### Sample item time intensity paramaters ###

      # Hyper parameters
      SS <- b.lambda + sum((lambda.min1 - mean(lambda.min1))^2) + (K*n0.lambda*mean(lambda.min1))/(2*(K + n0.lambda))
      var0.lambda <- 1 / rgamma(1, (K + a.lambda)/2, SS/2)
      mu0.lambda <- rnorm(1, (K*var0.lambda*mean(lambda.min1))/(var0.lambda*(K + n0.lambda)), sqrt(1/(var0.lambda*(K + n0.lambda))))

      #var.lambda <- 1/(N/(var.gen[1]) + 1/var0.lambda)
      #mu.lambda <- var.lambda*((N*(colMeans(RT) + mean(zeta)))/(var.gen[1]) + mu0.lambda/var0.lambda)

      var.lambda <- 1/(N*(var.gen[1]) + 1/var0.lambda)
      mu.lambda <- var.lambda*((N*(colMeans(RT) + mean(zeta)))*(var.gen[1]) + mu0.lambda/var0.lambda)

      # Draw K time intensity parameters
      chain[[1]][ii, ] <- lambda <- rnorm(K, mu.lambda, sqrt(var.lambda))


      ## Across all groups

      ### Sample measurement error variances

      # Within sum of squares
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta)
      mean.person <- apply(errors,1,mean)
      SSw <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),1,sum)) / K
      SSwk <- apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),2,sum)

      # Draw total measurement error variance
      #chain[[4]][ii,,1] <- sig2 <- (SSw/(K-1))/rgamma(1,shape=N/2,rate=1/2)
      chain[[4]][ii,,1] <- sig2 <- 1 / rgamma(1, (N+a.sig2)/2, (SSw+b.sig2)/2)

      # Draw K individual measurement error variances
      chain[[3]][ii,,1] <- sig2k <- 1 / rgamma(K, (N+a.sig2)/2, (SSwk+b.sig2)/2)


      ### Sample covariance parameter

      # Helmert transformation
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta) #+ matrix(zeta_i, nrow = N, ncol = K)
      #errors <- errors - mean(errors)
      tmat <- errors %*% t(hmat)

      # Between sum of squares
      mean.person <- apply(errors,1,mean)
      #print(cor(mean.person, zeta_i))
      #SSb1 <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K
      SSb <- sum((mean.person - mean(errors))^2)
      #cat(SSb, ",", SSb1, "\n")
      #SSb <- sum((mean.person - (mean(lambda) - mean(zeta_i)))^2)
      #mean.RT <- apply(RT,1,mean)
      #SSb <- sum((mean.RT - (mean(lambda) - zeta_i))^2)


      # Draw covariance parameter
      #chain[[5]][ii,,1] <- delta <- (SSb)/rgamma(1,shape=N/2,rate=1/2) - sig2/K
      chain[[5]][ii,,1] <- delta <- 1 / rgamma(1, N/2, SSb/2) - sig2/K

      ## Fractional Bayes Factor
      # delta != 0 vs delta = 0

      errors <- RT
      tmat <- errors %*% t(hmat)
      mean.person <- apply(errors,1,mean)
      SSb <- sum((mean.person - mean(errors))^2)
      SSw <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),1,sum)) / K
      SSwk <- apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),2,sum)

      m0[ii-1, cc] <- gmp::asNumeric((sig2/K)^((1-N)/2) * exp(((1-N)/(2*N)) * (sum(SSwk[2:K]/sig2k[2:K]) + SSb/(sig2/K))))
      m1[ii-1, cc] <- gmp::asNumeric(exp(-(1/2) * sum(SSwk[2:K]/sig2k[2:K])) * (gamma(N/2) * (SSb/2)^(-(N/2))) / (gamma(1/2) * (SSb/(2*N))^(-(1/2))))

      #FBF[ii-1, cc] <- gmp::asNumeric(m0[ii-1, cc] / m1[ii-1, cc])
      FBF[ii-1, cc] <- gmp::asNumeric(log(m0[ii-1, cc]) - log(m1[ii-1, cc]))

      ##

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
        SSw <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N.g))*(errors - matrix(mean.person,ncol=K,nrow=N.g)),1,sum)) / K
        SSwk <- apply((errors - matrix(mean.person,ncol=K,nrow=N.g))*(errors - matrix(mean.person,ncol=K,nrow=N.g)),2,sum)

        # Draw total measurement error variance
        chain[[4]][ii,,gg+1] <- sig2 <- 1 / rgamma(1, (N.g+a.sig2)/2, (SSw+b.sig2)/2)

        # Draw K individual measurement error variances
        chain[[3]][ii,,gg+1] <- sig2k <- 1 / rgamma(K, (N.g+a.sig2)/2, (SSwk+b.sig2)/2)

        # Draw covariance parameter
        chain[[5]][ii,,gg+1] <- delta <- 1 / rgamma(1, N.g/2, SSb/2) - sig2/K
      }


      if ((!silent) && (ii%%100 == 0))
        cat("Iteration ", ii, " ", "\n")
      flush.console()
    }

    chains.list[[cc]] <- chain
  }

  # Number of burnin iterations
  XG.burnin <- floor(XG * burnin)

  # Merge chains
  chain.1 <- chains.list[[1]]
  chain.2 <- chains.list[[2]]
  post.lambda <- colMeans((chain.1[[1]][XG.burnin:XG, ] + chain.2[[1]][XG.burnin:XG, ]) / 2)
  if(G > 1) {
    post.zeta <- colMeans((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
  }
  else {
    post.zeta <- mean((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
  }
  post.sig2k <- colMeans((chain.1[[3]][XG.burnin:XG,,] + chain.2[[3]][XG.burnin:XG,,]) / 2)
  post.sig2 <- colMeans((chain.1[[4]][XG.burnin:XG,,] + chain.2[[4]][XG.burnin:XG,,]) / 2)
  post.delta <- colMeans((chain.1[[5]][XG.burnin:XG,,] + chain.2[[5]][XG.burnin:XG,,]) / 2)

  #browser()
  # FBF
  raw.m0 <- (m0[XG.burnin:(XG-1), 1] + m0[XG.burnin:(XG-1), 2]) / 2
  raw.m1 <- (m1[XG.burnin:(XG-1), 1] + m1[XG.burnin:(XG-1), 2]) / 2
  raw.FBF <- (FBF[XG.burnin:(XG-1), 1] + FBF[XG.burnin:(XG-1), 2]) / 2
  valid.m0 <- !is.nan(raw.m0)
  valid.m1 <- !is.nan(raw.m1)
  valid.FBF <- valid.m0 * valid.m1
  m.m0 <- mean(raw.m0[valid.FBF])
  m.m1 <- mean(raw.m1[valid.FBF])
  m.FBF <- mean(raw.FBF[valid.FBF])
  FBF.perc.valid <- mean(valid.FBF)

  ### Sample person speed parameter

  if(est.person)
  {
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
          N.g <- nrow(RT.group[[gg]])
          sig2k.g <- post.sig2k[,gg+1] #chain[[3]][ii,,gg+1]
          delta.g <- post.delta[gg+1]
          zeta.g <- post.zeta[gg]
          lambda <- post.lambda
          zeta_i.min1 <- chain[[6]][ii-1,,gg]
          ones <- rep(1, K)
          #varinv <- diag(1/(sig2k.g[1:K] + delta.g))
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
  post.zeta_i <- colMeans((chain.1[[6]][XG.burnin:XG,,] + chain.2[[6]][XG.burnin:XG,,]) / 2)

  return(list(lambda = post.lambda, zeta = post.zeta, sig2k = post.sig2k, sig2 = post.sig2, delta = post.delta, zeta_i = post.zeta_i,
         m0 = m.m0, m1 = m.m1, FBF = m.FBF, FBF.perc.valid = FBF.perc.valid))
}
