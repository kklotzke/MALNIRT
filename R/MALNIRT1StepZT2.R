#' @importFrom truncnorm rtruncnorm
#' @importFrom invgamma dinvgamma
#' @export
MALNIRT.1StepZT2 <- function(Y, RT, Group = NULL, data, XG = 1000, burnin = 0.10, inits.1 = NULL, inits.2 = NULL, hyper.priors = NULL, est.person = FALSE, silent = FALSE) {

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
  chains.list[[1]][[13]] <- matrix(0, ncol = K, nrow = N) # Z|RT matrix
  chains.list[[1]][[14]] <- vector("list", G) # List to contain Z|RT matrices per group

  # Conditional response times
  chains.list[[1]][[15]] <- matrix(0, ncol = K, nrow = N) # RT|Z matrix
  chains.list[[1]][[16]] <- vector("list", G) # List to contain RT|Z matrices per group

  # SAT
  chains.list[[1]][[17]] <- matrix(NA, nrow = XG, ncol = K) # nu

  chains.list[[2]] <- chains.list[[1]]
  if (is.null(inits.1)) {
    inits.1 <- vector("list", 10)
    inits.1[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.1[[2]] <- rnorm(K, 5, 1) # Time intensity
    inits.1[[3]] <- rnorm(G, 0, 1) # Abiliy group means
    inits.1[[4]] <- rnorm(G, 10, 5) # Speed group means
    inits.1[[5]] <- rnorm(N, 0, 1) # Person ability parameter
    inits.1[[6]] <- rnorm(N, 10, 5) # Person speed parameter
    inits.1[[7]] <- runif(K, 0.5, 1.5) # Measurement variance for RT per item
    inits.1[[8]] <- mean(inits.1[[3]]) # Average measurement variance for RT
    inits.1[[9]] <- runif(1, 0, 1) # Response covariance
    inits.1[[10]] <- runif(1, 0, 1) # RT covariance
  }
  if (is.null(inits.2)) {
    inits.2 <- vector("list", 10)
    inits.2[[1]] <- rnorm(K, 0, 1) # Item difficulty
    inits.2[[2]] <- rnorm(K, 5, 1) # Time intensity
    inits.2[[3]] <- rnorm(G, 0, 1) # Abiliy group means
    inits.2[[4]] <- rnorm(G, 10, 5) # Speed group means
    inits.2[[5]] <- rnorm(N, 0, 1) # Person ability parameter
    inits.2[[6]] <- rnorm(N, 10, 5) # Person speed parameter
    inits.2[[7]] <- runif(K, 0.5, 1.5) # Measurement variance for RT per item
    inits.2[[8]] <- mean(inits.1[[3]]) # Average measurement variance for RT
    inits.2[[9]] <- runif(1, 0, 1) # Response covariance
    inits.2[[10]] <- runif(1, 0, 1) # RT covariance
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
    chains.list[[1]][[14]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
    chains.list[[2]][[14]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
    chains.list[[1]][[16]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))
    chains.list[[2]][[16]][[gg]] <- matrix(0, ncol = K, nrow = nrow(Y.group[[gg]]))

  }

  chains.list[[1]][[17]][1, ] <- -0.3
  chains.list[[2]][[17]][1, ] <- -0.3

  # Create helmert matrix
  hmat <- helmert(K)
  tmp1 <- hmat[K,1]
  tmp2 <- hmat[K,K]
  hmat2 <- matrix(tmp1, nrow = K, ncol = K)
  diag(hmat2) <- rep(tmp2, K)

  # Fractional Bayes factor
  m0 <- m1 <- FBF <- matrix(NA, nrow = XG-1, ncol = 2)

  ###### Run MCMC ######
  # For each chain
  for (cc in 1:2) {
    chain <- chains.list[[cc]] # Read the cc-th chain

    #ZT <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Z|T
    #ZT.cand <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Z|T
    ZT2 <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Zk|Z.mink, T1..p
    ZT2.cand <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Zk|Z.mink, T1..p
    mu_ZT2 <- matrix(0, ncol = K, nrow = N)
    mu_ZT <- mu_ZT.cand <- matrix(0, ncol = K, nrow = N)
    Sigma_ZT <- matrix(1, ncol = K, nrow = K)
    Z <- matrix(0, ncol = K, nrow = N)
    tau <- tau.cand <- 0
    delta <- delta.cand <- 0
    nu.cand <- numeric(K)
    sig2k <- sig2k.cand <- runif(K, 0.5, 1.5)
    sig2 <- sig2.cand <- mean(sig2k.cand)

    mh.accept <- 0

    mu.b <- var.b <- numeric(K)
    SSb.tau <- runif(1, 1, 100)
    SSb.delta <- runif(1, 1, 100)
    SSw.sig2 <- runif(1, 1, 100)
    SSwk.sig2k <- runif(K, 1, 100)

    #Z <- chain[[11]]
    #Z.group <- chain[[12]]
    #Z.RT <- chain[[13]]
    #Z.RT.group <- chain[[14]]
    #RT.Z <- chain[[15]]
    #RT.Z.group <- chain[[16]]

    # For each iteration
    ii <- 2
    while (ii <= XG) {
      #for (ii in 2:XG) {

      # Read values from former chain (t-1)
      beta.s0 <- chain[[1]][ii-1, ]
      lambda.s0 <- chain[[2]][ii-1, ]
      theta.s0 <- chain[[3]][ii-1, ]
      zeta.s0 <- chain[[4]][ii-1, ]
      sig2k.s0 <- chain[[7]][ii-1,,1]
      sig2.s0 <- chain[[8]][ii-1,,1]
      tau.s0 <- chain[[9]][ii-1,,]
      delta.s0 <- chain[[10]][ii-1,,]
      nu.s0 <- chain[[17]][ii-1, ]

      # Generalize sigma^2 (fixed to 1) and tau
      # ones <- rep(1, K)
      # varinv <- diag(1/(rep(1,K) + tau.s0[1]))
      # var.gen.Z <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          varinv <- diag(1/(rep(1,K) + tau.s0[gg+1]))
          tmp <- (t(ones) %*% varinv) %*% ones
          var.gen.Z <- c(var.gen.Z, tmp)
        }
      }

      # Generalize sigma^2 and delta
      ones <- rep(1, K)
      varinv <- diag(1/(sig2k.s0[1:K] + delta.s0[1]))
      var.gen.RT <- (t(ones) %*% varinv) %*% ones
      if (G > 1) {
        for (gg in 1:G) {
          sig2k.s0.g <- chain[[7]][ii-1,,gg+1]
          varinv <- diag(1/(sig2k.s0.g[1:K] + delta.s0[gg+1]))
          tmp <- (t(ones) %*% varinv) %*% ones
          var.gen.RT <- c(var.gen.RT, tmp)
        }
      }


      #nu.cand <- nu.s0

      ##### Responses #####

      ### Sample latent response variable ###

      ## Obtain latent respones from simulated Y's ##

      I <- diag(K)
      J <- matrix(1, nrow = K, ncol = K)
      a <- 1/delta.s0[1] + sum(1/sig2k.s0)
      a.cand <- 1/delta.cand + sum(1/sig2k.cand)

      # Means latent responses, response times
      mu_Z <- -beta.s0 # Assuming mu_theta = 0
      mu_T <- lambda.s0 # Assuming mu_zeta = 0

      # # Sigma_11, Sigma_22 and Sigma_12 in joint-model covariance matrix
      Sigma_Z <- I + tau.s0[1]*J
      # Sigma_T <- sig2k.s0*I + delta.s0[1]*J
      # Sigma_Z_T <- nu.s0*I
      # Sigma_Z.cand <- I + tau.cand*J
      # Sigma_T.cand <- sig2k.cand*I + delta.cand*J
      # Sigma_Z_T.cand <- nu.cand*I
      #
      # # Sigma_11.inv, Sigma_22.inv
      # Sigma_Z.inv <- I - J / (1/tau.s0[1] + K) #solve(Sigma_Z)
      # Sigma_T.inv <- (1/sig2k.s0)*I - ((1/sig2k.s0)%*%t(1/sig2k.s0)) / a #solve(Sigma_T)
      # Sigma_Z.inv.cand <- I - J / (1/tau.cand + K)
      # Sigma_T.inv.cand <- (1/sig2k.cand)*I - ((1/sig2k.cand)%*%t(1/sig2k.cand)) / a.cand
      #
      # # Sigma Z|T
      # Sigma_ZT <- Sigma_Z - Sigma_Z_T %*% Sigma_T.inv %*% Sigma_Z_T
      # Sigma_ZT.cand <- Sigma_Z.cand - Sigma_Z_T.cand %*% Sigma_T.inv.cand %*% Sigma_Z_T.cand
      #browser()
      #tmp.sigma <- Sigma_Z_T %*% Sigma_T.inv %*% Sigma_Z_T
      #print(Sigma_ZT + Sigma_Z_T %*% Sigma_T.inv %*% Sigma_Z_T)

      # Conditional mean of Z|T for each item
      #mu_ZT <- matrix(0, ncol = K, nrow = N)
      x <- x.cand <- matrix(NA, ncol = K, nrow = N)
      for (k in 1:K)
      {
        x[, k] <- (1/sig2k.s0[k]) * ((RT[, k] - mu_T[k]) - (RT[, k] - mu_T[k]) / (a * sig2k.s0[k]))
        x.cand[, k] <- (1/sig2k.cand[k]) * ((RT[, k] - mu_T[k]) - (RT[, k] - mu_T[k]) / (a.cand * sig2k.cand[k]))

        mu_ZT[, k] <- mu_Z[k] + nu.s0[k] * x[, k]
        mu_ZT.cand[, k] <- mu_Z[k] + nu.cand[k] * x.cand[, k]
      }

      # Sigma Z|T
      #Sigma_ZT <- Sigma_Z - Sigma_Z_T %*% Sigma_T.inv %*% Sigma_Z_T

      #lik <- sum(mvtnorm::dmvnorm(ZT, mean = mu_ZT[i, ], sigma = Sigma_ZT, log = TRUE))
      #lik.cand <- sum(mvtnorm::dmvnorm(ZT, mean = mu_ZT[i, ], sigma = Sigma_ZT.cand, log = TRUE))


      # # Sample Z|T
      #       for(i in 1:N)
      #       {
      #          ZT[i, ] <- mvtnorm::rmvnorm(1, mean = mu_ZT[i, ], sigma = Sigma_ZT)
      #          ZT.cand[i, ] <- mvtnorm::rmvnorm(1, mean = mu_ZT[i, ], sigma = Sigma_ZT.cand)
      #      }
      #       lik <- sum(mvtnorm::dmvnorm(ZT, mean = mu_ZT.s0[i, ], sigma = Sigma_ZT.s0, log = TRUE))
      #       lik.cand <- sum(mvtnorm::dmvnorm(ZT.cand, mean = mu_ZT.s0[i, ], sigma = Sigma_ZT.s0, log = TRUE))
      # #browser()
      #       u <- runif (1, 0, 1)
      #       ar <- exp(lik.cand - lik)
      #       if(is.nan(ar))
      #         ar <- 0
      #
      #       if (u <= ar) {
      #         cat("Accept tau: ", tau.cand, "\n")
      #         chain[[9]][ii,,1] <- tau <- tau.cand
      #       }
      #       else {
      #         #cat("Reject tau: ", tau.cand, "\n")
      #         chain[[9]][ii,,1] <- tau <- tau.s0[1]
      #       }
      #browser()
      #print(tau.s0[1])
      #print(tau.cand)



      # Conditional mean of Zk|Z.mink, T1..p for each item
      I.min1 <- diag(K-1)
      J.min1 <- matrix(1, nrow = K-1, ncol = K-1)
      ones.min1 <- rep(1, K-1)

      var_ZT2.cand <- var_ZT2.plus <- numeric(K)
      r <- matrix(NA, ncol = K, nrow = N)
      s <- matrix(NA, ncol = K, nrow = N)
      mu_ZT2.cand <- mu_ZT2.plus <- mu_ZT2 #matrix(0, ncol = K, nrow = N)
      ZT2.cand <- ZT2.plus <- ZT2
      lik_k <- lik_k.cand <- numeric(K)
      reset <- FALSE
      for (p in 1:2)
      {
        for (k in 1:K)
        {
          w.min1 <- nu.s0[-k]/sig2k.s0[-k]
          b.min1 <- sig2k.s0[-k] / (sig2k.s0[-k] - nu.s0[-k]^2)
          c.min1 <- 1/tau.s0[1] + sum(b.min1)
          A.min1.inv <- b.min1*I.min1 - (b.min1 %*% t(b.min1)) / c.min1
          d.min1 <- a + t(w.min1) %*% A.min1.inv %*% w.min1
          g.min1 <- sum(nu.s0[-k] / (sig2k.s0[-k] - nu.s0[-k]^2)) / c.min1
          A.min1.inv_w <- (nu.s0[-k] / (sig2k.s0[-k] - nu.s0[-k]^2) - g.min1*b.min1)

          B11 <- 1 + tau.s0[1] - nu.s0[k]^2 * (1/sig2k.s0[k] - ((1/sig2k.s0[k])^2) / a)
          B22 <- (1 - (nu.s0[-k]^2) / sig2k.s0[-k])*I.min1 + tau.s0[1]*J.min1 + w.min1 %*% t(w.min1) / a
          B12 <- tau.s0[1]*ones.min1 + ((nu.s0[k] / sig2k.s0[k]) %*% t(w.min1)) / a
          B21 <- t(B12)
          B22.inv <- A.min1.inv - (A.min1.inv_w %*% t(A.min1.inv_w)) / d.min1[1,1]

          #tmp <- B12 %*% B22.inv %*% t(ZT[, -k] - mu_ZT[, -k])
          tmp <- B12 %*% B22.inv %*% t(ZT2.plus[, -k] - mu_ZT2.plus[, -k])
          mu_ZT2.plus[, k] <-  mu_ZT[, k] + tmp
          var_ZT2.plus[k] <- B11 - B12 %*% B22.inv %*% B21


          ### MH candidates ###
          w.min1.cand <- nu.cand[-k]/sig2k.cand[-k]
          b.min1.cand <- sig2k.cand[-k] / (sig2k.cand[-k] - nu.cand[-k]^2)
          c.min1.cand <- 1/tau.cand + sum(b.min1.cand)
          A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
          d.min1.cand <- a.cand + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
          g.min1.cand <- sum(nu.cand[-k] / (sig2k.cand[-k] - nu.cand[-k]^2)) / c.min1.cand
          A.min1.inv_w.cand <- (nu.cand[-k] / (sig2k.cand[-k] - nu.cand[-k]^2) - g.min1.cand*b.min1.cand)

          B11.cand <- 1 + tau.cand - nu.cand[k]^2 * (1/sig2k.cand[k] - ((1/sig2k.cand[k])^2) / a.cand)
          B22.cand <- (1 - (nu.cand[-k]^2) / sig2k.cand[-k])*I.min1 + tau.cand*J.min1 + w.min1.cand %*% t(w.min1.cand) / a.cand
          B12.cand <- tau.cand*ones.min1 + ((nu.cand[k] / sig2k.cand[k]) %*% t(w.min1.cand)) / a.cand
          B21.cand <- t(B12.cand)
          B22.inv.cand <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]

          #tmp <- B12 %*% B22.inv %*% t(ZT[, -k] - mu_ZT[, -k])
          tmp.cand <- B12.cand %*% B22.inv.cand %*% t(ZT2.cand[, -k] - mu_ZT2.cand[, -k])
          mu_ZT2.cand[, k] <-  mu_ZT.cand[, k] + tmp.cand
          var_ZT2.cand[k] <- B11.cand - B12.cand %*% B22.inv.cand %*% B21.cand
          #######################

          if(!is.nan(var_ZT2.plus[k]) && !is.nan(var_ZT2.cand[k])) {
            if((var_ZT2.plus[k] < 0) || (var_ZT2.cand[k] < 0))
              reset <- TRUE
          }
          else if (is.nan(var_ZT2.plus[k]) || is.nan(var_ZT2.cand[k]))
            reset <- TRUE


          if(!reset) {
          # # Sample Zk|Z.mink, T1..p
          ZT2.plus[Y[,k]==0, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2[Y[, k]==0, k]), mean = mu_ZT2.plus[Y[, k]==0, k], sd = sqrt(var_ZT2.plus[k]), a = -Inf, b = 0)
          ZT2.plus[Y[,k]==1, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2[Y[, k]==1, k]), mean = mu_ZT2.plus[Y[, k]==1, k], sd = sqrt(var_ZT2.plus[k]), a = 0, b = Inf)
          ZT2.cand[Y[,k]==0, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2.cand[Y[, k]==0, k]), mean = mu_ZT2.cand[Y[, k]==0, k], sd = sqrt(var_ZT2.cand[k]), a = -Inf, b = 0)
          ZT2.cand[Y[,k]==1, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2.cand[Y[, k]==1, k]), mean = mu_ZT2.cand[Y[, k]==1, k], sd = sqrt(var_ZT2.cand[k]), a = 0, b = Inf)
}
          # browser()
          if(p == 2) {
            lik0 <- sum(log(truncnorm::dtruncnorm(x = ZT2[Y[,k]==0, k], mean = mu_ZT2.plus[Y[, k]==0, k], sd = sqrt(var_ZT2.plus[k]), a = -Inf, b = 0)))
            lik1 <- sum(log(truncnorm::dtruncnorm(x = ZT2[Y[,k]==1, k], mean = mu_ZT2.plus[Y[, k]==1, k], sd = sqrt(var_ZT2.plus[k]), a = 0, b = Inf)))
            lik_k[k] <- lik0 + lik1 + dnorm(nu.cand[k], mu.b[k], sqrt(var.b[k]), TRUE) + invgamma::dinvgamma(sig2k.cand[k], a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk.sig2k[k]/2, log = TRUE)

            # lik0.cand <- sum(log(truncnorm::dtruncnorm(x = ZT2.cand[Y[,k]==0, k], mean = mu_ZT2.cand[Y[, k]==0, k], sd = sqrt(var_ZT2.cand[k]), a = -Inf, b = 0)))
            # lik1.cand <- sum(log(truncnorm::dtruncnorm(x = ZT2.cand[Y[,k]==1, k], mean = mu_ZT2.cand[Y[, k]==1, k], sd = sqrt(var_ZT2.cand[k]), a = 0, b = Inf)))
            lik0.cand <- sum(log(truncnorm::dtruncnorm(x = ZT2[Y[,k]==0, k], mean = mu_ZT2.cand[Y[, k]==0, k], sd = sqrt(var_ZT2.cand[k]), a = -Inf, b = 0)))
            lik1.cand <- sum(log(truncnorm::dtruncnorm(x = ZT2[Y[,k]==1, k], mean = mu_ZT2.cand[Y[, k]==1, k], sd = sqrt(var_ZT2.cand[k]), a = 0, b = Inf)))
            lik_k.cand[k] <- lik0.cand + lik1.cand + dnorm(nu.s0[k], mu.b[k], sqrt(var.b[k]), TRUE) + invgamma::dinvgamma(sig2k.s0[k], a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk.sig2k[k]/2, log = TRUE)
          }
          # Sample Zk|Z.mink, T1..p
          #ZT2[Y[,k]==0, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2[Y[, k]==0, k]), mean = mu_ZT2[Y[, k]==0, k], sd = sqrt(var_ZT2[k]), a = -Inf, b = 0)
          #ZT2[Y[,k]==1, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2[Y[, k]==1, k]), mean = mu_ZT2[Y[, k]==1, k], sd = sqrt(var_ZT2[k]), a = 0, b = Inf)
          #ZT2.cand[Y[,k]==0, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2.cand[Y[, k]==0, k]), mean = mu_ZT2.cand[Y[, k]==0, k], sd = sqrt(var_ZT2.cand[k]), a = -Inf, b = 0)
          #ZT2.cand[Y[,k]==1, k] <- truncnorm::rtruncnorm(n = length(mu_ZT2.cand[Y[, k]==1, k]), mean = mu_ZT2.cand[Y[, k]==1, k], sd = sqrt(var_ZT2.cand[k]), a = 0, b = Inf)


          r[, k] <- (ZT2.plus[, k] - mu_Z[k]) - tmp[1, ]
          #s[, k] <- ZT2[, k] - tmp[1, ] - nu.s0[k]*x[,k]
        }
      }

      #if(anyNA(ZT2.plus) || anyNA(ZT2.cand))
      #  browser()
      #cat(sum(lik_k.cand), " ", sum(lik_k), "\n")
      #print(sum(lik_k))

      sum.lik <- sum(lik_k) + invgamma::dinvgamma(tau.cand + 1/K, N/2, SSb.tau/2, log = TRUE) + invgamma::dinvgamma(delta.cand + sig2.cand/K, N/2, SSb.delta/2, log = TRUE)
      sum.lik.cand <- sum(lik_k.cand) + invgamma::dinvgamma(tau.s0[1] + 1/K, N/2, SSb.tau/2, log = TRUE) + invgamma::dinvgamma(delta.s0[1] + sig2.cand/K, N/2, SSb.delta/2, log = TRUE)

      #if(is.nan(sum.lik) || is.nan(sum.lik.cand))
      #  browser()

      u <- runif (1, 0, 1)
      ar <- exp(sum.lik.cand - sum.lik)
      if(is.nan(ar))
        ar <- 0

      #if(ii <= 20)
      #  ar <- 1

      #
      # if (!is.nan(sum.lik.cand) && !is.nan(sum.lik)) {
      #   if(((sum.lik.cand == -Inf) && (sum.lik == -Inf)) || (sum.lik.cand == Inf))
      #     ar <- 1
      #   else
      #     ar <- exp(sum.lik.cand - sum.lik)
      # }
      # else if (!is.nan(sum.lik.cand) && is.nan(sum.lik)) {
      #   ar <- 1
      # }
      # else {
      #   ar <- 0
      # }
      #
      # if(is.nan(ar)) {
      #   browser()
      #   ar <- 0
      # }

      if (u <= ar) {
        #cat("Accept tau|delta|nu: ", tau.cand, "|", delta.cand, "|", nu.cand, "\n")
        chain[[7]][ii,,1] <- sig2k <- sig2k.cand
        chain[[8]][ii,,1] <- sig2 <- sig2.cand
        chain[[9]][ii,,1] <- tau <- tau.cand
        chain[[10]][ii,,1] <- delta <- delta.cand
        chain[[17]][ii, ] <- nu <- nu.cand
        mu_ZT2 <- mu_ZT2.cand
        var_ZT2 <- var_ZT2.cand
        ZT2 <- ZT2.cand
        mh.accept <- mh.accept + 1
      }
      else {
        #print("Reject tau")
        #cat("Reject tau|delta|nu: ", tau.cand, "|", delta.cand, "|", nu.cand, "\n")
        chain[[7]][ii,,1] <- sig2k <- sig2k.s0
        chain[[8]][ii,,1] <- sig2 <- sig2.s0
        chain[[9]][ii,,1] <- tau <- tau.s0[1]
        chain[[10]][ii,,1] <- delta <- delta.s0[1]
        chain[[17]][ii, ] <- nu <- nu.s0
        mu_ZT2 <- mu_ZT2.plus
        var_ZT2 <- var_ZT2.plus
        ZT2 <- ZT2.plus
      }


      ### Sample covariance responses-response times ###
      var0.b <- 10^10
      mu0.b <- 0

      tmp1.r <- tmp2.r <- rep(0, K)
      b <- var.b <- mu.b <- numeric(K)
      for (k in 1:K) {
        for (i in 1:N) {
          tmp1.r[k] <- tmp1.r[k] + x[i, k] * 1/var_ZT2[k] * x[i, k]
          tmp2.r[k] <- tmp2.r[k] + x[i, k] * 1/var_ZT2[k] * r[i, k] #ZT2[i, k] #
        }

        var.b[k] <- (1/(1/var0.b + tmp1.r[k]))
        mu.b[k] <- (var.b[k] * (mu0.b/var0.b + tmp2.r[k]))
        b[k] <- rnorm(1, mean = mu.b[k], sd = sqrt(var.b[k]))
      }

      nu.cand <- b #runif(K, nu.cand - 0.01, nu.cand + 0.01) #b
      #print(b)
      chain[[17]][ii, ] <- nu <- b


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
          mu.theta[gg] <- var.theta[gg]*((N.g*(mean(beta.s0) + mean(ZT2)))/(var.gen.Z[gg+1]) + mu0.theta/var0.theta)
        }
      }

      # Draw G speed parameter group means
      theta <- rnorm(G, mu.theta, sqrt(var.theta))
      chain[[3]][ii, ] <- theta


      ### Sample item difficulty paramaters ###

      ones <- rep(1, K)
      varinv <- diag(1/(rep(1,K) + tau))
      var.gen.Z <- (t(ones) %*% varinv) %*% ones

      # # Hyper parameters
      SS <- b.beta + sum((beta.s0 - mean(beta.s0))^2) + (K*n0.beta*mean(beta.s0))/(2*(K + n0.beta))
      var0.beta <- 1 / rgamma(1, (K + a.beta)/2, SS/2)
      mu0.beta <- rnorm(1, (-K*var0.beta*mean(beta.s0))/(var0.beta*(K + n0.beta)), sqrt(1/(var0.beta*(K + n0.beta))))

      #var0.beta <- 10^10
      #mu0.beta <- 0

      var.beta <- 1/(N*(var.gen.Z[1]) + 1/var0.beta)
      mu.beta <- var.beta*((N*(- colMeans(ZT2)))*(var.gen.Z[1]) + mu0.beta/var0.beta)

      # Draw K item difficulty parameters
      chain[[1]][ii, ] <- beta <- rnorm(K, mu.beta, sqrt(var.beta))
      # #print(beta)

      #
      # tmp1.s <- tmp2.s <- rep(0, K)
      # b2 <- var.b2 <- mu.b2 <- numeric(K)
      # for (k in 1:K) {
      #   for (i in 1:N) {
      #     #s <- r[i, k] - nu[k]*x[i, k] + mu_Z[k] # mu_Z = -beta
      #     tmp1.s[k] <- tmp1.s[k] + 1 * 1/var_ZT2[k] * 1
      #     tmp2.s[k] <- tmp2.s[k] + 1 * 1/var_ZT2[k] * s[i, k]
      #   }
      #
      #   var.b2[k] <- (1/(1/var0.beta + tmp1.s[k]))
      #   mu.b2[k] <- (var.b2[k] * (mu0.beta/var0.beta + tmp2.s[k]))
      #   b2[k] <- rnorm(1, mean = mu.b2[k], sd = sqrt(var.b2[k]))
      # }
      # chain[[1]][ii, ] <- beta <- -b2



      ### Sample covariance parameter

      # Generate proposal for tau from marginalized ability model
      mu0.beta <- 0
      var0.beta <- 10^10

      ones <- rep(1, K)
      varinv <- diag(1/(rep(1,K) + tau.cand))
      var.gen.mar <- (t(ones) %*% varinv) %*% ones
      var.beta.mar <- 1/(N*(var.gen.mar) + 1/var0.beta)
      mu.beta.mar <- var.beta.mar*((N*(- colMeans(Z)))*(var.gen.mar) + mu0.beta/var0.beta)
      beta.mar <- rnorm(K, mu.beta.mar, sqrt(var.beta.mar))

      Sjc <- matrix(tau.cand / (1 + (K - 1) * tau.cand), ncol = 1, nrow = K - 1)
      var.Z.mar <- (1 + K * tau.cand)/(1 + (K - 1) * tau.cand)
      theta1 <- matrix(0, ncol = K - 1, nrow = N)
      #beta <- data$beta
      for (kk in 1:K){
        beta1 <- beta.mar[-kk]
        Z1 <- Z[, -kk] # Latent responses to all but the current item
        mu.Z.mar <- (- beta.mar[kk]) + (Z1 - (theta1 - beta1)) %*% Sjc
        Z[Y[, kk]==0, kk] <- qnorm(runif(N, 0, pnorm(0, mu.Z.mar, sqrt(var.Z.mar))), mu.Z.mar, sqrt(var.Z.mar))[Y[, kk] == 0]
        Z[Y[, kk]==1, kk] <- qnorm(runif(N, pnorm(0, mu.Z.mar, sqrt(var.Z.mar)),1), mu.Z.mar, sqrt(var.Z.mar))[Y[, kk] == 1]
      }
      mean.person <- apply(Z,1,mean)
      SSb.tau <- sum((mean.person - mean(Z))^2)

      # Draw covariance parameter
      tau.cand <- 1 / rgamma(1, N/2, SSb.tau/2) - 1/K



      ## For each group
      for(gg in 1:G) {
        # N.g <- nrow(Z.group[[gg]])
        #
        # # Helmert transformation
        # errors <- Z.group[[gg]] + matrix(beta, nrow = N.g, ncol = K, byrow = TRUE) + theta[gg]
        # tmat <- errors %*% t(hmat)
        #
        # # Between sum of squares
        # SSb <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K

        # Draw covariance parameter
        chain[[9]][ii,,gg+1] <- data$tau #1 / rgamma(1, N.g/2, SSb/2) - 1/K
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
          mu.zeta[gg] <- var.zeta[gg]*((N.g*(mean(lambda.s0) - mean(RT.group[[gg]])))/(var.gen.RT[gg+1]) + mu0.zeta/var0.zeta)
        }
      }

      # Draw G speed parameter group means
      zeta <- rnorm(G, mu.zeta, sqrt(var.zeta))
      chain[[4]][ii, ] <- zeta


      ### Sample item time intensity paramaters ###
#browser()
      # Generalize sigma^2 and delta
      ones <- rep(1, K)
      varinv <- diag(1/(sig2k[1:K] + delta))
      var.gen.RT <- (t(ones) %*% varinv) %*% ones

      # Hyper parameters
      SS <- b.lambda + sum((lambda.s0 - mean(lambda.s0))^2) + (K*n0.lambda*mean(lambda.s0))/(2*(K + n0.lambda))
      var0.lambda <- 1 / rgamma(1, (K + a.lambda)/2, SS/2)
      mu0.lambda <- rnorm(1, (K*var0.lambda*mean(lambda.s0))/(var0.lambda*(K + n0.lambda)), sqrt(1/(var0.lambda*(K + n0.lambda))))

      var.lambda <- 1/(N*(var.gen.RT[1]) + 1/var0.lambda)
      mu.lambda <- var.lambda*((N*(colMeans(RT) + mean(zeta)))*(var.gen.RT[1]) + mu0.lambda/var0.lambda)

      # Draw K time intensity parameters
      chain[[2]][ii, ] <- lambda <- rnorm(K, mu.lambda, sqrt(var.lambda))


      ### Sample measurement error variances

#       ones <- rep(1, K)
#       varinv <- diag(1/(sig2k.cand[1:K] + delta.cand))
#       var.gen.RT <- (t(ones) %*% varinv) %*% ones
# #browser()
#       mu0.lambda <- 0
#       var0.lambda <- 10^10
#       var.lambda.mar <- 1/(N*(var.gen.RT[1]) + 1/var0.lambda)
#       mu.lambda.mar <- var.lambda.mar*(N*(colMeans(RT))*(var.gen.RT[1]) + mu0.lambda/var0.lambda)
#
#       # Draw K time intensity parameters
#       lambda.mar <- rnorm(K, mu.lambda.mar, sqrt(var.lambda.mar))


      # Within sum of squares
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta)
      # tmat <- errors %*% t(hmat2)
      mean.person <- apply(errors,1,mean)
      SSw.sig2 <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),1,sum))
      SSwk.sig2k <- apply((errors - matrix(mean.person,ncol=K,nrow=N))*(errors - matrix(mean.person,ncol=K,nrow=N)),2,sum)

      # Draw total measurement error variance
      #chain[[8]][ii,,1] <- sig2 <- mean(data$sig2k)#( 1 / rgamma(1, a.sig2 + (N*(K-1))/2, b.sig2 + SSw/2) )
      sig2.cand <- ( 1 / rgamma(1, a.sig2 + (N*(K-1))/2, b.sig2 + SSw.sig2/2) )

      # Draw K individual measurement error variances
      #chain[[7]][ii,,1] <- sig2k <- data$sig2k# (1 / rgamma(K, a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk/2))
      sig2k.cand <- (1 / rgamma(K, a.sig2 + (N*(K-1)/K)/2, b.sig2 + SSwk.sig2k/2))

      this.sig2k <- sig2k.cand #<- data.lnirt$sig2k.lnrt

#print(sig2k.cand[1:2])
      ### Sample covariance parameter

      # Helmert transformation
      errors <- RT - matrix(lambda, nrow = N, ncol = K, byrow = TRUE) #+ mean(zeta) #+ matrix(zeta_i, nrow = N, ncol = K)
      tmat <- errors %*% t(hmat)

      # Between sum of squares
      mean.person <- apply(errors,1,mean)
      SSb.delta <- sum((mean.person - mean(errors))^2)

      #print(mean(1 / rgamma(10000, N/2, SSb/2) - sig2/K))
      delta.cand <- 1 / rgamma(1, N/2, SSb.delta/2) - sig2.cand/K

      # Draw covariance parameter
      #chain[[10]][ii,,1] <- delta <- data$delta# 1 / rgamma(1, N/2, SSb/2) - sig2/K

      ## For each group
      for(gg in 1:G) {
        # N.g <- nrow(RT.group[[gg]])
        #
        # # Helmert transformation
        # errors <- RT.group[[gg]] - matrix(lambda, nrow = N.g, ncol = K, byrow = TRUE) + zeta[gg]
        # tmat <- errors %*% t(hmat)
        #
        # # Between sum of squares
        # SSb <- sum((tmat[, 1] - sqrt(K)*(mean(errors)))^2) / K
        #
        # # Within sum of squares
        # mean.person <- apply(errors,1,mean)
        # SSw <- sum(apply((errors - matrix(mean.person,ncol=K,nrow=N.g))*(errors - matrix(mean.person,ncol=K,nrow=N.g)),1,sum))
        # SSwk <- apply((errors - matrix(mean.person,ncol=K,nrow=N.g))*(errors - matrix(mean.person,ncol=K,nrow=N.g)),2,sum)

        # Draw total measurement error variance
        chain[[8]][ii,,gg+1] <- mean(data$sig2k) #1 / rgamma(1, a.sig2 + (N.g*(K-1))/2, b.sig2 + SSw/2)

        # Draw K individual measurement error variances
        chain[[7]][ii,,gg+1] <- data$sig2k #1 / rgamma(K, a.sig2 + ((N.g*(K-1)/K))/2, b.sig2 + SSwk/2)

        # Draw covariance parameter
        chain[[10]][ii,,gg+1] <- data$delta # 1 / rgamma(1, N.g/2, SSb/2) - sig2/K
      }





      if ((!silent) && (ii%%100 == 0))
        cat("Iteration ", ii, " | MH acceptance rate ", mh.accept/ii, "\n")

      if(((ii%%100 == 0) && ((mh.accept/ii) < 0.1)) || reset) {
        ZT2 <- matrix(beta.mar, nrow = N, ncol = K, byrow = TRUE) #+ rnorm(N*K, 0, 0.2) #data$Z + rnorm(N*K, 0, 0.5)# Z# <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Zk|Z.mink, T1..p
        #ZT2.cand <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N) # Zk|Z.mink, T1..p
        mu_ZT2 <- matrix(beta.mar, nrow = N, ncol = K, byrow = TRUE) #matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N)
        mu_ZT <- mu_ZT2 #+ rnorm(N*K, 0, 0.2)# <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N)
        Sigma_ZT <- Sigma_Z #toeplitz((runif(5, 1, 5))/5) #matrix(1, ncol = K, nrow = K)
        #Z <- matrix(rnorm(N*K, 0, 1), ncol = K, nrow = N)
        #tau.cand <- runif(1,0,1)
        #delta.cand <- runif(1,0,1)
        #nu.cand <- runif(K,-0.3,0.3)
        #SSb.tau <- runif(1, 1, 100)
        #SSb.delta <- runif(1, 1, 100)

        if(reset) {
          chain[[1]][ii-1, ] <- beta.mar #data$beta + rnorm(K, 0, 0.5) #rnorm(K, 0, 1)
          #chain[[2]][1, ] <- rnorm(K, 5, 1)
          chain[[7]][ii-1,,1] <- sig2k <- sig2k.cand
          chain[[8]][ii-1,,1] <- sig2 <- mean(sig2k) #sig2.cand
          chain[[9]][ii-1,,] <- rep(tau.cand, 2) #runif(2,0,1)
          chain[[10]][ii-1,,] <- rep(delta.cand, 2) #runif(2,0,1)
          chain[[17]][ii-1, ] <- diag(cov(RT, Z)) #nu.cand #data$nu + rnorm(K, 0, 0.3) #nu.cand #runif(K,-0.3,0.3)
          tau <- tau.cand
          delta <- delta.cand

          ii <- ii - 1
        }
        else {
          chain[[1]][1, ] <- beta.mar #data$beta + rnorm(K, 0, 0.5) #rnorm(K, 0, 1)
          #chain[[2]][1, ] <- rnorm(K, 5, 1)
          chain[[7]][1,,1] <- sig2k <- sig2k.cand
          chain[[8]][1,,1] <- sig2 <- sig2.cand
          chain[[9]][1,,] <- rep(tau.cand, 2) #runif(2,0,1)
          chain[[10]][1,,] <- rep(delta.cand, 2) #runif(2,0,1)
          chain[[17]][1, ] <- diag(cov(RT, Z)) # data$nu + rnorm(K, 0, 0.3) #nu.cand #runif(K,-0.3,0.3)
          tau <- tau.cand
          delta <- delta.cand

          mh.accept <- 0
          ii <- 1
        }



        print(reset)
      }

      ii <- ii + 1

      flush.console()
    }

    #chain[[11]] <- Z
    #chain[[12]] <- Z.group
    #chain[[13]] <- Z.RT
    #chain[[14]] <- Z.RT.group
    #chain[[15]] <- RT.Z
    #chain[[16]] <- RT.Z.group
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
  post.Z.RT <- (chain.1[[13]] + chain.2[[13]]) / 2
  post.Z.RT.group <- vector("list", G)
  #post.RT.Z <- (chain.1[[15]] + chain.2[[15]]) / 2
  #post.RT.Z.group <- vector("list", G)
  for (gg in 1:G) {
    post.Z.group[[gg]] <- (chain.1[[12]][[gg]] + chain.2[[12]][[gg]]) / 2
    post.Z.RT.group[[gg]] <- (chain.1[[14]][[gg]] + chain.2[[14]][[gg]]) / 2
    #post.RT.Z.group[[gg]] <- (chain.1[[16]][[gg]] + chain.2[[16]][[gg]]) / 2
  }
  post.nu <- colMeans((chain.1[[17]][XG.burnin:XG, ] + chain.2[[17]][XG.burnin:XG, ]) / 2)
  sd.nu <- apply((chain.1[[17]][XG.burnin:XG, ] + chain.2[[17]][XG.burnin:XG, ]) / 2, FUN = sd, MARGIN = 2)
  mce.nu <- sd.nu / sqrt(2*(XG - XG.burnin))

  data.chain1 <- data.frame("tau" = chain.1[[9]][,,], "delta" = chain.1[[10]], "nu" = chain.1[[17]], "sig2" = chain.1[[8]], "sig2k" = chain.1[[7]])
  data.chain2 <- data.frame("tau" = chain.2[[9]][,,], "delta" = chain.2[[10]], "nu" = chain.2[[17]], "sig2" = chain.2[[8]], "sig2k" = chain.2[[7]])


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
          theta_i.s0 <- chain[[5]][ii-1,,gg]
          ones <- rep(1, K)
          varinv <- diag(1/(rep(1, K)))
          var.g <- (t(ones) %*% varinv) %*% ones

          # Hyper parameters
          SS <- b.theta_i + sum((theta_i.s0 - theta.g)^2) + (N.g*n0.theta_i*theta.g)/(2*(N.g + n0.theta_i))
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
          zeta_i.s0 <- chain[[6]][ii-1,,gg]
          ones <- rep(1, K)
          varinv <- diag(1/(sig2k.g[1:K]))
          var.g <- (t(ones) %*% varinv) %*% ones

          # Hyper parameters
          SS <- b.zeta_i + sum((zeta_i.s0 - zeta.g)^2) + (N.g*n0.zeta_i*zeta.g)/(2*(N.g + n0.zeta_i))
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

  return(list(beta = post.beta, lambda = post.lambda, theta = post.theta, zeta = post.zeta, sig2k = post.sig2k, sig2 = post.sig2,
              tau = post.tau, delta = post.delta, theta_i = post.theta_i, zeta_i = post.zeta_i, nu = post.nu, sd.nu = sd.nu, mce.nu = mce.nu,
              ZT2 = ZT2, data.chain1 = data.chain1, data.chain2 = data.chain2, burnin = XG.burnin))
}
