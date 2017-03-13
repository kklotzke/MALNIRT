#' @importFrom truncnorm rtruncnorm
#' @importFrom invgamma dinvgamma
#' @export
MALNIRT <- function(Y, RT, Group = NULL, data, XG = 1000, XG.init = 100, burnin = 0.10, est.person = FALSE, silent = FALSE) {

  ###### Initialize ######
  if (!missing(data)) {
    # Try to find Y, RT and Group in the data set first
    tryCatch(Y <- eval(substitute(Y), data), error=function(e) NULL)
    tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
    if (!is.null(Group)) {
      tryCatch(Group <- eval(substitute(Group), data), error=function(e) NULL)
    }
    else {
      Group <- rep(1, nrow(Y))
    }
  } else {
    data = NULL
  }

  K <- ncol(Y) # Number of items
  N <- nrow(Y) # Number of persons

  # Split data for groups
  groups <- unique(Group)
  G <- length(groups)
  Yg <- vector("list", G)
  RTg <- vector("list", G)
  param <- vector("list", G)

  for (g in 1:G)
  {
    Yg[[g]] <- subset(Y, Group == groups[g])
    RTg[[g]] <- subset(RT, Group == groups[g])
    param[[g]]$Z.mar <- matrix(0, ncol = K, nrow = N)
    param[[g]]$beta.mar <- numeric(K)
    param[[g]]$theta.mar <- 0
    param[[g]]$Z <- matrix(0, ncol = K, nrow = N)
    param[[g]]$mu.Z <- matrix(0, ncol = K, nrow = N)
    param[[g]]$var.Z <- numeric(K) + 0.01
    param[[g]]$r <- matrix(0, ncol = K, nrow = N)
    param[[g]]$x <- matrix(0, ncol = K, nrow = N)
    param[[g]]$tau.cand <- 0
    param[[g]]$delta.cand <- 0
    param[[g]]$sig2k.cand <- numeric(K)
    param[[g]]$sig2.cand <- 0
    param[[g]]$nu.cand <- numeric(K)
    param[[g]]$q.tau <- vector("list", 0)
    param[[g]]$q.delta <- vector("list", 0)
    param[[g]]$q.sig2k <- vector("list", 0)
  }


  # Create chains
  chains <- vector("list", 2) # Two chains
  chains[[1]] <- vector("list", G)
  for (g in 1:G) {
    chains[[1]][[g]] <- vector("list", 9)

    # Item parameters
    chains[[1]][[g]][[1]] <- matrix(NA, nrow = XG, ncol = K) # Item difficulty
    chains[[1]][[g]][[2]] <- matrix(NA, nrow = XG, ncol = K) # Time intensity

    # Group parameters
    chains[[1]][[g]][[3]] <- matrix(NA, nrow = XG, ncol = 1) # Ability group mean
    chains[[1]][[g]][[4]] <- matrix(NA, nrow = XG, ncol = 1) # Speed group mean

    # Covariance matrix Z
    chains[[1]][[g]][[5]] <- matrix(NA, nrow = XG, ncol = 1) # tau

    # Covariance matrix RT
    chains[[1]][[g]][[6]] <- matrix(NA, nrow = XG, ncol = 1) # delta
    chains[[1]][[g]][[7]] <- matrix(NA, nrow = XG, ncol = K) # sig2 per item
    chains[[1]][[g]][[8]] <- matrix(NA, nrow = XG, ncol = 1) # mean sig2

    # Cross-covariance matrix
    chains[[1]][[g]][[9]] <- matrix(NA, nrow = XG, ncol = K) # nu per item

    chains[[2]][[g]] <- chains[[1]][[g]]
  }


  # Run the marginalized ability model to obtain sane starting values for the joint-model
  initMar <- function(e) {
    firstGroup <- TRUE
    for (g in 1:e$G) {
      if(g > 1)
        firstGroup <- FALSE

      for (xg in 1:XG.init) {
        out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[g]]$beta.mar,
                                     theta.mar = e$param[[g]]$theta.mar, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup)
        e$param[[g]]$Z.mar <- out.ab$Z.mar
        e$param[[g]]$beta.mar <- out.ab$beta.mar
        e$param[[g]]$theta.mar <- out.ab$theta.mar
        e$param[[g]]$tau.cand <- out.ab$tau.mar
        e$param[[g]]$q.tau <- out.ab$q.tau
      }

      # Use latent responses from marginalized model as sane starting values
      e$param[[g]]$Z <- e$param[[g]]$mu.Z <- e$param[[g]]$r <- e$param[[g]]$Z.mar
      e$param[[g]]$x <- e$param[[g]]$r + rnorm(N*K, 0, 1)
    }
  }


  # Initialize chains
  initChains <- function(reinit = FALSE, e) {
    for (c in 1:2) {
      for (g in 1:e$G) {
        # Item parameters
        e$chains[[c]][[g]][[1]][1, ] <- e$param[[g]]$beta.mar  # Item difficulty
        e$chains[[c]][[g]][[2]][1, ] <- rnorm(e$K, 5, 1) # Time intensity

        # Group parameters
        e$chains[[c]][[g]][[3]][1, ] <- rnorm(1, 0, 1) # Ability group mean
        e$chains[[c]][[g]][[4]][1, ] <-  rnorm(1, 0, 1) # Speed group mean

        # Covariance matrix Z
        e$chains[[c]][[g]][[5]][1, ] <- e$param[[g]]$tau.cand # tau

        # Covariance matrix RT
        e$chains[[c]][[g]][[6]][1, ] <- runif(1, 0, 1) # delta
        e$chains[[c]][[g]][[7]][1, ] <- runif(e$K, 0.5, 1.5) # sig2 per item
        e$chains[[c]][[g]][[8]][1, ] <- mean(e$chains[[c]][[g]][[7]][1, ]) # mean sig2

        # Cross-covariance matrix
        e$chains[[c]][[g]][[9]][1, ] <- diag(cov(e$param[[g]]$Z.mar, e$RTg[[g]])) # nu per item

        if(reinit) {
          e$chains[[c]][[g]][[6]][1, ] <- e$param[[g]]$delta.cand # delta
          e$chains[[c]][[g]][[7]][1, ] <- e$param[[g]]$sig2k.cand  # sig2 per item
          e$chains[[c]][[g]][[8]][1, ] <- e$param[[g]]$sig2.cand # mean sig2
        }
      }
    }
  }

  ################################################# Proposals ################################################
  sampleProposals <- function(lambda, zeta, init = FALSE, g, e) {

    ### Sample proposal for tau from marginalized ability model ###
    out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[g]]$beta.mar,
                                    theta.mar = e$param[[g]]$theta.mar, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup)
    e$param[[g]]$Z.mar <- out.ab$Z.mar
    e$param[[g]]$beta.mar <- out.ab$beta.mar
    e$param[[g]]$theta.mar <- out.ab$theta.mar
    e$param[[g]]$tau.cand <- out.ab$tau.mar
    e$param[[g]]$q.tau <- out.ab$q.tau

    ### Sample proposal for delta, sig2k and sig2 from marginalized speed model ###
    out.sp <- sampleMarSpeedModel(RT = e$RTg[[g]], lambda = lambda, zeta = zeta, delta.mar = e$param[[g]]$delta.cand, sig2k.mar = e$param[[g]]$sig2k.cand)

    e$param[[g]]$delta.cand <- out.sp$delta.mar
    e$param[[g]]$sig2k.cand <- out.sp$sig2k.mar
    e$param[[g]]$sig2.cand <- out.sp$sig2.mar
    e$param[[g]]$q.delta <- out.sp$q.delta
    e$param[[g]]$q.sig2k <- out.sp$q.sig2k
    e$param[[g]]$q.sig2 <- out.sp$q.sig2

    ### Sample proposal for nu as regression parameter ###
    out.cor <- sampleCorrelationNu (r = e$param[[g]]$r, x = e$param[[g]]$x, var.Z = e$param[[g]]$var.Z)
    e$param[[g]]$nu.cand <- out.cor$nu
    e$param[[g]]$q.nu <- out.cor$q.nu

    if(init)
      e$param[[g]]$nu.cand <- diag(cov(e$param[[g]]$Z.mar, e$RTg[[g]]))
  }

  ############################################################################################################


  ########################################### Metropolis-Hastings ############################################

  doMH <- function(e)
  {
    I <- diag(e$K)
    J <- matrix(1, nrow = e$K, ncol = e$K)
    a <- 1/e$delta.s0[1] + sum(1/e$sig2k.s0)
    a.cand <- 1/e$param[[g]]$delta.cand + sum(1/e$param[[g]]$sig2k.cand)

    muZ <- e$theta - e$beta
    muT <- e$lambda - e$zeta

    # Conditional mean of Z|T for each item
    x <- x.cand <- mu.ZT <- mu.ZT.cand <- matrix(NA, ncol = K, nrow = Ng)
    for (k in 1:K)
    {
      x[, k] <- (1/e$sig2k.s0[k]) * ((e$RTg[[g]][, k] - muT[k]) - (e$RTg[[g]][, k] - muT[k]) / (a * e$sig2k.s0[k]))
      x.cand[, k] <- (1/e$param[[g]]$sig2k.cand[k]) * ((e$RTg[[g]][, k] - muT[k]) - (e$RTg[[g]][, k] - muT[k]) / (a.cand * e$param[[g]]$sig2k.cand[k]))

      mu.ZT[, k] <- muZ[k] + e$nu.s0[k] * x[, k]
      mu.ZT.cand[, k] <- muZ[k] + e$param[[g]]$nu.cand[k] * x.cand[, k]
    }

    I.min1 <- diag(e$K-1)
    J.min1 <- matrix(1, nrow = e$K-1, ncol = e$K-1)
    ones.min1 <- rep(1, e$K-1)

    # Partioned matrix components for each of the K items
    partMatrix <- partMatrix.cand <- vector("list", e$K)
    for (k in 1:e$K)
    {
      w.min1 <- e$nu.s0[-k]/e$sig2k.s0[-k]
      b.min1 <- e$sig2k.s0[-k] / (e$sig2k.s0[-k] - e$nu.s0[-k]^2)
      c.min1 <- 1/e$tau.s0[1] + sum(b.min1)
      A.min1.inv <- b.min1*I.min1 - (b.min1 %*% t(b.min1)) / c.min1
      d.min1 <- a + t(w.min1) %*% A.min1.inv %*% w.min1
      g.min1 <- sum(e$nu.s0[-k] / (e$sig2k.s0[-k] - e$nu.s0[-k]^2)) / c.min1
      A.min1.inv_w <- (e$nu.s0[-k] / (e$sig2k.s0[-k] - e$nu.s0[-k]^2) - g.min1*b.min1)

      partMatrix[[k]]$B11 <- 1 + e$tau.s0[1] - e$nu.s0[k]^2 * (1/e$sig2k.s0[k] - ((1/e$sig2k.s0[k])^2) / a)
      partMatrix[[k]]$B22 <- (1 - (e$nu.s0[-k]^2) / e$sig2k.s0[-k])*I.min1 + e$tau.s0[1]*J.min1 + w.min1 %*% t(w.min1) / a
      partMatrix[[k]]$B12 <- e$tau.s0[1]*ones.min1 + ((e$nu.s0[k] / e$sig2k.s0[k]) %*% t(w.min1)) / a
      partMatrix[[k]]$B21 <- t(partMatrix[[k]]$B12)
      partMatrix[[k]]$B22.inv <- A.min1.inv - (A.min1.inv_w %*% t(A.min1.inv_w)) / d.min1[1,1]

      ### MH candidates ###
      w.min1.cand <- e$param[[g]]$nu.cand[-k]/e$param[[g]]$sig2k.cand[-k]
      b.min1.cand <- e$param[[g]]$sig2k.cand[-k] / (e$param[[g]]$sig2k.cand[-k] - e$param[[g]]$nu.cand[-k]^2)
      c.min1.cand <- 1/e$param[[g]]$tau.cand + sum(b.min1.cand)
      A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
      d.min1.cand <- a.cand + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
      g.min1.cand <- sum(e$param[[g]]$nu.cand[-k] / (e$param[[g]]$sig2k.cand[-k] - e$param[[g]]$nu.cand[-k]^2)) / c.min1.cand
      A.min1.inv_w.cand <- (e$param[[g]]$nu.cand[-k] / (e$param[[g]]$sig2k.cand[-k] - e$param[[g]]$nu.cand[-k]^2) - g.min1.cand*b.min1.cand)

      partMatrix.cand[[k]]$B11 <- 1 + e$param[[g]]$tau.cand - e$param[[g]]$nu.cand[k]^2 * (1/e$param[[g]]$sig2k.cand[k] - ((1/e$param[[g]]$sig2k.cand[k])^2) / a.cand)
      partMatrix.cand[[k]]$B22 <- (1 - (e$param[[g]]$nu.cand[-k]^2) / e$param[[g]]$sig2k.cand[-k])*I.min1 + e$param[[g]]$tau.cand*J.min1 + w.min1.cand %*% t(w.min1.cand) / a.cand
      partMatrix.cand[[k]]$B12 <- e$param[[g]]$tau.cand*ones.min1 + ((e$param[[g]]$nu.cand[k] / e$param[[g]]$sig2k.cand[k]) %*% t(w.min1.cand)) / a.cand
      partMatrix.cand[[k]]$B21 <- t(partMatrix.cand[[k]]$B12)
      partMatrix.cand[[k]]$B22.inv <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]
    }

    # Update mu.Z
    out.Z <- tryCatch({sampleZ(Y = e$Yg[[g]], e$param[[g]]$Z, e$param[[g]]$mu.Z, mu.ZT, partMatrix, likelihood = FALSE)},
      error = function(er) { return(NULL) })

    # Compute likelihood
    out.Z <- tryCatch({sampleZ(Y = e$Yg[[g]], out.Z$Z, out.Z$mu.Z, mu.ZT, partMatrix, likelihood = TRUE)},
                      error = function(er) { return(NULL) })

    # If current state is faulty, re-initialized chain
    reset <- FALSE
    if(is.null(out.Z))
      reset <- TRUE

    validProposals <- FALSE
    if(!reset) {
      # Update mu.Z for proposals
      out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], out.Z$Z, out.Z$mu.Z, mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
                        error = function(er) { return(NULL) })

      # Compute likelihood for proposals
      out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], out.Z.cand$Z, out.Z.cand$mu.Z, mu.ZT.cand, partMatrix.cand, likelihood = TRUE)},
                        error = function(er) { return(NULL) })

      if(!is.null(out.Z.cand)) {
        validProposals <- TRUE

        lik <- sum(out.Z$lik_k) + invgamma::dinvgamma(e$param[[g]]$tau.cand + 1/K, e$param[[g]]$q.tau$shape, e$param[[g]]$q.tau$rate, log = TRUE) + invgamma::dinvgamma(e$param[[g]]$delta.cand + 1/K, e$param[[g]]$q.delta$shape, e$param[[g]]$q.delta$rate, log = TRUE)
        lik.cand <- sum(out.Z.cand$lik_k) + invgamma::dinvgamma(e$tau.s0 + 1/K, e$param[[g]]$q.tau$shape, e$param[[g]]$q.tau$rate, log = TRUE) + invgamma::dinvgamma(e$delta.s0 + 1/K, e$param[[g]]$q.delta$shape, e$param[[g]]$q.delta$rate, log = TRUE)
        for (k in 1:K) {
          lik <- lik + dnorm(e$param[[g]]$nu.cand[k], e$param[[g]]$q.nu$mu[k], e$param[[g]]$q.nu$sd[k], TRUE) + invgamma::dinvgamma(e$param[[g]]$sig2k.cand[k], e$param[[g]]$q.sig2k$shape, e$param[[g]]$q.sig2k$rate[k], log = TRUE)
          lik.cand <- lik.cand + dnorm(e$nu.s0[k], e$param[[g]]$q.nu$mu[k], e$param[[g]]$q.nu$sd[k], TRUE) + invgamma::dinvgamma(e$sig2k.s0[k], e$param[[g]]$q.sig2k$shape, e$param[[g]]$q.sig2k$rate[k], log = TRUE)
        }

        ar <- exp(lik.cand - lik)
        if(is.nan(ar))
          ar <- 0
      }
      else
        ar <- 0
      u <- runif (1, 0, 1)

      if (u <= ar) {
        #cat("Accept tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[5]][xg, ] <- e$param[[g]]$tau.cand
        e$chains[[c]][[g]][[6]][xg, ] <- e$param[[g]]$delta.cand
        e$chains[[c]][[g]][[7]][xg, ] <- e$param[[g]]$sig2k.cand
        e$chains[[c]][[g]][[8]][xg, ] <- e$param[[g]]$sig2.cand
        e$chains[[c]][[g]][[9]][xg, ] <- e$param[[g]]$nu.cand
        e$param[[g]]$mu.Z <- out.Z.cand$mu.Z
        e$param[[g]]$var.Z <- out.Z.cand$var.Z
        e$param[[g]]$Z <- out.Z.cand$Z
        for (k in 1:K)
          e$param[[g]]$r[, k] <- (out.Z.cand$Z[, k]- muZ[k]) - out.Z.cand$tmp[[k]][1, ]
        e$param[[g]]$x <- x.cand

        e$mh.accept <- e$mh.accept + 1
      }
      else {
        #print("Reject tau")
        #cat("Reject tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[5]][xg, ] <- e$tau.s0
        e$chains[[c]][[g]][[6]][xg, ] <- e$delta.s0
        e$chains[[c]][[g]][[7]][xg, ] <- e$sig2k.s0
        e$chains[[c]][[g]][[8]][xg, ] <- e$sig2.s0
        e$chains[[c]][[g]][[9]][xg, ] <- e$nu.s0
        e$param[[g]]$mu.Z <- out.Z$mu.Z
        e$param[[g]]$var.Z <- out.Z$var.Z
        e$param[[g]]$Z <- out.Z$Z
        for (k in 1:K)
          e$param[[g]]$r[, k] <- (out.Z$Z[, k]- muZ[k]) - out.Z$tmp[[k]][1, ]
        e$param[[g]]$x <- x
      }
    }

    #cat(reset, " | " , validProposals, "\n")
    #print(e$param[[g]]$nu.cand)
    #print(e$param[[g]]$delta.cand)
    #print(e$param[[g]]$sig2k.cand)
    return(list(reset = reset, validProposals = validProposals))
  }

  ############################################################################################################

  e <- environment()
  e$initMar(e)
  e$initChains(reinit = FALSE, e)

  ###### Run MCMC ######
  # For each chain
  for (c in 1:2) {

    chain <- chains[[c]] # Read the cc-th chain
    mh.accept <- 0 # MH acceptance rate

    xg <- 2
    while (xg <= XG) {


      ############################################## Gibbs sampling ##############################################

      ### beta and lambda do not differ between groups ###
      ### Assume equally sized groups
      beta.s0 <- chains[[c]][[g]][[1]][xg-1, ]
      lambda.s0 <- chains[[c]][[g]][[2]][xg-1, ]
      Z.all <- param[[1]]$Z
      RT.all <- RTg[[1]]
      theta.all <- chains[[c]][[1]][[3]][xg-1, ]
      zeta.all <- chains[[c]][[1]][[4]][xg-1, ]
      tau.all <- chains[[c]][[1]][[5]][xg-1, ]
      delta.all <- chains[[c]][[1]][[6]][xg-1, ]
      sig2k.all <- chains[[c]][[1]][[7]][xg-1, ]
      g <- 2
      while(g <= G) {
        Z.all <- rbind(Z.all, param[[g]]$Z)
        RT.all <- rbind(RT.all, RTg[[g]])
        theta.all <- theta.all + chains[[c]][[g]][[3]][xg-1, ]
        zeta.all <- zeta.all + chains[[c]][[g]][[4]][xg-1, ]
        tau.all <- tau.all + chains[[c]][[g]][[5]][xg-1, ]
        delta.all <- delta.all + chains[[c]][[g]][[6]][xg-1, ]
        sig2k.all <- sig2k.all + chains[[c]][[g]][[7]][xg-1, ]
        g <- g + 1
      }
      theta.all <- theta.all / G
      zeta.all <- zeta.all / G
      tau.all <- tau.all / G
      delta.all <- delta.all / G
      sig2k.all <- sig2k.all / G

      ### Sample item difficulty paramaters ###
      chains[[c]][[1]][[1]][xg, ] <- beta <- sampleBeta(Z = Z.all, beta = beta.s0, theta = theta.all, tau = tau.all)

      ### Sample item time intensity paramaters ###
      chains[[c]][[1]][[2]][xg, ] <- lambda <- sampleLambda(RT = RT.all, lambda = lambda.s0, zeta = zeta.all, sig2k = sig2k.all, delta = delta.all)

      g <- 2
      while(g <= G) {
        chains[[c]][[g]][[1]][xg, ] <- beta
        chains[[c]][[g]][[2]][xg, ] <- lambda
        g <- g + 1
      }

      firstGroup <- TRUE
      for(g in 1:G) {

        Ng <- nrow(RTg[[g]])

        if(g > 1)
          firstGroup <- FALSE

        #### Read values from former chain (t-1) ###
        theta.s0 <- chains[[c]][[g]][[3]][xg-1, ]
        zeta.s0 <- chains[[c]][[g]][[4]][xg-1, ]
        tau.s0 <- chains[[c]][[g]][[5]][xg-1, ]
        delta.s0 <- chains[[c]][[g]][[6]][xg-1, ]
        sig2k.s0 <- chains[[c]][[g]][[7]][xg-1, ]
        sig2.s0 <- chains[[c]][[g]][[8]][xg-1, ]
        nu.s0 <- chains[[c]][[g]][[9]][xg-1, ]

        ### Sample group ability parameter ###
        if (g == 1)
          chains[[c]][[g]][[3]][xg, ] <- theta <- 0
        else
          chains[[c]][[g]][[3]][xg, ] <- theta <- sampleTheta(Z = param[[g]]$Z, beta = beta, tau = tau.s0)

        ### Sample group speed parameter ###
        if (g == 1)
          chains[[c]][[g]][[4]][xg, ] <- zeta <- 0
        else
          chains[[c]][[g]][[4]][xg, ] <- zeta <- sampleZeta(RT = RTg[[g]], lambda = lambda, sig2k = sig2k.s0, delta = delta.s0)

      ############################################################################################################


      ################################################# Proposals ################################################

        init <- FALSE
        if (xg < 3)
          init <- TRUE
        sampleProposals(lambda = lambda, zeta = zeta, init = init, g = g, e = e)

      ############################################################################################################


      ########################################### Metropolis-Hastings ############################################

        out.mh <- doMH(e)
        while (!out.mh$reset && !out.mh$validProposals) {
          sampleProposals(lambda = lambda, zeta = zeta, init = init, g = g, e = e)
          out.mh <- doMH(e)
        }

      ############################################################################################################

        if(out.mh$reset) {
          e$initMar(e)
          e$initChains(reinit = TRUE, e)
          mh.accept <- 0
          xg <- 2
        }
        else {
          if ((!silent) && (xg%%100 == 0))
            cat("Iteration ", xg, " | MH acceptance rate ", mh.accept/xg, "\n")

          if((mh.accept/xg < 0.05) && (xg%%100 == 0)) {
            e$initMar(e)
            e$initChains(reinit = TRUE, e)
            mh.accept <- 0
            xg <- 2
            print("Reset")
          }
          else {
            xg <- xg + 1
          }


        }

      }

    }
  }
}
