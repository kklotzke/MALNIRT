# @importFrom truncnorm rtruncnorm
#' @importFrom invgamma dinvgamma
#' @importFrom extraDistr rinvgamma dtnorm rtnorm
#' @export
MALNIRT3Steps <- function(Y, RT, group = NULL, data, XG = 1000, XG.init = 100, burnin = 0.10, est.person = FALSE,
                    doBIC.zeta = FALSE, doBIC.theta = FALSE, doBIC.nu = FALSE, silent = FALSE) {

  ###### Initialize ######
  if (!missing(data)) {
    # Try to find Y, RT and Group in the data set first
    tryCatch(Y <- eval(substitute(Y), data), error=function(e) NULL)
    tryCatch(RT <- eval(substitute(RT), data), error=function(e) NULL)
    if (!is.null(group)) {
      tryCatch(group <- eval(substitute(group), data), error=function(e) NULL)
    }
    else {
      group <- rep(1, nrow(Y))
    }
  } else {
    data = NULL
  }

  K <- ncol(Y) # Number of items
  N <- nrow(Y) # Number of persons

  # Split data for groups
  groups <- unique(group)
  G <- length(groups)
  Yg <- vector("list", G)
  RTg <- vector("list", G)
  param <- vector("list", G)

  for (g in 1:G)
  {
    Yg[[g]] <- subset(Y, group == groups[g])
    RTg[[g]] <- subset(RT, group == groups[g])
    Ng <- nrow(Yg[[g]])

    param[[g]]$Z.mar <- matrix(0, ncol = K, nrow = Ng)
    param[[g]]$beta.mar <- numeric(K)
    param[[g]]$theta.mar <- 0
    param[[g]]$lambda.mar <- numeric(K)
    param[[g]]$zeta.mar <- 0
    param[[g]]$sig2k.mar <- runif(K, 0.5, 1.5)
    param[[g]]$Z <- matrix(0, ncol = K, nrow = Ng)
    param[[g]]$mu.Z <- matrix(0, ncol = K, nrow = Ng)
    param[[g]]$var.Z <- numeric(K) + 0.01
    param[[g]]$r <- matrix(0, ncol = K, nrow = Ng)
    param[[g]]$x <- matrix(0, ncol = K, nrow = Ng)
    param[[g]]$tau.cand <- 0
    param[[g]]$delta.cand <- 0
    param[[g]]$nu.cand <- numeric(K)
    param[[g]]$uz <- matrix(0, nrow = Ng, ncol = K)
    param[[g]]$ut <- matrix(0, nrow = Ng, ncol = K)
    param[[g]]$mu.ZT <- matrix(0, ncol = K, nrow = Ng)
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
    chains[[1]][[g]][[7]] <- matrix(runif(K, 0.5, 1.5), nrow = XG, ncol = K) # sig2 per item
    chains[[1]][[g]][[8]] <- matrix(1, nrow = XG, ncol = 1) # mean sig2

    # Cross-covariance matrix
    chains[[1]][[g]][[9]] <- matrix(NA, nrow = XG, ncol = K) # nu per item

    chains[[2]][[g]] <- chains[[1]][[g]]
  }


  # Run the marginalized ability model to obtain sane starting values for the joint-model
  initMar <- function(e) {

    #cat("#### Initialize marginalized models #### \n")
    firstGroup <- TRUE
    for (g in 1:e$G) {
      if(g > 1)
        firstGroup <- FALSE

      for (xg in 1:XG.init) {

        if(firstGroup) {
          out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[g]]$beta.mar,
                                          theta.mar = e$param[[g]]$theta.mar, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup)
          out.sp <- sampleMarSpeedModel(RT = e$RTg[[g]], lambda = e$param[[g]]$lambda.mar, zeta = e$param[[g]]$zeta.mar,
                                        delta.mar = e$param[[g]]$delta.cand, sig2k.mar = e$param[[g]]$sig2k.mar)
          lambda.mar <- sampleLambda(RT = e$RTg[[g]], lambda = e$param[[g]]$lambda.mar, zeta = e$param[[g]]$zeta.mar,
                                     sig2k = out.sp$sig2k.mar, delta = out.sp$delta.mar)
          zeta.mar <- 0
        }
        else {
          out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[1]]$beta.mar,
                                          theta.mar = e$param[[g]]$theta.mar, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup) # Use beta's of first group, where theta = 0
          out.sp <- sampleMarSpeedModel(RT = e$RTg[[g]], lambda = e$param[[1]]$lambda.mar, zeta = e$param[[g]]$zeta.mar,
                                        delta.mar = e$param[[g]]$delta.cand, sig2k.mar = e$param[[g]]$sig2k.mar)

          zeta.mar <- sampleZeta(RT = e$RTg[[g]], lambda = e$param[[1]]$lambda.mar, sig2k = e$param[[g]]$sig2k.mar, delta = e$param[[g]]$delta.cand)
        }

        e$param[[g]]$Z.mar <- out.ab$Z.mar
        e$param[[g]]$beta.mar <- out.ab$beta.mar
        e$param[[g]]$theta.mar <- out.ab$theta.mar
        e$param[[g]]$tau.cand <- out.ab$tau.mar

        if(firstGroup) {
          e$param[[g]]$lambda.mar <- lambda.mar
        }
        else {
          e$param[[g]]$lambda.mar <- e$param[[1]]$lambda.mar
        }
        e$param[[g]]$zeta.mar <- zeta.mar
        e$param[[g]]$delta.cand <- out.sp$delta.mar
        e$param[[g]]$sig2k.mar <- out.sp$sig2k.mar
      }

      # Use latent responses from marginalized model as sane starting values
      e$param[[g]]$Z <- e$param[[g]]$mu.Z <- e$param[[g]]$r <- e$param[[g]]$Z.mar
      e$param[[g]]$x <- e$param[[g]]$r + rnorm(nrow(Yg[[g]])*K, 0, 1)
    }

    #print("Done initializing")
    #browser()
  }


  # Initialize chains
  initChains <- function(reinit = FALSE, which, e) {
    #cat("#### Initialize chains #### \n")

    for (c in which) {
      for (g in 1:e$G) {
        # Item parameters
        e$chains[[c]][[g]][[1]][1, ] <- e$param[[1]]$beta.mar  # Item difficulty
        e$chains[[c]][[g]][[2]][1, ] <- e$param[[1]]$lambda.mar #norm(e$K, 5, 1) # Time intensity

        # Group parameters
        e$chains[[c]][[g]][[3]][1, ] <- e$param[[g]]$theta.mar #rnorm(1, 0, 1) # Ability group mean
        e$chains[[c]][[g]][[4]][1, ] <- e$param[[g]]$zeta.mar #rnorm(1, 0, 1) # Speed group mean

        # Covariance matrix Z
        e$chains[[c]][[g]][[5]][1, ] <- e$param[[g]]$tau.cand # tau

        # Covariance matrix RT
        e$chains[[c]][[g]][[6]][1, ] <- e$param[[g]]$delta.cand # runif(1, 0, 1) # delta
        e$chains[[c]][[g]][[7]][1, ] <- e$param[[g]]$sig2k.mar # runif(e$K, 0.5, 1.5) # sig2 per item
        e$chains[[c]][[g]][[8]][1, ] <- mean(e$param[[g]]$sig2k.mar) #mean(e$chains[[c]][[g]][[7]][1, ]) # mean sig2

        # Cross-covariance matrix
        e$chains[[c]][[g]][[9]][1, ] <- diag(cov(e$param[[g]]$Z.mar, e$RTg[[g]])) # nu per item

        # if(reinit) {
        #   e$chains[[c]][[g]][[6]][1, ] <- e$param[[g]]$delta.cand # delta
        #   e$chains[[c]][[g]][[7]][1, ] <- e$param[[g]]$sig2k.cand  # sig2 per item
        #   e$chains[[c]][[g]][[8]][1, ] <- e$param[[g]]$sig2.cand # mean sig2
        # }
      }
    }
  }

  ################################################# Proposals ################################################
  sampleProposals <- function(reinit = FALSE, g, e) {

    ### Sample proposal for tau from marginalized ability model ###
    #out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[g]]$beta.mar,
    #                                theta.mar = e$param[[g]]$theta.mar, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup)
    #out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$param[[g]]$beta.mar,
    #                                theta.mar = 0, tau.mar <- e$param[[g]]$tau.cand, firstGroup = firstGroup, init = FALSE)
    out.ab <- sampleMarAbilityModel(Y = e$Yg[[g]], Z.mar = e$param[[g]]$Z.mar, beta.mar = e$beta,
                                     theta.mar = 0, tau.mar = e$param[[g]]$tau.cand, firstGroup = e$firstGroup, init = FALSE)
    e$param[[g]]$Z.mar <- out.ab$Z.mar
    e$param[[g]]$beta.mar <- out.ab$beta.mar
    e$param[[g]]$theta.mar <- out.ab$theta.mar
    e$param[[g]]$tau.cand <- out.ab$tau.mar

    ### Sample proposal for delta, sig2k and sig2 from marginalized speed model ###
    #out.sp <- sampleMarSpeedModel(RT = e$RTg[[g]], lambda = e$lambda, zeta = e$zeta, delta.mar = e$param[[g]]$delta.cand, sig2k.mar = e$param[[g]]$sig2k.cand)
    out.sp <- sampleMarSpeedModel(RT = e$RTg[[g]], lambda = e$lambda, zeta = 0, delta.mar = e$param[[g]]$delta.cand, sig2k.mar = e$param[[g]]$sig2k.mar)

    e$param[[g]]$delta.cand <- out.sp$delta.mar
    e$param[[g]]$sig2k.mar <- out.sp$sig2k.mar

    ### Sample proposal for nu as regression parameter ###

    # if(g==1) {
    out.cor <- sampleCorrelationNu (r = e$param[[g]]$r, x = e$param[[g]]$x, var.Z = e$param[[g]]$var.Z, Z = e$param[[g]]$Z)
    # if(g==1) {
    #   print(out.cor$nu)
    # }

    e$param[[g]]$nu.cand <- out.cor$nu #+ rnorm(K, 0, 0.2) # diag(cov(e$param[[g]]$Z, e$RTg[[g]])) #rep(0, K) out.cor$nu #
    #e$param[[g]]$q.nu <- out.cor$q.nu

    # if(g == 1)
    #   e$param[[g]]$nu.cand <- dat1$nu
    # else
    #   e$param[[g]]$nu.cand <- dat2$nu


    if(reinit)
      e$param[[g]]$nu.cand <- diag(cov(e$param[[g]]$Z.mar, e$RTg[[g]]))

    #if(g == 2)
    #  print(e$param[[g]]$tau.cand)
  }

  ############################################################################################################


  ########################################### Metropolis-Hastings ############################################

  doMH <- function(e)
  {
    ### Step 1 ####
    reset <- FALSE

    I <- diag(e$K)
    J <- matrix(1, nrow = e$K, ncol = e$K)
    a <- 1/e$delta.s0 + sum(1/e$sig2k)

    muZ <- e$theta - e$beta
    muT <- e$lambda - e$zeta

    Sigma.ZT <- ( diag(K) + J*e$tau.s0[1] ) - diag(e$nu.s0) %*% solve( diag(e$sig2k) + J*e$delta.s0[1] ) %*% t(diag(e$nu.s0))
    Sigma.ZT.cand <- ( diag(K) + J*e$param[[g]]$tau.cand ) - diag(e$nu.s0) %*% solve( diag(e$sig2k) + J*e$delta.s0[1] ) %*% t(diag(e$nu.s0))

    # Independent proposals may lead to a matrix that is not (semi-)positive-definite
    chol.ZT.cand <- tryCatch({chol(Sigma.ZT.cand)}, error = function(er) { return(NULL) })
    if(is.null(chol.ZT.cand))
      reset <- TRUE
    chol.ZT <- tryCatch({chol(Sigma.ZT)}, error = function(er) { return(NULL) })
    if(is.null(chol.ZT))
      reset <- TRUE

    if(!reset) {
      # Conditional mean of Z|T for each item
      x <- mu.ZT <- matrix(NA, ncol = K, nrow = Ng)
      for (k in 1:K)
      {
        x[, k] <- t( (1/e$sig2k[k]) * ( t(e$RTg[[g]][, k] - muT[k]) - (t(1/e$sig2k) %*% (t(e$RTg[[g]]) - matrix(muT, nrow = K, ncol = Ng))) / a ) )
        mu.ZT[, k] <- muZ[k] + e$nu.s0[k] * x[, k]
      }

      I.min1 <- diag(e$K-1)
      J.min1 <- matrix(1, nrow = e$K-1, ncol = e$K-1)
      ones.min1 <- rep(1, e$K-1)

      # Partioned matrix components for each of the K items
      partMatrix <- partMatrix.cand <- vector("list", e$K)
      for (k in 1:e$K)
      {
        w.min1 <- e$nu.s0[-k]/e$sig2k[-k]
        b.min1 <- e$sig2k[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2)
        c.min1 <- 1/e$tau.s0[1] + sum(b.min1)
        A.min1.inv <- b.min1*I.min1 - (b.min1 %*% t(b.min1)) / c.min1
        d.min1 <- a + t(w.min1) %*% A.min1.inv %*% w.min1
        g.min1 <- sum(e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2)) / c.min1
        A.min1.inv_w <- (e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2) - g.min1*b.min1)

        partMatrix[[k]]$B11 <- 1 + e$tau.s0[1] - e$nu.s0[k]^2 * (1/e$sig2k[k] - ((1/e$sig2k[k])^2) / a)
        partMatrix[[k]]$B22 <- (1 - (e$nu.s0[-k]^2) / e$sig2k[-k])*I.min1 + e$tau.s0[1]*J.min1 + w.min1 %*% t(w.min1) / a
        partMatrix[[k]]$B12 <- e$tau.s0[1]*ones.min1 + ((e$nu.s0[k] / e$sig2k[k]) %*% t(w.min1)) / a
        partMatrix[[k]]$B21 <- t(partMatrix[[k]]$B12)
        partMatrix[[k]]$B22.inv <- A.min1.inv - (A.min1.inv_w %*% t(A.min1.inv_w)) / d.min1[1,1]

        ### MH candidates ###
        w.min1.cand <- w.min1
        b.min1.cand <- b.min1
        c.min1.cand <- 1/e$param[[g]]$tau.cand + sum(b.min1.cand)
        A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
        d.min1.cand <- a + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
        g.min1.cand <- sum(e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2)) / c.min1.cand
        A.min1.inv_w.cand <- (e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2) - g.min1.cand*b.min1.cand)

        partMatrix.cand[[k]]$B11 <- 1 + e$param[[g]]$tau.cand - e$nu.s0[k]^2 * (1/e$sig2k[k] - ((1/e$sig2k[k])^2) / a)
        partMatrix.cand[[k]]$B22 <- (1 - (e$nu.s0[-k]^2) / e$sig2k[-k])*I.min1 + e$param[[g]]$tau.cand*J.min1 + w.min1.cand %*% t(w.min1.cand) / a
        partMatrix.cand[[k]]$B12 <- e$param[[g]]$tau.cand*ones.min1 + ((e$nu.s0[k] / e$sig2k[k]) %*% t(w.min1.cand)) / a
        partMatrix.cand[[k]]$B21 <- t(partMatrix.cand[[k]]$B12)
        partMatrix.cand[[k]]$B22.inv <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]
#
        # if(k == 10) {
        #   browser()
        #   upper <- cbind(partMatrix.cand[[k]]$B11[1], partMatrix.cand[[k]]$B12)
        #   lower <- cbind(partMatrix.cand[[k]]$B21, partMatrix.cand[[k]]$B22)
        #   Sigma_ZT2 <- rbind(upper, lower)
        #   print(all.equal(Sigma.ZT.cand, Sigma_ZT2))
        # }
      }

      # Sample Z
      #out.Z <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = e$param[[g]]$Z, mu.Z = e$param[[g]]$mu.Z, mu.ZT = mu.ZT, partMatrix, likelihood = FALSE)},
      #                    error = function(er) { return(NULL) })

      # If current state is faulty, re-initialized chain
      # if(is.null(out.Z)) {
      #   reset <- TRUE
      # }
    }

    validProposals <- FALSE
    if(!reset) {

      #out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = out.Z$Z, mu.Z = out.Z$mu.Z, mu.ZT = mu.ZT, partMatrix.cand, likelihood = FALSE)},
      #                         error = function(er) {  return(NULL) })
      #out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = e$param[[g]]$Z, mu.Z = e$param[[g]]$mu.Z, mu.ZT = mu.ZT, partMatrix.cand, likelihood = FALSE)},
     #                                              error = function(er) {  return(NULL) })

      # if(!is.null(out.Z.cand)) {
        validProposals <- TRUE

        # Evaluate likelihood 1 vs 2 in Z|T
        lik <- lik.cand <- 0
        for (i in 1:Ng) {
          #lik <- lik + sum(mvnfast::dmvn(X = out.Z$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT, log = TRUE, isChol = TRUE))
          #lik.cand <- lik.cand + sum(mvnfast::dmvn(X = out.Z.cand$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT.cand, log = TRUE, isChol = TRUE))
          lik <- lik + sum(mvnfast::dmvn(X = e$param[[g]]$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT, log = TRUE, isChol = TRUE))
          #lik <- lik + sum(mvnfast::dmvn(X = out.Z$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT, log = TRUE, isChol = TRUE))
          lik.cand <- lik.cand + sum(mvnfast::dmvn(X = e$param[[g]]$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT.cand, log = TRUE, isChol = TRUE))
          #lik.cand <- lik.cand + sum(mvnfast::dmvn(X = out.Z.cand$Z[i, ], mu = mu.ZT[i, ], sigma = chol.ZT.cand, log = TRUE, isChol = TRUE))

        }
        #browser()

        ar <- exp(lik.cand - lik)
        if(is.nan(ar) || is.na(ar))
          ar <- 0
        else
          ar <- min(1, ar)
      # }
      # else {
      #   ar <- 0
      #   #browser()
      # }
      u <- runif (1, 0, 1)

      if (u <= ar) {
        #cat("Accept tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[5]][xg, ] <- tau <- e$param[[g]]$tau.cand
        #e$param[[g]]$mu.Z <- out.Z.cand$mu.Z
        #e$param[[g]]$var.Z <- out.Z.cand$var.Z
        #e$param[[g]]$Z <- out.Z.cand$Z
        e$param[[g]]$mu.ZT <- mu.ZT

        #for (k in 1:K)
        #  e$param[[g]]$r[, k] <- (out.Z.cand$Z[, k]- muZ[k]) - out.Z.cand$tmp[[k]][1, ]
        #e$param[[g]]$x <- x

        e$mh.accept[g, 1] <- e$mh.accept[g, 1] + 1
        e$mh.accept.last[g, 1] <- e$mh.accept.last[g, 1] + 1

      }
      else {
        #print("Reject tau")
        #cat("Reject tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[5]][xg, ] <- tau <- e$tau.s0
        #e$param[[g]]$mu.Z <- out.Z$mu.Z
        #e$param[[g]]$var.Z <- out.Z$var.Z
        #e$param[[g]]$Z <- out.Z$Z
        e$param[[g]]$mu.ZT <- mu.ZT

        #for (k in 1:K)
        #  e$param[[g]]$r[, k] <- (out.Z$Z[, k]- muZ[k]) - out.Z$tmp[[k]][1, ]
        #e$param[[g]]$x <- x

      }
    }

    if (reset || !validProposals)
      return(list(reset = reset, validProposals = validProposals))


    ### Step 2 ###
    if (!reset) {
      a.cand <- 1/e$param[[g]]$delta.cand + sum(1/e$sig2k)
      Sigma.TZ <- diag(e$sig2k - e$nu.s0^2) + J*e$delta.s0 + (e$nu.s0 %*% t(e$nu.s0)) / (1 / tau + K)
      Sigma.TZ.cand <- diag(e$sig2k - e$nu.s0^2) + J*e$param[[g]]$delta.cand + (e$nu.s0 %*% t(e$nu.s0)) / (1 / tau + K)

      # Independent proposals may lead to a matrix that is not (semi-)positive-definite
      chol.TZ.cand <- tryCatch({chol(Sigma.TZ.cand)}, error = function(er) { return(NULL) })
      if(is.null(chol.TZ.cand))
        reset <- TRUE
      chol.TZ <- tryCatch({chol(Sigma.TZ)}, error = function(er) { return(NULL) })
      if(is.null(chol.TZ))
        reset <- TRUE


      if(!reset) {

        # Partioned matrix components for each of the K items
        partMatrix <- partMatrix.cand <- vector("list", e$K)
        for (k in 1:e$K)
        {
          ### MH candidates ###
          w.min1.cand <- e$nu.s0[-k]/e$sig2k[-k]
          b.min1.cand <- e$sig2k[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2)
          c.min1.cand <- 1/tau + sum(b.min1.cand)
          A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
          d.min1.cand <- a.cand + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
          g.min1.cand <- sum(e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2)) / c.min1.cand
          A.min1.inv_w.cand <- (e$nu.s0[-k] / (e$sig2k[-k] - e$nu.s0[-k]^2) - g.min1.cand*b.min1.cand)

          partMatrix.cand[[k]]$B11 <- 1 + tau - e$nu.s0[k]^2 * (1/e$sig2k[k] - ((1/e$sig2k[k])^2) / a.cand)
          partMatrix.cand[[k]]$B22 <- (1 - (e$nu.s0[-k]^2) / e$sig2k[-k])*I.min1 + tau*J.min1 + w.min1.cand %*% t(w.min1.cand) / a.cand
          partMatrix.cand[[k]]$B12 <- tau*ones.min1 + ((e$nu.s0[k] / e$sig2k[k]) %*% t(w.min1.cand)) / a.cand
          partMatrix.cand[[k]]$B21 <- t(partMatrix.cand[[k]]$B12)
          partMatrix.cand[[k]]$B22.inv <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]
        }

      # Conditional mean of Z|T for each item
      x.cand <- mu.ZT.cand <- matrix(NA, ncol = K, nrow = Ng)
      for (k in 1:K)
      {
        x.cand[, k] <- t( (1/e$sig2k[k]) * ( t(e$RTg[[g]][, k] - muT[k]) - (t(1/e$sig2k) %*% (t(e$RTg[[g]]) - matrix(muT, nrow = K, ncol = Ng))) / a.cand ) )
        mu.ZT.cand[, k] <- muZ[k] + e$nu.s0[k] * x.cand[, k]
        e$param[[g]]$uz[, k] <- e$nu.s0[k] * x.cand[, k]
      }

      #out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = e$param[[g]]$Z, mu.Z = e$param[[g]]$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
      #                       error = function(er) {  return(NULL) })

      validProposals <- FALSE
      # if(!is.null(out.Z.cand)) {
         validProposals <- TRUE

        mu.TZ <- t(muT + diag(e$nu.s0) %*% (diag(K) - J/(1/tau + K)) %*% t(e$param[[g]]$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))
        #mu.TZ.cand <- t(muT + diag(e$nu.s0) %*% (diag(K) - J/(1/tau + K)) %*% t(out.Z.cand$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))
        mu.TZ.cand <- t(muT + diag(e$nu.s0) %*% (diag(K) - J/(1/tau + K)) %*% t(e$param[[g]]$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))

        # Evaluate likelihood 1 vs 2 in Z|T
        lik <- lik.cand <- 0
        for (i in 1:Ng) {
          lik <- lik + sum(mvnfast::dmvn(X = e$RTg[[g]][i, ], mu = mu.TZ[i, ], sigma = chol.TZ, log = TRUE, isChol = TRUE))
          lik.cand <- lik.cand + sum(mvnfast::dmvn(X = e$RTg[[g]][i, ], mu = mu.TZ.cand[i, ], sigma = chol.TZ.cand, log = TRUE, isChol = TRUE))
        }

        ar <- exp(lik.cand - lik)
        if(is.nan(ar) || is.na(ar))
          ar <- 0
        else
          ar <- min(1, ar)
      #}
      # else {
      #   ar <- 0
      #   #browser()
      # }
      u <- runif (1, 0, 1)

      if (u <= ar) {
        #cat("Accept tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[6]][xg, ] <- delta <- e$param[[g]]$delta.cand
        #e$param[[g]]$mu.Z <- out.Z.cand$mu.Z
        #e$param[[g]]$var.Z <- out.Z.cand$var.Z
        #e$param[[g]]$Z <- out.Z.cand$Z
        e$param[[g]]$mu.ZT <- mu.ZT.cand

        #for(k in 1:K)
        #  e$param[[g]]$ut[, k] <- mu.TZ.cand[, k] - muT[k]

        #for (k in 1:K)
        #  e$param[[g]]$r[, k] <- (out.Z.cand$Z[, k]- muZ[k]) - out.Z.cand$tmp[[k]][1, ]
        #e$param[[g]]$x <- x.cand

        e$mh.accept[g, 2] <- e$mh.accept[g, 2] + 1
        e$mh.accept.last[g, 2] <- e$mh.accept.last[g, 2] + 1

      }
      else {
        #print("Reject tau")
        #cat("Reject tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
        e$chains[[c]][[g]][[6]][xg, ] <- delta <- e$delta.s0
       # for(k in 1:K)
         # e$param[[g]]$ut[, k] <- mu.TZ[, k] - muT[k]
      }
      }
    }

    if (reset || !validProposals)
      return(list(reset = reset, validProposals = validProposals))

    ### Step 3 ###
    if (!reset) {
      a.cand <- 1/delta + sum(1/e$sig2k)
      Sigma.TZ <- diag(e$sig2k - e$nu.s0^2) + J*delta + (e$nu.s0 %*% t(e$nu.s0)) / (1 / tau + K)
      Sigma.TZ.cand <- diag(e$sig2k - e$param[[g]]$nu.cand^2) + J*delta + (e$param[[g]]$nu.cand %*% t(e$param[[g]]$nu.cand)) / (1 / tau + K)

      # Independent proposals may lead to a matrix that is not (semi-)positive-definite
      chol.TZ.cand <- tryCatch({chol(Sigma.TZ.cand)}, error = function(er) { return(NULL) })
      if(is.null(chol.TZ.cand))
        reset <- TRUE
      chol.TZ <- tryCatch({chol(Sigma.TZ)}, error = function(er) { return(NULL) })
      if(is.null(chol.TZ))
        reset <- TRUE


      Sigma.ZT <- ( diag(K) + J*tau ) - diag(e$nu.s0) %*% solve( diag(e$sig2k) + J*delta ) %*% t(diag(e$nu.s0))
      Sigma.ZT.cand <- ( diag(K) + J*tau ) - diag(e$param[[g]]$nu.cand) %*% solve( diag(e$sig2k) + J*delta ) %*% t(diag(e$param[[g]]$nu.cand))

      #Independent proposals may lead to a matrix that is not (semi-)positive-definite
      chol.ZT.cand <- tryCatch({chol(Sigma.ZT.cand)}, error = function(er) { return(NULL) })
      if(is.null(chol.ZT.cand))
        reset <- TRUE
      chol.ZT <- tryCatch({chol(Sigma.ZT)}, error = function(er) { return(NULL) })
      if(is.null(chol.ZT))
        reset <- TRUE


      if(!reset) {

        # Partioned matrix components for each of the K items
        partMatrix.cand <- vector("list", e$K)
        for (k in 1:e$K)
        {
          ### MH candidates ###
          w.min1.cand <- e$param[[g]]$nu.cand[-k]/e$sig2k[-k]
          b.min1.cand <- e$sig2k[-k] / (e$sig2k[-k] - e$param[[g]]$nu.cand[-k]^2)
          c.min1.cand <- 1/tau + sum(b.min1.cand)
          A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
          d.min1.cand <- a.cand + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
          g.min1.cand <- sum(e$param[[g]]$nu.cand[-k] / (e$sig2k[-k] - e$param[[g]]$nu.cand[-k]^2)) / c.min1.cand
          A.min1.inv_w.cand <- (e$param[[g]]$nu.cand[-k] / (e$sig2k[-k] - e$param[[g]]$nu.cand[-k]^2) - g.min1.cand*b.min1.cand)

          partMatrix.cand[[k]]$B11 <- 1 + tau - e$param[[g]]$nu.cand[k]^2 * (1/e$sig2k[k] - ((1/e$sig2k[k])^2) / a.cand)
          partMatrix.cand[[k]]$B22 <- (1 - (e$param[[g]]$nu.cand[-k]^2) / e$sig2k[-k])*I.min1 + tau*J.min1 + w.min1.cand %*% t(w.min1.cand) / a.cand
          partMatrix.cand[[k]]$B12 <- tau*ones.min1 + ((e$param[[g]]$nu.cand[k] / e$sig2k[k]) %*% t(w.min1.cand)) / a.cand
          partMatrix.cand[[k]]$B21 <- t(partMatrix.cand[[k]]$B12)
          partMatrix.cand[[k]]$B22.inv <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]
        }

        # Conditional mean of Z|T for each item
        # x.cand <- mu.ZT <- mu.ZT.cand <- matrix(NA, ncol = K, nrow = Ng)
        # for (k in 1:K)
        # {
        #   x.cand[, k] <- t( (1/e$sig2k[k]) * ( t(e$RTg[[g]][, k] - muT[k]) - (t(1/e$sig2k) %*% (t(e$RTg[[g]]) - matrix(muT, nrow = K, ncol = Ng))) / a.cand ) )
        #   mu.ZT.cand[, k] <- muZ[k] + e$param[[g]]$nu.cand[k] * x.cand[, k]
        #   mu.ZT[, k] <- muZ[k] + e$nu.s0 * x.cand[, k]
        #   e$param[[g]]$uz[, k] <- e$param[[g]]$nu.cand[k] * x.cand[, k]
        # }

        #out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = e$param[[g]]$Z, mu.Z = e$param[[g]]$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
        #                       error = function(er) {  return(NULL) })

        # out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = out.Z.cand$Z, mu.Z = out.Z.cand$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
        #                        error = function(er) {  return(NULL) })
        #
        # out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = out.Z.cand$Z, mu.Z = out.Z.cand$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
        #                        error = function(er) {  return(NULL) })

        validProposals <- FALSE
        # if(!is.null(out.Z.cand)) {
          validProposals <- TRUE

          mu.TZ <- t(muT + diag(e$nu.s0) %*% (diag(K) - J/(1/tau + K)) %*% t(e$param[[g]]$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))
          #mu.TZ.cand <- t(muT + diag(e$param[[g]]$nu.cand) %*% (diag(K) - J/(1/tau + K)) %*% t(out.Z.cand$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))
          mu.TZ.cand <- t(muT + diag(e$param[[g]]$nu.cand) %*% (diag(K) - J/(1/tau + K)) %*% t(e$param[[g]]$Z - matrix(muZ, nrow = Ng, ncol = K, byrow = TRUE)))

          # Evaluate likelihood 1 vs 2 in Z|T
          lik <- lik.cand <- 0
          for (i in 1:Ng) {
            lik <- lik + sum(mvnfast::dmvn(X = e$RTg[[g]][i, ], mu = mu.TZ[i, ], sigma = chol.TZ, log = TRUE, isChol = TRUE))
            lik.cand <- lik.cand + sum(mvnfast::dmvn(X = e$RTg[[g]][i, ], mu = mu.TZ.cand[i, ], sigma = chol.TZ.cand, log = TRUE, isChol = TRUE))
            #lik <- lik + sum(mvnfast::dmvn(X = e$param[[g]]$Z[i, ], mu = e$param[[g]]$mu.ZT[i, ], sigma = chol.ZT, log = TRUE, isChol = TRUE))
            #lik.cand <- lik.cand + sum(mvnfast::dmvn(X = e$param[[g]]$Z[i, ], mu = mu.ZT.cand[i, ], sigma = chol.ZT.cand, log = TRUE, isChol = TRUE))
          }

          #cat(lik, "", lik.cand, "\n")

          ar <- exp(lik.cand - lik)
          if(is.nan(ar) || is.na(ar))
            ar <- 0
          else
            ar <- min(1, ar)
        # }
        # else {
        #   ar <- 0
        #   #browser()
        # }
        u <- runif (1, 0, 1)

        # print(e$nu.s0)
        # print(e$param[[g]]$nu.cand)
        # browser()


        if (u <= ar) {
          #cat("Accept tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
          e$chains[[c]][[g]][[9]][xg, ] <- nu <- e$param[[g]]$nu.cand
          #e$param[[g]]$mu.Z <- out.Z.cand$mu.Z
          #e$param[[g]]$var.Z <- out.Z.cand$var.Z
          #e$param[[g]]$Z <- out.Z.cand$Z
          #e$param[[g]]$mu.ZT <- mu.ZT.cand

          #or(k in 1:K)
          #  e$param[[g]]$ut[, k] <- mu.TZ.cand[, k] - muT[k]

          #for (k in 1:K)
          #  e$param[[g]]$r[, k] <- (out.Z.cand$Z[, k]- muZ[k]) - out.Z.cand$tmp[[k]][1, ]
          #e$param[[g]]$x <- x.cand

          e$mh.accept[g, 3] <- e$mh.accept[g, 3] + 1
          e$mh.accept.last[g, 3] <- e$mh.accept.last[g, 3] + 1

        }
        else {
          #print("Reject tau")
          #cat("Reject tau|delta|nu: ", e$param[[g]]$tau.cand, "|", e$param[[g]]$delta.cand, "|", e$param[[g]]$nu.cand, "\n")
          e$chains[[c]][[g]][[9]][xg, ] <- nu <- e$nu.s0
          #for(k in 1:K)
          #  e$param[[g]]$ut[, k] <- mu.TZ[, k] - muT[k]
        }
      }
    }
#
#     if(g == 1) {
#       e$chains[[c]][[g]][[7]][xg, ] <- dat1$sig2k
#       e$chains[[c]][[g]][[8]][xg, ] <- mean(dat1$sig2)
#       #e$chains[[c]][[g]][[9]][xg, ] <- dat1$nu
#     }
#     else {
#       e$chains[[c]][[g]][[7]][xg, ] <- dat2$sig2k
#       e$chains[[c]][[g]][[8]][xg, ] <- mean(dat2$sig2)
#       #e$chains[[c]][[g]][[9]][xg, ] <- dat2$nu
#     }


    if (reset || !validProposals)
      return(list(reset = reset, validProposals = validProposals))

    # Sample latent responses based on new set of parameters
    a.cand <- 1/delta + sum(1/e$sig2k)
    # Conditional mean of Z|T for each item
    x.cand <- mu.ZT.cand <- matrix(NA, ncol = K, nrow = Ng)
    for (k in 1:K)
    {
      x.cand[, k] <- t( (1/e$sig2k[k]) * ( t(e$RTg[[g]][, k] - muT[k]) - (t(1/e$sig2k) %*% (t(e$RTg[[g]]) - matrix(muT, nrow = K, ncol = Ng))) / a.cand ) )
      mu.ZT.cand[, k] <- muZ[k] + nu[k] * x.cand[, k]
      e$param[[g]]$uz[, k] <- nu[k] * x.cand[, k]
    }

    w.min1.cand <- nu[-k]/e$sig2k[-k]
    b.min1.cand <- e$sig2k[-k] / (e$sig2k[-k] - nu[-k]^2)
    c.min1.cand <- 1/tau + sum(b.min1.cand)
    A.min1.inv.cand <- b.min1.cand*I.min1 - (b.min1.cand %*% t(b.min1.cand)) / c.min1.cand
    d.min1.cand <- a.cand + t(w.min1.cand) %*% A.min1.inv.cand %*% w.min1.cand
    g.min1.cand <- sum(nu[-k] / (e$sig2k[-k] - nu[-k]^2)) / c.min1.cand
    A.min1.inv_w.cand <- (nu[-k] / (e$sig2k[-k] - nu[-k]^2) - g.min1.cand*b.min1.cand)

    partMatrix.cand[[k]]$B11 <- 1 + tau - nu[k]^2 * (1/e$sig2k[k] - ((1/e$sig2k[k])^2) / a.cand)
    partMatrix.cand[[k]]$B22 <- (1 - (nu[-k]^2) / e$sig2k[-k])*I.min1 + tau*J.min1 + w.min1.cand %*% t(w.min1.cand) / a.cand
    partMatrix.cand[[k]]$B12 <- tau*ones.min1 + ((nu[k] / e$sig2k[k]) %*% t(w.min1.cand)) / a.cand
    partMatrix.cand[[k]]$B21 <- t(partMatrix.cand[[k]]$B12)
    partMatrix.cand[[k]]$B22.inv <- A.min1.inv.cand - (A.min1.inv_w.cand %*% t(A.min1.inv_w.cand)) / d.min1.cand[1,1]


    out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = e$param[[g]]$Z, mu.Z = e$param[[g]]$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
                           error = function(er) {  return(NULL) })
    # for (k in 1:(K-1)) {
    #   out.Z.cand <- tryCatch({sampleZ(Y = e$Yg[[g]], Z = out.Z.cand$Z, mu.Z = out.Z.cand$mu.Z, mu.ZT = mu.ZT.cand, partMatrix.cand, likelihood = FALSE)},
    #                          error = function(er) {  return(NULL) })
    # }

    if(is.null(out.Z.cand)) {
      validProposals <- FALSE
    }
    else {
      e$param[[g]]$mu.Z <- out.Z.cand$mu.Z
      e$param[[g]]$var.Z <- out.Z.cand$var.Z
      e$param[[g]]$Z <- out.Z.cand$Z
      e$param[[g]]$mu.ZT <- mu.ZT.cand

      for (k in 1:K)
        e$param[[g]]$r[, k] <- (out.Z.cand$Z[, k]- muZ[k]) - out.Z.cand$tmp[[k]][1, ]
      e$param[[g]]$x <- x.cand
    }

    #e$chains[[c]][[g]][[7]][xg, ] <- data$sig2k
    #e$chains[[c]][[g]][[8]][xg, ] <- mean(data$sig2)
    #e$chains[[c]][[g]][[9]][xg, ] <- data$nu

    return(list(reset = reset, validProposals = validProposals))

  }

  ############################################################################################################

  e <- environment()
  e$initMar(e)
  e$initChains(reinit = FALSE, which = c(1, 2), e)
  #browser()

  ###### Run MCMC ######
  # For each chain
  for (c in 1:2) {

    chain <- chains[[c]] # Read the cc-th chain
    mh.accept <- matrix(0, nrow = G, ncol = 3) # MH acceptance rate 3 steps
    mh.accept.last <- matrix(0, nrow = G, ncol = 3) # MH acceptance rate last 100 iterations
    reset.count <- 0
    reinit.count <- 0

    xg <- 2
    while (xg <= XG) {

      ############################################## Gibbs sampling ##############################################

      ### beta and lambda do not differ between groups ###
      ### Assume equally sized groups
      beta.s0 <- chains[[c]][[1]][[1]][xg-1, ]
      lambda.s0 <- chains[[c]][[1]][[2]][xg-1, ]
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


      g <- 1
      while(g <= G) {
        beta.s0 <- chains[[c]][[g]][[1]][xg-1, ]
        lambda.s0 <- chains[[c]][[g]][[2]][xg-1, ]
        chains[[c]][[g]][[1]][xg, ] <- sampleBeta(Z = param[[g]]$Z, beta = beta.s0, theta = chains[[c]][[g]][[3]][xg-1, ], tau = chains[[c]][[g]][[5]][xg-1, ])
        chains[[c]][[g]][[2]][xg, ] <- sampleLambda(RT = RTg[[g]], lambda = lambda.s0, zeta = chains[[c]][[g]][[4]][xg-1, ], sig2k = chains[[c]][[g]][[7]][xg-1, ], delta = chains[[c]][[g]][[6]][xg-1, ])
        g <- g + 1
      }
# - e$param[[g]]$uz
      firstGroup <- TRUE
      reset <- FALSE
      g <- 1
      while((g <= G) && !reset) {

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
        beta <- chains[[c]][[g]][[1]][xg, ]
        lambda <- chains[[c]][[g]][[2]][xg, ]

        ### Sample group ability parameter ###
        if (g == 1) # Ability group mean of first group is fixed to zero
          chains[[c]][[g]][[3]][xg, ] <- theta <- 0
        else
          chains[[c]][[g]][[3]][xg, ] <- theta <- sampleTheta(Z = param[[g]]$Z, beta = chains[[c]][[1]][[1]][xg, ], tau = tau.s0)

        ### Sample group speed parameter ###
        if (g == 1) # Speed group mean of first group is fixed to zero
          chains[[c]][[g]][[4]][xg, ] <- zeta <- 0
        else
          chains[[c]][[g]][[4]][xg, ] <- zeta <- sampleZeta(RT = RTg[[g]], lambda = chains[[c]][[1]][[2]][xg, ], sig2k = sig2k.s0, delta = delta.s0)

        ### Sample measurement error variance parameters ###
        out.sig2 <- sampleSig2(RT = RTg[[g]], lambda = lambda, zeta = zeta, delta = delta.s0)
        e$chains[[c]][[g]][[7]][xg, ] <- sig2k <- out.sig2$sig2k
        e$chains[[c]][[g]][[8]][xg, ] <- sig2 <- out.sig2$sig2


        ############################################################################################################


        ################################################# Proposals ################################################

        e$sampleProposals(reinit = FALSE, g = g, e = e)

        ############################################################################################################


        ########################################### Metropolis-Hastings ############################################

        out.mh <- e$doMH(e)
        while (!out.mh$reset && !out.mh$validProposals && !reset) {
          e$sampleProposals(reinit = FALSE, g = g, e = e)
          out.mh <- e$doMH(e)
          reinit.count <- reinit.count + 1
          if(reinit.count == 100) {
            reset <- TRUE
            reinit.count <- 0
            #browser()
          }
          #browser()
          #print("Invalid Proposals")
        }

        ############################################################################################################

        if(out.mh$reset) {
          reset <- TRUE
          #print("Reset")
        }

        g <- g + 1
      }

      if(((any(mh.accept/xg < 0.01) || any(mh.accept.last/100 < 0.00)) && (xg%%100 == 0)) || reset) {
        #print(round(mh.accept.last/100, digits=2))

        cat("Reset | MH acceptance rate ", round(mh.accept/xg, digits=2), " | ", round(mh.accept.last/100, digits=2), "\n")

        if(!reset) {
          reset <- TRUE
        }
        reset.count <- reset.count + 1
        if(reset.count == 30) {
          return(NULL)
        }
        e$initMar(e)
        e$initChains(reinit = TRUE, which = c, e)
        e$mh.accept <- matrix(0, nrow = G, ncol = 3) # MH acceptance rate 3 steps
        e$mh.accept.last <- matrix(0, nrow = G, ncol = 3) # MH acceptance rate last 100 iterations
        xg <- 2
      }

      if ((!silent) && (xg%%100 == 0) && !reset) {
        cat("Iteration ", xg, " | MH acceptance rate ", round(mh.accept/xg, digits=2), " | ", round(mh.accept.last/100, digits=2), "\n")
        e$mh.accept.last <- matrix(0, nrow = G, ncol = 3)
      }

      if(!reset) {
        xg <- xg + 1
      }
    }
  }

  # Number of burnin iterations
  XG.burnin <- floor(XG * burnin)

  post.means <- vector("list", G)
  samples <- vector("list", G)
  g <- 1
  while(g <= G) {
    post.means[[g]] <- list()

    # Merge chains
    chain.1 <- chains[[1]][[g]]
    chain.2 <- chains[[2]][[g]]
    post.means[[g]]$beta <- colMeans((chain.1[[1]][XG.burnin:XG, ] + chain.2[[1]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$lambda <- colMeans((chain.1[[2]][XG.burnin:XG, ] + chain.2[[2]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$theta <- mean((chain.1[[3]][XG.burnin:XG, ] + chain.2[[3]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$zeta <- mean((chain.1[[4]][XG.burnin:XG, ] + chain.2[[4]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$tau <- mean((chain.1[[5]][XG.burnin:XG, ] + chain.2[[5]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$delta <- mean((chain.1[[6]][XG.burnin:XG, ] + chain.2[[6]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$sig2k <- colMeans((chain.1[[7]][XG.burnin:XG, ] + chain.2[[7]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$sig2 <- mean((chain.1[[8]][XG.burnin:XG, ] + chain.2[[8]][XG.burnin:XG, ]) / 2)
    post.means[[g]]$nu <- colMeans((chain.1[[9]][XG.burnin:XG, ] + chain.2[[9]][XG.burnin:XG, ]) / 2)


    post.means[[g]]$beta.sd <- apply(rbind(chain.1[[1]][XG.burnin:XG, ], chain.2[[1]][XG.burnin:XG, ]), FUN = sd, MARGIN = 2)
    post.means[[g]]$lambda.sd <- apply(rbind(chain.1[[2]][XG.burnin:XG, ], chain.2[[2]][XG.burnin:XG, ]), FUN = sd, MARGIN = 2)
    post.means[[g]]$theta.sd <- sd(c(chain.1[[3]][XG.burnin:XG, ], chain.2[[3]][XG.burnin:XG, ]))
    post.means[[g]]$zeta.sd <- sd(c(chain.1[[4]][XG.burnin:XG, ], chain.2[[4]][XG.burnin:XG, ]))
    post.means[[g]]$tau.sd <- sd(c(chain.1[[5]][XG.burnin:XG, ], chain.2[[5]][XG.burnin:XG, ]))
    post.means[[g]]$delta.sd <- sd(c(chain.1[[6]][XG.burnin:XG, ], chain.2[[6]][XG.burnin:XG, ]))
    post.means[[g]]$sig2k.sd <- apply(rbind(chain.1[[7]][XG.burnin:XG, ], chain.2[[7]][XG.burnin:XG, ]), FUN = sd, MARGIN = 2)
    post.means[[g]]$sig2.sd <- sd(c(chain.1[[8]][XG.burnin:XG, ], chain.2[[8]][XG.burnin:XG, ]))
    post.means[[g]]$nu.sd <- apply(rbind(chain.1[[9]][XG.burnin:XG, ], chain.2[[9]][XG.burnin:XG, ]), FUN = sd, MARGIN = 2)

    post.means[[g]]$beta.mce <- post.means[[g]]$beta.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$lambda.mce <- post.means[[g]]$lambda.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$theta.mce <- post.means[[g]]$theta.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$zeta.mce <- post.means[[g]]$zeta.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$tau.mce <- post.means[[g]]$tau.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$delta.mce <- post.means[[g]]$delta.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$sig2k.mce <- post.means[[g]]$sig2k.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$sig2.mce <- post.means[[g]]$sig2.sd / sqrt(2*(XG - XG.burnin))
    post.means[[g]]$nu.mce <- post.means[[g]]$nu.sd / sqrt(2*(XG - XG.burnin))

    samples[[g]] <- data.frame("beta" = rbind(chain.1[[1]][, ], chain.2[[1]][, ]),
                               "lambda" = rbind(chain.1[[2]][, ], chain.2[[2]][, ]),
                               "theta" = c(chain.1[[3]][, ], chain.2[[3]][, ]),
                               "zeta" = c(chain.1[[4]][, ], chain.2[[4]][, ]),
                               "tau" = c(chain.1[[5]][, ], chain.2[[5]][, ]),
                               "delta" = c(chain.1[[6]][, ], chain.2[[6]][, ]),
                               "sig2k" = rbind(chain.1[[7]][, ], chain.2[[7]][, ]),
                               "sig2" = c(chain.1[[8]][, ], chain.2[[8]][, ]),
                               "nu" = rbind(chain.1[[9]][, ], chain.2[[9]][, ]),
                               "chain" = c(rep(1, XG), rep(2, XG))
    )

    g <- g + 1
  }

  ###### Model selection ######
  if (G > 1) {

    # Model 1: G speed group parameters vs Model 2: 1 speed group parameter
    zetag <- deltag <- numeric(G)
    sig2gk <- list()

    g <- 1
    while (g <= G) {
      zetag[g] <- post.means[[g]]$zeta
      sig2gk[[g]] <- post.means[[g]]$sig2k
      deltag[g] <- post.means[[g]]$delta
      g <- g + 1
    }
    out.BIC.zeta <- computeBIC.zeta(RTg, post.means[[1]]$lambda, zetag, sig2gk, deltag)
  }
  else {
    out.BIC.zeta <- NULL
  }


  #print(param[[1]]$beta.mar)
  return(list(post.means = post.means, samples = samples, XG.burnin = XG.burnin, BIC.zeta = out.BIC.zeta))
}
