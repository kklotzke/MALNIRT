#' @export
helmert <- function(p) {
  H <- matrix(0, p, p)
  diag(H) <- -(0:(p-1)) * (-((0:(p-1))*((0:(p-1))+1))^(-0.5))
  for (i in 2:p){
    H[i,1:(i-1)]<- -((i-1)*(i))^(-0.5)
  }
  H[1,]<-1/sqrt(p)

  return(H)
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
