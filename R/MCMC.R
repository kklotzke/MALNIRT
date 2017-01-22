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
