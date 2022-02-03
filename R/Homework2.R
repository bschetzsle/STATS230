#'Multiply three matrices
#'
#' @param N the number of samples to generate from the MVN
#' @param mu an nx1 mean vector
#' @param Sigma an nxn positive definite covariance matrix
#' @return The function returns an Nxn matrix of samples;
#' the rows correspond to individual samples
#' @examples
#' matrix_prod(N=10, mu, Sigma))

cholesky_mvn = function(N, mu, Sigma){

  n = dim(mu)[1]
  L = chol(Sigma)
  Z = sapply(1:N, function(x){rnorm(n,0,1)})
  X = apply(Z, 2, function(z){mu + t(L) %*% z})
  return(t(X))
}

