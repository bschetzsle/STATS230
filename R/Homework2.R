#'Simulate Multivariate Normal using Cholesky Decomposition of Covariance Matrix
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
  #get the cholesky decomposition of Sigma
  L = chol(Sigma)
  #generate N samples of multivariate normals
  Z = sapply(1:N, function(x){rnorm(n,0,1)})
  #transform Z to get desired multivariate normal samples
  n = dim(mu)[1]
  L = chol(Sigma)
  Z = sapply(1:N, function(x){rnorm(n,0,1)})
  X = apply(Z, 2, function(z){mu + t(L) %*% z})
  return(t(X))
}

#'Estimate OLS regression coefficients using QR decomposition
#'
#' @param Y the response vector
#' @param X the coveriate matrix, including intercept
#' @return Returns a vector of regression coefficients
QR_regression = function(Y, X){
  QR = qr(X)
  Q = qr.Q(QR)
  R = qr.R(QR)
  QR_beta = solve(R) %*% t(Q) %*% Y
  return(QR_beta)
}

#'Estimate OLS regression coefficients using SVD decomposition
#'
#' @param Y the response vector
#' @param X the coveriate matrix, including intercept
#' @return Returns a vector of regression coefficients
SVD_regression = function(Y, X){
  SVD = svd(X)
  SVD_beta = SVD$v %*% diag(1/SVD$d) %*% t(SVD$u) %*% Y
  return(SVD_beta)
}
