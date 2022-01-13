#'Multiply three matrices
#'
#' @param A a M x N matrix
#' @param B a N x P matrix
#' @param x a P x 1 vector/matrix
#' @param method specify the order of multiplication
#' @return If method=1, returns (AB)x, otherwise returns A(Bx)
#' @examples
#' matrix_prod(A, B, x, 1))
#' matrix_prod(A, B, x, 2))

matrix_prod = function(A, B, x, method=1){
  if(dim(A)[2] != dim(B)[1] || dim(B)[2] != dim(x)[1]){
    print("Non-conforming matrices!!")
    return(0)
  }
  if(method==1){
    return( (A %*% B) %*% x )
  }
  return( A %*% (B %*% x) )
}
