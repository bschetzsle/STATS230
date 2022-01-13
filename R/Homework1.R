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
