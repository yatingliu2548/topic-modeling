library(Matrix)

simplex_dist <- function(theta, V){
  VV <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))%*%V
  D <- VV%*%t(VV)
  d <- VV%*%(theta-V[dim(V)[1],])
  
  A <- cbind(diag(rep(1,dim(V)[1]-1)), -rep(1,dim(V)[1]-1))
  b0 <- c(rep(0,dim(V)[1]-1),-1)
  
  # D <- matrix(nearPD(D)$mat, nrow(D), ncol(D))
  # D <- nearPD(D)
  obj <- solve.QP(D, d, A, b0)
  return(sum((theta-V[dim(V)[1],]) ^2)+ 2*obj$value)
}
