library(combinat)
library(clue)

error1_A <- function(A, A_hat){
  K <- dim(A)[2]
  all_perm <- permn(1:K)
  error <- Inf
  
  for (i in 1:length(all_perm)){
    error <- min(error, mean(colSums(abs(A[,all_perm[[i]]]-A_hat))))
  }
  
  return(error)
}

error2_A <- function(A, A_hat){
  K <- dim(A)[2]
  used <- rep(1,K)
  A_perm <- matrix(0,dim(A)[1],dim(A)[2])
  
  for (k in 1:K){
    dis <- colSums(abs(A-A_hat[,k]))*used
    index <- which(dis == min(dis))
    index <- index[1]
    A_perm[,k] <- A[,index]
    used[index] <- Inf
  }
  
  return(mean(colSums(abs(A_perm-A_hat))))
}


matrix_lp_distance <- function(A, B, lp=2){
  K <- dim(A)[2]
  if (lp ==2){
    error_matrix <- outer(seq_len(ncol(A)), seq_len(ncol(B)), Vectorize(function(i, j) l2_error(A[, i], B[, j])))
  }else{
    error_matrix <- outer(seq_len(ncol(A)), seq_len(ncol(B)), Vectorize(function(i, j) l1_error(A[, i], B[, j])))
  }
  # Find the optimal column permutation using the Hungarian algorithm
  permutation <- solve_LSAP(error_matrix) 
  B_permuted <- B[, permutation]

  if (lp ==2){
    error <- l2_error(A, B_permuted)
  }else{
    error <- l1_error(A, B_permuted)
  }
  
  return(error)
}

l2_error <- function(A, B) {
  sqrt(sum((A - B)^2))
}

l1_error <- function(A, B) {
  sum(abs(A - B))
}
