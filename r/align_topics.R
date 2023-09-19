
library(rdist)

normalize_rows <- function(mat) {
  row_norms <- sqrt(rowSums(mat^2))
  return(sweep(mat, 1, row_norms, FUN="/"))
}


align_topics<- function(A, B, dist="cosine", do.plot=TRUE){
  if(dist=="cosine"){
    A_normalized <- normalize_rows(A)
    B_normalized <- normalize_rows(B)
    match = A_normalized %*% t(B_normalized)
    match=  data.frame(match)
    permutation <- solve_LSAP(as.matrix(match), maximum=TRUE)
  }else{
    match = cdist(A, B, metric = "euclidean", p = 1)
    match=  data.frame(match)
    permutation <- solve_LSAP(as.matrix(match), maximum=FALSE)
  }

  match_permuted <- match[, permutation]
  if (do.plot){
    par(mar=c(1,1,1,1))
    colnames(match_permuted)= 1:ncol(match_permuted)
    match_permuted["X"] = 1:nrow(match_permuted)
    print(ggplot(pivot_longer(match_permuted, cols=-c("X")))+
      geom_tile(aes(x=X, y=name, fill=value)))
    #match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
    
  }

  B_permuted=B[permutation,]
  return(list("B_permuted"=B_permuted,
              "match" = match_permuted))
}
