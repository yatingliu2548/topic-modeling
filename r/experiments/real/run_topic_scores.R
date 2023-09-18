source("r/score.R")

run_topic_models <- function(X, train_index, #test_index, 
                             list_params=1:20, threshold = FALSE,
                             normalize="norm", VHMethod = 'SP',
                             Mquantile = 0.05,
                             alpha = 0.005, max_K=150, K0=NULL){
  ####
  active_train = which(apply(X[train_index,], 2, sum)>0)
  x_train = t(diag(1/ apply(X[train_index, active_train],1, sum)) %*% X[train_index, active_train])
  topic_models <- vector("list", length(list_params))  
  it = 1
  for (k in list_params){

    tm <- score(as.matrix(x_train), k, scatterplot=FALSE, 
                K0=NULL, m=NULL, N=median(apply(X, 1, sum)), threshold=threshold,
                Mquantile=Mquantile, VHMethod = VHMethod, normalize=normalize,
                alpha=alpha, max_K=max_K, returnW=TRUE, estimateK=FALSE)
    A_hat = matrix(0, ncol(X), k)
    A_hat[active_train, ] = tm$A_hat
    topic_models[[it]] <- list(
      beta = log(t(A_hat)) %>% magrittr::set_colnames(colnames(X)),
      gamma =  t(tm$W_hat) %>% magrittr::set_rownames(rownames(X)[train_index])
    )
    it <- it + 1
  }
  names(topic_models) <- sapply(list_params, function(x){paste0("k", x)})
  #names(topic_models_test) <- sapply(list_params,function(x){paste0("k",x)})
  return(topic_models)
  
}