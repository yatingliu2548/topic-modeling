library(VGAM)
library(tidyverse)
#library(ggtern)
library(tm)
library(topicmodels)
library(R.matlab)
library(tidyverse)
library(reticulate)
library(tidytext)

source("r/experiments/semi_synthetic/synthetic_AP.R")
#### the purpose here is to generate synthetic examples

synthetic_dataset_creation <- function(n, K, p, alpha=0.5, n_max_zipf=5 * 1e5, a_zipf=1,
                                        n_anchors=0, delta_anchor=1, N=500, seed=123){
  
  set.seed(seed)
  W <- rdiric(n, rep(alpha, K)) 
  #ggtern(data = data.frame(W), aes(x = X1, y = X2, z = X3)) +
  #  geom_point(size = 3) +
  #  theme_bw()
  # 
  if (n_anchors >0){
    A = matrix(0, nrow=K, ncol=p)
    for (k in 1:K){
      A[k, ((k-1)*n_anchors +1) : (k * n_anchors)] = delta_anchor
    }
    A[,(K * n_anchors+1):p ] <- sapply(1/(1: (p-n_anchors * K) + 2.7)^a_zipf, function(u){rexp(K, u)})
  }else{
    A <- sapply(1/(1:p + 2.7)^a_zipf, function(u){rexp(K, u)}) 
  }
  A = t(A/apply(A, 1, sum))
  D0 = A %*% t(W)
  X <- sapply(1:n, function(i){rmultinom(1, N, D0[,i])})
  return(list(D=t(X[which(apply(X,1, sum) >0 ),]),
              A= A[which(apply(X,1, sum) >0 ),], 
              W = W, vocab =which(apply(X,1, sum) >0 ),
              D0=t(D0[which(apply(X,1, sum) >0 ),]) ))
}

#anchor_word+weak lq ball
synthetic_dataset_creation_2 <- function(n, K, p, alpha=0.5, n_max_zipf=5 * 1e5, a_zipf=1,
                                        n_anchors=0, delta_anchor=1, N=500, s=20, seed=123){
  
  set.seed(seed)
  W <- rdiric(n, rep(0.5, K)) #n * K matrix
  # ggtern(data = data.frame(W), aes(x = X1, y = X2, z = X3)) +
  #   geom_point(size = 3) +
  #   theme_bw()
  # 
  A <- matrix(0, nrow=K, ncol=p)
  if (s<p){
    A[,1:s] <- abs(rnorm(s * K, mean = 0, sd = 100))
    A[,(s + 1):p] <- t(rdiric(p-s,rep(alpha,K))) * 0.00001
  }else{
    A = t(rdiric(p,rep(eta,K))) 
  }
  if (n_anchors >0){ #for all 1<=k \neq j <=K, there exists at least a column of A that can be represented as x_{kj 1}*e_k+ x
    kj=1
    for (k in 1:K){
      for (j in (k + 1):K) {
        for (i in ((kj-1)*n_anchors +1) : (kj * n_anchors)){
          A[, i] <- runif(1, min = 0, max = 1 / K) * diag(K)[, k] + runif(1, min = 0, max = 1 / K) * diag(K)[, j]
        }
        kj=kj+1
      }  
    }
  }
   A = t(A/apply(A, 1, sum))
  D0 = A %*% t(W)
  X <- sapply(1:n, function(i){rmultinom(1, N, D0[,i])})
  #X = X/N
  return(list(D=t(X[which(apply(X,1, sum) >0 ),]),
              A= A[which(apply(X,1, sum) >0 ),], 
              W = W, vocab =which(apply(X,1, sum) >0 ),
              D0=t(D0[which(apply(X,1, sum) >0 ),]) ))
}


run_synthetic_experiment <- function(n, K, p, alpha=0.5, a_zipf=1,
                                     n_anchors=0, delta_anchor=1, N=500,
                                     noise_generation = "uniform", seed=1, VHMethod="SVS"){
  
  data = synthetic_dataset_creation(n, K, p, alpha=alpha, n_max_zipf=50000, a_zipf=a_zipf,
                                    n_anchors=n_anchors, delta_anchor=delta_anchor, 
                                    N=N, seed=seed)
  #### Run check
  print("here")
  print(sprintf("Dim of data D = %s, %s", dim(data$D)[1], dim(data$D)[2]))
  print(sprintf("Dim of data A = %s, %s", dim(data$A)[1], dim(data$A)[2]))
  
  # #### Step 1: Run the LDA pipeline
  lda <- LDA((data$D), k = K, control = list(seed = seed), method = 'VEM')
  ap_topics <- tidy(lda, matrix = "beta")
  Ahat_lda = exp(t(lda@beta))
  What_lda = lda@gamma
  
  # resultsA <- process_results(Ahat_lda, "LDA", data$vocab)
  # resultsW <- process_results(What_lda, "LDA", seq_len(n), processingA=FALSE)
  error <- update_error(Ahat_lda, t(What_lda), (data$A), t(data$W), method = "LDA", error=NULL,
                        thresholded = 0)
  
  #### Step 2: Run Tracy's method
  score_recovery <- score(t(data$D), K, normalize = "norm", max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod)
  Khat_tracy = select_K(score_recovery$eigenvalues, p, n, N, method="tracy")
  Khat_olga = select_K(svd(data$D)$d, p,n, N, method="olga")
  # resultsA <- rbind(resultsA, 
  #                   process_results(score_recovery$A_hat, "TopicScore", data$vocab))
  # resultsW <- rbind(resultsW,
  #                   process_results(score_recovery$W_hat, "TopicScore", seq_len(n), processingA=FALSE))
  
  error <- update_error(score_recovery$A_hat, (score_recovery$W_hat), data$A, t(data$W), method = "TopicScore", error=error,
                        thresholded = 0)
  
  ### Step 3: Run AWR
  Ahat_awr <- tryCatch(
    AWR(data),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
    error = function(err) {
      # Code to handle the error (e.g., print an error message, log the error, etc.)
      cat("Error occurred while running AWR:", conditionMessage(err), "\n")
      # Return a default value or NULL to continue with the rest of the code
      return(NULL)
    }
  )
  
  if(is.null(Ahat_awr) == FALSE){
    # resultsA <- rbind(resultsA, 
    #                   process_results(Ahat_awr, "AWR", data$vocab))
    What_awr <- compute_W_from_AD(Ahat_awr, t(data$D))
    # resultsW <- rbind(resultsW, 
    #                   process_results(What_awr, "AWR", seq_len(n), processingA=FALSE))
    error <- update_error(Ahat_awr, (What_awr), data$A, t(data$W), method = "AWR", error=error)
  }
  
  
  
  
  ### Step 4: Run TSVD
  
 # id = ceiling(runif(n=1, max=1e6))
 # resultTSVD <- tryCatch(
 #   TSVD(data, n, K, p, id, matlab_path=matlab_path),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
 #   error = function(err) {
      # Code to handle the error (e.g., print an error message, log the error, etc.)
 #     cat("Error occurred while running TSVD:", conditionMessage(err), "\n")
      # Return a default value or NULL to continue with the rest of the code
 #     return(NULL)
 #   }
 # )
 # if (is.null(resultTSVD) == FALSE){
    # resultsA <- rbind(resultsA, 
    #                   process_results(resultTSVD$Ahat$M.hat, "TSVD", data$vocab))
    # resultsW <- rbind(resultsW, 
    #                   process_results(resultTSVD$What, "TSVD", seq_len(n), processingA=FALSE))
 #   error <- update_error(resultTSVD$Ahat$M.hat, (resultTSVD$What), data$A, t(data$W), method = "TSVD", error=error)
 # }
  
  
  ### Step 5: Run Bing's method
  
  # Code that might throw an error
  bing_recovery <- tryCatch(
    Bing(data, C0=0.1, C1=1.1),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
    error = function(err) {
      # Code to handle the error (e.g., print an error message, log the error, etc.)
      cat("Error occurred while running Bing:", conditionMessage(err), "\n")
      # Return a default value or NULL to continue with the rest of the code
      return(NULL)
    }
  )
  
  if (is.null(bing_recovery) == FALSE){
    # resultsA <- rbind(resultsA, 
    #                   process_results(bing_recovery$A, "Bing", data$vocab))
    #### Have to cluster
    if (dim(t(bing_recovery$A))[1]>K){
      clustered_res <- kmeans(t(bing_recovery$A), centers = K) 
      What_bing <- compute_W_from_AD(t(clustered_res$centers), t(data$D))
      error <- update_error(t(clustered_res$centers), t(What_bing), data$A, (data$W), method = "Bing", error=error)
    }else{
      What_bing <- compute_W_from_AD(bing_recovery$A, t(data$D))
      What_bing <- rbind(What_bing, 
                         matrix(0, ncol=ncol(What_bing), nrow = K - nrow(What_bing)  ))
      error <- update_error((bing_recovery$A), (What_bing), data$A, t(data$W), method = "Bing", error=error)
      
    }
    # resultsW <- rbind(resultsW, 
    #                   process_results(What_bing, "Bing", seq_len(n), processingA=FALSE))
  }
  
  
  
  
  
  #### Step 6: Run method
  for (alpha in c(0.1, 0.5 , 1, 2, 4, 8)){
    score_ours <- tryCatch(
      score(D = t(data$D), K=K, normalize = 'huy', 
            threshold =TRUE, alpha = alpha, N=N, max_K = min(min(dim(data$D))-1, 150),
            VHMethod=VHMethod),
      error = function(err) {
        # Code to handle the error (e.g., print an error message, log the error, etc.)
        paste0("Error occurred while running Score ", alpha, " :", conditionMessage(err), "\n")
        # Return a default value or NULL to continue with the rest of the code
        return(NULL)
      }
    )
    if (is.null(score_ours) == FALSE){
      # resultsA <- rbind(resultsA, 
      #                   process_results(score_ours$A_hat, paste0("Ours_", alpha), data$vocab))
      # resultsW <- rbind(resultsW,
      #                   process_results(score_ours$W_hat, paste0("Ours_", alpha), seq_len(n), processingA=FALSE))
      error <- update_error(score_ours$A_hat, (score_ours$W_hat), data$A, t(data$W), method = paste0("Ours_", alpha), error=error,
                            thresholded=score_ours$thresholded)
      Khat_huy = select_K(score_ours$eigenvalues, p,n, N, method="huy")
    }
    
  }
  
  
  
  
  return(list(#resultsA=resultsA, resultsW=resultsW,
    error=error, A=data$A, W=data$W, Aoriginal=data$Aoriginal, 
    Woriginal=data$Woriginal, Epsilon=data$Epsilon, vocab=data$original_vocab,
    Khat_huy=Khat_huy$Khat,
    Khat_huy_thresh = Khat_huy$thresh,
    Khat_olga=Khat_olga$Khat,
    Khat_olga_thresh = Khat_olga$thresh,
    Khat_tracy=Khat_tracy$Khat,
    Khat_tracy_thresh = Khat_tracy$thresh))
  
}


