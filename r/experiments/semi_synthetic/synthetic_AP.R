library(tm)
library(topicmodels)
library(R.matlab)
library(tidyverse)
library(reticulate)
library(tidytext)

setwd("~/Documents/topic-modeling/")
source("r/vertex_hunting_functions.R")
source('r/score.R')
source('r/evaluation_metrics.r')

use_condaenv("r-reticulate")

synthetic_dataset_generation <- function(dataset, K, doc_length=100, n=100, seed = 1234,
                                         A = NULL, W = NULL, vocab=NULL){
  # p is the number of words in the dictionary.
  # n is the number of documents.
  # N is a vector of length n, with ith entry as the total number of words in the ith documen.
  
  # read from sparse matrix txt file 
  if (dataset=="AP"){
    data("AssociatedPress")
    selected_docs = sample(seq_len(dim(AssociatedPress)[1]), n)
  }else{
    print("Dataset not implemented yet")
    return()
  }
  

  if (is.null(A) || is.null(W) ||is.null(vocab)){
    if (dataset=="AP"){
      ap_td <- tidy(AssociatedPress)
      # ap_td %>%
      #   count(term, sort = TRUE) 
      data(stop_words)
      ap_td <- ap_td %>%
        anti_join(stop_words %>% rename(term=word)) ### removes the stopwords
      D = ap_td %>%
        cast_tdm(document, term, count)
      vocab = colnames(D)
    }else{
      print("Dataset not implemented yet")
    }
    
    ap_lda <- LDA(D, k = K, control = list(seed = seed))
    ap_topics <- tidy(ap_lda, matrix = "beta")
    D_sim <- ap_lda@gamma[selected_docs,]%*%exp(ap_lda@beta)
    A = exp(t(ap_lda@beta))
    W = ap_lda@gamma
  }else{
    D_sim <- W[selected_docs,] %*% A
  }

  if (length(doc_length)==1){
    N <- rep(doc_length, n)
  }else{
    N <- doc_length
  }
  # ap_top_terms <- ap_topics %>%
  #   group_by(topic) %>%
  #   slice_max(beta, n = 10) %>% 
  #   ungroup() %>%
  #   arrange(topic, -beta)
  # 
  # ap_top_terms %>%
  #   mutate(term = reorder_within(term, beta, topic)) %>%
  #   ggplot(aes(beta, term, fill = factor(topic))) +
  #   geom_col(show.legend = FALSE) +
  #   facet_wrap(~ topic, scales = "free") +
  #   scale_y_reordered()
  p = dim(D_sim)[2]
  D_synth <- matrix(0, n, p)
  for (i in 1:length(N)){
    D_synth[i,] <- rmultinom(1, N[i], D_sim[i,])
  }
  counts = apply(D_synth, 2, sum)
  words2remove = which(counts == 0)
  return(list(D=D_synth[, -words2remove], vocab=vocab[-words2remove], 
              A = A[-words2remove, ], W =W[selected_docs,],
              Aoriginal=A, Woriginal = W, original_vocab = vocab))
}


process_results <- function(Ahat, method, vocab, processingA=TRUE){
  K = dim(Ahat)[2]
  temp = as_tibble(t(Ahat))
  colnames(temp) = vocab
  temp["topic"] = 1:K
  if(processingA){
    temp = pivot_longer(temp, cols = -c("topic")) %>% rename(word=name)
  }else{
    temp = pivot_longer(temp, cols = -c("topic")) %>% rename(document=name)
  }
  
  temp["method"] = method
  return(temp)
}


update_error <- function(Ahat, What, A, W, method, error){
  error_temp <- data.frame(l1_A=matrix_lp_distance(Ahat, A, lp=1),
                           l2_A=matrix_lp_distance(Ahat, A, lp=2),
                           l1_W=matrix_lp_distance(What, W, lp=1),
                           l2_W=matrix_lp_distance(What, W, lp=2),
                           method = method)
  if (is.null(error)){
    return(error_temp)
  }else{
    return(rbind(error,
                 error_temp))
  }
  return(error)
}

select_K <- function(){
  
}

run_experiment <- function(dataset, K, N=500, n=100, seed = 1234,
                           A=NULL, W=NULL, vocab=NULL){
  
  data = synthetic_dataset_generation(dataset,  K, doc_length=N, n=n, seed = seed,
                                    A=A, W=W, vocab=vocab)
  p = dim(data$D)[2]
  
  
  #### Step 1: Run the LDA pipeline
  lda <- LDA(data$D, k = K, control = list(seed = seed), method = 'VEM')
  ap_topics <- tidy(lda, matrix = "beta")
  Ahat_lda = exp(t(lda@beta))
  What_lda = lda@gamma
  
  resultsA <- process_results(Ahat_lda, "LDA", data$vocab)
  resultsW <- process_results(What_lda, "LDA", seq_len(n), processingA=FALSE)
  error <- update_error(Ahat_lda, What_lda, data$A, (data$W), method = "LDA", error=NULL)
  
  #### Step 2: Run Tracy's method
  score_recovery <- score(t(data$D), K, normalize = "norm")
  resultsA <- rbind(resultsA, 
                    process_results(score_recovery$A_hat, "TopicScore", data$vocab))
  resultsW <- rbind(resultsW,
                    process_results(score_recovery$W_hat, "TopicScore", seq_len(n), processingA=FALSE))
  error <- update_error(score_recovery$A_hat, t(score_recovery$W_hat), data$A, data$W, method = "TopicScore", error=error)

  ### Step 3: Run AWR
  source_python("r/alternative/anchor-word-recovery/learn_topics_function.py")
  awr_recovery <- learn_topics(M=r_to_py(t(data$D)), vocab=data$vocab, 
                       K=K, loss='L2', save2txt=FALSE,
                       new_dim=1000, top_words=10,
                       vocab_file=NULL, outfile=NULL, anchor_thresh=0,
                       log_prefix='log', checkpoint_prefix='checkpoint',
                       max_threads=0, eps=1e-7)
  Ahat_awr = awr_recovery[[1]]
  resultsA <- rbind(resultsA, 
                    process_results(Ahat_awr, "AWR", data$vocab))
  What_awr <- compute_W_from_AD(Ahat_awr, t(data$D))
  error <- update_error(Ahat_awr, What_awr, data$A, t(data$W), method = "AWR", error=error)
  
  
  
  ### Step 4: Run TSVD
  
  id = ceiling(runif(n=1, max=1e6))
  inputPath <- paste0(getwd(), "/r/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr.mat")
  writeMat(inputPath, arr = t(data$D))
  outputPath <-paste0(getwd(),"/r/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr-result.mat")
  # Construct the MATLAB command
  command <- sprintf("/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r \"addpath('%s/r/alternative/TSVD'); TSVD('%s', '%s', %d); exit;\"", getwd(), inputPath, outputPath, K)
  system(command)
  Ahat_tsvd <- readMat(outputPath)
  resultsA <- rbind(resultsA, 
                    process_results(Ahat_tsvd$M.hat, "TSVD", data$vocab))
  file.remove(inputPath)
  file.remove(outputPath)
  What_tsvd <- compute_W_from_AD(Ahat_tsvd$M.hat, t(data$D))
  error <- update_error(Ahat_tsvd$M.hat, What_tsvd, data$A, t(data$W), method = "TSVD", error=error)
  
  
  ### Step 5: Run Bing's method

    # Code that might throw an error
    source_python("r/alternative/Sparse-topic-modelling/code/Sp_Top.py")
    bing_recovery <- Sp_Top(r_to_py(t(data$D)), anchor_group =  r_to_py(c()), C0 = 0.1, C1 = r_to_py(c(1.1)), cv_rep = 50)
    resultsA <- rbind(resultsA, 
                      process_results(bing_recovery$A, "Bing", data$vocab))
    #### Have to cluster
    if (dim(t(bing_recovery$A))[2]>K){
      clustered_res <- kmeans(t(bing_recovery$A), centers = K) 
      What_bing <- compute_W_from_AD(t(clustered_res$centers), t(data$D))
    }else{
      What_bing <- compute_W_from_AD(bing_recovery$A, t(data$D))
      What_bing <- cbind(What_bing, matrix(0, nrow=nrow(What_bing), ncol = K -ncol(What_bing)  ))
    }
    
    error <- update_error(t(clustered_res$centers), What_bing, data$A, t(data$W), method = "Bing", error=error)
    

  
  
  
  #### Step 6: Run method
  for (alpha in c(0.1, 1, 2, 4, 8)){
    score_ours <- score(D = t(data$D), K=K, normalize = 'huy', 
                        threshold =TRUE, alpha = alpha, N=500)
    resultsA <- rbind(resultsA, 
                      process_results(score_ours$A_hat, "Ours", data$vocab))
    resultsW <- rbind(resultsW,
                      process_results(score_ours$W_hat, "Ours", seq_len(n), processingA=FALSE))
    error <- update_error(score_ours$A_hat, score_ours$W_hat, data$A, t(data$W), method = paste0("Ours_", alpha), error=error)
  }

  
  return(list(resultsA=resultsA, resultsW=resultsW,
              error=error, A=data$A, W=data$W, Aoriginal=data$Aoriginal, 
              Woriginal=data$Woriginal, vocab=data$original_vocab))
  
}


error <- c()
for (exp in 1:100){
  for (K in c(3, 4, 5, 7, 10, 20, 30, 50, 100)){
    A = NULL
    W = NULL
    vocab = NULL
    for (n in c(seq(from=100, to =1000, length.out=5), 2000, 5000)){
      for (n_frac in c(0.5, 0.8, 1, 2, 5, 10)){
        N = ceiling(n_frac * n)
        test <- run_experiment("AP", K, N=N, n=n, seed = (exp * 1000 +  n  + K  + 0.01 * n_frac)*1000, 
                               A = A, W=W, vocab=vocab)
        error_temp = test$error
        error_temp["K"] = K
        error_temp["N"] = N
        error_temp["n"] = n
        error_temp["seed"] =  (exp * 1000 +  n  + K  + 0.01 * n_frac)*1000
        error_temp["N"] = N
        error_temp["n_frac"] = n_frac
        error_temp["exp"] = exp
        error <- rbind(error,
                       error_temp)
        write_csv(error, paste0(getwd(), "/r/experiments/semi_synthetic/results_semisynthetic_right_K.csv"))
        if (is.null(A)){
          A = test$Aoriginal
          W = test$Woriginal
          vocab=test$original_vocab
        }
      }
    }
  }
}


