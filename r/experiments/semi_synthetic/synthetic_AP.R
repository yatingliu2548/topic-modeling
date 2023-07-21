library(tm)
library(topicmodels)
library(R.matlab)
library(tidyverse)
library(reticulate)
library(tidytext)

setwd("~/Documents/topic-modeling/")
source("./vertex_hunting_functions.R")
source('./score.R')
source('./evaluation_metrics.r')

use_condaenv("r-reticulate")

synthetic_dataset_generation <- function(dataset, K, doc_length=100, n=100, seed = 1234){
  # p is the number of words in the dictionary.
  # n is the number of documents.
  # N is a vector of length n, with ith entry as the total number of words in the ith documen.
  
  # read from sparse matrix txt file 
  if (dataset=="AP"){
    data("AssociatedPress")
    ap_td <- tidy(AssociatedPress)
    # ap_td %>%
    #   count(term, sort = TRUE) 
    data(stop_words)
    ap_td <- ap_td %>%
      anti_join(stop_words %>% rename(term=word)) ### removes the stopwords
    D = ap_td %>%
      cast_tdm(document, term, count)
  }else{
    print("Dataset not implemented yet")
  }
  vocab = colnames(D)
  ap_lda <- LDA(D, k = K, control = list(seed = seed))
  ap_topics <- tidy(ap_lda, matrix = "beta")
  
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
  
  D_sim <- ap_lda@gamma%*%exp(ap_lda@beta)
  p = dim(D_sim)[2]
  D_synth <- matrix(0, n, p)
  for (i in 1:length(N)){
    D_synth[i,] <- rmultinom(1, N[i], D_sim[i,])
  }
  counts = apply(D_synth, 2, sum)
  words2remove = which(counts == 0)
  return(list(D=D_synth[, -words2remove], vocab=vocab[-words2remove], A=exp(t(ap_lda@beta))[-words2remove, ], W = ap_lda@gamma))
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

run_experiment <- function(dataset, K, N=500, n=100, seed = 1234){
  
  data = synthetic_dataset_generation(dataset,  K, doc_length=500, n=100, seed = 1234)
  p = dim(data$D)[2]
  
  
  #### Step 1: Run the LDA pipeline
  lda <- LDA(data$D, k = K, control = list(seed = seed), method = 'VEM')
  ap_topics <- tidy(lda, matrix = "beta")
  Ahat_lda = exp(t(lda@beta))
  What_lda = lda@gamma
  
  resultsA <- process_results(Ahat_lda, "LDA", data$vocab)
  resultsW <- process_results(What_lda, "LDA", seq_len(n), processingA=FALSE)
  error <- update_error(Ahat_lda, What_lda, data$A, data$W, method = "LDA", error=NULL)
  
  #### Step 2: Run Tracy's method
  score_recovery <- score(t(data$D), K, normalize = "norm")
  resultsA <- rbind(resultsA, 
                    process_results(score_recovery$A_hat, "TopicScore", data$vocab))
  resultsW <- rbind(resultsW,
                    process_results(score_recovery$W_hat, "TopicScore", seq_len(n), processingA=FALSE))

  ### Step 3: Run AWR
  source_python("alternative/anchor-word-recovery/learn_topics_function.py")
  awr_recovery <- learn_topics(M=r_to_py(t(data$D)), vocab=data$vocab, 
                       K=K, loss='L2', save2txt=FALSE,
                       new_dim=1000, top_words=10,
                       vocab_file=NULL, outfile=NULL, anchor_thresh=0,
                       log_prefix='log', checkpoint_prefix='checkpoint',
                       max_threads=0, eps=1e-7)
  Ahat_awr = awr_recovery[[1]]
  resultsA <- rbind(resultsA, 
                    process_results(Ahat_awr, "AWR", data$vocab))
  What_awr <- compute_W_from_AD(Ahat_awr, D)
  error <- update_error(Ahat_awr, What_awr, data$A, data$W, method = "AWR", error=error)
  
  
  
  ### Step 4: Run TSVD
  
  id = ceiling(runif(n=1, max=1e6))
  inputPath <- paste0("/users/cdonnat/Documents/TopicSCORE/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr.mat")
  writeMat(inputPath, arr = data$D)
  outputPath <-paste0("/users/cdonnat/Documents/TopicSCORE/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr-result.mat")
  # Construct the MATLAB command
  command <- sprintf("/Applications/MATLAB_R2023a.app/bin/matlab -nodisplay -nosplash -r \"addpath('/users/cdonnat/Documents/TopicSCORE/alternative/TSVD'); TSVD('%s', '%s', %d); exit;\"", inputPath, outputPath, K)
  system(command)
  Ahat_tsvd <- readMat(outputPath)
  resultsA <- rbind(resultsA, 
                    process_results(Ahat_tsvd$M.hat, "TSVD", data$vocab))
  file.remove(inputPath)
  file.remove(outputPath)
  What_tsvd <- compute_W_from_AD(Ahat_tsvd, D)
  error <- update_error(Ahat_tsvd, What_tsvd, data$A, data$W, method = "TSVD", error=error)
  
  
  ### Step 5: Run Bing's method
  source_python("alternative/Sparse-topic-modelling/code/Sp_Top.py")
  bing_recovery <- Sp_Top(r_to_py(data$D), anchor_group =  r_to_py(c()), C0 = 0.01, C1 = r_to_py(c(1.1)), cv_rep = 50)
  resultsA <- rbind(resultsA, 
                    process_results(bing_recovery, "Bing", data$vocab))
  What_bing <- compute_W_from_AD(Ahat_bing, D)
  error <- update_error(Ahat_bing, What_bing, data$A, data$W, method = "Bing", error=error)
  
  
  
  #### Step 6: Run method
  score_ours <- score(t(data$D), K, normalize = "threshold")
  resultsA <- rbind(resultsA, 
                    process_results(score_recovery$A_hat, "Ours", data$vocab))
  resultsW <- rbind(resultsW,
                    process_results(score_recovery$W_hat, "Ours", seq_len(n), processingA=FALSE))
  error <- update_error(score_ours$A_hat, score_ours$W_hat, data$A, data$W, method = "Ours", error=error)
  
  return(list(resultsA=resultsA, resultsW=resultsW,
              error=error))
  
}
