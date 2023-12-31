library(tm)
library(topicmodels)
library(R.matlab)
library(tidyverse)
library(reticulate)
library(tidytext)

#setwd("~/topic-modeling/")
source(paste0(getwd(),"/r/vertex_hunting_functions.R"))
source(paste0(getwd(),'/r/score.r'))
source(paste0(getwd(),'/r/evaluation_metrics.r'))
source(paste0(getwd(),'/r/select_K.R'))

use_condaenv("r-reticulate")

DEFAULT_MATLAB="/Applications/MATLAB_R2023a.app/bin/matlab"

TSVD <- function(data, n, K, p, id, matlab_path=DEFAULT_MATLAB,
                 returnW=TRUE){
  inputPath <- paste0(getwd(), "/r/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr.mat")
  writeMat(inputPath, arr = t(data$D))
  outputPath <-paste0(getwd(),"/r/experiments/temporary/exp_","n", n, "-K", K, "-p",p , "-id", id,   "arr-result.mat")
  # Construct the MATLAB command
  command <- sprintf("%s -nodisplay -nosplash -r \"addpath('%s/r/alternative/TSVD'); TSVD('%s', '%s', %d); exit;\"", matlab_path, getwd(), inputPath, outputPath, K)
  system(command)
  Ahat_tsvd <- readMat(outputPath)
  file.remove(inputPath)
  file.remove(outputPath)
  if(returnW){
    What_tsvd <- compute_W_from_AD(Ahat_tsvd$M.hat, t(data$D))
  }else{
    What_tsvd=NULL
  }
  return(list(Ahat = Ahat_tsvd, 
              What = What_tsvd))
}


Cai <- function(data, C0=0.1, C1=1.1){
  source_python("r/alternative/Sparse-topic-modelling/Sp_Top.py")
  bing_recovery <- Sp_Top(r_to_py(t(data$D)), anchor_group =  r_to_py(c()), C0 = C0, C1 = r_to_py(c(C1)), cv_rep = 50)
  return(bing_recovery)
}


Bing <- function(data, C0=0.1, C1=1.1){
  source_python("r/alternative/Sparse-TM-Bing/Sp_Top.py")
  bing_recovery <- Sp_Top(r_to_py(t(data$D)), anchor_group =  r_to_py(c()), C0 = C0, C1 = r_to_py(c(C1)), cv_rep = 50)
  return(bing_recovery)
}


AWR <- function(data){
  source_python("r/alternative/AWR/learn_topics_function.py")
  awr_recovery <- learn_topics(M=r_to_py(t(data$D)), vocab=data$vocab, 
                               K=K, loss='L2', save2txt=FALSE,
                               new_dim=1000, top_words=10,
                               vocab_file=NULL, outfile=NULL, anchor_thresh=0,
                               log_prefix='log', checkpoint_prefix='checkpoint',
                               max_threads=0, eps=1e-7)
  Ahat_awr = awr_recovery[[1]]
  return(Ahat_awr)
}

synthetic_dataset_generation <- function(dataset, K, doc_length=100, n=100, seed = 1234,
                                         A = NULL, W = NULL, vocab=NULL, 
                                         remove_stop_words=TRUE,
                                         normalize_counts=FALSE){
  set.seed(seed)
  # p is the number of words in the dictionary.
  # n is the number of documents.
  # N is a vector of length n, with ith entry as the total number of words in the ith documen.
  
  # read from sparse matrix txt file 
  if (dataset=="AP"){
    data("AssociatedPress")
    if (n < dim(AssociatedPress)[1]){
      selected_docs = sample(seq_len(dim(AssociatedPress)[1]), n)
    }else{
      selected_docs = seq_len(dim(AssociatedPress)[1])
    }
    
  }else{
    print("Dataset not implemented yet")
    return()
  }
  

  if (is.null(A) || is.null(W) ||is.null(vocab)){
    if (dataset=="AP"){
      ap_td <- tidy(AssociatedPress)
      # ap_td %>%
      #   count(term, sort = TRUE) 
      if(remove_stop_words){
        data(stop_words)
         ap_td <- ap_td %>%
            anti_join(stop_words %>% rename(term=word)) ### removes the stopwords
      }
      D = ap_td %>%
        cast_tdm(document, term, count)
      vocab = colnames(D)
    }else{
      print("Dataset not implemented yet")
    }
    
    ap_lda <- LDA(D, k = K, control = list(seed = seed), method = 'VEM')
    ap_topics <- tidy(ap_lda, matrix = "beta")
    D_sim <- ap_lda@gamma[selected_docs,]%*%exp(ap_lda@beta)
    A = exp(t(ap_lda@beta))
    W = ap_lda@gamma
  }else{
    D_sim <- as.matrix(W)[selected_docs,] %*% t(as.matrix(A))
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
  D_synth <- t(sapply(1:length(N), function(i){rmultinom(1, N[i], D_sim[i,])}))
  D_synth[which(D_synth<0)] = 0

  print(sprintf("Dim of data D_synth = %s, %s", dim(D_synth)[1], dim(D_synth)[2]))
  counts = apply(D_synth, 2, sum)
  words2remove = which(counts < 1)
  print(sprintf("Dim of words2remove = %s", length(words2remove)))
  if (normalize_counts){
    D_synth = D_synth/mean(N)
  }
  if (length(words2remove) >0){
    Atemp = A[-words2remove, ]
    return(list(D=D_synth[, -words2remove], vocab=vocab[-words2remove], 
                A = Atemp %*% diag(1/apply(Atemp,2, sum)), W =W[selected_docs,],
                Aoriginal=A, Woriginal = W, original_vocab = vocab))
  }else{
    return(list(D=D_synth, vocab=vocab, 
                A = A, W =W[selected_docs,],
                Aoriginal=A, Woriginal = W, original_vocab = vocab))
  }
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


update_error <- function(Ahat, What, A, W, method, time, error, thresholded = 0){
  error_temp <- data.frame(l1_A=matrix_lp_distance(Ahat, A, lp=1),
                           l2_A=matrix_lp_distance(Ahat, A, lp=2),
                           l1_W=ifelse(is.null(What), NA, matrix_lp_distance(What, W, lp=1)),
                           l2_W=ifelse(is.null(What), NA, matrix_lp_distance(What, W, lp=2)),
                           thresholded = thresholded,
                           K = dim(A)[2],
                           p = dim(A)[1],
                           time = time,
                           method = method)
  
  if (is.null(error)){
    return(error_temp)
  }else{
    return(rbind(error,
                 error_temp))
  }
  return(error)
}

run_experiment <- function(dataset, K, N=500, n=100, seed = 1234,
                           A=NULL, W=NULL, vocab=NULL, plot_data=FALSE,
                           matlab_path=DEFAULT_MATLAB, 
                           VHMethod="SVS",
                           remove_stop_words=FALSE,
                           evaluateW=FALSE,
                           normalize_counts=TRUE){
  
  data = synthetic_dataset_generation(dataset,  K, doc_length=N, n=n, seed = seed,
                                    A=A, W=W, vocab=vocab, 
                                    remove_stop_words = remove_stop_words,
                                    normalize_counts=FALSE)
  #### Run check
  print("here")
  if (plot_data){
    ggplot(data.frame(count= apply(data$D,2,sum),
                      word=data$vocab)) +
      geom_histogram(aes(x=count)) +
      scale_x_log10()
    
    # Vertices of the simplex
    # library(ggtern)
    # ggtern(data = data.frame(data$W), aes(x = X1, y = X2, z = X3)) +
    #        geom_point(size = 3) +
    #        theme_bw()
    
    ### check the words

    norm = sqrt(apply(data$A^2,1,sum))
    sorted = sort(norm, index.return=T, decreasing=TRUE) 
    test= sorted$x * (1:length(sorted$x))#### Doesn't verify the assumption really
    print(max(test))
    res = data.frame(norm = sorted$x, ix = sorted$ix, x = 1:length(norm))
    ggplot(res, aes(x=x, y=norm * x)) + geom_point()
    
    counts = sort(apply(data$D,2,sum), decreasing=TRUE)
    ggplot(data.frame(x=log(1:length(counts)), c =log(counts))) + geom_line(aes(x=x, y=c)) + 
      geom_abline(aes(intercept =counts[1]-1, slope = -1), color = "blue")
    
    plot(sapply(1:length(norm), function(i){ i * norm[i]} ))
  }
  p = dim(data$D)[2]
  print(sprintf("Dim of data D = %s, %s", dim(data$D)[1], dim(data$D)[2]))
  print(sprintf("Dim of data A = %s, %s", dim(data$A)[1], dim(data$A)[2]))
  
  
  # #### Step 1: Run the LDA pipeline
  elapsed_timeLDA <- system.time({
    lda <- LDA(data$D, k = K, control = list(seed = seed), method = 'VEM')
    #ap_topics <- tidy(lda, matrix = "beta")
    Ahat_lda = exp(t(lda@beta))
  })["elapsed"]
  
  
  if(evaluateW){
    What_lda = lda@gamma
   error <- update_error(Ahat_lda, t(What_lda), data$A, t(data$W), 
                         time = elapsed_timeLDA, method = "LDA", error=NULL,
                        thresholded = 0)
  }else{
    What_lda = NULL
    error <- update_error(Ahat_lda, NULL, data$A, t(data$W), 
                          time = elapsed_timeLDA, method = "LDA", error=NULL,
                        thresholded = 0)
  }
  

  
  #### Step 2: Run Tracy's method
  elapsed_timeTracy <- system.time({
        score_recovery <- score(t(data$D)/N, K, normalize = "norm", max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod,
                                returnW=evaluateW)
  })["elapsed"]
  Khat_tracy = select_K(score_recovery$eigenvalues, p,n, N, method="tracy")
  Khat_olga = select_K(svd(data$D)$d, p,n, N, method="olga")
  
  error <- update_error(score_recovery$A_hat, (score_recovery$W_hat), 
                        data$A, t(data$W),  time = elapsed_timeTracy, 
                        method = "TopicScore", error=error,
                        thresholded = 0)

  ### Step 3: Run AWR
  elapsed_timeAWR <- system.time({
      Ahat_awr <- tryCatch(
        AWR(data),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
        error = function(err) {
          # Code to handle the error (e.g., print an error message, log the error, etc.)
          cat("Error occurred while running AWR:", conditionMessage(err), "\n")
          # Return a default value or NULL to continue with the rest of the code
          return(NULL)
        }
      )
  })["elapsed"]
  
  if(is.null(Ahat_awr) == FALSE){
    # resultsA <- rbind(resultsA, 
    #                   process_results(Ahat_awr, "AWR", data$vocab))
    if(evaluateW){
      What_awr <- compute_W_from_AD(Ahat_awr, t(data$D))
    }else{
      What_awr = NULL
    }
    # resultsW <- rbind(resultsW, 
    #                   process_results(What_awr, "AWR", seq_len(n), processingA=FALSE))
    error <- update_error(Ahat_awr, (What_awr), data$A, t(data$W), 
                          time = elapsed_timeAWR, method = "AWR", error=error)
  }

  
  
  
  ### Step 4: Run TSVD
  
#   id = ceiling(runif(n=1, max=1e6))
#   elapsed_timetsvd<- system.time({
#   resultTSVD <- tryCatch(
#     TSVD(data, n, K, p, id, matlab_path=matlab_path),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
#     error = function(err) {
#       # Code to handle the error (e.g., print an error message, log the error, etc.)
#       cat("Error occurred while running TSVD:", conditionMessage(err), "\n")
#       # Return a default value or NULL to continue with the rest of the code
#       return(NULL)
#     }
#   )
#   })["elapsed"]
#   if (is.null(resultTSVD) == FALSE){
#     error <- update_error(resultTSVD$Ahat$M.hat, (resultTSVD$What), data$A, t(data$W),
#                           time=elapsed_timetsvd, method = "TSVD", error=error)
#   }
  
  
  ### Step 5: Run Bing's method

    # Code that might throw an error
  elapsed_timeBing <- system.time({ 
          bing_recovery <- tryCatch(
          Bing(data, C0=0.1, C1=1.1),  # Replace arg1, arg2, ... with the actual arguments required by tSVD
          error = function(err) {
            # Code to handle the error (e.g., print an error message, log the error, etc.)
            cat("Error occurred while running Bing:", conditionMessage(err), "\n")
            # Return a default value or NULL to continue with the rest of the code
            return(NULL)
          }
        )   
    })["elapsed"]
  
  if (is.null(bing_recovery) == FALSE){
    # resultsA <- rbind(resultsA, 
    #                   process_results(bing_recovery$A, "Bing", data$vocab))
    #### Have to cluster
    if(evaluateW){
      if (dim(t(bing_recovery$A))[1]>K){
        clustered_res <- kmeans(t(bing_recovery$A), centers = K) 
        What_bing <- compute_W_from_AD(t(clustered_res$centers), t(data$D))
        error <- update_error(t(clustered_res$centers), (What_bing), data$A, t(data$W), 
                              time=elapsed_timeBing, method = "Bing", error=error,thresholded=bing_recovery$thresholded)
      }else{
        What_bing <- compute_W_from_AD(bing_recovery$A, t(data$D))
        What_bing <- rbind(What_bing, 
                           matrix(0, ncol=ncol(What_bing), nrow = K - nrow(What_bing)  ))
        
        
      }
    }else{
      What_bing = NULL
      if (dim(t(bing_recovery$A))[1]>K){
        clustered_res <- kmeans(t(bing_recovery$A), centers = K) 
        error <- update_error(t(clustered_res$centers), (What_bing), data$A, t(data$W), 
                              time=elapsed_timeBing, method = "Bing", error=error,thresholded=bing_recovery$thresholded)
      }else{
        if (dim(t(bing_recovery$A))[1]==K){
          error <- update_error((bing_recovery$A), (What_bing), data$A, t(data$W),
				time=elapsed_timeBing, method = "Bing", error=error,thresholded=bing_recovery$thresholded)

  	}
      }
    }
    
    
    # resultsW <- rbind(resultsW, 
    #                   process_results(What_bing, "Bing", seq_len(n), processingA=FALSE))
  }
    


  
  
  #### Step 6: Run method
  for (alpha in c(0.001, 0.003, 0.005, 0.007,  0.01, 0.05,  0.1, 0.5 , 1)){
    elapsed_timeOurs <- system.time({
          score_ours <- tryCatch(
            score(D = t(data$D)/N, K=K, normalize = 'huy', 
                  threshold =TRUE, alpha = alpha, N=N, max_K = min(min(dim(data$D))-1, 150),
                  VHMethod=VHMethod, returnW=evaluateW),
            error = function(err) {
              # Code to handle the error (e.g., print an error message, log the error, etc.)
              paste0("Error occurred while running Score ", alpha, " :", conditionMessage(err), "\n")
              # Return a default value or NULL to continue with the rest of the code
              return(NULL)
            }
          )
    })["elapsed"]
    
    if (is.null(score_ours) == FALSE){
      # resultsA <- rbind(resultsA, 
      #                   process_results(score_ours$A_hat, paste0("Ours_", alpha), data$vocab))
      # resultsW <- rbind(resultsW,
      #                   process_results(score_ours$W_hat, paste0("Ours_", alpha), seq_len(n), processingA=FALSE))
      error <- update_error(score_ours$A_hat, (score_ours$W_hat), data$A, t(data$W), 
                            time = elapsed_timeOurs, method = paste0("Ours_", alpha), error=error,
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


# error <- c()
# for (exp in 1:100){
#   for (K in c(3, 4, 5, 7, 10, 20, 30, 50, 100)){
#     A = NULL
#     W = NULL
#     vocab = NULL
#     for (n in c(c(100, 250, 500, 250, 1000), 2000, 5000)){
#       for (n_frac in c(0.5, 0.8, 1, 2, 5, 10)){
#         N = ceiling(n_frac * n)
#         test <- run_experiment("AP", K, N=N, n=n, seed = (exp * 1000 +  n  + K  + 0.01 * n_frac)*1000, 
#                                A = A, W=W, vocab=vocab)
#         error_temp = test$error
#         error_temp["Khat_huy"]=test$Khat_huy
#         error_temp["Khat_huy_thresh"] = test$Khat_huy_thresh
#         error_temp["Khat_olga"]=test$Khat_olg
#         error_temp["Khat_olga_thresh"] = test$Khat_olga_thresh
#         error_temp["Khat_tracy"]=test$Khat_tracy
#         error_temp["Khat_tracy_thresh"] = test$Khat_tracy_thresh
#         error_temp["N"] = N
#         error_temp["n"] = n
#         error_temp["seed"] =  (exp * 1000 +  n  + K  + 0.01 * n_frac)*1000
#         error_temp["n_frac"] = n_frac
#         error_temp["exp"] = exp
#         error <- rbind(error,
#                        error_temp)
#         write_csv(error, paste0(getwd(), "/r/experiments/semi_synthetic/results_semisynthetic_right_K2.csv"))
#         if (is.null(A)){
#           A = test$Aoriginal
#           W = test$Woriginal
#           vocab=test$vocab
#         }
#       }
#     }
#   }
# }
# 
# 
# ggplot(error %>% filter(n==100)) +
#   geom_line(aes(x=N, y=l1_A, colour=method))


