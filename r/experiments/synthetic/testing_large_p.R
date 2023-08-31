
setwd("~/Documents/topic-modeling/")


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
source("r/experiments/synthetic/synthetic_dataset.R")

result_file = "experiment_large_p3_"
error <- c()


delta_anchor = 1e-3
n_anchors = 5
a_zipf = 1
normalize_counts = TRUE
#vary_by_topic = FALSE
zipf_offset = 2.7
estimateW = F
estimateK = F
alpha_thresh = 0.005
alpha_dirichlet = 1
N = 500
K = 5

for (seed in 1:50){
  for (n in c(c(1000, 5000, 10000))){
    for (p in c(20000, 50000)){
      for (vary_by_topic in c(TRUE, FALSE)){
        if (K <5){
          VHMethod = "SVS"
        }else{
          VHMethod = "SP"
        }
        
        
        
        data = synthetic_dataset_creation(n, K, p, alpha=alpha_dirichlet, 
                                          n_max_zipf=50000, 
                                          a_zipf=a_zipf,
                                          n_anchors=n_anchors, 
                                          delta_anchor=delta_anchor, 
                                          N=N, seed=seed,
                                          offset_zipf=zipf_offset,
                                          vary_by_topic=vary_by_topic)
        
        
        
        if (p<20001){
          elapsed_timeTracy <- system.time({
            if (normalize_counts){
              score_recovery <- score(t(data$D)/N, K, normalize = "norm", Mquantile = 0.00,
                                      max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod,
                                      returnW=estimateW, estimateK=FALSE)
              # score_recovery <- evalWithTimeout({
              #   # Some potentially long-running code here
              #   score(t(data$D)/N, K, normalize = "norm", 
              #         max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod,
              #         returnW=estimateW)
              #   print("Tracy's code finished!")
              #   Khat_tracy = select_K(score_recovery$eigenvalues, p, n, N, method="tracy")
              #   Khat_olga = select_K(svd(data$D)$d, p,n, N, method="olga")
              #   
              # }, timeout = 300, onTimeout = "warning") ### stops if longer than 5min
              
            }else{
              # score_recovery <- evalWithTimeout({
              #   # Some potentially long-running code here
              #   score(t(data$D), K, normalize = "norm", 
              #         max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod,
              #         returnW=estimateW)
              #   print("Tracy's code finished!")
              #   Khat_tracy = select_K(score_recovery$eigenvalues, p, n, N, method="tracy")
              #   Khat_olga = select_K(svd(data$D)$d, p,n, N, method="olga")
              # }, timeout = 300, onTimeout = "warning") ### stops if longer than 5min
              score_recovery <- score(t(data$D), K, normalize = "norm",
                                      max_K = min(150, min(dim(data$D)-1)), VHMethod=VHMethod,
                                      returnW=estimateW, estimateK=FALSE)
            }
          })["elapsed"]
          
          
          error_temp <- update_error(score_recovery$A_hat, NULL, data$A, t(data$W), 
                                time = elapsed_timeTracy, method = "TopicScore", error=NULL,
                                thresholded = 0)
          
        }
        
       if (p > 20001){
          error_temp = NULL
       } 
        

        elapsed_timeOurs <- system.time({
          score_ours <- tryCatch(
            score(D = t(data$D)/N, K=K, 
                  normalize = 'huy', 
                  threshold =TRUE, alpha = alpha_thresh, N=N, max_K = min(min(dim(data$D))-1, 150),
                  VHMethod=VHMethod, returnW = estimateW, estimateK = FALSE)
            ,
            error = function(err) {
              # Code to handle the error (e.g., print an error message, log the error, etc.)
              paste0("Error occurred while running Score ", alpha_thresh, " :", conditionMessage(err), "\n")
              # Return a default value or NULL to continue with the rest of the code
              return(NULL)
            }
          )
        })["elapsed"]
        
        error_temp <- update_error(score_ours$A_hat, NULL, data$A, NULL, 
                              time=elapsed_timeOurs, method = paste0("Ours_", alpha_thresh), error=error_temp,
                              thresholded=score_ours$thresholded)

        error_temp["N"] = N
        error_temp["n"] = n
        error_temp["p"] = p
        error_temp["seed"] =  seed
        error_temp["n_anchors"] =  n_anchors
        error_temp["alpha"] = a_zipf
        error_temp["alpha_dirichlet"] = alpha_dirichlet
        error_temp["zipf_offset"] = zipf_offset
        error_temp["delta_anchor"] = delta_anchor
        error_temp["exp"] = result_file
        error_temp["vary_by_topic"] = vary_by_topic
        error_temp["VHMethod"] = VHMethod
        error <- rbind(error,
                       error_temp)
        write_csv(error, paste0(getwd(), paste0("/r/experiments/synthetic/results/",paste0(result_file, '_delta_anchor', delta_anchor, '_K_' , K), ".csv")))
      }
      
    }
    
  }
  
}



