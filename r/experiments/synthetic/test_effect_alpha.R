setwd("~/Documents/topic-modeling/")
source("r/experiments/semi_synthetic/synthetic_AP.R")
source("r/experiments/synthetic/synthetic_dataset.R")


n=500
K=5
alpha=1
n_max_zipf=50000
a_zipf=1
n_anchors=2
delta_anchor=0.01 
N=500
seed=1234
offset_zipf=2.7
vary_by_topic=FALSE

error = c()
VHMethod  = "SP"
estimateW = FALSE

vary_by_topic = FALSE
normalize_counts = TRUE

for (exp in  1:50){
  for (p in c(5000)){
    data = synthetic_dataset_creation(n, K, p, alpha=1, 
                                      n_max_zipf=50000, 
                                      a_zipf=a_zipf,
                                      n_anchors=n_anchors, 
                                      delta_anchor=delta_anchor, 
                                      N=N, seed=seed + exp,
                                      offset_zipf=offset_zipf,
                                      vary_by_topic=vary_by_topic)
    
    elapsed_timeTracy <- system.time({
      if (normalize_counts){
        score_recovery <- score(t(data$D)/N, K, normalize = "norm", Mquantile = 0,
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
    
    
    error <- update_error(score_recovery$A_hat, NULL, data$A, t(data$W), 
                          time = elapsed_timeTracy, method = "TopicScore", error=error,
                          thresholded = 0)
    for (alpha_thresh in 10^(-seq(from = 0.5, to = 3, length.out = 25))){
      
      
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
      
      error <- update_error(score_ours$A_hat, NULL, data$A, NULL, 
                            time=elapsed_timeOurs, method = paste0("Ours_", alpha_thresh), error=error,
                            thresholded=score_ours$thresholded)
      
      print(paste0("Done with alpha ", alpha_thresh))
      print(error)
      print("Done with p")
      
    }
    write_csv(error, paste0(getwd(), paste0("/r/experiments/synthetic/results/effect_alpha", ".csv")))
  }
}
