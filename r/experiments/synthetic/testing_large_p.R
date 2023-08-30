#setwd("~/Documents/topic-modeling/")
source("r/experiments/semi_synthetic/synthetic_AP.R")
source("r/experiments/synthetic/synthetic_dataset.R")

result_file = "experiment_large_p"
error <- c()


delta_anchor = 1e-3
n_anchors = 5
a_zipf = 1
#vary_by_topic = FALSE
zipf_offset = 2.7
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
        test <- run_synthetic_experiment(n, K, p, alpha=alpha_dirichlet, a_zipf=a_zipf,
                                         n_anchors=n_anchors, delta_anchor=delta_anchor, N=N,
                                         seed=  seed, 
                                         VHMethod=VHMethod,
                                         data_generation_method=1,
                                         normalize_counts = TRUE, 
                                         offset_zipf=zipf_offset,
                                         vary_by_topic=vary_by_topic)
        error_temp = test$error
        error_temp["Khat_huy"]=test$Khat_huy
        error_temp["Khat_huy_thresh"] = test$Khat_huy_thresh
        error_temp["Khat_olga"]= test$Khat_olga
        error_temp["Khat_olga_thresh"] = test$Khat_olga_thresh
        error_temp["Khat_tracy"]=test$Khat_tracy
        error_temp["Khat_tracy_thresh"] = test$Khat_tracy_thresh
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



