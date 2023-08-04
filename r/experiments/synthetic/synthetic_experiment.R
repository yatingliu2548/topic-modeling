source("r/experiments/semi_synthetic/synthetic_AP.R")
source("r/experiments/synthetic/synthetic_dataset.R")

args = commandArgs(trailingOnly=TRUE)
seed = ceiling(as.numeric(args[1]))
result_file = args[2]
K = ceiling(as.numeric(args[3]))
matlab_path = args[4]
error <- c()



anchors = c(0, 1, 5, 10)
tot = sapply(anchors, function(x){x * K})
p = 10000

for (n in c(c(100, 250, 500, 250, 1000), 2000)){
  for (n_frac in c(0.5, 0.8, 1, 2, 5, 10)){
    for (noise_level in c(0, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5)){
      for (noise_generation in c("uniform", "zipf")){
        for (n_anchors in anchors){
          N = ceiling(n_frac * n)
          if (K <5){
            VHMethod = "SVS"
          }else{
            VHMethod = "SP"
          }
          test <- run_synthetic_experiment(n, K, p, alpha=0.5, a_zipf=1,
                                           n_anchors=0, delta_anchor=1, N=N,
                                           noise_generation = noise_generation, seed=seed, VHMethod=VHMethod)
          error_temp = test$error
          error_temp["Khat_huy"]=test$Khat_huy
          error_temp["Khat_huy_thresh"] = test$Khat_huy_thresh
          error_temp["Khat_olga"]=test$Khat_olga
          error_temp["Khat_olga_thresh"] = test$Khat_olga_thresh
          error_temp["Khat_tracy"]=test$Khat_tracy
          error_temp["Khat_tracy_thresh"] = test$Khat_tracy_thresh
          error_temp["N"] = N
          error_temp["n"] = n
          error_temp["seed"] =  seed
          error_temp["n_anchors"] =  n_anchors
          error_temp["noise_generation"] =  noise_generation
          error_temp["n_frac"] = n_frac
          error_temp["exp"] = result_file
          error_temp["VHMethod"] = VHMethod
          error_temp["noise_level"] = ifelse((noise_level)=="auto", test$Epsilon, as.numeric(noise_level))
          error <- rbind(error,
                         error_temp)
          write_csv(error, paste0(getwd(), paste0("/r/experiments/synthetic/results/",result_file)))
        }
      }
    }
  }
}


