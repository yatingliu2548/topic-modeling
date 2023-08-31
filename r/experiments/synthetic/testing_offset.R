source("r/experiments/semi_synthetic/synthetic_AP.R")
source("r/experiments/synthetic/synthetic_dataset.R")

args = commandArgs(trailingOnly=TRUE)
seed <-  ceiling(as.numeric(args[1]))
result_file <- args[2]
K <- ceiling(as.numeric(args[3]))
#matlab_path <-  args[4]
alpha_dirichlet <-  as.numeric(args[4])
n_anchors <-  as.numeric(args[5])
delta_anchor <- as.numeric(args[6])

#matlab_path = DEFAULT_MATLAB
error <- c()

n <- 1000
a_zipf= 2.7
print(paste0("Delta anchor is: ", delta_anchor))

for (vary_by_topic in c(TRUE, FALSE)){
for (b_zipf in c(c(2.7, 10, 50, 100, 1000), 10000)){
    for (N in sort(c(500), decreasing=TRUE)){
        for (p in c(1000, 5000, 10000, 15000)){
            if (K <5){
                VHMethod = "SVS"
            }else{
                VHMethod = "SP"
            }
            test <- run_synthetic_experiment(n, K, p, alpha_dirichlet=alpha_dirichlet, 
                                                a_zipf=a_zipf,   offset_zipf = b_zipf,
                                                n_anchors=n_anchors, delta_anchor=delta_anchor, N=N,
                                                seed=seed, VHMethod=VHMethod,
                                                data_generation_method=1, vary_by_topic = vary_by_topic,
                                                normalize_counts = TRUE,
                                                sparsity = TRUE)
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
            error_temp["delta_anchor"] =  delta_anchor
            error_temp["a_zipf"] = 1
            error_temp["offset_zipf"] = b_zipf
            error_temp["alpha_dirichlet"] = alpha_dirichlet
            error_temp["exp"] = result_file
            error_temp["VHMethod"] = VHMethod
            error_temp["delta_anchor"] = delta_anchor
            error_temp["vary_by_topic"] = vary_by_topic
	    error <- rbind(error,
                            error_temp)
            write_csv(error, paste0(getwd(), paste0("/r/experiments/synthetic/results/",paste0(result_file, '_delta_anchor', delta_anchor, '_K_' , K), ".csv")))

        }
    }
        
}
        
}



