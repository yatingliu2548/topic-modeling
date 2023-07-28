args = commandArgs(trailingOnly=TRUE)
seed = ceiling(as.numeric(args[2]))
result_file = args[3]
K = ceiling(as.numeric(args[4]))
matlab_path = args[5]
error <- c()


setwd(args[1])
source("r/experiments/semi_synthetic/synthetic_AP.R")




A = NULL
W = NULL
vocab = NULL
for (n in c(c(100, 250, 500, 250, 1000), 2000, 5000)){
    for (n_frac in c(0.5, 0.8, 1, 2, 5, 10)){
    N = ceiling(n_frac * n)
    test <- run_experiment("AP", K, N=N, n=n, seed = seed, 
                            A = A, W=W, vocab=vocab, matlab_path=matlab_path)
    error_temp = test$error
    error_temp["Khat_huy"]=test$Khat_huy
    error_temp["Khat_huy_thresh"] = test$Khat_huy_thresh
    error_temp["Khat_olga"]=test$Khat_olg
    error_temp["Khat_olga_thresh"] = test$Khat_olga_thresh
    error_temp["Khat_tracy"]=test$Khat_tracy
    error_temp["Khat_tracy_thresh"] = test$Khat_tracy_thresh
    error_temp["N"] = N
    error_temp["n"] = n
    error_temp["seed"] =  (exp * 1000 +  n  + K  + 0.01 * n_frac)*1000
    error_temp["n_frac"] = n_frac
    error_temp["exp"] = exp
    error <- rbind(error,
                    error_temp)
    write_csv(error, paste0(getwd(), result_file))
    if (is.null(A)){
        A = test$Aoriginal
        W = test$Woriginal
        vocab=test$vocab
    }
}
    