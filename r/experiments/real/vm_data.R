library(alto)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tibble)
library(stringr)
library(tidyr)

load("~/Documents/alto/data/vm_data.rda")
names(vm_data)
library(combinat)
library(clue)

setwd("~/Documents/topic-modeling")
source("r/score.R")

# we create the human-friendly ASV names
# by concatenating Genus and Species and an identifying number
# as a given species can be represented by several ASVs
vm_data$taxonomy <-
  vm_data$taxonomy %>%
  as.data.frame() %>%
  group_by(Genus, Species) %>%
  mutate(
    ASV_name =
      paste0(
        ifelse(is.na(Genus),"[unknown Genus]", Genus), " (",
        ifelse(is.na(Species),"?", Species), ") ",
        row_number()
      )
  ) %>%
  ungroup()

# then, we replace the colnames of the count matrix by these new names.
# Note that the tax table rows are ordered as the count matrix columns,
# allowing us to assign without matching first.
colnames(vm_data$counts) <-  vm_data$taxonomy$ASV_name



run_topic_models <- function(data, train_index, #test_index, 
                             list_params=1:20, threshold = FALSE,
                             normalize="none", VHMethod = 'SP',
                             Mquantile = 0.05,
                             alpha = 0.005, max_K=150){
  ####
  n_total = 1:ncol(data$counts)
  active_train = which(apply(data$counts[train_index,], 2, sum)>0)
  X = as.matrix(data$counts)
  x_train = t(diag(1/ apply(X[train_index,active_train],1, sum)) %*% X[train_index, active_train])
  #active_test = which(apply(data$counts[test_index,], 2, sum)>0)
  #x_test = t(diag(1/ apply(data$counts[test_index, active_test],1, sum)) %*% data$counts[test_index, active_test])
  
  topic_models <- vector("list", length(list_params))  
  #topic_models_test <- vector("list", length(list_params)) 
  it = 1
  for (k in list_params){
    
    tm <- score(as.matrix(x_train), k, scatterplot=FALSE, 
                K0=NULL, m=NULL, N=mean(apply(data$counts, 1, sum)), threshold=threshold,
                Mquantile=Mquantile, VHMethod = VHMethod, normalize=normalize,
                alpha=alpha, max_K=max_K, returnW=TRUE, estimateK=FALSE)
    A_hat = matrix(0, ncol(data$counts), k)
    A_hat[active_train, ] = tm$A_hat
    #A_hat = tm$A_hat
    topic_models[[it]] <- list(
      beta = t(A_hat) %>% magrittr::set_rownames(rownames(data)),
      gamma =  t(tm$W_hat) %>% magrittr::set_colnames(colnames(data))
    )
    it <- it + 1
  }
  names(topic_models) <- sapply(list_params,function(x){paste0("k",x)})
  #names(topic_models_test) <- sapply(list_params,function(x){paste0("k",x)})
  return(topic_models)
  
}

lda_varying_params_lists <-  list()


for (k in 1:20) {
  lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
}



res = c()
for (exp in 1:50){
  print("exp is")
  print(exp)
  train_index = sample(1:nrow(vm_data$counts), ceiling(nrow(vm_data$counts)/2))
  test_index = (1:nrow(vm_data$counts))[-c(train_index)]
  
  for (k in 1:20) {
    lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
  }
  
  lda_models  <- 
    run_lda_models(
      data = vm_data$counts[train_index,],
      lda_varying_params_lists = lda_varying_params_lists,
      lda_fixed_params_list = list(method = "VEM"),
      dir = "microbiome_lda_models_train/",
      reset = TRUE,
      verbose = TRUE
    )
  
  lda_models_test  <- 
    run_lda_models(
      data = vm_data$counts[test_index,],
      lda_varying_params_lists = lda_varying_params_lists,
      lda_fixed_params_list = list(method = "VEM"),
      dir = "microbiome_lda_models_test/",
      reset = TRUE,
      verbose = TRUE
    )
  
  topic_models <- run_topic_models(vm_data, train_index, list_params=3:20,
                                   threshold = FALSE,
                                   normalize="none", VHMethod = 'SP',
                                   Mquantile = 0.00,
                                   alpha = 0.005, max_K=150)
  
  topic_models_test <- run_topic_models(vm_data, test_index, list_params=3:20, threshold = FALSE,
                                        normalize="none", VHMethod = 'SP',
                                        Mquantile = 0.00,
                                        alpha = 0.005, max_K=150)
  
  
  
  
  alphas= c(0.0001, 0.001, 0.005, 0.01)
  topic_models_huy <- vector("list", length(alphas))
  topic_models_huy_test <- vector("list", length(alphas))
  it = 1
  for(alpha in alphas){
    topic_models_huy[[it]] <- run_topic_models(vm_data, train_index, list_params=3:20, threshold = TRUE,
                                               normalize="huy", VHMethod = 'SP',
                                               Mquantile = 0.00,
                                               alpha = alpha, max_K=150)
    
    
    
    topic_models_huy_test[[it]] <- run_topic_models(vm_data, test_index, list_params=3:20, threshold = TRUE,
                                                    normalize="huy", VHMethod = 'SP',
                                                    Mquantile = 0.00,
                                                    alpha = alpha, max_K=150)
    it = it + 1
  }
  
  thres = 0.1
  
  for (k in 3:20){
    it = 1
    for (alpha in alphas){
      match = data.frame(((topic_models_huy[[it]][[paste0("k",k)]]$beta))%*% t(((topic_models_huy_test[[it]][[paste0("k",k)]]$beta))))
      #match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
      permutation <- solve_LSAP(as.matrix(match), maximum=TRUE) 
      match_permuted <- match[, permutation]
      res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                            "max" = max(diag(as.matrix(match_permuted))),
                            "mean" = mean(diag(as.matrix(match_permuted))), 
                            "median" = median(diag(as.matrix(match_permuted))), 
                            "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                            "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                            "k" = k,
                            "method" = paste0("huy", alpha),
                            "res1" = sum(diag(as.matrix(match_permuted))>0.1),
                            "res15" = sum(diag(as.matrix(match_permuted))>0.15),
                            "res2" = sum(diag(as.matrix(match_permuted))>0.2),
                            "res2.5" = sum(diag(as.matrix(match_permuted))>0.25),
                            "res3" = sum(diag(as.matrix(match_permuted))>0.3),
                            "res3.5" = sum(diag(as.matrix(match_permuted))>0.35),
                            "res45" = sum(diag(as.matrix(match_permuted))>0.4),
                            "res4" = sum(diag(as.matrix(match_permuted))>0.45),
                            "res5" = sum(diag(as.matrix(match_permuted))>0.5),
                            "res6" = sum(diag(as.matrix(match_permuted))>0.6),
                            "res7" = sum(diag(as.matrix(match_permuted))>0.7),
                            "res8" = sum(diag(as.matrix(match_permuted))>0.8),
                            "res9" = sum(diag(as.matrix(match_permuted))>0.9),
                            "exp" = exp,
                            "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
      res = rbind(res, res_temp)
      it = it + 1
    }
    
    
    match = data.frame(((topic_models[[paste0("k",k)]]$beta))%*% t(((topic_models_test[[paste0("k",k)]]$beta))))
    #match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
    permutation <- solve_LSAP(as.matrix(match), maximum=TRUE) 
    match_permuted <- match[, permutation]
    res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                          "max" = max(diag(as.matrix(match_permuted))),
                          "mean" = mean(diag(as.matrix(match_permuted))), 
                          "median" = median(diag(as.matrix(match_permuted))), 
                          "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                          "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                          "k" = k,
                          "method" = "score",
                          "res1" = sum(diag(as.matrix(match_permuted))>0.1),
                          "res15" = sum(diag(as.matrix(match_permuted))>0.15),
                          "res2" = sum(diag(as.matrix(match_permuted))>0.2),
                          "res2.5" = sum(diag(as.matrix(match_permuted))>0.25),
                          "res3" = sum(diag(as.matrix(match_permuted))>0.3),
                          "res3.5" = sum(diag(as.matrix(match_permuted))>0.35),
                          "res45" = sum(diag(as.matrix(match_permuted))>0.4),
                          "res4" = sum(diag(as.matrix(match_permuted))>0.45),
                          "res5" = sum(diag(as.matrix(match_permuted))>0.5),
                          "res6" = sum(diag(as.matrix(match_permuted))>0.6),
                          "res7" = sum(diag(as.matrix(match_permuted))>0.7),
                          "res8" = sum(diag(as.matrix(match_permuted))>0.8),
                          "res9" = sum(diag(as.matrix(match_permuted))>0.9),
                          "exp" = exp,
                          "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
    res = rbind(res, res_temp)
    
    match = data.frame((exp(lda_models[[paste0("k",k)]]$beta))%*% t((exp(lda_models_test[[paste0("k",k)]]$beta))))
    #match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
    permutation <- solve_LSAP(as.matrix(match), maximum=TRUE) 
    match_permuted <- match[, permutation]
    res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                          "max" = max(diag(as.matrix(match_permuted))),
                          "mean" = mean(diag(as.matrix(match_permuted))), 
                          "median" = median(diag(as.matrix(match_permuted))), 
                          "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                          "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                          "k" = k,
                          "method" = "lda",
                          "res1" = sum(diag(as.matrix(match_permuted))>0.1),
                          "res15" = sum(diag(as.matrix(match_permuted))>0.15),
                          "res2" = sum(diag(as.matrix(match_permuted))>0.2),
                          "res2.5" = sum(diag(as.matrix(match_permuted))>0.25),
                          "res3" = sum(diag(as.matrix(match_permuted))>0.3),
                          "res3.5" = sum(diag(as.matrix(match_permuted))>0.35),
                          "res45" = sum(diag(as.matrix(match_permuted))>0.4),
                          "res4" = sum(diag(as.matrix(match_permuted))>0.45),
                          "res5" = sum(diag(as.matrix(match_permuted))>0.5),
                          "res6" = sum(diag(as.matrix(match_permuted))>0.6),
                          "res7" = sum(diag(as.matrix(match_permuted))>0.7),
                          "res8" = sum(diag(as.matrix(match_permuted))>0.8),
                          "res9" = sum(diag(as.matrix(match_permuted))>0.9),
                          "exp" = exp,
                          "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
    res = rbind(res, res_temp)
    
  }
  write.csv(res,"~/Documents/topic-modeling/r/experiments/real/results_vm_data_with_Mquantile.csv")
  
  
}


