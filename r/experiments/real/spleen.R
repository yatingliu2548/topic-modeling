library(tidyverse)
library(alto)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tibble)
library(stringr)
library(tidyr)
library(clue)
library(combinat)
library(combinat)
library(clue)

setwd("~/Documents/topic-modeling")
source("r/score.R")

data <- read_csv("~/Downloads/spleen_cells_features.csv")
setwd("~/Documents/topic-modeling")
source("r/score.R")
source("r/experiments/real/run_topic_scores.R")


lda_varying_params_lists <-  list()
for (k in 1:15) {
  lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
}
####

res = c()
for (exp in 1:50){
  for (m in c(1000, 5000)){
    train_index = sample(1:nrow(data), m)
    test_index = sample((1:nrow(data))[-c(train_index)], m)
    
    
    
    lda_models  <- 
      run_lda_models(
        data = data[train_index,],
        lda_varying_params_lists = lda_varying_params_lists,
        lda_fixed_params_list = list(method = "VEM"),
        dir = "spleen_lda_models_train/",
        reset = TRUE,
        verbose = TRUE
      )
    
    lda_models_test  <- 
      run_lda_models(
        data = data[test_index,],
        lda_varying_params_lists = lda_varying_params_lists,
        lda_fixed_params_list = list(method = "VEM"),
        dir = "spleen_models_test/",
        reset = TRUE,
        verbose = TRUE
      )
    
    
    topic_models <- run_topic_models(as.matrix(data), train_index, list_params=3:15,
                                     threshold = FALSE,
                                     normalize="none", VHMethod = 'SP',
                                     Mquantile = 0.00,
                                     alpha = 0.005, max_K=150)
    
    topic_models_test <- run_topic_models(as.matrix(data), test_index, list_params=3:15, 
                                          threshold = FALSE,
                                          normalize="none", VHMethod = 'SP',
                                          Mquantile = 0.00,
                                          alpha = 0.005, max_K=150)
    
    
    alphas= c(0.0001,0.005)
    topic_models_huy <- vector("list", length(alphas))
    topic_models_huy_test <- vector("list", length(alphas))
    topic_models_huyNC <- vector("list", length(alphas))
    topic_models_huyNC_test <- vector("list", length(alphas))
    it = 1
    for(alpha in alphas){
      topic_models_huy[[it]] <- run_topic_models(as.matrix(data), train_index, list_params=3:15,
                                                 threshold = TRUE,
                                                 normalize="huy", VHMethod = 'SP',
                                                 Mquantile = 0.00,
                                                 alpha = alpha, max_K=150)
      
      
      
      topic_models_huy_test[[it]] <- run_topic_models(as.matrix(data), test_index, list_params=3:15, 
                                                      threshold = TRUE,
                                                      normalize="huy", VHMethod = 'SP',
                                                      Mquantile = 0.00,
                                                      alpha = alpha, max_K=150)
      
      
      topic_models_huyNC[[it]] <- run_topic_models(as.matrix(data), train_index, list_params=3:15,
                                                 threshold = TRUE,
                                                 normalize="huy_not_centered", VHMethod = 'SP',
                                                 Mquantile = 0.00,
                                                 alpha = alpha, max_K=150)
      
      
      
      topic_models_huyNC_test[[it]] <- run_topic_models(as.matrix(data), test_index, list_params=3:15, 
                                                      threshold = TRUE,
                                                      normalize="huy_not_centered", VHMethod = 'SP',
                                                      Mquantile = 0.00,
                                                      alpha = alpha, max_K=150)
      
      it = it + 1
    }
    
    
    thres = 0.1
    for (k in 3:15){
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
                              "m" = m,
                              "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
        res = rbind(res, res_temp)
        
        match = data.frame(((topic_models_huyNC[[it]][[paste0("k",k)]]$beta))%*% t(((topic_models_huyNC_test[[it]][[paste0("k",k)]]$beta))))
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
                              "method" = paste0("huyNC", alpha),
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
                              "m" = m,
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
                            "m" = m,
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
                            "m" = m,
                            "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
      res = rbind(res, res_temp)
      write_csv(res, "~/Documents/topic-modeling/r/experiments/real/spleen_results4.csv")
      
    }
    
    
  }
  
}


ggplot(res[,1:ncol(res)] %>%
         group_by(method, k) %>%
         summarise_all(median) %>%
         select(-c("min", "max", "mean", "median",
                   "q25", "q75")) %>%
         pivot_longer(cols=-c("method", "k")) %>%
         mutate(name= as.numeric(str_extract(name, "(?<=res)[0-9.]+"))) %>%
         mutate(name=ifelse(name> 10, 0.1 * name, name)), 
       aes(x=0.1 * name, y=value/k, colour=method, fill=method))+
  #geom_smooth(alpha=0.2,se = FALSE)+
  geom_line()+
  facet_wrap(.~ k, scales = "free_x") + 
  ylab("Cumulative distribution") +
  xlab("Topic Resolution") 



topic_to_compare <- vector("list", 2)
topic_to_compare[[1]] = lda_models$k8
topic_to_compare[[2]] = lda_models_test$k8
names(topic_to_compare) <- c("Train", "Test")


aligned_topics_transport_comp_lda <- 
  align_topics(
    models = lda_models,
    method = "transport", reg=0.001) 
plot(aligned_topics_transport_comp_lda, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_huy[[1]],
    method = "transport", reg=0.001) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp2 <- 
  align_topics(
    models = topic_models,
    method = "transport", reg=0.01) 

plot(aligned_topics_transport_comp2, add_leaves = TRUE, label_topics = TRUE)

plot(aligned_topics_transport_comp, color_by = "coherence")





