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
source("r/experiments/real/run_topic_scores.R")


lda_varying_params_lists <-  list()
for (k in 1:18) {
  lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
}
####

ggplot(data.frame(x= 1:ncol(data), y=sort(apply(data, 2, sum))),
        aes(x=x, y=y))+geom_line()+
  scale_x_log10()+
  scale_y_log10()

res = c()
set.seed(19930108)
for (exp in 1:50){
  for (nb_samples in c(5000, 10000)){
    train_index = sample(1:nrow(data),nb_samples)
    test_index = sample((1:nrow(data))[-c(train_index)], nb_samples)
    
    
    ggplot(data.frame(x= 1:ncol(data), y=sort(apply(data[train_index, ], 2, sum))),
           aes(x=x, y=y))+geom_line()+
      scale_x_log10()+
      scale_y_log10()
    
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
    
    alpha= 0.005


    topic_models_ours <- run_topic_models(as.matrix(data), train_index, list_params=3:18,
                                                 threshold = TRUE,
                                                 normalize="TTS", VHMethod = 'SP',
                                                 Mquantile = 0.00,
                                                 alpha = alpha, max_K=150)
      
    topic_models_ours_test <- run_topic_models(as.matrix(data), train_index, list_params=1:18,
                                           threshold = TRUE,
                                           normalize="TTS", VHMethod = 'SP',
                                           Mquantile = 0.00,
                                           alpha = alpha, max_K=150)
      
    
      


    
    
    topic_models <- run_topic_models(as.matrix(data), train_index, list_params=3:18,
                                     threshold = FALSE,
                                     normalize="norm", VHMethod = 'SP',
                                     Mquantile = 0.0,
                                     alpha = 0.005, max_K=150)
    
    topic_models_test <- run_topic_models(as.matrix(data), test_index, list_params=3:18, 
                                          threshold = FALSE,
                                          normalize="norm", VHMethod = 'SP',
                                          Mquantile = 0.0,
                                          alpha = 0.005, max_K=150)
    
    

    for (k in 3:18){
      it = 1

        match = data.frame((exp(topic_models_ours[[paste0("k",k)]]$beta))%*% t((exp(topic_models_ours_test[[paste0("k",k)]]$beta))))
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
                              "method" = paste0("ours", alpha),
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
                             "m" = nb_samples,
                              "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
        res = rbind(res, res_temp)
        
        it = it + 1
      
      
      
      match = data.frame((exp(topic_models[[paste0("k",k)]]$beta))%*% t((exp(topic_models_test[[paste0("k",k)]]$beta))))
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
                           "m" = nb_samples,
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
                           "m" = nb_samples,
                            "off_diag" = (sum(match_permuted) - sum(diag(as.matrix(match_permuted))))/(ncol(match_permuted)^2 - ncol(match_permuted)) )
      res = rbind(res, res_temp)
      write_csv(res, "~/Documents/topic-modeling/r/experiments/real/spleen_results_final.csv")
      
    }
    
    
  }
  
}

res = read_csv( "~/Documents/topic-modeling/r/experiments/real/spleen_results5.csv")
res[,1:ncol(res)] %>%
  group_by(method, k,m) %>%
  summarise(c=n())

ggplot(res[,1:ncol(res)] %>%
         #filter(exp<13) %>%
         group_by(method, k,m) %>%
         summarise_all(median) %>%
         select(-c("min", "max", "mean", "median",
                   "q25", "q75")) %>%
         pivot_longer(cols=-c("method", "k", "m")) %>%
         mutate(name= as.numeric(str_extract(name, "(?<=res)[0-9.]+"))) %>%
         drop_na(name) %>%
         mutate(name=ifelse(name> 10, 0.1 * name, name)) %>%
         filter(m==5000, 
       aes(x=0.1 * name, y=value/k, colour=method, fill=method))+
  #geom_smooth(alpha=0.2, se = FALSE)+
  geom_line()+
  facet_wrap(.~ k, scales = "free_x") + 
  ylab("Cumulative distribution") +
  xlab("Topic Resolution") 

theme_set(theme_bw(base_size = 14))
my_colors <- c(
  "dodgerblue", "red","orchid3")

ggplot(res[,1:ncol(res)] %>%
         #filter(exp<13) %>%
         group_by(method, k,m) %>%
         summarise_all(mean) %>%
         filter(m==5000, method !="ours0.005"), 
       aes(x=k, y=median, colour=method, fill=method))+
  #geom_smooth(alpha=0.2, se = FALSE)+
  geom_line(position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5), size=3)+
  geom_errorbar(aes(ymin=q25, ymax=q75),
                position=position_dodge(width=0.5), alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = c( "ours0.005", "score", "lda"),
                   labels=c('TTM (this paper)', 
                            'Topic Score (Ke et al)', 
                            'LDA \n(Blei et al)')) +
  labs(colour="Method") + 
  #facet_wrap(.~ k, scales = "free_x") + 
  ylab("Average Topic Resolution") +
  xlab("Number of Topics") 



topic_to_compare <- vector("list", 2)
topic_to_compare[[1]] = topic_models$k15
topic_to_compare[[2]] = topic_models_test$k15
names(topic_to_compare) <- c("Train", "Test")
aligned_topics_transport_comp2 <- 
  align_topics(
    models = topic_to_compare,
    method = "transport", reg=0.001) 
plot(aligned_topics_transport_comp2, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp_lda <- 
  align_topics(
    models = lda_models,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp_lda, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_ours,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp2 <- 
  align_topics(
    models = topic_models,
    method = "transport", reg=0.01) 

plot(aligned_topics_transport_comp2, add_leaves = TRUE, label_topics = TRUE)

plot(aligned_topics_transport_comp, color_by = "coherence")


match = data.frame(((topic_models_ours$k10$beta))%*% t(((topic_models_ours_test$k10$beta))))
colnames(match)= 1:ncol(match)
match["X"] = 1:nrow(match)
ggplot(pivot_longer(match, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))


match = data.frame(((topic_models_oursC$k15$beta))%*% t(((topic_models_oursC_test$k15$beta))))
colnames(match)= 1:ncol(match)
match["X"] = 1:nrow(match)
ggplot(pivot_longer(match, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))


match = data.frame(((topic_models$k15$beta))%*% t(((topic_models_test$k15$beta))))
colnames(match)= 1:ncol(match)
match["X"] = 1:nrow(match)
ggplot(pivot_longer(match, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))


match = data.frame((exp(lda_models$k15$beta))%*% t((exp(lda_models_test$k15$beta))))
colnames(match)= 1:ncol(match)
match["X"] = 1:nrow(match)
ggplot(pivot_longer(match, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))

aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_ours,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp_score <- 
  align_topics(
    models = topic_models,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp_score, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_lda <- 
  align_topics(
    models = lda_models,
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp_lda, add_leaves = TRUE, label_topics = TRUE)


