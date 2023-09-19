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
source("r/align_topics.R")
source("r/experiments/real/run_topic_scores.R")

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


lda_varying_params_lists <-  list()


for (k in 1:15) {
  lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
}

colnames(vm_data$sample_info)

vm_data$sample_info %>%
  group_by(SubjectID) %>%
  summarise(count = n())

apply(vm_data$counts, 1, sum)

##### Select just a few samples
dim(vm_data$counts)
apply(vm_data$counts, 1, sum)

apply(vm_data$counts, 1, sum)


res = c()
m = "all"
selected_samples = 1:nrow(vm_data$counts)
X = asinh(X)
for (exp in 1:50){
   
  #for (m in c(5, 10, 20)){
    print("exp is")
    print(exp)
    #selected_subjects = sample(unique(vm_data$sample_info$SubjectID), m)
    #selected_samples = which(vm_data$sample_info$SubjectID %in% selected_subjects)
    #X = vm_data$counts[selected_samples, ]
    #train_index = sample(1:nrow(X), ceiling(length(selected_samples)/2))
    X = vm_data$counts
    train_index = sample(1:nrow(X), ceiling(nrow(X)/2))
    test_index = (1:nrow(X))[-c(train_index)]
    
    
    
    lda_models  <-
      run_lda_models(
        data = ceiling(X)[train_index,],
        lda_varying_params_lists = lda_varying_params_lists,
        lda_fixed_params_list = list(method = "VEM"),
        dir = "microbiome_lda_models_train/",
        reset = TRUE,
        verbose = TRUE
      )

    lda_models_test  <-
      run_lda_models(
        data = ceiling(X)[test_index,],
        lda_varying_params_lists = lda_varying_params_lists,
        lda_fixed_params_list = list(method = "VEM"),
        dir = "microbiome_lda_models_test/",
        reset = TRUE,
        verbose = TRUE
      )
    
    topic_models <- run_topic_models(X, train_index, list_params=3:15,
                                     threshold = TRUE,
                                     normalize="TTS", VHMethod = 'SP',
                                     Mquantile = 0.00,
                                     alpha = 0.005, max_K=150)
    
    topic_models_test <- run_topic_models(X, test_index, list_params=3:15, threshold = FALSE,
                                          normalize="norm", VHMethod = 'SP',
                                          Mquantile = 0.05,
                                          alpha = 0.005, max_K=150)
    
    
    
    
    alphas= c(0.005)
    topic_models_TTS <- vector("list", length(alphas))
    topic_models_TTS_test <- vector("list", length(alphas))
    topic_models_TTSNC <- vector("list", length(alphas))
    topic_models_TTSNC_test <- vector("list", length(alphas))
    it = 1
    for(alpha in alphas){
      topic_models_TTS[[it]] <- run_topic_models(X, train_index, list_params=2:15, threshold = TRUE,
                                                 normalize="TTS", VHMethod = 'SP',
                                                 Mquantile = 0.00,
                                                 alpha = alpha, max_K=150)
      
      
      
      topic_models_TTS_test[[it]] <- run_topic_models(X, test_index, list_params=2:15, threshold = TRUE,
                                                      normalize="TTS", VHMethod = 'SP',
                                                      Mquantile = 0.00,
                                                      alpha = alpha, max_K=150)
      
      topic_models_TTSNC[[it]] <- run_topic_models(X, train_index, list_params=2:15, threshold = TRUE,
                                                 normalize="TTS_not_centered", VHMethod = 'SP',
                                                 Mquantile = 0.00,
                                                 alpha = alpha, max_K=150)
      
      
      
      topic_models_TTSNC_test[[it]] <- run_topic_models(X, test_index, list_params=2:15, threshold = TRUE,
                                                      normalize="TTS_not_centered", VHMethod = 'SP',
                                                      Mquantile = 0.00,
                                                      alpha = alpha, max_K=150)
      it = it + 1
    }
    
    thres = 0.1
    for (k in 3:15){
      it = 1
      for (alpha in alphas){
        alignment <- align_topics(exp(topic_models_TTS[[it]][[paste0("k",k)]]$beta),
                                  exp(topic_models_TTS_test[[it]][[paste0("k",k)]]$beta), 
                                  dist="l1", do.plot=FALSE)
        match_permuted <- alignment$match
        res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                              "max" = max(diag(as.matrix(match_permuted))),
                              "mean" = mean(diag(as.matrix(match_permuted))), 
                              "median" = median(diag(as.matrix(match_permuted))), 
                              "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                              "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                              "k" = k,
                              "m" = m,
                              "n_tot" = length(selected_samples),
                              "method" = paste0("TTS", alpha),
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
        
        
        alignment <- align_topics(exp(topic_models_TTSNC[[it]][[paste0("k",k)]]$beta),
                                  exp(topic_models_TTSNC_test[[it]][[paste0("k",k)]]$beta), 
                                  dist="l1", do.plot=FALSE)
        match_permuted <- alignment$match
        res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                              "max" = max(diag(as.matrix(match_permuted))),
                              "mean" = mean(diag(as.matrix(match_permuted))), 
                              "median" = median(diag(as.matrix(match_permuted))), 
                              "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                              "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                              "k" = k,
                              "m" = m,
                              "n_tot" = length(selected_samples),
                              "method" = paste0("TTSNC", alpha),
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
      
      
      alignment <- align_topics(exp(topic_models[[paste0("k",k)]]$beta),
                                exp(topic_models_test[[paste0("k",k)]]$beta), 
                                dist="l1", do.plot=FALSE)
      match_permuted <- alignment$match
      res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                            "max" = max(diag(as.matrix(match_permuted))),
                            "mean" = mean(diag(as.matrix(match_permuted))), 
                            "median" = median(diag(as.matrix(match_permuted))), 
                            "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                            "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                            "k" = k,
                            "m" = m,
                            "n_tot" = length(selected_samples),
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
      
      alignment <- align_topics(exp(lda_models[[paste0("k",k)]]$beta),
                                exp(lda_models_test[[paste0("k",k)]]$beta), 
                                dist="l1", do.plot=FALSE)
      match_permuted <- alignment$match
      res_temp = data.frame("min"=min(diag(as.matrix(match_permuted))), 
                            "max" = max(diag(as.matrix(match_permuted))),
                            "mean" = mean(diag(as.matrix(match_permuted))), 
                            "median" = median(diag(as.matrix(match_permuted))), 
                            "q25" = quantile(diag(as.matrix(match_permuted)), 0.25),
                            "q75" = quantile(diag(as.matrix(match_permuted)), 0.75),
                            "k" = k,
                            "m" = m,
                            "n_tot" = length(selected_samples),
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
    write.csv(res,"~/Documents/topic-modeling/r/experiments/real/results_vm_data_NC_transformed_l1.csv")
    
  #}
  
  
}

STOP



aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models,
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp <- 
  align_topics(
    models = lda_models,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)



library(tidyverse)
theme_set(theme_bw(base_size = 14))
res <- read_csv("~/Documents/topic-modeling/r/experiments/real/results_vm_data_NC_monday.csv")

res <- read_csv("~/Documents/topic-modeling/r/experiments/real/spleen_results.csv")


res[,2:ncol(res)] %>%
  group_by(method, k) %>%
  summarise(c=n())


my_colors <- c(
  "dodgerblue", "red","orchid3")

ggplot(res[,2:ncol(res)] %>%
         #filter(exp<13) %>%
         group_by(method, k,m) %>%
         summarise_all(mean) %>%
         filter(method !="TTSNC0.005", m=="all"), 
       aes(x=k, y=mean, colour=method
           ))+
  #geom_smooth(alpha=0.2, se = FALSE)+
  geom_line(position=position_dodge(width=0.5))+
  geom_point(position=position_dodge(width=0.5), size=3)+
  geom_errorbar(aes(ymin=q25, ymax=q75),
                position=position_dodge(width=0.5), alpha=0.5)+
  scale_color_manual(values = my_colors, breaks = c( "TTS0.005", "score", "lda"),
                     labels=c('TTM (this paper)', 
                              'Topic Score (Ke et al)', 
                              'LDA \n(Blei et al)')) +
  labs(colour="Method") + 
  #facet_wrap(.~ k, scales = "free_x") + 
  ylab("Average Topic Resolution") +
  xlab("Number of Topics") 

                
ggplot(res[,1:ncol(res)] %>%
         group_by(method, k, m) %>%
         summarise_all(mean) %>%
         select(-c("min", "max", "mean", "median",
                   "q25", "q75")) %>%
         pivot_longer(cols=-c("method", "k", "m")) %>%
         mutate(name= as.numeric(str_extract(name, "(?<=res)[0-9.]+"))) %>%
         drop_na(name) %>%
         mutate(name=ifelse(name> 10, 0.1 * name, name)), 
       aes(x=0.1 * name, y=value/k, colour=method, fill=method))+
  #geom_smooth(alpha=0.2,se = FALSE)+
  geom_line()+
  facet_wrap(.~ k, scales = "free") + 
  ylab("Cumulative distribution") +
  xlab("Topic Resolution")   


ggplot(res[,2:ncol(res)] %>%
         group_by(method, k, m) %>%
         summarise_all(median), 
       aes(x=k, y=median, colour=method, fill=method))+
  #geom_smooth(alpha=0.2,se = FALSE)+
  geom_line()+
  geom_errorbar(aes(ymin=q25, ymax=q75))+
  ylab("Cumulative distribution") +
  xlab("Topic Resolution") 


#geom_line(aes(y=q25), linetype="dashed")+
#geom_line(aes(y=q75), linetype="dashed")+
#geom_errorbar(aes(ymin=q25, ymax=q75))+
#geom_point(alpha=0.1)


ggplot(res[,1:ncol(res)] %>%
         group_by(method, k, m) %>%
         summarise_all(mean)  %>%
         filter(m==1000), 
       aes(x=k, y=mean, colour=method, fill=method))+
  #geom_smooth(alpha=0.2,se = FALSE)+
  geom_line()+
  #facet_wrap(m~ k, scales = "free") + 
  ylab("Mean Resolution") +
  xlab("Topics")  



test = data.frame(x = 1:ncol(X),
                  y = sort(apply(asinh(X), 2, sum), decreasing=TRUE))
library(tidyverse)
ggplot(test, aes(x=x, y=y)) + 
  geom_point() +
  scale_x_log10()+
  scale_y_log10()

