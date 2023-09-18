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

theme_set(theme_bw(base_size = 14))


library(phyloseq)
library(tidyverse)
library(topicmodels);

library(tidytext);

meta_df <- curatedMetagenomicData::sampleMetadata

mydata <- meta_df %>%
  dplyr::filter(study_name == "YachidaS_2019") %>%
  dplyr::select(sample_id, subject_id, body_site, disease, age, gender, BMI) %>%
  dplyr::mutate(case_status = ifelse(disease == "CRC", "CRC", 
                                     ifelse(disease == "healthy", "Control", "Other"))) %>%
  dplyr::filter(case_status != "Other") %>%
  na.omit(.)

keep_id <- mydata$sample_id

crc_se_sp <- curatedMetagenomicData::sampleMetadata %>%
  dplyr::filter(sample_id %in% keep_id) %>%
  curatedMetagenomicData::returnSamples("relative_abundance", counts = TRUE)

crc_sp_df <- as.data.frame(as.matrix(assay(crc_se_sp)))

crc_sp_df <- crc_sp_df %>%
  rownames_to_column(var = "Species")

otu_tab <- crc_sp_df %>%
  tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::select(-Kingdom, -Phylum, -Class, -Order, -Family, -Genus) %>%
  dplyr::mutate(Species = gsub("s__", "", Species)) %>%
  tibble::column_to_rownames(var = "Species")

tax_tab <- crc_sp_df %>%
  dplyr::select(Species) %>%
  tidyr::separate(Species, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|") %>%
  dplyr::mutate(Kingdom = gsub("k__", "", Kingdom),
                Phylum = gsub("p__", "", Phylum),
                Class = gsub("c__", "", Class),
                Order = gsub("o__", "", Order),
                Family = gsub("f__", "", Family),
                Genus = gsub("g__", "", Genus),
                Species = gsub("s__", "", Species)) %>%
  dplyr::mutate(spec_row = Species) %>%
  tibble::column_to_rownames(var = "spec_row")

rownames(mydata) <- NULL
mydata <- mydata %>%
  tibble::column_to_rownames(var = "sample_id")

(ps <- phyloseq(sample_data(mydata),
                otu_table(otu_tab, taxa_are_rows = TRUE),
                tax_table(as.matrix(tax_tab))))


sample_data(ps)$read_count <- sample_sums(ps)
ps <- subset_samples(ps, read_count > 0)

minTotRelAbun <- 1e-6    
x <- taxa_sums(ps)
keepTaxa <- (x / sum(x)) > minTotRelAbun
(ps <- prune_taxa(keepTaxa, ps))
rm(list= ls()[!(ls() %in% c("ps"))])
(ps_g <- tax_glom(ps, taxrank = "Genus"))
taxa_names(ps_g) <- tax_table(ps_g)[, 6]
head(taxa_names(ps_g))

count_matrix <- data.frame(t(data.frame(otu_table(ps))))

data = count_matrix
print(dim(data))
lda_varying_params_lists <-  list()
for (k in 1:10) {
  lda_varying_params_lists[[paste0("k",k)]] <- list(k = k)
}
####
setwd("~/Documents/topic-modeling")
source("r/score.R")
source("r/experiments/real/run_topic_scores.R")


library(alto)
lda_models_full  <- 
  run_lda_models(
    data = data[1:nrow(data),],
    lda_varying_params_lists = lda_varying_params_lists,
    lda_fixed_params_list = list(method = "VEM"),
    dir = "metagenomics_lda_models_full/",
    reset = TRUE,
    verbose = TRUE
  )

topic_models_full <- run_topic_models(as.matrix(data), 1:nrow(data), list_params=3:30,
                                      threshold = FALSE,
                                      normalize="norm", VHMethod = 'SP',
                                      Mquantile = 0.00,
                                      alpha = 0.005, max_K=150)


alpha = 0.005

topic_models_huy_full <- run_topic_models(as.matrix(data),  1:nrow(data), list_params=3:4,
                                          threshold = FALSE,
                                          normalize="huy_not_centered", VHMethod = 'SP',
                                          Mquantile = 0.00,
                                          alpha = 0.005, max_K=150)

test = as.data.frame(exp(lda_models_full$k7$beta))
test["index"]=1:nrow(test)
ggplot(pivot_longer(test, cols=-c("index")), aes(x=index, y=name, fill=value)) + geom_tile()

sort(apply(exp(lda_models_full$k7$beta), 2, max))

ggplot(data=data.frame(x=apply(ceiling(asinh(X)), 1, sum) ))+
  geom_histogram(aes(x=x))


aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_full[1:10],
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp <- 
  align_topics(
    models = lda_models_full[1:10],
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)




aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_huy_full[1:10],
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp <- 
  align_topics(
    models = lda_models_full[1:7],
    method = "transport") 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)



topic_models_huy_TRANSF_full <- run_topic_models(ceiling(asinh(as.matrix(data))),
                                                 1:nrow(data), list_params=1:15,
                                          threshold = TRUE,
                                          normalize="huy_not_centered", VHMethod = 'SP',
                                          Mquantile = 0.00,
                                          alpha = alpha, max_K=150)


topic_models_TRANSF_full <- run_topic_models(ceiling(asinh(as.matrix(data))),
                                             1:nrow(data), list_params=3:50,
                                      threshold = FALSE,
                                      normalize="none", VHMethod = 'SP',
                                      Mquantile = 0.00,
                                      alpha = 0.005, max_K=150)
save.image("microbiome_full.RData")

for (exp in 1:10){
  train_index = sample(1:nrow(data), ceiling(nrow(data)/2))
  test_index = (1:nrow(data))[-train_index]
  
  lda_models  <- 
    run_lda_models(
      data = data[train_index,],
      lda_varying_params_lists = lda_varying_params_lists,
      lda_fixed_params_list = list(method = "VEM"),
      dir = "metagenomics_lda_models_train/",
      reset = TRUE,
      verbose = TRUE
    )
  
  lda_models_test  <- 
    run_lda_models(
      data = data[test_index,],
      lda_varying_params_lists = lda_varying_params_lists,
      lda_fixed_params_list = list(method = "VEM"),
      dir = "metagenomics_lda_models_test/",
      reset = TRUE,
      verbose = TRUE
    )
  
  
  topic_models <- run_topic_models(as.matrix(data), train_index, list_params=3:50,
                                   threshold = FALSE,
                                   normalize="none", VHMethod = 'SP',
                                   Mquantile = 0.00,
                                   alpha = 0.005, max_K=150)
  
  topic_models_test <- run_topic_models(as.matrix(data), test_index, list_params=3:50, 
                                        threshold = FALSE,
                                        normalize="none", VHMethod = 'SP',
                                        Mquantile = 0.00,
                                        alpha = 0.005, max_K=150)
  
  topic_models_huy <- run_topic_models(as.matrix(data),  train_index, list_params=3:50,
                                             threshold = TRUE,
                                             normalize="huy", VHMethod = 'SP',
                                             Mquantile = 0.00,
                                             alpha = alpha, max_K=150)
  
  
  
  topic_models_huy_test <- run_topic_models(as.matrix(data), test_index, list_params=3:50, 
                                                  threshold = TRUE,
                                                  normalize="huy", VHMethod = 'SP',
                                                  Mquantile = 0.00,
                                                  alpha = alpha, max_K=150)
  
  
  topic_models_huyNC <- run_topic_models(as.matrix(data),train_index, list_params=3:50,
                                               threshold = TRUE,
                                               normalize="huy_not_centered", VHMethod = 'SP',
                                               Mquantile = 0.00,
                                               alpha = alpha, max_K=150)
  
  
  
  topic_models_huyNC_test <- run_topic_models(as.matrix(data), test_index, list_params=3:50, 
                                                    threshold = TRUE,
                                                    normalize="huy_not_centered", VHMethod = 'SP',
                                                    Mquantile = 0.00,
                                                    alpha = alpha, max_K=150)
  
  for (k in 3:50){
    it = 1
      match = data.frame(((topic_models_huy[[paste0("k",k)]]$beta))%*% t(((topic_models_huy_test[[paste0("k",k)]]$beta))))
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
      
      match = data.frame(((topic_models_huyNC[[paste0("k",k)]]$beta))%*% t(((topic_models_huyNC_test[[paste0("k",k)]]$beta))))
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
    write_csv(res, "~/Documents/topic-modeling/r/experiments/real/microbiome_full.csv")
    save.image("microbiome_exp.RData")
  }
  

}











topic_to_compare <-list()
theme_set(theme_bw(base_size = 14))




topic_to_compare[[2]] = topic_models_full$k10
topic_to_compare[[1]] = lda_models_full$k10
topic_to_compare[[3]] = topic_models_huy_full$k10
aligned_topics_transport_comp <- 
  align_topics(
    models = topic_to_compare,
    method = "product", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)


healthy = which(mydata$case_status == "Control")
disease = which(mydata$case_status == "Disease")
test = data.frame(lda_models_full$k7$gamma)
sort(apply(as.matrix(test), 2, sum), index.return=TRUE, decreasing=TRUE)$ix
test["case_status"] = mydata$case_status
test = test %>% 
  arrange_at(vars(sapply(sort(apply(as.matrix(test[, 1:(ncol(test)-1)]), 2, sum), index.return=TRUE, decreasing=TRUE)$ix,
                                 function(x){paste0("X",x)})))
test["index"] = 1:nrow(test)

ggplot(pivot_longer(test, cols=-c("case_status", 
                                  "index")),
       aes(x=as.factor(index), y=value, fill=name)) +
  geom_bar(stat="identity", position="stack") +
  labs(y="Proportion", x="Sample", fill="Topic") +
  #scale_x_continuous(breaks=test$sample_num, labels=test$sample)+
  facet_wrap(.~case_status, scales = "free_x", ncol = 1)


ap_top_terms <-  pivot_longer(test, cols=-c("case_status", 
                                            "index")) %>%
  group_by(index, case_status) %>%
  slice_max(beta, n = 20) %>% 
  ungroup() %>%
  arrange(topic, -beta)


topic_models_huy_full$k5$gamma

aligned_topics_transport_comp <- 
     align_topics(
         models = topic_models_huy_full[1:10],
         method = "transport", reg=0.02) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE,
     color_by = "coherence")

plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE, 
     color_by = "coherence")

aligned_topics_transport_comp_score <- 
  align_topics(
    models = topic_models_full[1:10],
    method = "transport", reg=0.02) 
plot(aligned_topics_transport_comp_score, add_leaves = TRUE, label_topics = TRUE)


aligned_topics_transport_comp_lda <- 
  align_topics(
    models = lda_models_full[1:10],
    method = "transport", reg=0.02) 
plot(aligned_topics_transport_comp_lda, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_lda <- 
  align_topics(
    models = lda_models_full[1:8],
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp_lda, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_lda_trans <- 
  align_topics(
    models = lda_models_TRANS,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp_lda_trans, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_trans <- 
  align_topics(
    models = topic_models_huy_TRANSF_full[1:8],
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp_trans, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_score <- 
  align_topics(
    models = topic_models_full,
    method = "product", reg=0.0001) 
plot(aligned_topics_transport_comp_score, add_leaves = TRUE, label_topics = TRUE)

aligned_topics_transport_comp_trans <- 
  align_topics(
    models = topic_models_TRANSF_full[1:10],
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp_trans, add_leaves = TRUE, label_topics = TRUE)



aligned_topics_transport_comp <- 
  align_topics(
    models = topic_models_huyNC_full,
    method = "transport", reg=0.1) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)



aligned_topics_transport_comp <- 
  align_topics(
    models = lda_models,
    method = "transport", reg=0.01) 
plot(aligned_topics_transport_comp, add_leaves = TRUE, label_topics = TRUE)



