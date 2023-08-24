# Replace 'path/to/folder' with the actual path to your folder containing the files.
library(tidyverse)
theme_set(theme_bw(base_size = 14))
folder_path <- "~/Documents/topic-modeling/r/experiments/semi_synthetic/results/"
#folder_path <- "~/Documents/topic-modeling/r/experiments/synthetic/results/"
# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^last_final", full.names = TRUE)
#file_list <- c(file_list, list.files(folder_path, pattern = "^759", full.names = TRUE))
# Read all the files into a single data frame using map_dfr.
data <- map_dfr(file_list, read_csv)  # Assuming your files are CSV files. Adjust the function accordingly if they are different file types.

unique(data$n)
unique(data$N)
unique(data$p)



data = data %>% 
  mutate(p_over_N = p/N,
         regime0 = p^2/N^(1.5),
         regime1 = p^3,
         regime2 = N * p^2,
         regime3 = N^2 * p^5,
         rate = sqrt(p * log(n)/ (N *n))
  ) %>%
  mutate(
         regime_max = max(regime1, regime2, regime3))

ggplot(data %>% filter(method %in% c("TopicScore",
                                     "Ours_0.005"),
                       K ==20), aes(x=rate * (1 +regime0 ), y=l2_A, colour=method,
                                     shape=as.factor(n))) +
  geom_point() + scale_x_log10() +scale_y_log10() +
  xlab(expression((1 + p^2/N^(3/2)) * sqrt(p * log(n)/ (N * n)))) +
  ylab(expression(italic(l[2]) ~ "error")) +
  labs(colour="Method", shape = "Corpus length (n)") +
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) 

res = data %>% 
  group_by(method, K, n, N) %>%
  summarise(l1_A_mean = mean(l1_A),
            l1_A_q50 = quantile(l1_A, 0.5),
            l1_A_q25 = quantile(l1_A, 0.25),
            l1_A_q75 = quantile(l1_A, 0.75),
            l1_A_q025 = quantile(l1_A, 0.025),
            l1_A_q975 = quantile(l1_A, 0.975),
            l2_A_mean= mean(l2_A),
            l2_A_q50 = quantile(l2_A, 0.5),
            l2_A_q25 = quantile(l2_A, 0.25),
            l2_A_q75 = quantile(l2_A, 0.75),
            l2_A_q025 = quantile(l2_A, 0.025),
            l2_A_q975 = quantile(l2_A, 0.975),
            Khat_huy = median(Khat_huy, na.rm=TRUE),
            Khat_huy_mean = mean(Khat_huy, na.rm=TRUE),
            Khat_olga = median(Khat_olga, na.rm=TRUE),
            Khat_olga_mean = mean(Khat_olga, na.rm=TRUE),
            Khat_tracy = median(Khat_tracy, na.rm=TRUE),
            Khat_tracy_mean = mean(Khat_tracy, na.rm=TRUE),
            thresholded = mean(thresholded, na.rm=TRUE),
            regime1= mean(regime1),
            regime2 = mean(regime2),
            regime3 = mean(regime3),
            p = mean(p),
            time=mean(time)
  )


res_l2  = pivot_wider(
  res %>% select(method, K, n,N, l2_A_mean), id_cols = c(K, n, N),names_from = method,
  values_from = l2_A_mean) %>%
  mutate(delta = (Ours_0.005 - TopicScore),
         delta_rel =  -100 * (Ours_0.005 - TopicScore)/TopicScore)


ggplot(res_l2 %>% 
         filter(n %in% c(250,  500, 1000, 2000), K %in% c(5, 10, 15, 20)), 
       aes(x=N, y=delta_rel, color=as.factor(K))) +
  geom_line(linewidth=1.5) +# 
  geom_point(size=2.5) + scale_x_log10()+
  xlab("N (Document size)") + 
  ylab("% Improvement upon Topic Score") +
  labs(colour="Number of\nTopics") +
  theme_bw(base_size = 16) +
  facet_grid(~n, scales="free", labeller = as_labeller(custom_labeller2)) 

  
res_l2 %>% filter( K==5, n==500 ) %>% select( K, n,N, delta, delta_rel)

unique(data$K)


unique(data$n)


legend_order <- c( "Ours_0.01","Ours_0.005","Ours_0.001", "TopicScore", "Bing")#, "AWR", "LDA")
legend_order <- c( "Ours_0.005", "TopicScore", "Bing")#, "AWR", "LDA")
custom_labeller2 <- as_labeller(c(`30` = "K = 30", `10` = "K = 10",
                                  `5` = "K = 5",
                                  `15` = "K = 15",
                                  `20` = "K = 20",
                                  `3` = "K = 3",
                                  `7` = "K = 7",
                                  `15` = "K = 15",
                                  `100` = "n = 100",
                                  `1000` = "n = 1,000",
                                  `2000` = "n = 2,000",
                                  `750` = "n = 750",
                                  `500` = "n = 500",
                                  `250` = "n = 250"
                                  ))

legend_order <-  c( "Ours_0.005", "TopicScore",  "Bing")#, "AWR", "LDA")
my_colors <-  c( "dodgerblue",  "red","chartreuse2")
labels_n <-  c('TTM (this paper)', 'Topic Score\n(Ke et al.)', 'Sparse Topic \n Estimation (Bing et al.)')

for (loss in c(1, 2)){
  for (with_bing in c(TRUE, FALSE)){
    ylab <- ifelse(loss ==1, expression(italic(l[1]) ~ "error"),
                   expression(italic(l[2]) ~ "error"))
    if (with_bing){
      legend_order <-  c( "Ours_0.005", "TopicScore",  "Bing")#, "AWR", "LDA")
      my_colors <-  c( "dodgerblue",  "red","chartreuse2")
      labels_n <-  c('TTM (this paper)', 'Topic Score\n(Ke et al.)', 'Sparse Topic Est. \n(Bing et al.)')
      method_list = c( "Ours_0.005", "Bing",
                       "TopicScore")
    }else{
      legend_order <-  c( "Ours_0.005", "TopicScore")#, "AWR", "LDA")
      my_colors <-  c( "dodgerblue",  "red")
      labels_n <-  c('TTM (this paper)', 'Topic Score\n(Ke et al.)')
      method_list = c( "Ours_0.005", 
                       "TopicScore")
    }

    title <- paste0(ifelse(loss == 1, "l1_error", "l2_error"),
                    ifelse(with_bing, "_with_bing_", "_without_bing"),
                    "k3_10_20_various_n")
    p<-ggplot(res %>% 
             filter(K %in% c(3, 10, 20), n %in% c(100, 500,  1000, 2000),
                    method %in%  method_list), 
           aes(x=N, y=l2_A_q50, color=method)) +
      geom_line()+
      geom_errorbar(aes(ymin = l2_A_q25, ymax = l2_A_q75), width = 0.1, alpha=0.4) + 
      scale_color_manual(values = my_colors, breaks = legend_order,
                        labels = labels_n) +
      geom_point() + facet_grid(K~n, scales="free", labeller = as_labeller(custom_labeller2)) + 
      scale_y_log10() +
      scale_x_log10() +
      xlab("N (Document size)") + 
      ylab(ylab) +
      labs(colour="Method") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("~/Documents/topic-modeling/r/experiments/semi_synthetic/plots/", title, ".png"), p, width = 7.5, height = 6.5, dpi = 300)
    
    
  }
}


legend_order <-  c( "Ours_0.005", "TopicScore",  "Bing")#, "AWR", "LDA")
my_colors <-  c( "dodgerblue",  "red","chartreuse2")
labels_n <-  c('TTM (this paper)', 'Topic Score\n(Ke et al.)', 'Sparse Topic Est. \n(Bing et al.)')
method_list = c( "Ours_0.005", "Bing",
                 "TopicScore")

ggplot(res %>% 
         filter(K %in% c(3, 10, 20), n %in% c(500,  1000, 2000),
                method %in%  c("Tracy", "Ours_0.005", "Bing")), 
       aes(x=N, y=100 * thresholded, colour=method)) +
  geom_line(linewidth=1.5)+
  geom_point(size=3.0) + facet_grid(K~n, scales="free", labeller = as_labeller(custom_labeller2)) + 
  scale_x_log10() +
  xlab("N (Document size)") + 
  ylab("Percentage of Words Thresholded") +
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  scale_y_continuous(labels = function(y) paste0(y,"%" ))+
  labs(colour="Method") + 
  theme_set(theme_bw(base_size = 18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


ggplot(res %>% 
         filter(K %in% c(3, 10, 15,  20), n %in% c(2000),
                N < 1000,
                method %in%  c("TopicScore", "Ours_0.005", "Bing")), 
       aes(x=N, y=time, colour=method)) +
  geom_line(linewidth=1.5)+
  geom_point(size=3.0) + facet_grid(n~K, scales="free", labeller = as_labeller(custom_labeller2)) + 
  scale_x_log10() +
  xlab("N (Document size)") + 
  ylab("Execution Time") +
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  scale_y_continuous(labels = function(y) paste0(y,"s" ))+
  labs(colour="Method") + 
  theme_set(theme_bw(base_size = 18))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



ggplot(res %>% 
         filter(K %in% c(3, 10, 20),
                method %in%  method_list), 
       aes(x=(p * p)/(n *N), y=l2_A_q50, color=method)) +
  geom_line()+
  geom_errorbar(aes(ymin = l2_A_q25, ymax = l2_A_q75), width = 0.1, alpha=0.4) + 
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels = labels_n) +
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  xlab("N (Document size)") + 
  ylab(ylab) +
  labs(colour="Method") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###


ggplot(res %>% 
         filter(#method != "Bing",
            K ==15,
                method %in% c(  "TopicScore",  "Ours_0.001",  "Ours_0.002",
                                "Ours_0.003",  "Ours_0.004",
                                "Ours_0.005",
                                "Ours_0.007",  "Ours_0.008",
                                "Ours_0.01")), 
       aes(x=N, y=l2_A, color=method)) +
  geom_line() +  
  facet_grid(K~n, scales="free") + scale_x_log10() + scale_y_log10() 


ggplot(res %>%
         filter(K==5, n %in% c(100, 500, 1000, 2000),n<2001,
                method %in% c("Ours_4",  "Ours_0.1",
                                          "Ours_1", "Ours_2", "Ours_8",
                                                 "TopicScore")), 
       aes(x=N, y=l1_A, color=method, size=thresholded)) +
  geom_point() + facet_grid(~n, scales="free") + scale_y_log10()


ggplot(data %>% 
         group_by(method, K, n, n_frac)  
       %>% summarise_if(is.numeric,mean) %>% 
         filter(K==3, n %in% c(100, 250, 500, 1000, 2000),n==500, N>50,
                method %in% c("Ours_8","Ours_4", "Ours_0.1",  "Ours_0.5",
                              "Ours_1",
                              "TopicScore")), 
       aes(x=N, y=l1_A, color=method)) +
  geom_line() + facet_grid(~n, scales="free")

ggplot(data %>% 
         group_by(method, K, n, n_frac)  %>% 
         summarise_if(is.numeric,mean) %>% 
         filter(method %in% c("Ours_4", "Ours_0.1",  "Ours_0.5",
                                      "Ours_1", "Ours_2", "Ours_8")), 
       aes(x=N, y=thresholded, color=method, shape=as.factor(K))) +
  geom_point()  + 
  facet_wrap(~n, scales = "free")

  
ggplot(data %>% 
         group_by(method, K, n, n_frac)  %>% 
         summarise_if(is.numeric,mean) %>% 
         filter(method %in% c("Ours_4", "Ours_0.1",  "Ours_0.5",
                              "Ours_1", "Ours_2", "Ours_8")), 
       aes(x=N, y=thresholded, color=method, shape=as.factor(K))) +
  geom_point()  + 
  facet_wrap(~n)


ggplot(data %>% 
         group_by(method, K, n, n_frac)  
       %>% summarise_if(is.numeric,mean) %>% 
         arrange(K, n, n_frac, l1_A) %>%
         filter(n==100, method %in% c("Ours_0.1","TopicScore")), 
       aes(x=K, y=l2_A, color=method)) +
  geom_line() +facet_wrap(~N, scales="free")


ggplot(data %>% 
         group_by(method, K, n, n_frac)  
       %>% summarise_if(is.numeric,mean) %>% 
         arrange(K, n, n_frac, l2_A) %>%
         filter(N==1000, K<30, method %in% c("Ours_0.1","Ours_4","Ours_8", "TopicScore")), 
       aes(x=K, y=l1_A, color=method)) +
  geom_line() +facet_wrap(~n, scales="free")

res = data %>% 
  group_by(method, K, n, n_frac, noise_level)  %>% 
  summarise_if(is.numeric,mean) %>% 
  arrange(K, n, n_frac, l1_A)

res2 = data %>% 
  group_by(method, K, n, n_frac, noise_level)  %>% 
  summarise_if(is.numeric,mean) %>% 
  arrange(K, n, n_frac, l1_A) %>%
  filter(method %in% c("Ours_0.1",TopicScore"))


