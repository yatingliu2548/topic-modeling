# Replace 'path/to/folder' with the actual path to your folder containing the files.
library(tidyverse)
theme_set(theme_bw(base_size = 14))
folder_path <- "~/Documents/topic-modeling/r/experiments/synthetic/results/"
# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^final_synthetic_results775", full.names = TRUE)
#file_list <- c(file_list, list.files(folder_path, pattern = "^759", full.names = TRUE))
# Read all the files into a single data frame using map_dfr.
data <- map_dfr(file_list, read_csv)  # Assuming your files are CSV files. Adjust the function accordingly if they are different file types.

unique(data$n)
unique(data$N)
unique(data$p)


ggplot(data, aes(x=K, y=l2_A, color=method)) +
  geom_smooth()+
  facet_grid(N~n)+
  scale_y_log10() 


unique(data$n_anchors)
unique(data$n)
unique(data$N)
unique(data$alpha)

unique(data$K)
ggplot(data %>% 
         
         filter(K==3, n_anchors ==0, n==100,#method != "Bing",
                method %in% c( "Ours_0.001", "Ours_0.005", "Ours_0.01",
                               "Ours_0.1","TopicScore")), 
       aes(x=N, y=l2_A, color=method, fill=method)) +
  #geom_smooth() +# 
  geom_point(aes(ymin = lower, ymax = upper), width = 0.2) +# 
  facet_wrap(.~alpha, scales="free")


unique(data$n)

res = data %>% 
  group_by(method, K, n, N, p, n_anchors) %>%
  summarise(l1_A_mean = mean(l1_A),
            l1_A_q50 = quantile(l1_A, 0.5),
            l1_A_q25 = quantile(l1_A, 0.75),
            l1_A_q75 = quantile(l1_A, 0.25),
            l2_A_mean= mean(l2_A),
            l2_A_q50 = quantile(l2_A, 0.5),
            l2_A_q25 = quantile(l2_A, 0.75),
            l2_A_q75 = quantile(l2_A, 0.25),
            Khat_huy = median(Khat_huy)
            
  )



# Custom labeller function
custom_labeller <- as_labeller(c(`0` = "Nb anchor words = 0", `1` = "Nb anchor words = 1",
                                 `1000` = "p = 1,000",
                                 `10000` = "p = 10,000",
                                 `5000` = "p = 5,000"))

legend_order <- c( "Ours_0.01","Ours_0.005","Ours_0.001", "TopicScore", "Bing")#, "AWR", "LDA")
my_colors <- c(
  "dodgerblue", "cyan", "navy", "red","green2", "yellow", "plum2")



custom_labeller2 <- as_labeller(c(`30` = "K = 30", `10` = "K = 10",
                                  `3` = "K = 3",
                                  `7` = "K = 7",
                                  `15` = "K = 15",
                                 `1000` = "p = 1,000",
                                 `10000` = "p = 10,000",
                                 `5000` = "p = 5,000"))

for (k in sort(unique(data$K))){
  for (nn in unique(data$n)){
    p<-ggplot(res %>% 
             filter(K==k, n==nn,
                    method %in%  c("Bing", "Ours_0.01","Ours_0.005","Ours_0.001",
                                  "TopicScore")), 
           aes(x=N, y=l2_A_mean, color=method)) +
      geom_line()+
      geom_errorbar(aes(ymin = l2_A_q25, ymax = l2_A_q75), width = 0.1, alpha=0.4) + 
      #scale_color_manual(values = my_colors, breaks = legend_order,
      #                   labels=c('TTM (this paper) 0.01', 'Topic Score (Ke et al)', 'Sparse Topic Estimation\n(Bing et al)')) +
      geom_point() + facet_grid(n_anchors~p, scales="free", labeller = as_labeller(custom_labeller)) + 
      scale_y_log10() +
      scale_x_log10() +
      xlab("N (Document size)") + 
      ylab(expression(italic(l[2]) ~ "error")) +
      labs(colour="Method") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("~/Documents/topic-modeling/r/experiments/synthetic/plots/n", nn, "_K", k, ".png"),
           p)
  }
}


custom_labeller2 <- as_labeller(c(`30` = "K = 30", `10` = "K = 10",
                                  `3` = "K = 3",
                                  `7` = "K = 7",
                                  `15` = "K = 15",
                                  `1000` = "p = 1,000",
                                  `10000` = "p = 10,000",
                                  `5000` = "p = 5,000"))

for (n_anchor in sort(unique(data$n_anchors))){
  for (nn in unique(data$n)){
    p<-ggplot(res %>% 
                filter(n_anchors==n_anchor, n==nn,
                       method %in%  c("Bing", "Ours_0.01","Ours_0.005","Ours_0.001",
                                      "TopicScore")), 
              aes(x=N, y=l2_A_mean, color=method)) +
      geom_line()+
      geom_errorbar(aes(ymin = l2_A_q25, ymax = l2_A_q75), width = 0.1, alpha=0.4) + 
      #scale_color_manual(values = my_colors, breaks = legend_order,
      #                   labels=c('TTM (this paper) 0.01', 'Topic Score (Ke et al)', 'Sparse Topic Estimation\n(Bing et al)')) +
      geom_point() + facet_grid(K~p, scales="free", labeller = as_labeller(custom_labeller2)) + 
      scale_y_log10() +
      scale_x_log10() +
      xlab("N (Document size)") + 
      ylab(expression(italic(l[2]) ~ "error")) +
      labs(colour="Method") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(paste0("~/Documents/topic-modeling/r/experiments/synthetic/plots/n", nn, "_n_anchor", n_anchor, ".png"),
           p)
  }
}



ggplot(res, 
       aes(x=N, y=l2_A_mean, color=method)) +
  geom_line()+
  geom_point() + facet_grid(K~p, scales="free", labeller = as_labeller(custom_labeller2)) + )

legend_order <- c("Ours_0.01", "Ours_0.05", "TopicScore", "Bing", "AWR")
my_colors <- c(
    "dodgerblue", "cyan",  "red","green2", "yellow")
    
ggplot(data %>% 
         group_by(method, K, n, n_frac)  
       %>% summarise_if(is.numeric,mean) %>% 
         filter(K==3, n %in% c(100, 250, 500, 1000, 2000),n==500, N>50,
                method %in% c("Ours_8","Ours_4", "Ours_0.1",  "Ours_0.5",
                              "Ours_1",
                              "TopicScore")), 
       aes(x=N, y=l1_A, color=method)) +
  scale_color_manual(values = my_colors, breaks = legend_order,
                     labels=c('TTM', '')) +
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


