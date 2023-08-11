# Replace 'path/to/folder' with the actual path to your folder containing the files.
library(tidyverse)
#folder_path <- "~/Documents/topic-modeling/r/experiments/semi_synthetic/results/"
folder_path <- "~/Documents/topic-modeling/r/experiments/synthetic/results/"
# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^synthetic_resultstest", full.names = TRUE)
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
      filter(K==100, #method != "Bing",
               method %in% c("Ours_8","Ours_4", "Ours_0.1",  "Ours_0.5","Ours_1","TopicScore")), 
       aes(x=N, y=l1_A, color=method)) +
  #geom_smooth() +# 
  geom_point() +# 
  facet_grid(n_anchors~alpha, scales="free")


unique(data$n)

ggplot(data %>% 
         group_by(method, K, n, n_frac)  
         %>% summarise_if(is.numeric,mean) %>% 
         arrange(K, n, n_frac, l1_A) %>%
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


