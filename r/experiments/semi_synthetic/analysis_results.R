# Replace 'path/to/folder' with the actual path to your folder containing the files.
library(tidyverse)
folder_path <- "~/Documents/topic-modeling/r/experiments/semi_synthetic/results/"

# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^758", full.names = TRUE)
file_list <- c(file_list, list.files(folder_path, pattern = "^759", full.names = TRUE))
# Read all the files into a single data frame using map_dfr.
data <- map_dfr(file_list, read_csv)  # Assuming your files are CSV files. Adjust the function accordingly if they are different file types.

unique(data$n)
unique(data$N)
unique(data$p)
unique(data$n_frac)
unique(data$noise_level)

ggplot(data, aes(x=K, y=l2_A, color=method)) +
  geom_smooth()+
  facet_grid(n_frac~n)+
  scale_y_log10() 

unique(data$n)
ggplot(data %>% 
         group_by(method, K, n, n_frac)  
         %>% summarise_if(is.numeric,mean) %>% 
         arrange(K, n, n_frac, l1_A) %>%
         filter(n==250, n_frac==5, method %in% c("Ours_4","TopicScore")), aes(x=K, y=l2_A, color=method)) +
  geom_line() 


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
       aes(x=N, y=p, color=method, shape=as.factor(K), size=p)) +
  geom_point()  + 
  facet_wrap(~n)


ggplot(data %>% 
         group_by(method, K, n, n_frac)  
       %>% summarise_if(is.numeric,mean) %>% 
         arrange(K, n, n_frac, l1_A) %>%
         filter(n==100, n_frac==0.5, method %in% c("Ours_0.1","TopicScore")), aes(x=K, y=l1_A, color=method)) +
  geom_line() 

res = data %>% 
  group_by(method, K, n, n_frac, noise_level)  %>% 
  summarise_if(is.numeric,mean) %>% 
  arrange(K, n, n_frac, l1_A)

res2 = data %>% 
  group_by(method, K, n, n_frac, noise_level)  %>% 
  summarise_if(is.numeric,mean) %>% 
  arrange(K, n, n_frac, l1_A) %>%
  filter(method %in% c("Ours_0.1",TopicScore"))


