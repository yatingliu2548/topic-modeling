# Replace 'path/to/folder' with the actual path to your folder containing the files.
library(tidyverse)
folder_path <- "~/Documents/topic-modeling/r/experiments/semi_synthetic/results/"

# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^756", full.names = TRUE)

# Read all the files into a single data frame using map_dfr.
data <- map_dfr(file_list, read_csv)  # Assuming your files are CSV files. Adjust the function accordingly if they are different file types.

unique(data$n)
unique(data$N)
unique(data$p)
unique(data$n_frac)
unique(data$noise_level)

ggplot(data %>% 
         filter(n_frac==0.5), aes(x=K, y=l2_A, color=method)) +
  geom_smooth()+
  facet_grid(noise_level~n)+
  scale_y_log10() 


ggplot(data %>% group_by(method, K, n, n_frac)  %>% summarise_if(is.numeric,mean) %>% arrange(K, n, n_frac, l1_A) %>%
         filter(n==500, n_frac==0.5), aes(x=K, y=l2_A, color=method)) +
  geom_line() +
  scale_y_log10()


res = data %>% group_by(method, K, n, n_frac, noise_level)  %>% summarise_if(is.numeric,mean) %>% arrange(K, n, n_frac, l1_A)
