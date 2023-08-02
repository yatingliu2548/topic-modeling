# Replace 'path/to/folder' with the actual path to your folder containing the files.
folder_path <- "~/Documents/topic-modeling/r/experiments/semi_synthetic/results/"

# Get a list of filenames starting with "7547" in the folder.
file_list <- list.files(folder_path, pattern = "^755", full.names = TRUE)

# Read all the files into a single data frame using map_dfr.
data <- map_dfr(file_list, read_csv)  # Assuming your files are CSV files. Adjust the function accordingly if they are different file types.

unique(data$n)
unique(data$N)
unique(data$p)
unique(data$n_frac)

ggplot(data, aes(x=K, y=l1_A, color=method)) +
  geom_point()+
  facet_grid(n_frac~n)


ggplot(data %>% filter(n==250), aes(x=K, y=l1_A, color=method)) +
  geom_point() +
  scale_y_log10()


res = data %>% filter(n==250) %>% group_by(method, K, n, n_frac)  %>% summarise_if(is.numeric,mean) %>% arrange(K, n, n_frac, l1_A)
