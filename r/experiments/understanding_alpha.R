
setwd("~/Documents/topic-modeling/")
source("r/experiments/semi_synthetic/synthetic_AP.R")
####
dataset = "AP"
N=25
n=100
A = NULL
W= NULL
doc_length=N
seed = 1234
vocab=NULL
noise_level =0
remove_stop_words = FALSE
Epsilon=0
dataset="AP"
K=5

res = c()
alpha = 8
for (N in c(25, 100, 500, 750, 1000, 5000, 10000)){
  print(N)
  data = synthetic_dataset_generation(dataset,  K, doc_length=N, n=n, seed = seed,
                                      A=A, W=W, vocab=vocab,
                                      remove_stop_words = remove_stop_words)
  p = ncol(data$D)
  M = rowMeans(t(data$D))
  D =  diag(sqrt(M^(-1))) %*% t(data$D)
  M = rowMeans(D)
  temp <- data.frame("M"=M,
                     "N"  =rep(N, length(t(data$D))),
                     "threshold" = alpha * sqrt(log(min(p,n))/(N *n)),
                     "p" = p,
                     "threshold_absolute" = alpha * sqrt(log(min(10473, n))/(N *n)))
  res  = rbind(res, temp)
  A = data$Aoriginal
  W = data$Woriginal
  vocab = data$original_vocab
}


ggplot(res, aes(x=M, colour=as.factor(N)))+
  geom_density() +
  geom_vline(aes(xintercept= threshold, colour=as.factor(N)))+
  scale_x_log10() + 
  theme_bw()

summary = res %>%
  group_by(N) %>%
  summarise(mean_m = mean(M),
            threshold = mean(threshold),
            threshold_absolute8 = mean(threshold_absolute),
            threshold_absolute4 = mean(threshold_absolute) * 4/8,
            threshold_absolute2 = mean(threshold_absolute) * 2/8,
            threshold_absolute1 = mean(threshold_absolute) * 1/8,
            threshold_absolute0.5 = mean(threshold_absolute) * 0.5/8,
            q01 = quantile(M, 0.01),
            q05 = quantile(M, 0.05),
            q10 = quantile(M, 0.1),
            q15 = quantile(M, 0.15),
            q25 = quantile(M, 0.25),
            q50 = quantile(M, 0.5),
            q75 = quantile(M, 0.75),
            q90 = quantile(M, 0.9),
            q95 = quantile(M, 0.95),
            )

ggplot(summary %>% filter(N<5000)) +
  geom_line(aes(x=N, y = threshold_absolute8, colour="threshold_abs"))+
  geom_line(aes(x=N, y = threshold, colour="threshold"))+
  geom_line(aes(x=N, y=q05, colour="q05"))+
  geom_line(aes(x=N, y=q10, colour="q10"))+
  geom_line(aes(x=N, y=q50, colour="q50"))+
  geom_line(aes(x=N, y=q95, colour="q95"))


for (N in c(25, 100, 500, 750, 1000, 5000, 10000)){
  print(mean(res$M[which(res$N == N)] < mean(res$threshold[which(res$N == N)])  ))
}
