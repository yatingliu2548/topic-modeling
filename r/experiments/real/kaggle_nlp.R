library(tm)
library(tidyverse)
test = read_csv("~/Downloads/archive (4)/train.csv")

corpus <- Corpus(VectorSource(test$ABSTRACT))
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, stopwords("en"))

dtm <- DocumentTermMatrix(corpus)


term_freq <- rowSums(as.matrix(dtm))
term_freq <- subset(term_freq, term_freq >= 5)  # Filter terms; this is optional

df <- data.frame(term = names(term_freq), freq = unname(term_freq))
df <- df[order(-df$freq), ]
install.packages("ggplot2")
library(ggplot2)

ggplot(df, aes(x = reorder(term, freq), y = freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Term Frequencies", x = "Terms", y = "Frequency") +
  coord_flip()
The coord_flip() function is used to flip the x and y axes, making it easier to read the terms.

Remember to adjust the plotting parameters, such as the number of terms you want to display, the ordering, axis labels, and other aesthetics to fit your data and preferences.







df <- data.frame(term = names(term_freq), freq = unname(term_freq))
df <- df[order(-df$freq), ]

corpus <- Corpus(VectorSource(test$ABSTRACT))
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeNumbers)
corpus <- tm_map(corpus, removeWords, stopwords("en"))
dtm <- DocumentTermMatrix(corpus)

df <- data.frame(term = names(term_freq), freq = unname(term_freq))
df <- df[order(-df$freq), ]
ggplot(df, aes(x = reorder(term, freq), y = freq)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Term Frequencies", x = "Terms", y = "Frequency") +
  coord_flip() +
  scale_y_log10()



X = score(t(as.matrix(dtm)), 6, scatterplot=FALSE, K0=NULL, m=NULL, N=NULL, 
          threshold=TRUE,
                  Mquantile=0.00, VHMethod = 'SP', normalize="TTS",
                  alpha=0.005, max_K=150, returnW=FALSE, estimateK=FALSE)

