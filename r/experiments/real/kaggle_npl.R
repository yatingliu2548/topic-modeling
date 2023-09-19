library(tm)
library(topicmodels)
library(R.matlab)
library(tidyverse)
library(reticulate)
library(tidytext)
library(dplyr)
library(matrixStats)
library(Rcpp)
library(tidyverse)
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
library(slam)

setwd("C:/Users/arvinyfw/topic-modeling")
setwd("~/Documents/topic-modeling/")
source("./r/vertex_hunting_functions.R")
source('./r/score.R')
source('./r/evaluation_metrics.r')
source("./r/score.R")


train = read_csv("C:/Users/arvinyfw/Downloads/archive/train.csv")
#test = read_csv("C://Users//建新//OneDrive - The University of Chicago//NMF//archive//test.csv")

train <- train %>%
  mutate(text = paste(TITLE, ": ", ABSTRACT, sep = ""))


seed=1234
K=6

docs <- Corpus(VectorSource(train$text))
docs <- tm_map(docs, content_transformer(tolower))
docs <- tm_map(docs, removePunctuation)
docs <- tm_map(docs, removeNumbers)
docs <- tm_map(docs, stripWhitespace)
docs <- tm_map(docs, removeWords, stopwords("en"))
library(SnowballC)
docs <- tm_map(docs, content_transformer(wordStem), language = "english")

dtm <- DocumentTermMatrix(docs)
dtm2 = as(dtm, 'sparseMatrix')
meanN=median(rowSums(dtm2))

#dtm=read.csv("C://Users//建新//OneDrive - The University of Chicago//NMF//dtm.csv")

#train=read.csv("C://Users//建新//OneDrive - The University of Chicago//NMF//research_data.csv")
columns_to_retain <- which(train$`Quantitative Biology` == 0 & train$`Quantitative Finance` ==0)

train=train[columns_to_retain, ]
train=train[columns_to_retain, ]
dtm=dtm[columns_to_retain,]

<<<<<<< HEAD
meanN=100 #mean(apply(as.data(dtm, 2, sy)))



D=t(dtm2)/meanN
=======
train$doc_word_length <- sapply(strsplit(as.character(train$text), "\\s+"), length)
meanN=mean(train$doc_word_length ) #mean(apply(as.data(dtm, 2, sy)))
D=t(dtm)/meanN
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac

lda<- LDA(dtm, k = K, control = list(seed = seed), method = 'VEM')
#ap_topics <- tidy(lda, matrix = "beta")

Ahat_lda = exp(t(lda@beta))
What_lda= t(lda@gamma)

score_recovery <- score(D, K, normalize = "norm", threshold =FALSE, 
                        max_K = min(150, min(dim(D)-1)), 
<<<<<<< HEAD
                        VHMethod="SP", estimateK=FALSE, 
                        returnW=TRUE)
=======
                        VHMethod="SP", returnW=TRUE)

p=dim(D)[1]
n=dim(D)[2]
N=meanN

#Khat_tracy = select_K(score_recovery$eigenvalues, p,n, N, method="tracy")
#Khat_olga = select_K(svd(dtm)$d, p,n, N, method="olga")
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac

Ahat_tracy=score_recovery$A_hat
What_tracy=score_recovery$W_hat

alpha=0.005
ours=score(D = D, K=K, normalize = 'TTS',
<<<<<<< HEAD
           threshold =TRUE, alpha = alpha, 
           N=meanN, max_K = min(min(dim(D))-1, 150),
           VHMethod="SP", returnW=TRUE, estimateK=FALSE)
=======
           threshold =TRUE, alpha = alpha, N=N, max_K = min(min(dim(D))-1, 150),
           VHMethod="SP", returnW=TRUE)
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac


Ahat_ours=ours$A_hat
What_ours=ours$W_hat


<<<<<<< HEAD
y_train=train[,4:(5+K-2)]
=======
y_train=train[,4:(4+K-1)]
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac
y_train=apply(y_train, 1, function(x) x / sum(x))

y_real=apply(y_train, 1, which.max)





#match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
alignment <- align_topics(exp(lda_models$k12$beta), 
                                exp(lda_models_test$k12$beta), dist = "l1",
                          do.plot=TRUE)
match_permuted_lda <- alignment$match
What_lda_permuted=alignment$B_permuted


A_normalized <- normalize_rows(y_train)
B_normalized <- normalize_rows(What_tracy)
match_tracy = data.frame(A_normalized %*% t(B_normalized))
permutation_tracy <- solve_LSAP(as.matrix(match_tracy), maximum=TRUE)

colnames(match_tracy)= 1:ncol(match_tracy)

match_tracy["X"] = 1:nrow(match_tracy)
ggplot(pivot_longer(match_tracy, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))
#match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
match_permuted_tracy <- match_tracy[, permutation_tracy]
What_tracy_permuted=What_tracy[permutation_tracy,]


A_normalized <- normalize_rows(y_train)
B_normalized <- normalize_rows(What_ours)
match_ours=  data.frame(A_normalized %*% t(B_normalized))
permutation_ours <- solve_LSAP(as.matrix(match_ours), maximum=TRUE)
colnames(match_ours)= 1:ncol(match_ours)

match_ours["X"] = 1:nrow(match_ours)
ggplot(pivot_longer(match_ours, cols=-c("X")))+
  geom_tile(aes(x=X, y=name, fill=value))
#match = data.frame((exp(lda_models$k12$beta))%*% t((exp(lda_models_test$k12$beta))))
match_permuted_ours<- match_ours[, permutation_ours]
What_ours_permuted=What_ours[permutation_ours,]

# Define a function to compute cross-entropy
cross_entropy <- function(q, p) {
  # Ensure the values are greater than 0 to avoid NaNs in the log function
  q <- pmax(q, .Machine$double.eps)
  p <- pmax(p, .Machine$double.eps)

  # Compute cross entropy
  ce <- -sum(q * log(p))
  return(ce)
}



# Apply the cross_entropy fu3nction to each pair of rows from mat1 and mat2
cross_entropy_values<- mapply(cross_entropy, split(t(What_lda_permuted), row(t(What_lda_permuted))), split(t(y_train), row(t(y_train))))
# Print the cross entropy values
print(mean(cross_entropy_values))

cross_entropy_values <- mapply(cross_entropy, split(t(What_ours_permuted), row(t(What_ours_permuted))), split(t(y_train), row(t(y_train))))
# Print the cross entropy values
print(mean(cross_entropy_values))


cross_entropy_values <- mapply(cross_entropy, split(t(What_tracy_permuted), row(t(What_tracy_permuted))), split(t(y_train), row(t(y_train))))
# Print the cross entropy values
print(mean(cross_entropy_values))


#top words
Ahat_lda_permuted=Ahat_lda[,permutation_lda]
Ahat_ours_permuted=Ahat_ours[,permutation_ours]
Ahat_tracy_permuted=Ahat_tracy[,permutation_tracy]
dictionary=dtm$dimnames$Terms
<<<<<<< HEAD


=======
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac
display_word_weighting <- function(beta, dictionary, words_to_display = 5) {

  topic_amt <- nrow(beta)
  word_df <- matrix(0, words_to_display,topic_amt)

  for (k in 1:topic_amt) {
    # Get indices of top words
    indices <- sort(beta[k,], decreasing = TRUE, index.return=TRUE)


    for (i in 1:words_to_display) {
      word_df[i,k] <- paste0(dictionary[indices$ix[i]], " (",round( beta[k,indices$ix[i]], 4), ")")

    }

  }
  return(word_df)
}

top_words_lda=data.frame(display_word_weighting(t(Ahat_lda_permuted), dictionary,10))
colnames(top_words_lda)=names(y_real)
print(top_words_lda)

top_words_ours=data.frame(display_word_weighting(t(Ahat_ours_permuted), dictionary,20))
colnames(top_words_ours)=names(y_real)
print(top_words_ours)
xtable(top_words_ours)

top_words_tracy=data.frame(display_word_weighting(t(Ahat_tracy_permuted), dictionary,20))
colnames(top_words_tracy)=names(y_real)
print(top_words_tracy)
<<<<<<< HEAD
xtable(top_words_tracy
       )
=======


#save results

cross_entropy_K<- matrix(0,words_to_display,topic_amt)
>>>>>>> 23204b47e8c81111ece034dcf6bfd6e21d053bac
