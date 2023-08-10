n =100
n_frac = 5
N= floor(n_frac * n)
A = NULL
W= NULL
doc_length=N
seed = 1234
vocab=NULL
noise_level =0
remove_stop_words = FALSE
Epsilon=0
dataset="AP"
K=20


data = synthetic_dataset_generation(dataset,  K, doc_length=N, n=n, seed = seed,
                                    A=A, W=W, vocab=vocab, noise_level =noise_level,
                                    remove_stop_words = remove_stop_words,
                                    Epsilon=Epsilon)


dim(data$D)
norm = sqrt(apply(data$A^2,1,sum))
sorted = sort(norm, index.return=T, decreasing=TRUE) 
test= sorted$x * (1:length(sorted$x))#### Doesn't verify the assumption really
print(max(test))


data2 = synthetic_dataset_generation(dataset,  K=3, doc_length=N, n=n, seed = seed,
                                    A=A, W=W, vocab=vocab, noise_level =noise_level,
                                    remove_stop_words = remove_stop_words,
                                    Epsilon=Epsilon)

h1 = sqrt(apply(data2$A,1,sum))
norm = sqrt(apply(data2$A^2,1,sum))
sorted = sort(norm, index.return=T, decreasing=TRUE) 
test= sorted$x * (1:length(sorted$x))#### Doesn't verify the assumption really
print(max(test))

