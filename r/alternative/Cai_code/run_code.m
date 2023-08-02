rng(200)

p=1000;
n=3000;
N=5000;
K=10;
s=5;
clambda=150

% Estimation of A and W
[dataset,result]=calculation(p,n,N,K,s,clambda)

% Confidence intervals for A and W
[CI]=compute_CI(p,n,N,K)
