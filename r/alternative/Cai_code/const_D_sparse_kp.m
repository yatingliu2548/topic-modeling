function [D1,A,W,L]=const_D_sparse_kp(p,n,N,K,s)
% Generate D with sparse W where each column of W has exactly s nonzero
% elements. anchor words in A occur with probability K/p.
A=rand(p,K);
if p/100*K>=p
    L=5;
else
    L=p/100;
end
for k=1:K
   A(((k-1)*L+1):k*L,:)=0;
   A(((k-1)*L+1):k*L,k)=K/p;
   a=A(K*L+1:p,k);
   A(K*L+1:p,k)=a/sum(a)*(1-L*K/p);
end
%A=A*diag(1./vecnorm(A,1));

W=rand(K,n);
for i=1:n
    %index=randperm(K,s);
    index=randi([1,K],1,s);
    W(setdiff(1:K,index),i)=0;
end
W=W*diag(1./vecnorm(W,1));


D1=zeros(p,n);
for i=1:n
    t=mnrnd(N,W(:,i));
    for k=1:K
        num=t(k);
        D1(:,i)=D1(:,i)+mnrnd(num,A(:,k))';
    end
end

D1=D1*diag(1./vecnorm(D1,1));
