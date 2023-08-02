function [hatA,hat_A,W0,pure,tildeA]=find_A_MLE(D,K,beta,tol,stepsize,iters,clambda,sa)
if  nargin<8
    sa=K;
end


[p,num]=size(D);
d_theta = sqrt(K/num*sum(D,2));

S = find(d_theta>0);
H = D;
H(S,:) = bsxfun(@times, D(S,:), 1./d_theta(S));
[Xi,~]=svds(H,K);

[hat_tilde_A,pure]=SVM_cone(Xi*Xi',K,0,1,beta,clambda);
hat_A=bsxfun(@times, hat_tilde_A, d_theta);
C=normalize_row_l1_s(hat_A');
hat_A=C';

% pure covered is the anchor word set

pure=sort(pure);

D_NOR=1./sum(D,2).*D;
D_NOR_anc=D_NOR(pure,:);

hatW=zeros(size(D_NOR_anc));
for j=1:num
   y=D(:,j);
   tic;
   w=lsqnonneg(hat_A,y);
   w=w/vecnorm(w,1);
   hatW(:,j)=w;
end

PiW=sum(hatW,2);
tildeW=1./PiW.*hatW;

tildeA=zeros(size(hat_A));
tildeA0=zeros(size(hat_A));

%% Determine all the anchor words
tilde_A=normalize_row_l1_s(hat_A);
anchor=[];
for k=1:K
    anck=[];
    for anc=1:p
        idx=find(tilde_A(anc,:)>0.85);
        if (length(idx)==1)&& idx==k
           anck=[anck,anc];
        end
    end
    anchor=[anchor,anck];
    tildeA0(anck,k)=1;
    tildeA(anck,k)=1;
end
anchor

%% Compute all non-anchor rows.
GD_iters=zeros(1,p);
for j=setdiff(1:p,anchor)
   y=normalize_row_l1_s(D(j,:));
   a0=lsqnonneg(tildeW',y');
   a0=a0/vecnorm(a0,1);
   tildeA0(j,:)=a0';
   [a,Giters]=GradAces(D(j,:)',tildeW',a0,tol,stepsize,iters);
   tildeA(j,:)=a';
   GD_iters(j)=Giters;
end
GD_iters;

hatA=sum(D,2).*tildeA;
C1=normalize_row_l1_s(hatA');
hatA=C1';
W0=D_NOR_anc;