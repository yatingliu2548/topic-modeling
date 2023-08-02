function [LB, UB,hlf_wid,center,b]=ci(N,D,varAW,estAW,loc1,loc2,alpha)
% compute confidence interval for an entry
% estAW=hatA if we consider the confidence interval of A
% estAW=hatW if we consider the confidence interval of W

% varAW is the input required to compute the variance of estimator
% varAW=W if we aim to compute the confidence interval of A
% varAW=A if we aim to compute the confidence interval of W

dim=size(varAW,2);
i=loc1;
j=loc2;
center=estAW(i,j);

% construct the canonical basis vector
I=eye(dim);
ej=I(:,j);

zalpha=norminv(1-alpha/2);
g1=diag(1./D(:,i));
g1(find(g1==Inf))=0;

b=(ej'/(varAW'*g1*varAW)*ej);
hlf_wid=zalpha*sqrt(b/N);
LB=center-hlf_wid;
UB=center+hlf_wid;

end