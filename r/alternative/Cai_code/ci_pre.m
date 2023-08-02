function [hatA,tildeW]=ci_pre(p,K,anchor,D,W,A,tol,stepsize,iters)

% inputAW=W if we compute the confidence interval of A
% inputAW=A if we compute the confidence interval of W

tildeA = zeros(p,K);

tildeW = 1./sum(W,2).*W;

for i = 1:p
    a0 = normalize_row_l1_s(A(i,:));
    [a,~] = GradAces(D(i,:)',tildeW',a0',tol,stepsize,iters);
    tildeA(i,:) = a';
end
hatA = sum(D,2).*tildeA;
hatA = transpose(normalize_row_l1_s(hatA'));
end