function [C1,pure,Z,yp] = SVM_cone_topic(A,k,beta,sa)
%
% SVM_cone of topic model inference
% Input: A: Word-document matrix
%        k: number of topics
%        beta: regularization parameter that does not use the points whose
%        norm is almong the buttom beta quantile among all the points,
%        default as 0.1
% Output: theta: node-community matrix, theta_{ij} is the probability node i is in community j
%         B: community-community matrix, B_{ij} is the probability there is
%         an edge between a node in community i and a node in community j

% Author: Xueyu Mao 
% Email: maoxueyu@gmail.com
% Last Update: Jan 15, 2019
%
% Reference: Overlapping Clustering Models, and One (class) SVM to Bind
% Them All, in Neural Information Processing Systems, 2018. ArXiv: 1806.06945
if nargin <4 
    sa=k;
end

if nargin < 3
    beta = 0.1;
end

H = normalize_row_l1_s( A' );

H = H';

[v,~]=svds(H,k);

[M1, pure,Z,yp] = SVM_cone(v,k,0,1,beta);
%[M1,M2, pure,Z,yp] = SVM_cone(v,k,0,1,beta,sa);

C1 = normalize_row_l1_s( M1' );
C1 = C1';
%C1 = sparse(C1');

%C2 = normalize_row_l1_s( M2' );
%C2 = C2'
%C2 = sparse(C2');

end



