function [data,result]=calculation(p,n,N,K,s,clambda)

tol=10^(-5);
stepsize=10^(-6);
iters=1000
ite_num=5 %100
beta=0.1

for l=1:ite_num
    l
    %% Construct observed D
    [D,A,W]=const_D_sparse_kp(p,n,N,K,s); % use const_D_sparse_r() for random anchor probabilites.
    data{l}.p=p;
    data{l}.n=n;
    data{l}.K=K;
    data{l}.N=N;
    data{l}.s=s;
    data{l}.beta=beta;
    data{l}.clambda=clambda;
    data{l}.D=D;
    data{l}.W=W;
    data{l}.A=A;
    
    %% Compute hatA
    t1=tic;
    [hatA,hat_A,Wtilde,anchor,tildeA]=find_A_MLE(D,K,beta,tol,stepsize,iters,clambda);
    time1=toc(t1)
    E1=sum(sum(abs(hatA*hatA'-A*A')))

    %% Compute hatW
    t2=tic;
    hatW=zeros(K,n);
    ERROR3=0;
    for j=1:n
        y=D(:,j);
        tic;
        w=lsqnonneg(hatA,y);
        w=w/vecnorm(w,1);
        [hatw,~]=GradAces(y,hatA,w,tol,stepsize,iters);
        hatW(:,j)=hatw;
        error3=sum(abs(hatw'*hatw-(W(:,j))'*W(:,j)));
        ERROR3=ERROR3+error3;
    end
    ERROR3/n
    time2=toc(t2)
    
    result{l}.hatA=hatA;
    result{l}.tildeA=tildeA;
    result{l}.anchor=anchor;
    result{l}.hatW=hatW;
    result{l}.Atime=time1;
    result{l}.Wtime=time2;
    result{l}.Aerr=E1;
    result{l}.Wcolerr=ERROR3/n;

end

save([pwd,strcat('/result/dataset_n',num2str(n),'_p',num2str(p),'_N',num2str(N),'_K',num2str(K),'.mat')], 'data')
save([pwd,strcat('/result/result_n',num2str(n),'_p',num2str(p),'_N',num2str(N),'_K',num2str(K),'.mat')], 'result')

