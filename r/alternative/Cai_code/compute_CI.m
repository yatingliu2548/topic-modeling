function [CI]=compute_CI(p,n,N,K)

load(strcat('result/dataset_n',num2str(n),'_p',num2str(p),'_N',num2str(N),'_K',num2str(K),'.mat'))
load(strcat('result/result_n',num2str(n),'_p',num2str(p),'_N',num2str(N),'_K',num2str(K),'.mat'))

ite_num=5 %100
alpha=0.05

for l=1:ite_num
    l
    dat=data{l};
    res=result{l};

    beta=dat.beta;
    clambda=dat.clambda;
    D=dat.D;
    W=dat.W;
    A=dat.A;

    hatA=res.hatA;
    anchor=res.anchor;
    hatW=res.hatW;
    E1=res.Aerr;
    E3=res.Wcolerr;

    tic;
    P=perms(1:K);
    Error1=Inf;
    for i = 1:size(P,1)
        colpmt=P(i,:);
        hatAnew=hatA(:,colpmt);
        E=sum(sum(abs(hatAnew-A)));
        if E<Error1
            hatAfinal=hatAnew;
            Error1=E;
            pmt=i
        end
    end
    t=toc
    hatAA=hatAfinal;
    Pmt=P(pmt,:)
    Error1

    hatWW=hatW(Pmt,:);
    tildeA=normalize_row_l1_s(hatAA);


    %% Compute CI
    t3=tic;

    time3=toc(t3)
    
    PiD=sum(D,2);
    PiW=sum(hatWW,2);
    
    Ainc1=0;
    Ainc2=0;
    Winc=0;
    Alength1=0;
    Alength2=0;
    Wlength=0;
    
    t4=tic;
    Acount=0;
    ACI=zeros(p,K);
    ACIUB=zeros(p,K);
    ACILB=zeros(p,K);
    for i=1:p
        j=randi(K);
        while A(i,j)==0
            j=randi(K);
        end
            [LBA2, UBA2,hlf_wid,center,b1]=ci(N,D',hatWW',hatAA,i,j,alpha);
            if hlf_wid<=0.5
                b2=PiD(i)*(tildeA(i,j)/PiW(j))^2;
                B=sqrt((b1+b2)/N);
                Ainc1=Ainc1+double((A(i,j)<=UBA2)&&(A(i,j)>=LBA2));
                Ainc2=Ainc2+double((A(i,j)<=center+norminv(1-alpha/2)*B)&&(A(i,j)>=center-norminv(1-alpha/2)*B));
                Alength1=Alength1+(UBA2-LBA2);
                Alength2=Alength2+2*norminv(1-alpha/2)*B;
                ACI(i,j)=UBA2-LBA2;
                ACIUB(i,j)=UBA2;
                ACILB(i,j)=LBA2;
                Acount=Acount+1;
            end
    end
    ACP1=Ainc1/Acount
    ALen1_ave=Alength1/Acount;
    ACP2=Ainc2/Acount
    ALen2_ave=Alength2/Acount;
    time4=toc(t4)
    
    t5=tic;
    Wcount=0;
    WCI=zeros(K,n);
    WCIUB=zeros(K,n);
    WCILB=zeros(K,n);
    for j=1:n
        i=randi(K);
        while W(i,j)==0
            i=randi(K);
        end
            [LBW, UBW,hlf_wid,~]=ci(N,D,hatAA,hatWW',j,i,alpha);
            if hlf_wid<=0.5
                Winc=Winc+double((W(i,j)<=UBW)&&(W(i,j)>=LBW));
                Wlength=Wlength+(UBW-LBW);
                WCI(i,j)=UBW-LBW;
                WCIUB(i,j)=UBW;
                WCILB(i,j)=LBW;
                Wcount=Wcount+1;
            end
    end
    WCP=Winc/Wcount
    WLen_ave=Wlength/Wcount;
    time5=toc(t5)
    
    CI{l}.centerA=hatAA;
    CI{l}.centerW=hatWW;
    CI{l}.tildeA=tildeA;
    CI{l}.alpha=alpha;
    CI{l}.CItime1=time3;
    CI{l}.CItime2=time4;
    CI{l}.CItime3=time5;
    
    CI{l}.ACI=ACI;
    CI{l}.ACIUB=ACIUB;
    CI{l}.ACILB=ACILB;
    CI{l}.Ainc1=Ainc1;
    CI{l}.Ainc2=Ainc2;
    CI{l}.ACP1=ACP1;
    CI{l}.ACP2=ACP2;
    CI{l}.Acount=Acount;
    CI{l}.LengthA1=ALen1_ave;
    CI{l}.LengthA2=ALen2_ave;
    CI{l}.WCI=WCI;
    CI{l}.WCIUB=WCIUB;
    CI{l}.WCILB=WCILB;
    CI{l}.Winc=Winc;
    CI{l}.WCP=WCP;
    CI{l}.Wcount=Wcount;
    CI{l}.LengthW=WLen_ave;
end

save([pwd,strcat('/result/CI_n',num2str(n),'_p',num2str(p),'_N',num2str(N),'_K',num2str(K),'.mat')], 'CI')


