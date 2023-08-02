function [estb,GD_iters]=GradAces(y,X,b,tol,stepsize,iters)
[p,~]=size(X);
GD_iters=0;
b0=b;

for ite=1:iters
    grad=zeros(size(b0));
    for i=1:p
        if (X(i,:)*b0)==0
            factor=0;
        else
            factor=y(i)/(X(i,:)*b0);
        end
        grad=grad+factor*X(i,:)';
    end
    b1=b0+stepsize*grad;
    
    b1((b1<=0))=0;
    b1=b1/vecnorm(b1,1);
    
    if max(abs(b1-b0))<=tol
        break
    else
       b0=b1;
       GD_iters=GD_iters+1; 
    end
end
estb=b0;
GD_iters;
