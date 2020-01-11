function d=HDISTANCE(X,Y)
% gives the Hausdorff distance between 2 finite sets X and Y,
% X, Y are arrays of column vectors
m=size(X,2);
n=size(Y,2);
if m==1 && n==1
    d=norm(X-Y);
else
    a=zeros(1,m);
    for i=1:m
        e1=zeros(1,n);
        for j=1:n
            e1(j)=norm(X(:,i)-Y(:,j));
        end
        a(i)=min(e1);
    end
    d1=max(a);
    b=zeros(1,n);
    for i=1:n
        e2=zeros(1,m);
        for j=1:m
            e2(j)=norm(Y(:,i)-X(:,j));
        end
        b(i)=min(e2);
    end
    d2=max(b);
    d=max([d1,d2]);
end
