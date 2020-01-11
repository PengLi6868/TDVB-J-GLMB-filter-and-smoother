function d=WDISTANCE(X,Y)
% gives the Wasserstein distance between 2 finite sets X and Y,
% X, Y are arrays of column vectors
n = size(X,2);
m = size(Y,2);
if n==1 && m==1
    d= norm(X-Y);
else
    A = ones(1,n);
    B = eye(n);
    for j = 2:m
       A = [A zeros(size(A,1),n);  zeros(1,size(A,2)) ones(1,n)];
       B = [B eye(n)];
    end
    f =[];
    for j=1:m
        for i=1:n
            f=[f;(norm(X(:,i)-Y(:,j)))^2];
        end
    end
    [C,d] = linprog(f,[],[],[A;B],[ones(m,1)/m; ones(n,1)/n],zeros(m*n,1),inf*ones(m*n,1));
    d= sqrt(d);
end