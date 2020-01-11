function d_hun=OSPA(X,Y,c,p)

m=size(Y,2);
n=size(X,2);
D=Normp(X,Y,c,p);
[assign_h,fval,flag]=Hungary(D);
N=max(m,n);
d_hun=(fval/N)^(1/p);
if d_hun>c
    d_hun=c;
end