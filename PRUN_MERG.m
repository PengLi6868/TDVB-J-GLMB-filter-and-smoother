function x_prun=PRUN_MERG(x_filter,T_prun,U_merg,J_max)
I.m=[];
I.P=[];
I.w=[];
I.n=0;
m=[];
P=[];
w=[];
l=[];
%==================修剪====================
for i=1:x_filter.J
    if x_filter.w(i)>T_prun
        I.m=[I.m,x_filter.m(:,i)];
        I.P=cat(3,I.P,x_filter.P(:,:,i));
        I.w=[I.w,x_filter.w(i)];
        I.n=I.n+1;
    end
end
%================合并======================
l=0;
while I.n>0
    l=l+1;
    [a,j]=max(I.w);
    L.m=[];
    L.P=[];
    L.w=[];
    L.n=0;
    L.b=[];
    for i=1:I.n
        merg=(I.m(:,i)-I.m(:,j))'*inv(I.P(:,:,i))*...
            (I.m(:,i)-I.m(:,j));
        if merg<=U_merg
            L.m=[L.m,I.m(:,i)];
            L.P=cat(3,L.P,I.P(:,:,i));
            L.w=[L.w,I.w(i)];
            L.n=L.n+1;
            L.b=[L.b,i];
        end
    end
    %===========剔除合并项===========
    for i=1:I.n
        n=1;
        s=1;
        while n==1 && s<=L.n
            a=L.b(s);
            if i==a
                n=0;
            else
                s=s+1;
            end
        end
        if n==1
            I.m=[I.m,I.m(:,i)];
            I.P=cat(3,I.P,I.P(:,:,i));
            I.w=[I.w,I.w(i)];
        end
    end          
    I.m(:,1:I.n)=[];
    I.P(:,:,1:I.n)=[];
    I.w(1:I.n)=[];
    I.n=I.n-L.n;
    %==============合并===============
    w(l)=sum(L.w);
    m1=zeros(4,1);
    for i=1:L.n
        m1=m1+L.w(i)*L.m(:,i);
    end
    m(:,l)=m1/w(l);
    P1=zeros(4);
    for i=1:L.n
        P1=P1+L.w(i)*(L.P(:,:,i)+...
            (m(:,l)-L.m(:,i))*(m(:,l)-L.m(:,i))');
    end
    P(:,:,l)=P1/w(l);
end
%==========================================
%        判断高斯分量个数是否超过上限
%==========================================
if l>J_max
    [a,j]=max(w);
    x_prun.m=repmat(m(:,j),1,J_max);
    x_prun.P=repmat(P(:,:,j),[1,1,J_max]);
    x_prun.w=repmat(a,1,J_max);
    x_prun.J=J_max;
else
    x_prun.m=m;
    x_prun.P=P;
    x_prun.w=w;
    x_prun.J=l;
end
    
    
    
    
    