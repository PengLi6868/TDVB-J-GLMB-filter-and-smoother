%%% ======用最少的直线覆盖所有零元素
function [y,Row,Col]=ZeroCover(C,Row,Col,Cur)
[m,n]=size(C);
temp=C(Cur,:);
flag=0;
Cur=[];
for i=1:m
    if temp(i)==-inf & Col(i)~=0
       Col(i)=0;
       Cur=[Cur,i];
       flag=1;
    end
end
if flag==0
   y=0;
else
   for i=Cur
       temp=C(:,i);
       for i=1:m
           if temp(i)==inf & Row(i)==0
              Row(i)=1;
              y=i;
              flag=0;
              [y,Row,Col]=ZeroCover(C,Row,Col,y);
              break;
           end
       end
       if flag==1
          y=0;
       end
    end
end


