%%用匈牙利法解指派问题
%%C是效益矩阵
function [y,fval,flag]=Hungary(C)
y = [];
[m,n]=size(C);
tempC=C;

%%%以下两个循环使每行每列都出现零元素
for i=1:m
    tempC(i,:)=tempC(i,:)-min(tempC(i,:));
end
for i=1:n
    tempC(:,i)=tempC(:,i)-min(tempC(:,i));
end

AssignMatrix=zeros(m,n);
tempC=TryAssign(tempC);
OneNumber=0;
for i=1:m
    for j=1:n
        if tempC(i,j)==inf
           OneNumber=OneNumber+1;
           break;
        end
    end
end
rum_time = 0;
while OneNumber<m && rum_time < 1000
      Row=zeros(m,1);
      Col=ones(1,n);
      Line=[];
      for i=1:m
          if IsInMatrix(inf,tempC(i,:))==0
             Line=[Line,i];
             Row(i)=1;
          end
      end
      for i=Line
          Cur=i;
          while Cur~=0
               [Cur,Row,Col]=ZeroCover(tempC,Row,Col,Cur);
          end
      end
      temp=inf;
      for i=1:m
          for j=1:n
              if Row(i)==1&&Col(j)==1&&tempC(i,j)<temp
                 temp=tempC(i,j);
              end
          end
      end
      for i=1:m
          for j=1:n
              if tempC(i,j)==inf|tempC(i,j)==-inf
                 tempC(i,j)=0;
              end
          end
      end
      for i=1:m
          if Row(i)==1
             tempC(i,:)=tempC(i,:)-temp;
          end
      end
      for j=1:n
          if Col(j)==0
             tempC(:,j)=tempC(:,j)+temp;
          end
      end
      tempC=TryAssign(tempC);
      OneNumber=0;
      for i=1:m
          for j=1:n
              if tempC(i,j)==inf
                 OneNumber=OneNumber+1;
                 break;
              end
          end
      end
      rum_time = rum_time + 1;
end


for i=1:m
    for j=1:n
        if tempC(i,j)==inf
           AssignMatrix(i,j)=1;
        end
    end
end
% y=AssignMatrix;
for i=1:m
    for j=1:m
        if AssignMatrix(i,j)==1
           y(i,1)=i;
           y(i,2)=j;
           break;
        end
    end
end
col=zeros(1,m);
row=zeros(1,m);
for i=1:m
    col=sum(AssignMatrix(i,:));
    row=sum(AssignMatrix(:,i));
end
if col==ones(1,m)&row==ones(1,m)
   flag=1;
else
   flag=0;
end
temp=C.*AssignMatrix;
fval=sum(temp(:));



