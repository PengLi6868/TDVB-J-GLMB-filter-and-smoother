% 对处理过的每行每列都有零元素的矩阵C进行试指派
function y=TryAssign(C)
[m,n]=size(C);
while IsInMatrix(0,C)==1
      flag=0;
      for i=1:m
          if ZeroNumber(C(i,:))==1
             temp=C(i,:);
             ZeroRowPos=find(temp==0);
             C(i,ZeroRowPos)=inf;
             temp=C(:,ZeroRowPos);
             ZeroColPos=find(temp==0);
             for j=ZeroColPos
                 C(j,ZeroRowPos)=-inf;
             end
             flag=flag+1;
          end
      end
      for i=1:n
          if ZeroNumber(C(:,i))==1
             temp=C(:,i);
             ZeroColPos=find(temp==0);
             C(ZeroColPos,i)=inf;
             temp=C(ZeroColPos,:);
             ZeroRowPos=find(temp==0);
             for j=ZeroRowPos
                 C(ZeroColPos,j)=-inf;
             end
             flag=flag+1;
          end
      end
      if flag==0
         temp=inf;
         for i=1:m
             if (ZeroNumber(C(i,:))<temp) && (ZeroNumber(C(i,:))>0)
                 temp=ZeroNumber(C(i,:));
                 ZeroRow=i;
             end
         end
         temp1=find(C(ZeroRow,:)==0);
         temp2=inf;
         for i=temp1
             if (ZeroNumber(C(:,i))<temp2) && (ZeroNumber(C(:,i))>0)
                 temp2=ZeroNumber(C(:,i));
                 ZeroCol=i;
             end
         end
         C(ZeroRow,ZeroCol)=inf;
         temp=find(C(ZeroRow,:)==0);
         for i=temp
             C(ZeroRow,i)=-inf;
         end
         temp=find(C(:,ZeroCol)==0);
         for i=temp
             C(i,ZeroCol)=-inf;
         end
      end
end
y=C;






