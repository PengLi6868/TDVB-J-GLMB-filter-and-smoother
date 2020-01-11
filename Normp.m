function D=Normp(X,Y,c,p)
%============求两个集合对应元素的p阶距离============
n=size(X,2);
m=size(Y,2);
if n~=0&&m~=0
   for i=1:n
       for j=1:m
                D(i,j)=(abs(X(1,i)-Y(1,j)))^p+(abs(X(2,i)-Y(2,j)))^p;
       end
   end
end
if n<m
   for i=(n+1):m
       for j=1:m
           D(i,j)=c^p;
       end
   end
end
if n>m
   for i=1:n
       for j=(m+1):n
           D(i,j)=c^p;
       end
   end
end    
