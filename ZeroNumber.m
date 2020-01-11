%%计算向量A中的零元素的个数
function y=ZeroNumber(A)
y=0;
for i=1:length(A)
    if A(i)==0
       y=y+1;
    end
end