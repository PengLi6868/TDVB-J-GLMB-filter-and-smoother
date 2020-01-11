function x=ESTIMATE(x_filter,tt,shiyan)
x=[];
for i=1:x_filter.J
    if x_filter.w(i)>0.5
        for j=1:round(x_filter.w(i))
            x=[x,x_filter.m(:,i)];
        end
    end
end
%================»­Í¼================
if tt==1
a=size(x,2);
figure(1)
c={'b.','ro','ks'};
if a~=0
    plot(x(1,:),x(3,:),c{shiyan}),hold on
end
end