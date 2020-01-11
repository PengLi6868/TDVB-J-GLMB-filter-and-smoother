function R=get_R(r_noise,u_noise)

cishu=500;
R=zeros(2,2);
for i=1:cishu
    if rand>0.1
        w_noise=r_noise;
    else
        w_noise=u_noise;
    end
     y(:,i)=w_noise(:,1).*randn(2,1);
     y(:,i)=y(:,i)+w_noise(:,2);
end

m=sum(y,2)/size(y,2);

for i=1:cishu
    R=R+(y(:,i)-m)*(y(:,i)-m)';
end

R=R/cishu;
     
% figure(10)
% plot(y(1,:),y(2,:),'b.')