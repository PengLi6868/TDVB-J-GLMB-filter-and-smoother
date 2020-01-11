function Y = MEASUREMENT(H,r_noise,outliers_noise,target_position,N,pd,lamda,tt,Xtime)
%=======================get measurements=======================
a = 2000;
b = 2000;                                   %detection range [a b]

for i=1:N
    Y(i).m = []; %targets' measurements & clutters
    Y(i).n = 0;  % the number of targets'measurements & clutters
end

% clutters
for i=1:N
    clutter=poissrnd(lamda);
    for j=1:clutter
        y=[a b]'.*rand(2,1);
        Y(i).m=[Y(i).m,y];
    end
end

% targets' masurements
for i=1:N
    for j=1:size(Xtime,2)
        if target_position.n(j,i)==1
            s=rand;
            if s<pd
                if rand<0.1
                    w_noise=outliers_noise;
                else
                    w_noise=r_noise;
                end
                y = target_position.m(j).s(:,i-Xtime(1,j)+1)+w_noise(:,1).*randn(2,1);
                Y(i).m=[Y(i).m,y];
            end
        end
    end
    Y(i).n = size(Y(i).m,2);
end
