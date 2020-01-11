function target_position=TARGET_FORM(MC_max_runs,Xtime,zongtime,Scenario)

m_num = size(Xtime,2);
target_position.n = zeros(m_num,zongtime);
for i=1:m_num
    target_position.m(i).s=[];
    target_position.n(i,Xtime(1,i):Xtime(2,i))=1;
end

% target tracks
t(1).xy = [200,1700,800,1500]';
t(2).xy = [200,1700,1500,800]';
t(3).xy = [200,1300,1500,1700]';
t(4).xy = [200,1000,800,700]';

for i=1:m_num
    target_position.m(i).s=[t(i).xy(1):(t(i).xy(2)-t(i).xy(1))/(Xtime(2,i)-Xtime(1,i)):t(i).xy(2) ; t(i).xy(3):(t(i).xy(4)-t(i).xy(3))/(Xtime(2,i)-Xtime(1,i)):t(i).xy(4)];
end


% plot figure
if MC_max_runs == 1 && Scenario == 1
    figure(1)
    plot(target_position.m(1).s(1,:),target_position.m(1).s(2,:),'k-'),hold on
    plot(target_position.m(2).s(1,:),target_position.m(2).s(2,:),'k-'),hold on
    plot(target_position.m(3).s(1,:),target_position.m(3).s(2,:),'k-'),hold on
    plot(target_position.m(4).s(1,:),target_position.m(4).s(2,:),'k-'),hold on
    for i = 1:4
        plot(target_position.m(i).s(1,1),target_position.m(i).s(2,1),'ko',target_position.m(i).s(1,size(target_position.m(i).s,2)),target_position.m(i).s(2,size(target_position.m(i).s,2)),'k>'),hold on
    end
    
    xlabel('X (m)'),ylabel('Y (m)')
    axis([0 2000 0 2000]);
end


