function meas = turn_plor_ekf(Meas)

meas.K = Meas.K;
for i = 1:Meas.K
    for j = 1:size(Meas.Z{i},2)
        y = Meas.Z{i}(:,j);
        meas.Z{i}(1,j)= atan2(y(1),y(2));
        meas.Z{i}(2,j)= ((y(1))^2 + (y(2))^2)^0.5;
    end
end
