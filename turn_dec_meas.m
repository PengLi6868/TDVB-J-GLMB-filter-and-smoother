function meas = turn_dec_meas(Meas)


meas.K = Meas.K;
for i = 1:Meas.K
    for j = 1:size(Meas.Z{i},2)
        y = Meas.Z{i}(:,j);
        meas.Z{i}(1,j) = y(2)*sin(y(1));
        meas.Z{i}(2,j) = y(2)*cos(y(1));
    end
end


