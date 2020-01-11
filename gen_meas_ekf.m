function [meas1 meas2]= gen_meas_ekf(model,truth)


%variables
meas1.K= truth.K;
meas2.K= truth.K;
meas1.Z= cell(truth.K,1);
meas2.Z= cell(truth.K,1);


%generate measurements
for k=1:truth.K
    if truth.N(k) > 0
        idx= find( rand(truth.N(k),1) <= model.P_D );                                            %detected target indices
        [meas2.Z{k} meas1.Z{k}]= gen_observation_fn_ekf(model,truth.X{k}(:,idx),'noise');        %single target observations if detected 
    end
    N_c= poissrnd(model.lambda_c);                                                               %number of clutter points
    C1 = repmat([4000;2000],[1 N_c]).*rand(2,N_c) - repmat([2000;0],[1 N_c]);
    C2 = zeros(2,N_c);
    C2(1,:)= atan2(C1(1,:),C1(2,:));   
    C2(2,:)= sqrt(sum(C1.^2));  
    meas1.Z{k}= [ meas1.Z{k} C1 ];                                                                %measurement is union of detections and clutter
    meas2.Z{k}= [ meas2.Z{k} C2 ];
end
    