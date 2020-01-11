function [Z2 Z1]= gen_observation_fn_ekf(model,X,W)

%r/t observation equation
if model.outliers_flg == 1
    if rand<0.1
        DD = model.D_outliers;
    else
        DD = model.D;
    end
else
    DD = model.D;
end

if ~isnumeric(W)
    if strcmp(W,'noise')
        W= DD*randn(2,size(X,2));
    elseif strcmp(W,'noiseless')
        W= zeros(2,size(X,2));
    end
end

Z2= [];
if isempty(X)
    Z1= [];
    Z2= [];
else %modify below here for user specified measurement model
    P= X([1 3],:); 
    Z1(1,:) = P(1,:) + W(1,:);
    Z1(2,:) = P(2,:) + W(2,:);

end
