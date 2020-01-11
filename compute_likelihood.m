function gz_vals= compute_likelihood(model,z,X)

% compute likelihood vector g= [ log_g(z|x_1), ... , log_g(z|x_M) ] -
% this is for bearings and range case with additive Gaussian noise
sa = 0;
sb = 0;
M= size(X,2);
y= X([1 3],:);
Phi= zeros(2,M);
Phi(1,:)= atan2(y(2)-sb,y(1)-sa);
Phi(2,:)= ((y(1)-sa).^2 + (y(2)-sb).^2).^0.5;
e_sq= sum( (diag(1./diag(model.D))*(repmat(z,[1 M])- Phi)).^2 );
gz_vals= exp(-e_sq/2 - log(2*pi*prod(diag(model.D))));
