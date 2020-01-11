function [qz_update, m_update, P_update] = kft_step_noliner(Xm, XPm, Z, C, R, Nu, Nvb,model)

% Input:
%  x - prior mean of the state
%  P - prior covariance matrix of the state
%  y - measurement vector
%  A - state transition matrix
%  Q - process noise covariance matrix
%  C - measurement model matrix
%  R - measurement noise covariance matrix
%  Nu - vector of degrees of freedom
%  Nvb - number of VB iterations (30 by default)
% Output:
%  x - the mean of the posterior
%  P - the covariance matrix of the posterior

plength = size(Xm,2);
zlength = size(Z,2);
qz_update = zeros(plength,zlength);
m_update = zeros(model.x_dim,plength,zlength);
P_update(1).P = zeros(model.x_dim,model.x_dim,plength);




for k = 1:zlength
    for j = 1:plength
        %-------------------------------------------------------------
        xm = Xm(:,j);
        Pm = XPm(:,:,j);
        y = Z(:,k);
        
        if nargin<9 || isempty(Nvb)
            Nvb = 30;
        end      
        
        n_y = numel(y); 
        eta = gen_observation_fn_ekf(model,xm,'noiseless');
        [H_ekf,U_ekf]= ekf_update_mat(model,xm);                 % user specified function for application
        C = [];C = H_ekf;
        
        % Inverse of the shape matrix
        Rinv = eye(size(R))/R;
        % Intializing the kurtosis variable lambda
        mlambda = ones(n_y,1);
        
        % Variational Bayes iterations
        for i = 1:Nvb
            Lambdainv = diag(1./mlambda);
            Sx = C*Pm*C' + Lambdainv*R;
            Kx = Pm*C'/Sx;
            x = xm + Kx*(y-eta);
            P = Pm - Kx*Sx*Kx';
            etax = gen_observation_fn_ekf(model,x,'noiseless');
            dd = y-etax;
            Psi = Rinv * (dd*dd' + C*P*C');
            alpha = Nu + 2;
            beta = Nu + diag(Psi);
            mlambda = alpha ./ beta;
            
        end
        
        
        t1 = trace(Psi);
        t2 = (Nu+2)/(Nu+t1);
        
        
        newRinv = eye(size(Lambdainv*R))/(Lambdainv*R);
        newPinv = eye(size(Pm))/Pm;
        Pinv = eye(size(P))/P;
        vp = t1;
        lq = exp(-0.5*t2*vp)*exp((Nu/2)*log(t2))*exp(-Nu*t2/2);
        lq2 = exp(-0.5*t2*(y-etax)'*Rinv*(y-etax))*exp(-0.5*(x-xm)'*newPinv*(x-xm));
        qx = lq*lq2;
        
        
        qz_update(j,k) = qx;
        m_update(:,j,k) = x;
        P_update(k).P(:,:,j) = P;
    end
end


    