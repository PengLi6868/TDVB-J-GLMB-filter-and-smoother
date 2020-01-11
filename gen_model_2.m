function model= gen_model_2(q_noise, R,G,filter_n,r,pd,ourliers_noise)


if filter_n == 7 || filter_n == 5
    % basic parameters
    model.x_dim= 5;   %dimension of state vector
    model.z_dim= 2;   %dimension of observation vector
    model.v_dim= 3;   %dimension of process noise
    model.w_dim= 2;   %dimension of observation noise
    
    % dynamical model parameters (CT model)
    % state transformation given by gen_newstate_fn, transition matrix is N/A in non-linear case
    model.T= 1;                         %sampling period
    model.sigma_vel= 5;
    model.sigma_turn= (pi/180);   %std. of turn rate variation (rad/s)
    model.bt= model.sigma_vel*[ (model.T^2)/2; model.T ];
    model.B2= [ model.bt zeros(2,2); zeros(2,1) model.bt zeros(2,1); zeros(1,2) model.T*model.sigma_turn ];
    model.B= eye(model.v_dim);
    model.Q= model.B*model.B';
    
    % survival/death parameters
    model.P_S= .99;
    model.Q_S= 1-model.P_S;
    
    % birth parameters (Poisson birth model, multiple Gaussian components)
    model.L_birth= 4;                                                     %no. of Gaussian birth terms
    model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
    model.m_birth= zeros(model.x_dim,model.L_birth);                      %means of Gaussian birth terms
    model.B_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %std of Gaussian birth terms
    model.P_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %cov of Gaussian birth terms
    
    model.w_birth(1)= 0.01;                                              %birth term 1
    model.m_birth(:,1)= [ -1500; 0; 250; 0; 0 ];
    model.B_birth(:,:,1)= diag([ 5; 30; 5; 30; 6*(pi/180) ]);
    model.P_birth(:,:,1)= model.B_birth(:,:,1)*model.B_birth(:,:,1)';
    
    model.w_birth(2)= 0.01;                                              %birth term 2
    model.m_birth(:,2)= [ -250; 0; 1000; 0; 0 ];
    model.B_birth(:,:,2)= diag([ 5; 30; 5; 30; 6*(pi/180) ]);
    model.P_birth(:,:,2)= model.B_birth(:,:,2)*model.B_birth(:,:,2)';
    
    model.w_birth(3)= 0.01;                                              %birth term 3
    model.m_birth(:,3)= [ 250; 0; 750; 0; 0 ];
    model.B_birth(:,:,3)= diag([5; 30; 5; 30; 6*(pi/180) ]);
    model.P_birth(:,:,3)= model.B_birth(:,:,3)*model.B_birth(:,:,3)';
    
    model.w_birth(4)= 0.01;                                              %birth term 4
    model.m_birth(:,4)= [ 1000; 0; 1500; 0; 0 ];
    model.B_birth(:,:,4)= diag([ 5; 30; 5; 30; 6*(pi/180) ]);
    model.P_birth(:,:,4)= model.B_birth(:,:,4)*model.B_birth(:,:,4)';
    
    % observation model parameters (noisy r/theta only)
    % measurement transformation given by gen_observation_fn, observation matrix is N/A in non-linear case
    model.D= diag([ (pi/90); 10 ]);      %std for angle and range noise
    model.R= model.D*model.D';              %covariance for observation noise
    
    
    % detection parameters
    model.P_D= pd;   %probability of detection in measurements
    model.Q_D= 1-model.P_D; %probability of missed detection in measurements
    
    % clutter parameters
    model.lambda_c= r;                             %poisson average rate of uniform clutter (per scan)
    model.range_c= [ -pi/2 pi/2; 0 2000 ];          %uniform clutter on r/theta
    model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
    
    if filter_n == 5
        model.range_c= [ -2000 2000; 0 2000];      %uniform clutter region
        model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
        model.H= [1,0,0,0,0
            0,0,1,0,0];
    end
    
elseif filter_n == 1 || filter_n == 2 || filter_n == 4 || filter_n == 6
    % basic parameters
    model.x_dim= 4;   %dimension of state vector
    model.z_dim= 2;   %dimension of observation vector
    
    % dynamical model parameters (CV model)
    model.T= 1;                                     %sampling period
    model.A0= [ 1 model.T; 0 1 ];                         %transition matrix
    model.F= [1,model.T,0,0
        0,1,0,0
        0,0,1,model.T
        0,0,0,1];
    model.B0= [ (model.T^2)/2; model.T ];
    model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
    model.sigma_v = 5;
    model.Q= G*diag([q_noise(1),q_noise(2)])*G';%(model.sigma_v)^2* model.B*model.B';   %process noise covariance
    
    % survival/death parameters
    model.P_S= .99;
    model.Q_S= 1-model.P_S;
    model.L_birth= 4;                                                     %no. of Gaussian birth terms
    model.w_birth= zeros(model.L_birth,1);                                %weights of Gaussian birth terms (per scan) [sum gives average rate of target birth per scan]
    model.m_birth= zeros(model.x_dim,model.L_birth);                      %means of Gaussian birth terms
    model.B_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %std of Gaussian birth terms
    model.P_birth= zeros(model.x_dim,model.x_dim,model.L_birth);          %cov of Gaussian birth terms
    
    model.w_birth(1)= 0.01;                                              %birth term 1
    model.m_birth(:,1)= [ -1500; 0; 250; 0];
    model.B_birth(:,:,1)= diag([ 5; 30; 5; 30]);
    model.P_birth(:,:,1)= model.B_birth(:,:,1)*model.B_birth(:,:,1)';
    
    model.w_birth(2)= 0.01;                                              %birth term 2
    model.m_birth(:,2)= [ -250; 0; 1000; 0 ];
    model.B_birth(:,:,2)= diag([ 5; 30; 5; 30 ]);
    model.P_birth(:,:,2)= model.B_birth(:,:,2)*model.B_birth(:,:,2)';
    
    model.w_birth(3)= 0.01;                                              %birth term 3
    model.m_birth(:,3)= [ 250; 0; 750; 0];
    model.B_birth(:,:,3)= diag([ 5; 30; 5; 30]);
    model.P_birth(:,:,3)= model.B_birth(:,:,3)*model.B_birth(:,:,3)';
    
    model.w_birth(4)= 0.01;                                              %birth term 4
    model.m_birth(:,4)= [ 1000; 0; 1500; 0];
    model.B_birth(:,:,4)= diag([ 5; 30; 5; 30 ]);
    model.P_birth(:,:,4)= model.B_birth(:,:,4)*model.B_birth(:,:,4)';
    
    model.H= [1,0,0,0
        0,0,1,0];
    %observation matrix
    
    model.D = diag([ R(1,1) R(2,2) ]);
    model.R = R;
    model.D_nonliner = [ (pi/90); 10 ];
    
    
    % detection parameters
    model.P_D= pd;   %probability of detection in measurements
    model.Q_D= 1-model.P_D; %probability of missed detection in measurements
    
    % clutter parameters
    model.lambda_c= r;                             %poisson average rate of uniform clutter (per scan)
    model.range_c= [ -2000 2000; 0 2000];      %uniform clutter region
    model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
elseif  filter_n == 3 
    % basic parameters
    model.x_dim= 4;   %dimension of state vector
    model.z_dim= 2;   %dimension of observation vector
    
    % dynamical model parameters (CV model)
    model.T= 1;                                     %sampling period
    model.A0= [ 1 model.T; 0 1 ];                         %transition matrix
    model.F= [1,model.T,0,0
        0,1,0,0
        0,0,1,model.T
        0,0,0,1];
    model.B0= [ (model.T^2)/2; model.T ];
    model.B= [ model.B0 zeros(2,1); zeros(2,1) model.B0 ];
    model.sigma_v = 5;
    model.Q= G*diag([q_noise(1),q_noise(2)])*G';%(model.sigma_v)^2* model.B*model.B';   %process noise covariance
    
    % survival/death parameters
    model.P_S= .99;
    model.Q_S= 1-model.P_S;
    
    model.T_birth= 4;         %no. of LMB birth terms
    model.L_birth= zeros(model.T_birth,1);                                          %no of Gaussians in each LMB birth term
    model.r_birth= zeros(model.T_birth,1);                                          %prob of birth for each LMB birth term
    model.w_birth= cell(model.T_birth,1);                                           %weights of GM for each LMB birth term
    model.m_birth= cell(model.T_birth,1);                                           %means of GM for each LMB birth term
    model.B_birth= cell(model.T_birth,1);                                           %std of GM for each LMB birth term
    model.P_birth= cell(model.T_birth,1);                                           %cov of GM for each LMB birth term
    
    model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
    model.r_birth(1)=0.01;                                                          %prob of birth
    model.w_birth{1}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
    model.m_birth{1}(:,1)= [ -1500; 0; 250; 0];                             %mean of Gaussians
    model.B_birth{1}(:,:,1)= diag([ 5; 30; 5; 30 ]);                           %std of Gaussians
    model.P_birth{1}(:,:,1)= model.B_birth{1}(:,:,1)*model.B_birth{1}(:,:,1)';   %cov of Gaussians
    
    model.L_birth(2)=1;                                                             %no of Gaussians in birth term 2
    model.r_birth(2)=0.01;                                                          %prob of birth
    model.w_birth{2}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
    model.m_birth{2}(:,1)= [ -250; 0; 1000; 0 ];                                  %mean of Gaussians
    model.B_birth{2}(:,:,1)= diag([ 5; 30; 5; 30 ]);                              %std of Gaussians
    model.P_birth{2}(:,:,1)= model.B_birth{2}(:,:,1)*model.B_birth{2}(:,:,1)';     %cov of Gaussians
    
    model.L_birth(3)=1;                                                             %no of Gaussians in birth term 2
    model.r_birth(3)=0.01;                                                          %prob of birth
    model.w_birth{3}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
    model.m_birth{3}(:,1)= [ 250; 0; 750; 0];                              %mean of Gaussians
    model.B_birth{3}(:,:,1)= diag([ 5; 30; 5; 30 ]);                            %std of Gaussians
    model.P_birth{3}(:,:,1)= model.B_birth{3}(:,:,1)*model.B_birth{3}(:,:,1)';     %cov of Gaussians
    
    model.L_birth(4)=1;                                                             %no of Gaussians in birth term 2
    model.r_birth(4)=0.01;                                                          %prob of birth
    model.w_birth{4}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
    model.m_birth{4}(:,1)= [ 1000; 0; 1500; 0];                              %mean of Gaussians
    model.B_birth{4}(:,:,1)= diag([ 5; 30; 5; 30 ]);                            %std of Gaussians
    model.P_birth{4}(:,:,1)= model.B_birth{4}(:,:,1)*model.B_birth{4}(:,:,1)';      %cov of Gaussians
    
    model.H= [1,0,0,0
        0,0,1,0];
    %observation matrix
    
    model.D = diag([ R(1,1) R(2,2) ]);
    model.R = R;
    
    
    
    % detection parameters
    model.P_D= pd;   %probability of detection in measurements
    model.Q_D= 1-model.P_D; %probability of missed detection in measurements
    
    % clutter parameters
    model.lambda_c= r;                             %poisson average rate of uniform clutter (per scan)
    model.range_c= [ -2000 2000; 0 2000];      %uniform clutter region
    model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
end
model.if_linear = 2;
model.D = R(1,1)^0.5;
model.D_outliers = ourliers_noise(1,1); %use to get the outliers