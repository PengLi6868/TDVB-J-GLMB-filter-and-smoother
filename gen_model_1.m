function model= gen_model_1(q_noise, R,G,filtrt_n,r,pd)

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

% birth parameters (Poisson birth model, multiple Gaussian components)
model.T_birth= 4;         %no. of GLMB birth terms
model.L_birth= zeros(model.T_birth,1);                                          %no of Gaussians in each GLMB birth term
model.r_birth= zeros(model.T_birth,1);                                          %prob of birth for each GLMB birth term
model.w_birth= cell(model.T_birth,1);                                           %weights of GM for each GLMB birth term
model.m_birth= cell(model.T_birth,1);                                           %means of GM for each GLMB birth term
model.B_birth= cell(model.T_birth,1);                                           %std of GM for each GLMB birth term
model.P_birth= cell(model.T_birth,1);                                           %cov of GM for each GLMB birth term

model.L_birth(1)=1;                                                             %no of Gaussians in birth term 1
model.r_birth(1)=0.03;                                                          %prob of birth
model.w_birth{1}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{1}(:,1)= [ -1500; 0; 250; 0 ];                                    %mean of Gaussians
model.B_birth{1}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{1}(:,:,1)= diag([20 100 20 100]);    %cov of Gaussians

model.L_birth(2)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(2)=0.03;                                                          %prob of birth
model.w_birth{2}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{2}(:,1)= [ -250; 0; 1000; 0 ];                                  %mean of Gaussians
model.B_birth{2}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{2}(:,:,1)= diag([20 100 20 100]);     %cov of Gaussians

model.L_birth(3)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(3)=0.03;                                                          %prob of birth
model.w_birth{3}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{3}(:,1)= [ 250; 0; 750; 0];                                  %mean of Gaussians
model.B_birth{3}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{3}(:,:,1)= diag([20 100 20 100]);     %cov of Gaussians

model.L_birth(4)=1;                                                             %no of Gaussians in birth term 2
model.r_birth(4)=0.03;                                                          %prob of birth
model.w_birth{4}(1,1)= 1;                                                       %weight of Gaussians - must be column_vector
model.m_birth{4}(:,1)= [ 1000; 0; 1500; 0];                                  %mean of Gaussians
model.B_birth{4}(:,:,1)= diag([ 10; 10; 10; 10 ]);                              %std of Gaussians
model.P_birth{4}(:,:,1)= diag([20 100 20 100]);     %cov of Gaussians

    

% observation model parameters (noisy x/y only)

model.H= [1,0,0,0
        0,0,1,0]; 
          %observation matrix

model.D = diag([ R(1,1)^0.5 R(2,2)^0.5 ]); 
model.R = R;



% detection parameters
model.P_D= pd;   %probability of detection in measurements
model.Q_D= 1-model.P_D; %probability of missed detection in measurements

% clutter parameters
model.lambda_c= r;                             %poisson average rate of uniform clutter (per scan)
model.range_c= [ 0 2000; 0 2000];      %uniform clutter region
model.pdf_c= 1/prod(model.range_c(:,2)-model.range_c(:,1)); %uniform clutter density
model.if_linear = 1;
model.outliers_flg = 0;