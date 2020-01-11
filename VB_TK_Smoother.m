function est_out = VB_TK_Smoother(est,meas,model, N)

Z_max_num = 1; % Pseudo target threshold
Time = meas.K;
[The_track,J] = Analys_track(est,meas.K); % Obtain the target trajectory from the state est.
Rinv = eye( size(model.R) ) / model.R;
for i = 1:meas.K
    X{i,1} = [];
end

% J trajectories
for j = 1:J
    m = [];
    P = [];
    ms = [];
    Ps = [];
    
    % Initialization
    %--------------------------------------
    meas.K = size(The_track(j).m,2);
    for i = 1:meas.K
        if sum( The_track(j).m(:,i) , 1 ) ~= 0
           the_first_m = The_track(j).m(:,i);
           break;
        end
    end
    if i == 1
        m0 = the_first_m;
    else
        for h = 1:size(model.m_birth,2)
            d(h) = (model.m_birth([1,3])' - the_first_m([1,3]))'*(model.m_birth([1,3])' - the_first_m([1,3]));
        end
        [~,d_i] = min(d);
        m0 = model.m_birth(:,d_i);
    end
    
    m(:,1) = m0;
    P(:,:,1) = diag([10 100 10 100 ]);
    lanmda = ones(1,meas.K);
    %---------------------------------
    % N times VB iteration
    for n = 1:N

        % Consider the pseudo target
        if size(The_track(j).Z,2) <= Z_max_num
             break;
        end
        % Estimate the state of trajectory
        for k=2:meas.K
            tic;
            if sum(The_track(j).Z(:,k),1) ~= 0
                m_pre = model.F*m(:,k-1);
                P_pre = model.F*P(:,:,k-1)*model.F' + model.Q;
                freedom_v = 8;
                [qz_temp, m_temp, P_temp] = kft_step(m_pre, P_pre, The_track(j).Z(:,k), model.H, model.R, freedom_v, 1, model,lanmda(k));
                m(:,k) = m_temp;
                P(:,:,k)= P_temp(1).P;  
            else
                m(:,k) = model.F*m(:,k-1);
                P(:,:,k) = model.F*P(:,:,k-1)*model.F'  + model.Q;
            end
            est.t(k) = est.t(k) + toc;
        end
        
        % The initial state of the smoother
        ms(:,meas.K) = m(:,meas.K);
        Ps(:,:,meas.K) = P(:,:,meas.K);
        
        % Smoothing
        for k = (meas.K - 1) :-1: 1
            tic;
            m_pre = model.F*m(:,k);
            P_pre = model.F*P(:,:,k)*model.F' + model.Q;
            G = P(:,:,k)*model.F'*inv(P_pre);
            ms(:,k) = m(:,k) + G*( ms(:,k+1) - m_pre);
            Ps(:,:,k) = P(:,:,k) - G*(  P_pre - Ps(:,:,k+1) )*G';
            est.t(k) = est.t(k) + toc;
        end
        
        % Update the auxiliary random variable
        for k = 2 : meas.K
            tic;
            if sum(The_track(j).Z(:,k),1) ~= 0
                dd = The_track(j).Z(:,k)  - model.H*ms(:,k);
                Psi = Rinv * (dd*dd' + model.H*Ps(:,:,k)*model.H');
                vp = trace(Psi);
                lanmda(k) = (freedom_v+2)/(freedom_v+vp);
            else
                lanmda(k) = lanmda(k-1);
            end
            est.t(k) = est.t(k) + toc;
        end
        
    end
%------------------------------
    % Save
    if size(The_track(j).Z,2) > Z_max_num
        for i = The_track(j).K(1):The_track(j).K(2)
            X{i,1} = [X{i,1} ms(:,i-The_track(j).K(1)+1)];
        end
    end
    
end

est_out = est;
est_out.X = X;
for k = 1:Time
    est_out.N(k,1) = size(est_out.X{k,1},2);
end
        
    
        