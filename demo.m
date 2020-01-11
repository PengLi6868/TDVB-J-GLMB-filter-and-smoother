%TDVB-J-GLMB Filter
%author:Peng Li
%data:2019/12/16
clear all
clc,close all

% This code is based on the code framework written by Professor BA-NGU VO: 
% http://ba-ngu.vo-au.com/publications.html

% the _common folder needs to be added to the MATLAB path
% MC_max_runs determines the Monte Carlo runs

%=============================
% Scenario  

MC_max_runs = 1;  % Monte Carlo runs
Outliers = 0;       % NO outliers, Outliers = 0, otherwise, Outliers = 1
lamda = 80;  % Poisson average rate of uniform clutter (per scan)
pd = 0.95;   % Detection probability

%=================================================
% Parameters of environment
ospa_c = 100;                                     % OSPA parameters
track_time=100;                                   % tracking time
filter_max_number = 4;                            % the number of fitlers
T=1;                                              % detection interval
q_noise=[ 1 ; 1 ];                                 % Process noise
G=[T^2/2,0
    T,0
    0,T^2/2
    0,T];                                         % Noise transfer matrix
r_noise = [10;10];                                % measurement noise
ourliers_noise = [50;50];                         % measurement noise of outliers
R = diag([r_noise(1,1)^2,r_noise(2,1)^2]);        % parametier R of the filters
range_K= [ 0 2000; 0 2000];                       % sensor detection region
pdf_K= lamda/prod(range_K(:,2)-range_K(:,1));     % uniform clutter density
ps=0.99;                                          % Survival probability
for i=1:filter_max_number
    plot_data(i).d=zeros(1,track_time);           % OSPA errors
    plot_data(i).n=zeros(1,track_time);           % the number of targets erros
    plot_data(i).t=zeros(1,track_time);           % time cost
end

%=======================================================

% start Monte Carlo runs
for MC_runs=1:MC_max_runs
    MC_runs
                                                      
    
    % generate model
    %====================================
    model = gen_model_2(q_noise, R,G,7,lamda,pd,ourliers_noise);
    if Outliers == 0
        model.outliers_flg = 0; % outliers flag off
    elseif Outliers == 1
        model.outliers_flg = 1; % outliers flag on
    end
    truth = gen_truth_ekf(model); % the the ture track of targets
    [meas1 meas2]=  gen_meas_ekf(model,truth);
    
    % plot figure when MC_max_runs =1
    %========================================

    
    % plot figure when MC_max_runs =1
    if MC_max_runs==1
        figure(2)
        hold on
        for i=1:track_time
            if size(meas1.Z{i},2)>0
                plot(meas1.Z{i}(1,:),meas1.Z{i}(2,:),'b.'),hold on
                axis([-2000 2000 0 2000]);
            end
        end
        title('measurement')
        xlabel('x position'),ylabel('y position')
    end
    
    % save target true positions for OSPA
    for i=1:track_time
        x(i).m=[];
        x(i).m=[x(i).m truth.X{i}([1,3],:)];
    end
    
    %==========================================================
    %                    start filtering
    %==========================================================
    
    % Prepare memory
    for huancun_i = 1:filter_max_number
        save_est(huancun_i). d = zeros(1,track_time); % OSPA errors
        save_est(huancun_i). n = zeros(1,track_time); % the number of targets erros
        save_est(huancun_i). t = zeros(1,track_time); % time cost
    end
    
    % filtering the data
    for filter_n= 1:filter_max_number 
        filter_n
        
        % Prepare initial model
        model = gen_model_1(q_noise, R,G,filter_n,lamda,pd);
        
        % filtering
        if filter_n == 1
            est = run_filter_j_glmb(model,meas1); % J-GLMB
            est3 = est;
        elseif filter_n == 2
            est = run_filter_TK_j_glmb(model,meas1); % STVB-J-GLMB
            est4 = est;
        elseif filter_n == 3
            est = VB_TK_Smoother( est3, meas1, model, 10 ); % J-GLMB smoother
        elseif filter_n == 4
            est = VB_TK_Smoother( est4, meas1, model, 10 ); % STVB-J-GLMB smoother
        end
        save_est(filter_n).t=est.t; % save time cost
        
        % plot figure if MC_max_runs = 1
        if MC_max_runs == 1
            for i = 1:track_time
                figure(1)
                hold on
                fuhao = {'b+','ro','g>','ks'};
                if size(est.X{i},2) > 0
                    plot(est.X{i}(1,:),est.X{i}(3,:),fuhao{filter_n});
                end
            end
        end
        
        % calculate OSPA errors
        for i = 1:track_time
            if size(est.X{i},2)>0
                d(i) = OSPA(est.X{i}([1,3],:),x(i).m,ospa_c,2);
            else
                d(i) = ospa_c;
            end
        end
        
        % save OSPA data
        save_est(filter_n).d=d;
        save_est(filter_n).n=est.N';
        
    end
    
    % save data for plot section
    for p_i = 1:filter_n
        plot_data(p_i).d=plot_data(p_i).d + save_est(p_i).d;
        plot_data(p_i).n=plot_data(p_i).n + save_est(p_i).n;
        plot_data(p_i).t=plot_data(p_i).t + save_est(p_i).t;
    end

end

%==========================================================
%                    plot figures
%==========================================================

for i=1:filter_max_number
    plot_data(i).d = plot_data(i).d / MC_max_runs;
    plot_data(i).n = plot_data(i).n / MC_max_runs;
    plot_data(i).t = plot_data(i).t / MC_max_runs;
end
c = {'-b.','-ro','-g>','-k+'};
number = truth.N;

%-----------------------------------------------------

figure(3)
plot(number),hold on
for i = 1:filter_max_number
    plot(plot_data(i).n,c{i})
end
legend('True number','J-GLMB','TDVB-J-GLMB','J-GLMB with smoothing','TDVB-J-GLMB with smoothing')
xlabel('Time (s)'),ylabel('Number of targets')
axis auto;

figure(4)
hold on
c={'b','r','g','k','m','c'};
for i = 1:filter_max_number
    plot(plot_data(i).d(1:track_time),c{i})
end
xlabel('Time (s)'),ylabel('OSPA, p=2, c=100')
legend('J-GLMB','TDVB-J-GLMB','J-GLMB with smoothing','TDVB-J-GLMB with smoothing')

figure(5)
hold on
c={'-b.','-ro','-g>','-k+','-m','-c'};
for i = 1:filter_max_number
    plot(plot_data(i).t(1:track_time),c{i})
end
xlabel('Time (s)'),ylabel('Time cost')
legend('J-GLMB','TDVB-J-GLMB','J-GLMB with smoothing','TDVB-J-GLMB with smoothing')

% average performance per scan per run
c = {'                    J-GLMB','               TDVB-J-GLMB','     J-GLMB with smoothing','TDVB-J-GLMB with smoothing'};
for i=1:filter_max_number
    the_OSPA(i) = sum(plot_data(i).d)/size(plot_data(i).d,2);
    the_time(i) = sum(plot_data(i).t/size(plot_data(i).t,2));
    disp([c{i},' : ','     OSPA is ',num2str(the_OSPA(i),4),...
        '      Time cost is ', num2str(the_time(i),4)])
end

