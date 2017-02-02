%% MC codes for DCFP model considering degradation-shock dependence
% Algorithm 2: Generate the cdf of failure times
clear; close all; clc;

%% Parameter definition
% Soft failures
phi = 0;
beta_mu = 1e-4;
beta_sigma = 1e-5;
D = 5;

% Hard failures
Y_mu = 1.2;
Y_sigma = 0.2;
alpha = 1;
H = 7.5;
C_l = 1.5;
C_u = 1.6;
lambda_0 = 2.5e-5;

% Dependence factor
gamma_set = [0;1e-5;1e-4;1e-3];

% Sampling
n = 1e4;
T = zeros(n,4); 
R = zeros(n,4);

%% Simulation
for j = 1:4
    gamma_a = gamma_set(j,1);
    for i = 1:n
        disp([num2str(i) '/' num2str(n) ', ' num2str(j) '/4'])
        beta_s = normrnd(beta_mu,beta_sigma);
        T_1 = (D - phi)/beta_s;
        % Generate the arrival time in [0,T_1]
        handle_lambda = @(t) (lambda_0 + gamma_a*(phi + beta_s*t));
        t_arrival = NHPP(handle_lambda,T_1);
        n_shock_total = length(t_arrival);
        if isnan(t_arrival) % if the first shock arrives later than T_1
            T(i,j) = T_1;
        else
            tau = 0;
            sum_W = 0;        
        	for k = 1:n_shock_total
                tau = t_arrival(k);
                Y = normrnd(Y_mu,Y_sigma);
                if Y < C_l
                    sum_W = sum_W + alpha*Y;
                    if sum_W >= H
                        T(i,j) = tau;
                        break;
                    end
                else
                    if Y < C_u
                        T(i,j) = tau;
                        break;
                    end                    
                end
                if k == n_shock_total
                    if T(i,j) == 0
                        T(i,j) = T_1;
                    end
                end
            end
        end
    end
T(:,j) = sort(T(:,j));
end
%% Plotting
R = 1 - (1:n)./n;
plot(T(:,1),R,'k-',T(:,2),R,'k-',...
    T(:,3),R,'k-',T(:,4),R,'k-');

xlabel('t/s')
ylabel('R(t)')
axis([0 6e4 0 1])
save('ft.mat')