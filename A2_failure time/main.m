%% MC codes for DCFP model considering degradation-shock dependence
% Algorithm 2: Generate the cdf of failure times

clear all; close; clc;

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
H = 17.5;
C_l = 1.5;
C_u = 1.6;
lambda_0 = 2.5e-5;

% Dependence factor
gamma_set = [0;1e-5;1e-4;1e-3];

% Sampling
n = 1000;
M = 200;
T = zeros(n,4); 
R = zeros(n,4);

%% Simulation
for j = 1:4
    gamma_a = gamma_set(j,1);
    for i = 1:n
        beta_s = normrnd(beta_mu,beta_sigma);
        T_1 = (D - phi)/beta_s;

        tau = 0;
        sum_W = 0;
        while 1
        handle_lambda = @(t) (lambda_0 + gamma_a*(phi + beta_s*t));
        delta_tau = NHPP(handle_lambda);
        tau = tau + delta_tau;
        
        if tau >= T_1
            T(i,j) = T_1;
            break
        end
        
        Y = normrnd(Y_mu,Y_sigma);
        if Y < C_l
            sum_W = sum_W + alpha*Y;
            if sum_W >= H
                T(i,j) = tau;
                break
            end
        elseif Y < C_u
            T(i,j) = tau;
            break
        end
        
        end
        R(i,j) = 1 - (i+1)./(n+3);
    end
T(:,j) = sort(T(:,j));
end

save('ft.mat');

%% Plotting
plot(T(:,1),R(:,1),'k-',T(:,2),R(:,2),'k-',...
    T(:,3),R(:,3),'k-',T(:,4),R(:,4),'k-');

xlabel('t/s')
ylabel('R(t)')
axis([0 6e4 0 1])