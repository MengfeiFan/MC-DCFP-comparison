%% MC codes for DCFP model considering degradation-shock dependence
% Algorithm 1: Calculate the R function

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
H = 7.5;
C_l = 1.5;
C_u = 1.6;
lambda_0 = 2.5e-5;
P_1 = normcdf(C_l,Y_mu,Y_sigma);
P_2 = normcdf(C_u,Y_mu,Y_sigma)-normcdf(C_l,Y_mu,Y_sigma);

% Dependence factor
gamma_0 = 0;
gamma_a_1 = 1e-5;
gamma_a_2 = 1e-4;
gamma_a_3 = 1e-3;

% Plotting 
num_t = 1000; 
t = linspace(0,6e4,num_t);
t = t';

% Sampling
n = 1000; 
beta_array = normrnd(beta_mu,beta_sigma,n,1);
R_s0 = zeros(n,num_t);
R_s1 = zeros(n,num_t);
R_s2 = zeros(n,num_t);
R_s3 = zeros(n,num_t);
M = 200; 
m = 150;
p_0_1 = zeros(m,1);
p_0_2 = zeros(m,1);
p_0_3 = zeros(m,1);
p_0_0 = zeros(m,1);
p_2 = zeros(m,1);

%% Simulation
for k = 1:num_t
    t_s = t(k);
    disp('k=');disp(k);
    for j = 1:n
        beta_s = beta_array(j);
        x = phi + beta_s*t_s;
        
        if x < D
            p_1 = 1; % p(E1/AB)
        else
            p_1 = 0;
        end
        
        lambda_integral_1 = lambda_0*t_s+gamma_a_1*phi*t_s+...
            .5*gamma_a_1*beta_s*(t_s^2); % dependent case 1
        lambda_integral_2 = lambda_0*t_s+gamma_a_2*phi*t_s+...
            .5*gamma_a_2*beta_s*(t_s^2); % dependent case 2
        lambda_integral_3 = lambda_0*t_s+gamma_a_3*phi*t_s+...
            .5*gamma_a_3*beta_s*(t_s^2); % dependent case 3
        lambda_integral_0 = lambda_0*t_s; % independent case
         
        p_3_1 = exp(-P_2*lambda_integral_1); % p(E3/AB) dependent case
        p_3_2 = exp(-P_2*lambda_integral_2); % p(E3/AB) dependent case
        p_3_3 = exp(-P_2*lambda_integral_3); % p(E3/AB) dependent case
        p_3_0 = exp(-P_2*lambda_integral_0); % p(E3/AB) independent case
        
        [a,Y] = damage_shock(M,Y_mu,Y_sigma,H);
        
        for i = 1:m
            p_0_1(i) = exp(-P_1*lambda_integral_1)*(P_1*lambda_integral_1)^i./factorial(i); % p(A/B) dependent case 1
            p_0_2(i) = exp(-P_1*lambda_integral_2)*(P_1*lambda_integral_2)^i./factorial(i); % p(A/B) dependent case 2
            p_0_3(i) = exp(-P_1*lambda_integral_3)*(P_1*lambda_integral_3)^i./factorial(i); % p(A/B) dependent case 3
            p_0_0(i) = exp(-P_1*lambda_integral_0)*(P_1*lambda_integral_0)^i./factorial(i); % p(A/B) independent case
            
            W_sum = 0;
            for q = 1:i
                W_sum = W_sum + alpha*Y(q);
            end
            if W_sum < H
                p_2(i) = 1;
            else
                p_2(i) = 0;
            end
                      
            R_s1(j,k) = R_s1(j,k) + p_1*p_2(i)*p_3_1*p_0_1(i);
            R_s2(j,k) = R_s2(j,k) + p_1*p_2(i)*p_3_2*p_0_2(i);
            R_s3(j,k) = R_s3(j,k) + p_1*p_2(i)*p_3_3*p_0_3(i);
            R_s0(j,k) = R_s0(j,k) + p_1*p_2(i)*p_3_0*p_0_0(i); 
        end 
        R_s1(j,k) = R_s1(j,k) + p_1*p_3_1*exp(-P_1*lambda_integral_1);
        R_s2(j,k) = R_s2(j,k) + p_1*p_3_2*exp(-P_1*lambda_integral_2);
        R_s3(j,k) = R_s3(j,k) + p_1*p_3_3*exp(-P_1*lambda_integral_3);
        R_s0(j,k) = R_s0(j,k) + p_1*p_3_0*exp(-P_1*lambda_integral_0);        
    end   
end

R_1 = mean(R_s1);
R_1 = R_1';
R_2 = mean(R_s2);
R_2 = R_2';
R_3 = mean(R_s3);
R_3 = R_3';
R_0 = mean(R_s0);
R_0 = R_0';

%% Plotting
figure
plot(t,R_0,'k-',...
    t,R_1,'k--',...
    t,R_2,'k:',...
    t,R_3,'k-.','Linewidth',1)
legend('\gamma_0=0',...
    '\gamma_1=10^-^5',...
    '\gamma_2=10^-^4',...
    '\gamma_3=10^-^3',...
    'Location','SouthWest')
xlabel('t/s')
ylabel('R(t)')
