clear; clc;
lambda_0 = 2.5e-5;
gamma_a = 1e-5;
beta_s = 1e-4;
phi = 0;
handle_lambda = @(t) (lambda_0 + gamma_a*(phi + beta_s*t));
T = 6e4;
NS = 1e5;
n_E = zeros(NS,1);
for i = 1:NS
    disp([num2str(i) '/' num2str(NS)])
    t_arrival = NHPP(handle_lambda,T);
    if isnan(t_arrival)
        n_E(i) = 0;
    else
        n_E(i) = length(t_arrival);
    end
end
Lambda = lambda_0*T+gamma_a*phi*T +...
    .5*gamma_a*beta_s*(T^2);
handle_P_i = @(i) exp(-Lambda)*(Lambda).^i./factorial(i);
n = 0:30;
P_real = handle_P_i(n);

n_E = sort(n_E);
P_simulate = zeros(length(P_real),1);
for i = n
    P_simulate(i+1) = length(find(n_E==i))/NS;
end

plot(n,P_real,'-o',n,P_simulate,'r--d');