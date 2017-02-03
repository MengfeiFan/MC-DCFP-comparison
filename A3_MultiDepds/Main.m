clear; clc;
%% Step 1
samplesize = 1e3; i = 0; x = 0; sum_W = 0;
tau = 0; delta_tau = 0;
t_f = zeros(samplesize,4);
phi = 0;
mu_beta = 1e-4; sigma_beta = 1e-5;
D = 5;
mu_Y_0 = 1.2;
sigma_Y = 0.2;
alpha = 1;
H = 7.5;
C = 1.5;
delta_C = .1;
lambda_0 = 2.5e-5;
number_ds = 0;
gamma_a = 1e-5;

gamma_set = [0 0 0;1e-5 0 0;1e-4 0 0;1e-3 0 0];
% gamma_set = [0 0 0;0 1e-5 0;0 1e-4 0;0 1e-3 0];
% gamma_set = [0 0 0;0 0 1e-5;0 0 1e-4;0 0 1e-3];

for k=1:4
    gamma_b = gamma_set(k,1);
    gamma_c = gamma_set(k,2);
    gamma_d = gamma_set(k,3);
    i = 0; 
    while i < samplesize
        disp([num2str(i) '/' num2str(samplesize) 'th interations'])
        x = 0; sum_W = 0; tau = 0; delta_tau = 0;mu_Y = mu_Y_0;
    %% Step 2
        beta_s = normrnd(mu_beta,sigma_beta);
        while 1
            handle_lambda = @(t) (lambda_0 + gamma_a*(phi + (beta_s + gamma_d*sum_W)*t) + gamma_b*sum_W);
            T_1 = 200/handle_lambda(tau); % let T_1 be 1 times of the current mean arrival time
            delta_tau = NHPP(handle_lambda,T_1,tau);
            tau = tau + delta_tau;
        %% Step 3
            x = phi + (beta_s + gamma_d*sum_W)*tau;
            if x > D
                i = i+1;
                t_f(i,k) = (D-phi)/(beta_s + gamma_d*sum_W);
                break;
            end
        %% Step 4
            mu_Y = mu_Y_0 + gamma_c*x;
            Y = normrnd(mu_Y,sigma_Y);
            if Y < C 
               %% Step 5
                number_ds = number_ds + 1;
                sum_W = sum_W + alpha*Y;
                if sum_W > H
                    i = i+1;
                    t_f(i,k) = tau;
                    break;
                end
            else if Y < C + delta_C
                    i = i+1;
                    t_f(i,k) = tau;
                    break;
                end
            end
        end
    end
t_f(:,k) = sort(t_f(:,k));
end
% figure;
R = zeros(samplesize,1);
for i = 1:samplesize
    R(i) = 1 - (i+1)./(samplesize+3);
end

load('Rfunction.mat','R_1','t');
plot(t,R_1,'r-');
hold on;

plot(t_f(:,1),R,'k-',t_f(:,2),R,'k--',t_f(:,3),R,'k:',t_f(:,4),R,'k-.');
xlabel('t/s')
ylabel('R(t)')
legend('\gamma_b=0','\gamma_b=0','\gamma_b=10^-^5','\gamma_b=10^-^4','\gamma_b=10^-^3','Location','Southwest')
% legend('\gamma_c=0','\gamma_c=10^-^5','\gamma_c=10^-^4','\gamma_c=10^-^3','Location','Southwest')
% legend('\gamma_d=0','\gamma_d=10^-^5','\gamma_d=10^-^4','\gamma_d=10^-^3','Location','Southwest')
axis([0,6e4,0,1])