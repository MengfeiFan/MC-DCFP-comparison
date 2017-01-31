clear; clc; close all;
% Load data set 1
load('Rfunction.mat','R_s0','R_s1','R_s2','R_s3','t');

% Nonparametric error estimations
R_mu_1 = mean(R_s0);
R_mu_2 = mean(R_s1);
R_mu_3 = mean(R_s2);
R_mu_4 = mean(R_s3);

plot(t,R_mu_1,'r-',t,R_mu_2,'g-',t,R_mu_3,'b-',t,R_mu_4,'c-');
hold on;

R_ci_1 = DataSeriesNonParamErr(R_s0, 2, .05);
R_ci_1 = R_ci_1';
R_1_u = R_ci_1(:,1);
R_1_l = R_ci_1(:,2);
plot(t,R_1_u,'k--',t,R_1_l,'k--')
hold on;

R_ci_2 = DataSeriesNonParamErr(R_s1, 2, .05);
R_ci_2 = R_ci_2';
R_2_u = R_ci_2(:,1);
R_2_l = R_ci_2(:,2);
plot(t,R_2_u,'k--',t,R_2_l,'k--')
hold on;

R_ci_3 = DataSeriesNonParamErr(R_s2, 2, .05);
R_ci_3 = R_ci_3';
R_3_u = R_ci_3(:,1);
R_3_l = R_ci_3(:,2);
plot(t,R_3_u,'k--',t,R_3_l,'k--')
hold on;

R_ci_4 = DataSeriesNonParamErr(R_s3, 2, .05);
R_ci_4 = R_ci_4';
R_4_u = R_ci_4(:,1);
R_4_l = R_ci_4(:,2);
plot(t,R_4_u,'k--',t,R_4_l,'k--')
hold on;

legend('\gamma_a=0',...
    '\gamma_a1=10^-^5',...
    '\gamma_a=10^-^4',...
    '\gamma_a=10^-^3',...
'95% confidence interval',...
    'Location','SouthEastOutside')

xlabel('t');
ylabel('R(t)');
axis([0 6e4 0 1])