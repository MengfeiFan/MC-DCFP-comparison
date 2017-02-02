% compare the results obtained by A1 & A2
clear; clc; close all;
addpath('C:\Users\Zhiguo ZENG\Documents\MATLAB\MengfeiFan_RESS\MC-DCFP-comparison\A1_R function')
addpath('C:\Users\Zhiguo ZENG\Documents\MATLAB\MengfeiFan_RESS\MC-DCFP-comparison\A2_failure time')
% Results obtained by A1: R function
load('Rfunction.mat','R_0','R_1','R_2','R_3','t');

plot(t,R_0,'r-',t,R_1,'r-',t,R_2,'r-',t,R_3,'r-');
hold on;

% Results obtained by A2: failure time
load('ft.mat','R','T');
plot(T(:,1),R,'k-',T(:,2),R,'k-',T(:,3),R,'k-',T(:,4),R,'k-');

xlabel('t/s')
ylabel('R(t)')
axis([0 7e4 0 1])
