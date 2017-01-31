% compare the results obtained by A1 & A2
clear all; clc; close all;
% Results obtained by A1: R function
load('Rfunction.mat','R_0','R_1','R_2','R_3','t');

R_fun_1 = mean(R_0);
R_fun_2 = mean(R_1);
R_fun_3 = mean(R_2);
R_fun_4 = mean(R_3);

plot(t,R_fun_1,'r-',t,R_fun_2,'r-',t,R_fun_3,'r-',t,R_fun_4,'r-');
hold on;

% Results obtained by A2: failure time
load('ft.mat','R','T');
plot(T(:,1),R(:,1),'k-',T(:,2),R(:,2),'k-',T(:,3),R(:,3),'k-',T(:,4),R(:,4),'k-');

xlabel('t/s')
ylabel('R(t)')
axis([0 7e4 0 1])
