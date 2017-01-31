function [ m,Y_array ] = damage_shock( M,Y_mu,Y_sigma,H )
% Generate an array of damage_shock_magnitudes
Y = normrnd(Y_mu,Y_sigma,M,1);
Y( Y>H )=[];
Y_array = Y;
m = length(Y_array);



