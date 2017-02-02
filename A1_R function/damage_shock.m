function Y = damage_shock(Y_mu,Y_sigma,C_l)
% Generate an array of damage_shock_magnitudes
% Y = normrnd(Y_mu,Y_sigma,M,1);
% Y( Y>H )=[];
% Y_array = Y;
% m = length(Y_array);

while 1
    Y = normrnd(Y_mu,Y_sigma);
    if Y < C_l % belong to damage zone
        return
    else
        continue;
    end
end


