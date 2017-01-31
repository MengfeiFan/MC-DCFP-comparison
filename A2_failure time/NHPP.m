function tau = NHPP(handle_lambda)
% NHPP simulates the first arrival time of an NHPP
%   Detailed explanation goes here
T = 1e3;
while 1
    lambda = handle_lambda(T);
    u = rand;
    t = -(1./lambda).*log(u);
    if t>T
        T = T + 1./lambda;
        continue;
    else
        v = rand;
        if v < handle_lambda(t)./lambda;
            tau = t;
            return;
        else
            continue;
        end
    end
end

