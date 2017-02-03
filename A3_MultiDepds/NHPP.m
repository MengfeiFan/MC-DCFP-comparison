function t_arrival = NHPP(handle_lambda,T,tau)
% Simulates the first arrival time of an NHPP in [0,T]
t = tau; % current time
lambda = handle_lambda(t + T); % intensity at t
while 1
    % Generate arrival time of the next event of HPP with lambda
    t = t + exprnd(1/lambda);
    % Reject or accept
    v = rand;
    if v < handle_lambda(t)/lambda % accept
        if t < T % verify if T has been reached
            t_arrival = t - tau;
            break;
        else
            T = T + 1/lambda; % generate a larger T
            lambda = handle_lambda(t + T);
            t = tau;
        end   
    end
end
