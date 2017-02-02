function t_arrival = NHPP(handle_lambda,T)
% Simulates the first arrival time of an NHPP in [0,T]
i = 0; % counter
t = 0; % current time
lambda = handle_lambda(T); % intensity at t
while 1
    % Generate arrival time of the next event of HPP with lambda
    t = t + exprnd(1/lambda);
    % Reject or accept
    v = rand;
    if v < handle_lambda(t)/lambda % accept
        if t < T % verify if T has been reached
            i = i+1;
            t_arrival(i) = t;
        else
            if i == 0
                t_arrival = nan;
            end
            break;
        end   
    end
end

