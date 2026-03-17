% function [R, h] = intcor(u, y)
% 
% N = length(u);
% 
% R = zeros(N, 1); % Initialize the correlation result vector
% h = zeros(N, 1); % Initialize the lag vector
% 
% for lag = 0:N-1
%     y_shifted = circshift(y, -lag);
%     R(lag+1) = 1/N * sum(u .* y_shifted);
%     h(lag+1) = lag; 
% end

function [R,h] = intcor(u,y)
    u = u(:);
    y = y(:);

    N = length(u);

    if length(y) ~= N
        error('u and y must have the same length');
    end

    R = zeros(N,1);
    h = (0:N-1)';

    for lag = 0:N-1
        s = 0;
        for k = 1:N
            idx = mod(k-1 - lag, N) + 1; 
            s = s + u(k) * y(idx);
        end
        R(lag+1) = s / N;
    end
end