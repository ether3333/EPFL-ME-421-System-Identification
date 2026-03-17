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
            idx = mod(k-1 + lag, N) + 1;
            s = s + u(k) * y(idx);
        end
        R(lag+1) = s / N;
    end
end