clearvars -except u y

% u, y are column vectors
u = u(:);
y = y(:);

% FIR model settings
m = 200;   % Number of FIR parameters
d = 1;     % Input delay

% Basic size check
N = length(u);
if length(y) ~= N
    error('u and y must have the same length.');
end
