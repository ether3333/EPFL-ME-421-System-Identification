clearvars -except u y

%%2.1.1
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
y_used = y(m+1:N);

%Make each row to 200 inputs
first_col = u(m:N-1);
first_row = u(m:-1:1);
Phi = toeplitz(first_col, first_row);

%define Theta hat
theta_hat = Phi \ y_used;

% Compute the predicted output
y_hat = Phi * theta_hat;

% Plot measured output and predicted output
figure;

% Compute the loss function
J = sum(prediction_error.^2);

% Display the loss function
fprintf('Loss function J(theta_hat) = %.6f\n', J);

