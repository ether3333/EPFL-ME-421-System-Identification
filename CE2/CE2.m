clearvars -except u y

%% 2.1.1 FIR model identification
[y, u, Ts] = GetExperimentData('data.m');

% u, y are column vectors
u = u(:);
y = y(:);

% FIR model settings
m = 200;
d = 1;

% Basic size check
N = length(u);
if length(y) ~= N
    error('u and y must have the same length.');
end

%Make each row to 200 inputs
first_col = u(m:N-1);
first_row = u(m:-1:1);
Phi = toeplitz(first_col, first_row);

y_used = y(m+1:N); % follow conditions of FIR
%1
theta_hat = Phi \ y_used; % Least-square solution
%2
y_hat = Phi * theta_hat;
prediction_error = y_used - y_hat;
J = sum(prediction_error.^2);
% Plot measured output and predicted output
figure;
plot(k, y_used, 'b', 'LineWidth', 1.2);
hold on;
plot(k, y_hat, 'r--', 'LineWidth', 1.2);
hold off;

xlabel('k');
ylabel('Output');
title('Measured output and predicted output');
legend('Measured output y(k)', 'Predicted output \ity\rm_hat(k,\theta_hat)', 'Location', 'best');
grid on;
% Display the loss function
fprintf('Loss function J(theta_hat) = %.6f\n', J);

%3
% Estimate noise variance
sigma_sq_hat = J / ((N-m)-m); %N-m is #of data actually used
cov_theta_hat = sigma_sq_hat * inv(Phi'*Phi);
