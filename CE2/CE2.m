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
k=(m+1:N); %defined for only used indexes; toeplitz
figure;
plot(k, y_used, 'b');
hold on;
plot(k, y_hat, 'r');
legend('True y', 'Estimated y')
hold off;

xlabel('index m+1 to N');
ylabel('output Y');
title('Measured output and predicted output');
grid on;
% Display the loss function
fprintf('Loss function J(theta_hat) = %.6f\n', J);

%3
% Estimate noise variance
sigma_sq_hat = J / ((N-m)-m); %N-m is #of data actually used <<Not really sure
cov_theta_hat = sigma_sq_hat * inv(Phi'*Phi);

theta_std = sqrt(diag(cov_theta_hat)); %compute the lenght of interval (σ)
%plot
figure;
errorbar(1:m, theta_hat, 2*theta_std, 'o');
xlabel('Coefficient index(m)');
ylabel('Estimated FIR');
title('FIR with ±2σ confidence interval');
grid on;


%% 2.1.2 ARX model identification

N = length(u);
Phi = [-y(2:N-1), -y(1:N-2), u(2:N-1), u(1:N-2)]; %k=3,4,...N
y_used = y(3:N);
%1
theta_hat = Phi\ y_used;
%2
y_hat = Phi * theta_hat;
prediction_error = y_used - y_hat;
J = sum(prediction_error.^2);
% Plot measured output and predicted output
k = (3:N)';
figure;
plot(k, y_used, 'b', 'LineWidth',1.5);
hold on;
plot(k, y_hat, 'r');
hold off;

%3
t = (0:N-1)' * Ts; %time vector; generated for lsim
sys_hat = tf([0 theta_hat(3) theta_hat(4)], [1 theta_hat(1) theta_hat(2)], Ts);
y_m = lsim(sys_hat, u, t);
y_m_used = y_m(3:N); %can only compare from k=3

error = y_used - y_m_used;
error_norm = norm(error, 2);

figure;
plot(k, y_used, 'b');
hold on;
plot(k, y_m_used, 'r');
hold off;
title('Measured output y and y_m');
grid on;

fprintf('Two-norm of the error = %.6f\n', error_norm);
%4 Instrumental Varibale method

% Phi_iv = [-y_m(k-1), -y_m(k-2), u(k-1), u(k-2)];
% theta_iv = (Phi_iv'*Phi)\(Phi_iv'*y_used);
% 
% y_hat_iv = Phi * theta_iv;
% J_iv = sum((y_used - y_m_iv_used))
