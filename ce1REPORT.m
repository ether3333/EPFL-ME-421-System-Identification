Ts = 0.1;
A = 0.9;  % Step Magnitude

s = tf('s');
Gs = (1.5 - s) / (s^2 + 0.85*s + 3);
bode(Gs);

%%
s = tf('s');
Gs = (1.5 - s) / (s^2 + 0.85*s + 3);

% Simulate with white noise input
Ts = 0.1;
t = 0:Ts:100;
u = randn(size(t));          % white noise input
y = lsim(Gs, u, t);          % simulate output

% FFT of output
N = length(y);
Y = fft(y);
freq = (0:N/2-1) / (N*Ts);           % frequency axis in Hz
magnitude = abs(Y(1:N/2)) / N;       % single-sided spectrum
omega = 2*pi * freq;   % rad/s
% Plot
plot(omega, magnitude);
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
title('Output spectrum');
grid on;

%% 

G = tf([-1 1.5], [1 0.85 3]);
bode(G)

[~, wn]   = damp(G);          % natural frequencies
wb        = bandwidth(G);      % -3 dB bandwidth (rad/s)
ws_used   = 2*pi / Ts;          % sampling frequency (rad/s)
 
fprintf('\n--- Sampling-time justification ---\n');
fprintf('Bandwidth  wb  = %.4f rad/s\n', wb);
fprintf('Sampling freq ws = 2pi/Ts = %.4f rad/s\n', ws_used);
fprintf('ws / wb = %.2f  (should be 20..30 for good ID)\n', ws_used / wb);
fprintf('Ts chosen  = %.2f s\n', Ts);
fprintf('Ts*wb      = %.4f  (should be in [2pi/30, 2pi/20] ≈ [0.21, 0.31])\n', Ts*wb);


%% 1.1 Step response
stepStartTime = 1.0;
t = (0:Ts:20)';

u_step = zeros(length(t),1);
u_step(t>=1) = A;

simin = struct;
simin.time = t;
simin.signals.values = u_step;

simout = sim("CE1simulink.slx");

y_step_noisy = simout.simout.Data;

G = tf([-1 1.5], [1 0.85 3]);
t_step = simout.tout;
y_step_true = step(G, t_step);



figure;
plot(t_step, y_step_true, 'LineWidth', 1.5)
hold on;
plot(t_step, y_step_noisy, '--', 'LineWidth', 1.2)
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noiseless step response', 'Noisy step response');
title('Unit Step Response: Noiseless vs Noisy');
grid on;
hold off


%% 1.1 Step response CORRECT
stepStartTime = 1.0;
t = (0:Ts:20)';
u_step = zeros(length(t),1);
u_step(t>=1) = A;
simin.time = t;
simin.signals.values = u_step;

% Run NOISY simulation and save everything immediately
simout_noisy = sim('CE1simulink');
y_step_noisy = simout_noisy.simout.Data / A;  % normalize to unit step
t_step = simout_noisy.tout;

% Set noise to 0 and run NOISELESS simulation
set_param('CE1simulink/Random Number', 'Variance', '0');
simout_noiseless = sim('CE1simulink');
y_step_true = simout_noiseless.simout.Data / A;  % normalize to unit step

% Restore noise
set_param('CE1simulink/Random Number', 'Variance', '0.015');

figure;
hold on;
plot(t_step, y_step_true, 'LineWidth', 1.5)
plot(t_step, y_step_noisy, 'LineWidth', 1.2)
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noiseless step response', 'Noisy step response');
title('Unit Step Response: Noiseless vs Noisy');
grid on;

%% 1.1  Impulse response
impulseStartTime = 1.0;
t = (0:Ts:20)';

u_impulse = zeros(length(t),1);
[~, idx_imp] = min(abs(t - impulseStartTime));
u_impulse(idx_imp) = 1/Ts;

simin = struct;
simin.time = (0:Ts:20)';
simin.signals.values = u_impulse;

simout = sim("CE1simulink.slx");

y_imp_noisy = simout.simout.Data;
t_imp = simout.tout;

% G = tf([-1 1.5], [1 0.85 3]);
% G_discrete = c2d(G, Ts, 'zoh');
% y_imp_true = impulse(G_discrete, t_imp) * Ts;


% Set noise variance to 0
set_param('CE1simulink/Random Number', 'Variance', '0');

% Run noiseless simulation
sim('CE1simulink');
simout_noiseless = simout;
y_imp_true = simout_noiseless.simout.Data;

% Restore original variance
set_param('CE1simulink/Random Number', 'Variance', '0.015');

figure;

hold on;


plot(t_imp, y_imp_true, 'LineWidth', 1.5)
plot(t_imp, y_imp_noisy, '--', 'LineWidth', 1.2)
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noisy impulse response', 'Noiseless impulse response');
title('Impulse Response: Noiseless vs Noisy');
grid on;
hold off

%% 1.1 Impulse response CORRECT
impulseStartTime = 1.0;
t = (0:Ts:20)';
u_impulse = zeros(length(t),1);
[~, idx_imp] = min(abs(t - impulseStartTime));
u_impulse(idx_imp) = 1/Ts;
simin.time = t;
simin.signals.values = u_impulse;

% Run NOISY simulation and save everything immediately
simout_noisy = sim('CE1simulink');
y_imp_noisy = simout_noisy.simout.Data;
t_imp = simout_noisy.tout;

% Set noise to 0 and run NOISELESS simulation
set_param('CE1simulink/Random Number', 'Variance', '0');
simout_noiseless = sim('CE1simulink');
y_imp_true = simout_noiseless.simout.Data;

% Restore noise
set_param('CE1simulink/Random Number', 'Variance', '0.015');

figure;
hold on;
plot(t_imp, y_imp_true, 'LineWidth', 1.5)
plot(t_imp, y_imp_noisy, '--', 'LineWidth', 1.2)
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noiseless impulse response', 'Noisy impulse response');
title('Impulse Response: Noiseless vs Noisy');
grid on;


%% 1.2 Auto Correlation of a PRBS signal
u_prbs = prbs(6,4); %not sure why you did (5,3) but I changed it as what GPT said
[R_uu, h] = intcor(u_prbs, u_prbs);

figure;
stem(h, R_uu, 'filled')
xlabel('Lag')
ylabel('Autocorrelation')
title('Autocorrelation of PRBS(6,4)')
grid on


%% 1.3 Impulse response by deconvolution method
Ts = 0.1;
Tf = 50;
N = (Tf/Ts);
t = 0:Ts:(N-1)*Ts;

u_rand = 1.8*rand(N,1)-0.9;
input_mat = toeplitz(u_rand, [u_rand(1), zeros(1, N-1)]);

simin = struct;
simin.time = t;
simin.signals.values = u_rand;

simout = sim("CE1simulink.slx");
Y = simout.simout.Data(1:N,:);

% Finite impulse response
K = 20/Ts; %find index for 20 seconds where the response dies out
input_mat_trunc = input_mat(:, 1:K);
u_K = u_rand(1:K,1);

Theta_fir = (input_mat_trunc \ Y);
t_fir = (0:K-1)' * Ts;


% Regularisation
lamda = 40;
Theta_reg = ((input_mat' * input_mat + lamda * eye(N)) \ (input_mat' * Y));

% True impulse response
G = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G, Ts, 'zoh');

[y_true, tOut_true] = impulse(G_discrete, t); %gives g(t) - multiply by Ts
y_true = y_true .* Ts;

% Comparison Plot
figure;
plot(t_fir, Theta_fir, 'o-', 'DisplayName', 'Deconvolution');
hold on;
plot(t, Theta_reg, 'x-', 'DisplayName', 'Regularization');
plot(tOut_true, y_true, 'k', 'LineWidth', 2, 'DisplayName', 'True Response');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
title('Comparison of Impulse Responses');
grid on;

% Claculate error
err_finite = norm(Theta_fir - y_true, 2);
err_reg = norm(Theta_reg - y_true, 2);

fprintf('2-norm error using finite impulse response: %.6f\n', err_finite);
fprintf('2-norm error using regularisation: %.6f\n', err_reg);

%% 1.4 Impulse response by correlation approach
Uprbs = prbs(8,3);
N = length(Uprbs);

t = 0:Ts:(N-1)*Ts;

simin = struct;
simin.time = t;
simin.signals.values = Uprbs;

simout = sim("CE1simulink.slx");
Y = simout.simout.Data(1:N,:);

[R_uu_hat, h_uu] = intcor(Uprbs, Uprbs);
[R_yu_hat, h_yu] = intcor(Y, Uprbs);

R_uu_hat_trunc = R_uu_hat(1:K, :);
R_yu_hat_trunc = R_yu_hat(1:K,:);

R_uu_matrix_trunc = toeplitz(R_uu_hat_trunc);

g_intcor = R_uu_matrix_trunc \ R_yu_hat_trunc; 

[R_yu_xcorr, lags] = xcorr(Y, Uprbs, 'unbiased');
[R_uu_xcorr, lags_uu] = xcorr(Uprbs, Uprbs, 'unbiased');
zero_index_uu = find(lags_uu == 0);
zero_index_yu = find(lags == 0);

R_uu_xcorr_trunc = R_uu_xcorr(zero_index_uu:zero_index_uu+K-1);
R_yu_xcorr_trunc = R_yu_xcorr(zero_index_yu:zero_index_yu+K-1);

R_uu_xcorr_matrix = toeplitz(R_uu_xcorr_trunc);
g_xcorr = R_uu_xcorr_matrix \ R_yu_xcorr_trunc;

t_corr = (0:K-1)' * Ts;

err_intcor = norm(g_intcor - y_true_K, 2);
err_xcorr = norm(g_xcorr - y_true_K, 2);

figure;
plot(t_corr, y_true_K, 'k', 'LineWidth', 2, 'DisplayName', 'True response')
hold on;
plot(t_corr, g_intcor, 'o','LineWidth', 1, 'DisplayName', 'Impulse response using intcor')
plot(t_corr, g_xcorr,'bx','LineWidth', 1, 'DisplayName', 'Impulse response using xcorr')

xlabel('Time (s)')
ylabel('Amplitude')
legend
title('Impulse Responses Identified by Correlation Methods')
grid on
hold off;

fprintf('2-norm error using intcor: %.6f\n', err_intcor);
fprintf('2-norm error using xcorr : %.6f\n', err_xcorr);


%% 1.5 Frequency Domain Identification (Periodic signal)
Ts = 0.1;
Uprbs = prbs(8,8); %gives a length of 2040
N = length(Uprbs);
t = 0:Ts:(N-1)*Ts;

simin = struct;
simin.time = t;
simin.signals.values = Uprbs;

simout = sim("CE1simulink.slx");
Y = simout.simout.Data(1:N,:);

periods = 8;
length_of_period = N/periods;
U_fft_mat = zeros(periods, length_of_period);
Y_fft_mat = zeros(periods, length_of_period);

for p = 1:periods
    U_fft_mat(p, :) = fft(Uprbs((p-1)*length_of_period + 1:p*length_of_period));
    Y_fft_mat(p, :) = fft(Y((p-1)*length_of_period + 1:p*length_of_period));
end

U_fft_avg = mean(U_fft_mat,1);
Y_fft_avg = mean(Y_fft_mat, 1);

G = Y_fft_avg ./ U_fft_avg;

omega_s = 2*pi/Ts;
frequencies = 0:omega_s/length_of_period:(length_of_period-1)*omega_s/length_of_period;

% Number of points to keep (0 to Nyquist frequency)
n_half = floor(length_of_period/2) + 1;

G_half = G(1:n_half);
freq_half = frequencies(1:n_half);

model = frd(G_half, freq_half);

G_true = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G_true, Ts, 'zoh');

figure;
bode(model,G_discrete, freq_half)
legend('Averaged (8 periods)', 'True response');
title('Frequency Response Identification via PRBS');
grid on;

%% 1.6 Frequency domain Identification (Random signal)
Ts = 0.1;
fs = 1/Ts;
A = 0.9;
N = 2000;

u_rand = A*sign(randn(1, N));

t = 0:Ts:(N-1)*Ts;

simin = struct;
simin.time = t;
simin.signals.values = u_rand';

simout = sim("CE1simulink.slx");
y = simout.simout.Data(1:N,:);

U = fft(u);
Y = fft(y);

Suu = (abs(U).^2) / (N * fs);               % auto PSD of input
Syu = (Y .* conj(U)) / (N * fs);            % cross PSD output/input

H = Syu ./ Suu;                             % FRF estimate

f = (0:N-1) * (fs/N);
f_half = f(1:N/2+1);
H_half = H(1:N/2+1);

figure;
subplot(2,1,1)
plot(f_half, 20*log10(abs(H_half)))
ylabel('Magnitude (dB)')
title('Frequency Response Function')
grid on

subplot(2,1,2)
plot(f_half, angle(H_half) * 180/pi)
ylabel('Phase (degrees)')
xlabel('Frequency (Hz)')
grid on


%% 1.6 Frequency domain Identification (Random signal) - Spectral analysis
Ts = 0.1;
A  = 0.9;
N  = 2000;
 
u_rand = A*sign(randn(N,1)); % Binary random signal (+/-A)
 
t = (0:Ts:(N-1)*Ts)';
 
simin = struct;
simin.time = t;
simin.signals.values = u_rand;
 
simout = sim("CE1simulink.slx");
y = simout.simout.Data(1:N,:);
 
% Basic spectral analysis: G = Phi_yu / Phi_uu
[R_uu_full, lags_uu] = xcorr(u_rand, u_rand, 'biased');
[R_yu_full, lags_yu] = xcorr(y, u_rand, 'biased');
 
R_uu_pos = R_uu_full(lags_uu >= 0);% Keep only the non-negative lags
R_yu_pos = R_yu_full(lags_yu >= 0); 

 
Phi_uu = fft(R_uu_pos);
Phi_yu = fft(R_yu_pos);
G_basic = Phi_yu ./ Phi_uu;
 
% Windowed 
M = 200; % number of lags in window 
w = hann(2*M+1);
w = w(M+1:end); % only the positive half
f_win = zeros(length(R_uu_pos),1);
f_win(1:M+1) = w; % zero-padded to full length
 
Phi_uu_w = fft(R_uu_pos .* f_win);
Phi_yu_w = fft(R_yu_pos .* f_win);
G_window = Phi_yu_w ./ Phi_uu_w;
 
% Averaging
m = 6;
Nseg = floor(N/m);
Phi_uu_avg = zeros(Nseg,1);
Phi_yu_avg = zeros(Nseg,1);
 
for i = 1:m
    idx = (i-1)*Nseg + 1 : i*Nseg;
    u_i = u_rand(idx);
    y_i = y(idx);
 
    [R_uu_i, lags_i] = xcorr(u_i, u_i, 'biased');
    [R_yu_i, ~] = xcorr(y_i, u_i, 'biased');
 
    R_uu_i_pos = R_uu_i(lags_i >= 0);
    R_yu_i_pos = R_yu_i(lags_i >= 0);
 
    Phi_uu_avg = Phi_uu_avg + fft(R_uu_i_pos);
    Phi_yu_avg = Phi_yu_avg + fft(R_yu_i_pos);
end
Phi_uu_avg = Phi_uu_avg / m;
Phi_yu_avg = Phi_yu_avg / m;
G_avg = Phi_yu_avg ./ Phi_uu_avg;
 
% Frequency vectors
omega_s = 2*pi/Ts;
w_full  = (0:N-1)'    * omega_s/N;
w_seg   = (0:Nseg-1)' * omega_s/Nseg;
 
n_half    = floor(N/2)    + 1;
n_half_s  = floor(Nseg/2) + 1;
 
model_basic  = frd(G_basic(1:n_half),    w_full(1:n_half));
model_window = frd(G_window(1:n_half),   w_full(1:n_half));
model_avg    = frd(G_avg(1:n_half_s),    w_seg(1:n_half_s));
 
% True discrete-time model for comparison
G_true     = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G_true, Ts, 'zoh');
 

figure;
bode(model_basic, G_discrete)
legend('Basic spectral analysis', 'True response');
title('Basic spectral analysis');
grid on;
 
figure;
bode(model_window, G_discrete)
legend('Hann window', 'True response');
title('Spectral analysis with Hann window');
grid on;
 
figure;
bode(model_avg, G_discrete)
legend(sprintf('Averaging (m = %d)', m), 'True response');
title('Spectral analysis with averaging');
grid on;