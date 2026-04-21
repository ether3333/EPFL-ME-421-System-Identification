Ts = 0.1;
A = 0.9;  % Step Magnitude

%% Step function
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

%% Impulse
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

G = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G, Ts, 'zoh');
y_imp_true = impulse(G_discrete, t_imp) * Ts;

figure;
plot(t_imp, y_imp_true, 'LineWidth', 1.5)
hold on;
plot(t_imp, y_imp_noisy, '--', 'LineWidth', 1.2)
xlabel('Time (s)');
ylabel('Amplitude');
legend('Noiseless impulse response', 'Noisy impulse response');
title('Impulse Response: Noiseless vs Noisy');
grid on;
hold off



%% 1.2 Auto Correlation of a PRBS signal
u_prbs = prbs(6,4); %can't remember why you did (5,3) but I changed it as what report requires
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

u_rand = 1.8*rand(N,1)-0.9;

input_mat = toeplitz(u_rand, [u_rand(1), zeros(1, N-1)]);

t = 0:Ts:(N-1)*Ts;

simin = struct;
simin.time = t;
simin.signals.values = u_rand;

simout = sim("CE1simulink.slx");
Y = simout.simout.Data(1:N,:);

% Finite impulse response
K = 20*1/Ts; %find index for 20 seconds where the response dies out
input_mat_trunc = input_mat(:, 1:K);
u_K = u_rand(1:K,1);

Theta_K = (input_mat_trunc \ Y);

t_ir = (0:K-1)' * Ts;


% Regularisation
lamda = 8;
Theta_reg = ((input_mat' * input_mat + lamda * eye(N)) \ (input_mat' * Y));

G = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G, Ts, 'zoh');

[y_true, tOut_true] = impulse(G_discrete, t_ir); %gives g(t) - multiply by Ts
y_true = y_true .* Ts;

% Comparison Plot
figure;
plot(t_ir, Theta_K, 'o-', 'DisplayName', 'Deconvolution');
hold on;
plot(t_ir, Theta_reg(1:K), 'x-', 'DisplayName', 'Regularization');
plot(tOut_true, y_true, 'k', 'LineWidth', 2, 'DisplayName', 'True Response');
xlabel('Time (s)');
ylabel('Amplitude');
legend;
title('Comparison of Impulse Responses');
grid on;


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

g = R_uu_matrix_trunc \ R_yu_hat_trunc; %%ASK IF WE NEED TO DO THIS

[R_yu_xcorr, lags] = xcorr(Y, Uprbs, 'unbiased');
[R_uu_xcorr, lags_uu] = xcorr(Uprbs, Uprbs, 'unbiased');
zero_index_uu = find(lags_uu == 0);
zero_index_yu = find(lags == 0);

R_uu_xcorr_trunc = R_uu_xcorr(zero_index_uu:zero_index_uu+K-1);
R_yu_xcorr_trunc = R_yu_xcorr(zero_index_yu:zero_index_yu+K-1);

R_uu_xcorr_matrix = toeplitz(R_uu_xcorr_trunc);
g_xcorr = R_uu_xcorr_matrix \ R_yu_xcorr_trunc;

t_corr = (0:K-1)' * Ts;

err_intcor = norm(g_intcor - y_true, 2);
err_xcorr = norm(g_xcorr - y_true, 2);

figure;
plot(t_corr, g_intcor, 'o-', 'DisplayName', 'Impulse response using intcor')
hold on;
plot(t_corr, g_xcorr, 'x-', 'DisplayName', 'Impulse response using xcorr')
plot(tOut_true, y_true, 'k', 'LineWidth', 2, 'DisplayName', 'True response')
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

Uprbs = prbs(8,8);

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


omega_s = 2*pi*(1/Ts);
frequencies = 0:omega_s/length_of_period:(length_of_period-1)*omega_s/length_of_period;

% Number of points to keep (0 to Nyquist)
n_half = floor(length_of_period/2) + 1;

G_half = G(1:n_half);
freq_half = frequencies(1:n_half);

% Create the model using only the meaningful half
model = frd(G_half, freq_half);

%model = frd(G, frequencies);

one_period_G = Y_fft_mat(8,:) ./ U_fft_mat(8,:);
one_period = frd(one_period_G, frequencies);

% Calculate the frequency response from the model
[mag, phase, w] = bode(model);

G = tf([-1 1.5], [1 0.85 3]);
G_discrete = c2d(G, Ts, 'zoh');

[mag_true, phase_true, w_true] = bode(G_discrete);

figure;
hold on;
% semilogx(w, 20*log10(squeeze(mag)));
% semilogx(w_true, 20*log10(squeeze(mag_true)));
bode(model,G_discrete)
% xlabel('Frequency (rad/s)');
% ylabel('Magnitude (dB)');
title('Frequency Response');
grid on;
hold off

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
