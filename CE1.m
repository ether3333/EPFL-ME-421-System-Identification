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

unitStepResponse = simout.simout.Data;

figure;
plot(simout.tout, unitStepResponse)
hold on;
plot(t,u_step)
xlabel('Time (s)'); ylabel('Amplitude');
legend('System Output', 'Step Input');
title('Step response')
hold off

%% Impulse
impulseStartTime = 1.0;
t = (0:Ts:20)';

u_impulse = zeros(length(t),1);
u_impulse(t == 1) = 1/Ts;

simin = struct;
simin.time = (0:Ts:20)';
simin.signals.values = u_impulse;

simout = sim("CE1simulink.slx");

figure;
plot(simout.tout, simout.simout.Data);
hold on;
%plot(t,u_impulse)
xlabel('Time (s)'); ylabel('Amplitude');
%legend('System Output', 'Step Input');
title('Impulse response')
hold off



%% 1.2 Auto Correlation of a PRBS signal
u_prbs = prbs(5,3);
[R_uu, h] = intcor(u_prbs, u_prbs);

figure;
stem(h, R_uu)


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
zero_index = find(lags == 0);
lags = lags .*Ts;
h_yu = h_yu .*Ts;

figure;
plot(h_yu, R_yu_hat)
hold on;
plot(lags((zero_index):(zero_index +600)), R_yu_xcorr((zero_index):(zero_index +600)))
plot(tOut_true, g)
plot(tOut_true, y_true, 'k', 'LineWidth', 2, 'DisplayName', 'True Response');
legend('Intcor', 'Xcor', 'True')
hold off;


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

model = frd(G, frequencies);

%one_period_G = Y_fft_mat(1,:) ./ U_fft_mat(1,:);
%one_period = frd(one_period_G, frequencies);

% Calculate the frequency response from the model
[mag, phase, w] = bode(model);
[mag_true, phase_true, w_true] = bode(G_discrete);

figure;
subplot(2,1,1);
hold on;
semilogx(w, 20*log10(squeeze(mag)));
semilogx(w_true, 20*log10(squeeze(mag_true)));
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');
title('Frequency Response');
grid on;
hold off
