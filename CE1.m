Ts = 0.1;

%Step Response
t = (0:Ts:20)';          % 0~20 sec
u = zeros(size(t));
u(t >= 1) = 0.9;         % 0.9 step after 1 sec

simin.time = t;
simin.signals.values = u;
simin.signals.dimensions = 1;

simres = sim("sim1.slx");


    % figure;
    plot(simres.tout, simres.simout.Data);
    grid on;
    xlabel('Time (s)');
    ylabel('Output');
    title('Step response');

%impulse response
impulseStartTime = 1.0;

simin = struct;
simin.time = (0:Ts:20)';
simin.signals.values = (zeros(length(simin.time),1));
simin.signals.dimensions = 1; %%Idk

k = round(impulseStartTime/Ts) + 1;
simin.signals.values(k) = 1/Ts;

simres = sim("sim1.slx");


% figure;
plot(simres.tout, simres.simout.Data);
grid on;
xlabel('Time (s)');
ylabel('Output');
title('Impulse response');

%%1.2
% u = prbs(n,p);

u = prbs(5,3);
u = u(:);

[R,h] = intcor(u,u);

% figure;
stem(h,R);
grid on;
title('PRBS autocorrelation');
xlabel('lag');
ylabel('R');

%% 1.3 Impulse response by deconvolution method

% t = (0:Ts:(N-1)*Ts)';
% u_rand = 1.8*rand(N,1) - 0.9;

%% 1.3 Impulse response by deconvolution method
% random input generation
N = 400;
m = 50;
lambda = 0.1;

t = (0:Ts:(N-1)*Ts)';
u_rand = 1.8 * rand(N,1) - 0.9;

%%check
disp('u_rand size:')
disp(size(u_rand))

% simlink output
simin.time = t;
simin.signals.values = u_rand;
simin.signals.dimensions = 1;

simres = sim("sim1.slx");

%%check
disp('sim done')

y = simres.simout.Data;
y = y(:);

%check
disp('y size:')
disp(size(y))

L = min(length(y), length(u_rand));
y = y(1:L);
u_rand = u_rand(1:L);
t = t(1:L);

% Construct Toeplitz input matrix */
first_col = u_rand;
first_row = [u_rand(1) zeros(1,m-1)];
U = toeplitz(first_col, first_row);

% Least-squares estimate */
g_hat = U \ y;

% Regularized estimate */
g_reg = (U' * U + lambda * eye(m)) \ (U' * y);

% True discrete-time impulse response */
Gs = tf([-1 1.5], [1 0.85 3]);
Gd = c2d(Gs, Ts, 'zoh');

tg = (0:m-1)' * Ts;
[g_true, t_true] = impulse(Gd, tg);
g_true = squeeze(g_true);

% Plots */
figure;
plot(t, u_rand);
grid on;
xlabel('Time (s)');
ylabel('Input');
title('Random input signal');

figure;
plot(t, y);
grid on;
xlabel('Time (s)');
ylabel('Output');
title('Measured output');

figure;
plot(tg, g_hat);
grid on;
xlabel('Time (s)');
ylabel('g');
title('Estimated impulse response by least squares');

figure;
plot(tg, g_reg);
grid on;
xlabel('Time (s)');
ylabel('g');
title('Estimated impulse response by regularization');

figure;
plot(t_true, g_true, 'k-', 'LineWidth', 1.5);
hold on;
plot(tg, g_hat, 'o-');
plot(tg, g_reg, 'x-');
grid on;
xlabel('Time (s)');
ylabel('g');
title('Comparison of impulse responses');
legend('True', 'Least squares', 'Regularized');