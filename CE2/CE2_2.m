%% 2.2 Parametric Identification of an Active Suspension System
load("ASSdata.mat")
data = iddata(y,u,Ts);
data = detrend(data);
%% 2.2.1 Order estimation
nk = 1;

orders = zeros(1,10);
losses = zeros(1,10);

for delta = 1:10
    na = delta;
    nb = delta;
    model = arx(data, [na nb nk]);
    orders(delta) = delta;
    losses(delta)= model.EstimationInfo.LossFcn;
end  

figure;
plot(orders, losses)
xlabel('Order of the Model');
ylabel('Loss Function Value');
title('Order Estimation for Active Suspension System');
grid on;


% Zero/pole cancellation
for delta = 2:6
    na = delta;
    nb = delta;
    nc = delta;
    model = armax(data, [na nb nc nk]);
    figure;
    h = iopzplot(model); %4 poles, cannot cancel the right ones that are mirrored
    showConfidence(h, 2);

end  

% Estimate delay
model = arx(data, [8 8 1]); %use estimated order, delay = 1 (first value is 0), the others dont matter
disp('B coefficients:'); disp(model.B)
disp('Standard deviations of B:'); disp(model.db)


NN = struc(1:10, 1:10, 1:5);  
V = arxstruc(data(1:500), data(501:end), NN);
nn = selstruc(V);  

% NN = struc(1:3,1:2,2:4);
% V = arxstruc(data(1:500), data(501:1000), NN);
% nn = selstruc(V)


% bode plot for order estimation - resonance peaks (1 - 2/3 order), then
% zero-pole etc was what the TA usually woukd do (start with Bode)









%% 2.2 Parametric Identification of an Active Suspension System
load("ASSdata.mat")
data = iddata(y,u,Ts);
data = detrend(data);

%% 2.2.1 Order estimation
nk = 1;
orders = zeros(1,10);
losses = zeros(1,10);
for delta = 1:10
    na = delta;
    nb = delta;
    model = arx(data, [na nb nk]);
    orders(delta) = delta;
    losses(delta) = model.EstimationInfo.LossFcn;
end
figure;
plot(orders, losses, '-o')
xlabel('Order of the Model');
ylabel('Loss Function Value');
title('Order Estimation for Active Suspension System');
grid on;

% Zero/pole cancellation
for delta = 2:6
    na = delta;
    nb = delta;
    nc = delta;
    model = armax(data, [na nb nc nk]);
    figure;
    h = iopzplot(model); % 4 poles, cannot cancel the right ones that are mirrored
    showConfidence(h, 2);
    title(['ARMAX zero/pole plot, order = ', num2str(delta)]);
end

% Estimate delay
model = arx(data, [8 8 1]); % use estimated order, delay = 1 (first value is 0)
disp('B coefficients:'); disp(model.B)
disp('Standard deviations of B:'); disp(model.db)

% Compare with struc/arxstruc/selstruc
NN = struc(1:10, 1:10, 1:5);
V = arxstruc(data(1:500), data(501:end), NN);
nn = selstruc(V);
disp('Order suggested by selstruc [na nb nk]:'); disp(nn)

%% 2.2.2 Parametric identification
% Choose the estimated orders from 2.2.1
na = 4;     % adjust based on your order-estimation conclusion
nb = 4;
nk = 1;
nc = na;
nd = na;
nf = na;
nx = na;    % number of states for n4sid (global order = delta)

% Split data: first half for identification, second half for validation
N = length(data.y);
Nid = floor(N/2);
data_id  = data(1:Nid);
data_val = data(Nid+1:end);

% Identify models with different structures
m_arx   = arx(data_id,   [na nb nk]);
m_iv4   = iv4(data_id,   [na nb nk]);
m_armax = armax(data_id, [na nb nc nk]);
m_oe    = oe(data_id,    [nb nf nk]);
m_bj    = bj(data_id,    [nb nc nd nf nk]);
m_n4sid = n4sid(data_id, nx);

%% 2.2.3 Model validation
% 1) Compare simulated output to measured output on the validation set
figure;
compare(data_val, m_arx, m_iv4, m_armax, m_oe, m_bj, m_n4sid);
title('Measured vs. simulated output (validation data)');

% 2) Compare frequency response with a nonparametric spectral model
g_spa = spa(data_val);   % nonparametric (spectral analysis) reference

figure;
bode(m_arx, m_iv4, m_armax, m_oe, m_bj, m_n4sid, g_spa);
legend('ARX','IV4','ARMAX','OE','BJ','N4SID','SPA');
title('Frequency response: parametric models vs. nonparametric (SPA)');
grid on;

% 3) Residual analysis: whiteness of residuals + cross-correlation with past inputs
figure; resid(data_val, m_arx);   title('Residuals - ARX');
figure; resid(data_val, m_iv4);   title('Residuals - IV4');
figure; resid(data_val, m_armax); title('Residuals - ARMAX');
figure; resid(data_val, m_oe);    title('Residuals - OE');
figure; resid(data_val, m_bj);    title('Residuals - BJ');
figure; resid(data_val, m_n4sid); title('Residuals - N4SID');

% Print fit percentages for quick comparison
[~, fit_arx]   = compare(data_val, m_arx);
[~, fit_iv4]   = compare(data_val, m_iv4);
[~, fit_armax] = compare(data_val, m_armax);
[~, fit_oe]    = compare(data_val, m_oe);
[~, fit_bj]    = compare(data_val, m_bj);
[~, fit_n4sid] = compare(data_val, m_n4sid);

fprintf('\nFit on validation data (%%):\n');
fprintf('  ARX  : %6.2f\n', fit_arx);
fprintf('  IV4  : %6.2f\n', fit_iv4);
fprintf('  ARMAX: %6.2f\n', fit_armax);
fprintf('  OE   : %6.2f\n', fit_oe);
fprintf('  BJ   : %6.2f\n', fit_bj);
fprintf('  N4SID: %6.2f\n', fit_n4sid);