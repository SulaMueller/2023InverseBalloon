clear all; close all; clc;

%% ===============  test inverse on noisy data  ===================

nT_model = 3000;
noise_level = 0.05;

%% define data structure
p = HDM_getParameters();
dt_data = p.seq.TR;
nT_data = floor(nT_model/dt_data*p.dt);
indice_data = [dt_data/p.dt:dt_data/p.dt:nT_data*dt_data/p.dt];  % indice in model

%% define stimulus; get clean signal
[stimulus_clean, t0, T] = generateTestStimulus(p, nT_model);
signal_clean = HDM_solveForward(p, stimulus_clean);
signal_clean_lo = signal_clean(:, indice_data);

%% create noisy data
noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
disp(['signal = [', num2str(max(signal_clean_lo(:))), ', ', num2str(min(signal_clean_lo(:))), ']']);
disp(['noise = [', num2str(max(noise(:))), ', ', num2str(min(noise(:))), ']']);
signal_noisy = signal_clean_lo + noise;

%% fit model to data -> get stimulus
[stimulus_estimate, t0_estimate, T_estimate] = HDM_solveInverse(p, signal_noisy);

%% plot result
signal_estimate = HDM_solveForward(p, stimulus_estimate);

figure;
t_hi = [p.dt:p.dt:size(signal_estimate,2)*p.dt];
t_lo = [dt_data:dt_data:nT_data*dt_data];
for d = 1:p.D
    subplot(p.D,1,d);
    plot(t_hi, signal_clean(d,1:size(signal_estimate,2))); hold on;
    plot(t_lo, signal_noisy(d,:));
    plot(t_hi, signal_estimate(d,:));
    if d==1
        legend('clean', 'noisy', 'estimate');
    end
end

% plot coefficients
figure;
subplot(2,1,1); plot(t0, 'o'); hold on; plot(t0_estimate, 'x'); legend('original', 'estimate'); title('t0');
subplot(2,1,2); plot(T , 'o'); hold on; plot(T_estimate , 'x'); legend('original', 'estimate'); title('T');

