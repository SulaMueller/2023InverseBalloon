
clear all; close all; clc;

%% ===============  test decontamination on noisy data  ===================

nT_model = 3000;
noise_level = 0.05;

%% define data structure
p = HDM_getParameters();
dt_data = p.seq.TR;
nT_data = floor(nT_model/dt_data*p.dt);
indice_data = [dt_data/p.dt:dt_data/p.dt:nT_data*dt_data/p.dt];  % indice in model

%% define stimulus; get clean signal
stimulus_clean = zeros([p.D,nT_model]);
t0 = [3.50, 3.60, 3.70, 5.00, 6.00, 1];
T = [1.50, 1.50, 1.00, 5.00, 6.00, 0];
i0 = round(t0/p.dt);
i1 = round((t0+T)/p.dt)-1;
for d = 1:p.D
    stimulus_clean(d, i0(d):i1(d)) = 1;
end

signal_clean = HDM_solveForward(stimulus_clean);
signal_clean_lo = signal_clean(:, indice_data);

%% create noisy data
noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
disp(['signal = [', num2str(max(signal_clean_lo(:))), ', ', num2str(min(signal_clean_lo(:))), ']']);
disp(['noise = [', num2str(max(noise(:))), ', ', num2str(min(noise(:))), ']']);
signal_noisy = signal_clean_lo + noise;

%% test decontaminate
[signal_decontaminated, stimulus_estimate, noise_estimate] = HDM_decontaminate(signal_noisy, 1, 'decontamination test');

HDM_plotD('stimuli', stimulus_clean, stimulus_estimate, 'original', 'estimate');
signal_clean_decontaminated = HDM_solveForward(stimulus_clean, 'clean');
HDM_plotD('decontaminated signals', signal_clean_decontaminated, signal_decontaminated, 'original', 'estimate');


