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
[stimulus_clean, t0, T] = generateTestStimulus(p, nT_model);
signal_clean = HDM_solveForward(p, stimulus_clean);
signal_clean_lo = signal_clean(:, indice_data);

%% create noisy data
noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
disp(['signal = [', num2str(max(signal_clean_lo(:))), ', ', num2str(min(signal_clean_lo(:))), ']']);
disp(['noise = [', num2str(max(noise(:))), ', ', num2str(min(noise(:))), ']']);
signal_noisy = signal_clean_lo + noise;

%% test decontaminate
[signal_decontaminated, stimulus_estimate, noise_estimate] = HDM_decontaminate(p, signal_noisy, 1, 'decontamination demo');

%% plot difference to clean test data
HDM_plotD(p.D, 'stimuli', stimulus_clean, stimulus_estimate, 'original', 'estimate');
p.contaminated = 0;
signal_clean_decontaminated = HDM_solveForward(p, stimulus_clean);
HDM_plotD(p.D, 'decontaminated signals', signal_clean_decontaminated, signal_decontaminated, 'original', 'estimate');



