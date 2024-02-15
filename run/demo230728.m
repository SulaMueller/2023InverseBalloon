clear all; close all; clc;

nT_model = 3000;
p = HDM_getParameters();

%% define stimulus; get clean signal
[stimulus_clean, t0, T] = generateTestStimulus(p, nT_model);
[signal_clean, y] = HDM_solveForward(p, stimulus_clean); 

%% plot and save figures
d = 3;
h = figure; 
plot(stimulus_clean(d,:),'LineWidth',8); title('stimulus');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', 'stimulus', '.png']);
titles = {'n_{excitation}', 'n_{inhibition}', 'vas', 'f_{arteriole}', 'f', '-', 'v', '-', 'q', '-', 'BOLD', 'VASO'};
for f = 1:length(titles)
    if strcmpi(titles{f},'-')
        continue;
    end
    h = figure; 
    plot(y(d+(f-1)*p.D,:),'LineWidth',8); title(titles{f});
    saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', titles{f}, '.png']);
end

%% downsample clean signal to size of data
dt_data = p.seq.TR;
nT_data = floor(nT_model/dt_data*p.dt);
indice_data = [dt_data/p.dt:dt_data/p.dt:nT_data*dt_data/p.dt];  % indice in model
signal_clean_lo = signal_clean(:, indice_data);

%% plot principle of estimation
noise_level = 0.1;
noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
signal_noisy = signal_clean_lo + noise;
h = figure; 
plot(signal_noisy(d,:),'LineWidth',8); title('data');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', 'data', '.png']);
[stimulus_estimate, t0_estimate, T_estimate] = HDM_solveInverse(p, signal_noisy);
h = figure; 
plot(stimulus_estimate(d,:),'LineWidth',8); title('estimate');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', 'estimate', '.png']);

%% plot comparison of noise levels
% start with existing 0.1
[signal_estimate, ~] = HDM_solveForward(p, stimulus_estimate);
t_hi = [p.dt:p.dt:size(signal_estimate,2)*p.dt];
t_lo = [dt_data:dt_data:nT_data*dt_data];
h = figure; 
plot(stimulus_clean(d,:),'LineWidth',8); title(['stimuli, ', num2str(noise_level)]); hold on;
plot(stimulus_estimate(d,:),'LineWidth',8); legend('original', 'estimate');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', ['stimuli, ', num2str(noise_level*100)], '.png']);
h = figure; 
plot(t_lo,signal_noisy(d,:),'LineWidth',8); title(['signals, ', num2str(noise_level)]); hold on;
plot(t_hi,signal_estimate(d,:),'LineWidth',8); legend('original', 'estimate');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', ['signals, ', num2str(noise_level*100)], '.png']);

%% compare more noise levels
noise_levels = [0.05, 0.2];
for i=1:length(noise_levels)
    noise_level = noise_levels(i);
    noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
    signal_noisy = signal_clean_lo + noise;

    [stimulus_estimate, t0_estimate, T_estimate] = HDM_solveInverse(p, signal_noisy);
    [signal_estimate, ~] = HDM_solveForward(p, stimulus_estimate);

    h = figure; 
    plot(stimulus_clean(d,:),'LineWidth',8); title(['stimuli, ', num2str(noise_level)]); hold on;
    plot(stimulus_estimate(d,:),'LineWidth',8); legend('original', 'estimate');
    saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', ['stimuli, ', num2str(noise_level*100)], '.png']);
    h = figure; 
    plot(t_lo,signal_noisy(d,:),'LineWidth',8); title(['signals, ', num2str(noise_level)]); hold on;
    plot(t_hi,signal_estimate(d,:),'LineWidth',8); legend('original', 'estimate');
    saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', ['signals, ', num2str(noise_level*100)], '.png']);
end

%% test decontaminate
noise_level = 0.1;
noise = random('Normal', 0, 1, size(signal_clean_lo)) * ( max(signal_clean_lo(:))-min(signal_clean_lo(:))) * noise_level ;
signal_noisy = signal_clean_lo + noise;
[signal_decontaminated, stimulus_estimate, ~] = HDM_decontaminate(p, signal_noisy);
T_total = nT_data*dt_data;
t_data = [dt_data:dt_data:T_total];
t_signal = [p.dt:p.dt:size(signal_decontaminated,2)*p.dt];
[signal_contaminated] = HDM_solveForward(p, stimulus_estimate);
% plot
h = figure;
plot(t_data, signal_noisy(d,:),'LineWidth',8); hold on;
plot(t_signal, signal_contaminated(d,:),'LineWidth',8);
plot(t_signal, signal_decontaminated(d,:),'LineWidth',8);
title('decontamination signal');
legend('data', 'contaminated estimate', 'decontaminated');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', 'decontamination', '.png']);

p.contaminated = 0;
signal_clean_decontaminated = HDM_solveForward(p, stimulus_clean);
h = figure;
plot(signal_clean_decontaminated(d,:),'LineWidth',8); hold on;
plot(signal_decontaminated(d,:),'LineWidth',8);
title('decontaminated comparison');
legend('original', 'estimate');
saveas(h,['C:\Users\someone\Desktop\work_stuff\BOLD-VASO model\img\', 'contaminated estimate', '.png']);

