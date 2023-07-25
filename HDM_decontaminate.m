function [signal_decontaminated, stimulus_estimate, noise] = HDM_decontaminate(data, plotFlag, plotTitle)

%% get decontaminated signal
[stimulus_estimate, ~, ~] = HDM_solveInverse(data);
[signal_decontaminated] = HDM_solveForward(stimulus_estimate, 'clean');

%% get noise
% get data structure
p = HDM_getParameters();
dt_data = p.seq.TR;
s = size(data); nT_data = s(end); clear 's';
T_total = nT_data*dt_data;
% match model to data
indice_data = round([dt_data/p.dt:dt_data/p.dt:T_total/p.dt]);  % corresponding indice in model
signal_decontaminated_lo = signal_decontaminated(:, indice_data);
% get noise by substraction
noise = data - signal_decontaminated_lo;

%% plot
if exist('plotFlag', 'var')
    if plotFlag
        t_data = [dt_data:dt_data:T_total];
        t_signal = [p.dt:p.dt:size(signal_decontaminated,2)*p.dt];
        [signal_contaminated] = HDM_solveForward(stimulus_estimate);
        figure;
        for d = 1:p.D
            subplot(p.D,1,d);
            plot(t_data, data(d,:)); hold on;
            plot(t_signal, signal_contaminated(d,:));
            plot(t_signal, signal_decontaminated(d,:));
            if d==1
                if exist('plotTitle', 'var')
                    title(plotTitle);
                end
                legend('data', 'contaminated estimate', 'decontaminated');
            end
        end
    end
end


end

