function [stimulus_estimate, t0, T] = HDM_solveInverse(p, data)
%% function [stimulus_estimate, t0, T] = HDM_solveInverse(p, data)
% data should be [D, nT]
% note: HDM_solveForward_rewritten4fit calls HDM_getParameters directly

%% get data structure
dt_data = p.seq.TR;
s = size(data); nT_data = s(end); clear 's';  % get number of time points
T_total = nT_data*dt_data;
% turn data into 1D array for fit function
data1D = reshape(data',[nT_data*p.D,1]);
t_axis1D = [dt_data:dt_data:T_total*p.D]';

%% fit stimulus to data
disp('fitting model to data......... ');
t1 = datetime('now');
fo = fitoptions('Method','NonlinearLeastSquares', ...
                'Lower', zeros([1,2*p.D]), ...
                'Upper', T_total * ones([1, 2*p.D]), ...
                'StartPoint', ones([1, 2*p.D]), ...
                'DiffMinChange', p.dt, ...
                'DiffMaxChange', dt_data);
ft = fittype('HDM_solveForward_rewritten4fit(x, T1, T2, T3, T4, T5, T6, t01, t02, t03, t04, t05, t06)', 'options', fo);
f = fit( t_axis1D, data1D, ft );
t2 = datetime('now');
disp(strcat("Fit complete. Required time: ", string(between(t1,t2))));

%% display result
% figure; plot(f, t_axis, data(1,:)');
t0 = [f.t01, f.t02, f.t03, f.t04, f.t05, f.t06]; 
T = [f.T1, f.T2, f.T3, f.T4, f.T5, f.T6];
disp(['Fit result (stimulus): t0 = ', mat2str(t0), ', T = [', num2str(T), ']']);

%% return stimulus
nT_model = T_total/p.dt;
i0 = max(1,    min(nT_model, round(t0/p.dt)));
i1 = max(i0-1, min(nT_model, round(t0/p.dt + T/p.dt)-1));
stimulus_estimate = zeros([p.D, nT_model]);
for d = 1:p.D
    stimulus_estimate(d, i0(d):i1(d)) = 1;
end

end

