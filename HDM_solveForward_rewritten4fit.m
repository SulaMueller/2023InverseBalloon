function [signal_lo_1D] = HDM_solveForward_rewritten4fit(x, T1, T2, T3, T4, T5, T6, t01, t02, t03, t04, t05, t06)
%% function [signal_lo_1D] = HDM_solveForward_rewritten4fit(x, T1, T2, T3, T4, T5, T6, t01, t02, t03, t04, t05, t06)
% wrapper for HDM_solveForward to match fitting procedure
% result gives info about stimulus:
%       * most of stimulus is all zeros
%       * active sections are defined by start point t and length T
% x: input created by fitter (array of a few increasing monospaced input values; corresponds to
% time axis; eg x = [0 1 2 3 4 5 6 7]

p = HDM_getParameters();

%% crop x to time axis
% (output only accepts same length as x)
% todo: probably not the best fitting procedure
nPerDepth = ceil(size(x,1)/p.D);  % number of time points: at least the first layer gets nPD points
t = x(1:nPerDepth);  % time axis
T_total = t(end);
i_model = ceil(t/p.dt);  % indice of all time points in model

%% define stimulus
nT_model = T_total/p.dt;
S = zeros([p.D, nT_model]);
nT_active = round([T1, T2, T3, T4, T5, T6]/p.dt);
i0_active = round([t01, t02, t03, t04, t05, t06]/p.dt);
for d = 1:p.D
    i_start = max(1,       min(nT_model, i0_active(d)));
    i_end   = max(i_start, min(nT_model, i0_active(d)+nT_active(d)-1));
    if i_start == i_end
        continue;
    end
    S(d, i_start:i_end) = 1;
end

%% get signal for all time points in x
signal = HDM_solveForward(p, S);
signal_lo_1D = reshape(signal(:, i_model)', [length(i_model)*p.D, 1]);
signal_lo_1D = signal_lo_1D(1:length(x));

end