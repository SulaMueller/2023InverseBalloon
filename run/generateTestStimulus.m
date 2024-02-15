function [stimulus_clean, t0, T] = generateTestStimulus(p, nT)
% function [stimulus_clean] = generateTestStimulus(p, nT)

stimulus_clean = zeros([p.D, nT]);
t0 = [3.50, 3.60, 3.70, 5.00, 6.00, 1];  % start of S on each depth-level
T = [1.50, 1.50, 1.00, 5.00, 6.00, 0];  % length of S on each depth-level
i0 = round(t0/p.dt);  % indice of S in model
i1 = round((t0+T)/p.dt)-1;
for d = 1:p.D
    stimulus_clean(d, i0(d):i1(d)) = 1;
end

end

