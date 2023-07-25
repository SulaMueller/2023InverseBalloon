clear all; close all; clc;

%% ================== compare Havlicek to mine (for performance) =============
t1 = datetime('now');
Example_dynamic_LBR;
disp(strcat("Havlicek complete. Required time: ", string(between(t1,datetime('now')))));

t1 = datetime('now');
P.N.T     = 30;               % Total lenght of the response (in seconds)
dur       = 2/P.N.dt;         % Stimulus duration (in second, e.g. 2 sec) ... dt - refers to integration step
onset     = 3/P.N.dt;         % Stimulus onset time (in seconds) 
offset    = onset + dur;      % Stimulus offset time (in seconds) 
U.u       = zeros(P.N.T/P.N.dt,K);   % Matrix with input vectors to the neuronal model (one column per depth)
U.u(onset:offset,:) = 1;             % Set one during stimulus window
signal_clean = HDM_solveForward(U.u');
disp(strcat("Forward complete. Required time: ", string(between(t1,datetime('now')))));

HDM_plotD('', signal_clean, LBR','mine', 'H');