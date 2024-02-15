clear all; close all; clc;

%% ===============  test generation of forward data  ===================

%% define data structure
nT = 3000;
p = HDM_getParameters();

%% define stimulus 
S = generateTestStimulus(p, nT);

%% get contaminated and clean signal
[signal_contaminated, y_contaminated] = HDM_solveForward(p, S);
p.contaminated = 0;
[signal_clean, y_clean] = HDM_solveForward(p, S);

%% plot result
HDM_plotY(p.D, y_contaminated);
HDM_plotY(p.D, y_clean);


