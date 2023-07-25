function [p] = HDM_getParameters(isContaminated)

if ~exist('isContaminated', 'var')
    isContaminated = 1;
end
p.contaminated = isContaminated;

%% basic parameters
p.D = 6;  % number of depth levels
p.K = 3;  % number of compartments [arteriole, venule, vein]
p.dt = 0.01;  % [s], dt_model
p.ARTERIOLE = 1;
p.VENULE = 2;
p.VEIN = 3;

%% sequence parameters
p.seq.TR = 1.4;  % [s]
p.seq.TE = 0.028;  % [s]

%% BOLD parameters (physical and physiological)
p.bold.B0 = 7;  % [T]
p.bold.dXi = 0.000000264;  % susceptibility difference between oxy/deoxy blood
p.bold.E0 = [0.35, 0.35, 0.35];  % oxygen extraction fraction [arteriole, venule, vein]
p.bold.eps = [0, 0.2706, 0.2334];  % intra/extra-vascular signal [arteriole, venule, vein]
p.bold.Hct = [0, 0.35, 0.38];  % hematocrit [arteriole, venule, vein]
p.bold.r0 = [0, 228, 232];  % signal relaxation/oxygen saturation [arteriole, venule, vein]
p.bold.gamma0 = 2*pi*42.58*power(10,6);

%% balloon parameters
p.balloon.F0 = [0.2083, 0.2083, 0.2083, 0.2083, 0.2083, 0.2083;  % blood flow in resting conditions
                0.2083, 0.2083, 0.2083, 0.2083, 0.2083, 0.2083;  % all depths; [arteriole, venule, vein]
                1.2498, 1.0415, 0.8332, 0.6249, 0.4166, 0.2083];  % copied from Python implementation
            
p.balloon.V0 = [0     , 0     , 0     , 0     , 0     , 0     ;  % blood volume in resting conditions
                0.2083, 0.2083, 0.2083, 0.2083, 0.2083, 0.2083;  % all depths; [arteriole, venule, vein]
                0.3125, 0.2708, 0.2292, 0.1875, 0.1458, 0.1042];  % copied from Python implementation
            
p.balloon.tau0 = [0.        , 0.        , 0.        , 0.        , 0.        ,0.        ;  % time constant in resting conditions
                  1.        , 1.        , 1.        , 1.        , 1.        ,1.        ;  % all depths; [arteriole, venule, vein]
                  0.25004001, 0.2600096 , 0.27508401, 0.30004801, 0.349976  ,0.50024004];  % copied from Python implementation

p.balloon.alpha = [0, 0.35, 0.2];  % flow/volume coupling [arteriole, venule, vein]
p.balloon.vet = [0, 2, 2; 
                 0, 2, 30];  % visco-elastic time constant [arteriole, venule, vein] [in; out]
p.balloon.n = [0, 4, 4];  % n-ratio (blood flow/oxygen consumption) [arteriole, venule, vein]

%% neural parameters [Havlicek2015]
p.neuro.C = 1;  % response weighting of stimulus
p.neuro.sigma = -3;  % excitatory self-connection 
p.neuro.mu = 1.5;  % inhibition of excitatory activation
p.neuro.lambda = 0.2;  % inhibitory self-connection
p.neuro.theta = 0.6;  % decay of vaso-active signal
p.neuro.phi = 1.5;  % gain of vaso-active signal
p.neuro.xi = 0.6;  % decay of blood inflow signal

end

