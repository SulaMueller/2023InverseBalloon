function [signal] = HDM_solveForward(S, modelType, y0)

if ~exist('modelType', 'var')
    modelType = 'contaminated';
end

% prepare
p = HDM_getParameters();
s = size(S); p.nT = s(end); clear 's';
p.contaminated = strcmpi(modelType, 'contaminated');
c = getConstants(p);
if ~exist('y0', 'var')
    y0 = getYinit(p);
end

% get signal by iteration
signal = zeros([p.D, p.nT]);
for t = 1:p.nT
    y0 = nextTimePoint(p, c, y0, S(:,t));
    signal(:,t) = y0(end-p.D+1:end);
end


end


function c = getConstants(p)
% store a few constants for hemodynamic model calculations (to only calculate them once) 
    nDir = 2;
    
    % prepare constants
    V0_ = p.balloon.V0 * p.D / 100;  % transformation from total volume to blood volume fraction
    sV0 = sum(V0_, 1);  % sum of V0_ over all compartments
    seV = sum(V0_ .* repmat(p.bold.eps', [1,p.D]), 1);  % sum of V0_*eps over all compartments
    H0 = 1 ./ (1 - sV0 + seV);  % scaling for BOLD-signal
    b1 = 4.3 * p.bold.dXi * p.bold.Hct * p.bold.gamma0 * p.bold.B0 .* p.bold.E0 * p.seq.TE;
    b2 = p.bold.eps .* p.bold.r0 .* p.bold.E0 * p.seq.TE;
    b3 = 1 - p.bold.eps;

    vet_3D = repmat(reshape(p.balloon.vet',[p.K,1,nDir]),[1,p.D,1]);
    F0_3D = repmat(p.balloon.F0,[1,1,nDir]);
    V0_3D = repmat(p.balloon.V0,[1,1,nDir]);
    F_denom = vet_3D .* F0_3D + V0_3D;
    F0_venule = repmat(F0_3D(p.VENULE,:,:),[p.K,1,1]);

    % init deeper layer constants
    c.fvq.f.deeperLayer = zeros([p.K, p.D, nDir]);
    c.fvq.vq.deeperLayer = zeros([p.K, p.D]);
    % define constants
    c.fvq.f.same = V0_3D ./ F_denom;
    c.fvq.f.prevComp = vet_3D .* F0_venule ./ F_denom;
    c.fvq.f.deeperLayer(p.VEIN,1:p.D-1,:) = vet_3D(p.VEIN,1:p.D-1,:) .* F0_3D(p.VEIN,2:p.D,:) ./ F_denom(p.VEIN,1:p.D-1,:);
    c.fvq.vq.same = - ones([p.K,p.D]) ./ p.balloon.tau0;
    c.fvq.vq.prevComp = [zeros([1,p.D]); p.balloon.F0(p.ARTERIOLE:p.VENULE,:) ./ p.balloon.F0(p.VENULE:p.VEIN,:) ./ p.balloon.tau0(p.VENULE:p.VEIN,:)];
    c.fvq.vq.deeperLayer(p.VEIN,1:p.D-1) = p.balloon.F0(p.VEIN, 2:p.D) ./ p.balloon.F0(p.VEIN, 1:p.D-1) ./ p.balloon.tau0(p.VEIN, 1:p.D-1);
    c.BOLD.q = repmat(H0,[p.K,1]) .* repmat(b1',[1,p.D]) .* V0_ .* (1 - repmat(sV0,[p.K,1]));
    c.BOLD.qv = repmat(H0,[p.K,1]) .* repmat(b2',[1,p.D]) .* V0_;
    c.BOLD.v = repmat(H0,[p.K,1]) .* repmat(b3',[1,p.D]) .* V0_;
end

function [i_start] = getYindex(p, funID)
    % funnames = {'n_excitation', 'n_inhibition', 'vaso_active', 'f_arteriole', 'f', 'v', 'q', 'BOLD'};
    % find is slow, lookup is faster
    if funID==1 %strcmpi(funname, 'n_excitation')
        i = 0;
    elseif funID==2 %strcmpi(funname, 'n_inhibition')
        i = 1;
    elseif funID==3 %strcmpi(funname, 'vaso_active')
        i = 2;
    elseif funID==4 %strcmpi(funname, 'f_arteriole')
        i = 3;
    elseif funID==5 %strcmpi(funname, 'f')
        i = 4;
    elseif funID==6 %strcmpi(funname, 'v')
        i = 6;
    elseif funID==7 %strcmpi(funname, 'q')
        i = 8;
    elseif funID==8 %strcmpi(funname, 'BOLD')
        i = 10;
    end
    i_start = i*p.D + 1;
end

function i_end = i_end(i_start, p)
    i_end = i_start+p.D-1;
end
function [y] = getVec(y, p, i_start)
    y = y(i_start:i_end(i_start, p));
end

function y0 = getYinit(p)
% get init values for each function
    funnames = {'n_excitation', 'n_inhibition', 'vaso_active', 'f_arteriole', 'f', 'v', 'q', 'BOLD'};
    y0 = zeros([1, length(funnames)*p.D+3*p.D]);  % f,v,q need 2*D values each
    % q, f_arteriole
    a_start = getYindex(p, 4);  %'f_arteriole');
    q_start = getYindex(p, 7);  %'q');
    y0(a_start:i_end(a_start, p)) = 1;
    y0(q_start:i_end(q_start, p)+p.D) = 1;
    % f
    F0 = p.balloon.F0;  % readability
    f_venule_start = getYindex(p, 5);  %'f');
    f_vein_start = f_venule_start+p.D;
    y0(f_venule_start:i_end(f_venule_start, p)) = F0(p.ARTERIOLE,:) ./ F0(p.VENULE,:);  % venule
    y0(i_end(f_vein_start, p)) = F0(p.ARTERIOLE,p.D)/F0(p.VEIN,p.D);  % deepest layer vein
    for d = p.D-1:-1:1  % higher layer veins
        y0(f_vein_start+d-1) = F0(p.ARTERIOLE,d)/F0(p.VEIN,d) + F0(p.VEIN,d+1)/F0(p.VEIN,d) * y0(f_vein_start+d) * p.contaminated;
    end
    % v
    v_venule_start = getYindex(p, 6);  %'v');
    v_vein_start = v_venule_start+p.D;
    y0(v_venule_start:i_end(v_venule_start, p)) = power( getVec(y0, p, f_venule_start), p.balloon.alpha(p.VENULE) );
    y0(v_vein_start  :i_end(v_vein_start  , p)) = power( getVec(y0, p, f_vein_start  ), p.balloon.alpha(p.VEIN  ) );
end

function dL = deeperLayer(y, p, i_start, defaultVal)
    dL = circshift( getVec(y, p, i_start), p.D-1 );
    dL(p.D) = defaultVal;
end

function flowdir = getFlowDir(dx_v)
    flowdir = squeeze(dx_v) < 0;  % inflation: 0, deflation: 1
end
function c_part = getC(p, c_part, flowdir)
% get the part of C that corresponds to flowdir on each depth level
    c_part = c_part(:);  % turn into 1D array (stack inflation/deflation after another)
    flowdir = linspace(1,p.D,p.D) + flowdir * p.D;  % get indice of stacked array
    c_part = c_part(flowdir)';
end

function y1 = updateVec(y0, y1, p, i_start, dx, type)
% types: 0 -> same, 1 -> Euler, 2 -> log-normal
    if type == 0  % same
        y1(i_start:i_end(i_start, p)) = dx;
    elseif type == 1  % Euler
        y1(i_start:i_end(i_start, p)) = getVec(y0, p, i_start) + dx * p.dt;
    elseif type == 2  % log-normal
        y1(i_start:i_end(i_start, p)) = getVec(y0, p, i_start) .* exp(dx * p.dt ./ getVec(y0, p, i_start));
    end
end

function y1 = nextTimePoint(p, c, y0, S_t)
    % funnames = {'n_excitation', 'n_inhibition', 'vaso_active', 'f_arteriole', 'f', 'v', 'q', 'BOLD'};
    y1 = zeros(size(y0));
    
    % neuro
    i_n_ex = getYindex(p, 1);  %'n_excitation');
    i_n_in = getYindex(p, 2);  %'n_inhibition');
    i_vas = getYindex(p, 3);  %'vaso_active');
    i_fa = getYindex(p, 4);  %'f_arteriole');
    dx_n_ex = p.neuro.sigma * getVec(y0, p, i_n_ex) - p.neuro.mu * getVec(y0, p, i_n_in) + p.neuro.C * S_t';
    dx_n_in = p.neuro.lambda * ( getVec(y0, p, i_n_ex) - getVec(y0, p, i_n_in) );
    dx_vas = getVec(y0, p, i_n_ex) - p.neuro.theta * getVec(y0, p, i_vas);
    dx_fa = p.neuro.phi * getVec(y0, p, i_vas) - p.neuro.xi * (getVec(y0, p, i_fa) - 1);
    y1 = updateVec(y0, y1, p, i_n_ex, dx_n_ex, 1);
    y1 = updateVec(y0, y1, p, i_n_in, dx_n_in, 1);
    y1 = updateVec(y0, y1, p, i_vas, dx_vas, 1);
    y1 = updateVec(y0, y1, p, i_fa, dx_fa, 2);
    % balloon
    i_f = getYindex(p, 5);  %'f');
    i_v = getYindex(p, 6);  %'v');
    i_q = getYindex(p, 7);  %'q');
    deeperLayer_f = deeperLayer(y0, p, i_f+p.D, 0);
    deeperLayer_v = deeperLayer(y0, p, i_v+p.D, 1);
    deeperLayer_q = deeperLayer(y0, p, i_q+p.D, 0);
    dx_v_venule = c.fvq.vq.same(p.VENULE,:)     .* getVec(y0, p, i_f) ...
                + c.fvq.vq.prevComp(p.VENULE,:) .* getVec(y0, p, i_fa);
    dx_v_vein = c.fvq.vq.same(p.VEIN,:)         .* getVec(y0, p, i_f+p.D) ...
              + c.fvq.vq.prevComp(p.VEIN,:)     .* getVec(y0, p, i_f) ...
              + c.fvq.vq.deeperLayer(p.VEIN,:)  .* deeperLayer_f * p.contaminated;
    flowdir_venule = getFlowDir(dx_v_venule);
    flowdir_vein   = getFlowDir(dx_v_vein);
    f_venule = getC(p, c.fvq.f.same(    p.VENULE,:,:), flowdir_venule) .* power( getVec(y0, p, i_v),     1/p.balloon.alpha(p.VENULE) ) ...
             + getC(p, c.fvq.f.prevComp(p.VENULE,:,:), flowdir_venule) .* getVec(y0, p, i_fa);
    f_vein   = getC(p, c.fvq.f.same(       p.VEIN,:,:), flowdir_vein) .* power( getVec(y0, p, i_v+p.D), 1/p.balloon.alpha(p.VEIN) ) ...
             + getC(p, c.fvq.f.prevComp(   p.VEIN,:,:), flowdir_vein) .* getVec(y0, p, i_f) ...
             + getC(p, c.fvq.f.deeperLayer(p.VEIN,:,:), flowdir_vein) .* deeperLayer_f * p.contaminated;
    n = p.balloon.n(p.VENULE);
    dx_q_venule = c.fvq.vq.same(p.VENULE,:) .* getVec(y0, p, i_f) .* getVec(y0, p, i_q) ./ getVec(y0, p, i_v) ...
                + c.fvq.vq.prevComp(p.VENULE,:) .* getVec(y0, p, i_fa) .* ( 1/n + 1./getVec(y0, p, i_fa) - 1/n./getVec(y0, p, i_fa) ) ;
    dx_q_vein = c.fvq.vq.same(p.VEIN,:) .* getVec(y0, p, i_f+p.D) .* getVec(y0, p, i_q+p.D) ./ getVec(y0, p, i_v+p.D) ...
              + c.fvq.vq.prevComp(p.VEIN,:) .* getVec(y0, p, i_f) .* getVec(y0, p, i_q) ./ getVec(y0, p, i_v) ...
              + c.fvq.vq.deeperLayer(p.VEIN,:) .* deeperLayer_f .* deeperLayer_q ./ deeperLayer_v * p.contaminated;
    y1 = updateVec(y0, y1, p, i_f, f_venule, 0);
    y1 = updateVec(y0, y1, p, i_f+p.D, f_vein, 0);
    y1 = updateVec(y0, y1, p, i_v, dx_v_venule, 2);
    y1 = updateVec(y0, y1, p, i_v+p.D, dx_v_vein, 2);
    y1 = updateVec(y0, y1, p, i_q, dx_q_venule, 2);
    y1 = updateVec(y0, y1, p, i_q+p.D, dx_q_vein, 2);
    % BOLD
    i_b = getYindex(p, 8);  %'BOLD');
    b_venule = c.BOLD.q(p.VENULE,:) .* ( 1 - getVec(y1, p, i_q) ) ...
             + c.BOLD.qv(p.VENULE,:) .* ( 1 - getVec(y1, p, i_q) ./ getVec(y1, p, i_v) ) ...
             + c.BOLD.v(p.VENULE,:) .* ( 1 - getVec(y1, p, i_v) );
    b_vein   = c.BOLD.q(p.VEIN,:) .* ( 1 - getVec(y1, p, i_q+p.D) ) ...
             + c.BOLD.qv(p.VEIN,:) .* ( 1 - getVec(y1, p, i_q+p.D) ./ getVec(y1, p, i_v+p.D) ) ...
             + c.BOLD.v(p.VEIN,:) .* ( 1 - getVec(y1, p, i_v+p.D) );
    signal = b_venule + b_vein;  % sum over compartments
    y1 = updateVec(y0, y1, p, i_b, signal, 0);
    
end
