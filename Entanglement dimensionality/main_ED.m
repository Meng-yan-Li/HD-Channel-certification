%% =========================
%  TASK 1
%  See-saw optimization for high-dimensional prepare-and-measure task
% =========================

clear;
clc;

% -------------------- General settings --------------------
repetition = 20;          % Number of random initializations (if used)
opts.maxIter = 200;       % Maximum see-saw iterations
opts.tol = 1e-8;          % Convergence tolerance
opts.verbose = true;      % Print iteration details
opts.solver = 'mosek';    % SDP solver
opts.initSeed = 'shuffle';% Random seed control

% -------------------- Problem parameters --------------------
n = 2;                    % Number of measurements
d = 5;                    % Local dimension
r = 5;                    % Schmidt number / rank parameter

% -------------------- Single run of see-saw optimization --------------------
% Note that if you want to implement the generalized DPS hierarchy, 
% you can make changes to it in the file 'optimize_phi.m'.

[Pbest, rho_cell, M_cell, Phi, history] = see_saw(n, d, r, opts);

% -------------------------------------------------------------------------
% The following blocks (commented out) are used for batch simulations
% over different dimensions d, Schmidt numbers r, and repetitions.
% -------------------------------------------------------------------------

% n = 2;
% for d = 2:5
%     ASP = arrayfun(@(r) arrayfun(@(i) see_saw(2,d,r,opts), ...
%         1:repetition), 1:d, 'UniformOutput', false);
%     save(sprintf('./data/ASP_n2_d%d.mat', d), 'ASP');
% end

%% =========================
%  TASK 2
%  Noisy scenario:
%  See-saw optimization under a communication channel
% =========================

clear;
clc;

% -------------------- Problem parameters --------------------
n = 2; 
d = 5; 
channel = 'deph';         % Channel type: 'depo' (depolarizing) or 'deph' (dephasing)
repetition = 20;          % Number of repetitions for averaging

% -------------------- Solver options --------------------
opts.maxIter = 200;
opts.tol = 1e-8;
opts.verbose = true;
opts.solver = 'mosek';
opts.initSeed = 'shuffle';    % Or set to a fixed integer (e.g., 0)

% -------------------- Storage for results --------------------
% ASP_v(j,i): j-th repetition at noise parameter v = (i-1)/100
ASP_v = zeros(repetition, 101);

% -------------------- Sweep over noise parameter v --------------------
i = 1;
for v = 0:0.01:1
    for j = 1:repetition
        ASP_v(j, i) = see_saw_noisy(n, d, v, channel, opts);
    end
    i = i + 1;
end

% -------------------- Save data --------------------
save('./data/deph_ASP_v_n2_d5.mat', 'ASP_v');
