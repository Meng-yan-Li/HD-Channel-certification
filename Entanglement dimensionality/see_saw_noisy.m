function [Pbest, rho_cell, M_cell, history] = see_saw_noisy(n, d, v, channel, opts)
% SEE_SAW_NOISY
% See-saw optimization in the presence of a noisy communication channel.
%
% This function implements a see-saw (alternating) optimization over:
%   (1) preparation states rho_x,
%   (2) measurement POVMs M_{b|y},
% while all operators are affected by a fixed noisy channel.
%
% The noise is incorporated via the map commu_channel(·), which represents
% either the channel itself or its dual action, depending on context.
%
% Inputs:
%   n, d, v, channel : problem parameters
%       n       - number of measurements
%       d       - system dimension
%       v       - noise parameter (e.g. visibility)
%       channel - identifier of the noisy channel model
%   opts        : struct of optional parameters
%       .maxIter   - maximum number of see-saw iterations (default: 200)
%       .tol       - convergence tolerance (default: 1e-6)
%       .verbose   - display iteration information (default: false)
%       .solver    - SDP solver ('scs','mosek','sdpt3')
%       .initSeed  - RNG seed (default: 'shuffle')
%
% Outputs:
%   Pbest     : best objective value achieved
%   rho_cell  : cell array of preparation states rho_x (each d x d)
%   M_cell    : cell array of POVMs, M_cell{y}{b} is d x d
%   history   : struct containing objective values and timing

    % -------------------- Default options --------------------
    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'maxIter'), opts.maxIter = 200; end
    if ~isfield(opts,'tol'), opts.tol = 1e-6; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end
    if ~isfield(opts,'solver'), opts.solver = 'mosek'; end
    if ~isfield(opts,'initSeed'), opts.initSeed = 'shuffle'; end
    rng(opts.initSeed);

    % -------------------- Basic parameters --------------------
    % Total number of preparation settings x ∈ {1,...,d}^n
    N = d^n;

    if opts.verbose
        fprintf('See-Saw: n=%d d=%d v=%d  (N = d^n = %d)\n', n, d, v, N);
    end

    % Overall normalization factor in the objective
    scale = 1 / (n * N);

    % Pre-generate all x strings in {1,...,d}^n (row-wise)
    Xlist = index2code_array(1:N, n, d);  % N x n matrix

    % -------------------- Initialization: rho_x --------------------
    % Initialize preparation states as random pure states
    rho_cell = cell(N,1);
    for x = 1:N
        vec = randn(d,1) + 1i*randn(d,1);
        vec = vec / norm(vec);
        rho_cell{x} = (vec*vec');  % rank-1 density matrix
    end

    % -------------------- Initialization: POVMs --------------------
    % Initialize POVMs as random projective measurements
    M_cell = cell(n, 1);
    for y = 1:n
        M_cell{y} = cell(d, 1);

        % Generate a random unitary U ∈ U(d)
        X = randn(d) + 1i * randn(d);
        [U, ~] = qr(X);

        % Define M_b = |u_b><u_b|
        for b = 1:d
            ub = U(:, b);
            M_cell{y}{b} = ub * ub';  % rank-1 projector
        end
    end

    % -------------------- History tracking --------------------
    history.obj = [];
    history.time = [];

    % -------------------- Main see-saw loop --------------------
    prevObj = - inf;
    tic;
    for iter = 1:opts.maxIter

        % ---------- (1) Optimize preparation states rho_x ----------
        for x = 1:N
            Cx = zeros(d);

            % Accumulate measurement operators for fixed x
            for y = 1:n
                Mb = M_cell{y}{Xlist(x, y)};
                Cx = Cx + Mb;
            end

            % Apply dual (or effective) noisy channel map
            Cx = commu_channel(Cx, v, d, channel);

            % Optimal rho_x is projector onto largest eigenvector of Cx
            [V, D] = eig((Cx + Cx')/2);  % ensure Hermiticity
            [~, idx] = max(real(diag(D)));
            vec = V(:, idx);
            rho_cell{x} = vec * vec';
        end

        % ---------- (2) Optimize POVMs M_{b|y} ----------
        for y = 1:n
            G_y = cell(d, 1);
            for b = 1:d
                x_indices = find(Xlist(:, y) == b);
                rho_yb = zeros(d);

                % Sum all rho_x consistent with outcome b at position y
                for idx = 1:length(x_indices)
                    rho_yb = rho_yb + rho_cell{x_indices(idx)};
                end

                % Apply noisy channel to effective operator
                G_y{b} = commu_channel(rho_yb, v, d, channel);
            end

            % Optimize POVM for fixed G_y
            M_cell{y} = optimize_POVM(G_y, d, opts.solver);
        end

        % ---------- Evaluate objective ----------
        obj = 0;
        for x = 1:N
            for y = 1:n
                Mb = M_cell{y}{Xlist(x, y)};
                obj = obj + real(trace( ...
                    commu_channel(rho_cell{x}, v, d, channel) * Mb ));
            end
        end
        obj = scale * obj;

        history.obj(end+1) = obj;
        history.time(end+1) = toc;

        if opts.verbose
            fprintf('iter %3d: obj = %.12f, delta = %.3e, time = %.2fs\n', ...
                iter, obj, abs(obj-prevObj), toc);
        end

        % Convergence check
        if abs(obj - prevObj) < opts.tol
            break;
        end
        prevObj = obj;
    end

    % Return best value
    Pbest = prevObj;
end
