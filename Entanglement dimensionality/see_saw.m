function [Pbest, rho_cell, M_cell, Phi, history] = see_saw(n, d, r, opts)
% SEE_SAW
% See-saw optimization algorithm for the target prepare-and-measure problem.
%
% This routine alternates between optimizing:
%   (1) preparation states rho_x,
%   (2) measurement POVMs M_{b|y},
%   (3) shared bipartite state Phi (with Schmidt number <= r),
% until convergence of the objective value.
%
% Inputs:
%   n, d, r : problem parameters
%       n   - number of measurements
%       d   - local dimension
%       r   - Schmidt number / rank constraint
%   opts    : struct of optional parameters
%       .maxIter     - maximum number of see-saw iterations (default: 200)
%       .tol         - convergence tolerance (default: 1e-6)
%       .verbose     - display iteration info (default: false)
%       .solver      - SDP solver ('scs','mosek','sdpt3')
%       .initSeed    - RNG seed (default: 'shuffle')
%
% Outputs:
%   Pbest     : best objective value achieved
%   rho_cell  : cell array of preparation states rho_x (each d x d)
%   M_cell    : cell array of POVMs, M_cell{y}{b} is d x d
%   Phi       : optimized bipartite state (d^2 x d^2)
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
        fprintf('See-Saw: n=%d d=%d r=%d  (N = d^n = %d)\n', n, d, r, N);
    end

    % Overall normalization factor in the objective
    scale = d / (n * N);

    % Pre-generate all x strings in {1,...,d}^n (row-wise)
    Xlist = index2code_array(1:N, n, d);  % N x n matrix

    % -------------------- Initialization: rho_x --------------------
    % Initialize each preparation as a random pure state
    rho_cell = cell(N,1);
    for x = 1:N
        v = randn(d,1) + 1i*randn(d,1);
        v = v / norm(v);
        rho_cell{x} = (v*v');  % rank-1 density matrix
    end

    % -------------------- Initialization: POVMs --------------------
    % Initialize POVMs as random projective measurements
    M_cell = cell(n, 1);
    for y = 1:n
        M_cell{y} = cell(d, 1);

        % Generate a random unitary U ∈ U(d)
        X = randn(d) + 1i * randn(d);
        [U, ~] = qr(X);

        % Define M_b = |u_b><u_b| for each outcome b
        for b = 1:d
            ub = U(:, b);
            M_cell{y}{b} = ub * ub';  % rank-1 projector
        end
    end

    % -------------------- Initialization: Phi --------------------
    % Initialize bipartite state Phi with Schmidt number <= r
    [Phi, ~] = random_SN_channel(d, r, opts, d, true);

    % -------------------- History tracking --------------------
    history.obj = [];
    history.time = [];

    % -------------------- Main see-saw loop --------------------
    prevObj = -inf;
    tic;
    for iter = 1:opts.maxIter

        % ---------- (1) Optimize preparation states rho_x ----------
        for x = 1:N
            Cx = zeros(d, d);
            I = eye(d);

            % Accumulate contributions from all measurements y
            for y = 1:n
                Mb = M_cell{y}{Xlist(x, y)};
                % Cx += Tr_B[(I ⊗ M_b) Phi]
                Cx = Cx + PartialTrace(Tensor(I, Mb) * Phi, 2, [d,d]);
            end

            % Optimal rho_x is projector onto largest eigenvector of Cx
            [V, D] = eig((Cx + Cx')/2);  % ensure Hermiticity
            [~, idx] = max(real(diag(D)));
            v = V(:, idx);
            rho_cell{x} = v * v';
        end

        % ---------- (2) Optimize POVMs M_{b|y} ----------
        % For each measurement y, construct effective operators G_y{b}
        for y = 1:n
            G_y = cell(d, 1);
            for b = 1:d
                x_indices = find(Xlist(:, y) == b);
                rho_yb = zeros(d);

                % Sum all rho_x consistent with outcome b at position y
                for idx = 1:length(x_indices)
                    rho_yb = rho_yb + rho_cell{x_indices(idx)};
                end

                % G_y{b} = Tr_A[(rho_yb ⊗ I) Phi]
                G_y{b} = PartialTrace(Tensor(rho_yb, eye(d)) * Phi, 1, [d,d]);
            end

            % Optimize POVM for fixed G_y
            M_cell{y} = optimize_POVM(G_y, d, opts.solver);
        end

        % ---------- (3) Optimize bipartite state Phi ----------
        % Construct operator A = sum_{x,y} rho_x ⊗ M_{x_y|y}
        A = sparse(d^2, d^2);
        for x = 1:N
            for y = 1:n
                Mb = M_cell{y}{Xlist(x, y)};
                A = A + Tensor(rho_cell{x}, Mb);
            end
        end

        % Optimize Phi under Schmidt number constraint
        Phi = optimize_phi(A, d, r, opts.solver);

        % ---------- Evaluate objective ----------
        obj = 0;
        for x = 1:N
            for y = 1:n
                Mb = M_cell{y}{Xlist(x, y)};
                obj = obj + real(trace(Tensor(rho_cell{x}, Mb) * Phi));
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
