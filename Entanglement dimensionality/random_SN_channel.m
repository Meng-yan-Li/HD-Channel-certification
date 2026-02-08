function [sigma, KrausOps] = random_SN_channel(d, r, opts, numKraus, exactRank)
% RANDOM_SN_CHANNEL
% Randomly generate a CPTP quantum channel whose Choi state has
% Schmidt number equal to r (if exactRank = true) or at most r.
%
% This function constructs a random set of Kraus operators with
% controlled rank, enforces the CPTP condition, and then builds
% the corresponding Choi state via the CJ isomorphism.
%
% Inputs:
%   d          : system dimension
%   r          : target Schmidt number
%   opts       : option structure (used here for RNG seed)
%   numKraus   : number of Kraus operators (default: 5)
%   exactRank  : if true, each Kraus operator has rank exactly r;
%                if false, rank is at most r
%
% Outputs:
%   sigma      : Choi state (d^2 × d^2 Hermitian matrix, trace 1)
%   KrausOps   : cell array containing the Kraus operators

    % -------------------- Default parameters --------------------
    if nargin < 5, exactRank = false; end
    if nargin < 4, numKraus = 5; end

    % Initialize RNG
    rng(opts.initSeed);

    % -------------------- 1. Generate Kraus operators --------------------
    KrausOps = cell(1, numKraus);
    for k = 1:numKraus
        if exactRank
            % Construct a rank-exactly-r operator via truncated SVD
            [U,~] = qr(randn(d) + 1i*randn(d));
            [V,~] = qr(randn(d) + 1i*randn(d));
            S = diag(1 + 0.1*rand(1,r));  % nonzero singular values
            A = U(:,1:r) * S * V(:,1:r)';
        else
            % Construct an operator with rank at most r
            X = randn(d,r) + 1i*randn(d,r);
            Y = randn(r,d) + 1i*randn(r,d);
            A = X * Y;
        end
        KrausOps{k} = A;
    end

    % -------------------- 2. CPTP normalization --------------------
    % Enforce sum_k A_k^† A_k = I via inverse square root normalization
    T = zeros(d);
    for k = 1:numKraus
        T = T + KrausOps{k}' * KrausOps{k};
    end

    % Compute inverse square root of T (numerically stabilized)
    [Vt, Dt] = eig(T);
    Dt = real(diag(Dt));
    Dt(Dt < 1e-14) = 1e-14;      % avoid numerical singularities
    T_inv_sqrt = Vt * diag(1./sqrt(Dt)) * Vt';

    for k = 1:numKraus
        KrausOps{k} = KrausOps{k} * T_inv_sqrt;
    end

    % -------------------- 3. Construct maximally entangled state |phi+> --------------------
    % |phi+> = (1/sqrt(d)) sum_i |i> ⊗ |i>
    phi = zeros(d^2, 1);
    for i = 1:d
        phi((i-1)*d + i) = 1;
    end
    phi = phi / sqrt(d);

    % Density matrix |phi+><phi+|
    rho_phi = phi * phi';

    % -------------------- 4. Construct Choi state --------------------
    % sigma = (I ⊗ Λ)(|phi+><phi+|)
    sigma = zeros(d^2);
    for k = 1:numKraus
        Ak = KrausOps{k};
        K = kron(eye(d), Ak);
        sigma = sigma + K * rho_phi * K';
    end

    % Numerical symmetrization and trace normalization
    sigma = (sigma + sigma')/2;    % enforce Hermiticity
    sigma = sigma / trace(sigma);  % enforce trace 1
end