function [Phi] = optimize_phi(A, d, r, solver)
% OPTIMIZE_PHI
% Optimize a bipartite quantum state Phi to maximize Tr(A * Phi)
% under a Schmidt number (or related entanglement) constraint.
%
% This function formulates and solves a semidefinite program (SDP)
% to find a two-party quantum state Phi acting on C^d ⊗ C^d that
% maximizes the linear objective Tr(A * Phi), subject to:
%   - positivity and normalization,
%   - an entanglement constraint controlled by the parameter r.
%
% Inputs:
%   A      : objective matrix (d^2 × d^2)
%   d      : local Hilbert space dimension
%   r      : upper bound on the Schmidt number
%   solver : SDP solver used by CVX (e.g. 'mosek', 'sdpt3', 'scs')
%
% Outputs:
%   Phi    : optimized bipartite quantum state (d^2 × d^2)

    % -------------------- Input dimension check --------------------
    if size(A,1) ~= d^2 || size(A,2) ~= d^2
        error('The dimension of matrix A should be d^2 × d^2');
    end

    % -------------------- CVX optimization problem --------------------
    cvx_begin sdp quiet
        cvx_solver(solver)

        % Optimization variables
        variable Phi(d^2,d^2) hermitian semidefinite
        variable omega(d^2*r^2, d^2*r^2) hermitian semidefinite

        % Objective: maximize linear functional Tr(A * Phi)
        maximize real(trace(A*Phi))

        subject to
            % Normalization of the quantum state
            trace(Phi) == 1

            % Optional marginal constraint (currently inactive)
            PartialTrace(Phi,2,[d,d]) == eye(d)/d

            % ---------------------------------------------------------
            % % Entanglement constraints (commented alternatives)
            % % DPS hierarchy (k = 1) with PPT constraint
            % Pi = Tensor(eye(d), Tensor(sqrt(r)*MaxEntangled(r), eye(d))); 
            % Pi' * omega * Pi == Phi; 
            % trace(omega) == r; 
            % PartialTranspose(omega/r,2,[d*r,d*r]) == hermitian_semidefinite(d^2*r^2)
            % ---------------------------------------------------------

            % Reduction map constraint (active constraint)
            % Ensures Schmidt number of Phi is bounded by r
            rhoA = PartialTrace(Phi,2,[d,d]);
            R = Tensor(rhoA, eye(d)) - (1/r) * Phi;
            R == hermitian_semidefinite(d^2);
    cvx_end

    % -------------------- Solver status check --------------------
    % If CVX does not report a successful solution, issue a warning
    if ~strcmp(cvx_status, 'Solved')
        warning('optimize_phi: CVX status = %s', cvx_status);
        if strcmp(cvx_status, 'Failed')
            Phi = [];
        end
    end
end
