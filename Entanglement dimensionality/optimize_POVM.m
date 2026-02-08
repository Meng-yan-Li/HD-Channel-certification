function Mlist = optimize_POVM(G_y, d, solver)
% OPTIMIZE_POVM
% Optimize a POVM {M_b} to maximize a linear objective of the form
%   sum_b Tr(G_b * M_b),
% subject to the standard POVM constraints.
%
% This function solves a small semidefinite program (SDP) where each
% POVM element M_b is a d × d positive semidefinite matrix and the
% completeness relation sum_b M_b = I is enforced.
%
% Inputs:
%   G_y    : cell(d,1) of fixed d × d matrices (treated as constants in CVX)
%   d      : local Hilbert space dimension
%   solver : SDP solver used by CVX (e.g. 'mosek', 'sdpt3', 'scs')
%
% Outputs:
%   Mlist  : cell(d,1), where each entry is a d × d POVM element

    % Cell array to store CVX variables corresponding to POVM elements
    Mcell = cell(d,1);

    % -------------------- CVX optimization problem --------------------
    % Declare one PSD matrix variable for each outcome b
    cvx_begin sdp quiet
        cvx_solver(solver);

        % Dynamically declare POVM elements M_b as CVX variables
        for b = 1:d
            eval(sprintf('variable Mb%d(%d,%d) hermitian semidefinite', b, d, d));
            % Store variable handle in cell array for convenient access
            eval(sprintf('Mcell{%d} = Mb%d;', b, b));
        end

        % -------------------- Objective function --------------------
        % Maximize sum_b Tr(G_b * M_b)
        obj = 0;
        for b = 1:d
            Gb = G_y{b};   % constant matrix
            obj = obj + trace(Gb * Mcell{b});
        end
        maximize(real(obj));

        % -------------------- POVM completeness constraint --------------------
        % Enforce sum_b M_b = I
        S = zeros(d);
        for b = 1:d
            S = S + Mcell{b};
        end
        S == eye(d);
    cvx_end

    % -------------------- Extract optimized POVM --------------------
    Mlist = cell(d,1);
    for b = 1:d
        Mlist{b} = Mcell{b};
    end
end
