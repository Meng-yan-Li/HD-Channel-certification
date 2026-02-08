function [GammaPHI, GammaIZ, GammaPHI_char, GammaIZ_char] = generateMomentMatrix(kOperatorList, d)
% GENERATEMOMENTMATRIX
% Construct the moment matrices Γ^Φ and Γ^IZ associated with the operator
% list up to hierarchy level k.
%
% Each entry of the moment matrix corresponds to an operator product of the
% form:
%       ⟨ O_i^† O_j ⟩
% where O_i, O_j are symbolic operators from kOperatorList.
%
% Inputs:
%   kOperatorList : cell array (K x 2), symbolic operator list {v1, v2}
%   d             : local Hilbert space dimension
%
% Outputs:
%   GammaPHI      : K x K cell array of simplified operators for Φ-matrix
%   GammaIZ       : K x K cell array of simplified operators for IZ-matrix
%   GammaPHI_char : string-encoded version of GammaPHI (for indexing)
%   GammaIZ_char  : string-encoded version of GammaIZ (for indexing)

    klistLength = size(kOperatorList, 1);

    % Initialize moment matrices
    GammaPHI = cell(klistLength);
    GammaIZ  = cell(klistLength);

    % Character-string representations (used for indexing / constraints)
    GammaPHI_char = cell(klistLength);
    GammaIZ_char  = cell(klistLength);

    %% Construct moment matrix entries
    for j = 1:klistLength          % Column index
        % Conjugate operator O_j^† :
        % reverse operator order and fix index convention
        conjOperator = { flipud(kOperatorList{j,1}), ...
                         flipud(kOperatorList{j,2}) };

        conjOperator = { fix_negative_indices(conjOperator{1,1}, d), ...
                         fix_negative_indices(conjOperator{1,2}, d) };

        for i = 1:klistLength      % Row index
            % Form symbolic product O_j^† O_i
            GammaPHI{i,j} = { [conjOperator{1,1}; kOperatorList{i,1}], ...
                              [conjOperator{1,2}; kOperatorList{i,2}] };
            GammaIZ{i,j}  = { [conjOperator{1,1}; kOperatorList{i,1}], ...
                              [conjOperator{1,2}; kOperatorList{i,2}] };

            % Apply algebraic simplification rules
            GammaPHI{i,j} = simplifyOperator(GammaPHI{i,j}, d, 1);
            GammaIZ{i,j}  = simplifyOperator(GammaIZ{i,j},  d, 2);

            % Store string representation for indexing
            if isempty(GammaPHI{i,j})
                GammaPHI_char{i,j} = '[]';
            else
                GammaPHI_char{i,j} = ...
                    [mat2str(GammaPHI{i,j}{1}) '_' mat2str(GammaPHI{i,j}{2})];
            end

            if isempty(GammaIZ{i,j})
                GammaIZ_char{i,j} = '[]';
            else
                GammaIZ_char{i,j} = ...
                    [mat2str(GammaIZ{i,j}{1}) '_' mat2str(GammaIZ{i,j}{2})];
            end
        end
    end
end


function v_new = fix_negative_indices(v, d)
% FIX_NEGATIVE_INDICES
% Convert negative operator indices into a canonical encoding using a
% base-d mapping rule.
%
% Negative entries correspond to matrix-unit operators E_{i,j}.
% This function enforces a consistent index ordering under conjugation.
%
% Mapping rule:
%   For each e_i < 0:
%       E = mod(floor((-e_i-1) ./ d.^(0:1)), d) + 1;
%       e_i -> -(E(2) + (E(1)-1)*d)
%
% Inputs:
%   v : vector of operator indices
%   d : local Hilbert space dimension
%
% Output:
%   v_new : vector with negative indices fixed

    v_new = v;                 % Copy input vector
    neg_idx = v < 0;           % Locate negative entries

    if any(neg_idx)
        e = -v(neg_idx);
        e = e - 1;
        E1 = mod(floor(e ./ d.^0), d) + 1;
        E2 = mod(floor(e ./ d.^1), d) + 1;

        % Re-encode indices in canonical form
        v_new(neg_idx) = -(E2 + (E1 - 1) * d);
    end
end
