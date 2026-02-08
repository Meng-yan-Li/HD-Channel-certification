function [kOperatorList] = generateOperatorList(n, d, k)
% GENERATEOPERATORLIST
% Generate the operator list used in the SDP hierarchy, determined by the
% experimental configuration (n, d) and the hierarchy level k.
%
% Each operator is represented symbolically by a pair {A, B}, corresponding
% to operators of the form:
%     A ? B
% where A acts on the preparation (Alice) space and B on the measurement
% (Bob) space.
%
% Encoding convention:
%   {0, 0}     : Identity operator I_A ? I_B
%   {-e, 0}    : Basis operator E_{i,j} ? I_B, where
%                E_{i,j} = |i><j| and
%                (i,j) = mod(floor((e-1)./d.^(0:1)), d) + 1
%   {x, 0}     : Preparation state ¦Ñ_x ? I_B
%   {0, by}    : Measurement operator I_A ? M_{b|y}, with
%                by = (y-1)(d-1) + b
%
% Index decoding:
%   y = floor((by-1)/(d-1)) + 1
%   b = mod(by-1, d-1) + 1
%
% Assumptions:
%   - Preparation states ¦Ñ_x are pure
%   - Measurements M_{b|y} are projective
%
% Inputs:
%   n : number of preparation settings
%   d : local Hilbert space dimension
%   k : hierarchy level
%
% Output:
%   kOperatorList : cell array containing all operators up to level k

    %% Level-0 operators: identity and matrix basis
    L1 = {};                        % Basic operator list
    L1(1,1:2) = {0, 0};             % Identity I_A ? I_B
    index = 2;
    for e = 1:d^2
        L1(index,1:2) = {-e, 0};    % E_{i,j} ? I_B
        index = index + 1;
    end

    %% Level-1 operators: preparations and measurements
    L2 = {};
    index = 1;

    % Preparation operators ¦Ñ_x ? I_B
    for x = 1:d^n
        L2(index,1:2) = {x, 0};
        index = index + 1;
    end

    % Measurement operators I_A ? M_{b|y}
    for y = 1:n
        for b = 1:d-1               % Last outcome fixed by completeness
            L2(index,1:2) = {0, (y-1)*(d-1) + b};
            index = index + 1;
        end
    end

    % Mixed preparation-measurement terms (k = 1 + AB)
    if k == 1.5
        for x = 1:d^n
            for y = 1:n
                for b = 1:d-1
                    L2(index,1:2) = {x, (y-1)*(d-1) + b};
                    index = index + 1;
                end
            end
        end
    end

    listLength = length(L2);

    %% Higher-order operators: Cartesian products up to level k
    kOperatorList = L2;             % Initialize hierarchy list
    for i = 2:k
        if i == 2
            O_k = L2;               % First Cartesian product
        else
            O_k = midList;
        end

        midList = {};               % Intermediate operator list
        midLength = length(O_k);
        index = 1;

        % Cartesian product O_k ¡Á L2
        for j = 1:midLength
            for l = 1:listLength
                midOperator = cell(1,2);
                midOperator{1,1} = [O_k{j,1}; L2{l,1}];
                midOperator{1,2} = [O_k{j,2}; L2{l,2}];

                % Simplify using algebraic identities
                midOperator = simplifyOperator(midOperator, d, 0);
                if ~isempty(midOperator)
                    midList(index,:) = midOperator;
                    index = index + 1;
                end
            end
        end

        % Remove duplicate operators
        midList = deleteDuplicates(midList);
        kOperatorList = [kOperatorList; midList];
    end

    %% Append level-0 operators and final cleanup
    kOperatorList = [L1; kOperatorList];
    if k >= 2
        kOperatorList = deleteDuplicates(kOperatorList);
    end
end
