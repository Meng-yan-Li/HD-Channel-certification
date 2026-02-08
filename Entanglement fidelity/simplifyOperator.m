function [Operator] = simplifyOperator(midOperator, d, GammaType)
% SIMPLIFYOPERATOR
% Simplify and validate composite operators generated in the hierarchy
% construction by enforcing algebraic consistency rules and equivalence
% relations.
%
% Each operator is represented as a pair {v1, v2}, corresponding to
% symbolic products of operators acting on Alice (v1) and Bob (v2).
%
% This function:
%   (i)   removes trivial identity contributions,
%   (ii)  checks algebraic compatibility of adjacent operators,
%   (iii) merges compatible basis operators,
%   (iv)  removes redundant repetitions,
%   (v)   enforces cyclic equivalence (if required).
%
% Invalid operator products are discarded by returning an empty cell {}.
%
% Inputs:
%   midOperator : 1x2 cell {v1, v2}, each a vector encoding operator indices
%   d           : local Hilbert space dimension
%   GammaType   : hierarchy type selector
%                 = 0 ：no cyclic trace
%                 = 1 : cyclic trace equivalence for Φ
%                 = 2 : cyclic trace equivalence for IZ
%
% Output:
%   Operator    : simplified {v1, v2} operator, or {} if invalid

    %% ================= Step 1: unpack operator components =================
    v1 = midOperator{1,1};   % Alice-side operator string
    v2 = midOperator{1,2};   % Bob-side operator string

    %% ================= Step 2: remove explicit identities =================
    % Replace [0;0] by scalar 0 (identity) and remove zero entries
    if isequal(v1, [0;0]), v1 = 0; end
    if numel(v1) >= 2, v1 = v1(v1 ~= 0); end

    if isequal(v2, [0;0]), v2 = 0; end
    if numel(v2) >= 2, v2 = v2(v2 ~= 0); end

    %% ================= Step 3: validate adjacent Alice operators ===========
    % Enforce algebraic consistency for products of E_{i,j} operators
    if GammaType == 1 && isequal(v2,0) || GammaType == 2
        len = numel(v1);          % cyclic trace
    else
        len = numel(v1) - 1;      % no cyclic trace
    end

    for i = 1:len
        a = v1(i);
        b = v1(mod(i, numel(v1)) + 1);

        % Both entries correspond to matrix units E_{i,j}
        if a < 0 && b < 0
            Ea = mod(floor((-a-1) ./ d.^(0:1)), d) + 1;
            Eb = mod(floor((-b-1) ./ d.^(0:1)), d) + 1;

            % Invalid contraction: inner indices do not match
            if Ea(2) ~= Eb(1)
                Operator = {};
                return;
            else
                % Merge E_{i,j} E_{j,k} -> E_{i,k}
                v1(i) = -(Ea(1) + (Eb(2)-1)*d);
                v1(mod(i, numel(v1)) + 1) = v1(i);
            end
        end
    end

    %% ================= Step 4: validate adjacent Bob operators ==============
    % Enforce projective measurement constraints:
    % incompatible outcomes of the same measurement setting vanish
    len = numel(v2) - 1;
    for i = 1:len
        a = v2(i);
        b = v2(i+1);
        if a > 0 && b > 0
            if floor((a-1)/(d-1)) + 1 == floor((b-1)/(d-1)) + 1 && ...
               mod(a-1, d-1) + 1 ~= mod(b-1, d-1) + 1
                Operator = {};
                return;
            end
        end
    end

    %% ================= Step 5: enforce cyclic equivalence ====================
    % For trace-like operators, rotate to lexicographically minimal form
    if GammaType == 1 && isequal(v2,0) || GammaType == 2
        v1 = rotate_min_lex(v1(:));
    end

    %% ================= Step 6: remove adjacent duplicates ===================
    % Idempotency: AA = A, M_b^2 = M_b
    if numel(v1) >= 2
        keep = [true; ~(v1(2:end) == v1(1:end-1))];
        v1 = v1(keep);
    end

    if numel(v2) >= 2
        keep = [true; ~(v2(2:end) == v2(1:end-1))];
        v2 = v2(keep);
    end

    %% ================= Final simplified operator ============================
    Operator = {v1, v2};
end

%% Helper: find lexicographically minimal cyclic rotation
function vmin = rotate_min_lex(v)
% Rotate vector v over all cyclic shifts and return the lexicographically
% smallest representative (used to quotient out cyclic trace equivalence).

    if isempty(v) || numel(v) <= 1
        vmin = v;
        return;
    end
    vmin = v;
    L = numel(v);
    for r = 1:(L-1)
        cand = circshift(v, r);
        if lexicographicLess(cand, vmin)
            vmin = cand;
        end
    end
end

%% Helper: lexicographic comparison (a < b ?)
function isLess = lexicographicLess(a, b)
% Compare two vectors lexicographically.

    na = numel(a); nb = numel(b);
    n = min(na, nb);
    if n == 0
        isLess = na < nb;
        return;
    end

    idx = find(a(1:n) ~= b(1:n), 1);
    if isempty(idx)
        isLess = na < nb;     % shorter prefix is smaller
    else
        isLess = a(idx) < b(idx);
    end
end
