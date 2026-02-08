function codes = index2code_array(indices, n, d)
% INDEX2CODE_ARRAY
% Convert linear indices into base-d codes of fixed length n.
%
% Each index in {1,...,d^n} is mapped to an n-dimensional vector whose
% entries take values in {1,...,d}. The output corresponds to the
% base-d representation of (index - 1), using 1-based indexing.
%
% Inputs:
%   indices : vector of integers in the range 1..d^n
%   n       : length of the code (number of digits)
%   d       : base / alphabet size
%
% Outputs:
%   codes   : N x n matrix, where each row is the base-d code associated
%             with the corresponding entry in indices

    % Number of indices to be converted
    N = numel(indices);

    % Preallocate output array
    codes = zeros(N, n);

    % Extract each digit of the base-d representation (most significant first)
    for k = 1:n
        pow = d^(n-k);
        codes(:,k) = mod(floor((indices-1)/pow), d) + 1;
    end
end
