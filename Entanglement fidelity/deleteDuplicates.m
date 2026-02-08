function [list_unique] = deleteDuplicates(list)
% DELETEDUPLICATES
% Remove duplicate operator entries from a symbolic operator list while
% preserving the original order.
%
% Each row of `list` represents one composite operator encoded as a pair
% {v1, v2}, where v1 and v2 are vectors describing Alice- and Bob-side
% operator strings, respectively.
%
% Two operators are considered identical if both v1 and v2 are exactly
% equal. This function converts each operator pair into a unique string
% representation and removes duplicate rows accordingly.
%
% Input:
%   list        : cell array (N x 2), each row is {v1, v2}
%
% Output:
%   list_unique : cell array with duplicate operators removed (stable order)

    % Convert each operator pair {v1, v2} into a unique string key
    % The combination "mat2str(v1) _ mat2str(v2)" uniquely labels the operator
    cstr = cellfun(@(a,b) [mat2str(a) '_' mat2str(b)], ...
                   list(:,1), list(:,2), 'UniformOutput', false);

    % Identify unique entries while preserving the original order
    [~, ia] = unique(cstr, 'stable');

    % Extract the unique operators
    list_unique = list(ia,:);
end
