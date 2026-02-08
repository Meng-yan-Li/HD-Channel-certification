function B = commu_channel(A, v, d, channel)
% COMMU_CHANNEL
% Apply a noisy communication channel (or its effective action) to an operator.
%
% This function models the effect of a simple noise process on a d-dimensional
% operator A. Depending on the specified channel type, the map corresponds to:
%   - Depolarizing channel
%   - Dephasing channel
%
% The noise strength is controlled by the parameter v ∈ [0,1].
%
% Inputs:
%   A        : d x d operator (density matrix or observable)
%   v        : noise parameter (e.g. visibility)
%   d        : system dimension
%   channel  : string specifying the channel type
%              'depo'  - depolarizing channel
%              otherwise - dephasing channel
%
% Outputs:
%   B        : d x d operator after the action of the noisy channel

    if isequal(channel,'depo')  % Depolarizing channel
        % B = v A + (1 - v) Tr(A) I / d
        B = v * A + (1 - v) * trace(A) / d * eye(d);  
    else                        % Dephasing channel
        % B = v A + (1 - v) Diag(A)
        B = v * A + (1 - v) * diag(diag(A));
    end
end
