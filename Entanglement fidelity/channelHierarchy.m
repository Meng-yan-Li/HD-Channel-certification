function [Entanglement_fidelity, Pvar, IZvar] = channelHierarchy(ASP, VPHI, VIZ, icPHI, icIZ, d, n)
% CHANNELHIERARCHY
% Semidefinite programming (SDP) hierarchy for bounding the entanglement
% fidelity of a quantum channel from observed ASP values.
%
% This function formulates and solves a CVX-based SDP corresponding to a
% completely positive (CP) map constraint, together with linear relations
% induced by the prepare-and-measure scenario and moment-matrix positivity.
%
% Inputs:
%   ASP    : observed average success probability
%   VPHI   : cell array encoding Φ-type operator monomials and constraints
%   VIZ    : cell array encoding I⊗Z-type operator monomials and constraints
%   icPHI  : index list specifying the principal Φ moment submatrix
%   icIZ   : index list specifying the principal IZ moment submatrix
%   d      : local Hilbert space dimension
%   n      : number of preparation settings
%
% Outputs:
%   Entanglement_fidelity : optimal value of the SDP (entanglement fidelity bound)
%   Pvar                 : moment variables associated with Φ operators
%   IZvar                : moment variables associated with I⊗Z operators

    %% Normalization factors and problem size
    scale = d / (n * d^n);
    numPHIvar = length(VPHI);
    numIZvar  = length(VIZ);
    klistLength = sqrt(length(icPHI));

    %% CVX: Completely positive (CP) map relaxation
    cvx_begin SDP quiet
        % Decision variables: moment vectors for Φ and I⊗Z operator families
        variable Pvar(numPHIvar,1) complex   % χ[Φ]
        variable IZvar(numIZvar,1) complex   % χ[I⊗Z]

        % Objective: minimize entanglement infidelity (normalized)
        minimize(real(IZvar(1)) / (d^2))

        subject to
            %% Trace normalization
            Pvar(1) == 1;    % Tr(Φ) = 1

            %% Linear constraints induced by operator identities (Φ-sector)
            indices = zeros(1, numPHIvar);
            for i = 1:numPHIvar
                if isempty(VPHI{i})
                    Pvar(i) == 0;   % Orthogonality constraint

                elseif isscalar(VPHI{i}{1})
                    if VPHI{i}{1} > 0
                        % Classical input-output consistency constraints
                        if isscalar(VPHI{i}{2}) && VPHI{i}{2} == 0
                            Pvar(i) == 1/d;   % Tr(Φ · ρ ⊗ I) = 1/d

                        elseif isscalar(VPHI{i}{2}) && VPHI{i}{2} ~= 0
                            % Contribution to the ASP linear constraint
                            x_vec = mod(floor((VPHI{i}{1}-1) ./ d.^(0:n-1)), d) + 1;
                            y = floor((VPHI{i}{2}-1)/(d-1)) + 1;
                            b = mod(VPHI{i}{2}-1, d-1) + 1;

                            if x_vec(y) == b
                                indices(i) = indices(i) + 1;
                            elseif x_vec(y) == d
                                indices(i) = indices(i) - 1;
                            end
                        end

                    elseif VPHI{i}{1} < 0
                        % Reduced-state constraints
                        if isscalar(VPHI{i}{2}) && VPHI{i}{2} == 0
                            E = mod(floor((-VPHI{i}{1}-1) ./ d.^(0:1)), d) + 1;
                            if E(1) == E(2)
                                Pvar(i) == 1/d;   % Tr_B(ρ) = I/d
                            else
                                Pvar(i) == 0;
                            end
                        end
                    end
                end
            end

            % Observed ASP constraint
            indices * Pvar == ASP/scale - n*d^(n-2);

            %% Linear constraints for I⊗Z sector
            for i = 1:numIZvar
                if isempty(VIZ{i})
                    IZvar(i) == 0;   % Orthogonality constraint

                elseif isscalar(VIZ{i}{1})
                    if VIZ{i}{1} < 0
                        E = mod(floor((-VIZ{i}{1}-1) ./ d.^(0:1)), d) + 1;
                        if E(1) ~= E(2)
                            IZvar(i) == 0;
                        end

                    elseif VIZ{i}{1} == 0
                        % Consistency between identity and reduced operators
                        for j = 1:numIZvar
                            if ~isempty(VIZ{j}) && ...
                               isscalar(VIZ{j}{1}) && VIZ{j}{1} > 0 && ...
                               isequal(VIZ{i}{2}, VIZ{j}{2})
                                IZvar(i) == d * IZvar(j);   % Tr(I⊗M)

                            end
                            if ~isempty(VIZ{j}) && ...
                               isscalar(VIZ{j}{1}) && VIZ{j}{1} < 0 && ...
                               isequal(VIZ{i}{2}, VIZ{j}{2})
                                E = mod(floor((-VIZ{j}{1}-1) ./ d.^(0:1)), d) + 1;
                                if E(1) == E(2)
                                    IZvar(i) == d * IZvar(j);
                                end
                            end
                        end
                    end
                end
            end

            %% Moment-matrix positivity constraints
            reshape(Pvar(icPHI), [klistLength, klistLength]) ...
                == hermitian_semidefinite(klistLength);

            reshape(IZvar(icIZ), [klistLength, klistLength]) ...
              - reshape(Pvar(icPHI), [klistLength, klistLength]) ...
                == hermitian_semidefinite(klistLength);

    cvx_end

    %% Output: entanglement fidelity bound
    Entanglement_fidelity = cvx_optval;
end
