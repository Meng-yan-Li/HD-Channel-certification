%% Entanglement fidelity
% Main script for computing the entanglement fidelity bound using
% the channel hierarchy SDP relaxation.

clear
clc
tic

%% ====================== Basic parameters ======================
n = 2;   % Number of encoding states (number of preparations)
d = 2;   % Local Hilbert space dimension
k = 3;   % Hierarchy level in the SDP relaxation
         % Note: k = 1.5 corresponds to the mixed level k = 1 + AB

%% ================= Step 1: Operator & moment matrix generation =================
% Generate the symbolic operator list up to hierarchy level k
kOperatorList = generateOperatorList(n, d, k);

% Construct moment matrices Γ^Φ and Γ^IZ
[GammaPHI, GammaIZ, GammaPHI_char, GammaIZ_char] = ...
    generateMomentMatrix(kOperatorList, d);

%% ================= Step 2: Extract unique operator variables =================
% Identify unique moment variables appearing in Γ^Φ
[~, iaPHI, icPHI] = unique(GammaPHI_char(:), 'stable');
VPHI = GammaPHI(iaPHI);

% Identify unique moment variables appearing in Γ^IZ
[~, iaIZ, icIZ] = unique(GammaIZ_char(:), 'stable');
VIZ = GammaIZ(iaIZ);

%% ================= Step 3: Channel hierarchy SDP =================
% Classical and quantum bounds for the ASP witness
classical_bound = 1/2 * (1 + 1/d);
quantum_bound   = 1/2 * (1 + 1/sqrt(d));

% Sweep the ASP value within the physically relevant range
ASP_range = [classical_bound : 0.001 : quantum_bound, quantum_bound];

% Solve the channel hierarchy SDP for each ASP value
Entanglement_fidelity = arrayfun( ...
    @(ASP) channelHierarchy(ASP, VPHI, VIZ, icPHI, icIZ, d, n), ...
    ASP_range);

toc


% [Entanglement_fidelity,~,~] = channelHierarchy(0.8,VPHI,VIZ,icPHI,icIZ,d,n);
