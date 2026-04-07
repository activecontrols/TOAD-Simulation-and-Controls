%% Realistic Monte Carlo Disturbances (Deltas Only)
num_sims = 100;

% =======================================================
% 1. Moment of Inertia Disturbances (Delta J)
% =======================================================
% We define the nominals here ONLY to calculate the 10% proportional 
% variance for the diagonals. They are not added to the final output.
J_nom = constantsTOAD.J;

% Preallocate cell array
J_d_vals = cell(1, num_sims);

for i = 1:num_sims
    % 1. Diagonals: Base MATLAB randn() scaled by 10% of nominal
    dI_xx = (0.10 * J_nom(1,1)) * randn();
    dI_yy = (0.10 * J_nom(2,2)) * randn();
    dI_zz = (0.10 * J_nom(3,3)) * randn();
    
    % 2. Off-Diagonals: Base MATLAB randn() scaled by 1.0 variance
    dI_xy = 1.0 * randn();
    dI_xz = 1.0 * randn();
    dI_yz = 1.0 * randn();
    
    % 3. Assemble the Delta matrix (enforcing symmetry)
    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];
end

% =======================================================
% 2. Lever Arm Disturbances (Delta Lever Arm)
% =======================================================
% Define uncertainty (Standard Deviation) for each axis
sigma_lever = [0.05; 0.02; 0.02]; 

% Preallocate
TB_d_vals = cell(1, num_sims);

for i = 1:num_sims
    % Generate zero-mean 3D noise using randn(3,1) and scale by axis uncertainties.
    % This outputs ONLY the disturbance vector [dX; dY; dZ]
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;
end