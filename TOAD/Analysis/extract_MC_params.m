function extract_MC_params(seed, runNum, constantsTOAD)
rng(seed)
J_nom = constantsTOAD.J;

for i = 1:runNum

        %% 1. Moment of Inertia Disturbances
    dI_xx = (0.1 * J_nom(1,1)) * rand();
    dI_yy = (0.1 * J_nom(2,2)) * rand();
    dI_zz = (0.1 * J_nom(3,3)) * rand();

    % Keep products diagonal to match the V3 disturbance model.
    dI_xy = 0;
    dI_xz = 0;
    dI_yz = 0;

    J_d_vals = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];

    %% 2. Center-of-Mass / Lever-Arm Disturbance
    % V3 used TB_d as the 3-axis lever-arm / off-center-CG disturbance.
    % Preserve that model-facing variable and keep its radial/axial
    % components available for sensitivity analysis.
    sigma_lever = [0.01; 0.01; 0.01];
    TB_d_vals = randn(3, 1) .* sigma_lever;

    %% 3. Wind Disturbances
    % V3 wind gain distribution and covariance range are retained.
    a = 0.35;
    b = 2;
    Wind_Gain_vals = a * (-log(rand(1,1))).^(1/b);
    Wind_Covar_vals = 8 * rand();
end

disp("Disturbances in Inertia: ")
disp(J_d_vals)
disp("Disturbances in TB: ")
disp(TB_d_vals)
disp("Wind gain: ")
disp(Wind_Gain_vals)
disp("Wind Covar: ")
disp(Wind_Covar_vals)
end