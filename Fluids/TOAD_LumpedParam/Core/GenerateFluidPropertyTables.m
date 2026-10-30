function GenerateFluidPropertyTables(outDir, N_pts)
% GENERATEFLUIDPROPERTYTABLES Generates high-density (u, rho) thermodynamic
% property tables for Nitrogen and Oxygen using CoolProp and exports them as CSV.
%
% Inputs:
%   outDir - Output directory for CSV files (default: 'sandbox/experiments/Fluid Properties')
%   N_pts  - Number of grid points per dimension (default: 150)

    if nargin < 1 || isempty(outDir)
        outDir = fullfile(pwd, 'sandbox', 'experiments', 'Fluid Properties');
    end
    if nargin < 2 || isempty(N_pts)
        N_pts = 150; % 150x150 = 22,500 grid points per fluid
    end

    if ~exist(outDir, 'dir')
        mkdir(outDir);
    end

    fprintf('====================================================\n');
    fprintf('  Generating High-Density Fluid Property Tables     \n');
    fprintf('  Grid Resolution: %dx%d (%d points per fluid)      \n', N_pts, N_pts, N_pts^2);
    fprintf('  Target Directory: %s\n', outDir);
    fprintf('====================================================\n');

    %% 1. Nitrogen (N2) Table Generation
    fprintf('\n[1/2] Computing Nitrogen Table...\n');
    tic;
    % Ranges: Covering 14.7 psi to 6000 psi, 70 K to 350 K
    rho_min_n2 = 0.5;   rho_max_n2 = 820.0;   % kg/m^3
    u_min_n2   = 60e3;  u_max_n2   = 350e3;   % J/kg

    rho_vec_n2 = linspace(rho_min_n2, rho_max_n2, N_pts);
    u_vec_n2   = linspace(u_min_n2, u_max_n2, N_pts);

    [U_grid_n2, RHO_grid_n2] = meshgrid(u_vec_n2, rho_vec_n2);
    P_grid_n2     = zeros(N_pts, N_pts);
    T_grid_n2     = zeros(N_pts, N_pts);
    H_grid_n2     = zeros(N_pts, N_pts);
    GAMMA_grid_n2 = zeros(N_pts, N_pts);
    MU_grid_n2    = zeros(N_pts, N_pts);
    K_grid_n2     = zeros(N_pts, N_pts);

    totalPts = N_pts * N_pts;
    validCount = 0;
    for r = 1:N_pts
        rho_val = rho_vec_n2(r);
        for c = 1:N_pts
            u_val = u_vec_n2(c);
            try
                p_val = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                t_val = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                h_val = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                g_val = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                try
                    mu_val = double(py.CoolProp.CoolProp.PropsSI('V', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                    k_val  = double(py.CoolProp.CoolProp.PropsSI('L', 'Umass', u_val, 'Dmass', rho_val, 'Nitrogen'));
                catch
                    mu_val = 1.8e-5;
                    k_val  = 0.026;
                end
                validCount = validCount + 1;
            catch
                % Fallback ideal gas near edge of EOS domain
                R_n2 = 296.8;
                t_val = max(u_val / 743.0, 70);
                p_val = rho_val * R_n2 * t_val;
                h_val = u_val + p_val / rho_val;
                g_val = 1.40;
                mu_val = 1.8e-5;
                k_val  = 0.026;
            end
            P_grid_n2(r, c)     = p_val;
            T_grid_n2(r, c)     = t_val;
            H_grid_n2(r, c)     = h_val;
            GAMMA_grid_n2(r, c) = g_val;
            MU_grid_n2(r, c)    = mu_val;
            K_grid_n2(r, c)     = k_val;
        end
        if mod(r, round(N_pts/5)) == 0
            fprintf('  Progress: %d%% (%d / %d rows)\n', round(100*r/N_pts), r, N_pts);
        end
    end
    t_n2 = toc;
    fprintf('  Nitrogen table complete in %.1f s (%d/%d valid EOS points)\n', t_n2, validCount, totalPts);

    % Save Nitrogen Table to CSV
    csv_n2 = fullfile(outDir, 'Nitrogen_Props.csv');
    n2_data = [U_grid_n2(:), RHO_grid_n2(:), P_grid_n2(:), T_grid_n2(:), H_grid_n2(:), GAMMA_grid_n2(:), MU_grid_n2(:), K_grid_n2(:)];
    T_n2_table = array2table(n2_data, 'VariableNames', {'u', 'rho', 'P', 'T', 'h', 'gamma', 'mu', 'k'});
    writetable(T_n2_table, csv_n2);
    fprintf('  Saved: %s\n', csv_n2);

    %% 2. Oxygen (O2) Table Generation
    fprintf('\n[2/2] Computing Oxygen Table...\n');
    tic;
    % Ranges: Subcooled LOX (70 K, 1141+ kg/m^3) to ambient gas (300 K)
    rho_min_ox = 0.5;    rho_max_ox = 1250.0;   % kg/m^3
    u_min_ox   = -170e3; u_max_ox   = 220e3;    % J/kg

    rho_vec_ox = linspace(rho_min_ox, rho_max_ox, N_pts);
    u_vec_ox   = linspace(u_min_ox, u_max_ox, N_pts);

    [U_grid_ox, RHO_grid_ox] = meshgrid(u_vec_ox, rho_vec_ox);
    P_grid_ox     = zeros(N_pts, N_pts);
    T_grid_ox     = zeros(N_pts, N_pts);
    H_grid_ox     = zeros(N_pts, N_pts);
    GAMMA_grid_ox = zeros(N_pts, N_pts);
    MU_grid_ox    = zeros(N_pts, N_pts);
    K_grid_ox     = zeros(N_pts, N_pts);

    validCount = 0;
    for r = 1:N_pts
        rho_val = rho_vec_ox(r);
        for c = 1:N_pts
            u_val = u_vec_ox(c);
            try
                p_val = double(py.CoolProp.CoolProp.PropsSI('P', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                t_val = double(py.CoolProp.CoolProp.PropsSI('T', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                h_val = double(py.CoolProp.CoolProp.PropsSI('Hmass', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                g_val = double(py.CoolProp.CoolProp.PropsSI('isentropic_expansion_coefficient', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                try
                    mu_val = double(py.CoolProp.CoolProp.PropsSI('V', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                    k_val  = double(py.CoolProp.CoolProp.PropsSI('L', 'Umass', u_val, 'Dmass', rho_val, 'Oxygen'));
                catch
                    mu_val = 1.9e-4;
                    k_val  = 0.15;
                end
                validCount = validCount + 1;
            catch
                % Fallback ideal gas / incompressible liquid
                if rho_val > 500
                    t_val = 90.0;
                    p_val = 101325;
                    h_val = u_val + p_val / rho_val;
                    g_val = 1.10;
                    mu_val = 1.9e-4;
                    k_val  = 0.15;
                else
                    R_ox = 259.8;
                    t_val = max(u_val / 657.0, 70);
                    p_val = rho_val * R_ox * t_val;
                    h_val = u_val + p_val / rho_val;
                    g_val = 1.40;
                    mu_val = 2.0e-5;
                    k_val  = 0.026;
                end
            end
            P_grid_ox(r, c)     = p_val;
            T_grid_ox(r, c)     = t_val;
            H_grid_ox(r, c)     = h_val;
            GAMMA_grid_ox(r, c) = g_val;
            MU_grid_ox(r, c)    = mu_val;
            K_grid_ox(r, c)     = k_val;
        end
        if mod(r, round(N_pts/5)) == 0
            fprintf('  Progress: %d%% (%d / %d rows)\n', round(100*r/N_pts), r, N_pts);
        end
    end
    t_ox = toc;
    fprintf('  Oxygen table complete in %.1f s (%d/%d valid EOS points)\n', t_ox, validCount, totalPts);

    % Save Oxygen Table to CSV
    csv_ox = fullfile(outDir, 'Oxygen_Props.csv');
    ox_data = [U_grid_ox(:), RHO_grid_ox(:), P_grid_ox(:), T_grid_ox(:), H_grid_ox(:), GAMMA_grid_ox(:), MU_grid_ox(:), K_grid_ox(:)];
    T_ox_table = array2table(ox_data, 'VariableNames', {'u', 'rho', 'P', 'T', 'h', 'gamma', 'mu', 'k'});
    writetable(T_ox_table, csv_ox);
    fprintf('  Saved: %s\n', csv_ox);

    % Also save a fast .mat cache for instant loading
    matFile = fullfile(outDir, 'FluidPropertyTables.mat');
    save(matFile, 'u_vec_n2', 'rho_vec_n2', 'P_grid_n2', 'T_grid_n2', 'H_grid_n2', 'GAMMA_grid_n2', 'MU_grid_n2', 'K_grid_n2', ...
                  'u_vec_ox', 'rho_vec_ox', 'P_grid_ox', 'T_grid_ox', 'H_grid_ox', 'GAMMA_grid_ox', 'MU_grid_ox', 'K_grid_ox', '-v7.3');
    fprintf('\nAll tables successfully generated and saved to:\n  CSV: %s\n  MAT: %s\n', outDir, matFile);
end
