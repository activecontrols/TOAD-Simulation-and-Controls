function PlotMonteCarlo(filename)
    % PlotMonteCarlo Generates robustness and performance plots for MC runs.
    
    % Close all previous figures
    close all

    % Set background color
    darkMode = 0;
    if darkMode
        bkgColor = 'k';
    else
        bkgColor = 'w';
    end

   %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('No file provided. Pulling data from the base workspace...');
        
        reqVars = {'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                   'RMSE_Controls_all', 'RMSE_Filter_all', ...
                   'GyroNoisePower_vals', 'Kg2_vals', 'G_RMAX_vals', ...
                   'kGrom_vals', 'bGrom_vals', 'G'}; % Added G here
        
        for i = 1:length(reqVars)
            varName = reqVars{i};
            try
                eval([varName ' = evalin(''base'', ''' varName ''');']);
            catch
                error('Variable "%s" not found in base workspace. Run the simulation first.', varName);
            end
        end
    else
        load(filename, 'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                       'RMSE_Controls_all', 'RMSE_Filter_all', ...
                       'GyroNoisePower_vals', 'Kg2_vals', 'G_RMAX_vals', ...
                       'kGrom_vals', 'bGrom_vals', 'G'); % Added G here
    end

    %% 2. Setup Plotting Parameters
    disp('Generating plots...');
    num_bins = 50; 
    clipUpper = @(x) min(x, prctile(x, 99)); % Clip top 1% for cleaner histograms

    % Safely extract arrays (in case they were preallocated as cells in the MC script)
    gyroVals = cell2mat(GyroNoisePower_vals);
    
    if iscell(Kg2_vals)
        kg2Vals = cell2mat(Kg2_vals);
    else
        kg2Vals = Kg2_vals;
    end
    
    if iscell(G_RMAX_vals)
        grmaxVals = cell2mat(G_RMAX_vals);
    else
        grmaxVals = G_RMAX_vals;
    end

    %% Plot 1: Bivariate Histograms (Disturbance Densities)
    figure('Name', 'MC Disturbance Densities', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    % Lever Arm Distributions (Stacked 1D)
    nexttile; hold on; grid on;
    histogram(Lever_Radial, num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Radial Shift');
    histogram(Lever_Axial, num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Axial Shift');
    xlabel('Shift (m)', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontWeight', 'bold');
    title('Lever Arm Disturbances');
    legend('show');

    % Inertia Scale Distributions (Stacked 1D)
    nexttile; hold on; grid on;
    histogram(J_Trans_Scale, num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Transverse Scale');
    histogram(J_Axial_Scale, num_bins, 'FaceColor', '#7E2F8E', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Axial Scale');
    xlabel('Scale Error', 'FontWeight', 'bold');
    ylabel('Frequency', 'FontWeight', 'bold');
    title('Inertia Scale Disturbances');
    legend('show');

    %% Plot 2: Controller Performance Histograms
    figure('Name', 'Controller RMSE Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % --- Attitude Histogram (Yaw, Pitch, Roll) ---
    nexttile; hold on; grid on;
    att_yaw   = RMSE_Controls_all(1, :);
    att_pitch = RMSE_Controls_all(2, :);
    att_roll  = RMSE_Controls_all(3, :);

    histogram(clipUpper(att_yaw), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Yaw');
    histogram(clipUpper(att_pitch), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Pitch');
    histogram(clipUpper(att_roll), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Roll');
    xlabel('Attitude RMSE (rad)'); ylabel('Frequency'); title('Attitude Error');
    legend('show', 'Location', 'northeast');

    % --- Lateral Position Histogram (X, Y Only) ---
    nexttile; hold on; grid on;
    pos_x = RMSE_Controls_all(4, :);
    pos_y = RMSE_Controls_all(5, :);

    histogram(clipUpper(pos_x), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'X');
    histogram(clipUpper(pos_y), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Y');
    xlabel('Lateral Position RMSE (m)'); title('Lateral Position Error');
    legend('show', 'Location', 'northeast');

    % --- Velocity Histogram ---
    nexttile; hold on; grid on;
    vel_x = RMSE_Controls_all(7, :);
    vel_y = RMSE_Controls_all(8, :);

    histogram(clipUpper(vel_x), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Vx');
    histogram(clipUpper(vel_y), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Vy');
    xlabel('Velocity RMSE (m/s)'); title('Velocity Error');
    legend('show', 'Location', 'northeast');

    %% Plot 3: Controller Sensitivities (Targeted Scatters)
    figure('Name', 'Targeted Controller Sensitivities', 'Color', bkgColor, 'WindowStyle', 'docked');
    
    % 'flow' automatically wraps tiles to fit the window size optimally
    tiledlayout('flow', 'TileSpacing', 'compact'); 
    
    lat_pos_err = sqrt(pos_x.^2 + pos_y.^2);
    att_trans_err = sqrt(att_pitch.^2 + att_roll.^2); % Pitch/Yaw comb

    % 1. Lateral Pos vs Lever Radial
    nexttile; grid on; hold on;
    scatter(Lever_Radial, lat_pos_err, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Lever Radial Shift (m)'); ylabel('Lateral Pos RMSE (m)');
    title('Pos vs. Off-Center CG');
    ylim([0 3]);

    % 2. Roll Error vs Axial Scale
    nexttile; grid on; hold on;
    scatter(J_Axial_Scale, att_roll, 30, [0 0.5 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Axial Scale Error (I_{zz})'); ylabel('Roll RMSE (rad)');
    title('Roll vs. Roll Inertia Error');

    % 3. Pitch/Yaw Error vs Transverse Scale (Ixx & Iyy)
    nexttile; grid on; hold on;
    scatter(J_Trans_Scale, att_trans_err, 30, [0.49 0.18 0.56], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Transverse Scale Error (I_{xx}, I_{yy})'); ylabel('Pitch/Yaw RMSE (rad)');
    title('Pitch/Yaw vs. Transverse Scale');

    % 4. Pitch/Yaw Error vs Lever Axial Shift
    nexttile; grid on; hold on;
    scatter(Lever_Axial, att_trans_err, 30, [0.30 0.75 0.93], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Lever Axial Shift (m)'); ylabel('Pitch/Yaw RMSE (rad)');
    title('Pitch/Yaw vs. CG Forward/Aft');

    %% Plot 4: Estimator Filter Performance Histogram
    figure('Name', 'Estimator RMSE Distribution', 'Color', bkgColor, 'WindowStyle', 'docked');
    hold on; grid on;

    filt_yaw   = RMSE_Filter_all(1, :);
    filt_pitch = RMSE_Filter_all(2, :);
    filt_roll  = RMSE_Filter_all(3, :);

    histogram(clipUpper(filt_yaw), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Yaw Error');
    histogram(clipUpper(filt_pitch), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Pitch Error');
    histogram(clipUpper(filt_roll), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Roll Error');
    
    xlabel('Filter RMSE (deg)'); ylabel('Frequency'); % Kept in degrees
    title('Estimator Component RMSE');
    legend('show', 'Location', 'northeast');

    disp('Plot generation complete.');

    %% Plot 5: Controller attitude with gimbal drift
    figure('Name', 'Controller attitude Err w/ Gyro Drift', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % 1. Yaw vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Yaw RMSE (rad)'); % Converted to Rad
    title('Yaw vs Gyro Drift');
    xscale log

    % 2. Pitch vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Pitch RMSE (rad)'); % Converted to Rad
    title('Pitch vs Gyro Drift');
    xscale log

    % 3. Roll Error vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Roll RMSE (rad)'); % Converted to Rad
    title('Roll vs. Gyro Drift');
    xscale log

    %% Plot 6: Filter attitude with gimbal drift combined
    figure('Name', 'Filter attitude Err w/ Gyro Drift combined', 'Color', bkgColor, 'WindowStyle', 'docked');

    scatter(gyroVals, filt_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7); hold on
    scatter(gyroVals, filt_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    scatter(gyroVals, filt_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    
    xlabel('Gyro Drift Power'); ylabel('Attitude RMSE (deg)'); % Kept in degrees
    title('Filter Attitude Error vs. Gyro Drift');
    legend('Yaw', 'Pitch', 'Roll', 'Location', 'best');
    xscale log
    grid on

    %% Plot 7: Bias & Noise Sampling Distributions
    figure('Name', 'Bias & Noise Sampling Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    
    % Gyro Drift Power Dist
    nexttile;
    histogram(gyroVals, num_bins, 'FaceColor', '#7E2F8E', 'FaceAlpha', 0.6);
    set(gca, 'XScale', 'log');
    xlabel('Gyro Drift Power'); ylabel('Frequency'); title('Gyro Noise Power Dist');
    grid on;

    % Kg2 Bias Dist
    nexttile;
    histogram(kg2Vals, num_bins, 'FaceColor', '#77AC30', 'FaceAlpha', 0.6);
    xlabel('Kg2 (g^2 bias)'); ylabel('Frequency'); title('Kg2 Bias Sensitivity Dist');
    grid on;
    
    % G_RMAX Dist
    nexttile;
    histogram(grmaxVals, num_bins, 'FaceColor', '#4DBEEE', 'FaceAlpha', 0.6);
    xlabel('G\_RMAX'); ylabel('Frequency'); title('G\_RMAX Distribution');
    grid on;

   %% Plot 8: Controller Lateral Pos Err w/ Kg2 Bias
    figure('Name', 'Controller Lat Pos Err w/ Kg2 Bias', 'Color', bkgColor, 'WindowStyle', 'docked');
    grid on; hold on;
    
    scatter(kg2Vals, pos_x, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'X Error');
    scatter(kg2Vals, pos_y, 30, 'filled', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Y Error');
    scatter(kg2Vals, lat_pos_err, 30, 'filled', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Magnitude');
    
    xlabel('Kg2 Bias Sensitivity'); ylabel('Lateral Pos RMSE (m)');
    title('Lateral Position Error vs Kg2 Bias (Log Scale)');
    legend('show', 'Location', 'best');

    %% Plot 9: Controller Lateral Pos Err w/ G_RMAX 
    figure('Name', 'Controller Lat Pos Err w/ G_RMAX', 'Color', bkgColor, 'WindowStyle', 'docked');
    grid on; hold on;
    
    scatter(grmaxVals, pos_x, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'X Error');
    scatter(grmaxVals, pos_y, 30, 'filled', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Y Error');
    scatter(grmaxVals, lat_pos_err, 30, 'filled', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Magnitude');
    
    xlabel('G\_RMAX'); ylabel('Lateral Pos RMSE (m)');
    title('Lateral Position Error vs G\_RMAX (Log Scale)');
    legend('show', 'Location', 'best');

    %% Plot 10: Controller Lateral Pos Err vs Filter Roll Error
    figure('Name', 'Controller Lat Pos Err vs Filter Roll', 'Color', bkgColor, 'WindowStyle', 'docked');
    grid on; hold on;
    
    scatter(filt_roll, pos_x, 30, 'filled', 'MarkerFaceColor', '#0072BD', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'X Error');
    scatter(filt_roll, pos_y, 30, 'filled', 'MarkerFaceColor', '#D95319', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Y Error');
    scatter(filt_roll, lat_pos_err, 30, 'filled', 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'DisplayName', 'Magnitude');
    
    set(gca, 'YScale', 'log'); % Set Y-axis to logarithmic
    max_y_val = max([max(pos_x(:)), max(pos_y(:)), max(lat_pos_err(:))]);
    ylim([prctile(lat_pos_err, 0.02), max_y_val * 1.2]); 
    xlabel('Filter Roll RMSE (deg)'); ylabel('Lateral Pos RMSE (m)');
    title('Lateral Position Error vs Filter Roll Error (Log Scale)');
    legend('show', 'Location', 'best');

    %% Plot 11: Stability Manifold (Surface)
    figure('Name', 'Stability Manifold: Kg2 vs G_RMAX', 'Color', bkgColor, 'WindowStyle', 'docked');
    
    % 1. Create a uniform grid for interpolation
    [KG, GR] = meshgrid(linspace(min(kg2Vals), max(kg2Vals), 100), ...
                        linspace(min(grmaxVals), max(grmaxVals), 100));
    
    % 2. Interpolate the scattered MC data onto the grid
    % We use 'natural' or 'linear' interpolation for the manifold
    F = scatteredInterpolant(kg2Vals(:), grmaxVals(:), lat_pos_err(:), 'linear', 'none');
    Z = F(KG, GR);
    
    % 3. Plot the surface
    surf(KG, GR, Z, 'EdgeColor', 'none');
    hold on;
    % Optional: Overlay the original MC points to see where data is dense
    scatter3(kg2Vals, grmaxVals, lat_pos_err, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
    
    view(2); % 3D View
    colormap turbo; colorbar;
    lowerLim = prctile(lat_pos_err, 15);
    lowerLim = min(lowerLim, 2.9);
    clim([lowerLim 3]);
    xlabel('Kg2 Bias Sensitivity'); ylabel('G\_RMAX'); zlabel('Lat Pos RMSE (m)');
    title('Performance Manifold & Stability Boundary');
    set(gca, 'ZScale', 'log'); % Often useful to see stability transitions

   %% Plot 12: Grommet Transmissibility with Confidence Bounds
    figure('Name', 'Grommet Transmissibility', 'Color', bkgColor, 'WindowStyle', 'docked');
    hold on; grid on;

    % 1. Setup Frequency Vector & Nominal Mass (from VRECalc context)
    f = logspace(0, log10(3000), 1000); % 1 Hz to 3000 Hz
    mBoard = 0.1; 

    % 2. Extract stiffness (k) and damping ratio (b) from the MC arrays
    k_vals = cell2mat(kGrom_vals);
    b_vals = cell2mat(bGrom_vals);
    
    % Calculate the true nominals from the generated distributions
    K_nom = mean(k_vals);
    B_nom = mean(b_vals);

    % 3. Preallocate and compute Transmissibility for all MC runs
    T_all = zeros(length(k_vals), length(f));
    for i = 1:length(k_vals)
        f_n_iter = 1 / (2 * pi) * sqrt(k_vals(i) / mBoard);
        r = f ./ f_n_iter;
        % 1-DoF Transmissibility Equation
        T_all(i, :) = sqrt((1 + (2 * b_vals(i) * r).^2) ./ ((1 - r.^2).^2 + (2 * b_vals(i) * r).^2));
    end

    % 4. Compute the Nominal Transmissibility Curve
    f_n_nom = 1 / (2 * pi) * sqrt(K_nom / mBoard);
    r_nom = f ./ f_n_nom;
    T_nom = sqrt((1 + (2 * B_nom * r_nom).^2) ./ ((1 - r_nom.^2).^2 + (2 * B_nom * r_nom).^2));

    % 5. Determine Confidence Region (1st and 99th percentiles)
    T_lower = prctile(T_all, 1, 1);
    T_upper = prctile(T_all, 99, 1);

    % 6. Plotting
    % Create combined x and y arrays to draw the shaded polygon 
    f_fill = [f, fliplr(f)];
    T_fill = [T_upper, fliplr(T_lower)];

    % Draw the shaded confidence region first using an RGB triplet
    fill(f_fill, T_fill, [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '98% Confidence Bound');
    
    % Draw the nominal curve
    plot(f, T_nom, 'k-', 'LineWidth', 1.5, 'DisplayName', 'Nominal Transmissibility');

    % At f = f_n_nom, the frequency ratio r is 1
    T_fn = sqrt((1 + (2 * B_nom * 1)^2) / ((1 - 1^2)^2 + (2 * B_nom * 1)^2));
    
    % Plot the point and add a text label next to it
    plot(f_n_nom, T_fn, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r', 'DisplayName', 'Natural Frequency (fn)');
    text(f_n_nom * 1.15, T_fn, sprintf('f_n = %.1f Hz', f_n_nom), 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    % 7. Formatting
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('Frequency (Hz)', 'FontWeight', 'bold');
    ylabel('Transmissibility (T)', 'FontWeight', 'bold');
    
    % Format dynamic title using G struct
    % Safely convert Name to string and force numbers to double to prevent ASCII recycling
    gName = string(G.Name); 
    duro_str = string(G.Durometer); % Keep as string for things like "OO"
    k_nom = str2double(string(G.K));
    c_nom = str2double(string(G.C));

    plotTitle = sprintf('Transmissibility Uncertainty: %s\n(Durometer: %s | K: %g | C: %g)', ...
                         gName, duro_str, k_nom, c_nom);
    title(plotTitle, 'Interpreter', 'none'); 
    
    legend('show', 'Location', 'southwest');
    
    % Clean up axes limits
    xlim([min(f) max(f)]);
    ylim([1e-2, max(T_upper)*2]);

    %% Plot 13: Controller Roll Error Manifold (Surface)
    figure('Name', 'Roll Error Manifold: Kg2 vs G_RMAX', 'Color', bkgColor, 'WindowStyle', 'docked');
    
    % 1. Create a uniform grid for interpolation (reusing previous grid logic)
    [KG, GR] = meshgrid(linspace(min(kg2Vals), max(kg2Vals), 100), ...
                        linspace(min(grmaxVals), max(grmaxVals), 100));
    
    % 2. Interpolate the scattered MC data onto the grid
    % Using att_roll (Controller Roll RMSE) calculated in Plot 2
    F_roll = scatteredInterpolant(kg2Vals(:), grmaxVals(:), filt_roll(:), 'linear', 'none');
    Z_roll = F_roll(KG, GR);
    
    % 3. Plot the surface
    surf(KG, GR, Z_roll, 'EdgeColor', 'none');
    hold on;
    
    % Overlay the original MC points
    scatter3(kg2Vals, grmaxVals, filt_roll, 10, 'k', 'filled', 'MarkerFaceAlpha', 0.2);
    
    % 4. Formatting
    view(2); % 2D Top-down view to match Plot 11
    colormap turbo; 
    colorbar;
    
    % Dynamically set color limits to prevent outliers from washing out the surface
    lowerLim_roll = prctile(filt_roll, 25);
    upperLim_roll = 16;
    if lowerLim_roll < upperLim_roll
        clim([lowerLim_roll upperLim_roll]);
    end
    
    xlabel('Kg2 Bias Sensitivity', 'FontWeight', 'bold'); 
    ylabel('G\_RMAX', 'FontWeight', 'bold'); 
    zlabel('MEKF Roll RMSE (deg)', 'FontWeight', 'bold');
    title('MEKF Roll Error Manifold');
    
    % Log scale for Z to easily spot stability transitions
    set(gca, 'ZScale', 'log');
end