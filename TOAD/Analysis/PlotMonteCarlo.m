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
        
        % Added Kg2_vals and G_RMAX_vals to requested variables
        reqVars = {'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                   'J_Wobble_Coup', 'J_Trans_Coup', 'RMSE_Controls_all', 'RMSE_Filter_all', ...
                   'GyroNoisePower_vals', 'Kg2_vals', 'G_RMAX_vals'};
        
        for i = 1:length(reqVars)
            varName = reqVars{i};
            try
                eval([varName ' = evalin(''base'', ''' varName ''');']);
            catch
                error('Variable "%s" not found in base workspace. Run the simulation first.', varName);
            end
        end
    else
        % Added Kg2_vals and G_RMAX_vals to file loader
        load(filename, 'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                       'J_Wobble_Coup', 'J_Trans_Coup', 'RMSE_Controls_all', 'RMSE_Filter_all', ...
                       'GyroNoisePower_vals', 'Kg2_vals', 'G_RMAX_vals');
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

    % Lever Arm Distributions
    nexttile;
    histogram2(Lever_Radial, Lever_Axial, 'DisplayStyle', 'bar3', 'FaceColor', 'flat');
    xlabel('Radial Shift (m)', 'FontWeight', 'bold');
    ylabel('Axial Shift (m)', 'FontWeight', 'bold');
    zlabel('MC Runs');
    title('Lever Arm Disturbances');
    view(45, 45); grid on;
    

    % Inertia Coupling Distributions
    nexttile;
    histogram2(J_Trans_Scale, J_Wobble_Coup, 'DisplayStyle', 'bar3', 'FaceColor', 'flat');
    xlabel('Transverse Scale Error', 'FontWeight', 'bold');
    ylabel('Wobble Coupling (XZ/YZ)', 'FontWeight', 'bold');
    zlabel('MC Runs');
    title('Inertia Disturbances');
    view(45, 45); grid on;

    %% Plot 2: Controller Performance Histograms
    figure('Name', 'Controller RMSE Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % --- Attitude Histogram (Yaw, Pitch, Roll) ---
    nexttile; hold on; grid on;
    att_yaw  = RMSE_Controls_all(1, :);
    att_pitch = RMSE_Controls_all(2, :);
    att_roll   = RMSE_Controls_all(3, :);

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
    
    lat_pos_err = sqrt(RMSE_Controls_all(4, :).^2 + RMSE_Controls_all(5, :).^2);
    att_trans_err = sqrt(RMSE_Controls_all(2, :).^2 + RMSE_Controls_all(3, :).^2); % Pitch/Yaw comb

    % 1. Lateral Pos vs Lever Radial
    nexttile; grid on; hold on;
    scatter(Lever_Radial, lat_pos_err, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Lever Radial Shift (m)'); ylabel('Lateral Pos RMSE (m)');
    title('Pos vs. Off-Center CG');
    ylim([0 3]);

    % 2. Lateral Pos vs Wobble Coupling
    nexttile; grid on; hold on;
    scatter(J_Wobble_Coup, lat_pos_err, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Wobble Coupling (XZ/YZ)'); ylabel('Lateral Pos RMSE (m)');
    title('Pos vs. Inertia Wobble');
    ylim([0 3]);

    % 3. Roll Error vs Axial Scale
    nexttile; grid on; hold on;
    scatter(J_Axial_Scale, att_roll, 30, [0 0.5 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Axial Scale Error (I_{zz})'); ylabel('Roll RMSE (rad)');
    title('Roll vs. Roll Inertia Error');

    % 4. Pitch/Yaw Error vs Transverse Coupling
    nexttile; grid on; hold on;
    scatter(J_Trans_Coup, att_trans_err, 30, [0.85 0.33 0.1], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Transverse Coupling (I_{xy})'); ylabel('Pitch/Yaw RMSE (rad)');
    title('Pitch/Yaw vs. Transverse Coupling');

    % 5. Pitch/Yaw Error vs Transverse Scale (Ixx & Iyy)
    nexttile; grid on; hold on;
    scatter(J_Trans_Scale, att_trans_err, 30, [0.49 0.18 0.56], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Transverse Scale Error (I_{xx}, I_{yy})'); ylabel('Pitch/Yaw RMSE (rad)');
    title('Pitch/Yaw vs. Transverse Scale');

    % 6. Pitch/Yaw Error vs Lever Axial Shift
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
    
    xlabel('Filter RMSE (rad)'); ylabel('Frequency'); 
    title('Estimator Component RMSE');
    legend('show', 'Location', 'northeast');

    disp('Plot generation complete.');

    %% Plot 5: Controller attitude with gimbal drift
    figure('Name', 'Controller attitude Err w/ Gyro Drift', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % 1. Yaw vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Yaw RMSE (deg)');
    title('Yaw vs Gyro Drift');
    xscale log

    % 2. Pitch vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Pitch RMSE (deg)');
    title('Pitch vs Gyro Drift');
    xscale log

    % 3. Roll Error vs Drift
    nexttile; grid on; hold on;
    scatter(gyroVals, att_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Roll RMSE (deg)');
    title('Roll vs. Gyro Drift');
    xscale log

    %% Plot 6: Filter attitude with gimbal drift combined
    figure('Name', 'Filter attitude Err w/ Gyro Drift combined', 'Color', bkgColor, 'WindowStyle', 'docked');

    scatter(gyroVals, filt_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7); hold on
    scatter(gyroVals, filt_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    scatter(gyroVals, filt_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    
    xlabel('Gyro Drift Power'); ylabel('Attitude RMSE (deg)');
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

    %% Plot 8: Controller attitude Err w/ Kg2 Bias
    figure('Name', 'Controller attitude Err w/ Kg2 Bias', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % 1. Yaw vs Kg2
    nexttile; grid on; hold on;
    scatter(kg2Vals, att_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Kg2 Bias Sensitivity'); ylabel('Yaw RMSE (deg)');
    title('Yaw vs Kg2 Bias');

    % 2. Pitch vs Kg2
    nexttile; grid on; hold on;
    scatter(kg2Vals, att_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Kg2 Bias Sensitivity'); ylabel('Pitch RMSE (deg)');
    title('Pitch vs Kg2 Bias');

    % 3. Roll vs Kg2
    nexttile; grid on; hold on;
    scatter(kg2Vals, att_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Kg2 Bias Sensitivity'); ylabel('Roll RMSE (deg)');
    title('Roll vs Kg2 Bias');

    %% Plot 9: Filter attitude Err w/ Kg2 Bias combined
    figure('Name', 'Filter attitude Err w/ Kg2 Bias combined', 'Color', bkgColor, 'WindowStyle', 'docked');

    scatter(kg2Vals, filt_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7); hold on
    scatter(kg2Vals, filt_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    scatter(kg2Vals, filt_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    
    xlabel('Kg2 Bias Sensitivity'); ylabel('Attitude RMSE (deg)');
    title('Filter Attitude Error vs Kg2 Bias');
    legend('Yaw', 'Pitch', 'Roll', 'Location', 'best');
    grid on

    %% Plot 10: Controller attitude Err w/ G_RMAX 
    figure('Name', 'Controller attitude Err w/ G_RMAX', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % 1. Yaw vs G_RMAX
    nexttile; grid on; hold on;
    scatter(grmaxVals, att_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('G\_RMAX'); ylabel('Yaw RMSE (deg)');
    title('Yaw vs G\_RMAX');

    % 2. Pitch vs G_RMAX
    nexttile; grid on; hold on;
    scatter(grmaxVals, att_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('G\_RMAX'); ylabel('Pitch RMSE (deg)');
    title('Pitch vs G\_RMAX');

    % 3. Roll vs G_RMAX
    nexttile; grid on; hold on;
    scatter(grmaxVals, att_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('G\_RMAX'); ylabel('Roll RMSE (deg)');
    title('Roll vs G\_RMAX');

    %% Plot 11: Filter attitude Err w/ G_RMAX combined
    figure('Name', 'Filter attitude Err w/ G_RMAX combined', 'Color', bkgColor, 'WindowStyle', 'docked');

    scatter(grmaxVals, filt_yaw, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7); hold on
    scatter(grmaxVals, filt_pitch, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    scatter(grmaxVals, filt_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    
    xlabel('G\_RMAX'); ylabel('Attitude RMSE (deg)');
    title('Filter Attitude Error vs G\_RMAX');
    legend('Yaw', 'Pitch', 'Roll', 'Location', 'best');
    grid on

end