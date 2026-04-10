function PlotMonteCarlo(filename)
    % PlotMonteCarlo Generates robustness and performance plots for MC runs.
    
    % Close all previous figures
    close all

    % Set background color
    darkMode = 1;
    if darkMode
        bkgColor = 'k';
    else
        bkgColor = 'w';
    end

    %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('No file provided. Pulling data from the base workspace...');
        
        reqVars = {'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                   'J_Wobble_Coup', 'J_Trans_Coup', 'RMSE_Controls_all', 'RMSE_Filter_all', 'GyroNoisePower_vals'};
        
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
                       'J_Wobble_Coup', 'J_Trans_Coup', 'RMSE_Controls_all', 'RMSE_Filter_all', 'GyroNoisePower_vals');
    end

    %% 2. Setup Plotting Parameters
    disp('Generating plots...');
    num_bins = 50; 
    clipUpper = @(x) min(x, prctile(x, 99)); % Clip top 1% for cleaner histograms

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

    % 2. Lateral Pos vs Wobble Coupling
    nexttile; grid on; hold on;
    scatter(J_Wobble_Coup, lat_pos_err, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Wobble Coupling (XZ/YZ)'); ylabel('Lateral Pos RMSE (m)');
    title('Pos vs. Inertia Wobble');

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

    %% Plot 5: Attituide with gimbal drift
    figure('Name', 'Attitude Err w/ Gyro Drift', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    gyroVals = cell2mat(GyroNoisePower_vals);

    % 1. Lateral Pos vs Lever Radial
    nexttile; grid on; hold on;
    scatter(gyroVals, att_pitch, 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Pitch RMSE (rad)');
    title('Pitch vs Gyro Drift');

    % 2. Lateral Pos vs Wobble Coupling
    nexttile; grid on; hold on;
    scatter(gyroVals, att_yaw, 30, [0 0.7 0], 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Yaw RMSE (rad)');
    title('Yaw vs Gyro Drift');

    % 3. Roll Error vs Axial Scale
    nexttile; grid on; hold on;
    scatter(gyroVals, att_roll, 30, 'b', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7);
    xlabel('Gyro Drift Power'); ylabel('Roll RMSE (rad)');
    title('Roll vs. Gyro Drift');

end