function PlotMonteCarlo(filename)
    % PlotMonteCarlo Generates robustness and performance plots for MC runs.
    %
    % Usage:
    %   PlotMonteCarlo()          - Uses data currently in the base workspace.
    %   PlotMonteCarlo('file.mat')- Loads data from the specified .mat file.

    %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('No file provided. Pulling data from the base workspace...');
        
        reqVars = {'J_mag', 'Lever_mag', 'RMSE_Controls_all', 'RMSE_Filter_all', 'metric_Filter'};
        
        for i = 1:length(reqVars)
            varName = reqVars{i};
            try
                eval([varName ' = evalin(''base'', ''' varName ''');']);
            catch
                error('Variable "%s" not found in base workspace. Run the simulation first or provide a valid .mat filename.', varName);
            end
        end
    else
        disp(['Loading data from: ', filename]);
        if ~isfile(filename)
            error('File "%s" does not exist. Please check the path and try again.', filename);
        end
        load(filename, 'J_mag', 'Lever_mag', 'RMSE_Controls_all', 'RMSE_Filter_all', 'metric_Filter');
    end

    %% 2. Setup Plotting Parameters
    disp('Generating plots...');
    num_bins = 50; 
    
    % RMSE is strictly positive. We only need to clip extreme upper outliers.
    clipUpper = @(x) min(x, prctile(x, 99));

    %% Plot 1: 3D Bivariate Histogram (Disturbance Density)
    figure('Name', 'MC Disturbance Density', 'Color', 'w', 'WindowStyle', 'docked');
    histogram2(J_mag, Lever_mag, 'DisplayStyle', 'bar3', 'FaceColor', 'flat');
    colorbar;
    xlabel('MoI Disturbance (Norm)', 'FontWeight', 'bold');
    ylabel('Lever Arm Disturbance (Norm)', 'FontWeight', 'bold');
    zlabel('Number of MC Runs', 'FontWeight', 'bold');
    title('3D Histogram: Distribution of Applied Disturbances');
    view(45, 45); grid on;

    %% Plot 2: Controller Performance Histograms (Stacked Subcomponents)
    figure('Name', 'Controller RMSE Distributions', 'Color', 'w', 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % --- Attitude Histogram (Yaw, Pitch, Roll) ---
    nexttile; hold on; grid on;
    att_roll  = RMSE_Controls_all(1, :);
    att_pitch = RMSE_Controls_all(2, :);
    att_yaw   = RMSE_Controls_all(3, :);

    % Plotted in order: Yaw, Pitch, Roll
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

    %% Plot 3: 3D Scatter (Lateral Position Error Baselined to Median)
    figure('Name', 'Robustness Scatter: Lateral Position Error', 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on;

    % Calculate Lateral Position Error (X and Y components only)
    lat_pos_err = sqrt(RMSE_Controls_all(4, :).^2 + RMSE_Controls_all(5, :).^2);
    
    Err_base_lat = median(lat_pos_err, 'omitnan');
    metric_Lat_c = lat_pos_err - Err_base_lat;

    p98_lat = prctile(metric_Lat_c, 98);
    p01_lat = prctile(metric_Lat_c, 1);
    idx_lat = (metric_Lat_c <= p98_lat) & (metric_Lat_c >= p01_lat);

    scatter3(J_mag(idx_lat), Lever_mag(idx_lat), metric_Lat_c(idx_lat), 75, 'b', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('MoI Disturbance (Norm)');
    ylabel('Lever Arm Disturbance (Norm)');
    zlabel('\Delta Lateral Pos RMSE (m) [From Median]');
    title('Robustness: Lateral Position Error (Extreme Outliers Omitted)');
    view(45, 30);

    %% Plot 4: Estimator Filter Performance Histogram (Overlaid Components)
    figure('Name', 'Estimator RMSE Distribution', 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on;

    % Extract individual raw RPY errors
    filt_yaw   = RMSE_Filter_all(1, :);
    filt_pitch = RMSE_Filter_all(2, :);
    filt_roll  = RMSE_Filter_all(3, :);

    histogram(clipUpper(filt_yaw), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Yaw Error');
    histogram(clipUpper(filt_pitch), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Pitch Error');
    histogram(clipUpper(filt_roll), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Roll Error');

    xlabel('Filter RMSE (degrees)'); 
    ylabel('Frequency'); 
    title('Estimator Component RMSE (Outliers Clipped)');
    legend('show', 'Location', 'northeast');

    %% Plot 5: 3D Scatter (Filter Error Baselined to Median)
    figure('Name', 'Robustness Scatter: Estimator Filter', 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on;

    Err_base_filt = median(metric_Filter, 'omitnan');
    metric_Filt_c = metric_Filter - Err_base_filt;

    p98_filt = prctile(metric_Filt_c, 98);
    p01_filt = prctile(metric_Filt_c, 1);
    idx_filt = (metric_Filt_c <= p98_filt) & (metric_Filt_c >= p01_filt);

    scatter3(J_mag(idx_filt), Lever_mag(idx_filt), metric_Filt_c(idx_filt), 75, [126 47 142]/255, 'filled', 'MarkerEdgeColor', 'k');
    xlabel('MoI Disturbance (Norm)');
    ylabel('Lever Arm Disturbance (Norm)');
    zlabel('\Delta Filter Total RMSE (degrees) [From Median]');
    title('Estimator Robustness (Extreme Outliers Omitted)');
    view(45, 30);
    
    disp('Plot generation complete.');
end