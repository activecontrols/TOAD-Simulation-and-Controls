function PlotMonteCarlo(filename)
    % PlotMonteCarlo Generates robustness and performance plots for MC runs.
    %
    % Usage:
    %   PlotMonteCarlo()          - Uses data currently in the base workspace.
    %   PlotMonteCarlo('file.mat')- Loads data from the specified .mat file.

    %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('No file provided. Pulling data from the base workspace...');
        
        % Added SSE_Filter_all to extract the raw RPY components
        reqVars = {'J_mag', 'Lever_mag', 'metric_Ctrl_Att', 'metric_Ctrl_Pos', ...
                   'metric_Ctrl_Vel', 'metric_Filter', 'SSE_Filter_all'};
        
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
        load(filename, 'J_mag', 'Lever_mag', 'metric_Ctrl_Att', 'metric_Ctrl_Pos', 'metric_Ctrl_Vel', 'metric_Filter', 'SSE_Filter_all');
    end

    %% 2. Setup Plotting Parameters
    disp('Generating plots...');
    
    % Recalculate means
    mean_att = mean(metric_Ctrl_Att, 'omitnan');
    mean_pos = mean(metric_Ctrl_Pos, 'omitnan');
    mean_vel = mean(metric_Ctrl_Vel, 'omitnan');
    num_bins = 50; 
    
    % Histogram outlier clipping
    clipBoth  = @(x) max(min(x, prctile(x, 99)), prctile(x, 1));
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

    %% Plot 2: Controller Performance Histograms
    figure('Name', 'Controller SSE Distributions', 'Color', 'w', 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    % Attitude Histogram
    nexttile; hold on; grid on;
    histogram(clipBoth(metric_Ctrl_Att - mean_att), num_bins, 'FaceColor', '#0072BD');
    xline(0, 'k-', 'Mean', 'LineWidth', 1.5);
    xlabel('\Delta Attitude SSE (rad)'); ylabel('Frequency'); title('Attitude Dev (Outliers Clipped)');

    % Position Histogram
    nexttile; hold on; grid on;
    histogram(clipBoth(metric_Ctrl_Pos - mean_pos), num_bins, 'FaceColor', '#D95319');
    xline(0, 'k-', 'Mean', 'LineWidth', 1.5);
    xlabel('\Delta Position SSE (m)'); title('Position Dev (Outliers Clipped)');

    % Velocity Histogram
    nexttile; hold on; grid on;
    histogram(clipBoth(metric_Ctrl_Vel - mean_vel), num_bins, 'FaceColor', '#EDB120');
    xline(0, 'k-', 'Mean', 'LineWidth', 1.5);
    xlabel('\Delta Velocity SSE (m/s)'); title('Velocity Dev (Outliers Clipped)');

    %% Plot 3: 3D Scatter (Position Error Baselined to Median)
    figure('Name', 'Robustness Scatter: Position Error', 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on;

    Err_base_pos = median(metric_Ctrl_Pos, 'omitnan');
    metric_Pos_c = metric_Ctrl_Pos - Err_base_pos;

    p98_pos = prctile(metric_Pos_c, 98);
    p01_pos = prctile(metric_Pos_c, 1);
    idx_pos = (metric_Pos_c <= p98_pos) & (metric_Pos_c >= p01_pos);

    scatter3(J_mag(idx_pos), Lever_mag(idx_pos), metric_Pos_c(idx_pos), 75, 'b', 'filled', 'MarkerEdgeColor', 'k');
    xlabel('MoI Disturbance (Norm)');
    ylabel('Lever Arm Disturbance (Norm)');
    zlabel('\Delta Position SSE (m) [From Median]');
    title('Robustness: Position Error (Extreme Outliers Omitted)');
    view(45, 30);

    %% Plot 4: Estimator Filter Performance Histogram (Overlaid Components)
    figure('Name', 'Estimator SSE Distribution', 'Color', 'w', 'WindowStyle', 'docked');
    hold on; grid on;

    % Extract individual raw RPY errors
    filt_roll  = SSE_Filter_all(1, :);
    filt_pitch = SSE_Filter_all(2, :);
    filt_yaw   = SSE_Filter_all(3, :);

    % Plot histograms with transparency (FaceAlpha) so they overlay cleanly
    % Using clipBoth since component errors can drift negative or positive
    histogram(clipBoth(filt_roll), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Roll Error');
    histogram(clipBoth(filt_pitch), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Pitch Error');
    histogram(clipBoth(filt_yaw), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Yaw Error');

    xline(0, 'k-', 'Zero Error', 'LineWidth', 1.5, 'HandleVisibility', 'off');

    xlabel('Filter Error (degrees)'); 
    ylabel('Frequency'); 
    title('Estimator Component Errors (Outliers Clipped)');
    legend('show', 'Location', 'best');

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
    zlabel('\Delta Filter RPY Error (degrees) [From Median]');
    title('Estimator Robustness (Extreme Outliers Omitted)');
    view(45, 30);
    
    disp('Plot generation complete.');
end