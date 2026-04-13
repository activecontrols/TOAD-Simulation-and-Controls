function PlotMCTrajectories(filename)
    % PlotMCTrajectories Generates trajectory overlays and metric histograms
    
    close all;
    bkgColor = 'w';
    alphaVal = 0.2; % Transparency for MC runs

   %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('Pulling data from the base workspace...');
        reqVars = {'pos_all', 'vel_all', 'pos_avg', 't_common', 'pos_targ_ref', ...
                   'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                   'RMSE_Controls_all', 'RMSE_Filter_all', ...
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
        load(filename);
    end

    Waypoints = TrajectoryBuilder();
    numWP = numel(Waypoints);
    wp_coords = zeros(3, numWP);
    for idx = 1:numWP
        wp_coords(:, idx) = Waypoints(idx).Position;
    end
    
    % Success/Failure Classification
    final_wp = Waypoints(end).Position;
    final_tol = Waypoints(end).PosTol;
    is_success = false(size(pos_all, 1), 1);
    
    for i = 1:size(pos_all, 1)
        % Find the last valid (non-NaN) index for this trajectory
        last_valid = find(~isnan(pos_all(i,:,1)), 1, 'last');
        
        if ~isempty(last_valid)
            % 1. Final Position Criterion
            final_pos = [pos_all(i, last_valid, 1), pos_all(i, last_valid, 2), pos_all(i, last_valid, 3)];
            pos_met = all(abs(final_pos(:) - final_wp(:)) <= final_tol);
            
            % 2. Approach Criteria (Evaluated at 2m crossing)
            appr_met = false;
            alts = pos_all(i, 1:last_valid, 3);
            idx_above = find(alts > 2.0, 1, 'last'); % Find last time above 2m
            
            if ~isempty(idx_above)
                appr_idx = min(idx_above + 1, last_valid);
                appr_alt = max(alts(appr_idx), 0.01); % Prevent division by zero
                
                % Lateral Error & Glideslope Angle
                appr_pos_xy = [pos_all(i, appr_idx, 1), pos_all(i, appr_idx, 2)];
                lat_err = norm(appr_pos_xy - final_wp(1:2));
                glideslope_deg = atan2d(lat_err, appr_alt); % Angle relative to vertical
                
                % Lateral Velocity
                appr_vel_xy = [vel_all(i, appr_idx, 1), vel_all(i, appr_idx, 2)];
                lat_vel = norm(appr_vel_xy);
                
                % Check if approach conditions are within bounds
                if glideslope_deg <= 10.0 && lat_vel <= 0.3
                    appr_met = true;
                end
            end
            
            % Classify as success if ALL conditions are met
            if pos_met && appr_met
                is_success(i) = true;
            end
        end
    end
    
    for i = 1:size(pos_all, 1)
        % Find the last valid (non-NaN) index for this trajectory
        last_valid = find(~isnan(pos_all(i,:,1)), 1, 'last');
        if ~isempty(last_valid)
            final_pos = [pos_all(i, last_valid, 1), pos_all(i, last_valid, 2), pos_all(i, last_valid, 3)];
            % Classify as success if final position is within the bounding box of the final WP
            if all(abs(final_pos(:) - final_wp(:)) <= final_tol)
                is_success(i) = true;
            end
        end
    end
    
    success_count = sum(is_success);
    disp(['Successful Runs: ', num2str(success_count), ' / ', num2str(size(pos_all, 1))]);
    
    %% Setup Plotting Parameters
    disp('Generating plots...');
    num_bins = 30; 
    clipUpper = @(x) min(x, prctile(x, 99)); 

    gyroVals = cell2mat(GyroNoisePower_vals);
    
    % Corrected MATLAB syntax for conditional cell-to-mat conversion
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

    %% Plot 1: 3D Mission Profile Overlay
    figure('Name', 'MC Mission Trajectories', 'Color', bkgColor, 'WindowStyle', 'docked');
    tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    axMain = nexttile(tl, 1, [3 3]); 
    hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
    xlabel(axMain, 'North [m]'); ylabel(axMain, 'West [m]'); zlabel(axMain, 'Altitude [m]');
    title(axMain, sprintf('3D Trajectory Envelope (%d Runs)', size(pos_all, 1)));

    % Draw Waypoint boxes
    for idx = 1:numWP
        wp = Waypoints(idx);
        if ~wp.IsPassAndGo && wp.Position(3) >= 0
            tol = wp.PosTol;
            v_cube = 2 * tol * [0.5 -0.5 -0.5; 0.5 0.5 -0.5; -0.5 0.5 -0.5; -0.5 -0.5 -0.5; ...
                                0.5 -0.5 0.5; 0.5 0.5 0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5] + wp.Position'; 
            [f_cube, ~] = convhull(v_cube);
            patch(axMain, 'Vertices', v_cube, 'Faces', f_cube, 'FaceColor', 'g', ...
                  'FaceAlpha', 0.1, 'EdgeColor', '#4DBEEE');
        end
    end

   % Draw Runs Transparently (Green = Success, Red = Failure)
    num_sims = size(pos_all,1);
    for i = 1:num_sims
        if is_success(i)
            plotColor = [0, 1, 0, alphaVal]; % Green
        else
            plotColor = [1, 0, 0, alphaVal]; % Red
        end
        plot3(axMain, pos_all(i,:,1), pos_all(i,:,2), pos_all(i,:,3), ...
              'Color', plotColor, 'LineWidth', 0.5);
    end
    
    % Draw Average Trajectory
    plot3(axMain, pos_avg(:,1), pos_avg(:,2), pos_avg(:,3), 'k', 'LineWidth', 2.5);
    
    % Draw Target Waypoints
    plot3(axMain, wp_coords(1,:), wp_coords(2,:), wp_coords(3,:), 'r--s', ...
        'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

    % Orthographic projections
    axTop = nexttile(tl, 4);  drawMCTraj(axTop, pos_all, pos_avg, wp_coords, Waypoints, 2, 1, 'West [m]', 'North [m]', alphaVal, is_success); title(axTop, 'Top View');
    axSide = nexttile(tl, 8); drawMCTraj(axSide, pos_all, pos_avg, wp_coords, Waypoints, 1, 3, 'North [m]', 'Altitude [m]', alphaVal, is_success); title(axSide, 'Side View');
    axFront = nexttile(tl, 12); drawMCTraj(axFront, pos_all, pos_avg, wp_coords, Waypoints, 2, 3, 'West [m]', 'Altitude [m]', alphaVal, is_success); title(axFront, 'Front View');

    %% Plot 2: Disturbance Densities
    figure('Name', 'Disturbance Densities', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact');

    nexttile; hold on; grid on;
    histogram(Lever_Radial, num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Radial Shift');
    histogram(Lever_Axial, num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Axial Shift');
    xlabel('Shift (m)', 'FontWeight', 'bold'); title('Lever Arm Disturbances'); legend('show');

    nexttile; hold on; grid on;
    histogram(J_Trans_Scale, num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Transverse');
    histogram(J_Axial_Scale, num_bins, 'FaceColor', '#7E2F8E', 'FaceAlpha', 0.6, 'EdgeColor', 'none', 'DisplayName', 'Axial');
    xlabel('Scale Error', 'FontWeight', 'bold'); title('Inertia Scale Disturbances'); legend('show');

    %% Plot 3: Controller Performance Histograms
    figure('Name', 'Controller RMSE Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');

    nexttile; hold on; grid on;
    histogram(clipUpper(RMSE_Controls_all(1, :)), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'DisplayName', 'Yaw');
    histogram(clipUpper(RMSE_Controls_all(2, :)), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'DisplayName', 'Pitch');
    histogram(clipUpper(RMSE_Controls_all(3, :)), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'DisplayName', 'Roll');
    xlabel('Attitude RMSE (rad)'); title('Attitude Error'); legend('show');

    nexttile; hold on; grid on;
    histogram(clipUpper(RMSE_Controls_all(4, :)), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'DisplayName', 'X');
    histogram(clipUpper(RMSE_Controls_all(5, :)), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'DisplayName', 'Y');
    xlabel('Lateral Position RMSE (m)'); title('Lateral Error'); legend('show');

    nexttile; hold on; grid on;
    histogram(clipUpper(RMSE_Controls_all(7, :)), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'DisplayName', 'Vx');
    histogram(clipUpper(RMSE_Controls_all(8, :)), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'DisplayName', 'Vy');
    xlabel('Velocity RMSE (m/s)'); title('Velocity Error'); legend('show');

    %% Plot 4: Estimator Filter Performance Histogram
    figure('Name', 'Estimator RMSE Distribution', 'Color', bkgColor, 'WindowStyle', 'docked');
    hold on; grid on;
    histogram(clipUpper(RMSE_Filter_all(1, :)), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'DisplayName', 'Yaw Error');
    histogram(clipUpper(RMSE_Filter_all(2, :)), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'DisplayName', 'Pitch Error');
    histogram(clipUpper(RMSE_Filter_all(3, :)), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'DisplayName', 'Roll Error');
    xlabel('Filter RMSE (deg)'); ylabel('Frequency'); title('Estimator Component RMSE'); legend('show');

    %% Plot 5: Bias & Noise Sampling Distributions
    figure('Name', 'Bias & Noise Sampling Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    
    nexttile; histogram(gyroVals, num_bins, 'FaceColor', '#7E2F8E', 'FaceAlpha', 0.6);
    set(gca, 'XScale', 'log'); xlabel('Gyro Drift Power'); title('Gyro Noise Power'); grid on;

    nexttile; histogram(kg2Vals, num_bins, 'FaceColor', '#77AC30', 'FaceAlpha', 0.6);
    xlabel('Kg2 (g^2 bias)'); title('Kg2 Bias Sensitivity'); grid on;
    
    nexttile; histogram(grmaxVals, num_bins, 'FaceColor', '#4DBEEE', 'FaceAlpha', 0.6);
    xlabel('G\_RMAX'); title('G\_RMAX Distribution'); grid on;

    disp('Plot generation complete.');

    %% Plot 6: Approach Conditions (Evaluated at 2m Altitude)
    figure('Name', 'Approach Conditions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    appr_lat_error = nan(num_sims, 1);
    appr_vel = nan(num_sims, 3);
    alt_threshold = 1; % Evaluate approach conditions at 1 meters
    
    final_wp = Waypoints(end).Position;
    
    for i = 1:num_sims
        % Find the last valid state of this simulation run
        last_valid = find(~isnan(pos_all(i,:,1)), 1, 'last');
        
        if ~isempty(last_valid)
            alts = pos_all(i, 1:last_valid, 3);
            
            % Find the LAST index where altitude was ABOVE the threshold.
            % This ensures we catch the final descent, ignoring the takeoff phase.
            idx_above = find(alts > alt_threshold, 1, 'last');
            
            if ~isempty(idx_above)
                % The approach point is the index right after it crossed below 2m
                appr_idx = min(idx_above + 1, last_valid);
            else
                % Fallback if it never went above 2m
                appr_idx = last_valid; 
            end
            
            % 1. Approach Lateral Error (Euclidean distance to target in X-Y plane)
            appr_pos_xy = [pos_all(i, appr_idx, 1), pos_all(i, appr_idx, 2)];
            appr_lat_error(i) = norm(appr_pos_xy - final_wp(1:2));
            
            % 2. Approach Velocity
            appr_vel(i, :) = vel_all(i, appr_idx, :);
        end
    end
    
    % --- Tile 1: Approach Lateral Error (Successful Runs Only) ---
    nexttile; hold on; grid on;
    
    % Filter to only include successful trajectories
    appr_lat_error_success = appr_lat_error(is_success); 
    
    if ~isempty(appr_lat_error_success)
        histogram(appr_lat_error_success, num_bins, 'FaceColor', '#7E2F8E', 'FaceAlpha', 0.6);
    end
    xlabel('Lateral Distance from Target (m)', 'FontWeight', 'bold');
    ylabel('Frequency');
    title(sprintf('Approach Lateral Error (Alt = %.1f m)', alt_threshold));
    
    % --- Tile 2: Approach Velocity (Overlay Histogram) ---
    nexttile; hold on; grid on;
    histogram(appr_vel(:, 1), num_bins, 'FaceColor', '#0072BD', 'FaceAlpha', 0.6, 'DisplayName', 'V_x');
    histogram(appr_vel(:, 2), num_bins, 'FaceColor', '#D95319', 'FaceAlpha', 0.6, 'DisplayName', 'V_y');
    histogram(appr_vel(:, 3), num_bins, 'FaceColor', '#EDB120', 'FaceAlpha', 0.6, 'DisplayName', 'V_z');
    xlabel('Velocity (m/s)', 'FontWeight', 'bold');
    ylabel('Frequency');
    title(sprintf('Approach Velocity (Alt = %.1f m)', alt_threshold));
    legend('show');
end

%% Helper function for drawing orthographic trajectory projections
function drawMCTraj(ax, pos_all, pos_avg, wp_coords, Waypoints, x_idx, y_idx, x_lbl, y_lbl, alphaVal, is_success)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

    % Draw Waypoint regions
    for idx = 1:numel(Waypoints)
        wp = Waypoints(idx);
        if ~wp.IsPassAndGo && wp.Position(3) >= 0
            tol = wp.PosTol;
            px = wp.Position(x_idx);
            py = wp.Position(y_idx);
            vx = [px-tol, px+tol, px+tol, px-tol];
            vy = [py-tol, py-tol, py+tol, py+tol];
            patch(ax, vx, vy, 'g', 'FaceAlpha', 0.1, 'EdgeColor', '#4DBEEE');
        end
    end

    % Draw all runs with Success/Fail colors
    num_sims = size(pos_all, 1);
    for i = 1:num_sims
        if is_success(i)
            plotColor = [0, 1, 0, alphaVal]; % Green for success
        else
            plotColor = [1, 0, 0, alphaVal]; % Red for failure
        end
        plot(ax, pos_all(i, :, x_idx), pos_all(i, :, y_idx), 'Color', plotColor, 'LineWidth', 0.5);
    end

    % Draw Average Trajectory
    plot(ax, pos_avg(:, x_idx), pos_avg(:, y_idx), 'k', 'LineWidth', 1.5);

    % Draw Target Points
    plot(ax, wp_coords(x_idx, :), wp_coords(y_idx, :), 'rs', ...
         'MarkerFaceColor', 'r', 'MarkerSize', 6);
         
    xlabel(ax, x_lbl); ylabel(ax, y_lbl);
end