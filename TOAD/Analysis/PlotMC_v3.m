function PlotMC_v3(filename)
    % PlotMC_v3 Generates combined trajectory and sensitivity overlays
    close all;
    bkgColor = 'w';
    alphaVal = 0.2; % Transparency for 3D lines

   %% 1. Handle Input Arguments & Data Loading
    if nargin < 1 || isempty(filename)
        disp('Pulling data from the base workspace...');
        reqVars = {'pos_all', 'vel_all', 'pos_avg', 't_common',...
                   'Lever_Radial', 'Lever_Axial', 'J_Trans_Scale', 'J_Axial_Scale', ...
                   'RMSE_Controls_all', 'RMSE_Filter_all', ...
                   'GyroNoisePower_vals', 'Kg2_vals', 'G_RMAX_vals', ...
                   'kGrom_vals', 'bGrom_vals', 'G', 'Gust_all'};
        
        for i = 1:length(reqVars)
            try
                % Corrected to pull directly into function workspace
                eval([reqVars{i} ' = evalin(''base'', ''' reqVars{i} ''');']);
            catch
                error('Variable "%s" not found in base workspace.', reqVars{i});
            end
        end
    else
        load(filename);
    end

    Waypoints = TrajectoryBuilder();
    numWP = numel(Waypoints);
    wp_coords = zeros(3, numWP);
    for idx = 1:numWP, wp_coords(:, idx) = Waypoints(idx).Position; end
    
    %% 2. Calculate Success Criteria (V2 Logic)
    final_wp = Waypoints(end).Position;
    final_tol = Waypoints(end).PosTol;
    num_sims = size(pos_all, 1);
    
    is_success = false(num_sims, 1);
    meets_final_pos = false(num_sims, 1); 
    
    % New: Pre-allocate for approach criteria storage
    gs_criteria_all = NaN(num_sims, 1);
    lv_criteria_all = NaN(num_sims, 1);

    for i = 1:num_sims
        last_valid = find(~isnan(pos_all(i,:,1)), 1, 'last');
        if ~isempty(last_valid)
            % Final Position Check
            final_pos = squeeze(pos_all(i, last_valid, :));
            pos_met = all(abs(final_pos(:) - final_wp(:)) <= final_tol(:));
            meets_final_pos(i) = pos_met; 
            
            % Approach Criteria Check (2m crossing)
            appr_met = false;
            alts = pos_all(i, 1:last_valid, 3);
            idx_above = find(alts > 1.0, 1, 'last'); 
            
            if ~isempty(idx_above)
                appr_idx = min(idx_above + 1, last_valid);
                appr_alt = max(alts(appr_idx), 0.01); 
                
                appr_pos_xy = [pos_all(i, appr_idx, 1), pos_all(i, appr_idx, 2)];
                lat_err = norm(appr_pos_xy - final_wp(1:2));
                glideslope_deg = atan2d(lat_err, appr_alt); 
                
                appr_vel_xy = [vel_all(i, appr_idx, 1), vel_all(i, appr_idx, 2)];
                lat_vel = norm(appr_vel_xy);
                
                % New: Save calculated values for *all* sims to plot later
                gs_criteria_all(i) = glideslope_deg;
                lv_criteria_all(i) = lat_vel;

                if glideslope_deg <= 15.0 && lat_vel <= 0.3
                    appr_met = true;
                end
            end
            
            if pos_met && appr_met, is_success(i) = true; end
        end
    end
    
    disp(['Successful Runs: ', num2str(sum(is_success)), ' / ', num2str(num_sims)]);
    
    %% 3. Setup Helper Variables
    gyroVals = cell2mat(GyroNoisePower_vals);
    if iscell(Kg2_vals), kg2Vals = cell2mat(Kg2_vals); else, kg2Vals = Kg2_vals; end
    if iscell(G_RMAX_vals), grmaxVals = cell2mat(G_RMAX_vals); else, grmaxVals = G_RMAX_vals; end

    % Indices for Scatter plots
    idx_pass = is_success;
    idx_fail = ~is_success;

    %% Plot 1: 3D Mission Profile Overlay
    figure('Name', 'MC Mission Trajectories', 'Color', bkgColor, 'WindowStyle', 'docked');
    tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    axMain = nexttile(tl, 1, [3 3]); 
    hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
    xlabel(axMain, 'North [m]'); ylabel(axMain, 'West [m]'); zlabel(axMain, 'Altitude [m]');
    
    pct_final_pos = (sum(meets_final_pos) / num_sims) * 100;
    pct_all_criteria = (sum(is_success) / num_sims) * 100;
    title(axMain, {sprintf('3D Trajectory Envelope (%d Runs)', num_sims), ...
                   sprintf('Meets Final Position: %.1f%% | Meets All Criteria: %.1f%%', pct_final_pos, pct_all_criteria)});

    % New: Draw properly defined 3D Waypoint Tolerance Boxes (not just 2D)
    for idx = 1:numel(Waypoints)
        wp = Waypoints(idx);
        if ~wp.IsPassAndGo && wp.Position(3) >= 0
            % Draw cuboid box of dimensions 2*tol
            c = wp.Position; % Centroid
            t = wp.PosTol; % Tolerance vector
            % Define Vertices for a cuboid centered at WP
            V = [c(1)-t(1) c(2)-t(2) c(3)-t(3); c(1)+t(1) c(2)-t(2) c(3)-t(3);
                 c(1)+t(1) c(2)+t(2) c(3)-t(3); c(1)-t(1) c(2)+t(2) c(3)-t(3);
                 c(1)-t(1) c(2)-t(2) c(3)+t(3); c(1)+t(1) c(2)-t(2) c(3)+t(3);
                 c(1)+t(1) c(2)+t(2) c(3)+t(3); c(1)-t(1) c(2)+t(2) c(3)+t(3)];
            % Define Faces
            F = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
            patch(axMain, 'Vertices', V, 'Faces', F, 'FaceColor', 'g', 'FaceAlpha', 0.1, 'EdgeColor', '#4DBEEE');
        end
    end

    % Draw Runs Transparently
    for i = 1:num_sims
        if meets_final_pos(i) 
            if is_success(i), plotColor = [0, 1, 0, alphaVal]; else, plotColor = [1, 0, 0, alphaVal]; end
            plot3(axMain, pos_all(i,:,1), pos_all(i,:,2), pos_all(i,:,3), 'Color', plotColor, 'LineWidth', 0.5);
        end
    end
    plot3(axMain, pos_avg(:,1), pos_avg(:,2), pos_avg(:,3), 'k', 'LineWidth', 2.5);
    
    axTop = nexttile(tl, 4); drawMCTraj(axTop, pos_all, pos_avg, wp_coords, Waypoints, 2, 1, 'West', 'North', alphaVal, is_success, meets_final_pos); title(axTop, 'Top View');
    axSide = nexttile(tl, 8); drawMCTraj(axSide, pos_all, pos_avg, wp_coords, Waypoints, 1, 3, 'North', 'Alt', alphaVal, is_success, meets_final_pos); title(axSide, 'Side View');
    axFront = nexttile(tl, 12); drawMCTraj(axFront, pos_all, pos_avg, wp_coords, Waypoints, 2, 3, 'West', 'Alt', alphaVal, is_success, meets_final_pos); title(axFront, 'Front View');

    %% New Plot 2: Monte Carlo Parameter Distributions (Added per request)
    figure('Name', 'Monte Carlo Parameter Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact'); 

    % Subplot 1: Lever Arm Disturbances
    nexttile; hold on; grid on;
    histogram(Lever_Radial * 1000, 30, 'FaceAlpha', 0.3, 'DisplayName', 'Radial Shift (mm)');
    histogram(Lever_Axial * 1000, 30, 'FaceAlpha', 0.3, 'DisplayName', 'Axial Shift (mm)');
    xlabel('Shift (mm)'); ylabel('Frequency'); title('Lever Arm Disturbances');
    legend('show', 'Location', 'best');

    % Subplot 2: Inertia Scale Disturbances
    nexttile; hold on; grid on;
    histogram(J_Trans_Scale, 30, 'FaceAlpha', 0.3, 'DisplayName', 'Transverse Scale (I_{xx}, I_{yy})');
    histogram(J_Axial_Scale, 30, 'FaceAlpha', 0.3, 'DisplayName', 'Axial Scale (I_{zz})');
    xlabel('Scale Error'); ylabel('Frequency'); title('Inertia Scale Disturbances');
    legend('show', 'Location', 'best');

    %% Plot 3: Targeted Sensitivities (Colorized by Success)
    figure('Name', 'Sensitivities by Pass/Fail', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout('flow', 'TileSpacing', 'compact'); 
    
    pos_x = RMSE_Controls_all(4, :); pos_y = RMSE_Controls_all(5, :);
    att_pitch = RMSE_Controls_all(2, :); att_roll = RMSE_Controls_all(3, :);
    
    lat_pos_err = sqrt(pos_x.^2 + pos_y.^2);
    att_trans_err = sqrt(att_pitch.^2 + att_roll.^2); 

    nexttile; plotSuccessScatter(Lever_Radial, lat_pos_err, 'Lever Radial (m)', 'Lat Pos RMSE (m)', 'Pos vs Off-Center CG', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(J_Axial_Scale, att_roll, 'Axial Scale Err (I_{zz})', 'Roll RMSE (rad)', 'Roll vs Roll Inertia', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(J_Trans_Scale, att_trans_err, 'Transverse Scale (I_{xx}, I_{yy})', 'Pitch/Yaw RMSE (rad)', 'Pitch/Yaw vs Transverse Scale', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(Lever_Axial, att_trans_err, 'Lever Axial Shift (m)', 'Pitch/Yaw RMSE (rad)', 'Pitch/Yaw vs Forward/Aft CG', idx_pass, idx_fail);
    legend('show', 'Location', 'best');

    %% Plot 4: Hardware Parameters & Error vs Drift (Colorized)
    figure('Name', 'Hardware Params vs Stability', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(2, 3, 'TileSpacing', 'compact');
    
    nexttile; plotSuccessScatter(gyroVals, att_roll, 'Gyro Drift Power', 'Roll RMSE (rad)', 'Roll vs Gyro Drift', idx_pass, idx_fail); set(gca, 'XScale', 'log');
    nexttile; plotSuccessScatter(gyroVals, lat_pos_err, 'Gyro Drift Power', 'Lat Pos RMSE (m)', 'Lat Pos vs Gyro Drift', idx_pass, idx_fail); set(gca, 'XScale', 'log');
    nexttile; plotSuccessScatter(kg2Vals, lat_pos_err, 'Kg2 Bias', 'Lat Pos RMSE (m)', 'Lat Pos vs Kg2 Bias', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(grmaxVals, lat_pos_err, 'G\_RMAX', 'Lat Pos RMSE (m)', 'Lat Pos vs G\_RMAX', idx_pass, idx_fail);
    
    filt_roll  = RMSE_Filter_all(3, :);
    nexttile([1 2]); plotSuccessScatter(filt_roll, lat_pos_err, 'Filter Roll RMSE (deg)', 'Lat Pos RMSE (m)', 'Control Lat Pos vs Estimator Roll Error', idx_pass, idx_fail);
    set(gca, 'YScale', 'log');

    %% Plot 5: Manifold Surfaces with Pass/Fail Overlay (Transparency & CLim updated)
    figure('Name', 'Stability Manifolds', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    
    [KG, GR] = meshgrid(linspace(min(kg2Vals), max(kg2Vals), 100), linspace(min(grmaxVals), max(grmaxVals), 100));
    
    % Surface 1: Lateral Pos Error
    nexttile; hold on;
    F_lat = scatteredInterpolant(kg2Vals(:), grmaxVals(:), lat_pos_err(:), 'linear', 'none');
    surf(KG, GR, F_lat(KG, GR), 'EdgeColor', 'none', 'FaceAlpha', 0.8); % Changed: Alpha to 0.8
    
    scatter3(kg2Vals(idx_fail), grmaxVals(idx_fail), lat_pos_err(idx_fail), 15, 'r', 'filled', 'MarkerEdgeColor', 'k');
    scatter3(kg2Vals(idx_pass), grmaxVals(idx_pass), lat_pos_err(idx_pass), 15, 'g', 'filled', 'MarkerEdgeColor', 'k');
    
    view(2); colormap turbo; colorbar; set(gca, 'ZScale', 'log');
    xlabel('Kg2 Bias Sensitivity'); ylabel('G\_RMAX'); title('Lat Pos Error Manifold');
    
    % Changed: CLim limits now dynamic based on success criteria (anchored to worst passing error)
    min_lat = min(lat_pos_err); max_pass_lat = max(lat_pos_err(idx_pass));
    if ~isempty(max_pass_lat) && ~isempty(min_lat), clim([min_lat, max_pass_lat]); end
    
    % Surface 2: Roll Error
    nexttile; hold on;
    F_roll = scatteredInterpolant(kg2Vals(:), grmaxVals(:), filt_roll(:), 'linear', 'none');
    surf(KG, GR, F_roll(KG, GR), 'EdgeColor', 'none', 'FaceAlpha', 0.8); % Changed: Alpha to 0.8
    
    scatter3(kg2Vals(idx_fail), grmaxVals(idx_fail), filt_roll(idx_fail), 15, 'r', 'filled', 'MarkerEdgeColor', 'k');
    scatter3(kg2Vals(idx_pass), grmaxVals(idx_pass), filt_roll(idx_pass), 15, 'g', 'filled', 'MarkerEdgeColor', 'k');
    
    view(2); colormap turbo; colorbar; set(gca, 'ZScale', 'log');
    xlabel('Kg2 Bias Sensitivity'); ylabel('G\_RMAX'); title('MEKF Roll Error Manifold');

    % Changed: CLim limits now dynamic based on success criteria (anchored to worst passing error)
    min_roll = min(filt_roll); max_pass_roll = max(filt_roll(idx_pass));
    if ~isempty(max_pass_roll) && ~isempty(min_roll), clim([min_roll, max_pass_roll]); end

    %% Plot 6: Grommet Transmissibility
    figure('Name', 'Grommet Transmissibility', 'Color', bkgColor, 'WindowStyle', 'docked');
    hold on; grid on;
    f = logspace(0, log10(3000), 1000); 
    mBoard = 0.1; 
    k_vals = cell2mat(kGrom_vals); b_vals = cell2mat(bGrom_vals);
    K_nom = mean(k_vals); B_nom = mean(b_vals);

    T_all = zeros(length(k_vals), length(f));
    for i = 1:length(k_vals)
        f_n_iter = 1 / (2 * pi) * sqrt(k_vals(i) / mBoard);
        r = f ./ f_n_iter;
        T_all(i, :) = sqrt((1 + (2 * b_vals(i) * r).^2) ./ ((1 - r.^2).^2 + (2 * b_vals(i) * r).^2));
    end

    T_lower = prctile(T_all, 1, 1); T_upper = prctile(T_all, 99, 1);
    f_n_nom = 1 / (2 * pi) * sqrt(K_nom / mBoard);
    r_nom = f ./ f_n_nom;
    T_nom = sqrt((1 + (2 * B_nom * r_nom).^2) ./ ((1 - r_nom.^2).^2 + (2 * B_nom * r_nom).^2));

    fill([f, fliplr(f)], [T_upper, fliplr(T_lower)], [0 0.4470 0.7410], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(f, T_nom, 'k-', 'LineWidth', 1.5);
    
    T_fn = sqrt((1 + (2 * B_nom * 1)^2) / ((1 - 1^2)^2 + (2 * B_nom * 1)^2));
    plot(f_n_nom, T_fn, 'ro', 'MarkerSize', 7, 'MarkerFaceColor', 'r');
    text(f_n_nom * 1.15, T_fn, sprintf('f_n = %.1f Hz', f_n_nom), 'Color', 'r', 'FontWeight', 'bold', 'VerticalAlignment', 'bottom');

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([min(f) max(f)]); ylim([1e-2, max(T_upper)*2]);
    xlabel('Frequency (Hz)'); ylabel('Transmissibility (T)');
    title(sprintf('Transmissibility Uncertainty: %s\n(Durometer: %s | K: %g | C: %g)', string(G.Name), string(G.Durometer), str2double(string(G.K)), str2double(string(G.C))), 'Interpreter', 'none'); 

   %% New Plot 7: Approach Criteria Distributions (Corrected Binning)
    figure('Name', 'Approach Criteria Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact'); 

    nexttile; hold on; grid on;
    % Plot All Sims and grab the object to extract its BinEdges
    h_gs = histogram(gs_criteria_all, 40, 'DisplayName', 'All Sims');
    % Force the Failed histogram to use the EXACT same bins
    histogram(gs_criteria_all(idx_fail), h_gs.BinEdges, 'DisplayName', 'Failed');
    
    xlabel('Glideslope Angle (deg)'); ylabel('Frequency'); title('Approach Glideslope Angle Distribution');
    xline(10, 'r--', 'LineWidth', 2, 'DisplayName', 'Criterion');
    xlim([0, min(prctile(gs_criteria_all, 99.5), 180)]);
    legend('show', 'Location', 'best');

    nexttile; hold on; grid on;
    % Plot All Sims and grab the object to extract its BinEdges
    h_lv = histogram(lv_criteria_all, 40, 'DisplayName', 'All Sims');
    % Force the Failed histogram to use the EXACT same bins
    histogram(lv_criteria_all(idx_fail), h_lv.BinEdges, 'DisplayName', 'Failed');
    
    xlabel('Lateral Velocity (m/s)'); ylabel('Frequency'); title('Approach Lateral Velocity Distribution');
    xline(0.3, 'r--', 'LineWidth', 2, 'DisplayName', 'Criterion');
    xlim([0, min(prctile(lv_criteria_all, 99.5), 2)]);
    legend('show', 'Location', 'best');

    %% Plot 8: Targeted Sensitivities (Colorized by Success)
    figure('Name', 'Pos Error w/ Wind', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout('horizontal', 'TileSpacing', 'compact');  

    nexttile; plotSuccessScatter(Gust_all(1,:), lat_pos_err, 'Max Gust X (m/s)', 'Lat Pos RMSE (m)', 'Pos vs X Wind', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(Gust_all(2,:), lat_pos_err, 'Max Gust Y (m/s)', 'Lat Pos RMSE (m)', 'Pos vs Y Wind', idx_pass, idx_fail);
    nexttile; plotSuccessScatter(Gust_all(3,:), lat_pos_err, 'Max Gust Z (m/s)', 'Lat Pos RMSE (m)', 'Pos vs Z Wind', idx_pass, idx_fail);
    legend('show', 'Location', 'best');
end

%% Local Helper function for Success/Fail Scatters
function plotSuccessScatter(x, y, xlbl, ylbl, t, idx_pass, idx_fail)
    hold on; grid on;
    scatter(x(idx_fail), y(idx_fail), 30, 'r', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.3, 'DisplayName', 'Failed');
    scatter(x(idx_pass), y(idx_pass), 30, 'g', 'filled', 'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8, 'DisplayName', 'Success');
    xlabel(xlbl); ylabel(ylbl); title(t);
end

%% Helper function for Orthographic Trajectory Projections
function drawMCTraj(ax, pos_all, pos_avg, wp_coords, Waypoints, x_idx, y_idx, x_lbl, y_lbl, alphaVal, is_success, meets_final_pos)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');

    % The orthographic boxes remain 2D planes ('vx', 'vy')
    for idx = 1:numel(Waypoints)
        wp = Waypoints(idx);
        if ~wp.IsPassAndGo && wp.Position(3) >= 0
            tol = wp.PosTol; px = wp.Position(x_idx); py = wp.Position(y_idx);
            vx = [px-tol(x_idx), px+tol(x_idx), px+tol(x_idx), px-tol(x_idx)];
            vy = [py-tol(y_idx), py-tol(y_idx), py+tol(y_idx), py+tol(y_idx)];
            patch(ax, vx, vy, 'g', 'FaceAlpha', 0.1, 'EdgeColor', '#4DBEEE');
        end
    end

    for i = 1:size(pos_all, 1)
        if meets_final_pos(i) 
            if is_success(i), c = [0, 1, 0, alphaVal]; else, c = [1, 0, 0, alphaVal]; end
            plot(ax, pos_all(i, :, x_idx), pos_all(i, :, y_idx), 'Color', c, 'LineWidth', 0.5);
        end
    end

    plot(ax, pos_avg(:, x_idx), pos_avg(:, y_idx), 'k', 'LineWidth', 1.5);
    plot(ax, wp_coords(x_idx, :), wp_coords(y_idx, :), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    xlabel(ax, x_lbl); ylabel(ax, y_lbl);
end