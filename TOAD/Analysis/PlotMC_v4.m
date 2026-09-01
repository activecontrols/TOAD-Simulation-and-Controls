function PlotMC_v4(filename)
% PlotMC_v4
% Focused Monte Carlo plotting for the V4 full-state target architecture.

close all;
bkgColor = 'w';
alphaVal = 0.12;

% Convention: Pitch (X), Yaw (Y), Roll (Z)
oranPitch = [0.85 0.60 0.45];
blueYaw   = [0.55 0.71 0.84];
yellRoll  = [0.93 0.81 0.54];
pyrColors = {oranPitch, blueYaw, yellRoll}; 

%% 1. Load Data
if nargin < 1 || isempty(filename)
    reqVars = {'out', 't_common', 'Lever_Radial', 'Lever_Axial', ...
               'J_Trans_Scale', 'J_Axial_Scale', 'J_Wobble_Coup', 'J_Trans_Coup', ...
               'Wind_Gain_all', 'Wind_Covar_all', 'RMSE_Wind_all', ...
               'RMSE_Controls_all', 'RMSE_Filter_all', 'MaxLESODist_all'};
    for i = 1:numel(reqVars)
        try evalin('base', [reqVars{i} ';']); catch, end
        eval([reqVars{i} ' = evalin(''base'', ''' reqVars{i} ''');']);
    end
else
    S = load(filename);
    vars = fieldnames(S);
    for i = 1:numel(vars), eval([vars{i} ' = S.' vars{i} ';']); end
end

num_sims = numel(out);

%% 2. Extract & Interpolate State Histories
firstValid = find(arrayfun(@(x) isempty(x.ErrorMessage), out), 1, 'first');
nState = size(extractLoggedState(out(firstValid).state_log), 2);
nTarget = size(extractLoggedState(out(firstValid).target_pos_log), 2);

actual_interp = nan(num_sims, numel(t_common), nState);
target_interp = nan(num_sims, numel(t_common), nTarget);

for i = 1:num_sims
    if ~isempty(out(i).ErrorMessage), continue; end 
    actual_interp(i,:,:) = interpolateLoggedState(out(i).state_log, t_common, nState);
    target_interp(i,:,:) = interpolateLoggedState(out(i).target_pos_log, t_common, nTarget);
end

quat_idx = 1:4; pos_idx = 5:7; vel_idx = 8:10; rate_idx = 11:13;

%% 3. Target-Relative Errors
pos_error_all = actual_interp(:,:,pos_idx) - target_interp(:,:,pos_idx);
vel_error_all = actual_interp(:,:,vel_idx) - target_interp(:,:,vel_idx);

pos_error_rmse = reshape(squeeze(sqrt(mean(pos_error_all.^2, 2, 'omitnan'))), [], 3);
vel_error_rmse = reshape(squeeze(sqrt(mean(vel_error_all.^2, 2, 'omitnan'))), [], 3);

Pos_RMSE_total = sqrt(sum(pos_error_rmse.^2, 2));
Vel_RMSE_total = sqrt(sum(vel_error_rmse.^2, 2));
target_pos_ref = squeeze(target_interp(firstValid,:,pos_idx));

%% 4. Kinematic Trajectory Overlay (3D + Projections)
figure('Name', 'MC 3D Trajectories', 'Color', bkgColor, 'WindowStyle', 'docked');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

% Main 3D Plot
axMain = nexttile(tl, 1, [3 3]); 
hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
xlabel(axMain, 'Pitch / X [m]'); ylabel(axMain, 'Yaw / Y [m]'); zlabel(axMain, 'Roll / Z [m]');
title(axMain, sprintf('MC Trajectories vs Full-State Target (%d Runs)', num_sims));

for i = 1:num_sims
    xyz = squeeze(actual_interp(i,:,pos_idx));
    if isempty(xyz), continue; end
    valid = all(isfinite(xyz),2);
    plot3(axMain, xyz(valid,1), xyz(valid,2), xyz(valid,3), 'Color', [0.25 0.25 0.25 alphaVal], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end
plot3(axMain, target_pos_ref(:,1), target_pos_ref(:,2), target_pos_ref(:,3), 'b--', 'LineWidth', 2.5, 'DisplayName', 'Target');

% Orthographic Projections
axTop = nexttile(tl, 4); 
drawMCTraj_v4(axTop, actual_interp(:,:,pos_idx), target_pos_ref, 1, 2, 'Pitch / X [m]', 'Yaw / Y [m]', alphaVal); 
title(axTop, 'Top View');

axSide = nexttile(tl, 8); 
drawMCTraj_v4(axSide, actual_interp(:,:,pos_idx), target_pos_ref, 1, 3, 'Pitch / X [m]', 'Roll / Z [m]', alphaVal); 
title(axSide, 'Side View');

axFront = nexttile(tl, 12); 
drawMCTraj_v4(axFront, actual_interp(:,:,pos_idx), target_pos_ref, 2, 3, 'Yaw / Y [m]', 'Roll / Z [m]', alphaVal); 
title(axFront, 'Front View');

%% 5. Kinematic Plot (3-Sigma State Distributions)
figure('Name', 'State Distributions (3-Sigma)', 'Color', bkgColor, 'WindowStyle', 'docked');
tiledlayout(2,3, 'TileSpacing', 'compact');
state_titles = {'Pitch (X)', 'Yaw (Y)', 'Roll (Z)'};
vars = {pos_idx, vel_idx};
ylbls = {'Pos', 'Vel'};
units = {'m', 'm/s'};

for r = 1:2
    for c = 1:3
        nexttile; hold on; grid on;
        data_block = squeeze(actual_interp(:,:,vars{r}(c)));
        target_block = squeeze(target_interp(firstValid,:,vars{r}(c))); % Extract specific target axis
        
        mu_val = mean(data_block, 1, 'omitnan');
        sig_val = std(data_block, 0, 1, 'omitnan');
        
        plot(t_common, data_block', 'Color', [0.45 0.68 0.88 alphaVal*2], 'HandleVisibility', 'off');
        plot(t_common, mu_val + 3*sig_val, 'r--', 'LineWidth', 1.5, 'DisplayName', '+3\sigma');
        plot(t_common, mu_val - 3*sig_val, 'r--', 'LineWidth', 1.5, 'DisplayName', '-3\sigma');
        plot(t_common, target_block, 'k', 'LineWidth', 2, 'DisplayName', 'Target');
        
        title(sprintf('%s - %s', ylbls{r}, state_titles{c}));
        xlabel('Time (s)'); ylabel(sprintf('%s (%s)', ylbls{r}, units{r}));
        if r == 1 && c == 1, legend('Location', 'best'); end
    end
end

%% 5b. Rotational Kinematics: Tilt Angle & Angular Rates (3-Sigma)
% Tilt angle (not Euler) is used for the attitude envelope since the
% flip maneuver passes through/near singular Euler configurations.
R33_actual = actual_interp(:,:,1).^2 - actual_interp(:,:,2).^2 - ...
             actual_interp(:,:,3).^2 + actual_interp(:,:,4).^2;
tilt_actual = acosd(max(min(R33_actual, 1), -1));

R33_target = target_interp(firstValid,:,1).^2 - target_interp(firstValid,:,2).^2 - ...
             target_interp(firstValid,:,3).^2 + target_interp(firstValid,:,4).^2;
tilt_target = acosd(max(min(R33_target, 1), -1));

figure('Name', 'Attitude & Angular Rates', 'Color', bkgColor, 'WindowStyle', 'docked');
tl_att = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

axTilt = nexttile(tl_att, 1, [1 3]); hold(axTilt, 'on'); grid(axTilt, 'on');
plot(axTilt, t_common, tilt_actual', 'Color', [0.45 0.68 0.88 alphaVal*2], 'HandleVisibility', 'off');
mu_tilt = mean(tilt_actual, 1, 'omitnan');
sig_tilt = std(tilt_actual, 0, 1, 'omitnan');
plot(axTilt, t_common, mu_tilt + 3*sig_tilt, 'r--', 'LineWidth', 1.5, 'DisplayName', '+3\sigma');
plot(axTilt, t_common, mu_tilt - 3*sig_tilt, 'r--', 'LineWidth', 1.5, 'DisplayName', '-3\sigma');
plot(axTilt, t_common, tilt_target, 'k', 'LineWidth', 2, 'DisplayName', 'Target');
title(axTilt, 'Tilt Angle from Vertical'); xlabel(axTilt, 'Time (s)'); ylabel(axTilt, 'Tilt (deg)');
legend(axTilt, 'Location', 'best');

rate_titles = {'Pitch Rate (X)', 'Yaw Rate (Y)', 'Roll Rate (Z)'};
for c = 1:3
    ax = nexttile(tl_att); hold(ax, 'on'); grid(ax, 'on');
    data_block   = squeeze(actual_interp(:,:,rate_idx(c))) * (180/pi);
    target_block = squeeze(target_interp(firstValid,:,rate_idx(c))) * (180/pi);
    mu_val  = mean(data_block, 1, 'omitnan');
    sig_val = std(data_block, 0, 1, 'omitnan');

    plot(ax, t_common, data_block', 'Color', [0.45 0.68 0.88 alphaVal*2], 'HandleVisibility', 'off');
    plot(ax, t_common, mu_val + 3*sig_val, 'r--', 'LineWidth', 1.5, 'DisplayName', '+3\sigma');
    plot(ax, t_common, mu_val - 3*sig_val, 'r--', 'LineWidth', 1.5, 'DisplayName', '-3\sigma');
    plot(ax, t_common, target_block, 'k', 'LineWidth', 2, 'DisplayName', 'Target');
    title(ax, rate_titles{c}); xlabel(ax, 'Time (s)'); ylabel(ax, 'Rate (deg/s)');
end

%% 6. Landing Angular-Rate Distribution
landing_rate = nan(num_sims,3);
for i = 1:num_sims
    rates = squeeze(actual_interp(i,:,rate_idx)) * (180/pi);
    pos = squeeze(actual_interp(i,:,pos_idx));
    if isempty(rates) || isempty(pos), continue; end
    
    validRows = all(isfinite(rates),2) & all(isfinite(pos),2);
    if any(validRows)
        above_thresh = pos(:,3) > 0.05;
        crossings = find(diff(above_thresh) == -1);
        if ~isempty(crossings)
            land_idx = crossings(end) + 1;
        else
            land_idx = find(validRows, 1, 'last');
        end
        landing_rate(i,:) = rates(land_idx,:);
    end
end

figure('Name', 'Landing Rates', 'Color', bkgColor, 'WindowStyle', 'docked'); 
tiledlayout(3,1, 'TileSpacing', 'compact');
rate_lbls = {'Pitch Rate', 'Yaw Rate', 'Roll Rate'};
for c = 1:3
    ax = nexttile;
    plotSmartHistogram(ax, landing_rate(:,c), pyrColors{c}, rate_lbls{c});
    title(sprintf('Landing %s', rate_lbls{c})); 
    xlabel('Rate (deg/s)'); grid on; legend('show');
end

%% 7. Control RMSE Distributions 
if size(RMSE_Controls_all,1) == 12
    labels = {'Attitude', 'Angular Rate', 'Position', 'Velocity'};
    axis_lbls = {'Pitch (X)','Yaw (Y)','Roll (Z)'};
    units = {'deg', 'deg/s', 'm', 'm/s'};
    multipliers = [(180/pi), (180/pi), 1, 1];
    
    for g = 1:4
        figure('Name', sprintf('%s Control RMSE', labels{g}), 'Color', bkgColor, 'WindowStyle', 'docked');
        tiledlayout(3,1, 'TileSpacing', 'compact');
        for c = 1:3
            idx = (g-1)*3 + c;
            ax = nexttile;
            data_raw = real(RMSE_Controls_all(idx,:)) * multipliers(g);
            plotSmartHistogram(ax, data_raw, pyrColors{c}, axis_lbls{c});
            title(sprintf('%s - %s', labels{g}, axis_lbls{c}));
            xlabel(sprintf('RMSE (%s)', units{g}));
            grid on; legend('show');
        end
    end
end

%% 8. Filter Attitude RMSE Plot (Degrees natively)
if exist('RMSE_Filter_all', 'var') && size(RMSE_Filter_all,1) >= 3
    figure('Name', 'Filter Attitude RMSE', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(3,1, 'TileSpacing', 'compact');
    filter_lbls = {'Pitch RMSE', 'Yaw RMSE', 'Roll RMSE'}; 
    
    for c = 1:3
        ax = nexttile;
        plotSmartHistogram(ax, real(RMSE_Filter_all(c,:)), pyrColors{c}, filter_lbls{c});
        title(sprintf('Filter %s', filter_lbls{c}));
        xlabel('RMSE (degrees)'); grid on; legend('show');
    end
end

%% 9. Curated Sensitivities
figure('Name', 'Focused Disturbance Sensitivities', 'Color', bkgColor, 'WindowStyle', 'docked');
tiledlayout(2,2, 'TileSpacing', 'compact');

nexttile; plotSensitivityScatter(Lever_Radial, Pos_RMSE_total, 'Radial Lever [m]', 'Total Pos RMSE [m]', 'Position Error vs Radial Lever');
nexttile; plotSensitivityScatter(Wind_Gain_all, RMSE_Wind_all, 'Wind Gain', 'Wind RMSE', 'Wind RMSE vs Wind Gain');
nexttile; plotSensitivityScatter(RMSE_Wind_all, Pos_RMSE_total, 'Wind RMSE', 'Total Pos RMSE [m]', 'Position Error vs Wind RMSE');
nexttile; plotSensitivityScatter(RMSE_Wind_all, Vel_RMSE_total, 'Wind RMSE', 'Total Vel RMSE [m/s]', 'Velocity Error vs Wind RMSE');

%% 10. Inertia Diagnostics
Attitude_RMSE_total = sqrt(sum(real(RMSE_Controls_all(1:3,:)).^2, 1))' * (180/pi);
inertiaColor = [0.65 0.55 0.75];

figure('Name', 'Inertia Diagnostics', 'Color', bkgColor, 'WindowStyle', 'docked');
tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

inertia_vals = {J_Trans_Scale, J_Axial_Scale, J_Wobble_Coup, J_Trans_Coup};
inertia_lbls = {'J Trans Scale', 'J Axial Scale', 'J Wobble Coup', 'J Trans Coup'};

for k = 1:4
    ax = nexttile;
    plotSmartHistogram(ax, inertia_vals{k}, inertiaColor, inertia_lbls{k});
    title(ax, inertia_lbls{k}); xlabel(ax, sprintf('%s (kg\\cdotm^2)', inertia_lbls{k}));
    grid(ax, 'on'); legend(ax, 'show');
end

for k = 1:4
    nexttile;
    plotSensitivityScatter(inertia_vals{k}, Attitude_RMSE_total, inertia_lbls{k}, ...
        'Attitude RMSE [deg]', sprintf('Attitude Error vs %s', inertia_lbls{k}));
end

%% 11. LESO Disturbance Sensitivities
if exist('MaxLESODist_all', 'var') && size(MaxLESODist_all, 1) == 6
    % Figure 1: Acceleration LESO vs Pos RMSE
    figure('Name', 'Acceleration LESO Sensitivities', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    accel_lbls = {'Accel X', 'Accel Y', 'Accel Z'};
    for c = 1:3
        nexttile;
        plotSensitivityScatter(MaxLESODist_all(c, :), Pos_RMSE_total, ...
            sprintf('Max %s', accel_lbls{c}), 'Total Pos RMSE [m]', ...
            sprintf('Pos RMSE vs %s', accel_lbls{c}));
    end

    % Figure 2: Attitude LESO vs Pos RMSE
    figure('Name', 'Attitude LESO Sensitivities', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact');
    att_lbls = {'Attitude Pitch', 'Attitude Yaw', 'Attitude Roll'};
    for c = 1:3
        nexttile;
        plotSensitivityScatter(MaxLESODist_all(c+3, :), Pos_RMSE_total, ...
            sprintf('Max %s', att_lbls{c}), 'Total Pos RMSE [m]', ...
            sprintf('Pos RMSE vs %s', att_lbls{c}));
    end
end

end

%% Local Helpers
function data = extractLoggedState(ts)
    if isempty(ts), data = []; return; end
    data = squeeze(ts.Data);
    if size(data,1) ~= numel(ts.Time), data = data'; end
end

function data_i = interpolateLoggedState(ts, t_common, nState)
    raw = extractLoggedState(ts);
    data_i = nan(numel(t_common), nState);
    for k = 1:nState
        data_i(:,k) = interp1(ts.Time, raw(:,k), t_common, 'linear', 'extrap');
    end
end

function plotSensitivityScatter(x, y, xlbl, ylbl, ttl)
    x = x(:); y = y(:); 
    good = isfinite(x) & isfinite(y);
    if ~any(good), return; end
    scatter(x(good), y(good), 32, 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.3);
    grid on; xlabel(xlbl); ylabel(ylbl); title(ttl);
end

function plotSmartHistogram(ax, data, clr, name)
    data = data(isfinite(data));
    if isempty(data), return; end
    hold(ax, 'on');

    % Bin the 2nd-98th percentile core at full resolution; collapse the
    % tails into single labeled overflow/underflow bins instead of
    % stretching the axis to fit a handful of outlier runs.
    loB = prctile(data, 2);
    hiB = prctile(data, 98);
    if hiB - loB < eps
        hiB = loB + max(abs(loB)*0.05, eps);
    end

    nBins = max(10, min(40, round(sqrt(numel(data)))));
    edges = linspace(loB, hiB, nBins+1);
    binW  = edges(2) - edges(1);
    centers = edges(1:end-1) + binW/2;

    counts    = histcounts(data, edges);
    underflow = sum(data < loB);
    overflow  = sum(data > hiB);

    bar(ax, centers, counts, 1, 'FaceColor', clr, 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'DisplayName', name);

    if underflow > 0
        bar(ax, loB - binW, underflow, binW, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, loB - binW, underflow, num2str(underflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    if overflow > 0
        bar(ax, hiB + binW, overflow, binW, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, hiB + binW, overflow, num2str(overflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end

    % Log-scale the y-axis when one bin dwarfs the rest, so smaller bars
    % stay visible instead of flattening to nothing next to the spike.
    allCounts = [counts, underflow, overflow];
    peak = max(allCounts);
    rest = allCounts(allCounts < peak & allCounts > 0);
    if ~isempty(rest) && peak > 4 * max(rest)
        set(ax, 'YScale', 'log');
        ylim(ax, [0.8, peak*1.3]);
    end
end

function drawMCTraj_v4(ax, pos_all, target_pos, x_idx, y_idx, x_lbl, y_lbl, alphaVal)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    
    % Draw run trajectories
    for i = 1:size(pos_all, 1)
        xy = squeeze(pos_all(i, :, [x_idx, y_idx]));
        if isempty(xy), continue; end
        valid = all(isfinite(xy), 2);
        plot(ax, xy(valid, 1), xy(valid, 2), 'Color', [0.25 0.25 0.25 alphaVal], 'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
    
    % Draw target reference
    plot(ax, target_pos(:, x_idx), target_pos(:, y_idx), 'b--', 'LineWidth', 2);
    xlabel(ax, x_lbl); ylabel(ax, y_lbl);
end