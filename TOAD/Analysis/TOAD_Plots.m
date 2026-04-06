% Load the Simulink model & Waypoints
LoadTOADSim;
Waypoints = TrajectoryBuilder(); % Load the new TOADWaypoint array
load_system("TOAD_Simulation.slx");

% Run the simulation
simOut = sim("TOAD_Simulation.slx");

% Extract state timeseries
stateTimeSeries = simOut.state_log;
MEKFsTimeSeries = simOut.MEKF_state;
MEKFpTimeSeries = simOut.MEKF_P;
inputTimeSeries = simOut.inputCMD;
QerrTimeSeries  = simOut.AttErrorControl;

if ~isprop(simOut, 'target_pos_log')
    error('target_pos_log not found in simOut. Please ensure trg is logged to the workspace.');
end

%% Post-Processing & Plotting Script
formatData = @(ts) correctDims(ts.Data, length(ts.Time));

% Extract Time & Data
t = stateTimeSeries.Time;
x_true    = formatData(stateTimeSeries); 
x_est     = formatData(MEKFsTimeSeries);  
P_diag    = formatData(MEKFpTimeSeries);  
u_cmd     = formatData(inputTimeSeries);  
q_err_raw = formatData(QerrTimeSeries);   
pos_targ  = formatData(simOut.target_pos_log); 

% --- Variable Mapping ---
q_true   = x_true(:, 1:4);
pos_true = x_true(:, 5:7);
vel_true = x_true(:, 8:10); % Extracted Velocity for Kinematics plot

q_est    = x_est(:, 1:4);
pos_est  = x_est(:, 5:7);

sig3_att  = 3 * sqrt(P_diag(:, 1:3));
sig3_pos  = 3 * sqrt(P_diag(:, 4:6));

gimbal_theta = u_cmd(:, 1);
gimbal_phi   = u_cmd(:, 2);
thrust_cmd   = u_cmd(:, 3);

% --- Calculations ---
% 1. Reconstruct Target Attitude & Errors
q_target = zeros(size(q_true));
att_err_angle = zeros(length(t), 3);
for k = 1:length(t)
    q_T = q_true(k, :)';
    H = HamiltonianProd(q_T); 
    q_target(k, :) = (H * q_err_raw(k, :)')';
    
    q_E_conj = [q_est(k,1); -q_est(k,2:4)']; 
    H_inv = HamiltonianProd(q_E_conj);
    q_diff = H_inv * q_T;
    
    if q_diff(1) < 0, q_diff = -q_diff; end
    q_diff = q_diff / norm(q_diff);
    att_err_angle(k, :) = 2 * q_diff(2:4)';
end

eul_true = quat2eul(q_true, 'ZYX');
eul_targ = quat2eul(q_target, 'ZYX');

% 2. Infer Active Waypoint Index & Extract Waypoint Data
numWP = numel(Waypoints);
wp_coords = zeros(3, numWP);
wp_max_vels = zeros(3, numWP); 

for idx = 1:numWP
    wp_coords(:, idx)   = Waypoints(idx).Position;
    wp_max_vels(:, idx) = Waypoints(idx).MaxVel; 
end

wp_index = ones(length(t), 1);
current_wp = 1;
for k = 1:length(t)
    if current_wp < numWP
        if norm(pos_targ(k, :)' - wp_coords(:, current_wp+1)) < 1e-3
            current_wp = current_wp + 1;
        end
    end
    wp_index(k) = current_wp;
end

% Find transition times for mission phase plotting
trans_idx = find(diff(wp_index) > 0);
trans_times = t(trans_idx);

% --- APPLY TIME SHIFT ---
% Shift all times relative to the end of the first waypoint (Liftoff)
if ~isempty(trans_times)
    t_liftoff = trans_times(1);
    t = t - t_liftoff;
    trans_times = trans_times - t_liftoff;
end


% =========================================================================
% --- Plot 1: 3D Mission Profile ---
% =========================================================================
figure('Name', 'Figure 5.1: Mission Trajectory', 'WindowStyle', 'docked');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

N = pos_true(:,1); W = pos_true(:,2); U = pos_true(:,3);

drawTraj = @(ax, x_idx, y_idx, x_lbl, y_lbl) ...
    plotProjection(ax, N, W, U, wp_coords, Waypoints, x_idx, y_idx, x_lbl, y_lbl);

axMain = nexttile(tl, 1, [3 3]); 
hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
xlabel(axMain, 'North [m]'); ylabel(axMain, 'West [m]'); zlabel(axMain, 'Altitude [m]');
title(axMain, '3D Mission Trajectory');

patch(axMain, [N; NaN], [W; NaN], [U; NaN], [t; NaN], ...
      'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5);
cb = colorbar(axMain); cb.Label.String = 'Time [s]'; colormap(axMain, 'turbo');

plot3(axMain, wp_coords(1,:), wp_coords(2,:), wp_coords(3,:), 'r--s', ...
    'LineWidth', 1.5, 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');

for idx = 1:numWP
    wp = Waypoints(idx);
    if ~wp.IsPassAndGo && wp.Position(3) >= 0
        tol = wp.PosTol;
        v_cube = 2 * tol * [0.5 -0.5 -0.5; 0.5 0.5 -0.5; -0.5 0.5 -0.5; -0.5 -0.5 -0.5; ...
                            0.5 -0.5 0.5; 0.5 0.5 0.5; -0.5 0.5 0.5; -0.5 -0.5 0.5] + wp.Position'; 
        [f_cube, ~] = convhull(v_cube);
        patch(axMain, 'Vertices', v_cube, 'Faces', f_cube, 'FaceColor', 'g', ...
              'FaceAlpha', 0.1, 'EdgeColor', '#4DBEEE');
        text(axMain, wp.Position(1), wp.Position(2), wp.Position(3) + tol + 1, ...
             strtrim(wp.Label), 'FontSize', 9, 'FontWeight', 'bold');
    end
end

axTop = nexttile(tl, 4);  drawTraj(axTop, 2, 1, 'West [m]', 'North [m]'); title(axTop, 'Top View (X-Y)');
axSide = nexttile(tl, 8); drawTraj(axSide, 1, 3, 'North [m]', 'Altitude [m]'); title(axSide, 'Side View (X-Z)');
axFront = nexttile(tl, 12); drawTraj(axFront, 2, 3, 'West [m]', 'Altitude [m]'); title(axFront, 'Front View (Y-Z)');


% =========================================================================
% --- Plot 2: Estimator Performance ---
% =========================================================================
figure('Name', 'Figure 5.2: Estimator Performance', 'WindowStyle', 'docked');
subplot(2,1,1); hold on; grid on;
plot(t, pos_est - pos_true, 'LineWidth', 1.5);
plot(t, sig3_pos, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(t, -sig3_pos, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
ylabel('Pos Error [m]'); legend('N', 'W', 'U'); title('Position Error bounds \pm 3\sigma');

subplot(2,1,2); hold on; grid on;
plot(t, rad2deg(att_err_angle), 'LineWidth', 1.5);
plot(t, rad2deg(sig3_att), 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(t, -rad2deg(sig3_att), 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');
ylabel('Att Error [deg]'); legend('\alpha_x', '\alpha_y', '\alpha_z'); title('Attitude Error bounds \pm 3\sigma');


% =========================================================================
% --- Plot 3: Controller Performance ---
% =========================================================================
figure('Name', 'Figure 5.3: Controller Performance', 'WindowStyle', 'docked');

subplot(3,1,1); hold on; grid on;
plot(t, rad2deg(eul_targ(:,3)), '--', 'Color', '#D95319', 'LineWidth', 1.5);
plot(t, rad2deg(eul_true(:,3)), '-', 'Color', '#D95319', 'LineWidth', 1);
plot(t, rad2deg(eul_targ(:,2)), '--', 'Color', '#0072BD', 'LineWidth', 1.5);
plot(t, rad2deg(eul_true(:,2)), '-', 'Color', '#0072BD', 'LineWidth', 1);
ylabel('Angle [deg]'); legend('Roll Cmd', 'Roll True', 'Pitch Cmd', 'Pitch True', 'Location', 'best');
title('Attitude Tracking');

subplot(3,1,2); hold on; grid on;
plot(t, rad2deg(gimbal_theta), 'Color', '#7E2F8E', 'LineWidth', 1.5);
plot(t, rad2deg(gimbal_phi), 'Color', '#77AC30', 'LineWidth', 1.5);
ylabel('Gimbal [deg]'); legend('\theta', '\phi');
title('Gimbal Deflection');

subplot(3,1,3); hold on; grid on;
area(t, thrust_cmd, 'FaceColor', '#EDB120', 'FaceAlpha', 0.4, 'EdgeColor', '#A27712', 'LineWidth', 1.5);
ylabel('Thrust [N]'); xlabel('Time [s]');
title('Thrust Command');


% =========================================================================
% --- Plot 4: Mission Kinematics (Per-Axis Triple Plot) ---
% =========================================================================
figure('Name', 'Figure 5.4: Mission Kinematics', 'WindowStyle', 'docked');
tl_kin = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_kin, 'Mission Kinematics: Per-Axis Position and Velocity Envelopes', 'FontWeight', 'bold');
xlabel(tl_kin, 'Time [s]', 'FontWeight', 'bold'); 

labels_pos = {'North [m]', 'West [m]', 'Altitude [m]'};
labels_vel = {'North Vel [m/s]', 'West Vel [m/s]', 'Vertical Vel [m/s]'};
colors = {'#0072BD', '#D95319', '#EDB120'};

for ax_idx = 1:3
    % --- Position Tracking Subplot ---
    ax_pos = nexttile(tl_kin, (ax_idx-1)*2 + 1); hold on; grid on;
    plot(t, pos_targ(:, ax_idx), 'k--', 'LineWidth', 1.5);
    plot(t, pos_true(:, ax_idx), 'Color', colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_pos{ax_idx}, 'FontWeight', 'bold');
    
    % Massive Top Padding (40%) to create a clear header for text
    y_limits = ylim(ax_pos);
    y_range = diff(y_limits);
    if y_range == 0, y_range = 1; end 
    ylim(ax_pos, [y_limits(1) - max(y_range*0.1, 0.5), y_limits(2) + max(y_range*0.40, 1.0)]);
    
    if ax_idx == 1
        legend('Target Position', 'True Position', 'Location', 'best'); 
    end
    
    % --- Velocity Tracking Subplot ---
    ax_vel = nexttile(tl_kin, (ax_idx-1)*2 + 2); hold on; grid on;
    max_vel_curve = wp_max_vels(ax_idx, wp_index)';
    
    plot(t, max_vel_curve, 'r--', 'LineWidth', 1.5);
    plot(t, -max_vel_curve, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
    plot(t, vel_true(:, ax_idx), 'Color', colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_vel{ax_idx}, 'FontWeight', 'bold');
    
    % Massive Top Padding (40%) to create a clear header for text
    max_v = max(max_vel_curve);
    if max_v == 0, max_v = 1; end
    ylim(ax_vel, [-(max_v * 1.15), max_v * 1.45]); 
    
    if ax_idx == 1
        legend('Max Velocity Limit', 'True Velocity', 'Location', 'best'); 
    end
    
    % --- Add Custom Lines & Floating Boxes ---
    for ax = [ax_pos, ax_vel]
        axes(ax);
        
        y_bot = min(ylim(ax));
        y_range = diff(ylim(ax));
        
        % The magic fix: Draw lines that STOP 10% before the top edge
        y_line_top = min(ylim(ax)) + y_range * 0.90; 
        
        for i = 1:length(trans_times)
            plot([trans_times(i), trans_times(i)], [y_bot, y_line_top], 'k:', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        
        % Only draw labels on the top row
        if ax_idx == 1
            phase_starts = [t(1); trans_times]; 
            phase_ends = [trans_times; t(end)];
            for i = 1:numel(phase_starts)
                mid_time = (phase_starts(i) + phase_ends(i)) / 2;
                
                % Place text exactly 8% down from the absolute top (Safely above the lines!)
                y_text = max(ylim(ax)) - y_range * 0.08; 
                
                text(mid_time, y_text, strtrim(Waypoints(i).Label), ...
                    'HorizontalAlignment', 'center', 'FontSize', 9, 'FontWeight', 'bold', ...
                    'BackgroundColor', 'w', 'EdgeColor', 'k', 'Margin', 3);
            end
        end
    end
end

% Link the X-axes for zooming across all 6 kinematic subplots
linkaxes(findobj(gcf, 'Type', 'axes'), 'x');

%% --- Helper Functions ---
function out = correctDims(data, timeLen)
    sqData = squeeze(data);
    sz = size(sqData);
    if sz(1) == timeLen
        out = sqData;         
    elseif sz(2) == timeLen
        out = sqData';        
    else
        error('Data dimension mismatch: No dimension matches time vector length.');
    end
end

function plotProjection(ax, N, W, U, wp_coords, Waypoints, x_idx, y_idx, x_lbl, y_lbl)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    coords = [N, W, U];
    
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
    
    plot(ax, coords(:, x_idx), coords(:, y_idx), 'b', 'LineWidth', 1.5);
    
    plot(ax, wp_coords(x_idx, :), wp_coords(y_idx, :), 'rs', ...
         'MarkerFaceColor', 'r', 'MarkerSize', 6);
         
    xlabel(ax, x_lbl); ylabel(ax, y_lbl);
end