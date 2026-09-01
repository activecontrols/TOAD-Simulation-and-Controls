%% Load in parameters
if ~exist("constantsTOAD", "var")
    LoadTOADSim;
end
trajectoryName = "Backflip_v1.csv";
% ReadGains now returns the cascaded gain set (4 outputs), not the old
% single [K_List, ~] pair.
[K_trans_List, K_rot_List, LA_List, LT_List] = ReadGains(trajectoryName);
m  = readmatrix(".\Guidance\Trajectories\" + trajectoryName);
dT = 1/500; %m(2,1) - m(1,1);
x  = m(:, 2:16)';   % 15 x N   [q(4); r(3); v(3); omega(3); m_lox; m_ipa]
u  = m(:, 17:20)';  % 4  x N   [theta; phi; thrust; roll]
N = size(K_trans_List, 1);
assert(size(K_rot_List, 1) == N, ...
    'Translational/rotational gain horizons do not match -- check ReadGains output.');
assert(size(x, 2) == N, ...
    'Trajectory length does not match gain horizon -- check trajectoryName / CSV.');

%% Build Jacobian handles
% Uses the same JacobianX/JacobianU functions RicattiRecursion calls.
A_fcn = str2func('JacobianX');
B_fcn = str2func('JacobianU');

%% Reference (single endpoint) linearization + discretization
% Sanity check only -- the per-step spectral radius below is the maneuver-wide check.
x_n  = x(:, end);
u_n1 = u(:, end-1);
A_lin = A_fcn(x_n, u_n1);
B_lin = B_fcn(x_n, u_n1);

q_ref = x_n(1:4);
T = zeros(15, 12);
T(1:4, 1:3)   = 0.5 * XiMat(q_ref);
T(5:13, 4:12) = eye(9);
A_lin = pinv(T) * A_lin * T;
B_lin = pinv(T) * B_lin;

% Translational model: fixed analytic double integrator (matches RicattiRecursion).
A_trans_c = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
B_trans_c = [zeros(3,3); eye(3,3)];
M_c_trans     = [A_trans_c, B_trans_c; zeros(3,6), zeros(3,3)];
M_d_trans     = expm(M_c_trans * dT);
A_trans_d_ref = M_d_trans(1:6, 1:6);
B_trans_d_ref = M_d_trans(1:6, 7:9);

% Rotational model: attitude+rate rows/cols; B uses theta/phi/roll cols only (thrust excluded).
A_rot_c = [A_lin(1:3,1:3),   A_lin(1:3,10:12);
           A_lin(10:12,1:3), A_lin(10:12,10:12)];
B_rot_c = [B_lin(1:3,  [1,2,4]);
           B_lin(10:12,[1,2,4])];
M_c_rot     = [A_rot_c, B_rot_c; zeros(3,6), zeros(3,3)];
M_d_rot     = expm(M_c_rot * dT);
A_rot_d_ref = M_d_rot(1:6, 1:6);
B_rot_d_ref = M_d_rot(1:6, 7:9);

Q_trans = constantsTOAD.Q_trans; R_trans = constantsTOAD.R_trans;
Q_rot   = constantsTOAD.Q_rot;   R_rot   = constantsTOAD.R_rot;

P_trans_ref = idare(A_trans_d_ref, B_trans_d_ref, Q_trans, R_trans);
P_rot_ref   = idare(A_rot_d_ref,   B_rot_d_ref,   Q_rot,   R_rot);
K_trans_ref = (B_trans_d_ref'*P_trans_ref*B_trans_d_ref + R_trans) \ (B_trans_d_ref'*P_trans_ref*A_trans_d_ref);
K_rot_ref   = (B_rot_d_ref'*P_rot_ref*B_rot_d_ref     + R_rot)   \ (B_rot_d_ref'*P_rot_ref*A_rot_d_ref);

%% Diagnostics
error_norm_trans = zeros(N,1);   % Frobenius error vs single-point DARE, per loop
error_norm_rot   = zeros(N,1);
K_fro_trans      = zeros(N,1);  % gain "magnitude" over the horizon, per loop
K_fro_rot        = zeros(N,1);
rho_trans        = nan(N,1);    % spectral radius, LOCAL linearization at each step
rho_rot          = nan(N,1);
bad_trans        = false(N,1);  % NaN/Inf flag per step, per loop
bad_rot          = false(N,1);

t_traj     = m(:, 1);
pos_traj   = x(5:7, :);
vel_traj   = x(8:10, :);

% LESO Bandwidth (using explicit constants to avoid reverse-calculation flaws)
omega_att_nom = constantsTOAD.OmegaAtt;
omega_thr_nom = constantsTOAD.OmegaThr;

% Crossover tracking (gain crossover frequency per actuator/channel)
wc_trans = nan(N, 3); % [X, Y, Z]
wc_rot   = nan(N, 3); % [Theta, Phi, Roll]

% Turn off bode/margin warnings for edge cases
warning('off', 'Control:analysis:MarginUnstable');

% Re-linearize/re-discretize at every step, exactly as RicattiRecursion does internally.
for n = 1:N
    Kt = squeeze(K_trans_List(n,:,:));   % 3x6
    Kr = squeeze(K_rot_List(n,:,:));     % 3x6
    bad_trans(n) = any(~isfinite(Kt(:)));
    bad_rot(n)   = any(~isfinite(Kr(:)));
    error_norm_trans(n) = norm(Kt - K_trans_ref, 'fro');
    error_norm_rot(n)   = norm(Kr - K_rot_ref,   'fro');
    K_fro_trans(n)      = norm(Kt, 'fro');
    K_fro_rot(n)        = norm(Kr, 'fro');
    
    if bad_trans(n) || bad_rot(n)
        continue % skip eig/bandwidth computations on garbage data
    end
    if n == 1
        continue % no x_{n-1} available
    end
    
    x_n_local  = x(:, n);
    u_n1_local = u(:, n-1);
    dT_local   = m(n,1) - m(n-1,1);
    A_lin_local = A_fcn(x_n_local, u_n1_local);
    B_lin_local = B_fcn(x_n_local, u_n1_local);
    
    q_local = x_n_local(1:4);
    T_local = zeros(15, 12);
    T_local(1:4, 1:3)   = 0.5 * XiMat(q_local);
    T_local(5:13, 4:12) = eye(9);
    A_lin_local = pinv(T_local) * A_lin_local * T_local;
    B_lin_local = pinv(T_local) * B_lin_local;
    
    % Translational: local dT re-discretization of the fixed double integrator.
    Mct = [A_trans_c, B_trans_c; zeros(3,6), zeros(3,3)];
    Mdt = expm(Mct * dT_local);
    A_trans_d_local = Mdt(1:6, 1:6);
    B_trans_d_local = Mdt(1:6, 7:9);
    
    % Rotational: local re-linearization + re-discretization.
    A_rot_c_local = [A_lin_local(1:3,1:3),   A_lin_local(1:3,10:12);
                     A_lin_local(10:12,1:3), A_lin_local(10:12,10:12)];
    B_rot_c_local = [B_lin_local(1:3,  [1,2,4]);
                     B_lin_local(10:12,[1,2,4])];
    Mcr = [A_rot_c_local, B_rot_c_local; zeros(3,6), zeros(3,3)];
    Mdr = expm(Mcr * dT_local);
    A_rot_d_local = Mdr(1:6, 1:6);
    B_rot_d_local = Mdr(1:6, 7:9);
    
    % Check Strict Stability
    eig_cl_t = eig(A_trans_d_local - B_trans_d_local * Kt);
    eig_cl_r = eig(A_rot_d_local   - B_rot_d_local   * Kr);
    rho_trans(n) = max(abs(eig_cl_t));
    rho_rot(n)   = max(abs(eig_cl_r));

    % Crossover Frequency via channel-specific loop transfer function L = K(sI-A)^-1B
    for ch = 1:3
        % Translational loop (1=X, 2=Y, 3=Z)
        % Break the loop at the specific actuator input
        sys_t = ss(A_trans_d_local, B_trans_d_local(:, ch), Kt(ch, :), 0, dT_local);
        [~, ~, ~, wcp_t] = margin(sys_t);
        if ~isempty(wcp_t) && isfinite(wcp_t(1))
            wc_trans(n, ch) = wcp_t(1);
        end
        
        % Rotational loop (1=Theta, 2=Phi, 3=Roll)
        sys_r = ss(A_rot_d_local, B_rot_d_local(:, ch), Kr(ch, :), 0, dT_local);
        [~, ~, ~, wcp_r] = margin(sys_r);
        if ~isempty(wcp_r) && isfinite(wcp_r(1))
            wc_rot(n, ch) = wcp_r(1);
        end
    end
end
warning('on', 'Control:analysis:MarginUnstable'); % Restore warnings

K_trans_converged = squeeze(K_trans_List(1,:,:));
K_rot_converged   = squeeze(K_rot_List(1,:,:));
first_step_zero_trans = all(K_trans_converged(:) == 0);
first_step_zero_rot   = all(K_rot_converged(:) == 0);

n_unstable_trans = sum(rho_trans(~isnan(rho_trans)) >= 1);
n_unstable_rot   = sum(rho_rot(~isnan(rho_rot)) >= 1);
n_bad = sum(bad_trans | bad_rot);

%% Bandwidth Separation Diagnostics (cascade-specific)
SEP_THRESHOLD = 3; % conventional minimum for a valid time-scale-separated cascade
disp('--- Cascaded LQR/LESO Diagnostics ---');
fprintf('Horizon length (N): %d steps\n', N);
fprintf('Translational rho range: [%.4f, %.4f]\n', min(rho_trans,[],'omitnan'), max(rho_trans,[],'omitnan'));
fprintf('Rotational rho range:    [%.4f, %.4f]\n', min(rho_rot,[],'omitnan'),   max(rho_rot,[],'omitnan'));

if n_unstable_trans > 0
    fprintf(2, 'WARNING: %d/%d translational-loop steps have rho >= 1 against local linearization\n', n_unstable_trans, N);
end
if n_unstable_rot > 0
    fprintf(2, 'WARNING: %d/%d rotational-loop steps have rho >= 1 against local linearization\n', n_unstable_rot, N);
end
if n_bad > 0
    fprintf(2, 'WARNING: %d/%d steps contain non-finite gain entries\n', n_bad, N);
end

% Median Separation Checks
med_wc_trans = median(wc_trans, 1, 'omitnan');
med_wc_rot   = median(wc_rot, 1, 'omitnan');

% X -> Theta & Thrust LESO
sep_x_theta = med_wc_rot(1) / med_wc_trans(1);
sep_x_tleso = omega_thr_nom / med_wc_trans(1);
fprintf('Median bandwidth separation (X-axis -> Theta): %.2fx\n', sep_x_theta);
if sep_x_theta < SEP_THRESHOLD, fprintf(2, '  -> WARNING: X/Theta separation < %gx\n', SEP_THRESHOLD); end
fprintf('Median bandwidth separation (X-axis -> Thrust LESO): %.2fx\n', sep_x_tleso);
if sep_x_tleso < SEP_THRESHOLD, fprintf(2, '  -> WARNING: X/Thrust LESO separation < %gx\n', SEP_THRESHOLD); end

% Y -> Phi & Thrust LESO
sep_y_phi = med_wc_rot(2) / med_wc_trans(2);
sep_y_tleso = omega_thr_nom / med_wc_trans(2);
fprintf('Median bandwidth separation (Y-axis -> Phi): %.2fx\n', sep_y_phi);
if sep_y_phi < SEP_THRESHOLD, fprintf(2, '  -> WARNING: Y/Phi separation < %gx\n', SEP_THRESHOLD); end
fprintf('Median bandwidth separation (Y-axis -> Thrust LESO): %.2fx\n', sep_y_tleso);
if sep_y_tleso < SEP_THRESHOLD, fprintf(2, '  -> WARNING: Y/Thrust LESO separation < %gx\n', SEP_THRESHOLD); end

% Z -> Thrust LESO
sep_z_leso = omega_thr_nom / med_wc_trans(3);
fprintf('Median bandwidth separation (Z-axis -> Thrust LESO): %.2fx\n', sep_z_leso);
if sep_z_leso < SEP_THRESHOLD, fprintf(2, '  -> WARNING: Z/Thrust LESO separation < %gx\n', SEP_THRESHOLD); end

disp('---------------------------------------');

%% Plots -- Translational (outer) loop gain heatmaps
state_names_trans = {'x','y','z','v_x','v_y','v_z'};
figure('Name', 'Translational Loop Diagnostics', 'WindowStyle', 'docked');
tl1 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Outer (Translational) Loop Gain Heatmaps');

nexttile; hold on; grid on;
K_slice_t1 = squeeze(K_trans_List(:,1,:));
imagesc(K_slice_t1); axis tight;
clim_val_t1 = max(abs(K_slice_t1(:)), [], 'omitnan');
if clim_val_t1 > 0, clim([-clim_val_t1, clim_val_t1]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_trans);
title('Gain Heatmap: a_x cmd'); xlabel('State'); ylabel('Horizon step n');

nexttile; hold on; grid on;
K_slice_t2 = squeeze(K_trans_List(:,2,:));
imagesc(K_slice_t2); axis tight;
clim_val_t2 = max(abs(K_slice_t2(:)), [], 'omitnan');
if clim_val_t2 > 0, clim([-clim_val_t2, clim_val_t2]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_trans);
title('Gain Heatmap: a_y cmd'); xlabel('State'); ylabel('Horizon step n');

nexttile; hold on; grid on;
K_slice_t3 = squeeze(K_trans_List(:,3,:));
imagesc(K_slice_t3); axis tight;
clim_val_t3 = max(abs(K_slice_t3(:)), [], 'omitnan');
if clim_val_t3 > 0, clim([-clim_val_t3, clim_val_t3]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_trans);
title('Gain Heatmap: a_z cmd'); xlabel('State'); ylabel('Horizon step n');

%% Plots -- Rotational (inner) loop gain heatmaps
state_names_rot = {'\delta\phi_1','\delta\phi_2','\delta\phi_3','\omega_1','\omega_2','\omega_3'};
figure('Name', 'Rotational Loop Diagnostics', 'WindowStyle', 'docked');
tl2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Inner (Rotational) Loop Gain Heatmaps');

nexttile; hold on; grid on;
K_slice_r1 = squeeze(K_rot_List(:,1,:));
imagesc(K_slice_r1); axis tight;
clim_val_r1 = max(abs(K_slice_r1(:)), [], 'omitnan');
if clim_val_r1 > 0, clim([-clim_val_r1, clim_val_r1]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_rot);
title('Gain Heatmap: \theta'); xlabel('State'); ylabel('Horizon step n');

nexttile; hold on; grid on;
K_slice_r2 = squeeze(K_rot_List(:,2,:));
imagesc(K_slice_r2); axis tight;
clim_val_r2 = max(abs(K_slice_r2(:)), [], 'omitnan');
if clim_val_r2 > 0, clim([-clim_val_r2, clim_val_r2]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_rot);
title('Gain Heatmap: \phi'); xlabel('State'); ylabel('Horizon step n');

nexttile; hold on; grid on;
K_slice_r3 = squeeze(K_rot_List(:,3,:));
imagesc(K_slice_r3); axis tight;
clim_val_r3 = max(abs(K_slice_r3(:)), [], 'omitnan');
if clim_val_r3 > 0, clim([-clim_val_r3, clim_val_r3]); end
colormap(gca, redblue_local()); colorbar;
set(gca, 'XTick', 1:6, 'XTickLabel', state_names_rot);
title('Gain Heatmap: Roll'); xlabel('State'); ylabel('Horizon step n');

%% Plot -- Actuator-Specific Crossover Separation
figure('Name', 'Actuator-Specific Crossover', 'WindowStyle', 'docked');
tl_cross = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_cross, 'Actuator-Specific Loop Crossover Across Horizon');

% Tile 1: X -> Theta separation
nexttile; hold on; grid on;
plot(1:N, wc_trans(:, 1), 'LineWidth', 2, 'DisplayName', 'Outer: X-Axis');
plot(1:N, wc_rot(:, 1), 'LineWidth', 2, 'DisplayName', 'Inner: Pitch (\theta)');
yline(omega_att_nom, '--', 'Attitude LESO', 'Color', [0.6 0.2 0.6], 'LineWidth', 1.5, 'DisplayName', 'Attitude LESO');
yline(omega_thr_nom, ':', 'Thrust LESO', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.5, 'DisplayName', 'Thrust LESO');
set(gca, 'YScale', 'log'); ylim([1e-1, 2e2]);
title('X-Axis / Pitch (\theta) Channel'); xlabel('Horizon step n'); ylabel('\omega_c [rad/s]'); legend('Location', 'best', 'FontSize', 8);

% Tile 2: Y -> Phi separation
nexttile; hold on; grid on;
plot(1:N, wc_trans(:, 2), 'LineWidth', 2, 'DisplayName', 'Outer: Y-Axis');
plot(1:N, wc_rot(:, 2), 'LineWidth', 2, 'DisplayName', 'Inner: Yaw (\phi)');
yline(omega_att_nom, '--', 'Attitude LESO', 'Color', [0.6 0.2 0.6], 'LineWidth', 1.5, 'DisplayName', 'Attitude LESO');
yline(omega_thr_nom, ':', 'Thrust LESO', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.5, 'DisplayName', 'Thrust LESO');
set(gca, 'YScale', 'log'); ylim([1e-1, 2e2]);
title('Y-Axis / Yaw (\phi) Channel'); xlabel('Horizon step n'); ylabel('\omega_c [rad/s]'); legend('Location', 'best', 'FontSize', 8);

% Tile 3: Z -> Thrust separation
nexttile; hold on; grid on;
plot(1:N, wc_trans(:, 3), 'LineWidth', 2, 'DisplayName', 'Outer: Z-Axis');
yline(omega_thr_nom, ':', 'Thrust LESO', 'Color', [0.2 0.4 0.6], 'LineWidth', 1.5, 'DisplayName', 'Thrust LESO');
set(gca, 'YScale', 'log'); ylim([1e-1, 2e2]);
title('Z-Axis / Thrust Channel'); xlabel('Horizon step n'); ylabel('\omega_c [rad/s]'); legend('Location', 'best', 'FontSize', 8);

% Tile 4: Roll separation
nexttile; hold on; grid on;
plot(1:N, wc_rot(:, 3), 'LineWidth', 2, 'DisplayName', 'Inner: Roll');
yline(omega_att_nom, '--', 'Attitude LESO', 'Color', [0.6 0.2 0.6], 'LineWidth', 1.5, 'DisplayName', 'Attitude LESO');
set(gca, 'YScale', 'log'); ylim([1e-1, 2e2]);
title('Roll Channel'); xlabel('Horizon step n'); ylabel('\omega_c [rad/s]'); legend('Location', 'best', 'FontSize', 8);

%% Trajectory Position/Velocity & 3D Profile (unchanged)
figure('Name', 'Trajectory: Mission Kinematics', 'WindowStyle', 'docked');
tl_kin = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl_kin, 'Position and Velocity Profiles', 'FontWeight', 'bold');
xlabel(tl_kin, 'Time [s]', 'FontWeight', 'bold');
labels_pos = {'X Pos [m]', 'Y Pos [m]', 'Z Pos [m]'};
labels_vel = {'X Vel [m/s]', 'Y Vel [m/s]', 'Z Vel [m/s]'};
kin_colors = {'#0072BD', '#D95319', '#EDB120'};
for ax_idx = 1:3
    nexttile(tl_kin, (ax_idx-1)*2 + 1); hold on; grid on;
    plot(t_traj, pos_traj(ax_idx, :), 'Color', kin_colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_pos{ax_idx}, 'FontWeight', 'bold');
    nexttile(tl_kin, (ax_idx-1)*2 + 2); hold on; grid on;
    plot(t_traj, vel_traj(ax_idx, :), 'Color', kin_colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_vel{ax_idx}, 'FontWeight', 'bold');
end

figure('Name', 'Trajectory: 3D Mission Profile', 'WindowStyle', 'docked');
tl_3d = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
axMain = nexttile(tl_3d, 1, [3 3]);
hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
xlabel(axMain, 'X [m]'); ylabel(axMain, 'Y [m]'); zlabel(axMain, 'Z (Alt) [m]');
title(axMain, '3D Mission Trajectory');
patch(axMain, [pos_traj(1,:), NaN], [pos_traj(2,:), NaN], [pos_traj(3,:), NaN], [t_traj', NaN], ...
      'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5);
cb = colorbar(axMain); cb.Label.String = 'Time [s]'; colormap(axMain, 'turbo');
plot3(axMain, pos_traj(1,1), pos_traj(2,1), pos_traj(3,1), 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8, 'DisplayName', 'Start');
plot3(axMain, pos_traj(1,end), pos_traj(2,end), pos_traj(3,end), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8, 'DisplayName', 'End');
axTop = nexttile(tl_3d, 4);
hold(axTop, 'on'); grid(axTop, 'on'); axis(axTop, 'equal');
plot(axTop, pos_traj(1,:), pos_traj(2,:), 'b', 'LineWidth', 1.5);
xlabel(axTop, 'X [m]'); ylabel(axTop, 'Y [m]'); title(axTop, 'Top View (X-Y)');
axSide = nexttile(tl_3d, 8);
hold(axSide, 'on'); grid(axSide, 'on'); axis(axSide, 'equal');
plot(axSide, pos_traj(1,:), pos_traj(3,:), 'b', 'LineWidth', 1.5);
xlabel(axSide, 'X [m]'); ylabel(axSide, 'Z (Alt) [m]'); title(axSide, 'Side View (X-Z)');
axFront = nexttile(tl_3d, 12);
hold(axFront, 'on'); grid(axFront, 'on'); axis(axFront, 'equal');
plot(axFront, pos_traj(2,:), pos_traj(3,:), 'b', 'LineWidth', 1.5);
xlabel(axFront, 'Y [m]'); ylabel(axFront, 'Z (Alt) [m]'); title(axFront, 'Front View (Y-Z)');

function cmap = redblue_local()
    n = 128;
    neg = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];
    pos = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];
    cmap = [neg; pos];
end