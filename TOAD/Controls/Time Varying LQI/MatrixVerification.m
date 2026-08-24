%% Load in parameters
if ~exist("constantsTOAD", "var")
    LoadTOADSim;
end
trajectoryName = "Circle_v1.csv";
% ReadGains now returns the cascaded gain set (4 outputs), not the old
% single [K_List, ~] pair.
[K_trans_List, K_rot_List, LA_List, LT_List] = ReadGains(trajectoryName);
m  = readmatrix(".\Guidance\Trajectories\" + trajectoryName);
dT = 1/1000; %m(2,1) - m(1,1);
x  = m(:, 2:16)';   % 15 x N   [q(4); r(3); v(3); omega(3); m_lox; m_ipa]
u  = m(:, 17:20)';  % 4  x N   [theta; phi; thrust; roll]
N = size(K_trans_List, 1);
assert(size(K_rot_List, 1) == N, ...
    'Translational/rotational gain horizons do not match -- check ReadGains output.');
assert(size(x, 2) == N, ...
    'Trajectory length does not match gain horizon -- check trajectoryName / CSV.');
%% Build Jacobian handles
% RicattiRecursion evaluates the nonlinear dynamics Jacobians via the
% generated JacobianX/JacobianU functions, called as f(x_n, u_n-1) (TWO
% arguments). This replaces the old locally-defined build_AB_jacobians(),
% which built its own symbolic Jacobian taking THREE arguments
% (x_n, x_n-1, u_n-1) -- a parallel derivation that is no longer guaranteed
% to match what actually generates the gains. Using the same function
% handles RicattiRecursion uses means this script's "reference" gain is
% produced by the exact call the production recursion makes.
%
% NOTE: if JacobianX/JacobianU are not on the path, this will error --
% that's intentional, since a stale/legacy symbolic Jacobian silently
% substituted back in would defeat the point of this check.
A_fcn = str2func('JacobianX');
B_fcn = str2func('JacobianU');
%% Reference (single endpoint) linearization + discretization
% Mirrors RicattiRecursion's per-step block exactly, evaluated once at the
% trajectory's final sample. This is a "did the recursion converge to
% something sane" sanity check, NOT a claim that one fixed-point
% linearization is valid across the whole (highly dynamic) maneuver -- the
% local, per-step closed-loop spectral radius further down is the
% maneuver-wide check.
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
% Translational model: fixed analytic double integrator, exactly as in
% RicattiRecursion (NOT taken from the linearization -- inertial-frame
% position/velocity kinematics are linear by construction).
A_trans_c = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
B_trans_c = [zeros(3,3); eye(3,3)];
M_c_trans     = [A_trans_c, B_trans_c; zeros(3,6), zeros(3,3)];
M_d_trans     = expm(M_c_trans * dT);
A_trans_d_ref = M_d_trans(1:6, 1:6);
B_trans_d_ref = M_d_trans(1:6, 7:9);
% Rotational model: attitude+rate rows/cols of the linearization; B uses
% only the theta/phi/roll columns of u (thrust, column 3, is excluded --
% it only actuates the translational loop).
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
error_norm_rot    = zeros(N,1);
K_fro_trans       = zeros(N,1);  % gain "magnitude" over the horizon, per loop
K_fro_rot         = zeros(N,1);
rho_trans         = nan(N,1);    % spectral radius, LOCAL linearization at each step
rho_rot           = nan(N,1);
bw_trans          = nan(N,1);    % approx closed-loop bandwidth [rad/s]
bw_rot            = nan(N,1);
bad_trans         = false(N,1);  % NaN/Inf flag per step, per loop
bad_rot           = false(N,1);
t_traj     = m(:, 1);
quat_traj  = x(1:4, :);
pos_traj   = x(5:7, :);
vel_traj   = x(8:10, :);
theta_cmd  = u(1, :);
phi_cmd    = u(2, :);
thrust_cmd = u(3, :);
R33 = quat_traj(1,:).^2 - quat_traj(2,:).^2 - quat_traj(3,:).^2 + quat_traj(4,:).^2;
tilt_angle = acosd(max(min(R33, 1), -1));

% Dynamically estimate LESO bandwidths from the gain matrices 
LA_mat0 = squeeze(LA_List(1,:,:));
LT_mat0 = squeeze(LT_List(1,:,:));
dT_flight = 1/1000;
omega_att_nom = max(abs(LA_mat0(:))) / (3 * dT_flight);
omega_thr_nom = max(abs(LT_mat0(:))) / (3 * dT_flight);

% The trajectory includes highly dynamic maneuvers, so a single fixed-point
% linearization is not representative of local closed-loop stability
% throughout the horizon. Re-linearize/re-discretize at every step exactly
% as RicattiRecursion does internally, then check each loop's gain against
% ITS OWN local A_d/B_d.
for n = 1:N
    Kt = squeeze(K_trans_List(n,:,:));   % 3x6
    Kr = squeeze(K_rot_List(n,:,:));     % 3x6
    bad_trans(n) = any(~isfinite(Kt(:)));
    bad_rot(n)   = any(~isfinite(Kr(:)));
    error_norm_trans(n) = norm(Kt - K_trans_ref, 'fro');
    error_norm_rot(n)   = norm(Kr - K_rot_ref,   'fro');
    K_fro_trans(n)       = norm(Kt, 'fro');
    K_fro_rot(n)         = norm(Kr, 'fro');
    if bad_trans(n) || bad_rot(n)
        continue % skip eig/bandwidth computations on garbage data
    end
    if n == 1
        continue % no x_{n-1} available, same gap as RicattiRecursion's loop bounds
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
    % Translational: local dT re-discretization of the fixed double
    % integrator (no re-linearization needed -- matches RicattiRecursion).
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
    eig_cl_t = eig(A_trans_d_local - B_trans_d_local * Kt);
    eig_cl_r = eig(A_rot_d_local   - B_rot_d_local   * Kr);
    rho_trans(n) = max(abs(eig_cl_t));
    rho_rot(n)   = max(abs(eig_cl_r));
    % Continuous-time equivalent bandwidth of each closed loop: map the
    % slowest (max-|z|) discrete pole back through s = log(z)/dT and take
    % its magnitude. Only meaningful while rho < 1.
    if rho_trans(n) < 1
        [~, idx_t] = max(abs(eig_cl_t));
        bw_trans(n) = abs(log(eig_cl_t(idx_t)) / dT_local);
    end
    if rho_rot(n) < 1
        [~, idx_r] = max(abs(eig_cl_r));
        bw_rot(n) = abs(log(eig_cl_r(idx_r)) / dT_local);
    end
end
K_trans_converged = squeeze(K_trans_List(1,:,:));
K_rot_converged    = squeeze(K_rot_List(1,:,:));
first_step_zero_trans = all(K_trans_converged(:) == 0);
first_step_zero_rot    = all(K_rot_converged(:) == 0);
n_unstable_trans = sum(rho_trans(~isnan(rho_trans)) >= 1);
n_unstable_rot    = sum(rho_rot(~isnan(rho_rot)) >= 1);
n_bad = sum(bad_trans | bad_rot);
%% Bandwidth Separation Diagnostics (cascade-specific)
% For a cascaded design to behave like its decoupled sub-loops, the inner
% (rotational) loop must be meaningfully faster than the outer
% (translational) loop it's nested inside, and the LESO observers must be
% faster still than the loop they're estimating disturbances for.
sep_rot_over_trans = bw_rot ./ bw_trans;              % inner/outer, want >> 1
sep_leso_att_over_rot = omega_att_nom ./ bw_rot;      % attitude LESO vs inner loop
sep_leso_thr_over_trans = omega_thr_nom ./ bw_trans;  % thrust LESO vs outer loop
SEP_THRESHOLD = 3; % conventional minimum for a valid time-scale-separated cascade
% LESO gains are placed against fixed kinematic observer models
% (A_LESO_A/C_LESO_A, A_LESO_T/C_LESO_T) and fixed target poles every
% iteration -- they do NOT depend on the trajectory linearization. So
% LA_List/LT_List should be (numerically) constant across the horizon.
% If they're not, that's a red flag that something trajectory-dependent
% leaked into the observer gain computation.
LA_flat = reshape(LA_List, N, []);
LT_flat = reshape(LT_List, N, []);
disp('--- Cascaded LQR/LESO Diagnostics ---');
fprintf('Horizon length (N): %d steps\n', N);
fprintf('Translational rho range: [%.4f, %.4f]\n', min(rho_trans,[],'omitnan'), max(rho_trans,[],'omitnan'));
fprintf('Rotational rho range:    [%.4f, %.4f]\n', min(rho_rot,[],'omitnan'),   max(rho_rot,[],'omitnan'));
if n_unstable_trans > 0
    fprintf(2, 'WARNING: %d/%d translational-loop steps have rho >= 1 against local linearization\n', n_unstable_trans, N);
else
    disp('Translational loop stabilizes its local linearization at every step (rho < 1).');
end
if n_unstable_rot > 0
    fprintf(2, 'WARNING: %d/%d rotational-loop steps have rho >= 1 against local linearization\n', n_unstable_rot, N);
else
    disp('Rotational loop stabilizes its local linearization at every step (rho < 1).');
end
if n_bad > 0
    fprintf(2, 'WARNING: %d/%d steps contain non-finite gain entries\n', n_bad, N);
end
if first_step_zero_trans || first_step_zero_rot
    fprintf(2, 'WARNING: first-step gain(s) are all-zero before the n=2 copy-down -- check RicattiRecursion loop bounds\n');
end
med_sep = median(sep_rot_over_trans, 'omitnan');
fprintf('Median inner/outer bandwidth separation (rot/trans): %.2fx\n', med_sep);
if med_sep < SEP_THRESHOLD
    fprintf(2, ['WARNING: inner/outer bandwidth separation < %gx -- cascade/', ...
        'singular-perturbation assumption is weak here\n'], SEP_THRESHOLD);
else
    fprintf('Inner/outer loop bandwidth separation exceeds the %gx threshold.\n', SEP_THRESHOLD);
end
disp('---------------------------------------');
%% Plots -- Translational (outer) loop
state_names_trans = {'x','y','z','v_x','v_y','v_z'};
row_names_trans    = {'a_x cmd', 'a_y cmd', 'a_z cmd'};
row_colors_trans   = lines(3);
figure('Name', 'Translational Loop Diagnostics', 'WindowStyle', 'docked');
tl1 = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl1, 'Outer (Translational) Loop Diagnostics');
n_top = 4;
for row = 1:3
    nexttile;
    hold on; grid on;
    [~, sorted_idx] = sort(abs(K_trans_ref(row, :)), 'descend');
    top_idx = sorted_idx(1:min(n_top, numel(sorted_idx)));
    for i = 1:length(top_idx)
        col = top_idx(i);
        component_trajectory = squeeze(K_trans_List(:, row, col));
        plot(N:-1:1, flipud(component_trajectory), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('K_{%d,%s}', row, state_names_trans{col}));
        yline(K_trans_ref(row, col), '--k', 'HandleVisibility', 'off');
    end
    set(gca, 'XDir', 'reverse');
    title(sprintf('Top Gains: %s', row_names_trans{row}));
    if row == 3
        xlabel('Backward Recursion Step (N \rightarrow 1)');
    end
    ylabel('Gain Value');
    legend('Location', 'best', 'FontSize', 7);
    hold off;
end
nexttile; hold on; grid on;
for row = 1:3
    row_norm = nan(N,1);
    for n = 1:N
        if ~bad_trans(n)
            row_norm(n) = norm(squeeze(K_trans_List(n, row, :)));
        end
    end
    plot(1:N, row_norm, 'LineWidth', 1.5, 'Color', row_colors_trans(row,:), ...
         'DisplayName', row_names_trans{row});
end
set(gca, 'YScale', 'log');
title('Per-Channel Gain Magnitude (log scale)');
xlabel('Horizon step n'); ylabel('|| K_{trans}(row,:) ||_2');
legend('Location', 'best', 'FontSize', 7);
hold off;
nexttile; hold on; grid on;
plot(1:N, rho_trans, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
yline(1, '--r', 'Marginally Stable', 'LabelHorizontalAlignment', 'left');
title('Closed-Loop Spectral Radius (local A_d, B_d per step)');
xlabel('Horizon step n'); ylabel('max |eig(A_d(n) - B_d(n) K_n)|');
ylim([0, max(1.1, max(rho_trans, [], 'omitnan')*1.1)]);
hold off;

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
%% Plots -- Rotational (inner) loop
state_names_rot = {'\delta\phi_1','\delta\phi_2','\delta\phi_3','\omega_1','\omega_2','\omega_3'};
row_names_rot    = {'\theta (gimbal 1)', '\phi (gimbal 2)', 'Roll'};
row_colors_rot   = lines(3);
figure('Name', 'Rotational Loop Diagnostics', 'WindowStyle', 'docked');
tl2 = tiledlayout(2, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl2, 'Inner (Rotational) Loop Diagnostics');
for row = 1:3
    nexttile;
    hold on; grid on;
    [~, sorted_idx] = sort(abs(K_rot_ref(row, :)), 'descend');
    top_idx = sorted_idx(1:min(n_top, numel(sorted_idx)));
    for i = 1:length(top_idx)
        col = top_idx(i);
        component_trajectory = squeeze(K_rot_List(:, row, col));
        plot(N:-1:1, flipud(component_trajectory), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('K_{%d,%s}', row, state_names_rot{col}));
        yline(K_rot_ref(row, col), '--k', 'HandleVisibility', 'off');
    end
    set(gca, 'XDir', 'reverse');
    title(sprintf('Top Gains: %s', row_names_rot{row}));
    if row == 3
        xlabel('Backward Recursion Step (N \rightarrow 1)');
    end
    ylabel('Gain Value');
    legend('Location', 'best', 'FontSize', 7);
    hold off;
end
nexttile; hold on; grid on;
for row = 1:3
    row_norm = nan(N,1);
    for n = 1:N
        if ~bad_rot(n)
            row_norm(n) = norm(squeeze(K_rot_List(n, row, :)));
        end
    end
    plot(1:N, row_norm, 'LineWidth', 1.5, 'Color', row_colors_rot(row,:), ...
         'DisplayName', row_names_rot{row});
end
set(gca, 'YScale', 'log');
title('Per-Channel Gain Magnitude (log scale)');
xlabel('Horizon step n'); ylabel('|| K_{rot}(row,:) ||_2');
legend('Location', 'best', 'FontSize', 7);
hold off;
nexttile; hold on; grid on;
plot(1:N, rho_rot, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
yline(1, '--r', 'Marginally Stable', 'LabelHorizontalAlignment', 'left');
title('Closed-Loop Spectral Radius (local A_d, B_d per step)');
xlabel('Horizon step n'); ylabel('max |eig(A_d(n) - B_d(n) K_n)|');
ylim([0, max(1.1, max(rho_rot, [], 'omitnan')*1.1)]);
hold off;

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
%% Plots -- Bandwidth Separation (cascade validity)
figure('Name', 'Cascade Bandwidth Separation', 'WindowStyle', 'docked');
tl3 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl3, 'Time-Scale Separation Across the Cascade');
nexttile; hold on; grid on;
plot(1:N, bw_trans, 'LineWidth', 1.5, 'DisplayName', 'Translational (outer)');
plot(1:N, bw_rot, 'LineWidth', 1.5, 'DisplayName', 'Rotational (inner)');
yline(omega_att_nom, '--', 'Attitude LESO nominal \omega', 'Color', [0.6 0.2 0.6], 'LabelHorizontalAlignment', 'left');
yline(omega_thr_nom, ':', 'Thrust LESO nominal \omega', 'Color', [0.2 0.4 0.6], 'LabelHorizontalAlignment', 'left');
ylim([1e-1, 2e2]);
set(gca, 'YScale', 'log');
title('Closed-Loop Bandwidth Estimate');
xlabel('Horizon step n'); ylabel('\omega_{cl} [rad/s] (log)');
legend('Location', 'best', 'FontSize', 7);
hold off;
nexttile; hold on; grid on;
plot(1:N, sep_rot_over_trans, 'LineWidth', 1.5, 'Color', [0 0.4 0.8]);
yline(SEP_THRESHOLD, '--r', sprintf('%gx min separation', SEP_THRESHOLD), 'LabelHorizontalAlignment', 'left');
title('Inner/Outer Bandwidth Separation');
xlabel('Horizon step n'); ylabel('\omega_{rot} / \omega_{trans}');
hold off;
nexttile; hold on; grid on;
LA_norm = vecnorm(LA_flat, 2, 2);
LT_norm = vecnorm(LT_flat, 2, 2);
plot(1:N, LA_norm, 'LineWidth', 1.5, 'DisplayName', 'Attitude LESO ||L_{att}||_F');
plot(1:N, LT_norm, 'LineWidth', 1.5, 'DisplayName', 'Thrust LESO ||L_{thr}||_F');
title('LESO Gain Consistency (should be flat)');
xlabel('Horizon step n'); ylabel('Frobenius norm');
legend('Location', 'best', 'FontSize', 7);
hold off;
nexttile; hold on; grid on;
plot(1:N, error_norm_trans, 'LineWidth', 1.5, 'DisplayName', 'Translational');
plot(1:N, error_norm_rot, 'LineWidth', 1.5, 'DisplayName', 'Rotational');
title('Gain Error vs Single-Point DARE Reference');
xlabel('Horizon step n'); ylabel('|| K_n - K_{ref} ||_F');
legend('Location', 'best', 'FontSize', 7);
hold off;
%% Trajectory Context Figure (ported from TrajectoryGenerator.m, unchanged)
figure('Name', 'Trajectory: Control & Attitude Context', 'WindowStyle', 'docked');
subplot(3,1,1); hold on; grid on;
yyaxis left;
plot(t_traj, quat_traj(1,:), 'k', 'LineWidth', 1.5, 'DisplayName', 'q_w (Scalar)');
plot(t_traj, quat_traj(2,:), 'r', 'LineWidth', 1.5, 'DisplayName', 'q_x');
plot(t_traj, quat_traj(3,:), 'g', 'LineWidth', 1.5, 'DisplayName', 'q_y');
plot(t_traj, quat_traj(4,:), 'b', 'LineWidth', 1.5, 'DisplayName', 'q_z');
ylabel('Quaternion Component');
ylim([-1.1 1.1]);
yyaxis right;
plot(t_traj, tilt_angle, '--m', 'LineWidth', 2, 'DisplayName', 'Tilt Angle (deg)');
ylabel('Tilt Angle [deg]');
ylim([0 185]);
legend('Location', 'eastoutside');
title('Vehicle Attitude (Quaternions & Absolute Tilt)');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'm';
subplot(3,1,2); hold on; grid on;
plot(t_traj, rad2deg(theta_cmd), 'Color', '#7E2F8E', 'LineWidth', 1.5, 'DisplayName', '\theta');
plot(t_traj, rad2deg(phi_cmd), 'Color', '#77AC30', 'LineWidth', 1.5, 'DisplayName', '\phi');
ylabel('Gimbal [deg]'); legend('Location', 'best');
title('Gimbal Deflection Commands');
subplot(3,1,3); hold on; grid on;
area(t_traj, thrust_cmd, 'FaceColor', '#EDB120', 'FaceAlpha', 0.4, 'EdgeColor', '#A27712', 'LineWidth', 1.5);
ylabel('Thrust [N]'); xlabel('Time [s]');
title('Thrust Command');
%% Trajectory Position/Velocity & 3D Profile (ported from TrajectoryGenerator.m, unchanged)
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
% linkaxes(findobj(gcf, 'Type', 'axes'), 'x');
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
    % Simple diverging colormap so zero-centered gain heatmaps read
    % cleanly (blue = negative, white = zero, red = positive)
    n = 128;
    neg = [linspace(0,1,n)', linspace(0,1,n)', ones(n,1)];        % blue -> white
    pos = [ones(n,1), linspace(1,0,n)', linspace(1,0,n)'];        % white -> red
    cmap = [neg; pos];
end