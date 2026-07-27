%% Test Script for TVLQR Gain Generation
% Call and load in parameters
if ~exist("constantsTOAD")
    LoadTOADSim;
end

% Hand tuning for Q for now
a_weights = ones(12,1);
a_weights = a_weights / norm(a_weights);
max_x = [0.23, 0.23, 0.04, 4, 4, 4, 5, 5, 5, 0.72, 0.72, 1];
Q = eye(12) .* a_weights ./ max_x.^2;
R = diag([3.5, 3.5, 1/2400^2, 5]);

K_List=ReadGains("GainMatrixTrajectory1.csv");
m = readmatrix(".\Guidance\Trajectories\Trajectory1.csv");
dT = m(2,1)-m(1,1);
x = m(:, 2:16)';
u = m(:, 17:20)';

% Build the A/B Jacobian function handles once (expensive symbolic step)
% and reuse them both for the single-point DARE reference and for the
% per-step local linearization used below. This replaces the old
% generateAB() which only ever evaluated at the trajectory's last sample.
[A_fcn, B_fcn] = build_AB_jacobians(constantsTOAD);

x_n   = x(:, end);
x_n1  = x(:, end-1);
u_n1  = u(:, end-1);
A_lin = A_fcn(x_n, x_n1, u_n1);
B_lin = B_fcn(x_n, x_n1, u_n1);

% Kinematic Maps
q_ref = x_n(1:4);
T = zeros(15, 12);
T(1:4, 1:3) = 0.5 * XiMat(q_ref); 
T(5:13, 4:12) = eye(9);
A_lin = pinv(T) * A_lin * T;
B_lin = pinv(T) * B_lin;

SanityCheck = lqrd(A_lin, B_lin, Q, R, dT);

%% Diagnostics
N = size(K_List, 1);
nx = 12; nu = 4;
error_norm      = zeros(N, 1);   % Frobenius error vs single-point DARE
K_fro_norm      = zeros(N, 1);   % gain "magnitude" over the horizon
closed_loop_rho = nan(N, 1);     % spectral radius of A_d(n) - B_d(n)*K_t, LOCAL linearization at each step
has_bad_values  = false(N, 1);   % NaN/Inf flag per step

% Trajectory context (ported from TrajectoryGenerator.m) -- pulled from
% the same x/u/m arrays already loaded above, no re-solve needed. Used
% for the standalone context figures below.
t_traj     = m(:, 1);
quat_traj  = x(1:4, :);
pos_traj   = x(5:7, :);
vel_traj   = x(8:10, :);
theta_cmd  = u(1, :);
phi_cmd    = u(2, :);
thrust_cmd = u(3, :);
R33 = quat_traj(1,:).^2 - quat_traj(2,:).^2 - quat_traj(3,:).^2 + quat_traj(4,:).^2;
tilt_angle = acosd(max(min(R33, 1), -1));

% The trajectory includes highly dynamic maneuvers, so a single
% fixed-point linearization is not representative of local closed-loop
% stability throughout the horizon. Instead, re-linearize and
% re-discretize (ZOH) at every step exactly as RicattiRecursion does
% internally (using x_n, x_n-1, u_n-1 and the local dT), then check each
% gain against ITS OWN local A_d/B_d rather than the endpoint's.
for n = 1:N
    K_t = squeeze(K_List(n, :, :))';

    has_bad_values(n) = any(~isfinite(K_t(:)));
    if has_bad_values(n)
        continue % skip norm/eig computations on garbage data
    end

    if n == 1
        continue % no x_{n-1} available at the first sample, same gap as RicattiRecursion's loop bounds
    end

    x_n_local  = x(:, n);
    x_n1_local = x(:, n-1);
    u_n1_local = u(:, n-1);
    dT_local   = m(n,1) - m(n-1,1);

    A_lin_local = A_fcn(x_n_local, x_n1_local, u_n1_local);
    B_lin_local = B_fcn(x_n_local, x_n1_local, u_n1_local);

    q_ref_local = x_n_local(1:4);
    T_local = zeros(15, 12);
    T_local(1:4, 1:3) = 0.5 * XiMat(q_ref_local);
    T_local(5:13, 4:12) = eye(9);

    A_lin_local = pinv(T_local) * A_lin_local * T_local;
    B_lin_local = pinv(T_local) * B_lin_local;

    M_c_local = [A_lin_local, B_lin_local; zeros(nu, nx), zeros(nu, nu)];
    M_d_local = expm(M_c_local * dT_local);
    A_d_local = M_d_local(1:nx, 1:nx);
    B_d_local = M_d_local(1:nx, (nx+1):end);

    eig_cl = eig(A_d_local - B_d_local*K_t(:,1:12));
    closed_loop_rho(n) = max(abs(eig_cl));
end

K_converged = squeeze(K_List(1, :, :))';

% Row-by-row: is the all-zero first entry (see RicattiRecursion loop
% bounds n = N:-1:2) actually populated?
first_step_is_zero = all(K_converged(:) == 0);

n_unstable = sum(closed_loop_rho(~isnan(closed_loop_rho)) >= 1);
n_bad      = sum(has_bad_values);

% Command Window Summary (kept in the console only -- see plotting
% notes below for why the on-figure text/DARE-norm panels were dropped)
disp('--- LQR Recursion Diagnostics ---');
disp(['Horizon Length (N): ', num2str(N), ' steps']);
disp(['Closed-loop spectral radius range: [', num2str(min(closed_loop_rho)), ', ', num2str(max(closed_loop_rho)), ']']);
if n_unstable > 0
    fprintf(2, 'WARNING: %d of %d gains give spectral radius >= 1 against the reference linearization (possible instability)\n', n_unstable, N);
else
    disp('All gains stabilize the reference linearization (spectral radius < 1).');
end
disp('---------------------------------');

%% Plots
row_names = {'\theta (gimbal 1)', '\phi (gimbal 2)', 'Thrust', 'Roll'};
row_colors = lines(4);

figure('Name', 'Riccati Recursion Diagnostics', 'WindowStyle', 'docked');
tl = tiledlayout(3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, 'TVLQR Riccati Recursion Diagnostics');

% (Col 1, rows 1-3) Top gain components per control row: theta, phi, thrust
n_top = 4; % how many state-columns to show per control row
for row = 1:3
    nexttile;
    hold on; grid on;
    [~, sorted_idx] = sort(abs(SanityCheck(row, :)), 'descend');
    top_idx = sorted_idx(1:min(n_top, numel(sorted_idx)));
    for i = 1:length(top_idx)
        col = top_idx(i);
        component_trajectory = squeeze(K_List(:, col, row));
        plot(N:-1:1, flipud(component_trajectory), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('K_{%d,%d}', row, col));
        yline(SanityCheck(row, col), '--k', 'HandleVisibility', 'off');
    end
    set(gca, 'XDir', 'reverse');
    title(sprintf('Top Gains: %s', row_names{row}));
    if row == 3
        xlabel('Backward Recursion Step (N \rightarrow 1)');
    end
    ylabel('Gain Value');
    legend('Location', 'best', 'FontSize', 7);
    hold off;
end

% (Col 2, rows 1-3) Per-channel gain heatmap (state vs horizon step),
% each with its own color scale -- avoids thrust magnitude swamping the
% gimbal channels, and makes the sparsity pattern within each channel
% actually visible.
for row = 1:3
    nexttile;
    K_channel = squeeze(K_List(:, :, row)); % N x 12 (time x state)
    imagesc(K_channel);
    clim_val = max(abs(K_channel(:)), [], 'omitnan');
    if clim_val > 0
        clim([-clim_val, clim_val]);
    end
    colormap(gca, redblue_local());
    colorbar;
    title(sprintf('Gain Heatmap: %s', row_names{row}));
    if row == 3
        xlabel('State index (1-12)');
    end
    ylabel('Horizon step n');
end

% (Col 3, row 1) Roll channel heatmap, shown for completeness even though
% it wasn't part of the "top gains" request above
nexttile;
K_roll = squeeze(K_List(:, :, 4));
imagesc(K_roll);
clim_val = max(abs(K_roll(:)), [], 'omitnan');
if clim_val > 0
    clim([-clim_val, clim_val]);
end
colormap(gca, redblue_local());
colorbar;
title('Gain Heatmap: Roll');
xlabel('State index (1-12)');
ylabel('Horizon step n');

% (Col 3, row 2) Relative gain magnitude by control row, log-scaled so
% the ~10^3 thrust/gimbal scale gap doesn't hide the smaller channels
nexttile;
hold on; grid on;
for row = 1:4
    row_norm = nan(N, 1);
    for n = 1:N
        if ~has_bad_values(n)
            row_norm(n) = norm(squeeze(K_List(n, :, row)));
        end
    end
    plot(1:N, row_norm, 'LineWidth', 1.5, 'Color', row_colors(row, :), ...
         'DisplayName', row_names{row});
end
set(gca, 'YScale', 'log');
title('Per-Channel Gain Magnitude (log scale)');
xlabel('Horizon step n');
ylabel('|| K_n(row,:) ||_2');
legend('Location', 'best', 'FontSize', 7);
hold off;

% (Col 3, row 3) Closed-loop stability margin, evaluated against the
% LOCAL linearization at each step (not a single fixed operating point)
% -- this is the diagnostic that actually matters given how dynamic the
% maneuver is.
nexttile;
plot(1:N, closed_loop_rho, 'LineWidth', 1.5, 'Color', [0 0.5 0]);
hold on;
yline(1, '--r', 'Marginally Stable', 'LabelHorizontalAlignment', 'left');
grid on;
title('Closed-Loop Spectral Radius (local A_d, B_d per step)');
xlabel('Horizon step n');
ylabel('max |eig(A_d(n) - B_d(n) K_n)|');
ylim([0, max(1.1, max(closed_loop_rho, [], 'omitnan')*1.1)]);
hold off;

%% Trajectory Context Figure (ported from TrajectoryGenerator.m)
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

%% Trajectory Position/Velocity & 3D Profile (ported from TrajectoryGenerator.m)
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
linkaxes(findobj(gcf, 'Type', 'axes'), 'x');

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

function [A_fcn, B_fcn] = build_AB_jacobians(constantsTOAD)
% Builds the symbolic A/B Jacobian function handles once. Evaluate the
% returned handles as A_fcn(x_n, x_n-1, u_n-1) / B_fcn(x_n, x_n-1, u_n-1)
% for as many trajectory points as needed -- this is the expensive part,
% so callers should build these handles a single time and reuse them.

% state variables, x
syms q0 q1 q2 q3
syms r1 r2 r3
syms v1 v2 v3
syms omega1 omega2 omega3
syms m_lox m_ipa
% control for u
syms theta phi thrust roll

% state x_n-1
syms q01 q11 q21 q31
syms r11 r21 r31
syms v11 v21 v31
syms omega11 omega21 omega31
syms m_lox1 m_ipa1

q = [q0;q1;q2;q3];
r = [r1;r2;r3];
v = [v1;v2;v3];
m = [m_lox; m_ipa];
omegaB = [omega1; omega2; omega3];

xn = [q;r;v; omegaB; m]; % x_n

% U
un = [theta;phi;thrust;roll];

% X_n-1
qn1 = [q01;q11;q21;q31];
rn1 = [r11;r21;r31];
vn1 = [v11;v21;v31];
mn1 = [m_lox1; m_ipa1];
omegaBn1 = [omega11;omega21;omega31];
xn1 = [qn1;rn1;vn1;omegaBn1;mn1]; % x_n-1

%% Dynamics
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
m_dry = constantsTOAD.m_dry;
m = m_dry + m_lox + m_ipa;

% Propellant Fill height
OxFluidHeight = (m_lox / constantsTOAD.OxMass) * constantsTOAD.OxHeight * 0.9;
FuFluidHeight = (m_ipa / constantsTOAD.FuMass) * constantsTOAD.FuHeight * 0.9;

% Propellant inertias
J_xx = 1/12 * m_lox * (3 * constantsTOAD.OxRadius^2 + OxFluidHeight^2);
J_zz = 1/2 * m_lox * constantsTOAD.OxRadius^2;
J_lox = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

J_xx = 1/12 * m_ipa * (3 * constantsTOAD.FuRadius^2 + FuFluidHeight^2);
J_zz = 1/2 * m_ipa * constantsTOAD.FuRadius^2;
J_ipa = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

% Fluid fills & location of CoM
OxFluidLocation = constantsTOAD.Ox_Z + OxFluidHeight / 2;
FuFluidLocation = constantsTOAD.Fu_Z + FuFluidHeight / 2;
CGz = (constantsTOAD.m_dry * constantsTOAD.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

% Distances to CG
d_dry = constantsTOAD.rTB - CGz;
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Shifted Inertias
J_dry = constantsTOAD.J + m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);
% Bring everything to instantaneous CG (w.r.t Engine attachment frame)
J_tot = J_dry + J_lox + J_ipa;

% Angular dynamics
MB = zetaCross([0; 0; -CGz])*TB + (TB * roll)/thrust;
M = [q(1) -q(2) -q(3) -q(4);
     q(2)  q(1) -q(4)  q(3);
     q(3)  q(4)  q(1) -q(2);
     q(4) -q(3)  q(2)  q(1)];
qdot = 0.5 * M * [0; omegaB];
omegaBdot = J_tot \ (MB - zetaCross(omegaB) * J_tot * omegaB);

% Linear Dynamics
C_BI = quatRot(q);
C_IB = C_BI.';
FI = C_IB * TB + [0; 0;-m*constantsTOAD.g];

rdot = [v1;v2;v3];
vdot = FI/m;

% Mass Dynamics
mdot_lox = -thrust / constantsTOAD.MaxThrust * constantsTOAD.OF / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);
mdot_ipa = -thrust / constantsTOAD.MaxThrust * 1 / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);

%% State vec derivative
xdot = [qdot;rdot;vdot; omegaBdot;mdot_lox;mdot_ipa];
A = jacobian(xdot, xn);
B = jacobian(xdot, un);

% Turn the jacobians into matlab functions for evaluation speed
A_fcn = matlabFunction(A, 'Vars', {xn, xn1, un});
B_fcn = matlabFunction(B, 'Vars', {xn, xn1, un});

end