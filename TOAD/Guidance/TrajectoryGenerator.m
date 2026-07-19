%% The following script uses an advanced optimizer and algorithmic differentiation
%% algorithm like CasADi to solve for a complete vehicle trajectory given a set of 
%% constraints and loose waypoints. 
%
% Authors: Andrew Lulo & Pablo Plata
% 7/17/2026

%% Setup and initialization 
close all
import casadi.*
opti = casadi.Opti();
N = 400;              % Number of control intervals
T_total = 40; % opti.variable(); % Total flight time (can be a decision variable or fixed)
% opti.subject_to(T_total > 0);
% opti.set_initial(T_total, 10); % Initial guess, using as decision 
dt = T_total / N;     % Time step

% Decision Variables
% State vector x: [q(4); r(3); v(3); omegaB(3); m_lox(1); m_ipa(1)]
X = opti.variable(15, N+1);
% Control vector u: [theta(1); phi(1); thrust(1); roll(1)]
U = opti.variable(4, N);

% Force dispersion parameters to zero
MaxMdot_d_val = 0;
J_d_vec = zeros(9, 1);
TB_d_val = zeros(3, 1);

% System parameters
params_val = [
    constantsTOAD.m_dry;      % m_dry
    constantsTOAD.g;          % g
    constantsTOAD.rTB;        % rTB
    constantsTOAD.Ox_Z;       % Ox_Z
    constantsTOAD.OxMass;     % OxMassI
    constantsTOAD.OxHeight;   % OxHeight
    constantsTOAD.Fu_Z;       % Fu_Z
    constantsTOAD.FuMass;     % FuMassI
    constantsTOAD.FuHeight;   % FuHeight
    constantsTOAD.J(:);       % Inertia
    constantsTOAD.OxRadius;   % OxRadius
    constantsTOAD.FuRadius;   % FuRadius
    constantsTOAD.MaxThrust;  % MaxThrust
    constantsTOAD.OF;         % OF
    constantsTOAD.MaxMdot;    % MaxMdot
    MaxMdot_d_val;            % MaxMdot_d 
    J_d_vec;                  % J_d(:)
    TB_d_val                  % TB_d 
];
%% Integration Loop (RK4), Enforcing dynamics between nodes.
for k = 1:N
    x_k = X(:, k);
    u_k = U(:, k);
    
    % RK4 substeps using your dynamics_fnc
    k1 = dynamics_fnc(x_k,             u_k, params_val);
    k2 = dynamics_fnc(x_k + dt/2 * k1, u_k, params_val);
    k3 = dynamics_fnc(x_k + dt/2 * k2, u_k, params_val);
    k4 = dynamics_fnc(x_k + dt * k3,   u_k, params_val);
    
    x_next = x_k + dt/6 * (k1 + 2*k2 + 2*k3 + k4);
    
    % "Close the gap" constraint (Multiple Shooting)
    opti.subject_to(X(:, k+1) == x_next);
end

%% Path Constraints
for k = 1:N+1
    % Enforce unit for quaternion
    q_k = X(1:4, k);
    opti.subject_to(sum(q_k.^2) == 1);
end

%% Control Bounds 
MaxThrust_val = constantsTOAD.MaxThrust;
% opti.subject_to(-pi/24 <= U(1,:) <= pi/24); % Gimbal theta limits
% opti.subject_to(-pi/24 <= U(2,:) <= pi/24); % Gimbal phi limits 
opti.subject_to(-pi/15 <= U(1,:) <= pi/15); % Gimbal theta limits
opti.subject_to(-pi/15 <= U(2,:) <= pi/15); % Gimbal phi limits 
opti.subject_to(0.5 * MaxThrust_val <= U(3,:) <= MaxThrust_val); % Throttle
opti.subject_to(-10 <= U(4,:) <= 10);     % Roll torque limits

%% Rate constraints
max_gimbal_rate = deg2rad(20);   % deg/s, tune to lin act spec.
max_thrust_rate = 5000;          % N/s

for k = 1:N-1
    opti.subject_to(-max_gimbal_rate*dt <= U(1,k+1)-U(1,k) <= max_gimbal_rate*dt);
    opti.subject_to(-max_gimbal_rate*dt <= U(2,k+1)-U(2,k) <= max_gimbal_rate*dt);
    opti.subject_to(-max_thrust_rate*dt <= U(3,k+1)-U(3,k) <= max_thrust_rate*dt);
end

%% Trajectory Design (Boundaries)

    N_ascent = round(0.2*N);
    N_flip = round(0.5 * N);
    N_approach = round(0.8*N);
    
% Initial state (On the pad)
    q0 = [1; 0; 0; 0]; % Upright
    r0 = [0; 0; 0];
    v0 = [0; 0; 0];
    w0 = [0; 0; 0];
    m_lox0 = constantsTOAD.OxMass;
    m_ipa0 = constantsTOAD.FuMass;
    opti.subject_to(X(:, 1) == [q0; r0; v0; w0; m_lox0; m_ipa0]);

% Final state (On the landing zone)
    r_f = [5; 5; 0]; % Example downrange landing pad
    opti.subject_to(X(1:4, end) == q0); % Upright upon landing
    opti.subject_to(X(5:7, end) == r_f); 
    opti.subject_to(X(8:10, end) == [0;0;0]); % Zero velocity

%% Trajectory
    
% Ascent 
    for k = 1:N_ascent
        opti.subject_to(-1 <= X(5:6, k) <= 1);
        opti.subject_to(X(7, k) >= -1);
    end
    opti.subject_to(X(10, N_ascent) >= 5);
    
% Flip
    q_inverted = [0; 0; 1; 0];
    q_flip = X(1:4, N_flip);
    
    % Squared quaternion dot product >= 0.98 ensures the attitude is very close 
    % to inverted, but allows the solver a slight tolerance.
    dot_prod = q_inverted' * q_flip;
    opti.subject_to(dot_prod^2 >= 0.98);
    opti.subject_to(X(7, N_flip) >= 30);

    for k = 1:N+1
        if k ~= N_flip
            opti.subject_to(-50 <= X(5:6, k) <= 50);
            opti.subject_to(-1 <= X(7, k) <= 200);
        end
    end

% Descent (glideslope constrained)
    % glideslope_angle = deg2rad(15);
    % for k = N_approach:N+1
    %     horiz_dist_sq = (X(5,k) - r_f(1))^2 + (X(6,k) - r_f(2))^2;
    %     opti.subject_to(horiz_dist_sq <= (tan(glideslope_angle) * X(7,k))^2);
    %     opti.subject_to(X(7, k) >= 0);
    % end

%% Initial Guess
% Linearly interpolate positions from start to end
r_guess = [linspace(r0(1), r_f(1), N+1);
           linspace(r0(2), r_f(2), N+1);
           linspace(r0(3) + 10, r_f(3), N+1)];

% Set constant upright quaternion guess (avoids zero-norm singularity)
q_guess = repmat([1; 0; 0; 0], 1, N+1);

opti.set_initial(X(1:4, :), q_guess);
opti.set_initial(X(5:7, :), r_guess);
opti.set_initial(U(3, :), repmat(constantsTOAD.m_wet * constantsTOAD.g, 1, N));

%% Cost Function (needs lots of improvement)
% Add weightings for control usage
w_rate = [5e1; 5e1; 1e-4; 1e-1];
w_mag = [1e-2; 1e-2; 0; 1e-3];
w_omega = 1e-2;
w_pos = 1e-2;
w_vel = 1e-1;

dU = U(:, 2:end) - U(:, 1:end-1);
pos_cost = w_vel * sum(sum(X(5:7, :).^2));
vel_cost = w_vel * sum(sum(X(8:10, :).^2));
rate_cost = sum(sum(w_rate .* dU.^2));
reg_cost = sum(sum(w_mag .* U.^2));
omega_cost = w_omega * sum(sum(X(11:13, :).^2));

opti.minimize(-(X(14,end) + X(15,end)) + 5*T_total + reg_cost + rate_cost + omega_cost + vel_cost);

%% Solver Configuration
p_opts = struct('expand', true);
s_opts = struct('max_iter', 3000, 'tol', 1e-4, 'constr_viol_tol', 1e-4);
opti.solver('ipopt', p_opts, s_opts);

%% Solve
try
    sol = opti.solve();
    disp('Optimal trajectory found!');
catch
    disp('Solver failed. Retrieving debug values...');
    sol = opti.debug;
end

% Extract results for plotting
X_res = sol.value(X);
U_res = sol.value(U);
T_res = sol.value(T_total);

%% --- Post-Processing & Visualization ---
% Extract state and control arrays for plotting
t_state = linspace(0, T_res, N+1);
t_ctrl  = linspace(0, T_res - (T_res/N), N);

pos = X_res(5:7, :);
vel = X_res(8:10, :);
quat = X_res(1:4, :); % [w; x; y; z]

theta_cmd = U_res(1, :);
phi_cmd   = U_res(2, :);
thrust    = U_res(3, :);

% Convert quaternions to Euler angles (ZYX sequence)[cite: 7]
% quat2eul expects an Nx4 matrix in [w x y z] format
eul_angles = quat2eul(quat', 'ZYX'); 

% Plot 1: 3D Mission Trajectory 
figure('Name', 'Trajectory: 3D Mission Profile', 'WindowStyle', 'docked');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact'); %[cite: 6, 7]

% Main 3D View
axMain = nexttile(tl, 1, [3 3]); 
hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
xlabel(axMain, 'X [m]'); ylabel(axMain, 'Y [m]'); zlabel(axMain, 'Z (Alt) [m]');
title(axMain, '3D Mission Trajectory');

% Plot trajectory with a color gradient mapped to time using patch[cite: 7]
patch(axMain, [pos(1,:), NaN], [pos(2,:), NaN], [pos(3,:), NaN], [t_state, NaN], ...
      'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5);
cb = colorbar(axMain); cb.Label.String = 'Time [s]'; colormap(axMain, 'turbo'); %[cite: 7]

% Plot start and end boundary constraints
plot3(axMain, r0(1), r0(2), r0(3), 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
plot3(axMain, r_f(1), r_f(2), r_f(3), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);

% Orthographic Projections[cite: 6, 7]
axTop = nexttile(tl, 4);  
hold(axTop, 'on'); grid(axTop, 'on'); axis(axTop, 'equal');
plot(axTop, pos(1,:), pos(2,:), 'b', 'LineWidth', 1.5);
xlabel(axTop, 'X [m]'); ylabel(axTop, 'Y [m]'); title(axTop, 'Top View (X-Y)');

axSide = nexttile(tl, 8); 
hold(axSide, 'on'); grid(axSide, 'on'); axis(axSide, 'equal');
plot(axSide, pos(1,:), pos(3,:), 'b', 'LineWidth', 1.5);
xlabel(axSide, 'X [m]'); ylabel(axSide, 'Z (Alt) [m]'); title(axSide, 'Side View (X-Z)');

axFront = nexttile(tl, 12); 
hold(axFront, 'on'); grid(axFront, 'on'); axis(axFront, 'equal');
plot(axFront, pos(2,:), pos(3,:), 'b', 'LineWidth', 1.5);
xlabel(axFront, 'Y [m]'); ylabel(axFront, 'Z (Alt) [m]'); title(axFront, 'Front View (Y-Z)');

% Plot 2: Controller & Attitude Performance 
figure('Name', 'Trajectory: Control & Attitude', 'WindowStyle', 'docked');

% Attitude Tracking (Quaternions & Tilt Angle)
subplot(3,1,1); hold on; grid on;

% Left Axis: Raw Quaternions (Smooth and singularity-free)
yyaxis left;
plot(t_state, quat(1, :), 'k', 'LineWidth', 1.5, 'DisplayName', 'q_w (Scalar)');
plot(t_state, quat(2, :), 'r', 'LineWidth', 1.5, 'DisplayName', 'q_x');
plot(t_state, quat(3, :), 'g', 'LineWidth', 1.5, 'DisplayName', 'q_y');
plot(t_state, quat(4, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'q_z');
ylabel('Quaternion Component');
ylim([-1.1 1.1]);

% Right Axis: Absolute Tilt Angle from Vertical
% R33 element of the Direction Cosine Matrix gives the Z-Z projection
R33 = quat(1,:).^2 - quat(2,:).^2 - quat(3,:).^2 + quat(4,:).^2;
tilt_angle = acosd(max(min(R33, 1), -1)); % Clamp between [-1, 1] for numerical safety

yyaxis right; 
plot(t_state, tilt_angle, '--m', 'LineWidth', 2, 'DisplayName', 'Tilt Angle (deg)');
ylabel('Tilt Angle [deg]');
ylim([0 185]); 

% Add legend and format
legend('Location', 'eastoutside');
title('Vehicle Attitude (Quaternions & Absolute Tilt)');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = 'm';

% Gimbal Deflection Commands
subplot(3,1,2); hold on; grid on;
plot(t_ctrl, rad2deg(theta_cmd), 'Color', '#7E2F8E', 'LineWidth', 1.5);
plot(t_ctrl, rad2deg(phi_cmd), 'Color', '#77AC30', 'LineWidth', 1.5);
ylabel('Gimbal [deg]'); legend('\theta', '\phi', 'Location', 'best');
title('Gimbal Deflection Commands');

% Thrust Command
subplot(3,1,3); hold on; grid on;
area(t_ctrl, thrust, 'FaceColor', '#EDB120', 'FaceAlpha', 0.4, 'EdgeColor', '#A27712', 'LineWidth', 1.5);
ylabel('Thrust [N]'); xlabel('Time [s]');
title('Thrust Command');

% Plot 3: Mission Kinematics (Position & Velocity)
figure('Name', 'Trajectory: Mission Kinematics', 'WindowStyle', 'docked');
tl_kin = tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact'); %[cite: 7]
title(tl_kin, 'Position and Velocity Profiles', 'FontWeight', 'bold');
xlabel(tl_kin, 'Time [s]', 'FontWeight', 'bold'); 

labels_pos = {'X Pos [m]', 'Y Pos [m]', 'Z Pos [m]'};
labels_vel = {'X Vel [m/s]', 'Y Vel [m/s]', 'Z Vel [m/s]'};
colors = {'#0072BD', '#D95319', '#EDB120'};

for ax_idx = 1:3
    % Position Tracking Subplot[cite: 7]
    nexttile(tl_kin, (ax_idx-1)*2 + 1); hold on; grid on;
    plot(t_state, pos(ax_idx, :), 'Color', colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_pos{ax_idx}, 'FontWeight', 'bold');
    
    % Velocity Tracking Subplot[cite: 7]
    nexttile(tl_kin, (ax_idx-1)*2 + 2); hold on; grid on;
    plot(t_state, vel(ax_idx, :), 'Color', colors{ax_idx}, 'LineWidth', 2);
    ylabel(labels_vel{ax_idx}, 'FontWeight', 'bold');
end

% Link the X-axes for zooming across all 6 kinematic subplots[cite: 7]
linkaxes(findobj(gcf, 'Type', 'axes'), 'x');