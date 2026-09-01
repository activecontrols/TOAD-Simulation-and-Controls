%% The following script uses an advanced optimizer and algorithmic differentiation
%% algorithm like CasADi to solve for a complete vehicle trajectory given a set of 
%% constraints and loose waypoints. 
%
% Authors: Andrew Lulo & Pablo Plata
% 7/17/2026

%% Setup and initialization 
CasADiDynamics;
if ~exist("constantsTOAD")
    LoadTOADSim;
end
Filename = "Backflip_v2";
SaveFile = true;
import casadi.*
opti = casadi.Opti();
N = 200;              % Number of control intervals
T_total = opti.variable();
opti.subject_to(15 <= T_total <= 50);  
opti.set_initial(T_total, 35);         
dt = T_total / N;     % Time step

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
%% Problem Scaling
% Characteristic scales chosen so that Xhat, Uhat live near O(1)-O(10) at most.
L_c   = 50;                                   % position scale
V_c   = sqrt(constantsTOAD.g * L_c);          % velocity scale
W_c   = sqrt(constantsTOAD.g / L_c);          % angular rate scale
F_c   = constantsTOAD.MaxThrust;              % thrust scal
G_c   = pi/15;                                % gimbal angle scale
Roll_c = 10;                                  % roll torque scale
Mlox_c = constantsTOAD.OxMass;                % oxidizer mass scale
Mipa_c = constantsTOAD.FuMass;                % fuel mass scale

% State scale vector
Sx = [1; 1; 1; 1; L_c; L_c; L_c; V_c; V_c; V_c; W_c; W_c; W_c; Mlox_c; Mipa_c];

% Input scale vector
Su = [G_c; G_c; F_c; Roll_c];

% Decision Variables
Xhat = opti.variable(15, N+1);
Uhat = opti.variable(4, N);

% Physical-unit aliasesm, use in the script CasADi broadcasts the Nx1 scale vector
% against the NxM variable matrix elementwise.
X = Sx .* Xhat;
U = Su .* Uhat;

%% Integration (RK4), Enforcing dynamics between nodes.
% Build a single-step RK4 integrator as its own CasADi Function
xhat_sym = MX.sym('xhat', 15);
uhat_sym = MX.sym('uhat', 4);
dt_sym   = MX.sym('dt');

% Convert to physical units
x_sym = Sx .* xhat_sym; 
u_sym = Su .* uhat_sym;

k1 = dynamics_fnc(x_sym,                u_sym, params_val);
k2 = dynamics_fnc(x_sym + dt_sym/2*k1,  u_sym, params_val);
k3 = dynamics_fnc(x_sym + dt_sym/2*k2,  u_sym, params_val);
k4 = dynamics_fnc(x_sym + dt_sym*k3,    u_sym, params_val);
x_next_sym = x_sym + dt_sym/6*(k1 + 2*k2 + 2*k3 + k4);

% Back to nondimensional 
xhat_next_sym = x_next_sym ./ Sx;
F_step = Function('F_step', {xhat_sym, uhat_sym, dt_sym}, {xhat_next_sym});

% Map it across all N intervals at once
F_map = F_step.map(N);             
dt_row = repmat(dt, 1, N);          
Xhat_next_all = F_map(Xhat(:, 1:N), Uhat, dt_row);   % single vectorized call, shape 15 x N

% Replace the entire for-loop with one vectorized constraint
opti.subject_to(Xhat(:, 2:end) == Xhat_next_all); 
opti.subject_to(sum(X(1:4, :).^2, 1) == 1.00);
%% Boundaries

MaxThrust_val = constantsTOAD.MaxThrust;
max_gimbal_rate = deg2rad(30);   % deg/s, tune to lin act spec.
max_thrust_rate = 1000;          % N/s
max_roll_rate = 4;

% Initial state (On the pad)
q0 = [1; 0; 0; 0];               % Upright
r0 = [0; 0; 0];
v0 = [0; 0; 0];
w0 = [0; 0; 0];
m_lox0 = constantsTOAD.OxMass;
m_ipa0 = constantsTOAD.FuMass;

% Final state (On the landing zone)
r_f = [0; 0; 0];                 

% Sandbox constraint
opti.subject_to(-30 <= X(5:6, :) <= 30);
opti.subject_to(-1 <= X(7, :) <= 150);

% Initial state
    opti.subject_to(X(:, 1) == [q0; r0; v0; w0; m_lox0; m_ipa0]);

% Final state
    opti.subject_to(X(1:4, end) == q0);
    opti.subject_to(X(5:7, end) == r_f);
    opti.subject_to(sum(X(8:10, end).^2) <= 0.1^2);
    prop_margin_frac = 0.10;   % require >=10% of loaded propellant as reserve
    opti.subject_to(X(14, end) >= prop_margin_frac * m_lox0);
    opti.subject_to(X(15, end) >= prop_margin_frac * m_ipa0);
%% Control Bounds
    thrust_margin  = 0.15;   
    gimbal_margin  = 0.15;   
    
    opti.subject_to((0.25 + thrust_margin) * MaxThrust_val <= U(3,:) <= (1 - thrust_margin) * MaxThrust_val);
    opti.subject_to(-(1 - gimbal_margin) * pi/15 <= U(1,:) <= (1 - gimbal_margin) * pi/15);
    opti.subject_to(-(1 - gimbal_margin) * pi/15 <= U(2,:) <= (1 - gimbal_margin) * pi/15);
    opti.subject_to(-(1 - thrust_margin) * max_roll_rate <= U(4,:) <= (1 - thrust_margin) * max_roll_rate);

%% Rate constraints (physical rate limits — keep on U, not Uhat)
dU_phys = U(:, 2:end) - U(:, 1:end-1);
opti.subject_to(-max_gimbal_rate*dt <= dU_phys(1,:) <= max_gimbal_rate*dt);
opti.subject_to(-max_gimbal_rate*dt <= dU_phys(2,:) <= max_gimbal_rate*dt);
opti.subject_to(-max_thrust_rate*dt <= dU_phys(3,:) <= max_thrust_rate*dt);
opti.subject_to(-max_roll_rate*dt <= dU_phys(4,:) <= max_roll_rate*dt);

%% Trajectory 
N_ascent   = round(0.05*N);
N_flip     = round(0.5*N);
N_approach = round(0.95*N);

% Ascent
    opti.subject_to(X(1:4, N_ascent) == q0)
    pos_desc = X(5:6, 1:N_ascent);
    vel_desc = X(8:9, 1:N_ascent);
    opti.subject_to(pos_desc(1,:).^2 + pos_desc(2,:).^2 <= 1^2);
    opti.subject_to(vel_desc(1,:).^2 + vel_desc(2,:).^2 <= 1^2);
    
% Flip Maneuver
    theta_tol = deg2rad(25); % allow 15 degrees of rotational slack
    att_tol = cos(theta_tol/2); 
    
    % Flip target attitude
    q_inverted = [0; 0; 1; 0];
    q_flip = X(1:4, N_flip);
    
    % Attitude Slack (Quaternion Inner Product)
    opti.subject_to( (q_inverted' * q_flip)^2 >= att_tol^2 );
    opti.subject_to(X(7, N_flip) >= 40);
    
% Descent 
    opti.subject_to(X(1:4, N_approach) == q0)
    pos_desc = X(5:7, N_approach:end);
    vel_desc = X(8:10, N_approach:end);
    opti.subject_to( (pos_desc(1,:) - r_f(1)).^2 + (pos_desc(2,:) - r_f(2)).^2 <= 1^2 );
    opti.subject_to( vel_desc(1,:).^2 + vel_desc(2,:).^2 <= 1^2 );
    opti.subject_to(-2 <= vel_desc(3,:) <= 2)
        
%% Initial Guess
% Linearly interpolate positions from start to end
r_guess = [linspace(r0(1), r_f(1), N+1);
           linspace(r0(2), r_f(2), N+1);
           linspace(r0(3) + 10, r_f(3), N+1)];

% Set constant upright quaternion guess (avoids zero-norm singularity)
q_guess = repmat([1; 0; 0; 0], 1, N+1);
opti.set_initial(Xhat(1:4, :), q_guess);             
opti.set_initial(Xhat(5:7, :), r_guess / L_c);         
opti.set_initial(Xhat(14,:), linspace(1, 0.05, N+1)); 
opti.set_initial(Xhat(15,:), linspace(1, 0.05, N+1));  
opti.set_initial(Xhat(8:10,:), zeros(3, N+1));
opti.set_initial(Xhat(11:13,:), zeros(3, N+1));
opti.set_initial(Uhat(3, :), repmat(constantsTOAD.m_wet * constantsTOAD.g / F_c, 1, N));  

%% Cost Function — Controllability & Survivability Weighted. Backflip specific.
% Time in critical region
% R33 = cos(tilt angle from vertical): +1 upright, -1 fully inverted.
R33 = X(1,:).^2 - X(2,:).^2 - X(3,:).^2 + X(4,:).^2;

k_risk      = 6;    % logistic steepness; higher = sharper on/off transition
tilt_thresh = 0;    % R33 = 0 <-> 90 deg tilt. Shift to redefine "critical".
risk = 1 ./ (1 + exp(k_risk .* (R33 - tilt_thresh)));   
dt_frac = dt / T_total;                                  
J_critical = N * sum(risk(1:end-1) .* dt_frac);           

% Penalize body rates while inside the critical band
J_rate_flip = sum(risk(1:end-1) .* sum(Xhat(11:13, 1:end-1).^2, 1));

% Margin 
gimbal_bound   = (1 - gimbal_margin);                     
J_marginGimbal = sum(sum((Uhat(1:2, :) / gimbal_bound).^2));
Tmin_hat  = (0.25 + thrust_margin) * MaxThrust_val / F_c;
Tmax_hat  = (1 - thrust_margin)    * MaxThrust_val / F_c;
Tmid_hat  = (Tmin_hat + Tmax_hat) / 2;
Thalf_hat = (Tmax_hat - Tmin_hat) / 2;
J_marginThrust = sum(((Uhat(3, :) - Tmid_hat) / Thalf_hat).^2);

roll_bound   = (1 - thrust_margin) * max_roll_rate / Roll_c;
J_marginRoll = sum((Uhat(4, :) / roll_bound).^2);

%% Main costs
w_crit         = 1;
w_rate_flip    = 2;
w_marginGimbal = 1e-2;
w_marginThrust = 1e-4;
w_marginRoll   = 1e-2;

opti.minimize( ...
      w_crit         * J_critical      ...
    + w_rate_flip    * J_rate_flip     ...
    + w_marginGimbal * J_marginGimbal  ...
    + w_marginThrust * J_marginThrust  ...
    + w_marginRoll   * J_marginRoll);

fprintf('Starting Coarse Solve!\n');
p_opts_coarse = struct('expand', true);
s_opts_coarse = struct('max_iter', 1500, 'tol', 5e-3, 'constr_viol_tol', 1e-2);
opti.solver('ipopt', p_opts_coarse, s_opts_coarse);

try
    sol_coarse = opti.solve();
    fprintf('Coarse solve successful!\n');
    
    % Extract values from successful solve
    x_warm   = sol_coarse.value(opti.x);
    lam_warm = sol_coarse.value(opti.lam_g);
    
catch
    fprintf('Coarse solve failed to converge. Extracting debug values...\n');
    
    % Extract values from the last solver iteration using opti.debug
    x_warm   = opti.debug.value(opti.x);
    lam_warm = opti.debug.value(opti.lam_g);
end

% Set the initial guess for the next (fine) solve
opti.set_initial(opti.x, x_warm);
opti.set_initial(opti.lam_g, lam_warm);

% Stage B: full-precision solve, warm-started
p_opts = struct('expand', true);
s_opts = struct('max_iter', 5000, 'tol', 1e-4, 'constr_viol_tol', 1e-3, ...
                 'warm_start_init_point', 'yes', 'mu_strategy', 'adaptive');
opti.solver('ipopt', p_opts, s_opts);
fprintf('Starting Full Solve!\n');

%% Solve
try
    sol = opti.solve();
    disp('Optimal trajectory found!');
    
    % Extract results for success
    X_res = Sx .* sol.value(Xhat);
    U_res = Su .* sol.value(Uhat);
    T_res = sol.value(T_total);
    
catch e
    disp('Solver failed. Retrieving debug values...');
    
    % Extract results using the debug interface directly
    X_res = Sx .* opti.debug.value(Xhat);
    U_res = Su .* opti.debug.value(Uhat);
    T_res = opti.debug.value(T_total);
end

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

%% Save Trajectory
if SaveFile
    % Define the target directory
    save_dir = "Guidance\Trajectories\";
    
    % Extract additional variables from the results matrices
    ang_rate = X_res(11:13, :);
    m_lox  = X_res(14, :);
    m_fuel = X_res(15, :);
    roll_cmd = U_res(4, :);
    
    % Pad control arrays (length N) to match the length of state arrays (length N+1)
    % by holding the final control command for the last time step.
    theta_out  = [theta_cmd, theta_cmd(end)]';
    phi_out    = [phi_cmd, phi_cmd(end)]';
    thrust_out = [thrust, thrust(end)]';
    roll_out   = [roll_cmd, roll_cmd(end)]';
    
    % Compile data into a table matching the newly specified order:
    % quat, pos, vel, ang_rate, masslox, massfuel, gimbal theta, phi, thrust, roll
    trajectory_data = table(t_state',...
        quat(1,:)', quat(2,:)', quat(3,:)', quat(4,:)', ...
        pos(1,:)', pos(2,:)', pos(3,:)', ...
        vel(1,:)', vel(2,:)', vel(3,:)', ...
        ang_rate(1,:)', ang_rate(2,:)', ang_rate(3,:)', ...
        m_lox', m_fuel', ...
        theta_out, phi_out, thrust_out, roll_out, ...
        'VariableNames', {'Time', 'QuatW', 'QuatX', 'QuatY', 'QuatZ', ...
                          'PosX', 'PosY', 'PosZ', ...
                          'VelX', 'VelY', 'VelZ', ...
                          'AngRateX', 'AngRateY', 'AngRateZ', ...
                          'MassLox', 'MassFuel', ...
                          'GimbalTheta', 'GimbalPhi', 'Thrust', 'RollCmd'});
                          
    % Define the full file path and write the table to a .csv file
    full_path = save_dir + Filename + ".csv";
    writetable(trajectory_data, full_path);
    fprintf('Trajectory successfully saved to: %s\n', full_path);
end