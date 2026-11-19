function out = Plot6DoF(varargin)
% PLOT6DOF Visualizes 6-DoF trajectory, 2D traces, and real-time 3D flight animation for TOAD.
%
% Visual Features in Real-Time 3D Animation:
%   - Dynamic XY ground projection footprint & vertical altitude drop line
%   - Orange thrust vector line scaled proportionally with engine thrust command
%   - 3D body frame coordinate triad (+X Red, +Y Green, +Z Blue) centered at vehicle middle
%   - Real-time HUD banner displaying flight time, altitude, thrust %, speed, and tracking error
%
% Syntax:
%   Plot6DoF()
%   Plot6DoF(simOut)
%   Plot6DoF(..., 'TimeScale', 1.0, 'VehicleScale', 2.5, 'ShowWaypoints', false)
%   out = Plot6DoF(...)
%
% Inputs:
%   simOut        - (Optional) Simulink.SimulationOutput object from TOAD_Simulation.
%                   If omitted, LoadTOADSim is executed and TOAD_Simulation is run.
%
% Parameter Name-Value Pairs:
%   'TimeScale'     - Playback speed multiplier for animation (default: 1.0 = real-time).
%   'VehicleScale'  - Visual scale multiplier for vehicle 3D model (default: 1.0).
%   'Animate'       - Logical flag to run AnimateFlight (default: true).
%   'ShowWaypoints' - Logical flag to display waypoint/checkpoint markers (default: false).
%   'SampleTime'    - Logging sample time configured on "To Workspace" blocks (default: 0.02s = 50Hz).
%   'StopTime'      - Simulation stop time if running simulation (default: model setting).
%   'Vehicle'       - Vehicle type if running simulation: 1 for TOAD, 0 for ASTRA (default: from LoadTOADSim).
%
% Outputs:
%   out - Struct containing processed flight data (t, r, v, q, u, maxThrust, r_target, simOut).

%% 1. Parse Input Arguments
p = inputParser;
p.KeepUnmatched = true;
addOptional(p, 'simOut', [], @(x) isempty(x) || isa(x, 'Simulink.SimulationOutput') || isstruct(x));
addParameter(p, 'TimeScale', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'VehicleScale', 1, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'Animate', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'ShowWaypoints', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'SampleTime', 0.02, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'StopTime', '', @(x) ischar(x) || isstring(x) || isnumeric(x));
addParameter(p, 'Vehicle', [], @(x) isempty(x) || isnumeric(x));

parse(p, varargin{:});
simOut        = p.Results.simOut;
timeScale     = p.Results.TimeScale;
vehicleScale  = p.Results.VehicleScale;
doAnimate     = logical(p.Results.Animate);
showWaypoints = logical(p.Results.ShowWaypoints);
sampleTime    = p.Results.SampleTime;
stopTime      = p.Results.StopTime;
vehicleSel    = p.Results.Vehicle;

%% 2. Run Simulation if simOut Not Provided
modelName = 'TOAD_Simulation';
if isempty(simOut)
    fprintf('Running simulation setup for %s...\n', modelName);
    
    % If vehicle specified, set in base workspace
    if ~isempty(vehicleSel)
        assignin('base', 'Vehicle', vehicleSel);
    end
    
    % Run LoadTOADSim to populate workspace
    evalin('base', 'LoadTOADSim;');
    
    % Load Simulink model
    if ~bdIsLoaded(modelName)
        load_system(modelName);
    end
    
    % Configure "To Workspace" blocks for high-resolution logging
    ConfigureToWorkspaceBlocks(modelName, sampleTime);
    
    % Configure stop time if specified
    if ~isempty(stopTime)
        if isnumeric(stopTime), stopTime = num2str(stopTime); end
        set_param(modelName, 'StopTime', stopTime);
    end
    
    fprintf('Simulating %s...\n', modelName);
    simOut = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
    fprintf('Simulation completed.\n');
end

%% 3. Extract States and Target Trajectory
% Verify state_log presence
if ~isprop(simOut, 'state_log') && ~isfield(simOut, 'state_log')
    error('Plot6DoF:MissingStateLog', ...
        'state_log not found in simOut. Ensure TOAD_Simulation logs state_log.');
end

stateTimeSeries = simOut.state_log;
t = stateTimeSeries.Time;
x_data = formatData(stateTimeSeries.Data, length(t));

% Map states: 1:4 quaternion [qw qx qy qz], 5:7 position [x y z], 8:10 velocity
q = x_data(:, 1:4);
r = x_data(:, 5:7);
if size(x_data, 2) >= 10
    v = x_data(:, 8:10);
else
    v = gradient(r, mean(diff(t)));
end

% Normalize quaternions
q_norms = sqrt(sum(q.^2, 2));
q_norms(q_norms < 1e-8) = 1.0;
q = q ./ q_norms;

% Extract Target Trajectory from target_pos_log
r_target = [];
if isprop(simOut, 'target_pos_log') || isfield(simOut, 'target_pos_log')
    targetTS = simOut.target_pos_log;
    if isa(targetTS, 'timeseries')
        t_targ = targetTS.Time;
        targ_data = formatData(targetTS.Data, length(t_targ));
        if length(t_targ) == length(t) && max(abs(t_targ - t)) < 1e-6
            r_target = targ_data(:, 5:7);
        else
            r_target = interp1(t_targ, targ_data(:, 5:7), t, 'linear', 'extrap');
        end
    end
end

% Fallback target if target_pos_log was empty or missing
if isempty(r_target)
    try
        constants6DoF = evalin('base', 'constants6DoF');
        if isfield(constants6DoF, 'Traj') && isfield(constants6DoF.Traj, 'States')
            trajStates = constants6DoF.Traj.States;
            trajTime   = constants6DoF.Traj.Time;
            r_target = interp1(trajTime, trajStates(:, 5:7), t, 'linear', 'extrap');
        end
    catch
        % No reference trajectory available; use initial position
        r_target = repmat(r(1,:), length(t), 1);
    end
end

% Extract Waypoints only if explicitly requested
Waypoints = [];
if showWaypoints
    try
        Waypoints = evalin('base', 'Waypoints');
    catch
        try
            Waypoints = TrajectoryBuilder();
        catch
            Waypoints = [];
        end
    end
end

% Extract Control Inputs for Thrust & Gimbal Visualization
u_data = [];
if isprop(simOut, 'input_log') || isfield(simOut, 'input_log')
    inTS = simOut.input_log;
    if isa(inTS, 'timeseries')
        t_in = inTS.Time;
        raw_in = formatData(inTS.Data, length(t_in));
        if length(t_in) == length(t) && max(abs(t_in - t)) < 1e-6
            u_data = raw_in;
        else
            u_data = interp1(t_in, raw_in, t, 'linear', 'extrap');
        end
    end
elseif isprop(simOut, 'inputCMD') || isfield(simOut, 'inputCMD')
    inTS = simOut.inputCMD;
    if isa(inTS, 'timeseries')
        t_in = inTS.Time;
        raw_in = formatData(inTS.Data, length(t_in));
        if length(t_in) == length(t) && max(abs(t_in - t)) < 1e-6
            u_data = raw_in;
        else
            u_data = interp1(t_in, raw_in, t, 'linear', 'extrap');
        end
    end
end

% Fallback to constants6DoF.Traj.Inputs if simOut did not contain input logs
if isempty(u_data)
    try
        constants6DoF = evalin('base', 'constants6DoF');
        if isfield(constants6DoF, 'Traj') && isfield(constants6DoF.Traj, 'Inputs')
            trajIn   = constants6DoF.Traj.Inputs;
            trajTime = constants6DoF.Traj.Time;
            u_data = interp1(trajTime, trajIn, t, 'linear', 'extrap');
        end
    catch
        u_data = zeros(length(t), 4);
    end
end

% Extract MaxThrust for proportional thrust line scaling
maxThrust = 0;
try
    constants6DoF = evalin('base', 'constants6DoF');
    if isfield(constants6DoF, 'MaxThrust')
        maxThrust = constants6DoF.MaxThrust;
    end
catch
end
if maxThrust <= 0
    try
        constantsTOAD = evalin('base', 'constantsTOAD');
        if isfield(constantsTOAD, 'MaxThrust')
            maxThrust = constantsTOAD.MaxThrust;
        end
    catch
    end
end
if maxThrust <= 0
    if ~isempty(u_data) && size(u_data, 2) >= 3
        maxThrust = max(u_data(:, 3));
    end
    if maxThrust <= 0
        maxThrust = 16.25; % Default fallback (ASTRA nominal max thrust)
    end
end

%% 4. Figure 1: 3D Trajectory & Target Following
figure('Name', 'TOAD: 3D Flight Trajectory', 'Color', 'w');
colormap('turbo');

% Surface color-gradient for flight path
patch([r(:,1); NaN], [r(:,2); NaN], [r(:,3); NaN], [t; NaN], ...
     'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5, ...
     'DisplayName', 'Vehicle Trajectory');
hold on; grid on; axis equal; view(3);

cb = colorbar;
cb.Label.String = 'Flight Time (s)';
cb.Label.FontSize = 10;

% Plot Target Trajectory
if ~isempty(r_target)
    plot3(r_target(:,1), r_target(:,2), r_target(:,3), 'r--', ...
        'LineWidth', 1.8, 'DisplayName', 'Target Trajectory');
end

% Plot Waypoints only if showWaypoints is true
if showWaypoints && ~isempty(Waypoints)
    numWP = numel(Waypoints);
    wp_coords = zeros(3, numWP);
    for idx = 1:numWP
        wp_coords(:, idx) = Waypoints(idx).Position;
    end
    plot3(wp_coords(1,:), wp_coords(2,:), wp_coords(3,:), 'ks', ...
        'MarkerSize', 8, 'MarkerFaceColor', '#EDB120', 'LineWidth', 1.2, ...
        'DisplayName', 'Waypoints');
end

% Mark Start and End Positions
plot3(r(1,1), r(1,2), r(1,3), 'go', 'MarkerSize', 9, 'MarkerFaceColor', 'g', ...
    'DisplayName', 'Liftoff');
plot3(r(end,1), r(end,2), r(end,3), 'mo', 'MarkerSize', 9, 'MarkerFaceColor', 'm', ...
    'DisplayName', 'Touchdown / End');

xlabel('North [m]', 'FontWeight', 'bold');
ylabel('West [m]', 'FontWeight', 'bold');
zlabel('Altitude [m]', 'FontWeight', 'bold');
title('TOAD 6-DoF 3D Trajectory with Target Tracking', 'FontWeight', 'bold');
legend('show', 'Location', 'best');
hold off;

%% 5. Figure 2: 2D Trajectory Traces (Synchronized Multi-View)
figure('Name', 'TOAD: 2D Trajectory Projections', 'Color', 'w');
colormap('turbo');
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% XY Trace (Top View: West vs North)
nexttile; hold on; grid on; axis equal;
patch([r(:,2); NaN], [r(:,1); NaN], [t; NaN], ...
     'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2);
if ~isempty(r_target)
    plot(r_target(:,2), r_target(:,1), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Target');
end
if showWaypoints && ~isempty(Waypoints)
    plot(wp_coords(2,:), wp_coords(1,:), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', '#EDB120');
end
xlabel('West [m]', 'FontWeight', 'bold');
ylabel('North [m]', 'FontWeight', 'bold');
title('Top View (XY)', 'FontWeight', 'bold');

% XZ Trace (Side View: North vs Altitude)
nexttile; hold on; grid on; axis equal;
patch([r(:,1); NaN], [r(:,3); NaN], [t; NaN], ...
     'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2);
if ~isempty(r_target)
    plot(r_target(:,1), r_target(:,3), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Target');
end
if showWaypoints && ~isempty(Waypoints)
    plot(wp_coords(1,:), wp_coords(3,:), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', '#EDB120');
end
xlabel('North [m]', 'FontWeight', 'bold');
ylabel('Altitude [m]', 'FontWeight', 'bold');
title('Side View (XZ)', 'FontWeight', 'bold');

% YZ Trace (Front View: West vs Altitude)
nexttile; hold on; grid on; axis equal;
patch([r(:,2); NaN], [r(:,3); NaN], [t; NaN], ...
     'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2);
if ~isempty(r_target)
    plot(r_target(:,2), r_target(:,3), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Target');
end
if showWaypoints && ~isempty(Waypoints)
    plot(wp_coords(2,:), wp_coords(3,:), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', '#EDB120');
end
xlabel('West [m]', 'FontWeight', 'bold');
ylabel('Altitude [m]', 'FontWeight', 'bold');
title('Front View (YZ)', 'FontWeight', 'bold');

% Shared Colorbar
cb2 = colorbar('Position', [0.93 0.15 0.015 0.7]);
cb2.Label.String = 'Time [s]';
cb2.Label.FontSize = 10;
sgtitle('TOAD 2D Trajectory Traces (North-West-Up)', 'FontWeight', 'bold');

%% 6. Spectral Analysis / Spectrogram
if isprop(simOut, 'meas_log') || isfield(simOut, 'meas_log')
    measurementLog = simOut.meas_log;
    if isfield(measurementLog, 'Data') || isprop(measurementLog, 'Data')
        figure('Name', 'Sensor Spectral Analysis', 'Color', 'w');
        dt_mean = mean(diff(t));
        fs = 1 / dt_mean;
        windowSize = min(256, 2^nextpow2(floor(length(t)/4)));
        if windowSize >= 16
            overlap = floor(windowSize * 0.9);
            nfft = max(512, 2 * windowSize);
            measData = measurementLog.Data;
            if size(measData, 1) > 3
                sig = measData(4, :);
            else
                sig = measData(1, :);
            end
            spectrogram(sig, kaiser(windowSize, 5), overlap, nfft, fs, 'yaxis');
            title('IMU Measurement Spectrogram', 'FontWeight', 'bold');
        end
    end
end

%% 7. Real-Time Flight Animation
if doAnimate
    AnimateFlight(t, r, q, r_target, ...
        'TimeScale', timeScale, ...
        'VehicleScale', vehicleScale, ...
        'ShowWaypoints', showWaypoints, ...
        'Waypoints', Waypoints, ...
        'Velocity', v, ...
        'Inputs', u_data, ...
        'MaxThrust', maxThrust);
end

%% Package Output
if nargout > 0
    out.t = t;
    out.r = r;
    out.v = v;
    out.q = q;
    out.u = u_data;
    out.maxThrust = maxThrust;
    out.r_target = r_target;
    out.simOut = simOut;
end

end

%% ========================================================================
%% Helper Function: AnimateFlight (Real-Time 6-DoF 3D Flight Animation)
%% ========================================================================
function AnimateFlight(t, r, q, r_target, varargin)
% ANIMATEFLIGHT Animates 6-DoF vehicle motion following the target trajectory in real time.

p = inputParser;
addParameter(p, 'TimeScale', 1.0, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'VehicleScale', 2.5, @(x) isnumeric(x) && isscalar(x) && x > 0);
addParameter(p, 'ShowWaypoints', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'Waypoints', [], @(x) isempty(x) || isstruct(x) || isa(x, 'TOADWaypoint'));
addParameter(p, 'Velocity', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'Inputs', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'MaxThrust', 16.25, @(x) isnumeric(x) && isscalar(x));
parse(p, varargin{:});
timeScale     = p.Results.TimeScale;
scale         = p.Results.VehicleScale;
showWaypoints = logical(p.Results.ShowWaypoints);
Waypoints     = p.Results.Waypoints;
v             = p.Results.Velocity;
u_data        = p.Results.Inputs;
maxThrust     = p.Results.MaxThrust;

numPts = length(t);
if numPts < 2
    warning('AnimateFlight:NotEnoughPoints', 'Not enough points to animate.');
    return;
end

% Create Animation Figure with Dark Mode Theme
fAnim = figure('Name', 'TOAD 6-DoF Flight Animation (Real-Time)', ...
               'Color', [0.1 0.1 0.1], ...
               'NumberTitle', 'off');

ax = axes(fAnim);
set(ax, 'Color', 'k', ...
        'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'ZColor', [0.8 0.8 0.8], ...
        'GridColor', [0.35 0.35 0.35], 'GridAlpha', 0.5);

hold(ax, 'on');
grid(ax, 'on');
axis(ax, 'equal');
xlabel(ax, 'North [m]', 'Color', 'w', 'FontWeight', 'bold');
ylabel(ax, 'West [m]',  'Color', 'w', 'FontWeight', 'bold');
zlabel(ax, 'Altitude [m]', 'Color', 'w', 'FontWeight', 'bold');
view(ax, 45, 25);

% 1. Render Static Target Trajectory Curve
if ~isempty(r_target)
    plot3(ax, r_target(:,1), r_target(:,2), r_target(:,3), ...
          '--', 'Color', [1.0 0.4 0.3], 'LineWidth', 2.0, 'DisplayName', 'Target Trajectory');
end

% 2. Render Full Trajectory Outline (Ghost / History)
plot3(ax, r(:,1), r(:,2), r(:,3), ...
      ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'DisplayName', 'Flight Path');

% 3. Render Trailing Wake / Recent Trajectory History
hTrail = plot3(ax, r(1,1), r(1,2), r(1,3), ...
               '-', 'Color', [0.2 0.85 1.0], 'LineWidth', 2.2, 'DisplayName', 'Vehicle Trail');

% 4. Render Waypoints ONLY IF explicitly requested (prevents small checkpoints)
if showWaypoints && ~isempty(Waypoints)
    numWP = numel(Waypoints);
    wp_pos = zeros(3, numWP);
    for idx = 1:numWP
        wp_pos(:, idx) = Waypoints(idx).Position;
    end
    plot3(ax, wp_pos(1,:), wp_pos(2,:), wp_pos(3,:), ...
          's', 'MarkerSize', 7, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', 'w', ...
          'DisplayName', 'Waypoints');
end

% 5. Animated Moving Target Marker
if ~isempty(r_target)
    hTarget = plot3(ax, r_target(1,1), r_target(1,2), r_target(1,3), ...
                    'o', 'MarkerSize', 9, 'MarkerFaceColor', [1.0 0.25 0.2], ...
                    'MarkerEdgeColor', 'w', 'LineWidth', 1.5, 'DisplayName', 'Target Marker');
else
    hTarget = [];
end

% 6. Define Realistic 3D Rocket / Lander Geometry (No Double Cone)
% Fuselage cylinder with conical nosecone at +Z and 4 distinct landing legs at -Z.
nSides = 8;
theta_cyl = linspace(0, 2*pi, nSides + 1);
theta_cyl(end) = []; % 8 unique azimuthal angles

R_body = 0.40;   % Fuselage radius
Z_nose = 3.8;   % Nose tip (+Z)
Z_top  = 2.6;   % Upper body ring
Z_bot  = 0.8;  % Lower body ring
Z_leg  = 0;  % Leg footpad level (-Z)
R_leg  = 1.15;   % Radial leg spread

% Vertices:
% 1: Nose Tip (0, 0, Z_nose)
% 2..9: Top Ring (8 vertices)
% 10..17: Bottom Ring (8 vertices)
% 18..21: 4 Landing Leg Footpads (+X, -X, +Y, -Y)
% 22: Engine Nozzle Tip (0, 0, -1.2)

v_nose = [0, 0, Z_nose];
v_top  = [R_body * cos(theta_cyl)', R_body * sin(theta_cyl)', Z_top * ones(nSides, 1)];
v_bot  = [R_body * cos(theta_cyl)', R_body * sin(theta_cyl)', Z_bot * ones(nSides, 1)];
v_legs = [ R_leg,   0.0, Z_leg;   % Leg 1 (+X)
          -R_leg,   0.0, Z_leg;   % Leg 2 (-X)
            0.0,  R_leg, Z_leg;   % Leg 3 (+Y)
            0.0, -R_leg, Z_leg];  % Leg 4 (-Y)
v_engine = [0, 0, 0.6];

bodyVerts = [v_nose; v_top; v_bot; v_legs; v_engine];
bodyVerts = bodyVerts * scale;

% Faces Construction
faces = [];

% A. Nosecone facets (Tip to Top Ring)
for i = 1:nSides
    next_i = mod(i, nSides) + 1;
    faces = [faces; 1, 1 + i, 1 + next_i]; %#ok<AGROW>
end

% B. Fuselage Cylinder Walls (Top Ring to Bottom Ring)
for i = 1:nSides
    next_i = mod(i, nSides) + 1;
    t1 = 1 + i;
    t2 = 1 + next_i;
    b1 = 1 + nSides + i;
    b2 = 1 + nSides + next_i;
    faces = [faces; t1, b1, b2; t1, b2, t2]; %#ok<AGROW>
end

% C. Engine Nozzle (Bottom Ring to Engine Tip)
for i = 1:nSides
    next_i = mod(i, nSides) + 1;
    b1 = 1 + nSides + i;
    b2 = 1 + nSides + next_i;
    faces = [faces; 22, b1, b2]; %#ok<AGROW>
end

% D. Landing Legs (Trusses connecting Bottom Ring to Footpads)
% theta_cyl indices for +X, +Y, -X, -Y are approximately 1, 3, 5, 7
leg_joints = [1, 5, 3, 7]; % matching +X, -X, +Y, -Y
for l = 1:4
    j_idx = 1 + nSides + leg_joints(l);
    j_prev = 1 + nSides + mod(leg_joints(l)-2, nSides) + 1;
    j_next = 1 + nSides + mod(leg_joints(l), nSides) + 1;
    foot_idx = 1 + 2*nSides + l;
    faces = [faces; j_idx, foot_idx, j_prev; j_idx, foot_idx, j_next]; %#ok<AGROW>
end

% Vehicle patch
hVehicle = patch(ax, 'Vertices', bodyVerts, 'Faces', faces, ...
                 'FaceColor', [0.0 0.88 1.0], 'FaceAlpha', 0.90, ...
                 'EdgeColor', [1.0 1.0 1.0], 'LineWidth', 1.5);

% 7. Define 3D Body Axes (+X: Red, +Y: Green, +Z: Blue) Centered at Vehicle Middle
axisLen = 2.4 * scale;
vX = [axisLen, 0, 0];
vY = [0, axisLen, 0];
vZ = [0, 0, axisLen];

% Center body triad at geometric middle of vehicle body
v_mid_scaled = [0, 0, 1.9] * scale;
v_engine_scaled = [0, 0, 0.6] * scale;

% Vehicle height and max thrust vector length (half scale of vehicle at max thrust)
H_veh = (Z_nose - Z_leg) * scale;
L_max_thrust = 0.5 * H_veh; % 1.9 * scale

hAxisX = plot3(ax, [0 0], [0 0], [0 0], 'r-', 'LineWidth', 3.0, 'DisplayName', 'Body +X');
hAxisY = plot3(ax, [0 0], [0 0], [0 0], 'g-', 'LineWidth', 3.0, 'DisplayName', 'Body +Y');
hAxisZ = plot3(ax, [0 0], [0 0], [0 0], 'b-', 'LineWidth', 3.5, 'DisplayName', 'Body +Z (Nose)');

% 7b. Gimbaled Thrust Vector (Orange line, scaled with thrust command)
hThrust = plot3(ax, [0 0], [0 0], [0 0], '-', ...
                'Color', [1.0 0.45 0.0], 'LineWidth', 3.5, ...
                'DisplayName', 'Thrust Vector');

% 8. Axis Limits with Padding
allX = r(:,1); allY = r(:,2); allZ = [r(:,3); 0]; % Ensure ground z=0 is encompassed
if ~isempty(r_target)
    allX = [allX; r_target(:,1)];
    allY = [allY; r_target(:,2)];
    allZ = [allZ; r_target(:,3)];
end
spanX = max(allX) - min(allX);
spanY = max(allY) - min(allY);
spanZ = max(allZ) - min(allZ);
margin = max(2.0, 0.12 * max([spanX, spanY, spanZ]));
xlim(ax, [min(allX) - margin, max(allX) + margin]);
ylim(ax, [min(allY) - margin, max(allY) + margin]);
zlim(ax, [max(-0.5, min(allZ) - 1.0), max(allZ) + margin]);

% 9. Real-Time HUD Telemetry Banner (Translucent Overlay)
% Positioned at top to avoid collision with MATLAB figure/axes toolbar
hHUD = annotation(fAnim, 'textbox', [0.10, 0.91, 0.80, 0.07], ...
    'String', 'Initializing Real-Time Playback...', ...
    'Color', 'w', 'FontSize', 10.5, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'BackgroundColor', [0.15 0.15 0.15 0.85], 'EdgeColor', [0.4 0.4 0.4], ...
    'FitBoxToText', 'off');

title(ax, 'TOAD 6-DoF Real-Time Flight Simulation', 'Color', [0.8 0.8 0.8], 'FontSize', 11);
drawnow;
pause(0.3); % Brief settle pause

t0 = t(1);
tFinal = t(end);

currIdx = 1;
trailLen = max(10, round(50 / timeScale)); % Trailing points

startTime = tic;

while true
    % Check if window closed by user
    if ~ishghandle(fAnim) || ~isvalid(fAnim)
        break;
    end
    
    % Elapsed wall-clock time scaled by playback speed
    wallElapsed = toc(startTime) * timeScale;
    simTimeTarget = t0 + wallElapsed;
    
    if simTimeTarget >= tFinal
        currIdx = numPts;
    else
        % Fast monotonic index pointer advance (O(1) amortized)
        while currIdx < numPts && t(currIdx) < simTimeTarget
            currIdx = currIdx + 1;
        end
    end
    
    % Current vehicle state
    pos    = r(currIdx, :);
    q_curr = q(currIdx, :);
    
    % Rotate body vertices into world frame: v' = R(q) * v + pos
    rotVerts = RotateVector(bodyVerts, q_curr);
    worldVerts = rotVerts + pos;
    set(hVehicle, 'Vertices', worldVerts);
    
    % Rotate body axes centered at the geometric MIDDLE of the vehicle
    r_mid_w = pos + RotateVector(v_mid_scaled, q_curr);
    x_rot = RotateVector(vX, q_curr);
    y_rot = RotateVector(vY, q_curr);
    z_rot = RotateVector(vZ, q_curr);
    
    set(hAxisX, 'XData', [r_mid_w(1), r_mid_w(1) + x_rot(1)], ...
                'YData', [r_mid_w(2), r_mid_w(2) + x_rot(2)], ...
                'ZData', [r_mid_w(3), r_mid_w(3) + x_rot(3)]);
    set(hAxisY, 'XData', [r_mid_w(1), r_mid_w(1) + y_rot(1)], ...
                'YData', [r_mid_w(2), r_mid_w(2) + y_rot(2)], ...
                'ZData', [r_mid_w(3), r_mid_w(3) + y_rot(3)]);
    set(hAxisZ, 'XData', [r_mid_w(1), r_mid_w(1) + z_rot(1)], ...
                'YData', [r_mid_w(2), r_mid_w(2) + z_rot(2)], ...
                'ZData', [r_mid_w(3), r_mid_w(3) + z_rot(3)]);
            
    % Update Thrust Vector (Orange line from nozzle, scaled with thrust command)
    if ~isempty(u_data) && size(u_data, 1) >= currIdx
        u_k = u_data(currIdx, :);
        th_g = u_k(1);
        phi_g = u_k(2);
        T_val = u_k(3);
    else
        th_g = 0; phi_g = 0; T_val = 0;
    end
    thrustRatio = min(max(T_val / maxThrust, 0), 1.0);
    L_thrust = L_max_thrust * thrustRatio;
    
    % Engine nozzle in world coordinates
    r_nozzle = pos + RotateVector(v_engine_scaled, q_curr);
    
    % Plume direction (outward / downward from nozzle in body frame)
    d_plume_B = - [cos(th_g)*sin(phi_g), -sin(th_g), cos(th_g)*cos(phi_g)];
    d_norm = norm(d_plume_B);
    if d_norm > 1e-6
        d_plume_B = d_plume_B / d_norm;
    else
        d_plume_B = [0, 0, -1];
    end
    d_world = RotateVector(d_plume_B, q_curr);
    r_thrust_end = r_nozzle + d_world * L_thrust;
    
    set(hThrust, 'XData', [r_nozzle(1), r_thrust_end(1)], ...
                 'YData', [r_nozzle(2), r_thrust_end(2)], ...
                 'ZData', [r_nozzle(3), r_thrust_end(3)]);
            
  
    % Update Trailing Wake
    trailStart = max(1, currIdx - trailLen);
    set(hTrail, 'XData', r(trailStart:currIdx, 1), ...
                'YData', r(trailStart:currIdx, 2), ...
                'ZData', r(trailStart:currIdx, 3));
            
    % Update Target Marker
    trkErr = 0;
    if ~isempty(r_target) && ~isempty(hTarget)
        posTarg = r_target(currIdx, :);
        set(hTarget, 'XData', posTarg(1), 'YData', posTarg(2), 'ZData', posTarg(3));
        trkErr = norm(pos - posTarg);
    end
    
    % Update HUD Banner Text
    if ~isempty(v) && size(v, 1) >= currIdx
        spd = norm(v(currIdx, :));
    else
        spd = 0;
    end
    thrustPct = thrustRatio * 100;
    hudStr = sprintf('Time: %.2f / %.2f s   |   Alt: %.2f m   |   Thrust: %.0f%%   |   Speed: %.2f m/s   |   Error: %.2f m   |   [%.1fx RT]', ...
                     t(currIdx), tFinal, pos(3), thrustPct, spd, trkErr, timeScale);
    if ishghandle(hHUD) && isvalid(hHUD)
        hHUD.String = hudStr;
    end
    
    drawnow limitrate;
    
    if currIdx >= numPts
        break;
    end
end

if ishghandle(fAnim) && isvalid(fAnim)
    if ishghandle(hHUD) && isvalid(hHUD)
        hHUD.String = sprintf('Flight Playback Complete (Final Time: %.2f s | Final Error: %.2f m)', t(end), trkErr);
        hHUD.Color = [0.4 1.0 0.4];
    end
    drawnow;
end

end

%% ========================================================================
%% Helper Function: RotateVector (Quaternion Rotation)
%% ========================================================================
function v_rot = RotateVector(v, q)
% ROTATEVECTOR Rotates 3D vector(s) v (Nx3) by quaternion q [w, x, y, z].
% Uses the standard Rodrigues formula: v' = v + 2w(u x v) + 2(u x (u x v))

w = q(1);
u = q(2:4);

numV = size(v, 1);
u_rep = repmat(u, numV, 1);

uv = cross(u_rep, v, 2);
uuv = cross(u_rep, uv, 2);

v_rot = v + 2 * w * uv + 2 * uuv;
end

%% ========================================================================
%% Helper Function: formatData (Matrix Dimension Alignment)
%% ========================================================================
function out = formatData(data, timeLen)
% FORMATDATA Robustly aligns time series data to [timeLen x numStates].
sqData = squeeze(data);
sz = size(sqData);

if sz(1) == timeLen
    out = sqData;
elseif sz(2) == timeLen
    out = sqData';
else
    nCols = round(numel(data) / timeLen);
    if nCols * timeLen == numel(data)
        out = reshape(data, nCols, timeLen)';
    else
        error('formatData:DimensionMismatch', ...
            'Data dimension mismatch: Neither dimension matches time vector length %d.', timeLen);
    end
end
end

%% ========================================================================
%% Helper Function: ConfigureToWorkspaceBlocks
%% ========================================================================
function ConfigureToWorkspaceBlocks(modelName, sampleTime)
% CONFIGURETOWORKSPACEBLOCKS Sets sample time and limits on To Workspace blocks,
% and attaches input_log in memory if missing.

% 1. Ensure input_log exists to log commanded control inputs
try
    foundInp = find_system(modelName, 'SearchDepth', 1, 'Name', 'input_log');
    if isempty(foundInp)
        fromBlk = [modelName '/From_INP_Log'];
        toWsBlk = [modelName '/input_log'];
        add_block('simulink/Signal Routing/From', fromBlk, 'GotoTag', 'INP');
        add_block('simulink/Sinks/To Workspace', toWsBlk, ...
                  'VariableName', 'input_log', ...
                  'SaveFormat', 'Timeseries', ...
                  'SampleTime', num2str(sampleTime), ...
                  'MaxDataPoints', 'inf', ...
                  'Decimation', '1');
        add_line(modelName, 'From_INP_Log/1', 'input_log/1');
    end
catch ME
    warning('Plot6DoF:InputLogConfigWarning', 'Could not configure input_log block: %s', ME.message);
end

% 2. Configure sample times on target logging blocks
targetBlocks = {'state_log', 'target_pos_log', 'input_log'};
for i = 1:length(targetBlocks)
    blkName = targetBlocks{i};
    bPath = [modelName '/' blkName];
    try
        found = find_system(modelName, 'SearchDepth', 1, 'Name', blkName);
        if ~isempty(found)
            set_param(bPath, 'SampleTime', num2str(sampleTime));
            set_param(bPath, 'MaxDataPoints', 'inf');
            set_param(bPath, 'Decimation', '1');
            set_param(bPath, 'SaveFormat', 'Timeseries');
        end
    catch ME
        warning('Plot6DoF:ConfigWarning', 'Could not configure block %s: %s', bPath, ME.message);
    end
end
end
