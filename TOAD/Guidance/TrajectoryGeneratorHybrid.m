%% TRAJECTORYGENERATORHYBRID
% Master 6-DoF trajectory setup, generation, and TV-LQI gain export script.
% Powered by the Universal Two-Stage Hybrid Trajectory Engine:
%   Stage 1: Minimum-snap polynomial trajectory optimization (PolyTrajectoryQP)
%            with corridor enforcement & differential-flatness state inversion.
%   Stage 2: 6-DoF direct collocation refinement (IPOPT / CasADi) with
%            actuator chatter suppression and takeoff/landing flight corridors.
%
% Supported Vehicles:
%   - TOAD    (Vehicle = 1): High-thrust liquid bipropellant lander
%   - ASTRAv2 (Vehicle = 0): Precision electric ducted-fan / TVC testbed
%
% Supported Maneuver Presets:
%   - "Hop":      Point-to-point parabolic trajectory with vertical liftoff/touchdown
%   - "Circle":   Climb, horizontal orbit inspection loop, and flared descent
%   - "Backflip": High-altitude 360-degree pitch inversion maneuver
%   - "Waypoint": Multi-target 3D survey route with intermediate corridor gates
%   - "Custom":   Arbitrary user-defined 3D boundary conditions & waypoints
%
% Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

clear; clc; close all;

%% =========================================================================
%% 1. Mission Configuration
%% =========================================================================
Vehicle     = 1;            % 1 = TOAD, 0 = ASTRA
Maneuver    = "Backflip";     % "Hop", "Circle", "Backflip", "Waypoint", "Custom"
Version     = 1;            % Output version integer: formats as v001, v002, etc.

% Discretization & Mesh
N_nodes     = 60;           % Control intervals (60-80 recommended for Hybrid engine)
T_initial   = [];           % Optional duration guess [s] (leave empty [] for auto-physics)

% Position Boundaries [m] (East, North, Up)
r_launch    = [0; 0; 0];    % Launch pad position
r_target    = [0; 0; 0];    % Target touchdown position (set per maneuver below)

%% =========================================================================
%% 2. Vehicle Constants & Project Paths
%% =========================================================================
% Ensure CasADi optimization suite is available on search path
if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir') && isempty(which('casadi.Opti'))
    addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
end

% Load vehicle mass, inertia, and actuator parameters
constants6DoF = LoadTOADParams(Vehicle);

% Target directory for trajectory CSV export
save_dir = fullfile(pwd, 'Guidance', 'Trajectories');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

isTOAD = (Vehicle == 1 || strcmpi(string(Vehicle), "TOAD"));
if isTOAD
    veh_name = "TOAD";
else
    veh_name = "ASTRA";
end

fprintf('================================================================================\n');
fprintf('  PSP ACTIVE CONTROLS - UNIVERSAL 6-DoF TRAJECTORY GENERATOR (HYBRID ENGINE)     \n');
fprintf('  Vehicle: %s | Maneuver: %s | Version: v%03d | Nodes: %d                       \n', ...
        veh_name, Maneuver, Version, N_nodes);
fprintf('================================================================================\n\n');

%% =========================================================================
%% 3. Instantiate & Configure Optimizer
%% =========================================================================
opt = TrajectoryOptimizerHybrid(constants6DoF, ...
    'Vehicle',         Vehicle, ...
    'Maneuver',        Maneuver, ...
    'Version',         Version, ...
    'N',               N_nodes, ...
    'SaveDir',         save_dir, ...
    'GlideslopeAngle', 10, ...
    'FunnelCurvature', 0.015, ...
    'MaxIter',         150, ...
    'Tol',             1.5e-2, ...
    'PrintLevel',      0);

if ~isempty(T_initial)
    opt.T_initial = T_initial;
end

%% =========================================================================
%% 4. Configure Maneuver-Specific Waypoints & Goals
%% =========================================================================
switch Maneuver
    case "Hop"
        % Parabolic hop: vertical climb, lateral translation, flared descent
        if isTOAD
            r_target = [10; 0; 0];      % 10 m East touchdown (or [50; 0; 0] for long hop)
            apex_alt = 50.0;            % Target apex altitude [m]
        else
            r_target = [8; 0; 0];       % 8 m East touchdown
            apex_alt = 15.0;            % Target apex altitude [m]
        end
        opt.setBoundaries(r_launch, r_target);
        opt.setManeuver('Hop', 'apex_alt', apex_alt);

    case "Circle"
        % Orbit maneuver: vertical climb, circular survey loop, flared touchdown
        opt.setBoundaries(r_launch, r_target);
        if isTOAD
            opt.setManeuver('Circle', 'circle_radius', 15.0, 'circle_alt', 35.0, 'circle_center', [0; 0]);
        else
            opt.setManeuver('Circle', 'circle_radius', 5.0,  'circle_alt', 7.0,  'circle_center', [0; 0]);
        end

    case "Backflip"
        % Backflip maneuver: vertical climb, full 360-deg pitch flip, recovery landing
        opt.setBoundaries(r_launch, r_target);
        if isTOAD
            opt.setManeuver('Backflip', 'apex_alt', 45.0, 'flip_start_frac', 0.35, 'flip_end_frac', 0.65);
        else
            opt.setManeuver('Backflip', 'apex_alt', 20.0, 'flip_start_frac', 0.35, 'flip_end_frac', 0.65);
        end

    case "Waypoint"
        % Multi-waypoint 3D survey route
        if isTOAD
            wps = [ 0.0,  15.0,  30.0,  15.0,   0.0; ...
                    0.0,  15.0,   0.0, -15.0,   0.0; ...
                    0.0,  35.0,  50.0,  35.0,   0.0];
            t_wp = 28.0;
        else
            wps = [ 0.0,   4.0,   8.0,   4.0,   0.0; ...
                    0.0,   4.0,   0.0,  -4.0,   0.0; ...
                    0.0,   8.0,  12.0,   8.0,   0.0];
            t_wp = 18.0;
        end
        opt.setWaypoints(wps, 'T_total', t_wp, 'Tolerances', 0.80);

    case "Custom"
        % User custom boundaries and parameters
        opt.setBoundaries(r_launch, r_target);

    otherwise
        error('TrajectoryGeneratorHybrid:UnknownManeuver', 'Unrecognized maneuver: %s', Maneuver);
end

%% =========================================================================
%% 5. Solve Optimal Trajectory (Stage 1 QP -> Stage 2 Collocation)
%% =========================================================================
fprintf('Starting two-stage optimization pipeline...\n');
tic;
sol = opt.solve();
t_solve = toc;

%% =========================================================================
%% 6. Post-Processing, Figures, & Data Export
%% =========================================================================
if strcmp(sol.Status, 'Success')
    fprintf('\n>>> OPTIMAL TRAJECTORY FOUND in %.2f s! <<<\n', t_solve);
    fprintf('    Stage 1 (QP Polynomial): %.2f ms\n', sol.stats.t_stage1_ms);
    fprintf('    Stage 2 (NLP Refine):   %.2f s (%d iterations)\n', sol.stats.t_stage2_s, sol.stats.iter_count);
    fprintf('    Mission Duration:       %.2f s\n', sol.T_total);
    fprintf('    Flight Path Length:     %.2f m\n', sol.L_path);
    fprintf('    Discretization:         %d intervals (dt = %.3f s)\n', opt.N, sol.T_total / opt.N);

    % Generate publication dashboard figure with glideslope cones & landing funnels
    opt.plot();

    % Export trajectory table to CSV
    [~, csv_path] = opt.exportCSV();
    [~, traj_stem] = fileparts(csv_path);

    % Pre-generate corresponding TV-LQI tracking gains via SaveGains
    try
        fprintf('\nGenerating TV-LQI tracking gains for %s...\n', traj_stem);
        SaveGains(traj_stem, constants6DoF);
        fprintf('Gain matrices successfully generated and saved to Guidance/Trajectories/Gains/\n');
    catch ME_gains
        fprintf('Note on Gain Generation: %s\n', ME_gains.message);
    end

    fprintf('\nReady for flight simulation! In LoadTOADSim.m, set:\n');
    fprintf('  filename = "%s";\n\n', traj_stem);
else
    error('TrajectoryGeneratorHybrid:OptimizationFailed', ...
        'Solver did not achieve full optimal convergence (Status: %s). Check opt.Solution for details.', sol.Status);
end
