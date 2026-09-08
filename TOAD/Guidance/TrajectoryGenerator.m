%% TRAJECTORYGENERATOR
% Master trajectory setup, generation, and export script for TOAD & ASTRAv2.
% Uses the TrajectoryOptimizer engine to generate 6-DoF optimal trajectories
% and export them to Guidance/Trajectories/ for Simulink and TV-LQI tracking.
%
% Authors: PSP Active Controls (Pablo Plata, Andrew Lulo, & Antigravity)

clear; clc; close all;

%% Mission Configuration
Vehicle     = 1;                % Vehicle model: 1 for TOAD, 0 for ASTRA (legacy strings "TOAD", "ASTRAv2" also supported)
Maneuver    = "Hop";            % Maneuver preset: "Backflip", "Circle", "Hop", "Custom"
Version     = 1;                % Version integer: formats as v001, v002, etc.

% Discretization & Mesh
N_nodes     = 100;              % Number of control intervals (80 - 200 recommended)
T_initial   = 35;               % Initial duration guess [s]

% Position Boundaries [m] (East, North, Up)
r_launch    = [0; 0; 0];        % Launch pad position
r_target    = [0; 0; 0];        % Target landing touchdown position

%% Initialize Simulation Parameters & Project Paths
project_root = '';
p_scan = fileparts(mfilename('fullpath'));
while ~isempty(p_scan)
    if exist(fullfile(p_scan, 'LoadTOADSim.m'), 'file') && exist(fullfile(p_scan, 'Navigation'), 'dir')
        project_root = p_scan;
        break;
    end
    parent_scan = fileparts(p_scan);
    if strcmp(parent_scan, p_scan), break; end
    p_scan = parent_scan;
end
if isempty(project_root)
    project_root = pwd;
end
cd(project_root);
addpath(genpath(project_root));

% Pre-set vehicle in constants6DoF so LoadTOADSim loads matching parameters
constants6DoF.Vehicle = Vehicle;
LoadTOADSim;

% Target directory for trajectory CSV export
save_dir    = fullfile(project_root, 'Guidance', 'Trajectories');

%% Instantiate & Configure TrajectoryOptimizer
% Normalize vehicle for comparison
isTOAD = (Vehicle == 1 || strcmpi(string(Vehicle), "TOAD"));
if isTOAD
    veh_name = "TOAD";
else
    veh_name = "ASTRA";
end

fprintf('  PSP ACTIVE CONTROLS - 6-DoF TRAJECTORY GENERATOR\n\n');
fprintf('  Vehicle: %s | Maneuver: %s | Version: v%03d\n', veh_name, Maneuver, Version);

opt = TrajectoryOptimizer(constants6DoF, ...
    'Vehicle',   Vehicle, ...
    'Maneuver',  Maneuver, ...
    'Version',   Version, ...
    'N',         N_nodes, ...
    'T_initial', T_initial, ...
    'SaveDir',   save_dir, ...
    'max_iter',  1000);

opt.setBoundaries(r_launch, r_target);

%% Configure Maneuver-Specific Waypoints & Goals
switch Maneuver
    case "Backflip"
        opt.setBoundaries(r_launch, r_target);
        if isTOAD
            opt.setManeuver('Backflip', 'apex_alt', 50.0);
        else
            opt.setManeuver('Backflip', 'apex_alt', 20.0);
        end
        
    case "Circle"
        opt.setBoundaries(r_launch, r_target);
        if isTOAD
            opt.setManeuver('Circle', 'circle_radius', 15.0, 'circle_alt', 35.0);
        else
            opt.setManeuver('Circle', 'circle_radius', 5.0, 'circle_alt', 12.0);
        end
        
    case "Hop"
        opt.setBoundaries(r_launch, [10; 0; 0]);
        if isTOAD
            opt.setManeuver('Hop', 'apex_alt', 50.0);
        else
            opt.setManeuver('Hop', 'apex_alt', 15.0);
        end
end

%% Inspect Initial Guess Before Solving
% figure_guess = opt.plotInitialGuess();

%% Solve Optimal Trajectory
tic;
sol = opt.solve();
t_solve = toc;

%% Post-Processing, Figures, & Data Export
if strcmp(sol.Status, 'Success')
    fprintf('\n>>> OPTIMAL TRAJECTORY FOUND in %.2f s! <<<\n', t_solve);
    fprintf('    Mission Duration: %.2f s\n', sol.T_total);
    fprintf('    Discretization:   %d intervals (dt = %.3f s)\n', opt.N, sol.T_total / opt.N);
    
    % Generate 12-panel mission dashboard
    opt.plot();
    
    % Export CSV to Guidance/Trajectories/<Vehicle>_<Maneuver>_v###.csv
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
    
    fprintf('\nReady for simulation! In LoadTOADSim.m, set:\n');
    fprintf('  filename = "%s";\n\n', traj_stem);
else
    error('TrajectoryGenerator:OptimizationFailed', ...
        'Solver did not achieve full optimal convergence (Status: %s). Check opt.Solution for debug values.', sol.Status);
end
