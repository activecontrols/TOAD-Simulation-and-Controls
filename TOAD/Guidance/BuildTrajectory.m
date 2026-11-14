%% BUILDTRAJECTORY
% Top-level 6-DoF trajectory setup, closed-form Polynomial QP generation,
% 20-column CSV export, TV-LQI gain generation, and 4-figure docked dashboard.

clear; clc; close all;

%% Mission & Solver Config
Vehicle     = 0;            % 1 = TOAD, 0 = ASTRA
MissionName = "Waypoint";   % Mission identifier for CSV output naming
Version     = 1;            % Output version integer: formats as v001, v002, etc.

Opts = struct('R', 5);      % Minimum-jounce
ExportCSVAndGains = true;   
ShowPlots         = true;   

%% Load Vehicle Constants
constants6DoF = LoadTOADParams(Vehicle);
if Vehicle == 1
    vehName = "TOAD";
else
    vehName = "ASTRA";
end

fprintf('================================================================================\n');
fprintf('  PSP ACTIVE CONTROLS - 6-DoF POLYNOMIAL QP TRAJECTORY GENERATOR (PolyQP v2)    \n');
fprintf('  Vehicle: %-6s | Mission: %-10s | Version: v%03d | Order: R = %d            \n', ...
    vehName, MissionName, Version, Opts.R);
fprintf('================================================================================\n\n');

%% Define Mission Waypoints
dummyWP   = TOADWaypoint([0; 0; 0]);
Waypoints = repmat(dummyWP, 1, 6);

Waypoints(1) = TOADWaypoint([0; 0; 0], ...
    'Label', 'PAD HLD 1', 'HoldTime', 10.0);

Waypoints(2) = TOADWaypoint([0; 0; 2], ...
    'Label', 'FLT 1', 'MaxVel', [5; 5; 5], 'StraightLine', true);

Waypoints(3) = TOADWaypoint([10; 4; 15], ...
    'Label', 'FLT 2', 'MaxVel', [5; 5; 5], 'VelBC', [NaN, NaN, 0]);

Waypoints(4) = TOADWaypoint([10; 10; 10], ...
    'Label', 'FLT 3', 'MaxVel', [5; 5; 5]);

Waypoints(5) = TOADWaypoint([10; 10; 2], ...
    'Label', 'FLT 4', 'HoldTime', 1.0, 'MaxVel', [5; 5; 5]);

Waypoints(6) = TOADWaypoint([10; 10; 0], ...
    'Label', 'FLT END', 'HoldTime', 1.0, 'MaxVel', [3; 3; 1], 'StraightLine', true);

%% Polynomial Trajectory & Flatness Inversion
tic;
[InitialGuess, QPTraj, Report] = PolyQP(Waypoints, constants6DoF, Opts);
tSolveMs = toc * 1000;

dr = diff(InitialGuess.X(5:7, :), 1, 2);
pathLen = sum(vecnorm(dr, 2, 1));
mStart  = constants6DoF.m_dry + InitialGuess.X(14, 1)   + InitialGuess.X(15, 1);
mEnd    = constants6DoF.m_dry + InitialGuess.X(14, end) + InitialGuess.X(15, end);

fprintf('>>> POLYQP v2 TRAJECTORY SOLVED in %.2f ms (Feasible = %d) <<<\n', tSolveMs, Report.Feasible);
fprintf('    Segments:            %d (%d Hold, %d Cruise, %d Flight)\n', ...
    QPTraj.Count, sum(QPTraj.IsHold(1:QPTraj.Count)), sum(QPTraj.IsCruise(1:QPTraj.Count)), ...
    QPTraj.Count - sum(QPTraj.IsHold(1:QPTraj.Count)) - sum(QPTraj.IsCruise(1:QPTraj.Count)));
fprintf('    Outer Dilation Iter: %d\n', Report.OuterIterations);
fprintf('    Mission Duration:    %.2f s (%d samples @ dt = %.3f s)\n', ...
    QPTraj.T_total, numel(InitialGuess.Time), InitialGuess.dt_row(1));
fprintf('    Flight Path Length:  %.2f m\n', pathLen);
fprintf('    Vehicle Mass:        %.3f kg -> %.3f kg (Burned: %.3f kg)\n', mStart, mEnd, mStart - mEnd);
fprintf('    Peak Utilization:    %.1f%% (Vel %.1f%%, Acc %.1f%%, ThrHi %.1f%%, ThrLo %.1f%%, Gim %.1f%%, Tilt %.1f%%)\n\n', ...
    100 * Report.MaxOverallUtil, 100 * Report.UtilVel, 100 * Report.UtilAccel, ...
    100 * Report.UtilThrustHi, 100 * Report.UtilThrustLo, 100 * Report.UtilGimbal, 100 * Report.UtilTilt);

%% Export CSV & Pre-Generate TV-LQI Gains
if ExportCSVAndGains
    [~, csvPath, trajStem] = ExportTrajectoryCSV(InitialGuess, constants6DoF, ...
        'MissionName', MissionName, ...
        'Version',     Version, ...
        'RunSaveGains', true);

    fprintf('\nReady for flight simulation! In LoadTOADSim.m, set:\n');
    fprintf('    filename = "%s";\n\n', trajStem);
end

%% Plots
if ShowPlots
    PlotTrajectoryQP(InitialGuess, QPTraj, Report, Waypoints, constants6DoF);
end
