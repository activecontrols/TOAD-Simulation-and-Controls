function Waypoints = TrajectoryBuilder()
% TRAJECTORYBUILDER  Construct the ordered TOADWaypoint mission array.
%
% Returns a 1xN array of TOADWaypoint objects. Waypoints are evaluated
% in order by CheckpointCommands. Each waypoint transitions when its
% CriteriaMet() method returns true AND any required hold time elapses.

% Pre allocate
dummyWP = TOADWaypoint([0; 0; 0]);
Waypoints = repmat(dummyWP, 1, 5);

% Trajectory
Waypoints(1) = TOADWaypoint([0; 0; 0], 'Label', 'PAD HLD 1', 'HoldTime', 20.0, ...
    'PosTol',   1);

Waypoints(2) = TOADWaypoint([0; 0; 5], 'Label', 'PAD DIV 1', 'MaxVel', ...
    [3; 3; 6], 'PosTol', 1, 'VelFF', [0; 0; 6]);

Waypoints(3) = TOADWaypoint([0; 5; 50], 'Label', 'FLT HOV 1', 'PosTol', 1, ...
    'HoldTime', 10);

Waypoints(4) = TOADWaypoint([0; 5; 35], 'Label', 'FLT DIV 2', 'MaxVel', ...
    [3; 3; 6], 'PosTol', 1, 'VelFF', [0; 0; -6]);

Waypoints(5) = TOADWaypoint([0; 10; 0], 'Label', 'FLT END 1', 'PosTol', 1, ...
    'HoldTime', 10, 'MaxVel', [3; 3; 6]);