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
Waypoints(1) = TOADWaypoint([0; 0; 0], 'Label', 'PAD HLD 1', 'HoldTime', 10.0, ...
    'PosTol',  [1; 1; 1]);

Waypoints(2) = TOADWaypoint([0; 0; 5], 'Label', 'PAD DIV 1', 'MaxVel', ...
    [3; 3; 6], 'PosTol', [4; 4; 1], 'VelFF', [0; 0; 6], 'MaxTime', 5);

Waypoints(3) = TOADWaypoint([10; 10; 50], 'Label', 'FLT HOV 1', 'PosTol', [1; 1; 1], ...
    'HoldTime', 11, 'MaxTime', 25, 'HDGRef', [1; 0; 0]);

Waypoints(4) = TOADWaypoint([10; 10; 15], 'Label', 'FLT HOV 2', 'PosTol', [2; 2; 3], ...
    'MaxVel', [3; 3; 6], 'VelFF', [0; 0; -4]);

Waypoints(5) = TOADWaypoint([10; 10; 0], 'Label', 'FLT END 1', 'PosTol', [3; 3; 0.5], ...
    'HoldTime', 1, 'MaxVel', [3; 3; 4], 'VelFF', [0; 0; -0.1]);