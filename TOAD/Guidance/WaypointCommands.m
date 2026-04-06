function [trg, MaxVel, VelFF] = WaypointCommands(x, t)

    % Store full 13-state vector before slicing for CriteriaMet evaluation
    X_full = x; 
    ABORT = 0;

    % Persistents
    persistent Waypoints
    persistent timeFlag
    persistent i
    persistent timeCounter
    persistent prevTime
    persistent HoldMode
    persistent errorAccum

    if isempty(timeFlag)
        % Initialize the mission waypoints once at the start of the simulation
        Waypoints = TrajectoryBuilder();
        
        timeFlag = 999;
        i = 1;
        timeCounter = 0;
        prevTime = 0;
        HoldMode = 0;
        errorAccum = zeros(3, 1);
    end
    dt = t - prevTime;
    prevTime = t;

    % Ignores lateral position gain if time is past set abort value
    % ABORT MODE LOGIC: If artificial abort is triggered, enter abort mode
    % loop. If the position hold mode is disabled, set lateral velocity
    % references to zero by setting their gains to zero. If lateral
    % velocities are below threshold, pick current position as hold and
    % activate HoldMode.

    isABORT = ABORT > 0 && t >= ABORT;
    
    % Retrieve current waypoint from array
    wp = Waypoints(i);
    currentPosTarget = wp.Position;
    
    if isABORT
        if HoldMode == 0
            currentPosTarget = zeros(3,1);
        end
        if norm(x(7:8)) < 0.2 && HoldMode == 0
            currentPosTarget = [x(5:6); 0];
            HoldMode = 1;
        end
        % Build a temporary waypoint to safely bypass CriteriaMet during abort
        wp = TOADWaypoint(currentPosTarget, 'PosTol', 1.0);
    end
    
    PosError = currentPosTarget - x(5:7);
    MaxVel = wp.MaxVel;
    if isABORT == 0
        VelFF = wp.VelFF;
    else
        VelFF = zeros(3,1);
    end

    % Advance logic
    if isABORT == 0 
        % Timer only accumulates if position is within tolerance
        if CriteriaMet(wp, X_full)
            timeCounter = timeCounter + dt;
            errorAccum = errorAccum + PosError * dt;
        end
        
        % Advance based solely on position & hold time
        advanceNow = CriteriaMet(wp, X_full) && ...
                     (wp.IsPassAndGo || timeCounter >= wp.HoldTime);
        
        if advanceNow && i < numel(Waypoints)
            i = i + 1;
            timeCounter = 0;
        end
    end
    
    trg = currentPosTarget;
end
function met = CriteriaMet(wp, X)
    met = norm(X(5:7) - wp.Position) < wp.PosTol;
end