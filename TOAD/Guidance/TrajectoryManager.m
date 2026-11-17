function [X, U, K, LA, LT] = TrajectoryManager(t, X_est, constantsTOAD)
    persistent Time
    persistent States
    persistent Inputs
    persistent tAbort

    if isempty(Time)
        Time = constantsTOAD.Traj.Time;
        States = constantsTOAD.Traj.States;
        Inputs = constantsTOAD.Traj.Inputs;
        tAbort = 0;
    end

    ABORT = 0; 
    % Persist abort calculation state across function calls
    persistent abortCalculated
    if isempty(abortCalculated)
        abortCalculated = false;
    end

    if ABORT > 0 && t > ABORT && ~abortCalculated
        Opts = struct('R', 5);      % Minimum-jounce
        dummyWP   = TOADWaypoint([0; 0; 0]);
        Waypoints = repmat(dummyWP, 1, 3);

        Waypoints(1) = TOADWaypoint(X_est(5:7), 'VelBC', X_est(8:10), 'Label', 'ABT START');
        Waypoints(2) = TOADWaypoint([nan; nan; 2], 'Label', 'ABT GND', 'HoldTime', 0.0, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]);
        Waypoints(3) = TOADWaypoint([nan; nan; 0], 'Label', 'ABT END', 'HoldTime', 0.0, 'MaxVel', [3; 3; 1], 'StraightLine', true, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]);
        [InitialGuess, QPTraj, Report] = PolyQP(Waypoints, constantsTOAD, Opts);  

        Time = InitialGuess.Time(:);
        States = InitialGuess.X';
        Inputs = [InitialGuess.U'; InitialGuess.U(:, end)'];

        % Turn off abort to avoid recalculating and reference time from
        % when the abort was called
        abortCalculated = true;
        tAbort = t;
    end

    % FBGain maps to K_trans, FBCost maps to K_rot
    K_trans_Gain = constantsTOAD.Traj.KTGain;
    K_rot_Gain = constantsTOAD.Traj.KRGain; 
    LA_Gain = constantsTOAD.Traj.LAGain;
    LT_Gain = constantsTOAD.Traj.LTGain;

    X = zeros(15, 1);
    XTEST = zeros(3, 1);
    U = zeros(4, 1);
    K = zeros(6, 6);
    LA = zeros(6, 3);
    LT = zeros(9, 6);

    K(1:3, 1:6) = K_trans_Gain(end, :, :);
    K(4:6, 1:6) = K_rot_Gain(end, :, :);

    LA(:) = LA_Gain(end, :, :);
    LT(:) = LT_Gain(end, :, :);
    
    %% Handle Boundary Conditions
    if t - tAbort <= Time(1) || t - tAbort == 0
        X(:) = States(1,:);
        U(:) = Inputs(1,:);
        return;
    elseif t - tAbort >= Time(end)
        X(:) = States(end,:);
        X(1:4) = [1; 0; 0; 0];
        X(8:end) = 0;
        total_mass = constantsTOAD.m_dry + X_est(14,1) + X_est(15,1);
        U(:) = [0; 0; total_mass * constantsTOAD.g; 0];
        return;
    end
    
    n_low = 1;
    for i = 1:(length(Time) - 1)
        if Time(i) <= (t - tAbort) && Time(i+1) > (t - tAbort)
            n_low = i;
            break;
        end
    end    
    
    n_high = n_low + 1;
    T_low = Time(n_low);     
    T_high = Time(n_high);
    dt = T_high - T_low;
    ratio = (t - tAbort - T_low) / dt;
    
    %% Interpolation
    X_low = States(n_low,:);
    X_up  = States(n_high,:);
    X(:) = X_low + (X_up - X_low) .* ratio;
    
    % Quaternion normalization to eliminate norm shrinkage and DCM defect
    q_norm = norm(X(1:4));
    if q_norm > 1e-6
        X(1:4) = X(1:4) / q_norm;
    end
    
    U_low = Inputs(n_low,:);
    U_up  = Inputs(n_high,:);
    U(:) = U_low + (U_up - U_low) .* ratio;
end