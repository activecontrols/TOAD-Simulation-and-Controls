function S = AbortGuidance(t, X_est, constantsTOAD, ABORT, S, ABORT_HOVER_ALT, u_cmd, a_dist)
% Abort guidance: the abort block of Guidance/TrajectoryManager.m from branch abort-fixes, with the start-of-abort
% acceleration / jerk boundary condition from attitude and thrust (FIX-F02, ASTRA integration Phase 5).
% TrajectoryManager.m was later replaced by the interpolating trajectory manager (commit fd08005), which has
% no abort block, so the block lives here. It is the MATLAB reference that cpp/PolyQP/AbortGuidance.cpp transpiles.
% The block's persistents (Time, States, Inputs, tAbort, abortCalculated, abortStatus) are carried in S.
% ABORT / ABORT_HOVER_ALT default to the constants the block had (0, 2). u_cmd is the last command [gimbal pitch; gimbal yaw;
% thrust (N); roll] and a_dist the estimated disturbance acceleration (LESO_Position); they set the start acceleration and
% jerk (default: hover thrust with the gimbal centred, no disturbance).
    if nargin < 4
        ABORT = 0;
    end
    if nargin < 6
        ABORT_HOVER_ALT = 2; % DECISION-D5 (TRANSPILE-FLAG F01): abort waypoint 2 altitude (m); the lead's spec says 1 m above the pad
    end
    if nargin < 7 || isempty(u_cmd)
        u_cmd = [0; 0; constantsTOAD.m_dry * constantsTOAD.g; 0]; % FIX-F02: default command: hover thrust, gimbal centred, no roll
    end
    if nargin < 8 || isempty(a_dist)
        a_dist = [0; 0; 0]; % FIX-F02: no disturbance estimate
    end
    if nargin < 5 || isempty(S)
        S = struct('Time', [], 'States', [], 'Inputs', [], 'tAbort', 0, 'abortCalculated', false, 'abortStatus', []);
    end
    Time = S.Time; States = S.States; Inputs = S.Inputs; tAbort = S.tAbort;
    % Persist abort calculation state across function calls
    abortCalculated = S.abortCalculated;
    abortStatus = S.abortStatus; % FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    if isempty(abortStatus)
        abortStatus = 0; % FIX-F08: 0 = OK, 1 = INFEASIBLE_FLOWN (D2), 2 = FALLBACK (D1)
    end

    if ABORT > 0 && t > ABORT && ~abortCalculated
        % FIX-F29: plan from the estimated propellant
        Opts = struct('R', 5, 'InitOxMass', X_est(14), 'InitFuMass', X_est(15));      % Minimum-jounce
        dummyWP   = TOADWaypoint([0; 0; 0]);
        Waypoints = repmat(dummyWP, 1, 3);

        % FIX-F02: body thrust vector from the last command, as LESO_Position computes it: T * [cos(th) sin(ph); -sin(th); cos(th) cos(ph)]
        q_est = X_est(1:4) / norm(X_est(1:4));
        C_IB = quatRot(q_est)';
        T_B = u_cmd(3) * [cos(u_cmd(1)) * sin(u_cmd(2)); -sin(u_cmd(1)); cos(u_cmd(1)) * cos(u_cmd(2))];
        w_I = C_IB * X_est(11:13);
        % FIX-F02: start acceleration a0 = C_IB * T_B / m + g + a_dist (g = [0; 0; -g], z up; m = m_dry, ASTRA carries no propellant)
        a0 = C_IB * T_B / constantsTOAD.m_dry + [0; 0; -constantsTOAD.g] + a_dist;
        % FIX-F02: start jerk j0 = (w_I x C_IB * T_B) / m with w_I = C_IB * w_B: the thrust is held constant over the abort instant (DECISION-D6)
        j0 = cross(w_I, C_IB * T_B) / constantsTOAD.m_dry;
        % FIX-F03: thrust axis b3 = C_IB * T_B / |T_B| (the body z axis when the gimbal is centred) and db3/dt = w_I x b3, consistent with a0 and j0
        if norm(T_B) > 0
            b3_est = C_IB * T_B / norm(T_B);
        else
            b3_est = C_IB * [0; 0; 1];
        end
        b3dot_est = cross(w_I, b3_est);
        % DECISION-D9: the heading reference is the estimated heading, x_C = C_IB(:, 2) x e3, so the plan starts at the estimated attitude
        hdg_est = cross(C_IB(:, 2), [0; 0; 1]);
        if norm(hdg_est) < 1e-6
            hdg_est = [1; 0; 0];
        end
        Waypoints(1) = TOADWaypoint(X_est(5:7), 'VelBC', X_est(8:10), 'AccelBC', a0, 'JerkBC', j0, 'Label', 'ABT START', 'B3BC', b3_est, 'B3DotBC', b3dot_est, 'HDGRef', hdg_est); % FIX-F02/F03: knot-1 state, acceleration, jerk and attitude
        % FIX-F02: the end knots keep AccelBC = 0 and JerkBC = 0: the a0/j0 formula at hover (C_IB = I, T_B = m g e3, a_dist = 0, w = 0) gives exactly these
        Waypoints(2) = TOADWaypoint([nan; nan; ABORT_HOVER_ALT], 'Label', 'ABT GND', 'HoldTime', 0.0, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0], 'HDGRef', hdg_est);
        Waypoints(3) = TOADWaypoint([nan; nan; 0], 'Label', 'ABT END', 'HoldTime', 0.0, 'MaxVel', [3; 3; 1], 'StraightLine', true, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0], 'HDGRef', hdg_est);
        % FIX-F08: a non-finite start state cannot be planned; it goes straight to the fallback
        planOk = all(isfinite(X_est(5:10)));
        if planOk
            if coder.target('MATLAB') % FIX-F08: code generation has no try/catch, so a PolyQP error is caught only in interpreted MATLAB
                try
                    [InitialGuess, QPTraj, Report] = PolyQP(Waypoints, constantsTOAD, Opts);
                catch
                    planOk = false;
                end
            else
                [InitialGuess, QPTraj, Report] = PolyQP(Waypoints, constantsTOAD, Opts);
            end
            if planOk
                Time = InitialGuess.Time(:);
                States = InitialGuess.X';
                Inputs = [InitialGuess.U'; InitialGuess.U(:, end)'];
                % FIX-F08: an infeasible plan is flown anyway (D2) and flagged
                if Report.Feasible
                    abortStatus = 0;
                else
                    abortStatus = 1;
                end
            end
        end
        if ~planOk % FIX-F08: closed-form vertical descent that cannot fail (D1)
            [Time, States, Inputs] = abortFallback(X_est, constantsTOAD, Waypoints(3).MaxVel(3));
            abortStatus = 2;
        end

        % Turn off abort to avoid recalculating and reference time from
        % when the abort was called
        abortCalculated = true; % FIX-F32: set on every path (plan, infeasible plan, fallback): one attempt only
        tAbort = t;
    end

    S.Time = Time; S.States = States; S.Inputs = Inputs;
    S.tAbort = tAbort; S.abortCalculated = abortCalculated; S.abortStatus = abortStatus;
end
