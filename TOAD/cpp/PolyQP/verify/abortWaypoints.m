function Waypoints = abortWaypoints(X_est, constantsTOAD, u_cmd, a_dist, ABORT_HOVER_ALT)
% Verification only (never part of the flight build).
% Waypoint lines copied verbatim from the abort block in Guidance/AbortGuidance.m (shared by verify_polyqp, fix_regressions, diagnose_aborts).
% constantsTOAD defaults to LoadTOADParams(0); u_cmd / a_dist / ABORT_HOVER_ALT default as in the abort block (hover thrust with
% the gimbal centred, no disturbance, 2 m).
    if nargin < 2 || isempty(constantsTOAD)
        constantsTOAD = LoadTOADParams(0);
    end
    if nargin < 3 || isempty(u_cmd)
        u_cmd = [0; 0; constantsTOAD.m_dry * constantsTOAD.g; 0];
    end
    if nargin < 4 || isempty(a_dist)
        a_dist = [0; 0; 0];
    end
    if nargin < 5
        ABORT_HOVER_ALT = 2;
    end
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
end
