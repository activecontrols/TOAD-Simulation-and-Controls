function Waypoints = abortWaypoints(X_est, constantsTOAD, ABORT_PIN_THRUST, ABORT_HOVER_ALT)
% Verification only (never part of the flight build).
% Waypoint lines copied verbatim from the abort block in Guidance/AbortGuidance.m (shared by verify_polyqp, fix_regressions, diagnose_aborts).
% constantsTOAD is used only when ABORT_PIN_THRUST is true; ABORT_PIN_THRUST / ABORT_HOVER_ALT default to the file's constants (false, 2).
    if nargin < 2
        constantsTOAD = struct('g', NaN);
    end
    if nargin < 3
        ABORT_PIN_THRUST = false;
    end
    if nargin < 4
        ABORT_HOVER_ALT = 2;
    end
        dummyWP   = TOADWaypoint([0; 0; 0]);
        Waypoints = repmat(dummyWP, 1, 3);

        % FIX-F03: actual thrust axis b3 = C_IB * e3 (quatRot returns C_BI) and db3/dt = w_I x b3 with w_I = C_IB * w_B
        q_est = X_est(1:4) / norm(X_est(1:4));
        C_IB = quatRot(q_est)';
        b3_est = C_IB * [0; 0; 1];
        b3dot_est = cross(C_IB * X_est(11:13), b3_est);
        Waypoints(1) = TOADWaypoint(X_est(5:7), 'VelBC', X_est(8:10), 'Label', 'ABT START', 'B3BC', b3_est, 'B3DotBC', b3dot_est); % FIX-F03: knot-1 attitude
        if ABORT_PIN_THRUST % FIX-F02: a0 = (T_hover / m) * b3 - g * e3 and j0 = (T_hover / m) * db3/dt, T_hover / m = g
            Waypoints(1) = TOADWaypoint(X_est(5:7), 'VelBC', X_est(8:10), 'Label', 'ABT START', 'B3BC', b3_est, 'B3DotBC', b3dot_est, ...
                'AccelBC', constantsTOAD.g * b3_est - [0; 0; constantsTOAD.g], 'JerkBC', constantsTOAD.g * b3dot_est);
        end
        Waypoints(2) = TOADWaypoint([nan; nan; ABORT_HOVER_ALT], 'Label', 'ABT GND', 'HoldTime', 0.0, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]);
        Waypoints(3) = TOADWaypoint([nan; nan; 0], 'Label', 'ABT END', 'HoldTime', 0.0, 'MaxVel', [3; 3; 1], 'StraightLine', true, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]);
end
