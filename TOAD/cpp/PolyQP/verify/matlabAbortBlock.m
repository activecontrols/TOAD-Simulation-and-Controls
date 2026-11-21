function S = matlabAbortBlock(t, X_est, constantsTOAD, ABORT, S, ABORT_HOVER_ALT, ABORT_PIN_THRUST)
% Verification only (never part of the flight build).
% The TrajectoryManager.m abort block, verbatim, with the persistents carried in S
% (verify_polyqp.m checkAbortCopyIsVerbatim asserts every line still matches the file).
% ABORT_HOVER_ALT / ABORT_PIN_THRUST default to the TrajectoryManager.m constants (2, false).
    if nargin < 6
        ABORT_HOVER_ALT = 2;
    end
    if nargin < 7
        ABORT_PIN_THRUST = false;
    end
    Time = S.Time; States = S.States; Inputs = S.Inputs;
    tAbort = S.tAbort; abortCalculated = S.abortCalculated; abortStatus = S.abortStatus;

    if ABORT > 0 && t > ABORT && ~abortCalculated
        % FIX-F29: plan from the estimated propellant
        Opts = struct('R', 5, 'InitOxMass', X_est(14), 'InitFuMass', X_est(15));      % Minimum-jounce
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
        % FIX-F08: a non-finite start state cannot be planned; it goes straight to the fallback
        planOk = all(isfinite(X_est(5:10)));
        if planOk
            if coder.target('MATLAB') % FIX-F08: code generation has no try/catch, so a PolyQP error is caught only in interpreted MATLAB
                try
                    [InitialGuess, QPTraj, Report] = PolyQP(Waypoints, constantsTOAD, Opts); %#ok<ASGLU>
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
