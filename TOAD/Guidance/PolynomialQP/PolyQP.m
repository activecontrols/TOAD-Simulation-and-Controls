function [InitialGuess, QPTraj, Report] = PolyQP(Waypoints, Constants, Opts)
%% Waypoint-Driven Closed-Form Polynomial Trajectory Generator
    if nargin < 2
        Constants = [];
    end
    if nargin < 3 || isempty(Opts)
        Opts = struct('R', 5);
    end
    if ~isfield(Opts, 'R') || isempty(Opts.R)
        r = 5;
    else
        r = Opts.R;
    end

    % 1. Build vehicle limits and resolve waypoints
    Limits = buildVehicleLimits(Constants);
    Plan = resolveWaypoints(Waypoints, r, Limits);
    [A0, ~, H0] = assembleQP(r);

    % 2. Convert CG position BCs to Centre-of-Percussion P-space and solve all legs
    % Run Richter time redistribution only on legs where intermediate knots have no fixed VelBC
    [Plan, QPTraj] = solveAllLegs(Plan, A0, H0, Limits, true);

    % 3. If no vehicle constants provided, return trajectory only
    if ~Limits.HasVehicle
        InitialGuess = struct();
        Report = buildReport(true, zeros(9, 1), 0, QPTraj);
        return;
    end

    % 4. Unified Evaluate -> Dilation Loop
    MaxOuter = 40;
    stallCount = 0;
    outerIter = 0;

    [EvalOut, Feasible] = evaluateConstraints(QPTraj, Plan, Constants, Limits);
    bestUtil = EvalOut.MaxOverallUtil;
    bestSegT = Plan.SegT;
    prevSegUtil = max(EvalOut.SegUtil, [], 1);

    while ~Feasible && outerIter < MaxOuter
        outerIter = outerIter + 1;

        % Monotone per-segment dilation: T_m <- T_m * 1.02 * max_c (u_c)^(1/k_c)
        didStretch = false;
        M = Plan.SegCount;
        for m = 1:M
            if ~Plan.IsHold(m) && ~Plan.IsCruise(m) && EvalOut.SegDilate(m) > 1.0
                s_apply = min(2.0, 1.02 * EvalOut.SegDilate(m));
                Plan.SegT(m) = Plan.SegT(m) * s_apply;
                didStretch = true;
            end
        end

        if ~didStretch
            break;
        end

        % Re-solve legs with updated segment durations (no time redistribution inside loop)
        [Plan, QPTraj] = solveAllLegs(Plan, A0, H0, Limits, false);
        [EvalOut, Feasible] = evaluateConstraints(QPTraj, Plan, Constants, Limits);

        % Stall detection: continue as long as any active segment's max utilisation improved
        curSegUtil = max(EvalOut.SegUtil, [], 1);
        violMask = curSegUtil > 1.0 + 1e-3;
        madeProgress = (EvalOut.MaxOverallUtil < bestUtil - 1e-4) || ...
                       any(curSegUtil(violMask) < prevSegUtil(violMask) - 1e-4);
        prevSegUtil = curSegUtil;

        if EvalOut.MaxOverallUtil < bestUtil - 1e-4
            bestUtil = EvalOut.MaxOverallUtil;
            bestSegT = Plan.SegT;
        end

        if madeProgress
            stallCount = 0;
        else
            stallCount = stallCount + 1;
            if stallCount >= 4
                Plan.SegT = bestSegT;
                [Plan, QPTraj] = solveAllLegs(Plan, A0, H0, Limits, false);
                [EvalOut, Feasible] = evaluateConstraints(QPTraj, Plan, Constants, Limits);
                break;
            end
        end
    end

    % 5. Sample uniform output grid and build Report
    InitialGuess = sampleOutputs(QPTraj, Plan, Constants, Limits);
    Report = buildReport(Feasible, EvalOut.MaxUtilPerType, outerIter, QPTraj);
end

function [Plan, QPTraj] = solveAllLegs(Plan, A0, H0, Limits, doAllocateTimes)
%% Solve All Flight Legs with P-Space Fixed-Point Refinement
    MAX_SEG = 40;
    MAX_KNOTS = MAX_SEG + 1;
    r = Plan.R;
    nc = 2 * r;
    M = Plan.SegCount;
    K = Plan.KnotCount;
    rho = Limits.Rho;
    g = Limits.g;

    % Initialize P-space boundary conditions: P = r_cg + rho * b3 (init b3 = e3)
    FixP = Plan.FixCg(:, 1:r, 1:K);
    for k = 1:K
        for d = 1:3
            if ~isnan(FixP(d, 1, k))
                if d == 3
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho;
                else
                    FixP(d, 1, k) = Plan.FixCg(d, 1, k);
                end
            end
        end
    end

    CoefAll  = zeros(12, MAX_SEG, 3);
    DerivAll = FixP;
    totalCost = 0;

    % Optional one-time Richter time redistribution on multi-segment legs where intermediate knots have no pinned velocities
    if doAllocateTimes
        for leg = 1:Plan.LegCount
            m1 = Plan.LegStartSeg(leg);
            m2 = Plan.LegEndSeg(leg);
            if m2 > m1
                midVelBins = FixP(:, 2, (m1 + 1):m2);
                if ~any(~isnan(midVelBins(:)))
                    FixLeg = FixP(:, 1:r, m1:(m2 + 1));
                    T_leg = Plan.SegT(m1:m2);
                    [T_opt, ~] = allocateTimes(FixLeg, T_leg, A0, H0, Limits.MinSegTime);
                    Plan.SegT(m1:m2) = T_opt;
                end
            end
        end
    end

    % Fixed-point P-space BC refinement at pass-through knots (<= 3 passes)
    hasPassThrough = any(~Plan.KnotIsRest(1:K));
    nPasses = 1;
    if hasPassThrough && Limits.HasVehicle
        nPasses = 3;
    end

    for pass = 1:nPasses
        totalCost = 0;

        % Solve each decoupled flight leg
        for leg = 1:Plan.LegCount
            m1 = Plan.LegStartSeg(leg);
            m2 = Plan.LegEndSeg(leg);
            FixLeg = FixP(:, 1:r, m1:(m2 + 1));
            T_leg  = Plan.SegT(m1:m2);

            [CoefLeg, DerivLeg, costLeg] = solveLeg(FixLeg, T_leg, A0, H0);
            CoefAll(1:nc, m1:m2, :) = CoefLeg;
            DerivAll(:, 1:r, m1:(m2 + 1)) = DerivLeg;
            totalCost = totalCost + costLeg;
        end

        % Populate Hold and Cruise segment polynomial coefficients
        for m = 1:M
            if Plan.IsHold(m)
                for d = 1:3
                    CoefAll(1:nc, m, d) = 0;
                    CoefAll(1, m, d) = FixP(d, 1, m);
                end
            elseif Plan.IsCruise(m)
                for d = 1:3
                    CoefAll(1:nc, m, d) = 0;
                    CoefAll(1, m, d) = FixP(d, 1, m);
                    CoefAll(2, m, d) = FixP(d, 2, m);
                end
            end
        end

        % Update P = r_cg + rho * b3 at pass-through knots using solved acceleration
        if pass < nPasses
            for k = 2:(K - 1)
                if ~Plan.KnotIsRest(k)
                    a_k = DerivAll(:, 3, k);
                    e_k = a_k + [0; 0; g];
                    ne  = norm(e_k);
                    if ne > 1e-6
                        b3_k = e_k / ne;
                        for d = 1:3
                            if ~isnan(Plan.FixCg(d, 1, k))
                                FixP(d, 1, k) = Plan.FixCg(d, 1, k) + rho * b3_k(d);
                            end
                        end
                    end
                end
            end
        end
    end

    KeyTimes = zeros(1, MAX_KNOTS);
    KeyTimes(1:(M + 1)) = [0, cumsum(Plan.SegT(1:M))];

    FixOut = nan(3, 6, MAX_KNOTS);
    FixOut(:, 1:r, 1:K) = DerivAll;

    QPTraj = struct( ...
        'R',          r, ...
        'Count',      M, ...
        'T',          Plan.SegT, ...
        'KeyTimes',   KeyTimes, ...
        'Coef',       CoefAll(1:nc, :, :), ...
        'IsHold',     Plan.IsHold, ...
        'IsCruise',   Plan.IsCruise, ...
        'Fix',        FixOut, ...
        'Cost',       totalCost, ...
        'T_total',     sum(Plan.SegT(1:M)), ...
        'FlatOffset', rho ...
    );
end

function Report = buildReport(Feasible, MaxUtilPerType, outerIter, QPTraj)
%% Build Diagnostic Output Report
    M = QPTraj.Count;
    Report = struct( ...
        'Feasible',        Feasible, ...
        'MaxOverallUtil',  max(MaxUtilPerType), ...
        'UtilVel',         MaxUtilPerType(1), ...
        'UtilAccel',       MaxUtilPerType(2), ...
        'UtilThrustHi',    MaxUtilPerType(3), ...
        'UtilThrustLo',    MaxUtilPerType(4), ...
        'UtilGimbal',      MaxUtilPerType(5), ...
        'UtilGimbalRate',  MaxUtilPerType(6), ...
        'UtilThrustRate',  MaxUtilPerType(7), ...
        'UtilRoll',        MaxUtilPerType(8), ...
        'UtilTilt',        MaxUtilPerType(9), ...
        'OuterIterations', outerIter, ...
        'Cost',            QPTraj.Cost, ...
        'SegTimes',        QPTraj.T(1:M), ...
        'T_total',         QPTraj.T_total ...
    );
end
