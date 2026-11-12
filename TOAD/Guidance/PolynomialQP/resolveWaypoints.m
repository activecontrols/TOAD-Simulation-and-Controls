function Plan = resolveWaypoints(Waypoints, r, Limits)
%% Validate Waypoints, Expand Holds & Cruise Profiles, and Partition Legs
    MAX_WP  = 16;
    MAX_SEG = 40;
    MAX_KNOTS = MAX_SEG + 1;
    MAX_R   = 6;

    nWP = numel(Waypoints);
    if nWP < 2 || nWP > MAX_WP
        error('PolyQP:InvalidWaypointCount', ...
            'Waypoint count must be between 2 and %d (got %d).', MAX_WP, nWP);
    end
    if r < 3 || r > MAX_R
        error('PolyQP:InvalidOrder', 'Order r must be in 3..%d (got %d).', MAX_R, r);
    end

    % Validate heading reference consistency and rest knots
    hdg0 = Waypoints(1).HDGRef(:);
    if norm(hdg0) < 1e-9
        error('PolyQP:InvalidHDGRef', 'First waypoint HDGRef must be non-zero.');
    end
    hdg0 = hdg0 / norm(hdg0);

    for k = 1:nWP
        wp = Waypoints(k);
        if norm(wp.HDGRef(:) / max(1e-9, norm(wp.HDGRef(:))) - hdg0) > 1e-6
            error('PolyQP:InconsistentHDGRef', ...
                'Waypoint %d HDGRef must match first waypoint HDGRef.', k);
        end
        if any(wp.MaxVel <= 0) || any(wp.MaxAccel <= 0)
            error('PolyQP:InvalidLimits', ...
                'Waypoint %d MaxVel and MaxAccel must be strictly positive.', k);
        end
    end

    % Peak velocity factor c_v = (2r-1)! / ((r-1)!^2 * 4^(r-1)) for degree 2r-1 rest-to-rest profile
    if r == 3
        c_v = 15 / 8;
    elseif r == 4
        c_v = 35 / 16;
    elseif r == 5
        c_v = 315 / 128;
    else
        c_v = 693 / 256;
    end

    % Fixed-capacity knot and segment buffers
    FixCg       = nan(3, MAX_R, MAX_KNOTS);
    KnotIsRest  = false(1, MAX_KNOTS);
    SegT        = zeros(1, MAX_SEG);
    IsHold      = false(1, MAX_SEG);
    IsCruise    = false(1, MAX_SEG);
    SegMaxVel   = inf(3, MAX_SEG);
    SegMaxAccel = inf(3, MAX_SEG);

    % Helper to populate a knot's BCs from a waypoint
    nKnots = 1;
    nSegs  = 0;
    [FixCg(:, 1:r, 1), KnotIsRest(1)] = wpToFix(Waypoints(1), r);

    % If first waypoint has a positive HoldTime, insert initial hold segment
    if ~Waypoints(1).IsPassAndGo && Waypoints(1).HoldTime > 0
        nSegs = nSegs + 1;
        nKnots = nKnots + 1;
        FixCg(:, 1:r, nKnots) = FixCg(:, 1:r, nKnots - 1);
        KnotIsRest(nKnots) = true;
        SegT(nSegs) = Waypoints(1).HoldTime;
        IsHold(nSegs) = true;
        SegMaxVel(:, nSegs) = Waypoints(1).MaxVel;
        SegMaxAccel(:, nSegs) = min(Waypoints(1).MaxAccel, Limits.DefaultMaxAccel);
    end

    % Expand each waypoint transition m = 1 .. nWP-1
    for m = 1:(nWP - 1)
        wpB = Waypoints(m + 1);
        v_max = wpB.MaxVel;
        a_max = min(wpB.MaxAccel, Limits.DefaultMaxAccel);

        pA = FixCg(:, 1, nKnots);
        pB = wpB.PosBC;
        dp = pB - pA;

        % Check if StraightLine is requested on segment m -> m+1
        if wpB.StraightLine && all(~isnan(dp)) && norm(dp) > 1e-6
            FixCg(:, 1:r, nKnots) = pinTransverse(FixCg(:, 1:r, nKnots), dp, r);
        end

        % Check if segment is long enough to saturate MaxVel and split into [Ramp-Up, Cruise, Ramp-Down]
        didSplit = false;
        if all(~isnan(dp)) && norm(dp) > 1e-3
            T_ideal = max(abs(dp) ./ v_max);
            if T_ideal > 2 * Limits.MinSegTime
                v_cruise = dp / T_ideal;
                vA = FixCg(:, 2, nKnots);
                vA(isnan(vA)) = 0;
                vB = wpB.VelBC;
                vB(isnan(vB)) = 0;

                % Bound vertical ramp acceleration by vehicle thrust band [T_lo, T_hi]
                a_ramp = a_max;
                if Limits.HasVehicle
                    a_up_lim = max(0.30, 0.80 * (Limits.T_hi / Limits.m_init - Limits.g));
                    a_dn_lim = max(0.30, 0.80 * (Limits.g - Limits.T_lo / Limits.m_dry));
                    a_ramp(3) = min([a_ramp(3), a_up_lim, a_dn_lim]);
                end

                % Check whether endpoint velocities are compatible with a zero-accel Cruise segment:
                % Only split into [Ramp-Up, Cruise, Ramp-Down] when endpoints are either rest knots
                % (HoldTime >= 0) or have an explicitly specified VelBC aligned with dp.
                % Free pass-through knots (HoldTime = NaN, VelBC = NaN) stay in a coupled multi-segment leg
                % so their velocity & curvature are solved smoothly by the QP across the corner.
                hasExplicitVelA = KnotIsRest(nKnots) || all(~isnan(FixCg(:, 2, nKnots)));
                hasExplicitVelB = ~wpB.IsPassAndGo || all(~isnan(wpB.VelBC));
                u_dp = dp / norm(dp);
                vA_trans = norm(vA - dot(vA, u_dp) * u_dp);
                vB_trans = norm(vB - dot(vB, u_dp) * u_dp);

                if hasExplicitVelA && hasExplicitVelB && vA_trans < 1e-3 && vB_trans < 1e-3
                    T_up = max(Limits.MinSegTime, c_v * max(abs(v_cruise - vA) ./ a_ramp));
                    T_dn = max(Limits.MinSegTime, c_v * max(abs(vB - v_cruise) ./ a_ramp));

                    dp_up = 0.5 * (vA + v_cruise) * T_up;
                    dp_dn = 0.5 * (v_cruise + vB) * T_dn;
                    dp_rem = dp - dp_up - dp_dn;

                    % Require cruise segment to cover at least 25% of segment distance and >= 2*MinSegTime
                    proj_rem = dot(dp_rem, dp) / dot(dp, dp);
                    T_cruise = proj_rem * T_ideal;
                    if proj_rem >= 0.25 && T_cruise >= 2 * Limits.MinSegTime && norm(dp_rem - proj_rem * dp) < 1e-3 * norm(dp)
                        didSplit = true;
                        p_in  = pA + dp_up;
                        p_out = pB - dp_dn;

                        % 1. Ramp-Up segment (nKnots -> nKnots+1)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        FixCg(:, 1:r, nKnots) = makeCruiseKnot(p_in, v_cruise, r);
                        if wpB.StraightLine
                            FixCg(:, 1:r, nKnots - 1) = pinTransverse(FixCg(:, 1:r, nKnots - 1), dp, r);
                        end
                        KnotIsRest(nKnots) = true;
                        SegT(nSegs) = T_up;
                        SegMaxVel(:, nSegs) = v_max;
                        SegMaxAccel(:, nSegs) = a_max;

                        % 2. Constant-Velocity Cruise segment (nKnots -> nKnots+1)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        FixCg(:, 1:r, nKnots) = makeCruiseKnot(p_out, v_cruise, r);
                        KnotIsRest(nKnots) = true;
                        SegT(nSegs) = T_cruise;
                        IsCruise(nSegs) = true;
                        SegMaxVel(:, nSegs) = v_max;
                        SegMaxAccel(:, nSegs) = a_max;

                        % 3. Ramp-Down segment (nKnots -> nKnots+1 at wpB)
                        nSegs = nSegs + 1;
                        nKnots = nKnots + 1;
                        [FixCg(:, 1:r, nKnots), KnotIsRest(nKnots)] = wpToFix(wpB, r);
                        if wpB.StraightLine
                            FixCg(:, 1:r, nKnots) = pinTransverse(FixCg(:, 1:r, nKnots), dp, r);
                        end
                        SegT(nSegs) = T_dn;
                        SegMaxVel(:, nSegs) = v_max;
                        SegMaxAccel(:, nSegs) = a_max;
                    end
                end
            end
        end

        % Standard single polynomial flight segment if not split
        if ~didSplit
            vA_prev = FixCg(:, 2, nKnots);
            nSegs = nSegs + 1;
            nKnots = nKnots + 1;
            [FixCg(:, 1:r, nKnots), KnotIsRest(nKnots)] = wpToFix(wpB, r);
            if wpB.StraightLine && all(~isnan(dp)) && norm(dp) > 1e-6
                FixCg(:, 1:r, nKnots) = pinTransverse(FixCg(:, 1:r, nKnots), dp, r);
            end
            a_seg_init = a_max;
            if Limits.HasVehicle
                a_up_lim = max(0.30, 0.80 * (Limits.T_hi / Limits.m_init - Limits.g));
                a_dn_lim = max(0.30, 0.80 * (Limits.g - Limits.T_lo / Limits.m_dry));
                a_seg_init(3) = min([a_seg_init(3), a_up_lim, a_dn_lim]);
            end
            SegT(nSegs) = initSingleSegTime(dp, vA_prev, FixCg(:, 2, nKnots), v_max, a_seg_init, r, Limits.MinSegTime);
            SegMaxVel(:, nSegs) = v_max;
            SegMaxAccel(:, nSegs) = a_max;
        end

        % Append Hold segment if wpB has HoldTime > 0
        if ~wpB.IsPassAndGo && wpB.HoldTime > 0
            nSegs = nSegs + 1;
            nKnots = nKnots + 1;
            FixCg(:, 1:r, nKnots) = FixCg(:, 1:r, nKnots - 1);
            KnotIsRest(nKnots) = true;
            SegT(nSegs) = wpB.HoldTime;
            IsHold(nSegs) = true;
            SegMaxVel(:, nSegs) = v_max;
            SegMaxAccel(:, nSegs) = a_max;
        end
    end

    % Partition non-hold, non-cruise flight segments into decoupled legs
    % A knot decouples the QP if all 3xR derivatives are fixed
    KnotFullyPinned = false(1, nKnots);
    for k = 1:nKnots
        block = FixCg(:, 1:r, k);
        KnotFullyPinned(k) = all(~isnan(block(:)));
    end

    LegStartSeg = zeros(1, MAX_SEG);
    LegEndSeg   = zeros(1, MAX_SEG);
    nLegs = 0;
    m = 1;
    while m <= nSegs
        if IsHold(m) || IsCruise(m)
            m = m + 1;
            continue;
        end
        m_start = m;
        while m < nSegs && ~IsHold(m + 1) && ~IsCruise(m + 1) && ~KnotFullyPinned(m + 1)
            m = m + 1;
        end
        nLegs = nLegs + 1;
        LegStartSeg(nLegs) = m_start;
        LegEndSeg(nLegs)   = m;
        m = m + 1;
    end

    Plan = struct( ...
        'R',           r, ...
        'KnotCount',   nKnots, ...
        'SegCount',    nSegs, ...
        'FixCg',       FixCg, ...
        'KnotIsRest',  KnotIsRest, ...
        'SegT',        SegT, ...
        'IsHold',      IsHold, ...
        'IsCruise',    IsCruise, ...
        'SegMaxVel',   SegMaxVel, ...
        'SegMaxAccel', SegMaxAccel, ...
        'LegCount',    nLegs, ...
        'LegStartSeg', LegStartSeg, ...
        'LegEndSeg',   LegEndSeg, ...
        'HDGRef',      hdg0 ...
    );
end

function [slab, isRest] = wpToFix(wp, r)
%% Convert TOADWaypoint to 3 x r boundary condition slab
    slab = nan(3, r);
    slab(:, 1) = wp.PosBC;
    slab(:, 2) = wp.VelBC;
    slab(:, 3) = wp.AccelBC;
    if r >= 4
        slab(:, 4) = wp.JerkBC;
    end
    isRest = ~wp.IsPassAndGo;
    if isRest
        for k = 2:r
            for d = 1:3
                if isnan(slab(d, k))
                    slab(d, k) = 0;
                end
            end
        end
    end
end

function slab = makeCruiseKnot(pos, v_cruise, r)
%% Build fully-pinned constant-velocity transition knot
    slab = zeros(3, r);
    slab(:, 1) = pos;
    slab(:, 2) = v_cruise;
end

function slab = pinTransverse(slab, dp, r)
%% Pin derivatives >= 1 to zero on axes transverse to straight-line displacement
    for d = 1:3
        if abs(dp(d)) < 1e-6
            for k = 2:r
                if isnan(slab(d, k))
                    slab(d, k) = 0;
                end
            end
        end
    end
end

function T = initSingleSegTime(dp, vA, vB, v_max, a_max, r, MinSegTime)
%% Kinematic duration initializer for a single polynomial segment
    if r == 3
        c_v = 15 / 8;    c_a = 5.7735;
    elseif r == 4
        c_v = 35 / 16;   c_a = 7.5131;
    elseif r == 5
        c_v = 315 / 128; c_a = 9.3720;
    else
        c_v = 693 / 256; c_a = 11.2;
    end

    dp_abs = abs(dp);
    dp_abs(isnan(dp_abs)) = 0;
    dv_abs = max(abs(vA), abs(vB));
    dv_abs(isnan(dv_abs)) = 0;

    T_v  = max(c_v * (dp_abs ./ v_max));
    T_a  = max(sqrt(c_a * (dp_abs ./ a_max)));
    T_dv = max(c_v * (dv_abs ./ a_max));
    T = max([MinSegTime, T_v, T_a, T_dv]);
end
