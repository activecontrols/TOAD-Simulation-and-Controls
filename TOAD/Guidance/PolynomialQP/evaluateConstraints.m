function [EvalOut, Feasible] = evaluateConstraints(QPTraj, Plan, Constants, Limits)
%% Dense-Grid Constraint Utilisation & Per-Segment Dilation Ratios
    M = QPTraj.Count;
    T_total = QPTraj.T_total;
    dt_eval = 0.5 * Limits.DtSample;
    nGrid = min(3999, max(41, ceil(T_total / dt_eval) + 1));
    t_grid = linspace(0, T_total, nGrid);
    dt_step = t_grid(2) - t_grid(1);

    % Invert flat trajectory and allocate controls across grid
    [X, Wrench, AccelCg] = flatnessInvert(QPTraj, t_grid, Plan.HDGRef, Constants, Limits);
    U = allocateControls(Wrench, Limits);

    % Finite-difference actuator rates on dense grid
    U_rate = zeros(4, nGrid);
    if nGrid >= 2
        U_rate(:, 1:(nGrid - 1)) = diff(U, 1, 2) / dt_step;
        U_rate(:, nGrid) = U_rate(:, nGrid - 1);
    end

    % Map grid samples to segments
    kt = QPTraj.KeyTimes(1:(M + 1));
    seg_idx = ones(1, nGrid);
    for k = 1:nGrid
        tk = t_grid(k);
        m_found = 1;
        for m = 1:M
            if tk > kt(m)
                m_found = m;
            elseif tk == kt(m) && ~QPTraj.IsHold(m) && ~Plan.IsCruise(m)
                m_found = m;
            end
        end
        seg_idx(k) = m_found;
    end

    % Zero finite-difference actuator rates at hold/cruise segment boundaries
    for k = 1:(nGrid - 1)
        m1 = seg_idx(k);
        m2 = seg_idx(k + 1);
        if QPTraj.IsHold(m1) || QPTraj.IsHold(m2) || Plan.IsCruise(m1) || Plan.IsCruise(m2)
            U_rate(:, k) = 0;
        end
    end
    U_rate(:, nGrid) = U_rate(:, max(1, nGrid - 1));

    % Utilisation arrays per segment:
    % Rows: 1=Vel, 2=Accel, 3=ThrustHi, 4=ThrustLo, 5=Gimbal, 6=GimbalRate, 7=ThrustRate, 8=Roll, 9=Tilt
    SegUtil     = zeros(9, M);
    SegDilate   = ones(1, M);

    g = Limits.g;
    G_vec = [0; 0; g];

    for k = 1:nGrid
        m = seg_idx(k);
        if QPTraj.IsHold(m)
            continue;
        end

        vk = abs(X(8:10, k));
        ak = AccelCg(:, k);
        qk = X(1:4, k);
        uk = U(:, k);
        urk = abs(U_rate(:, k));
        mk = Limits.m_dry + X(14, k) + X(15, k);

        % 1. Velocity utilisation (k = 1)
        u_vel = max(vk ./ Plan.SegMaxVel(:, m));

        % 2. Acceleration utilisation (k = 2)
        u_acc = max(abs(ak) ./ Plan.SegMaxAccel(:, m));

        % 3. Thrust upper & lower utilisation
        T_k = uk(3);
        u_thi = T_k / Limits.T_hi;
        u_tlo = Limits.T_lo / max(1e-3, T_k);

        % 4. Gimbal angle (k = 4)
        u_gim = max(abs(uk(1:2))) / Limits.GimbalLimit;

        % 5. Gimbal rate (k = 5)
        u_grate = max(urk(1:2)) / Limits.GimbalRateLim;

        % 6. Thrust rate (k = 3)
        u_trate = urk(3) / Limits.ThrustRateLim;

        % 7. Roll torque & rate (k = 4)
        u_roll = max(abs(uk(4)) / Limits.RollLimit, urk(4) / Limits.RollRateLim);

        % 8. Attitude tilt angle from vertical (k = 2)
        b3_z = 1 - 2 * (qk(2)^2 + qk(3)^2);
        tilt_k = acos(max(-1, min(1, b3_z)));
        u_tilt = tilt_k / Limits.MaxTilt;

        % Skip recording finite-difference actuator rates on exact 1D cruise segments
        if Plan.IsCruise(m)
            u_grate = 0;
            u_trate = 0;
        end

        u_col = [u_vel; u_acc; u_thi; u_tlo; u_gim; u_grate; u_trate; u_roll; u_tilt];
        SegUtil(:, m) = max(SegUtil(:, m), u_col);

        % Compute required segment stretch ratio s >= 1 if any constraint is violated
        if ~Plan.IsCruise(m)
            s_req = 1.0;
            if u_vel > 1.0,   s_req = max(s_req, u_vel);        end
            if u_acc > 1.0,   s_req = max(s_req, u_acc^(1/2));  end
            if u_gim > 1.0,   s_req = max(s_req, u_gim^(1/4));  end
            if u_grate > 1.0, s_req = max(s_req, u_grate^(1/5)); end
            if u_trate > 1.0, s_req = max(s_req, u_trate^(1/3)); end
            if u_roll > 1.0,  s_req = max(s_req, u_roll^(1/4)); end
            if u_tilt > 1.0,  s_req = max(s_req, u_tilt^(1/2)); end

            % Closed-form quadratic thrust scaling for T > T_hi: |a*w + G| = T_allow / m, w = 1/s^2
            if u_thi > 1.0
                s_thi = solveThrustStretch(ak, G_vec, Limits.T_hi / mk, u_thi);
                s_req = max(s_req, s_thi);
            end
            if u_tlo > 1.0
                s_tlo = solveThrustFloorStretch(ak, G_vec, Limits.T_lo / mk, u_tlo);
                s_req = max(s_req, s_tlo);
            end

            SegDilate(m) = max(SegDilate(m), s_req);
        end
    end

    maxUtilPerType = max(SegUtil, [], 2);
    maxOverallUtil = max(maxUtilPerType);
    Feasible = (maxOverallUtil <= 1.0 + 1e-3);

    EvalOut = struct( ...
        'SegUtil',        SegUtil, ...
        'MaxUtilPerType', maxUtilPerType, ...
        'MaxOverallUtil', maxOverallUtil, ...
        'SegDilate',      SegDilate ...
    );
end

function s = solveThrustStretch(a, G_vec, e_allow, u_thi)
%% Quadratic root for time stretch s when specific force exceeds e_allow
    g = G_vec(3);
    a2 = dot(a, a);
    aG = dot(a, G_vec);
    c  = g^2 - e_allow^2;
    if a2 < 1e-8 || c >= 0
        s = sqrt(u_thi);
        return;
    end
    disc = max(0, aG^2 - a2 * c);
    w = (-aG + sqrt(disc)) / a2;
    if w > 0 && w < 1
        s = 1 / sqrt(w);
    else
        s = sqrt(u_thi);
    end
end

function s = solveThrustFloorStretch(a, G_vec, e_floor, ~)
%% Time stretch s when downward acceleration causes thrust to drop below T_lo
    g = G_vec(3);
    a2 = dot(a, a);
    aG = dot(a, G_vec);
    c  = g^2 - e_floor^2;
    if a2 < 1e-8 || c <= 0
        % Hover thrust m*g is already below T_lo: time dilation cannot raise thrust
        s = 1.0;
        return;
    end
    disc = max(0, aG^2 - a2 * c);
    w = (-aG - sqrt(disc)) / a2;
    if w > 0 && w < 1
        s = 1 / sqrt(w);
    else
        s = 1.0;
    end
end
