function [InitialGuess, QPTraj] = PolyQP(FixArray, SegTimes, Constants, Opts)
%% POLYQP  Unified Polynomial QP & Differential Flatness Guidance Engine.
%   Closed-form minimum-derivative trajectory optimization (Richter et al., IJRR 2016)
%   with corridor enforcement, time redistribution, and 6-DoF flatness inversion.
%
%   Inputs:
%       FixArray  : D x R x K matrix of keyframe boundary conditions (NaN = free)
%       SegTimes  : 1 x (K-1) vector of segment durations [s]
%       Constants : Vehicle parameters struct from LoadTOADParams (optional)
%       Opts      : Configuration options struct
%
%   Outputs:
%       InitialGuess : Struct containing (Time, X, U, T_total, dt_row)
%       QPTraj       : Struct containing polynomial coefficients, times, and Eval handle

%   ABORT MODE: PolyQP(AbortIn, [], Constants, struct('Mode',"Abort", ...))
%       AbortIn = flight state struct (r, v, a, q, w, m, ml, mi, lz)
%       Returns InitialGuess.X (14 x N+1 = x_target) and InitialGuess.U (4 x N = U_ff)

    %% 0. Abort-mode
    if nargin >= 4 && isstruct(Opts) && isfield(Opts, 'Mode') && strcmpi(string(Opts.Mode), "Abort")
        [InitialGuess, QPTraj] = abort_guidance(FixArray, Constants, Opts); return;
    end

    %% 1. Configuration & Options Handling
    if nargin < 3, Constants = []; end
    if nargin < 4 || isempty(Opts), Opts = struct(); end

    def = struct('Tag', {repmat({''}, 1, size(FixArray, 3))}, 'MinSegTime', 0.25, ...
        'OptimizeTimes', false, 'MaxIter', 10, 'Corridors', {{}}, 'MaxInsert', 20, ...
        'DtSample', 0.05, 'N', 60, 'Vehicle', 0, 'Maneuver', "Circle", ...
        'ManeuverParams', struct(), 'q0', [1; 0; 0; 0], 'w0', [0; 0; 0], ...
        'r0', [0; 0; 0], 'rf', [0; 0; 0], 'ThrustMargin', 0.05, 'GimbalMargin', 0.15, ...
        'MaxGimbalAngle', pi/15, 'MaxGimbalRate', deg2rad(35), 'MaxThrustRate', 1500, 'MaxRollRate', 5);
    f = fieldnames(def);
    for i = 1:numel(f), if ~isfield(Opts, f{i}), Opts.(f{i}) = def.(f{i}); end, end

    %% 2. Solve Closed-Form Polynomial QP (Richter Eq. 26)
    QPTraj = poly_solve(FixArray, SegTimes, Opts.Tag, Opts.MinSegTime);

    %% 3. Segment Duration & Corridor Optimization (Richter Sec. V & VI-B)
    if Opts.OptimizeTimes && numel(QPTraj.T) >= 2
        M_t = numel(QPTraj.T); Ttot = sum(QPTraj.T);
        fcost = @(a_log) time_cost(QPTraj.Fix, a_log, Ttot, QPTraj.MinSegTime);
        a_log = log(QPTraj.T(:)); J0 = fcost(a_log);

        for it = 1:Opts.MaxIter
            g = zeros(M_t, 1); h_step = 1e-3;
            for i = 1:M_t
                ap = a_log; ap(i) = ap(i) + h_step; g(i) = (fcost(ap) - J0) / h_step;
            end
            g = g - mean(g);
            if max(abs(g)) < 1e-9 * max(J0, 1e-12), break; end
            gn = g / max(abs(g));
            improved = false;
            for step = [0.4, 0.2, 0.1, 0.05, 0.025]
                an = a_log - step * gn;
                Jn = fcost(an);
                if Jn < J0 * (1 - 1e-6)
                    a_log = an; J0 = Jn; improved = true; break;
                end
            end
            if ~improved, break; end
        end
        w = exp(a_log - max(a_log));
        Tn = max(Ttot * w(:).' / sum(w), QPTraj.MinSegTime);
        QPTraj = poly_solve(QPTraj.Fix, Ttot * Tn / sum(Tn), QPTraj.Tag, QPTraj.MinSegTime);
    end

    if ~isempty(Opts.Corridors)
        for it = 1:Opts.MaxInsert
            Ttot = sum(QPTraj.T);
            ts = linspace(0, Ttot, max(100, ceil(Ttot / Opts.DtSample) + 1));
            rs = poly_eval(QPTraj, ts, 0);
            gap = min(abs(QPTraj.KeyTimes(:) - ts), [], 1);
            ins = zeros(1 + size(rs, 1), 0);

            for c = 1:numel(Opts.Corridors)
                [viol, tgt] = Opts.Corridors{c}.fn(ts, rs, QPTraj.TagTimes);
                viol(gap < QPTraj.MinSegTime) = -Inf;
                [vm, i] = max(viol);
                if vm > 0, ins = [ins, [ts(i); tgt(:, i)]]; end %#ok<AGROW>
            end
            if isempty(ins), break; end

            [~, ord] = sort(ins(1, :));
            ins = ins(:, ord);
            ins = ins(:, [true, diff(ins(1, :)) >= QPTraj.MinSegTime]);

            did_insert = false;
            kt = QPTraj.KeyTimes;
            for j = 1:size(ins, 2)
                t_ins = ins(1, j); pos = ins(2:end, j);
                m = find(t_ins > kt(1:end-1) & t_ins < kt(2:end), 1);
                if isempty(m) || (t_ins - kt(m) < QPTraj.MinSegTime) || (kt(m+1) - t_ins < QPTraj.MinSegTime)
                    continue;
                end
                slab = nan(size(QPTraj.Fix, 1), QPTraj.R, 1); slab(:, 1, 1) = pos(:);
                NewFix = cat(3, QPTraj.Fix(:, :, 1:m), slab, QPTraj.Fix(:, :, m+1:end));
                NewTag = [QPTraj.Tag(1:m), {''}, QPTraj.Tag(m+1:end)];
                NewT   = [QPTraj.T(1:m-1), t_ins - kt(m), kt(m+1) - t_ins, QPTraj.T(m+1:end)];
                QPTraj = poly_solve(NewFix, NewT, NewTag, QPTraj.MinSegTime);
                kt = QPTraj.KeyTimes; did_insert = true;
            end
            if ~did_insert, break; end
        end
    end

    %% 4. Differential Flatness Inversion (6-DoF State & Actuator Synthesis)
    if isempty(Constants)
        InitialGuess = struct();
        return;
    end

    g = Constants.g; MT = Constants.MaxThrust;
    thr_allow = (1 - Opts.ThrustMargin) * MT;

    % Analytical thrust sizing check: scale duration if required thrust exceeds capability
    m_chk = Constants.m_dry;
    if Opts.Vehicle == 1, m_chk = m_chk + Constants.OxMass + Constants.FuMass; end
    tv_c = linspace(0, sum(QPTraj.T), 40);
    a_c = poly_eval(QPTraj, tv_c, 2);
    kap = sqrt(max(m_chk * sqrt(sum((a_c + [0; 0; g]).^2, 1))) / thr_allow);
    if kap > 1.0
        QPTraj = poly_solve(QPTraj.Fix, QPTraj.T * kap, QPTraj.Tag, QPTraj.MinSegTime);
    end

    % Sample polynomial at collocation nodes
    Nc = Opts.N; Nn = Nc + 1;
    T_total = sum(QPTraj.T);
    tv = linspace(0, T_total, Nn); dt = T_total / Nc;
    r_s = poly_eval(QPTraj, tv, 0); v_s = poly_eval(QPTraj, tv, 1);
    a_s = poly_eval(QPTraj, tv, 2); j_s = poly_eval(QPTraj, tv, 3);

    % Mass depletion tracking
    m_dry = Constants.m_dry;
    ml = zeros(1, Nn); mi = zeros(1, Nn);
    if Opts.Vehicle == 1, ml(1) = Constants.OxMass; mi(1) = Constants.FuMass; end
    m_vec = (m_dry + ml(1) + mi(1)) * ones(1, Nn);

    q_s = zeros(4, Nn); w_s = zeros(3, Nn); U_s = zeros(4, Nc);
    plim = sin((1 - Opts.GimbalMargin) * Opts.MaxGimbalAngle);
    t_lo = (0.25 + Opts.ThrustMargin) * MT; t_hi = (1.00 - Opts.ThrustMargin) * MT;
    OF = Constants.OF; Mdot = Constants.MaxMdot;

    if Opts.Maneuver == "Backflip"
        mp = Opts.ManeuverParams;
        ks = max(3, round(mp.flip_start_frac * Nc)) + 1;
        ke = min(Nc - 2, round(mp.flip_end_frac * Nc)) + 1;
        Lf = max(2, ke - ks); Tf = tv(ke) - tv(ks);
        Jyy = Constants.J(2, 2); cgz = Constants.rTB;

        for k = 1:Nn
            if k <= ks
                q_s(:, k) = Opts.q0;
            elseif k <= ke
                tau = (k - ks) / Lf;
                th  = -2*pi * (10*tau^3 - 15*tau^4 + 6*tau^5);
                dth = -(2*pi/Tf) * (30*tau^2 - 60*tau^3 + 30*tau^4);
                q_s(:, k) = [cos(th/2); 0; sin(th/2); 0];
                w_s(:, k) = [0; dth; 0];
            else
                q_s(:, k) = -Opts.q0;
            end
        end

        for k = 1:Nc
            F_req = m_vec(k) * (a_s(:, k) + [0; 0; g]);
            U_s(3, k) = max(t_lo, min(t_hi, norm(F_req)));
            if k > ks && k <= ke
                tau = (k - ks) / Lf;
                th  = -2*pi * (10*tau^3 - 15*tau^4 + 6*tau^5);
                ddth = -(2*pi/Tf^2) * (60*tau - 180*tau^2 + 120*tau^3);
                if cos(th) < 0, U_s(3, k) = max(t_lo, U_s(3, k) * (1.0 + 0.20 * cos(th))); end
                sp = -(Jyy * ddth) / max(1.0, cgz * U_s(3, k));
                U_s(2, k) = asin(max(-plim, min(plim, sp)));
            end
            if Opts.Vehicle == 1
                ml(k+1) = max(0.02*Constants.OxMass, ml(k) - dt*(U_s(3, k)/MT)*(OF/(1+OF))*Mdot);
                mi(k+1) = max(0.02*Constants.FuMass, mi(k) - dt*(U_s(3, k)/MT)*(1/(1+OF))*Mdot);
                m_vec(k+1) = m_dry + ml(k+1) + mi(k+1);
            end
        end
    else
        for k = 1:Nc
            F_req = m_vec(k) * (a_s(:, k) + [0; 0; g]);
            U_s(3, k) = max(t_lo, min(t_hi, norm(F_req)));
            if Opts.Vehicle == 1
                ml(k+1) = max(0.02*Constants.OxMass, ml(k) - dt*(U_s(3, k)/MT)*(OF/(1+OF))*Mdot);
                mi(k+1) = max(0.02*Constants.FuMass, mi(k) - dt*(U_s(3, k)/MT)*(1/(1+OF))*Mdot);
                m_vec(k+1) = m_dry + ml(k+1) + mi(k+1);
            end
        end

        for k = 1:Nn
            F_req = m_vec(k) * (a_s(:, k) + [0; 0; g]);
            T_mag = norm(F_req);
            b3 = F_req / max(1e-3, T_mag);
            q_flat = vectorToQuat(b3);

            if k <= 4
                alpha = (k - 1) / 4;
                if dot(Opts.q0, q_flat) < 0, q_flat = -q_flat; end
                qb = (1 - alpha) * Opts.q0 + alpha * q_flat;
                q_s(:, k) = qb / norm(qb);
            elseif k >= Nn - 3
                alpha = (Nn - k) / 3;
                if dot(Opts.q0, q_flat) < 0, q_flat = -q_flat; end
                qb = (1 - alpha) * Opts.q0 + alpha * q_flat;
                q_s(:, k) = qb / norm(qb);
            else
                q_s(:, k) = q_flat;
            end

            b3_dot = (j_s(:, k) - b3 * (b3' * j_s(:, k))) * m_vec(k) / max(1e-3, T_mag);
            b3d_B = quatRot(q_s(:, k)) * b3_dot;
            w_s(:, k) = [-b3d_B(2); b3d_B(1); 0];
        end
        q_s(:, 1) = Opts.q0; q_s(:, end) = Opts.q0;
        w_s(:, 1) = Opts.w0; w_s(:, end) = [0; 0; 0];

        J_mat = Constants.J; rTB = Constants.rTB;
        for k = 1:Nc
            w_dot = (w_s(:, k+1) - w_s(:, k)) / dt;
            w_mid = 0.5 * (w_s(:, k) + w_s(:, k+1));
            tau_req = J_mat * w_dot + cross(w_mid, J_mat * w_mid);
            denom = max(1.0, U_s(3, k) * rTB);
            U_s(1, k) = asin(max(-plim, min(plim, -tau_req(2) / denom)));
            U_s(2, k) = asin(max(-plim, min(plim,  tau_req(1) / denom)));
            U_s(4, k) = 0;
        end
    end

    % Forward-backward actuator rate smoothing
    rate_limits = [Opts.MaxGimbalRate; Opts.MaxGimbalRate; Opts.MaxThrustRate; Opts.MaxRollRate];
    max_du = 0.90 * rate_limits * dt;
    for k = 2:Nc, U_s(:, k) = max(U_s(:, k-1) - max_du, min(U_s(:, k-1) + max_du, U_s(:, k))); end
    for k = Nc-1:-1:1, U_s(:, k) = max(U_s(:, k+1) - max_du, min(U_s(:, k+1) + max_du, U_s(:, k))); end

    r_s(:, 1) = Opts.r0; r_s(:, end) = Opts.rf; v_s(:, end) = [0; 0; -0.05];

    InitialGuess = struct('Time', tv, 'X', [q_s; r_s; v_s; w_s; ml; mi], ...
        'U', U_s, 'T_total', T_total, 'dt_row', dt * ones(1, Nc));
end

%% Section 1: Local Functions

function Traj = poly_solve(FixArray, SegTimes, KeyframeTags, MinSegTime)
%% POLY_SOLVE  Closed-form constrained QP solve (Richter Eq. 26).
    if nargin < 3 || isempty(KeyframeTags), KeyframeTags = repmat({''}, 1, size(FixArray, 3)); end
    if nargin < 4 || isempty(MinSegTime),   MinSegTime   = 0.25; end

    R = size(FixArray, 2); K = size(FixArray, 3); D = size(FixArray, 1);
    M = K - 1; nc = 2 * R; nv = R * K; SegTimes = SegTimes(:).';

    A0 = zeros(nc); Q0 = zeros(nc);
    for k = 0:R-1
        A0(k + 1, k + 1) = factorial(k);
        for i = k:(nc - 1), A0(R + k + 1, i + 1) = factorial(i) / factorial(i - k); end
    end
    for i = R:(nc - 1)
        for j = R:(nc - 1)
            Q0(i + 1, j + 1) = (factorial(i)/factorial(i-R)) * (factorial(j)/factorial(j-R)) / (i + j - 2*R + 1);
        end
    end
    A0inv = inv(A0);
    H0 = 0.5 * (A0inv.' * Q0 * A0inv + (A0inv.' * Q0 * A0inv).');

    H = zeros(nv);
    for m = 1:M
        s = repmat(SegTimes(m) .^ (0:R-1), 1, 2).';
        Hm = SegTimes(m)^(1 - 2*R) * (s * s.') .* H0;
        ix = (m - 1)*R + (1:nc);
        H(ix, ix) = H(ix, ix) + Hm;
    end

    Deriv = zeros(D, R, K); Cost = 0;
    for d = 1:D
        fx = reshape(FixArray(d, :, :), nv, 1);
        F = find(~isnan(fx)); P = find(isnan(fx));
        dv = fx;
        if ~isempty(P)
            HPP = 0.5 * (H(P, P) + H(P, P).');
            rc = rcond(HPP);
            if isnan(rc) || rc < 1e-12, dv(P) = -pinv(HPP) * (H(P, F) * fx(F));
            else, dv(P) = -(HPP \ (H(P, F) * fx(F))); end
        end
        Deriv(d, :, :) = reshape(dv, 1, R, K); Cost = Cost + dv.' * H * dv;
    end

    Coef = zeros(nc, M, D);
    for m = 1:M
        s = repmat(SegTimes(m) .^ (0:R-1), 1, 2).';
        for d = 1:D
            Coef(:, m, d) = A0inv * (s .* reshape(Deriv(d, :, m:m+1), nc, 1)); %#ok<MINV>
        end
    end

    kt = [0, cumsum(SegTimes)];
    tg = struct();
    for k = 1:numel(KeyframeTags)
        if ~isempty(KeyframeTags{k}), tg.(KeyframeTags{k}) = kt(k); end
    end
    Traj = struct('Fix', FixArray, 'T', SegTimes, 'Tag', {KeyframeTags}, ...
                  'MinSegTime', MinSegTime, 'R', R, 'Coef', Coef, 'Deriv', Deriv, ...
                  'Cost', Cost, 'T_total', sum(SegTimes), 'KeyTimes', kt, 'TagTimes', tg);
    Traj.Eval = @(t, deriv) poly_eval(Traj, t, deriv);
end

function Out = poly_eval(Traj, T_query, DerivOrder)
%% POLY_EVAL  Vectorized evaluation of the k-th time derivative.
    if nargin < 3, DerivOrder = 0; end
    T_query = T_query(:).';
    M = numel(Traj.T); D = size(Traj.Coef, 3); nc = size(Traj.Coef, 1);
    cumT = [0, cumsum(Traj.T)];
    idx = min(max(sum(bsxfun(@ge, T_query, cumT(1:M).'), 1), 1), M);
    tau = min(max((T_query - cumT(idx)) ./ Traj.T(idx), 0), 1);
    Out = zeros(D, numel(T_query));
    for i = DerivOrder:(nc - 1)
        f = factorial(i) / factorial(i - DerivOrder);
        Ci = reshape(Traj.Coef(i + 1, idx, :), numel(T_query), D).';
        Out = Out + f * Ci .* repmat(tau .^ (i - DerivOrder), D, 1);
    end
    Out = Out ./ repmat(Traj.T(idx) .^ DerivOrder, D, 1);
end

function J = time_cost(FixArray, LogDurations, TotalTime, MinSegTime)
%% TIME_COST  QP evaluation inside log-duration search space (Rule 5 compliant).
    w = exp(LogDurations - max(LogDurations));
    T = max(TotalTime * w(:).' / sum(w), MinSegTime);
    T_eval = poly_solve(FixArray, TotalTime * T / sum(T), {}, MinSegTime);
    J = T_eval.Cost;
end



function [IG, QPTraj] = abort_guidance(Ab, C, Opts)
%% ABORT_GUIDANCE  In-flight abort: waypoint over LZ -> min-energy touchdown trajectory.
%   Ab fields: r, v (required, 3x1 NWU); a (3x1 current accel, default 0); q (4x1 current
%   quat, default [1;0;0;0]); w (3x1, default 0); m (current mass); ml, mi (current
%   oxidizer / fuel mass, Vehicle==1); lz (3x1 landing point, default [0;0;0]).
    d = struct('N',60,'MinSegTime',0.25,'DtSample',0.05,'MaxInsert',20,'Vehicle',0, ...
        'ThrustMargin',0.05,'GimbalMargin',0.15,'MaxGimbalAngle',pi/15, ...
        'MaxGimbalRate',deg2rad(35),'MaxThrustRate',1500,'MaxRollRate',5, ...
        'AbortWpAlt',5,'TouchdownVz',0.3,'FunnelRadius',0.3,'FunnelAngle',deg2rad(15), ...
        'AbortTfRange',[3 30],'AbortNTf',14,'MaxTilt',deg2rad(20), ...
        'BlendNodes',4,'Yaw',[1;0;0],'AbortTopK',3);
    Opts = abort_fill(Opts, d);

    r0  = Ab.r(:);  
    v0 = Ab.v(:);
    a0  = abort_get(Ab,'a',[0;0;0]);
    q0  = abort_get(Ab,'q',[1;0;0;0]);  
    q0 = q0/norm(q0);
    w0  = abort_get(Ab,'w',[0;0;0]);
    r_lz = abort_get(Ab,'r_lz',[0;0;0]);
    ml0 = 0; 
    mi0 = 0;
    if Opts.Vehicle == 1
        ml0 = abort_get(Ab,'ml',C.OxMass); mi0 = abort_get(Ab,'mi',C.FuMass);
    end
    m0  = abort_get(Ab,'m', C.m_dry + ml0 + mi0);

    % --- 1. Keyframes: start / waypoint above LZ / touchdown (R=4 -> min-snap) ---
    z0 = r0(3) - r_lz(3);
    use_wp = z0 > 1.0;                          % skip waypoint if already very low. 1 meter is an arbitrary number, I don't know if this is a good height
    h_wp = min(Opts.AbortWpAlt, max(0.5*z0, 0.5));
    r_wp = [r_lz(1:2); r_lz(3) + h_wp];
    K = 2 + use_wp;
    Fix = nan(3, 4, K);
    Fix(:,1,1) = r0; Fix(:,2,1) = v0; Fix(:,3,1) = a0;        % jerk free
    if use_wp, Fix(:,1,2) = r_wp; end                         % waypoint: position only
    Fix(:,1,K) = r_lz; Fix(:,2,K) = [0; 0; -Opts.TouchdownVz]; Fix(:,3,K) = 0;
    if use_wp, Tag = {'','wp','land'}; else, Tag = {'','land'}; end

    qpo = struct(); qpo.Tag = Tag; qpo.MinSegTime = Opts.MinSegTime;
    qpo.DtSample = Opts.DtSample; qpo.MaxInsert = Opts.MaxInsert;

    % --- 2. Sweep final time Tf. pick least-energy feasible (no corridor yet, cheap) ---
    Tgrid = logspace(log10(Opts.AbortTfRange(1)), log10(Opts.AbortTfRange(2)), Opts.AbortNTf);
    Gc = inf(1, numel(Tgrid)); Jc = inf(1, numel(Tgrid));
    for i = 1:numel(Tgrid)
        T = abort_seg_times(Tgrid(i), use_wp, r0, r_wp, r_lz, Opts.MinSegTime);
        [~, Q] = PolyQP(Fix, T, [], qpo);
        [Gc(i), Jc(i)] = abort_metrics(Q, C, Opts, m0, r_lz);
    end
    key = Jc; key(Gc > 0) = Inf;
    if all(isinf(key)), [~, order] = sort(Gc); else, [~, order] = sort(key); end

    % --- 3. Re-solve best candidates with the landing-funnel corridor on the last leg ---
    cor = struct('fn', @(ts, rs, tg) abort_funnel(ts, rs, tg, r_lz, Opts));
    qpo.Corridors = {cor};
    bestScore = Inf; QPTraj = []; feas = false; Gm = NaN; En = NaN; Cv = NaN;
    for ii = order(1:min(Opts.AbortTopK, numel(order)))
        T = abort_seg_times(Tgrid(ii), use_wp, r0, r_wp, r_lz, Opts.MinSegTime);
        [~, Q] = PolyQP(Fix, T, [], qpo);
        [gm, E] = abort_metrics(Q, C, Opts, m0, r_lz);
        tsc = linspace(0, Q.T_total, 200);
        cv = max(abort_funnel(tsc, poly_eval(Q, tsc, 0), Q.TagTimes, r_lz, Opts));
        score = max(gm, cv);
        if score < bestScore
            bestScore = score; QPTraj = Q; Gm = gm; En = E; Cv = cv;
            feas = (gm <= 0) && (cv <= 1e-3);
        end
        if feas, break; end
    end

    % --- 4. Position trace -> accel -> quaternion -> rates -> U_ff ---
    IG = abort_synthesis(QPTraj, C, Opts, m0, ml0, mi0, q0, w0);
    IG.Abort = struct('Feasible', feas, 'ThrustAttMetric', Gm, 'CorridorViol', Cv, ...
                      'Energy', En, 'Tf', QPTraj.T_total, 'UseWaypoint', use_wp);
end

function T = abort_seg_times(Tf, use_wp, r0, r_wp, rlz, Tmin)
    if use_wp
        L1 = norm(r_wp - r0) + 1e-3; L2 = norm(r_wp - rlz) + 1e-3;
        T1 = max(Tmin, Tf * L1 / (L1 + L2)); T = [T1, max(Tmin, Tf - T1)];
    else
        T = Tf;
    end
end

function [gm, E] = abort_metrics(Q, C, Opts, m0, rlz)
%% Feasibility metric (<=0 means feasible) and energy (impulse, integral of |F| dt).
    MT = C.MaxThrust; t_lo = (0.25 + Opts.ThrustMargin) * MT; t_hi = (1 - Opts.ThrustMargin) * MT;
    ts = linspace(0, Q.T_total, 120);
    r = poly_eval(Q, ts, 0); a = poly_eval(Q, ts, 2);
    F = m0 * (a + [0; 0; C.g]);                 % constant m0 (conservative, see notes)
    Fn = sqrt(sum(F.^2, 1));
    tilt = acos(min(1, max(-1, F(3,:) ./ max(Fn, 1e-6))));
    gm = max([max(Fn)/t_hi - 1, (t_lo - min(Fn))/MT, max(tilt)/Opts.MaxTilt - 1, ...
              (rlz(3) - min(r(3,:))) - 0.02]);
    E = trapz(ts, Fn);
end

function [viol, tgt] = abort_funnel(ts, rs, tg, rlz, Opts)
%% Corridor for the last leg (after tag 'wp'): descent funnel + ground floor.
%   Same interface as Opts.Corridors{c}.fn in PolyQP: viol>0 means violated,
%   tgt = 3xN position to pin if a keyframe gets inserted there.
    t0 = 0; if isfield(tg, 'wp'), t0 = tg.wp; end
    dxy = rs(1:2,:) - rlz(1:2);
    rad = sqrt(sum(dxy.^2, 1));
    zh  = rs(3,:) - rlz(3);
    rmax = Opts.FunnelRadius + max(zh, 0) * tan(Opts.FunnelAngle);
    viol = max(rad - rmax, -(zh + 0.02));
    viol(ts < t0) = -Inf;
    sc  = min(1, rmax ./ max(rad, 1e-9)) * 0.9;
    tgt = [rlz(1:2) + dxy .* sc; rlz(3) + max(zh, 0.05)];
end

function IG = abort_synthesis(Q, C, Opts, m0, ml0, mi0, q_meas, w_meas)
%% Flatness inversion on the abort trajectory.
    Nc = Opts.N; Nn = Nc + 1; Tt = Q.T_total; tv = linspace(0, Tt, Nn); dt = Tt / Nc;
    r_s = poly_eval(Q, tv, 0); v_s = poly_eval(Q, tv, 1); a_s = poly_eval(Q, tv, 2);
    g = C.g; MT = C.MaxThrust; OF = C.OF; Mdot = C.MaxMdot;
    t_lo = (0.25 + Opts.ThrustMargin) * MT; t_hi = (1 - Opts.ThrustMargin) * MT;
    plim = sin((1 - Opts.GimbalMargin) * Opts.MaxGimbalAngle);

    ml = zeros(1, Nn); mi = zeros(1, Nn); ml(1) = ml0; mi(1) = mi0;
    m_vec = m0 * ones(1, Nn); U = zeros(4, Nc);
    for k = 1:Nc                                             % thrust + mass (as in PolyQP)
        U(3,k) = max(t_lo, min(t_hi, norm(m_vec(k) * (a_s(:,k) + [0;0;g]))));
        if Opts.Vehicle == 1
            ml(k+1) = max(0.02*C.OxMass, ml(k) - dt*(U(3,k)/MT)*(OF/(1+OF))*Mdot);
            mi(k+1) = max(0.02*C.FuMass, mi(k) - dt*(U(3,k)/MT)*(1/(1+OF))*Mdot);
            m_vec(k+1) = C.m_dry + ml(k+1) + mi(k+1);
        end
    end

    q_s = zeros(4, Nn); nb = max(1, Opts.BlendNodes);        % attitude from accel vector
    for k = 1:Nn
        qf = abort_acc2quat(m_vec(k) * (a_s(:,k) + [0;0;g]), Opts.Yaw);
        if k <= nb + 1
            if dot(qf, q_meas) < 0, qf = -qf; end
            al = (k - 1) / nb; qb = (1 - al) * q_meas + al * qf; qf = qb / norm(qb);
        end
        if k > 1 && dot(qf, q_s(:,k-1)) < 0, qf = -qf; end
        q_s(:,k) = qf;
    end
    w_s = abort_quat_rates(q_s, tv);                         % body rates from q trace
    w_s(:,1) = w_meas; w_s(:,end) = [0;0;0];

    J = C.J; rTB = C.rTB;                                    % torque -> gimbal (same map as PolyQP)
    for k = 1:Nc
        w_dot = (w_s(:,k+1) - w_s(:,k)) / dt; w_mid = 0.5 * (w_s(:,k) + w_s(:,k+1));
        tau = J * w_dot + cross(w_mid, J * w_mid);
        den = max(1.0, U(3,k) * rTB);
        U(1,k) = asin(max(-plim, min(plim, -tau(2) / den)));
        U(2,k) = asin(max(-plim, min(plim,  tau(1) / den)));
        U(4,k) = 0;
    end

    lim = [Opts.MaxGimbalRate; Opts.MaxGimbalRate; Opts.MaxThrustRate; Opts.MaxRollRate];
    md = 0.90 * lim * dt;                                    % forward-backward rate smoothing
    for k = 2:Nc, U(:,k) = max(U(:,k-1) - md, min(U(:,k-1) + md, U(:,k))); end
    for k = Nc-1:-1:1, U(:,k) = max(U(:,k+1) - md, min(U(:,k+1) + md, U(:,k))); end

    IG = struct('Time', tv, 'X', [q_s; r_s; v_s; w_s; ml; mi], 'U', U, ...
                'T_total', Tt, 'dt_row', dt * ones(1, Nc));
end

function q = abort_acc2quat(F, psi)
%% Thrust vector -> body-to-earth quaternion (scalar first). Zb = F/|F|, heading from psi.
%   Same construction as ASTRA doc Eq 4.10-4.13, but Yb is normalized and F includes gravity.
    n = norm(F); if n < 1e-6, F = [0;0;1]; n = 1; end
    zb = F / n; yb = cross(zb, psi);
    if norm(yb) < 1e-6, yb = cross(zb, [0;1;0]); end
    yb = yb / norm(yb); xb = cross(yb, zb);
    q = abort_rot2quat([xb, yb, zb]);
end

function q = abort_rot2quat(R)
%% Shepperd's method, all four branches (redundant with the controller's version).
    tr = R(1,1) + R(2,2) + R(3,3);
    [~, i] = max([tr, R(1,1), R(2,2), R(3,3)]);
    switch i
        case 1, s = 2*sqrt(1+tr);
            q = [s/4; (R(3,2)-R(2,3))/s; (R(1,3)-R(3,1))/s; (R(2,1)-R(1,2))/s];
        case 2, s = 2*sqrt(1+R(1,1)-R(2,2)-R(3,3));
            q = [(R(3,2)-R(2,3))/s; s/4; (R(1,2)+R(2,1))/s; (R(1,3)+R(3,1))/s];
        case 3, s = 2*sqrt(1+R(2,2)-R(1,1)-R(3,3));
            q = [(R(1,3)-R(3,1))/s; (R(1,2)+R(2,1))/s; s/4; (R(2,3)+R(3,2))/s];
        otherwise, s = 2*sqrt(1+R(3,3)-R(1,1)-R(2,2));
            q = [(R(2,1)-R(1,2))/s; (R(1,3)+R(3,1))/s; (R(2,3)+R(3,2))/s; s/4];
    end
    q = q / norm(q); if q(1) < 0, q = -q; end
end

function W = abort_quat_rates(Q, tv)
%% Body rates from q trace: [0;w] = 2 * conj(q) (x) qdot  (matches qdot = 1/2 M(q)[0;w]).
    n = size(Q, 2); W = zeros(3, n);
    for k = 1:n
        k1 = max(k-1, 1); k2 = min(k+1, n);
        qd = (Q(:,k2) - Q(:,k1)) / (tv(k2) - tv(k1));
        p = abort_qmul([Q(1,k); -Q(2:4,k)], qd);
        W(:,k) = 2 * p(2:4);
    end
end

function r = abort_qmul(p, q)
    r = [p(1)*q(1)-p(2)*q(2)-p(3)*q(3)-p(4)*q(4); ...
         p(1)*q(2)+p(2)*q(1)+p(3)*q(4)-p(4)*q(3); ...
         p(1)*q(3)-p(2)*q(4)+p(3)*q(1)+p(4)*q(2); ...
         p(1)*q(4)+p(2)*q(3)-p(3)*q(2)+p(4)*q(1)];
end

function S = abort_fill(S, d)
    f = fieldnames(d);
    for i = 1:numel(f), if ~isfield(S, f{i}), S.(f{i}) = d.(f{i}); end, end
end

function x = abort_get(S, n, def)
    if isfield(S, n) && ~isempty(S.(n)), x = S.(n)(:); else, x = def; end
end
