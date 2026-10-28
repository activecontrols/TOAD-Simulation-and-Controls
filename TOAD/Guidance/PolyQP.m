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
