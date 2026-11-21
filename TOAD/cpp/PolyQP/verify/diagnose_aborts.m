function diagnose_aborts(tag)
% Verification only (never part of the flight build).
% Phase-1 abort diagnostics on the MATLAB PolyQP currently on the path:
%   1.1 infeasibility root cause, 1.2 conditioning envelope, 1.3 roll units (J_zz, planned U(4)),
%   1.4 feedforward step at abort, 1.5 knot-1 P-offset magnitude.
% Cases are the 200 + 200 random abort states frozen in results/baseline_cases_slim.mat.
% Results go to results/diag_<tag>.mat; diagnosis_report.md is rewritten from every results/diag_*.mat
% present (baseline first), so a run on the fixed code adds the "after" columns.
if nargin < 1
    tag = 'baseline';
end
here    = fileparts(mfilename('fullpath'));
toadDir = fileparts(fileparts(fileparts(here)));
addpath(toadDir);
addpath(fullfile(toadDir, 'Guidance'));   % AbortGuidance.m (the abort block)
addpath(fullfile(toadDir, 'Guidance', 'PolynomialQP'));
addpath(fullfile(toadDir, 'Helper'));
warnState = warning('off', 'all');
cleanupWarn = onCleanup(@() warning(warnState));

C = {LoadTOADParams(0), LoadTOADParams(1)};
vehName = {'ASTRA', 'TOAD'};
Lim = {buildVehicleLimits(C{1}), buildVehicleLimits(C{2})};
B = load(fullfile(here, 'results', 'baseline_cases_slim.mat'));
utilNames = {'Vel', 'Accel', 'ThrustHi', 'ThrustLo', 'Gimbal', 'GimbalRate', 'ThrustRate', 'Roll', 'Tilt'};

D = struct();
D.tag = tag;
D.date = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
[~, rev] = system(sprintf('git -C "%s" rev-parse --short HEAD', toadDir));
[~, dirty] = system(sprintf('git -C "%s" status --porcelain -- Guidance Helper cpp/PolyQP', toadDir));
D.rev = strtrim(rev);
D.dirty = ~isempty(strtrim(dirty));
D.utilNames = utilNames;
D.vehName = vehName;
D.rho = [Lim{1}.Rho, Lim{2}.Rho];
D.f07 = contains(fileread(which('solveLeg')), 'FIX-F07');   % the code under test reparametrises shift-invariant segments

%% 1.1 / 1.4 / 1.5 (plan part): the frozen random abort cases
fprintf('[diag %s] abort cases\n', tag);
nC = numel(B.caseLog.abort);
cs = repmat(struct('veh', 0, 'q', 0, 'X_est', [], 'status', '', 'feasible', NaN, 'outerIter', NaN, 'T_total', NaN, ...
    'segT', [], 'maxUtil', [], 'domType', 0, 'domVal', NaN, 'domSeg', 0, 'domAtT0', false, 'util0', [], ...
    'overCap', false, 'overCapAxes', false(3, 1), 'v0Util', NaN, 'staticD3Feasible', NaN, 'staticD3Max', NaN, 'staticD3Type', 0, ...
    'dT', NaN, 'dTrel', NaN, 'dThrustDir', NaN, 'gimbal0', NaN, 'roll0', NaN, 'dAtt', NaN, ...
    'cgPosErr0', NaN, 'cgVelErr0', NaN, 'maxAbsRoll', NaN, 'maxRollRate', NaN, 'thrustSlew', NaN, 'gimbalSlew', NaN), 1, nC);
for i = 1:nC
    a = B.caseLog.abort{i};
    v = a.veh;
    c = C{v};
    L = Lim{v};
    X_est = a.X_est;
    cs(i).veh = v;
    cs(i).q = a.q;
    cs(i).X_est = X_est;
    W = abortWaypoints(X_est);
    Plan0 = resolveWaypoints(W, 5, L);
    cs(i).overCapAxes = abs(X_est(8:10)) > Plan0.SegMaxVel(:, 1);
    cs(i).overCap = any(cs(i).overCapAxes);
    cs(i).v0Util = max(abs(X_est(8:10)) ./ Plan0.SegMaxVel(:, 1));
    try
        [IG, QT, RP] = PolyQP(W, c, struct('R', 5));
    catch ME
        cs(i).status = ME.identifier;
        continue;
    end
    cs(i).status = 'ok';
    cs(i).feasible = RP.Feasible;
    cs(i).outerIter = RP.OuterIterations;
    cs(i).T_total = RP.T_total;
    cs(i).segT = RP.SegTimes;
    cs(i).maxUtil = [RP.UtilVel; RP.UtilAccel; RP.UtilThrustHi; RP.UtilThrustLo; RP.UtilGimbal; RP.UtilGimbalRate; RP.UtilThrustRate; RP.UtilRoll; RP.UtilTilt];

    PlanE = Plan0;
    M = QT.Count;
    PlanE.SegT(1:M) = QT.T(1:M);
    [EO, ~] = evaluateConstraints(QT, PlanE, c, L);
    [uk, segk] = perSampleUtil(QT, PlanE, c, L, PlanE.SegMaxVel);
    [cs(i).domVal, cs(i).domType] = max(EO.MaxUtilPerType);
    [~, cs(i).domSeg] = max(EO.SegUtil(cs(i).domType, :));
    cs(i).util0 = uk(:, 1);
    cs(i).domAtT0 = uk(cs(i).domType, 1) > 1.0 + 1e-3;
    % Static D3 counterfactual on the final plan: segment velocity cap raised to |v_BC| where the start knot pins a faster velocity
    capD3 = PlanE.SegMaxVel;
    capD3(:, 1) = max(capD3(:, 1), abs(X_est(8:10)));
    ukD3 = perSampleUtil(QT, PlanE, c, L, capD3);
    [mxD3, kD3] = max(max(ukD3, [], 2));
    cs(i).staticD3Max = mxD3;
    cs(i).staticD3Type = kD3;
    cs(i).staticD3Feasible = cs(i).staticD3Max <= 1.0 + 1e-3;

    % 1.4 feedforward / attitude step at t = 0 against the hover command implied by X_est
    m0 = c.m_dry + X_est(14) + X_est(15);
    T_hover = m0 * c.g;
    b3_est = quatRot(X_est(1:4))' * [0; 0; 1];
    u1 = IG.U(:, 1);
    tB = [cos(u1(1)) * sin(u1(2)); -sin(u1(1)); cos(u1(1)) * cos(u1(2))];
    tI = quatRot(IG.X(1:4, 1))' * tB;
    cs(i).dT = u1(3) - T_hover;
    cs(i).dTrel = (u1(3) - T_hover) / T_hover;
    cs(i).dThrustDir = rad2deg(acos(max(-1, min(1, dot(tI, b3_est)))));
    cs(i).gimbal0 = rad2deg(max(abs(u1(1:2))));
    cs(i).roll0 = u1(4);
    cs(i).dAtt = rad2deg(2 * acos(min(1, abs(dot(IG.X(1:4, 1), X_est(1:4))))));
    cs(i).thrustSlew = abs(cs(i).dT) / L.ThrustRateLim;
    cs(i).gimbalSlew = deg2rad(cs(i).gimbal0) / L.GimbalRateLim;
    % 1.5 (plan part): CG error of the plan at t = 0
    cs(i).cgPosErr0 = norm(IG.X(5:7, 1) - X_est(5:7));
    cs(i).cgVelErr0 = norm(IG.X(8:10, 1) - X_est(8:10));
    % 1.3: planned roll channel
    cs(i).maxAbsRoll = max(abs(IG.U(4, :)));
    cs(i).maxRollRate = max(abs(diff(IG.U(4, :)))) / (IG.Time(2) - IG.Time(1));
end
D.cases = cs;

%% 1.2 Conditioning envelope
fprintf('[diag %s] conditioning sweep\n', tag);
z0s = [0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 180, 200, 225, 250, 275, 300];
v0s = [0, 0, 0; 3, -2, -3; -3.5, 3.5, 4.5; 2, 1, -4.5].';
xy0 = [5; -3];
sw = repmat(struct('veh', 0, 'z0', NaN, 'v0', [], 'status', '', 'T0', [], 'Tf', [], 'condInit', [], 'rcondInit', [], ...
    'condFinal', [], 'rcondFinal', [], 'feasible', NaN), 1, 0);
for v = 1:2
    for iz = 1:numel(z0s)
        for iv = 1:size(v0s, 2)
            X_est = [1; 0; 0; 0; xy0; z0s(iz); v0s(:, iv); 0; 0; 0; C{v}.OxMass; C{v}.FuMass];
            W = abortWaypoints(X_est);
            P0 = resolveWaypoints(W, 5, Lim{v});
            s = struct('veh', v, 'z0', z0s(iz), 'v0', v0s(:, iv), 'status', '', 'T0', P0.SegT(1:P0.SegCount), 'Tf', [], ...
                'condInit', [], 'rcondInit', [], 'condFinal', [], 'rcondFinal', [], 'feasible', NaN);
            [s.condInit, s.rcondInit] = abortLegConditioning(P0, P0.SegT, Lim{v}, D.f07);
            try
                [~, QT, RP] = PolyQP(W, C{v}, struct('R', 5));
                s.status = 'ok';
                s.Tf = QT.T(1:QT.Count);
                s.feasible = RP.Feasible;
                [s.condFinal, s.rcondFinal] = abortLegConditioning(P0, QT.T, Lim{v}, D.f07);
            catch ME
                s.status = ME.identifier;
            end
            sw(end + 1) = s; %#ok<AGROW>
        end
    end
end
D.sweep = sw;
D.sweepZ0 = z0s;
D.sweepV0 = v0s;

%% 1.3 Roll units: J_zz per vehicle
for v = 1:2
    J_full = ComputeJtot(C{v}.OxMass, C{v}.FuMass, C{v});
    J_empty = ComputeJtot(0, 0, C{v});
    D.Jzz(v, :) = [J_full(3, 3), J_empty(3, 3)];
    D.RollLimit(v) = Lim{v}.RollLimit;
    D.RollRateLim(v) = Lim{v}.RollRateLim;
end

%% 1.5 Knot-1 P-offset magnitude for sampled tilts (<= 30 deg) and body rates (<= 1 rad/s)
rs = RandStream('mt19937ar', 'Seed', 15);
nS = 5000;
for v = 1:2
    rho = Lim{v}.Rho;
    ePos = zeros(1, nS);
    eVel = zeros(1, nS);
    for k = 1:nS
        tilt = deg2rad(30) * rand(rs);
        az = 2 * pi * rand(rs);
        ax = [cos(az); sin(az); 0];
        q = [cos(tilt / 2); sin(tilt / 2) * ax];
        w = randn(rs, 3, 1);
        w = w / norm(w) * rand(rs);
        C_IB = quatRot(q)';
        b3 = C_IB * [0; 0; 1];
        b3dot = cross(C_IB * w, b3);
        ePos(k) = rho * norm(b3 - [0; 0; 1]);
        eVel(k) = rho * norm(b3dot);
    end
    D.f03(v) = struct('ePos', ePos, 'eVel', eVel);
end

if ~isfolder(fullfile(here, 'results'))
    mkdir(fullfile(here, 'results'));
end
save(fullfile(here, 'results', ['diag_', tag, '.mat']), 'D', '-v7');
writeDiagnosisReport(here);
fprintf('[diag %s] done: %s\n', tag, fullfile(here, 'diagnosis_report.md'));
end

% ================================================================ helpers

function [cnd, rc] = abortLegConditioning(Plan, SegT, L, useReduced)
% Per-axis cond and rcond of R_PP_eq for the abort leg, FixP built as solveAllLegs pass 1. useReduced: the matrix the
% FIX-F07 solveLeg factorises; otherwise the pre-fix matrix (no reparametrisation).
    r = 5;
    K = Plan.KnotCount;
    FixP = Plan.FixCg(:, 1:r, 1:K);
    FixP(3, 1, ~isnan(FixP(3, 1, :))) = FixP(3, 1, ~isnan(FixP(3, 1, :))) + L.Rho;
    [~, ~, H0] = assembleQP(r);
    m1 = Plan.LegStartSeg(1);
    m2 = Plan.LegEndSeg(Plan.LegCount);
    FixLeg = FixP(:, :, m1:(m2 + 1));
    [Req, ~, ReqRaw] = rppEqAsSolved(FixLeg, SegT(m1:m2), H0);
    if ~useReduced
        Req = ReqRaw;
    end
    cnd = nan(3, 1);
    rc = nan(3, 1);
    for d = 1:3
        if ~isempty(Req{d})
            cnd(d) = cond(Req{d});
            rc(d) = rcond(Req{d});
        end
    end
end

function [UtilK, seg_idx] = perSampleUtil(QPTraj, Plan, Constants, Limits, SegMaxVel)
% Per-grid-sample utilisation (9 x nGrid), the per-k computation of evaluateConstraints.m (hold samples NaN).
    M = QPTraj.Count;
    T_total = QPTraj.T_total;
    dt_eval = 0.5 * Limits.DtSample;
    nGrid = min(3999, max(41, ceil(T_total / dt_eval) + 1));
    t_grid = linspace(0, T_total, nGrid);
    dt_step = t_grid(2) - t_grid(1);
    [X, Wrench, AccelCg] = flatnessInvert(QPTraj, t_grid, Plan.HDGRef, Constants, Limits);
    U = allocateControls(Wrench, Limits);
    U_rate = zeros(4, nGrid);
    U_rate(:, 1:(nGrid - 1)) = diff(U, 1, 2) / dt_step;
    U_rate(:, nGrid) = U_rate(:, nGrid - 1);
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
    for k = 1:(nGrid - 1)
        m1 = seg_idx(k);
        m2 = seg_idx(k + 1);
        if QPTraj.IsHold(m1) || QPTraj.IsHold(m2) || Plan.IsCruise(m1) || Plan.IsCruise(m2)
            U_rate(:, k) = 0;
        end
    end
    U_rate(:, nGrid) = U_rate(:, max(1, nGrid - 1));
    UtilK = nan(9, nGrid);
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
        u_vel = max(vk ./ SegMaxVel(:, m));
        u_acc = max(abs(ak) ./ Plan.SegMaxAccel(:, m));
        T_k = uk(3);
        u_thi = T_k / Limits.T_hi;
        u_tlo = Limits.T_lo / max(1e-3, T_k);
        u_gim = max(abs(uk(1:2))) / Limits.GimbalLimit;
        u_grate = max(urk(1:2)) / Limits.GimbalRateLim;
        u_trate = urk(3) / Limits.ThrustRateLim;
        u_roll = max(abs(uk(4)) / Limits.RollLimit, urk(4) / Limits.RollRateLim);
        b3_z = 1 - 2 * (qk(2)^2 + qk(3)^2);
        u_tilt = acos(max(-1, min(1, b3_z))) / Limits.MaxTilt;
        if Plan.IsCruise(m)
            u_grate = 0;
            u_trate = 0;
        end
        UtilK(:, k) = [u_vel; u_acc; u_thi; u_tlo; u_gim; u_grate; u_trate; u_roll; u_tilt];
    end
end

function s = pct(x, p)
    x = sort(x(~isnan(x)));
    if isempty(x)
        s = NaN;
        return;
    end
    s = x(max(1, ceil(p * numel(x))));
end

function t = dist3(x, fmt)
% "median / 95th / max"
    t = sprintf([fmt, ' / ', fmt, ' / ', fmt], pct(x, 0.5), pct(x, 0.95), max(x));
end

function writeDiagnosisReport(here)
    files = dir(fullfile(here, 'results', 'diag_*.mat'));
    tags = cellfun(@(f) f(6:end - 4), {files.name}, 'UniformOutput', false);
    order = [find(strcmp(tags, 'baseline')), find(~strcmp(tags, 'baseline'))];
    R = cell(1, numel(order));
    for i = 1:numel(order)
        S = load(fullfile(here, 'results', files(order(i)).name));
        R{i} = S.D;
    end
    fid = fopen(fullfile(here, 'diagnosis_report.md'), 'w');
    closer = onCleanup(@() fclose(fid));
    p = @(varargin) fprintf(fid, varargin{:});
    p('# Abort diagnosis (Phase 1)\n\n');
    p('Generated by `diagnose_aborts.m`. Runs: ');
    for i = 1:numel(R)
        p('`%s` (%s, commit %s%s)', R{i}.tag, R{i}.date, R{i}.rev, ternary(R{i}.dirty, ' + uncommitted changes', ''));
        if i < numel(R)
            p('; ');
        else
            p('.\n\n');
        end
    end
    p('All numbers are measured with the MATLAB code (the specification). The case set is the 200 ASTRA + 200 TOAD random abort states of `verify_polyqp.m` (seed 42), frozen in `results/baseline_cases_slim.mat`: x, y in [-20, 20] m, z in [0.5, 40] m, vx, vy in [-5, 5] m/s, vz in [-6, 4] m/s, upright, zero rates, full tanks.\n\n');
    vn = R{1}.vehName;
    un = R{1}.utilNames;

    %% 1.1
    p('## 1.1 Infeasibility root cause\n\n');
    p('Segment 1 (ABT START to ABT GND) takes its velocity cap from waypoint 2, the `TOADWaypoint` default MaxVel = [4; 4; 5] m/s. A case is *over-cap* if |v0| exceeds that cap on any axis. Predicted over-cap rate: 1 - 0.8 x 0.8 x 0.9 = 42.4%%.\n\n');
    p('| Run | Vehicle | Cases | PolyQP errors | Infeasible | Over-cap | Over-cap & infeasible | Over-cap & feasible | Within-cap & infeasible |\n|---|---|---|---|---|---|---|---|---|\n');
    for i = 1:numel(R)
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v);
            ok = strcmp({c.status}, 'ok');
            inf_ = ok & [c.feasible] == 0;
            oc = [c.overCap];
            p('| %s | %s | %d | %d | %d (%.1f%%) | %d (%.1f%%) | %d | %d | %d |\n', R{i}.tag, vn{v}, numel(c), nnz(~ok), nnz(inf_), 100 * nnz(inf_) / numel(c), ...
                nnz(oc), 100 * nnz(oc) / numel(c), nnz(oc & inf_), nnz(oc & ok & [c.feasible] == 1), nnz(~oc & inf_));
        end
    end
    p('\n');
    for i = 1:numel(R)
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v);
            ax = [c.overCapAxes];
            p('- `%s` %s over-cap axes: vx %d, vy %d, vz %d (cases may exceed on several).', R{i}.tag, vn{v}, nnz(ax(1, :)), nnz(ax(2, :)), nnz(ax(3, :)));
            cf = c([c.overCap] & strcmp({c.status}, 'ok') & [c.feasible] == 1);
            if numel(cf) > 10
                p(' Over-cap but feasible: %d cases.', numel(cf));
            elseif ~isempty(cf)
                p(' Over-cap but feasible: %s (v0 utilisation %s; the feasibility test allows 1 + 1e-3).', strjoin(arrayfun(@(x) sprintf('#%d', x.q), cf, 'UniformOutput', false), ', '), ...
                    strjoin(arrayfun(@(x) sprintf('%.4f', x.v0Util), cf, 'UniformOutput', false), ', '));
            end
            p('\n');
        end
    end
    p('\n### Infeasible cases by dominant utilisation (largest entry of `Report.MaxUtilPerType`)\n\n');
    p('"At t = 0" counts cases whose dominant utilisation already exceeds 1 + 1e-3 at the first grid sample. "Static D3" re-evaluates the final plan with segment 1''s velocity cap raised to |v_BC| (the D3 rule, no re-planning): it shows whether anything other than the pinned start velocity is violated.\n\n');
    for i = 1:numel(R)
        p('**Run `%s`**\n\n', R{i}.tag);
        p('| Vehicle | Group | Dominant type | Cases | Dominant value median / max | Segment 1 / 2 | Above 1 at t = 0 | Static D3 feasible | Static D3 still violated (dominant type: count, max) |\n|---|---|---|---|---|---|---|---|---|\n');
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v & strcmp({R{i}.cases.status}, 'ok'));
            c = c([c.feasible] == 0);
            grp = {'over-cap', [c.overCap]; 'within cap', ~[c.overCap]};
            for g = 1:2
                cg = c(grp{g, 2});
                if isempty(cg)
                    p('| %s | %s | - | 0 | - | - | - | - | - |\n', vn{v}, grp{g, 1});
                    continue;
                end
                for ty = unique([cg.domType])
                    ct = cg([cg.domType] == ty);
                    cr = ct([ct.staticD3Feasible] == 0);
                    rem = strjoin(arrayfun(@(t) sprintf('%s: %d, %.3f', un{t}, nnz([cr.staticD3Type] == t), max([cr([cr.staticD3Type] == t).staticD3Max])), unique([cr.staticD3Type]), 'UniformOutput', false), '; ');
                    p('| %s | %s | %s | %d | %.3f / %.3f | %d / %d | %d | %d | %s |\n', vn{v}, grp{g, 1}, un{ty}, numel(ct), median([ct.domVal]), max([ct.domVal]), ...
                        nnz([ct.domSeg] == 1), nnz([ct.domSeg] == 2), nnz([ct.domAtT0]), nnz([ct.staticD3Feasible] == 1), rem);
                end
            end
        end
        p('\n');
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v & strcmp({R{i}.cases.status}, 'ok'));
            ci = c([c.feasible] == 0 & ~[c.overCap]);
            if isempty(ci)
                continue;
            end
            p('Within-cap infeasible %s cases (%s):\n\n| # | X_est(5:10) | Dominant | Value | Seg | At t = 0 | Util at t = 0 (> 1 listed) | T_total (s) | OuterIter |\n|---|---|---|---|---|---|---|---|---|\n', vn{v}, R{i}.tag);
            for j = 1:numel(ci)
                hi = find(ci(j).util0 > 1 + 1e-3);
                hs = strjoin(arrayfun(@(t) sprintf('%s %.3f', un{t}, ci(j).util0(t)), hi(:).', 'UniformOutput', false), ', ');
                p('| %d | %s | %s | %.3f | %d | %s | %s | %.2f | %d |\n', ci(j).q, mat2str(ci(j).X_est(5:10).', 3), un{ci(j).domType}, ci(j).domVal, ...
                    ci(j).domSeg, ternary(ci(j).domAtT0, 'yes', 'no'), hs, ci(j).T_total, ci(j).outerIter);
            end
            p('\n');
        end
    end

    %% 1.2
    p('## 1.2 Conditioning envelope\n\n');
    p('Sweep: start (5, -3, z0), upright, z0 in %s m, start velocities (columns) %s m/s (all within the caps). T1, T2 are the abort leg''s two segment times; "init" is the resolveWaypoints initial guess (where pass 1 solves first), "final" the returned plan. cond and rcond are of the equilibrated R_PP_eq per axis (x, y; z is well conditioned) that `solveLeg.m` factorises: the full matrix before FIX-F07, the reduced (reparametrised) matrix after it (`rppEqAsSolved.m`). The `SingularRPP` guard trips at rcond < 1e-13.\n\n', mat2str(R{1}.sweepZ0), mat2str(R{1}.sweepV0.'));
    for i = 1:numel(R)
        sw = R{i}.sweep;
        p('**Run `%s`**\n\n', R{i}.tag);
        p('| Vehicle | z0 (m) | T1/T2 init (s) | T1/T2 init ratio | max cond x,y (init) | min rcond x,y (init) | PolyQP status (4 velocities) | T1/T2 final ratio (max) | max cond x,y (final) | max cond z |\n|---|---|---|---|---|---|---|---|---|---|\n');
        for v = 1:2
            for z0 = R{i}.sweepZ0
                s = sw([sw.veh] == v & [sw.z0] == z0);
                T0 = vertcat(s.T0);
                ci = arrayfun(@(x) max(x.condInit(1:2)), s);
                ri = arrayfun(@(x) min(x.rcondInit(1:2)), s);
                czi = arrayfun(@(x) x.condInit(3), s);
                st = strjoin(cellfun(@(x) strrep(x, 'PolyQP:', ''), {s.status}, 'UniformOutput', false), ',');
                okm = strcmp({s.status}, 'ok');
                if any(okm)
                    rf = max(arrayfun(@(x) x.Tf(1) / x.Tf(2), s(okm)));
                    cf = max(arrayfun(@(x) max(x.condFinal(1:2)), s(okm)));
                    czf = max(arrayfun(@(x) x.condFinal(3), s(okm)));
                    rfs = sprintf('%.1f', rf);
                    cfs = sprintf('%.2e', cf);
                else
                    rfs = '-'; cfs = '-'; czf = NaN;
                end
                p('| %s | %g | %.1f-%.1f / %.2f-%.2f | %.1f-%.1f | %.2e | %.2e | %s | %s | %s | %.2e |\n', vn{v}, z0, min(T0(:, 1)), max(T0(:, 1)), min(T0(:, 2)), max(T0(:, 2)), ...
                    min(T0(:, 1) ./ T0(:, 2)), max(T0(:, 1) ./ T0(:, 2)), max(ci), min(ri), st, rfs, cfs, max([czi, czf]));
            end
        end
        p('\n');
        for v = 1:2
            s = sw([sw.veh] == v);
            bad = ~strcmp({s.status}, 'ok');
            tripInit = arrayfun(@(x) min(x.rcondInit(1:2)) < 1e-13, s);
            p('- %s (`%s`): first z0 where PolyQP returns SingularRPP: %s; first z0 where rcond(R_PP_eq) at the initial times is below 1e-13: %s. PolyQP errors in the sweep: %d of %d. Max cond(R_PP_eq) over successful solves (final times, any axis): %.2e.\n', ...
                vn{v}, R{i}.tag, firstZ(s, bad & contains({s.status}, 'SingularRPP')), firstZ(s, tripInit), nnz(bad), numel(s), ...
                max(arrayfun(@(x) max([x.condFinal; -Inf]), s)));
        end
        allS = sw;
        lr = arrayfun(@(x) log10(x.T0(1) / x.T0(2)), allS);
        lc = arrayfun(@(x) log10(max(x.condInit(1:2))), allS);
        okf = isfinite(lr) & isfinite(lc) & lr > 0.5;
        if nnz(okf) > 3
            pf = polyfit(lr(okf), lc(okf), 1);
            p('- Fit over the sweep (T1/T2 > 3, initial times, both vehicles): log10 cond = %.2f x log10(T1/T2) + %.2f, i.e. cond grows as (T1/T2)^%.1f. Predicted rcond = 1e-13 (cond ~ 1e13) at T1/T2 = %.1f.\n', pf(1), pf(2), pf(1), 10 ^ ((13 - pf(2)) / pf(1)));
        end
        p('\n');
    end

    %% 1.3
    p('## 1.3 F30: roll units\n\n');
    p('| Vehicle | J_zz full tanks (kg m^2) | J_zz empty (kg m^2) | PolyQP RollLimit (on U(4)) | PolyQP RollRateLim (on dU(4)/dt) |\n|---|---|---|---|---|\n');
    for v = 1:2
        p('| %s | %.4g | %.4g | %.3g | %.3g |\n', vn{v}, R{1}.Jzz(v, 1), R{1}.Jzz(v, 2), R{1}.RollLimit(v), R{1}.RollRateLim(v));
    end
    p('\nPlanned roll channel over the abort cases (max over each plan''s samples):\n\n| Run | Vehicle | max abs U(4) median / 95th / max (N m) | implied roll accel U(4)/J_zz max (rad/s^2) | max abs dU(4)/dt (N m/s) |\n|---|---|---|---|---|\n');
    for i = 1:numel(R)
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v & strcmp({R{i}.cases.status}, 'ok'));
            p('| %s | %s | %s | %.3g | %.3g |\n', R{i}.tag, vn{v}, dist3([c.maxAbsRoll], '%.3g'), max([c.maxAbsRoll]) / R{1}.Jzz(v, 1), max([c.maxRollRate]));
        end
    end
    p('\n%s\n\n', f30Trace());

    %% 1.4
    p('## 1.4 F02: feedforward step at the abort instant\n\n');
    p('Plan t = 0 values against the command implied by `X_est` (hover thrust m g along the current b3 = C_IB e3, zero gimbal, zero roll). The thrust direction is the plan''s inertial thrust vector C_IB(q_plan(0)) t_B(theta, phi). Slew time = step / rate limit (`ThrustRateLim`, `GimbalRateLim`). Values are median / 95th / max.\n\n');
    p('| Run | Vehicle | Thrust step (N) | Thrust step / hover | Thrust slew time (s) | Thrust direction step (deg) | Gimbal step (deg) | Gimbal slew time (s) | Roll U(4) at t = 0 (N m) | Attitude step (deg) |\n|---|---|---|---|---|---|---|---|---|---|\n');
    for i = 1:numel(R)
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v & strcmp({R{i}.cases.status}, 'ok'));
            p('| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', R{i}.tag, vn{v}, dist3(abs([c.dT]), '%.3g'), dist3(abs([c.dTrel]), '%.3f'), dist3([c.thrustSlew], '%.3g'), ...
                dist3([c.dThrustDir], '%.2f'), dist3([c.gimbal0], '%.2f'), dist3([c.gimbalSlew], '%.3f'), dist3(abs([c.roll0]), '%.2g'), dist3([c.dAtt], '%.2f'));
        end
    end
    p('\nSigned thrust step (plan minus hover), min / max: ');
    for v = 1:2
        c = R{1}.cases([R{1}.cases.veh] == v & strcmp({R{1}.cases.status}, 'ok'));
        p('%s %.3g / %.3g N (%.1f%% / %.1f%% of hover)', vn{v}, min([c.dT]), max([c.dT]), 100 * min([c.dTrel]), 100 * max([c.dTrel]));
        if v == 1
            p('; ');
        else
            p('.\n\n');
        end
    end
    p('Limits for comparison: gimbal %.1f deg (planner %.1f deg), gimbal rate %.0f deg/s, thrust rate %g N/s (ASTRA) and %g N/s (TOAD), actuator lags tau = 0.15 s (thrust) and 0.075 s (roll), gimbal 2nd order (a = 800, b = 70) in `Actuator_Models_2.slx`; controller tilt cone 25 deg, planner MaxTilt 30 deg.\n\n', ...
        10, rad2deg(buildGimbal()), 30, 30, 1500);

    %% 1.5
    p('## 1.5 F03: knot-1 P-offset magnitude\n\n');
    p('No MC_v5 dispersion outputs are in the repo (`Analysis/MonteCarlo_v5.m` saves to `Analysis/Monte Carlo Runs/`, which does not exist), so tilts are sampled uniformly in [0, 30] deg about a random horizontal axis and body rates uniformly in magnitude up to 1 rad/s (5000 draws per vehicle). Position error = rho |b3 - e3|, velocity error = rho |omega_I x b3|.\n\n');
    p('| Vehicle | rho (m) | Position error median / 95th / max (m) | Velocity error median / 95th / max (m/s) |\n|---|---|---|---|\n');
    for v = 1:2
        p('| %s | %.4f | %s | %s |\n', vn{v}, R{1}.rho(v), dist3(R{1}.f03(v).ePos, '%.4f'), dist3(R{1}.f03(v).eVel, '%.4f'));
    end
    p('\nPlan CG at t = 0 against `X_est` for the (upright, zero-rate) abort cases, |X(5:7,1) - X_est(5:7)| and |X(8:10,1) - X_est(8:10)|. With b3_actual = e3 the P-space BC is exact, so this residual is rho |e3 - b3_plan(0)|: the plan''s own t = 0 tilt, set by its free initial acceleration.\n\n');
    p('| Run | Vehicle | CG position error median / 95th / max (m) | CG velocity error median / 95th / max (m/s) |\n|---|---|---|---|\n');
    for i = 1:numel(R)
        for v = 1:2
            c = R{i}.cases([R{i}.cases.veh] == v & strcmp({R{i}.cases.status}, 'ok'));
            p('| %s | %s | %s | %s |\n', R{i}.tag, vn{v}, dist3([c.cgPosErr0], '%.3g'), dist3([c.cgVelErr0], '%.3g'));
        end
    end
    p('\n');
end

function g = buildGimbal()
    g = 0.9 * deg2rad(10);
end

function s = firstZ(sw, mask)
    if any(mask)
        s = sprintf('%g m', min([sw(mask).z0]));
    else
        s = 'none up to 300 m';
    end
end

function out = ternary(c, a, b)
    if c
        out = a;
    else
        out = b;
    end
end

function t = f30Trace()
t = strjoin({
'Trace of U(4) from the planner to the plant (file:line at `abort-fixes`):'
''
'| Stage | File:line | Quantity assumed | Evidence |'
'|---|---|---|---|'
'| Planner output | `Guidance/PolynomialQP/allocateControls.m:24`, `:36`, `:44` | roll **torque**, N m | `roll_cmd = tauk(3) * (T_mag / Fz)` (ASTRA) / `roll_cmd = tauk(3)` (TOAD), with `tauk = J_tot * wdk + cross(wk, J_tot * wk)` (`flatnessInvert.m:130`): an Euler-equation moment. |'
'| Planner limit | `Guidance/PolynomialQP/buildVehicleLimits.m:8`, `:46`; `evaluateConstraints.m:107` | rate-named constants applied to a torque | `MaxRollRate = 5`, `RollLimit = 0.9 * MaxRollRate`; `u_roll = max(abs(uk(4)) / RollLimit, urk(4) / RollRateLim)`: 4.5 N m and 5 N m/s in effect. |'
'| Abort hand-off | `Guidance/TrajectoryManager.m:58`, `:139` | passes U(4) through unchanged | `Inputs = [InitialGuess.U''; ...]`, then linear interpolation into `U`. |'
'| Flight controller | `cpp/ASTRAv2_Controller.cpp:185-191` | label says rad/s^2; arithmetic is N m | `U_ff_chan(2) = U_ff(3)` is summed with `Trim(2) = Trim_FB(2) + U_dist_req(2)`, where `U_dist(2) = torque(2) = ang_accel_dist(2) * J(2)` (`:96`, `:109`) is a torque, and `Trim_FB = -K_rot * X_err_rot` with `K_rot` from an LQR whose roll input enters the moment balance. Output clamped to +-4 and labelled `// Roll (rad/s^2)`. |'
'| LQR design | `Controls/TOAD_Controller_Gen.m:8`, `:61`, `:67` | label says rad/s^2; model is N m | `syms roll % roll input (rad/sec^2)`, but `MB = ... + (T_vec * roll)` and `omegaBdot = J_tot \ (MB - ...)`: roll is a moment. Same in `Controls/Time Varying LQI/RicattiRecursion.m:36-38`. |'
'| Controller docs | `cpp/README.md:67`, `:99` | rad/s^2 | "`U_ff`: ... roll accel", "`roll_rad_sec_squared`: Roll acceleration command, rad/s^2". |'
'| Actuator model | `Flight Dynamics/Simulink Models/Actuator_Models_2.slx`, MATLAB Function `Actuators` | unit-agnostic first-order lag | `UDOT(6) = (1/(tau/2)) * (U_TRG(4) - U(6))`, tau = 0.15 s; no scaling. |'
'| Plant | `Flight Dynamics/FlightDynamicsGen.m:10`, `:107`, `:109`, `:114`; `Dynamics Files/TOAD_Dynamics.m:5` | roll **torque**, N m | `syms roll % roll input (Nm)`; `MB = zetaCross([0; 0; -CGz])*TB + roll * thrustDir` (ASTRA) / `+ [0; 0; roll]` (TOAD); `omegaBdot = J_tot \ netTau`. |'
''
'**Finding.** Every computation on the path (planner, LQR model, controller arithmetic, actuator, plant) treats the roll channel as a torque in N m. The planner is consistent with the plant. What is wrong is on the controller side and is documentation plus one limit: the `// Roll (rad/s^2)` label and the +-4 clamp in `ASTRAv2_Controller.cpp:191`, the `rad/sec^2` comment in `TOAD_Controller_Gen.m:8`, and `cpp/README.md:67`/`:99`. If the clamp was chosen as 4 rad/s^2, it is in fact 4 N m, which is 200 rad/s^2 on ASTRA (J_zz = 0.02) and 0.19 rad/s^2 on TOAD (J_zz = 20.6). Inside PolyQP, `RollLimit`/`RollRateLim` are built from a constant named `MaxRollRate` but are applied to the torque; their values (4.5 N m, 5 N m/s) are not derived from J_zz and the planner''s 4.5 exceeds the controller''s 4 clamp. This is the controller/plant owners'' call, so it is reported as BLOCKER B1 and not edited.'
}, newline);
end
