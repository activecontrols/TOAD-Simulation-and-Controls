function verify_polyqp(runTag)
% Verification only (never part of the flight build).
% runTag (optional): if given, every abort case and solveLeg check (inputs, MATLAB and C++ outputs)
% is saved to results/<runTag>_cases.mat for before/after comparisons.
% Per-function I/O checks of the C++ PolyQP transpile (polyqp_mex)
% against the MATLAB .m files. Run from TOAD/cpp/PolyQP/verify after
%   build_polyqp_mex;
% Writes verify_polyqp_report.md next to this file and errors if any double-precision check fails.

if nargin < 1
    runTag = '';
end
here    = fileparts(mfilename('fullpath'));
toadDir = fileparts(fileparts(fileparts(here)));
addpath(toadDir);
addpath(fullfile(toadDir, 'Guidance'));   % AbortGuidance.m (the abort block)
addpath(fullfile(toadDir, 'Guidance', 'PolynomialQP'));
addpath(fullfile(toadDir, 'Helper'));
addpath(here);
clear polyqp_mex polyqp_mex_wide polyqp_mex_p3;

assert(exist('polyqp_mex', 'file') == 3, 'polyqp_mex not built: run build_polyqp_mex first.');
assert(exist('polyqp_mex_wide', 'file') == 3, 'polyqp_mex_wide not built: run build_polyqp_mex(''wide'') first.');
% STATIC: polyqp_mex has the flight (abort-sized) capacities and runs every abort case; polyqp_mex_wide is the same source with
% larger capacities (-D) for what exceeds the abort envelope by design: the missions, fuzzed legs, and grid-sized I/O checks.
mxF = @polyqp_mex;
capF = polyqp_mex('capacities');   % [MAX_WP, MAX_SEG, MAX_KNOTS, MAX_FREE, MAX_NODES] of the flight build
POLYQP_FLIGHT_KNOTS = capF(3);
mxW = @polyqp_mex_wide;
haveP3 = exist('polyqp_mex_p3', 'file') == 3;   % the last variable-size build (Phase 3), for the identity-padding check

warnState = warning('off', 'all');
cleanupWarn = onCleanup(@() warning(warnState));

checkAbortCopyIsVerbatim(toadDir, here);

rng(42);
tStart = tic;

% ---------------------------------------------------------------- tolerances
TOL_ALG = 1e-14;   % closed-form algebra (rel)
TOL_LIN = 1e-9;    % linear solves (normwise rel)
TOL_NL  = 1e-9;    % nonlinear chains (rel or abs)
TOL_E2E = 1e-6;    % end-to-end (abs)

rows = struct('name', {}, 'kind', {}, 'tol', {}, 'cases', {}, 'fails', {}, 'roundoff', {}, 'maxErr', {}, 'maxRatio', {}, ...
    'notes', {}, 'failDetail', {}, 'roundoffDetail', {});
ROUNDOFF_FACTOR = 10;   % a strict miss is roundoff-explained if err <= ROUNDOFF_FACTOR * (MATLAB self-sensitivity floor)
exactRef = {};          % vpa reference rows for solveLeg strict misses
statusCounts = struct('ok', 0, 'IllPosedQP', 0, 'SingularRPP', 0);
atMoved = 0;
caseLog = struct('abort', {{}}, 'solveLeg', {{}});

    function addRow(name, kind, tol)
        rows(end + 1) = struct('name', name, 'kind', kind, 'tol', tol, 'cases', 0, 'fails', 0, 'roundoff', 0, 'maxErr', 0, 'maxRatio', 0, ...
            'notes', '', 'failDetail', {{}}, 'roundoffDetail', {{}});
    end
    function rec(name, ok, err, detail, floorFn)
        % floorFn (optional): MATLAB self-sensitivity of the same metric to ~2 ulp input perturbations.
        i = find(strcmp({rows.name}, name), 1);
        assert(~isempty(i), 'unknown row %s', name);
        rows(i).cases = rows(i).cases + 1;
        if ~isnan(err)
            rows(i).maxErr = max(rows(i).maxErr, err);
        end
        if ~ok
            flo = NaN;
            if nargin >= 5 && isfinite(err)
                flo = floorFn();
            end
            if isfinite(flo) && err <= ROUNDOFF_FACTOR * flo
                rows(i).roundoff = rows(i).roundoff + 1;
                rows(i).maxRatio = max(rows(i).maxRatio, err / flo);
                if numel(rows(i).roundoffDetail) < 6
                    rows(i).roundoffDetail{end + 1} = sprintf('%s; floor=%.3g, err/floor=%.2f', detail, flo, err / flo);
                end
            else
                rows(i).fails = rows(i).fails + 1;
                if numel(rows(i).failDetail) < 8
                    rows(i).failDetail{end + 1} = sprintf('%s; floor=%.3g', detail, flo);
                end
            end
        end
    end
    function note(name, txt)
        i = find(strcmp({rows.name}, name), 1);
        rows(i).notes = txt;
    end

addRow('buildVehicleLimits',                 'closed-form (rel)',            TOL_ALG);
addRow('ComputeJtot',                        'closed-form (rel)',            TOL_ALG);
addRow('DCM_Quat_Conversion',                'closed-form (rel)',            TOL_ALG);
addRow('quatRot',                            'closed-form (rel)',            TOL_ALG);
addRow('TOADWaypoint',                       'closed-form (rel) + exact',    TOL_ALG);
addRow('assembleQP A0, Q0',                  'closed-form (rel)',            TOL_ALG);
addRow('assembleQP H0',                      'linear solve (rel)',           TOL_LIN);
addRow('resolveWaypoints structural',        'exact',                        0);
addRow('resolveWaypoints numeric',           'closed-form (rel)',            TOL_ALG);
addRow('resolveWaypoints error IDs',         'exact',                        0);
addRow('solveLeg (from solveAllLegs)',       'linear solve (rel)',           TOL_LIN);
addRow('solveLeg (fuzzed)',                  'linear solve (rel)',           TOL_LIN);
addRow('solveLeg error IDs',                 'exact',                        0);
addRow('allocateTimes',                      'linear solve (rel)',           TOL_LIN);
addRow('PolyEval',                           'closed-form (rel)',            TOL_ALG);
addRow('flatnessInvert',                     'nonlinear (rel or abs)',       TOL_NL);
addRow('allocateControls',                   'nonlinear (rel or abs)',       TOL_NL);
addRow('evaluateConstraints structural',     'exact',                        0);
addRow('evaluateConstraints numeric',        'nonlinear (rel or abs)',       TOL_NL);
addRow('sampleOutputs',                      'nonlinear (rel or abs)',       TOL_NL);
addRow('PolyQP structural / status',         'exact',                        0);
addRow('PolyQP end-to-end',                  'end-to-end (abs)',             TOL_E2E);
addRow('PolyQP error IDs',                   'exact',                        0);
addRow('AbortGuidance structural',           'exact',                        0);
addRow('AbortGuidance end-to-end',           'end-to-end (abs)',             TOL_E2E);
addRow('Missions structural / status',       'exact',                        0);
addRow('Missions end-to-end',                'end-to-end (abs)',             TOL_E2E);
addRow('Padded solve vs variable-size solve', 'linear solve (rel)',           TOL_LIN);
addRow('rcond guard on unpadded block',      'closed-form (rel)',            1e-12);

% ---------------------------------------------------------------- 1. constants
% ASTRA only: the C++ solver is ASTRA-only (vehicle 0, compiled-in constants, r = 5); the TOAD and r ~= 5 cases were dropped.
C = {LoadTOADParams(0)};
vehName = {'ASTRA'};
LIMITS_REMOVED = {'Vehicle', 'OxMass', 'FuMass', 'MaxMdot', 'OF'};   % Limits fields the ASTRA-only C++ does not carry

% ---------------------------------------------------------------- closed-form leaves
fprintf('[leaves] buildVehicleLimits, ComputeJtot, DCM_Quat_Conversion, TOADWaypoint, assembleQP\n');
LimML = cell(1, 2);
for v = 1:numel(C)
    LimML{v} = buildVehicleLimits(C{v});
    [Lm, st] = polyqp_mex('buildVehicleLimits', C{v});
    [ok, err, msg] = cmpStruct(rmfield(LimML{v}, LIMITS_REMOVED), Lm, 'rel', TOL_ALG, {'HasVehicle'});
    ok = ok && isempty(st);
    rec('buildVehicleLimits', ok, err, sprintf('%s: %s %s', vehName{v}, msg, st));
end

for v = 1:numel(C)
    c = C{v};
    massCases = [0, 0; c.OxMass, c.FuMass; 0.5 * c.OxMass, 0.25 * c.FuMass];   % ASTRA carries no propellant: every case is (0, 0)
    for q = 1:size(massCases, 1)
        [J_ml, CGz_ml] = ComputeJtot(massCases(q, 1), massCases(q, 2), c);
        [J_cx, CGz_cx, st] = polyqp_mex('ComputeJtot', massCases(q, 1), massCases(q, 2), c);
        [ok1, e1] = cmpNum(J_ml, J_cx, 'rel', TOL_ALG);
        [ok2, e2] = cmpNum(CGz_ml, CGz_cx, 'rel', TOL_ALG);
        rec('ComputeJtot', ok1 && ok2 && isempty(st), max(e1, e2), sprintf('%s m=(%g,%g) eJ=%.3g eCG=%.3g', vehName{v}, massCases(q, :), e1, e2));
    end
end

dcmCases = {};
dcmCases{end + 1} = eye(3);
for ang = [pi, pi - 0.1, pi - 1e-3, 2.5, 2.0944]
    dcmCases{end + 1} = rotx(ang); %#ok<AGROW>
    dcmCases{end + 1} = roty(ang); %#ok<AGROW>
    dcmCases{end + 1} = rotz(ang); %#ok<AGROW>
end
for q = 1:200
    qq = randn(4, 1);
    qq = qq / norm(qq);
    dcmCases{end + 1} = quatToDCM(qq); %#ok<AGROW>
end
nTraceNeg = 0;
for q = 1:numel(dcmCases)
    Rm = dcmCases{q};
    if trace(Rm) <= 0
        nTraceNeg = nTraceNeg + 1;
    end
    q_ml = DCM_Quat_Conversion(Rm);
    [q_cx, st] = polyqp_mex('DCM_Quat_Conversion', Rm);
    [ok, err] = cmpNum(q_ml, q_cx, 'rel', TOL_ALG);
    rec('DCM_Quat_Conversion', ok && isempty(st), err, sprintf('case %d trace=%.3g err=%.3g', q, trace(Rm), err));
end
note('DCM_Quat_Conversion', sprintf('%d cases with trace <= 0 (all three dominant-diagonal branches)', nTraceNeg));
qrs = RandStream('mt19937ar', 'Seed', 303);
for q = 1:100
    qq = randn(qrs, 4, 1);
    if q > 50
        qq = qq / norm(qq);
    end
    [Rc, st] = polyqp_mex('quatRot', qq);
    [ok, err] = cmpNum(quatRot(qq), Rc, 'rel', TOL_ALG);
    rec('quatRot', ok && isempty(st), err, sprintf('case %d err=%.3g', q, err));
end
note('quatRot', '50 unit and 50 non-unit quaternions (FIX-F03 uses it for b3 = C_IB e3)');

wpCases = {
    {[1; 2; 3]}
    {[1; 2; 3], 'Label', 'ABT START'}
    {[1; 2; 3], 'Label', 'ABT GND', 'HoldTime', 0.0}
    {[1; 2; 3], 'Label', 'A VERY LONG LABEL', 'HoldTime', 1.5, 'VelBC', [NaN; 1; NaN]}
    {[nan; nan; 2], 'Label', 'ABT GND', 'HoldTime', 0.0, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]}
    {[nan; nan; 0], 'Label', 'ABT END', 'HoldTime', 0.0, 'MaxVel', [3; 3; 1], 'StraightLine', true, 'VelBC', [0; 0; 0], 'AccelBC', [0;0;0], 'JerkBC', [0;0;0]}
    {[5; 0; 7], 'Label', 'FLT 2', 'VelBC', [NaN, NaN, 0]}
    {[0; 0; 0], 'HoldTime', NaN, 'AccelBC', [1; NaN; 2], 'JerkBC', [NaN; 3; NaN], 'MaxAccel', [1; 2; 3], 'HDGRef', [0; 1; 0], 'StraightLine', true}
    {[0; 0; 0], 'Label', '123456789'}
    {[0; 0; 0], 'Label', '1234567890'}
    {[1; 2; 3], 'VelBC', [1; 0; -1], 'Label', 'ABT START', 'B3BC', [0.1; -0.2; 0.97], 'B3DotBC', [0.3; 0.1; 0]}
    };
for q = 1:numel(wpCases)
    args = wpCases{q};
    w_ml = TOADWaypoint(args{:});
    [w_cx, st] = polyqp_mex('TOADWaypoint', args{:});
    [ok, err, msg] = cmpStruct(w_ml, w_cx, 'rel', TOL_ALG, {'StraightLine', 'Label', 'IsPassAndGo'});
    rec('TOADWaypoint', ok && isempty(st), err, sprintf('case %d: %s', q, msg));
end

A0ML = cell(1, 6);
H0ML = cell(1, 6);
for r = 5
    [A0, Q0, H0] = assembleQP(r);
    A0ML{r} = A0;
    H0ML{r} = H0;
    [A0c, Q0c, H0c, st] = polyqp_mex('assembleQP', r);
    [ok1, e1] = cmpNum(A0, A0c, 'rel', TOL_ALG);
    [ok2, e2] = cmpNum(Q0, Q0c, 'rel', TOL_ALG);
    [ok3, e3] = cmpNum(H0, H0c, 'rel', TOL_LIN);
    rec('assembleQP A0, Q0', ok1 && ok2 && isempty(st), max(e1, e2), sprintf('r=%d eA0=%.3g eQ0=%.3g', r, e1, e2));
    rec('assembleQP H0', ok3 && isempty(st), e3, sprintf('r=%d eH0=%.3g', r, e3));
end
note('assembleQP H0', sprintf('r = 5 only; log10 cond(A11) = %.1f', log10(cond(A0ML{5}(6:10, 6:10)))));

% ---------------------------------------------------------------- resolveWaypoints error IDs
fprintf('[resolveWaypoints] error-ID cases\n');
c0 = C{1};
Wok = abortWaypoints([1; 0; 0; 0; 3; -2; 10; 1; 0.5; -1; 0; 0; 0; 0; 0]);
errCases = {};
errCases{end + 1} = {Wok(1), 5};
% (the r = 2 and r = 7 InvalidOrder cases were dropped: the C++ hardcodes r = 5)
Wbad = Wok; Wbad(1).HDGRef = [0; 0; 0];   errCases{end + 1} = {Wbad, 5};
Wbad = Wok; Wbad(2).HDGRef = [0; 1; 0];   errCases{end + 1} = {Wbad, 5};
Wbad = Wok; Wbad(2).MaxVel = [1; 0; 1];   errCases{end + 1} = {Wbad, 5};
Wbad = Wok; Wbad(3).MaxAccel = [1; 1; -1]; errCases{end + 1} = {Wbad, 5};
Wbad = Wok; Wbad(3).HDGRef = [2; 0; 0];   errCases{end + 1} = {Wbad, 5};
for q = 1:numel(errCases)
    a = errCases{q};
    idML = '';
    try
        resolveWaypoints(a{1}, a{2}, LimML{1});
    catch ME
        idML = ME.identifier;
    end
    [~, st] = polyqp_mex('resolveWaypoints', a{1}, a{2}, LimML{1});
    rec('resolveWaypoints error IDs', strcmp(idML, st), NaN, sprintf('case %d: MATLAB "%s" vs C++ "%s"', q, idML, st));
    idML = '';
    try
        PolyQP(a{1}, c0, struct('R', a{2}));
    catch ME
        idML = ME.identifier;
    end
    [~, ~, ~, st] = polyqp_mex('PolyQP', a{1}, c0, struct('R', a{2}));
    rec('PolyQP error IDs', strcmp(idML, st), NaN, sprintf('case %d: MATLAB "%s" vs C++ "%s"', q, idML, st));
end

% ---------------------------------------------------------------- 2. abort cases
nAbort = 200;
abortStats = struct('veh', {}, 'n', {}, 'feasible', {}, 'outerHist', {}, 'T_total', {}, 'maxLogCond', {}, 'nOver', {}, 'feasOver', {});
maxLogCondAll = -Inf;
for v = 1:numel(C)
    c = C{v};
    Lim = LimML{v};
    fprintf('[abort] %s: %d cases\n', vehName{v}, nAbort);
    feas = 0;
    nOver = 0;
    feasOver = 0;
    outerIters = zeros(1, nAbort);
    Ttot = zeros(1, nAbort);
    for q = 1:nAbort
        X_est = zeros(15, 1);
        X_est(1:4) = [1; 0; 0; 0];
        X_est(5:7) = [-20 + 40 * rand; -20 + 40 * rand; 0.5 + 39.5 * rand];
        X_est(8:10) = [-5 + 10 * rand; -5 + 10 * rand; -6 + 10 * rand];
        X_est(14) = c.OxMass;
        X_est(15) = c.FuMass;
        tag = sprintf('%s abort #%d X=%s', vehName{v}, q, mat2str(X_est(5:10).', 4));

        Waypoints = abortWaypoints(X_est);

        % ---- PolyQP end-to-end
        idML = '';
        try
            [IG, QT, RP] = PolyQP(Waypoints, c, struct('R', 5));
        catch ME
            idML = ME.identifier;
        end
        [IGc, QTc, RPc, st] = polyqp_mex('PolyQP', Waypoints, c, struct('R', 5));
        okS = strcmp(idML, st);
        if okS && isempty(idML)
            [okR, ~, mR] = cmpStruct(RP, RPc, 'exactonly', 0, {'Feasible', 'OuterIterations'});
            [okQ, ~, mQ] = cmpStruct(QT, QTc, 'exactonly', 0, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okG, ~, mG] = sameShapes(IG, IGc);
            okS = okR && okQ && okG;
            rec('PolyQP structural / status', okS, NaN, sprintf('%s: %s %s %s', tag, mR, mQ, mG));
            [okE1, e1, m1] = cmpStruct(IG, IGc, 'abs', TOL_E2E, {});
            [okE2, e2, m2] = cmpStruct(QT, QTc, 'abs', TOL_E2E, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okE3, e3, m3] = cmpStruct(RP, RPc, 'abs', TOL_E2E, {'Feasible', 'OuterIterations'});
            rec('PolyQP end-to-end', okE1 && okE2 && okE3, max([e1, e2, e3]), sprintf('%s: IG %s | QT %s | RP %s', tag, m1, m2, m3), ...
                @() floorPolyQP(Waypoints, c, 5, IG, QT, RP));
            feas = feas + RP.Feasible;
            isOver = any(abs(X_est(8:10)) > [4; 4; 5]);
            nOver = nOver + isOver;
            feasOver = feasOver + (isOver && RP.Feasible);
            outerIters(q) = RP.OuterIterations;
            Ttot(q) = RP.T_total;
        else
            rec('PolyQP structural / status', okS, NaN, sprintf('%s: MATLAB "%s" vs C++ "%s"', tag, idML, st));
        end
        if isempty(idML)
            caseLog.abort{end + 1} = struct('veh', v, 'q', q, 'X_est', X_est, 'idML', idML, 'st', st, 'IG', IG, 'QT', QT, 'RP', RP, ...
                'IGc', IGc, 'QTc', QTc, 'RPc', RPc);
        else
            caseLog.abort{end + 1} = struct('veh', v, 'q', q, 'X_est', X_est, 'idML', idML, 'st', st, 'IG', [], 'QT', [], 'RP', [], ...
                'IGc', [], 'QTc', [], 'RPc', []);
        end

        % ---- AbortGuidance vs the MATLAB abort block (persistent semantics as explicit state)
        ABORT = 1;
        tFire = 2;
        S0 = struct('abortCalculated', false, 'tAbort', 0, 'abortStatus', 0, 'Time', [], 'States', [], 'Inputs', []);
        idML = '';
        try
            S1 = AbortGuidance(tFire, X_est, c, ABORT, S0);
        catch ME
            idML = ME.identifier;
        end
        st0 = struct('abortCalculated', false, 'tAbort', 0, 'abortStatus', 0);
        [Tm, Sm, Im, tAb, aC, ~, aSt, st] = polyqp_mex('AbortGuidance', tFire, X_est, c, st0, ABORT);
        okS = strcmp(idML, st);
        if okS && isempty(idML)
            okS = isequal(size(S1.Time), size(Tm)) && isequal(size(S1.States), size(Sm)) && isequal(size(S1.Inputs), size(Im)) ...
                && tAb == S1.tAbort && aC == S1.abortCalculated && aSt == S1.abortStatus;
            % Second call: abortCalculated persists, nothing is recomputed
            S2 = AbortGuidance(tFire + 1, X_est, c, ABORT, S1);
            [Tm2, ~, ~, tAb2, aC2, ~, aSt2, st2] = polyqp_mex('AbortGuidance', tFire + 1, X_est, c, struct('abortCalculated', aC, 'tAbort', tAb, 'abortStatus', aSt), ABORT);
            okS = okS && isempty(st2) && isempty(Tm2) && tAb2 == S2.tAbort && aC2 == S2.abortCalculated && aSt2 == S2.abortStatus;
            % Default trigger (MATLAB ABORT = 0): never fires
            S3 = AbortGuidance(tFire, X_est, c, 0, S0);
            [Tm3, ~, ~, tAb3, aC3, ~, ~, st3] = polyqp_mex('AbortGuidance', tFire, X_est, c, st0);
            okS = okS && isempty(st3) && isempty(Tm3) && tAb3 == S3.tAbort && aC3 == S3.abortCalculated;
            rec('AbortGuidance structural', okS, NaN, sprintf('%s: shapes/tAbort/abortCalculated/abortStatus/persistence mismatch', tag));
            [ok1, e1] = cmpNum(S1.Time, Tm, 'abs', TOL_E2E);
            [ok2, e2] = cmpNum(S1.States, Sm, 'abs', TOL_E2E);
            [ok3, e3] = cmpNum(S1.Inputs, Im, 'abs', TOL_E2E);
            rec('AbortGuidance end-to-end', ok1 && ok2 && ok3, max([e1, e2, e3]), sprintf('%s: eTime=%.3g eStates=%.3g eInputs=%.3g', tag, e1, e2, e3), ...
                @() max(floorAbort(tFire, X_est, c, ABORT, S0, S1), exactFinalSolveFloor(Waypoints, c, 5, QT)));
        else
            rec('AbortGuidance structural', okS, NaN, sprintf('%s: MATLAB "%s" vs C++ "%s"', tag, idML, st));
        end

        % ---- per-function checks on inputs from this run
        if isempty(idML)
            lc = perFunctionChecks(Waypoints, c, Lim, QT, 5, tag, mxF);
            maxLogCondAll = max(maxLogCondAll, lc);
        end
    end
    abortStats(end + 1) = struct('veh', vehName{v}, 'n', nAbort, 'feasible', feas, 'outerHist', outerIters, 'T_total', Ttot, 'maxLogCond', maxLogCondAll, 'nOver', nOver, 'feasOver', feasOver); %#ok<AGROW>
end

% ---- Tilted, rotating starts (FIX-F03 code path); separate RNG stream
fprintf('[abort] tilted / depleted starts\n');
trs = RandStream('mt19937ar', 'Seed', 4242);
nTilt = 40;
tiltStats = zeros(numel(C), 2);
for v = 1:numel(C)
    c = C{v};
    for q = 1:nTilt
        X_est = tiltedAbortState(trs, c);
        OptsQ = struct('R', 5, 'InitOxMass', X_est(14), 'InitFuMass', X_est(15));
        LimQ = buildVehicleLimits(c, OptsQ);
        Waypoints = abortWaypoints(X_est);
        tag = sprintf('%s tilted #%d tilt=%.1fdeg |w|=%.2f m=(%.3g,%.3g) X=%s', vehName{v}, q, rad2deg(acos(1 - 2 * (X_est(2)^2 + X_est(3)^2))), ...
            norm(X_est(11:13)), X_est(14), X_est(15), mat2str(X_est(5:10).', 4));
        idML = '';
        try
            [IG, QT, RP] = PolyQP(Waypoints, c, OptsQ);
        catch ME
            idML = ME.identifier;
        end
        [IGc, QTc, RPc, st] = polyqp_mex('PolyQP', Waypoints, c, OptsQ);
        okS = strcmp(idML, st);
        if okS && isempty(idML)
            [okR, ~, mR] = cmpStruct(RP, RPc, 'exactonly', 0, {'Feasible', 'OuterIterations'});
            [okQ, ~, mQ] = cmpStruct(QT, QTc, 'exactonly', 0, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okG, ~, mG] = sameShapes(IG, IGc);
            rec('PolyQP structural / status', okR && okQ && okG, NaN, sprintf('%s: %s %s %s', tag, mR, mQ, mG));
            [okE1, e1, m1] = cmpStruct(IG, IGc, 'abs', TOL_E2E, {});
            [okE2, e2, m2] = cmpStruct(QT, QTc, 'abs', TOL_E2E, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okE3, e3, m3] = cmpStruct(RP, RPc, 'abs', TOL_E2E, {'Feasible', 'OuterIterations'});
            rec('PolyQP end-to-end', okE1 && okE2 && okE3, max([e1, e2, e3]), sprintf('%s: IG %s | QT %s | RP %s', tag, m1, m2, m3), ...
                @() floorPolyQP(Waypoints, c, 5, IG, QT, RP, OptsQ));
            tiltStats(v, :) = tiltStats(v, :) + [1, RP.Feasible];
            lc = perFunctionChecks(Waypoints, c, LimQ, QT, 5, tag, mxF);
            maxLogCondAll = max(maxLogCondAll, lc);
        else
            rec('PolyQP structural / status', okS, NaN, sprintf('%s: MATLAB "%s" vs C++ "%s"', tag, idML, st));
        end
        S0 = struct('abortCalculated', false, 'tAbort', 0, 'abortStatus', 0, 'Time', [], 'States', [], 'Inputs', []);
        S1 = AbortGuidance(2, X_est, c, 1, S0);
        [Tm, Sm, Im, tAb, aC, ~, aSt] = polyqp_mex('AbortGuidance', 2, X_est, c, struct('abortCalculated', false, 'tAbort', 0, 'abortStatus', 0), 1);
        okS = isequal(size(S1.Time), size(Tm)) && isequal(size(S1.States), size(Sm)) && isequal(size(S1.Inputs), size(Im)) ...
            && tAb == S1.tAbort && aC == S1.abortCalculated && aSt == S1.abortStatus;
        rec('AbortGuidance structural', okS, NaN, sprintf('%s: shapes/tAbort/abortCalculated/abortStatus mismatch', tag));
        if okS
            [ok1, e1] = cmpNum(S1.Time, Tm, 'abs', TOL_E2E);
            [ok2, e2] = cmpNum(S1.States, Sm, 'abs', TOL_E2E);
            [ok3, e3] = cmpNum(S1.Inputs, Im, 'abs', TOL_E2E);
            rec('AbortGuidance end-to-end', ok1 && ok2 && ok3, max([e1, e2, e3]), sprintf('%s: eTime=%.3g eStates=%.3g eInputs=%.3g', tag, e1, e2, e3), ...
                @() floorAbort(2, X_est, c, 1, S0, S1));
        end
    end
end

% ---------------------------------------------------------------- 4. missions
fprintf('[missions]\n');
missions = {};
missions{end + 1} = {'Circle (BuildTrajectory.m)', circleMissionWaypoints()};
missions{end + 1} = {'Synthetic: hold + cruise split', cruiseMissionWaypoints()};
missions{end + 1} = {'Synthetic: free pass-through corners', cornerMissionWaypoints()};
missionLog = {};
for mi = 1:numel(missions)
    for v = 1:numel(C)
        Waypoints = missions{mi}{2};
        tag = sprintf('%s / %s', missions{mi}{1}, vehName{v});
        idML = '';
        try
            [IG, QT, RP] = PolyQP(Waypoints, C{v}, struct('R', 5));
        catch ME
            idML = ME.identifier;
        end
        [IGc, QTc, RPc, st] = polyqp_mex_wide('PolyQP', Waypoints, C{v}, struct('R', 5));
        okS = strcmp(idML, st);
        if okS && isempty(idML)
            [okR, ~, mR] = cmpStruct(RP, RPc, 'exactonly', 0, {'Feasible', 'OuterIterations'});
            [okQ, ~, mQ] = cmpStruct(QT, QTc, 'exactonly', 0, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okG, ~, mG] = sameShapes(IG, IGc);
            rec('Missions structural / status', okR && okQ && okG, NaN, sprintf('%s: %s %s %s', tag, mR, mQ, mG));
            [okE1, e1, m1] = cmpStruct(IG, IGc, 'abs', TOL_E2E, {});
            [okE2, e2, m2] = cmpStruct(QT, QTc, 'abs', TOL_E2E, {'R', 'Count', 'IsHold', 'IsCruise'});
            [okE3, e3, m3] = cmpStruct(RP, RPc, 'abs', TOL_E2E, {'Feasible', 'OuterIterations'});
            rec('Missions end-to-end', okE1 && okE2 && okE3, max([e1, e2, e3]), sprintf('%s: IG %s | QT %s | RP %s', tag, m1, m2, m3), ...
                @() floorPolyQP(Waypoints, C{v}, 5, IG, QT, RP));
            Plan = resolveWaypoints(Waypoints, 5, LimML{v});
            missionLog{end + 1} = sprintf('| %s | %d | %d | %d | %d | %d | %d | %.3f | %d | %.3g |', tag, Plan.SegCount, Plan.LegCount, ...
                nnz(Plan.IsHold), nnz(Plan.IsCruise), RP.Feasible, RP.OuterIterations, RP.T_total, RPc.OuterIterations, max([e1, e2, e3])); %#ok<AGROW>
            lc = perFunctionChecks(Waypoints, C{v}, LimML{v}, QT, 5, tag, mxW);
            maxLogCondAll = max(maxLogCondAll, lc);
        else
            rec('Missions structural / status', okS, NaN, sprintf('%s: MATLAB "%s" vs C++ "%s"', tag, idML, st));
            missionLog{end + 1} = sprintf('| %s | error: MATLAB "%s", C++ "%s" | | | | | | | | |', tag, idML, st); %#ok<AGROW>
        end
    end
end

% ---------------------------------------------------------------- 3. solveLeg fuzz (500)
fprintf('[solveLeg] 500 fuzzed cases\n');
nFuzz = 500;
fuzzKinds = zeros(1, 4);
fuzzLogCond = -Inf;
for q = 1:nFuzz
    if q <= 300
        kind = 1;   % random NaN patterns
        r = 5;
        K = randi([2, 5]);
        FixLeg = 3 * randn(3, r, K);
        p = 0.2 + 0.5 * rand;
        FixLeg(rand(3, r, K) < p) = NaN;
        SegTimes = 0.2 + 9.8 * rand(1, K - 1);
    elseif q <= 380
        kind = 2;   % fewer than r fixed values on at least one axis
        r = 5;
        K = randi([2, 4]);
        FixLeg = 3 * randn(3, r, K);
        d = randi(3);
        slice = FixLeg(d, :, :);
        keep = randperm(numel(slice), randi([0, r - 1]));
        tmp = nan(size(slice));
        tmp(keep) = slice(keep);
        FixLeg(d, :, :) = tmp;
        SegTimes = 0.2 + 9.8 * rand(1, K - 1);
    elseif q <= 460
        kind = 3;   % near-singular: abort-like leg with unbalanced segment times
        r = 5;
        K = 3;
        FixLeg = abortLikeFixLeg(r);
        T1 = 4 + 6 * rand;
        T2 = 0.2 + 0.8 * rand;
        if rand < 0.5
            SegTimes = [T1, T2];
        else
            SegTimes = [T2, T1];
        end
    else
        kind = 4;   % fixed boundary-pinned legs with free interior knots
        r = 5;
        K = randi([3, 6]);
        FixLeg = nan(3, r, K);
        FixLeg(:, 1, :) = 5 * randn(3, 1, K);
        FixLeg(:, 2:r, 1) = 0;
        FixLeg(:, 2:r, K) = 0;
        SegTimes = 0.2 + 9.8 * rand(1, K - 1);
    end
    if q == 461
        kind = 3; r = 5; K = 3; FixLeg = abortLikeFixLeg(r); SegTimes = [10, 0.6];
    elseif q == 462
        kind = 3; r = 5; K = 3; FixLeg = abortLikeFixLeg(r); SegTimes = [6, 2];
    end
    fuzzKinds(kind) = fuzzKinds(kind) + 1;
    lc = checkSolveLeg(FixLeg, SegTimes, A0ML{r}, H0ML{r}, mxW, 'solveLeg (fuzzed)', sprintf('fuzz #%d kind %d r=%d K=%d T=%s', q, kind, r, K, mat2str(SegTimes, 4)));
    fuzzLogCond = max(fuzzLogCond, lc);
end
note('solveLeg (fuzzed)', sprintf('300 random NaN patterns, 80 with < r fixed, 82 near-singular (unbalanced T incl. (10,0.6) and (6,2)), 38 boundary-pinned; max log10 cond(R_PP_eq) among solved cases = %.1f; MATLAB statuses: %d ok, %d IllPosedQP, %d SingularRPP', ...
    fuzzLogCond, statusCounts.ok, statusCounts.IllPosedQP, statusCounts.SingularRPP));
note('solveLeg (from solveAllLegs)', sprintf('FixP built as solveAllLegs pass 1, at initial and final SegT; max log10 cond(R_PP_eq) = %.1f', maxLogCondAll));

% ---------------------------------------------------------------- allocateTimes fuzz
fprintf('[allocateTimes] fuzzed multi-segment legs\n');
nAT = 150;
for q = 1:nAT
    r = 5;
    if q <= 10
        K = 2;
    else
        K = randi([3, 6]);
    end
    FixLeg = nan(3, r, K);
    FixLeg(:, 1, :) = 6 * randn(3, 1, K);
    FixLeg(:, 2:r, 1) = 0;
    FixLeg(:, 2:r, K) = 0;
    if mod(q, 4) == 0 && K > 2
        FixLeg(3, 2, 2:(K - 1)) = 0;   % one axis with pinned interior velocity, others free
    end
    T0 = 0.3 + 3.7 * rand(1, K - 1);
    checkAllocateTimes(FixLeg, T0, A0ML{r}, H0ML{r}, 0.2, mxW, sprintf('allocTimes #%d r=%d K=%d T0=%s', q, r, K, mat2str(T0, 4)));
end
note('allocateTimes', sprintf('%d fuzzed legs (10 single-segment for the M < 2 return), pinned rest endpoints and free interior derivatives; %d legs where the line search moved the times; also run on every qualifying mission leg', nAT, atMoved));

% ---------------------------------------------------------------- STATIC: identity-padded R_PP solve (Phase 4)
% Every solveLeg case of this run that fits the flight capacities, through polyqp_mex (fixed POLYQP_MAX_FREE x POLYQP_MAX_FREE
% system, identity-padded) and polyqp_mex_p3 (the Phase-3 build: variable-size, fixed-max storage, the system solved at its
% own size). Same status required; outputs within the linear-solve tier or 10 x cond(R_PP_eq) * eps. Separately, the rcond the
% padded system gives against the rcond of its unpadded leading nP x nP block (fixed-size Eigen in the MEX), and whether the
% 1e-13 guard decides the same.
fprintf('[padding] identity-padded solve vs variable-size solve\n');
padStats = struct('n', 0, 'beyond', 0, 'statusDiff', 0, 'rcN', 0, 'rcDecisionDiff', 0, 'maxRcRel', 0, 'maxPadEstRel', 0, 'nPadEstDiff', 0, 'padDecisionDiff', 0);
for i = 1:numel(caseLog.solveLeg)
    e = caseLog.solveLeg{i};
    if size(e.FixLeg, 3) > POLYQP_FLIGHT_KNOTS
        padStats.beyond = padStats.beyond + 1;
        continue;
    end
    [Cf, Df, Jf, stf] = polyqp_mex('solveLeg', e.FixLeg, e.SegTimes, A0ML{5}, H0ML{5});
    if strcmp(stf, 'PolyQP:CapacityExceeded')
        padStats.beyond = padStats.beyond + 1;
        continue;
    end
    padStats.n = padStats.n + 1;
    [rcG, nFree, rcU, rcP] = polyqp_mex('solveLegRcond', e.FixLeg, e.SegTimes, A0ML{5}, H0ML{5});
    for d = 1:3
        if nFree(d) > 0 && isfinite(rcU(d))
            padStats.rcN = padStats.rcN + 1;
            relRc = abs(rcG(d) - rcU(d)) / max(abs(rcU(d)), realmin);
            padStats.maxRcRel = max(padStats.maxRcRel, relRc);
            padStats.rcDecisionDiff = padStats.rcDecisionDiff + ((rcG(d) >= 1e-13) ~= (rcU(d) >= 1e-13));
            rec('rcond guard on unpadded block', relRc <= 1e-12 && ((rcG(d) >= 1e-13) == (rcU(d) >= 1e-13)), relRc, sprintf('%s axis %d nFree=%d rc guard %.3g, unpadded %.3g', e.tag, d, nFree(d), rcG(d), rcU(d)));
            % informational: Eigen's estimate on the padded matrix itself (same exact 1-norm condition, different estimate)
            relPad = abs(rcP(d) - rcU(d)) / max(abs(rcU(d)), realmin);
            padStats.maxPadEstRel = max(padStats.maxPadEstRel, relPad);
            padStats.nPadEstDiff = padStats.nPadEstDiff + (relPad > 1e-12);
            padStats.padDecisionDiff = padStats.padDecisionDiff + ((rcP(d) >= 1e-13) ~= (rcU(d) >= 1e-13));
        end
    end
    if haveP3
        [Cp, Dp, Jp, stp] = polyqp_mex_p3('solveLeg', e.FixLeg, e.SegTimes, A0ML{5}, H0ML{5});
        if ~strcmp(stp, stf)
            padStats.statusDiff = padStats.statusDiff + 1;
            rec('Padded solve vs variable-size solve', false, NaN, sprintf('%s: status padded "%s" vs variable-size "%s"', e.tag, stf, stp));
            continue;
        end
        if isempty(stf)
            [ok1, e1] = cmpNum(Cp, Cf, 'rel', TOL_LIN);
            [ok2, e2] = cmpNum(Dp, Df, 'rel', TOL_LIN);
            [ok3, e3] = cmpNum(Jp, Jf, 'rel', TOL_LIN);
            lcp = logCondRPP(e.FixLeg, e.SegTimes, H0ML{5});
            rec('Padded solve vs variable-size solve', ok1 && ok2 && ok3, max([e1, e2, e3]), sprintf('%s: eCoef=%.3g eDeriv=%.3g eCost=%.3g log10cond=%.1f', e.tag, e1, e2, e3, lcp), ...
                @() 10 ^ lcp * eps);
        else
            rec('Padded solve vs variable-size solve', true, 0, '');
        end
    end
end
note('Padded solve vs variable-size solve', sprintf('%d solveLeg cases within the flight capacities (%d beyond them skipped); %s; status mismatches %d', ...
    padStats.n, padStats.beyond, ternary(haveP3, 'against polyqp_mex_p3 (Phase 3 build)', 'polyqp_mex_p3 NOT available: outputs not compared'), padStats.statusDiff));
note('rcond guard on unpadded block', sprintf('%d axis solves; guard vs independent unpadded rcond: max rel diff %.2e, decisions that differ %d. Informational: Eigen''s estimate on the padded matrix differs from the unpadded one in %d solves (max rel %.2f; guard decisions it would change: %d), which is why the guard uses the unpadded block', ...
    padStats.rcN, padStats.maxRcRel, padStats.rcDecisionDiff, padStats.nPadEstDiff, padStats.maxPadEstRel, padStats.padDecisionDiff));

% ---------------------------------------------------------------- self-review (step 9): comments + flight-source scan
selfReview = selfReviewRows(toadDir, fileparts(here));
if any(~[selfReview.ok])
    addRow('Self-review (comments, flight sources)', 'exact', 0);
    for i = 1:numel(selfReview)
        rec('Self-review (comments, flight sources)', selfReview(i).ok, NaN, sprintf('%s: %s', selfReview(i).cpp, selfReview(i).msg));
    end
end

% ---------------------------------------------------------------- 7. report
elapsed = toc(tStart);
if ~isempty(runTag)
    if ~isfolder(fullfile(here, 'results'))
        mkdir(fullfile(here, 'results'));
    end
    save(fullfile(here, 'results', [runTag, '_cases.mat']), 'caseLog', 'rows', 'abortStats', '-v7');
end
% ---------------------------------------------------------------- 8. abort-fix regressions (one section per FIX, ABORT_FIXES.md)
fprintf('[fix regressions]\n');
fixR = fix_regressions();
anyFail = any([rows.fails] > 0) || any(~[fixR.pass]);   % unexplained failures only (strict misses within the roundoff floor are reported separately)
writeReport(fullfile(here, 'verify_polyqp_report.md'));

fprintf('\n%-36s %7s %9s %6s %11s %9s\n', 'Check', 'Cases', 'Roundoff', 'Fails', 'Max err', 'err/floor');
for i = 1:numel(rows)
    fprintf('%-36s %7d %9d %6d %11.3g %9.2f\n', rows(i).name, rows(i).cases, rows(i).roundoff, rows(i).fails, rows(i).maxErr, rows(i).maxRatio);
end
fprintf('Report: %s (%.0f s)\n', fullfile(here, 'verify_polyqp_report.md'), elapsed);
if anyFail
    for i = 1:numel(rows)
        for j = 1:numel(rows(i).failDetail)
            fprintf('  FAIL %s: %s\n', rows(i).name, rows(i).failDetail{j});
        end
    end
    error('verify_polyqp:failed', 'One or more double-precision checks failed (see report).');
end
fprintf('ALL DOUBLE-PRECISION CHECKS PASSED\n');

% ================================================================ nested helpers
    function lcMax = perFunctionChecks(Waypoints, c, Lim, QT, r, tag, mx)
        % mx: the MEX build for this case (flight for aborts, wide for missions); grid-sized I/O always uses the wide build
        % resolveWaypoints
        Plan = resolveWaypoints(Waypoints, r, Lim);
        [Pc, st] = mx('resolveWaypoints', Waypoints, r, Lim);
        structF = {'R', 'KnotCount', 'SegCount', 'KnotIsRest', 'IsHold', 'IsCruise', 'LegCount', 'LegStartSeg', 'LegEndSeg'};
        okS = isempty(st);
        if okS
            [okS, ~, mS] = cmpStruct(Plan, Pc, 'exactonly', 0, structF);
            [okN, eN, mN] = cmpStruct(Plan, Pc, 'rel', TOL_ALG, structF);
            rec('resolveWaypoints numeric', okN, eN, sprintf('%s: %s', tag, mN));
        else
            mS = st;
        end
        rec('resolveWaypoints structural', okS, NaN, sprintf('%s: %s', tag, mS));

        % solveLeg / allocateTimes with inputs reconstructed as solveAllLegs builds them (pass 1)
        lcMax = -Inf;
        K = Plan.KnotCount;
        FixP = Plan.FixCg(:, 1:r, 1:K);
        for k = 1:K
            hasB3 = all(~isnan(Plan.B3BC(:, k)));
            hasB3Dot = all(~isnan(Plan.B3DotBC(:, k)));
            for d = 1:3
                if ~isnan(FixP(d, 1, k))
                    if hasB3
                        FixP(d, 1, k) = Plan.FixCg(d, 1, k) + Lim.Rho * Plan.B3BC(d, k);
                    elseif d == 3
                        FixP(d, 1, k) = Plan.FixCg(d, 1, k) + Lim.Rho;
                    else
                        FixP(d, 1, k) = Plan.FixCg(d, 1, k);
                    end
                end
                if hasB3Dot && ~isnan(FixP(d, 2, k))
                    FixP(d, 2, k) = Plan.FixCg(d, 2, k) + Lim.Rho * Plan.B3DotBC(d, k);
                end
            end
        end
        for leg = 1:Plan.LegCount
            m1 = Plan.LegStartSeg(leg);
            m2 = Plan.LegEndSeg(leg);
            FixLeg = FixP(:, 1:r, m1:(m2 + 1));
            lcMax = max(lcMax, checkSolveLeg(FixLeg, Plan.SegT(m1:m2), A0ML{r}, H0ML{r}, mx, 'solveLeg (from solveAllLegs)', sprintf('%s leg %d initial T', tag, leg)));
            lcMax = max(lcMax, checkSolveLeg(FixLeg, QT.T(m1:m2), A0ML{r}, H0ML{r}, mx, 'solveLeg (from solveAllLegs)', sprintf('%s leg %d final T', tag, leg)));
            if m2 > m1
                midVelBins = FixP(:, 2, (m1 + 1):m2);
                if all(any(isnan(midVelBins), 1))
                    checkAllocateTimes(FixLeg, Plan.SegT(m1:m2), A0ML{r}, H0ML{r}, Lim.MinSegTime, mx, sprintf('%s leg %d', tag, leg));
                end
            end
        end

        % PolyEval: k = 0..5, the k >= nc branch, the nargin default, and query times exactly at knots
        M = QT.Count;
        T_total = QT.T_total;
        tq = [linspace(0, T_total, 57), QT.KeyTimes(1:(M + 1)), -0.5, T_total + 0.5, 0.5 * (QT.KeyTimes(1:M) + QT.KeyTimes(2:(M + 1)))];
        for k = [0:5, 2 * r]
            Out = PolyEval(QT, tq, k);
            [Oc, st] = mx('PolyEval', QT, tq, k);
            [ok, err] = cmpNum(Out, Oc, 'rel', TOL_ALG);
            rec('PolyEval', ok && isempty(st), err, sprintf('%s k=%d err=%.3g', tag, k, err), @() floorPolyEval(QT, tq, k, Out));
        end
        Out = PolyEval(QT, tq);
        [Oc, st] = mx('PolyEval', QT, tq);
        [ok, err] = cmpNum(Out, Oc, 'rel', TOL_ALG);
        rec('PolyEval', ok && isempty(st), err, sprintf('%s default k err=%.3g', tag, err), @() floorPolyEval(QT, tq, 0, Out));

        % flatnessInvert / allocateControls / evaluateConstraints / sampleOutputs on MATLAB's QPTraj
        PlanE = Plan;
        PlanE.SegT(1:M) = QT.T(1:M);
        dt_eval = 0.5 * Lim.DtSample;
        nGrid = min(3999, max(41, ceil(T_total / dt_eval) + 1));
        t_grid = linspace(0, T_total, nGrid);
        [X, Wr, AccelCg] = flatnessInvert(QT, t_grid, PlanE.HDGRef, c, Lim);
        [Xc, Wrc, Ac, st] = polyqp_mex_wide('flatnessInvert', QT, t_grid, PlanE.HDGRef, c, Lim);   % grid I/O (up to 3999 samples)
        [ok1, e1] = cmpNum(X, Xc, 'relabs', TOL_NL);
        [ok2, e2, m2] = cmpStruct(Wr, Wrc, 'relabs', TOL_NL, {});
        [ok3, e3] = cmpNum(AccelCg, Ac, 'relabs', TOL_NL);
        rec('flatnessInvert', ok1 && ok2 && ok3 && isempty(st), max([e1, e2, e3]), sprintf('%s eX=%.3g W:%s eA=%.3g', tag, e1, m2, e3));

        U = allocateControls(Wr, Lim);
        [Uc, st] = polyqp_mex_wide('allocateControls', Wr, Lim);
        [ok, err] = cmpNum(U, Uc, 'relabs', TOL_NL);
        rec('allocateControls', ok && isempty(st), err, sprintf('%s eU=%.3g', tag, err));

        [EO, Fe] = evaluateConstraints(QT, PlanE, c, Lim);
        [EOc, Fec, st] = mx('evaluateConstraints', QT, PlanE, c, Lim);
        rec('evaluateConstraints structural', isempty(st) && isequal(Fe, Fec) && isequal(size(EO.SegUtil), size(EOc.SegUtil)), NaN, ...
            sprintf('%s Feasible %d vs %d', tag, Fe, Fec));
        [ok, err, msg] = cmpStruct(EO, EOc, 'relabs', TOL_NL, {});
        rec('evaluateConstraints numeric', ok && isempty(st), err, sprintf('%s %s', tag, msg));

        IG = sampleOutputs(QT, PlanE, c, Lim);
        [IGc, st] = mx('sampleOutputs', QT, PlanE, c, Lim);
        [ok, err, msg] = cmpStruct(IG, IGc, 'relabs', TOL_NL, {});
        rec('sampleOutputs', ok && isempty(st), err, sprintf('%s %s', tag, msg));
    end

    function lc = checkSolveLeg(FixLeg, SegTimes, A0, H0, mx, rowName, tag)
        idML = '';
        try
            [Cml, Dml, Jml] = solveLeg(FixLeg, SegTimes, A0, H0);
        catch ME
            idML = ME.identifier;
        end
        [Ccx, Dcx, Jcx, st] = mx('solveLeg', FixLeg, SegTimes, A0, H0);
        if isempty(idML)
            caseLog.solveLeg{end + 1} = struct('row', rowName, 'tag', tag, 'FixLeg', FixLeg, 'SegTimes', SegTimes, 'idML', idML, 'st', st, ...
                'Cml', Cml, 'Dml', Dml, 'Jml', Jml, 'Ccx', Ccx, 'Dcx', Dcx, 'Jcx', Jcx);
        else
            caseLog.solveLeg{end + 1} = struct('row', rowName, 'tag', tag, 'FixLeg', FixLeg, 'SegTimes', SegTimes, 'idML', idML, 'st', st, ...
                'Cml', [], 'Dml', [], 'Jml', [], 'Ccx', Ccx, 'Dcx', Dcx, 'Jcx', Jcx);
        end
        if strcmp(rowName, 'solveLeg (fuzzed)')
            if isempty(idML)
                statusCounts.ok = statusCounts.ok + 1;
            else
                f = matlab.lang.makeValidName(strrep(idML, 'PolyQP:', ''));
                if ~isfield(statusCounts, f)
                    statusCounts.(f) = 0;
                end
                statusCounts.(f) = statusCounts.(f) + 1;
            end
        end
        lc = NaN;
        if isempty(idML)
            lc = logCondRPP(FixLeg, SegTimes, H0);
        end
        rec('solveLeg error IDs', strcmp(idML, st), NaN, sprintf('%s: MATLAB "%s" vs C++ "%s" (log10 cond R_PP_eq = %.1f, MATLAB min rcond = %.3g)', ...
            tag, idML, st, logCondRPP(FixLeg, SegTimes, H0), minRcondRPP(FixLeg, SegTimes, H0)));
        if isempty(idML) && isempty(st)
            [ok1, e1] = cmpNum(Cml, Ccx, 'rel', TOL_LIN);
            [ok2, e2] = cmpNum(Dml, Dcx, 'rel', TOL_LIN);
            [ok3, e3] = cmpNum(Jml, Jcx, 'rel', TOL_LIN);
            okAll = ok1 && ok2 && ok3;
            floorVal = NaN;
            if ~okAll
                % Roundoff floor = MATLAB's own error against the exact-arithmetic solve (same metric)
                [Cex, Dex, Jex] = exactSolveLeg(FixLeg, SegTimes, A0, H0);
                [~, f1] = cmpNum(Cex, Cml, 'rel', Inf);
                [~, f2] = cmpNum(Dex, Dml, 'rel', Inf);
                [~, f3] = cmpNum(Jex, Jml, 'rel', Inf);
                floorVal = max([f1, f2, f3]);
                if ~isnan(lc)
                    % Textbook forward-error bound of a backward-stable solve, cond(R_PP_eq) * eps (no n factor)
                    floorVal = max(floorVal, 10 ^ lc * eps / ROUNDOFF_FACTOR);
                end
                % Cost only: the rounding bound of evaluating the quadratic form d' R d, nv * eps * (|d|' |R| |d|) / |Cost|. A near-
                % optimal cost can be orders of magnitude smaller than its terms (kappa_cost up to 6e10 on fuzzed legs), so the
                % cost's own floor is the larger of floorVal and this bound; Coef and Deriv keep floorVal.
                costFloor = max(floorVal, quadFormFloor(FixLeg, SegTimes, H0, Dml, Jml) / ROUNDOFF_FACTOR);
                e3 = e3 * floorVal / costFloor;   % the cost error expressed against floorVal
                if numel(exactRef) < 25
                    [~, g1] = cmpNum(Cex, Ccx, 'rel', Inf);
                    [~, g2] = cmpNum(Dex, Dcx, 'rel', Inf);
                    [~, g3] = cmpNum(Jex, Jcx, 'rel', Inf);
                    exactRef{end + 1} = sprintf('| %s | %.1f | %.2e | %.2e | %.2e | %.2e |', tag, lc, max([f1, f2, f3]), max([g1, g2, g3]), max([e1, e2, e3]), 10 ^ lc * eps);
                end
            end
            rec(rowName, okAll, max([e1, e2, e3]), sprintf('%s: eCoef=%.3g eDeriv=%.3g eCost=%.3g log10cond=%.1f', tag, e1, e2, e3, lc), ...
                @() floorVal);
        end
    end

    function checkAllocateTimes(FixLeg, T0, A0, H0, MinSegTime, mx, tag)
        idML = '';
        try
            [Tml, Jml] = allocateTimes(FixLeg, T0, A0, H0, MinSegTime);
        catch ME
            idML = ME.identifier;
        end
        [Tcx, Jcx, st] = mx('allocateTimes', FixLeg, T0, A0, H0, MinSegTime);
        if ~strcmp(idML, st)
            rec('allocateTimes', false, NaN, sprintf('%s: MATLAB "%s" vs C++ "%s"', tag, idML, st));
            return;
        end
        if ~isempty(idML)
            rec('allocateTimes', true, NaN, '');
            return;
        end
        if max(abs(Tml - T0)) > 1e-12
            atMoved = atMoved + 1;
        end
        [ok1, e1] = cmpNum(Tml, Tcx, 'rel', TOL_LIN);
        [ok2, e2] = cmpNum(Jml, Jcx, 'rel', TOL_LIN);
        rec('allocateTimes', ok1 && ok2, max(e1, e2), sprintf('%s: eT=%.3g eCost=%.3g', tag, e1, e2), ...
            @() floorAllocateTimes(FixLeg, T0, A0, H0, MinSegTime, Tml, Jml));
    end

    function writeReport(path)
        fid = fopen(path, 'w');
        closer = onCleanup(@() fclose(fid));
        fprintf(fid, '# PolyQP C++ transpile: verification report\n\n');
        fprintf(fid, 'Generated by `verify_polyqp.m` (seed `rng(42)`) on %s, MATLAB %s, %.0f s.\n\n', ...
            char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm')), version, elapsed);
        if anyFail
            fprintf(fid, '**Result: FAIL**: at least one double-precision check failed.\n\n');
        else
            fprintf(fid, '**Result: PASS**: every double-precision check passed.\n\n');
        end
        fprintf(fid, 'Relative error is `max|a - b| / max(max|a|, 1e-300)`; NaN and Inf positions must match and array sizes must be identical. ');
        fprintf(fid, '"rel or abs" passes if either is within tolerance. Exact rows report no error magnitude.\n\n');
        fprintf(fid, '## Per-function results (double precision, `polyqp_mex`)\n\n');
        fprintf(fid, '| Check | Criterion | Tol | Cases | Pass (strict) | Strict miss, roundoff-explained | Fail | Max error | Max err / floor | Notes |\n|---|---|---|---|---|---|---|---|---|---|\n');
        for i = 1:numel(rows)
            if strcmp(rows(i).kind, 'exact')
                es = '-';
                ts = 'exact';
            else
                es = sprintf('%.2e', rows(i).maxErr);
                ts = sprintf('%.0e', rows(i).tol);
            end
            if rows(i).roundoff > 0
                rs = sprintf('%.2f', rows(i).maxRatio);
            else
                rs = '-';
            end
            fprintf(fid, '| %s | %s | %s | %d | %d | %d | %d | %s | %s | %s |\n', rows(i).name, rows(i).kind, ts, rows(i).cases, ...
                rows(i).cases - rows(i).fails - rows(i).roundoff, rows(i).roundoff, rows(i).fails, es, rs, rows(i).notes);
        end
        fprintf(fid, '\n### How strict misses are classified\n\n');
        fprintf(fid, 'The tolerances above are applied as written in the transpile spec and are not loosened. A case that misses its tolerance is classed as *roundoff-explained* only if the C++-vs-MATLAB difference is at most %d x the case''s **roundoff floor**, a measurement of MATLAB''s own rounding error for that case:\n\n', ROUNDOFF_FACTOR);
        fprintf(fid, '- `solveLeg`: the C++-vs-MATLAB difference must be at most max(%d x MATLAB''s own error against the same solve in 60-digit `vpa` arithmetic, cond(R_PP_eq) * eps). The second term is the textbook first-order forward-error bound of a backward-stable (Cholesky/LU) solve, without its dimension factor; it covers cases where MATLAB''s own rounding happens to land unusually close to the exact answer. The `Cost` output alone may also use the rounding bound of its quadratic form, nv * eps * (|d|''|R||d|) / |Cost| (Phase 4: a near-optimal cost can be ~1e10 times smaller than its terms on fuzzed legs).\n', ROUNDOFF_FACTOR);
        fprintf(fid, '- `PolyQP`, `AbortGuidance`, missions: the larger of (a) MATLAB''s abs error on its final QP solve (every leg rebuilt from `QPTraj.Fix` with the `Plan.FixCg` free pattern at `QPTraj.T`, solved in `vpa`) and (b) the largest change in MATLAB''s outputs when the waypoint positions and velocities are multiplied by `1 + 2*eps*randn` (3 draws, separate RNG stream; draws that change `OuterIterations` or the sample count are discarded). `AbortGuidance` uses the floor of the identical `PolyQP` call.\n');
        fprintf(fid, '- `allocateTimes`: (b) applied to `FixLeg` and the initial times. `PolyEval`: (b) applied to `QPTraj.Coef` only (query times exactly at knots must not move).\n\n');
        fprintf(fid, 'Input perturbation alone underestimates the floor of the linear solves: it perturbs `R_PP` only along structured directions, while Cholesky''s backward error is unstructured, so its forward error is set by cond(R_PP_eq). The `vpa` table below shows MATLAB and the C++ are equally far from the exact solution. Any other miss counts as a failure and makes the script exit with an error.\n\n');
        fprintf(fid, 'Total strict misses: %d of %d numeric cases, all roundoff-explained: %s.\n\n', sum([rows.roundoff]) + sum([rows(~strcmp({rows.kind}, 'exact')).fails]), ...
            sum([rows(~strcmp({rows.kind}, 'exact')).cases]), mat2str(~anyFail));
        if ~isempty(exactRef)
            fprintf(fid, '### Exact-arithmetic reference for solveLeg strict misses (first %d)\n\n', numel(exactRef));
            fprintf(fid, 'The same solve (`Deriv`, `CoefPhys`, `Cost`) in 60-digit `vpa` arithmetic from the identical double inputs. Errors are the max over the three outputs of the relative error.\n\n');
            fprintf(fid, '| Case | log10 cond(R_PP_eq) | MATLAB vs exact | C++ vs exact | C++ vs MATLAB | cond * eps |\n|---|---|---|---|---|---|\n');
            fprintf(fid, '%s\n', exactRef{:});
            fprintf(fid, '\n');
        end
        fprintf(fid, '\n');
        if anyFail
            fprintf(fid, '### Failures (first few per check)\n\n');
            for i = 1:numel(rows)
                for j = 1:numel(rows(i).failDetail)
                    fprintf(fid, '- **%s**: %s\n', rows(i).name, rows(i).failDetail{j});
                end
            end
            fprintf(fid, '\n');
        end
        if any([rows.roundoff] > 0)
            fprintf(fid, '### Roundoff-explained strict misses (first few per check)\n\n');
            for i = 1:numel(rows)
                for j = 1:numel(rows(i).roundoffDetail)
                    fprintf(fid, '- **%s**: %s\n', rows(i).name, rows(i).roundoffDetail{j});
                end
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '## Self-review: each .m beside its .cpp\n\n');
        fprintf(fid, 'The comment columns and the flight-source scan are checked automatically on every run: every MATLAB comment (`%%`/`%%%%`, trailing included) must appear verbatim as a `//` comment in the same order, and the `.cpp` may carry no other comment except the file header and `TRANSPILE-FLAG`s. The scan rejects MEX, stream/file I/O, printing and timing calls in flight sources (the only `snprintf`/`abort` are inside the `POLYQP_ERROR_MESSAGES`/`POLYQP_DEBUG_BOUNDS` blocks). Statement order and deviations come from a manual side-by-side read of each file (2026-11-18, `abort-cpp`); the `abort-fixes` edits were written side by side with their `.m` change (2026-11-19) and their `FIX-Fxx` comments are covered by the automated comment-order check, but the unchanged parts were not re-read.\n\n');
        fprintf(fid, '| .cpp | .m | MATLAB comments | C++ comments (excl. header/flags) | Verbatim, in order | Flight-source scan | Statement order (manual) | Unflagged deviations (manual) |\n|---|---|---|---|---|---|---|---|\n');
        for i = 1:numel(selfReview)
            fprintf(fid, '| `%s` | `%s` | %d | %d | %s | %s | matches | none |\n', selfReview(i).cpp, selfReview(i).m, selfReview(i).nM, selfReview(i).nC, ...
                passStr(selfReview(i).commentsOk), passStr(selfReview(i).scanOk));
        end
        fprintf(fid, '\n');
        fprintf(fid, '## Coverage notes\n\n');
        fprintf(fid, '- Local functions are covered through their parents: `solveAllLegs`, `buildReport` (via `PolyQP`); `wpToFix`, `makeCruiseKnot`, `pinTransverse`, `initSingleSegTime` (via `resolveWaypoints`); `softmaxTimes` (via `allocateTimes`); `solveThrustStretch`, `solveThrustFloorStretch` (via `evaluateConstraints`).\n');
        fprintf(fid, '- ASTRA only (vehicle 0, r = 5): the C++ solver is ASTRA-only with compiled-in constants, so the TOAD and r ~= 5 cases of the transpile verification were dropped. Abort cases: %d; x, y in [-20, 20] m, z in [0.5, 40] m, vx, vy in [-5, 5] m/s, vz in [-6, 4] m/s, quaternion [1 0 0 0], zero propellant. Waypoints from `abortWaypoints` (lines copied verbatim from `Guidance/AbortGuidance.m`; the script asserts they still match the file). Plus %d tilted, rotating starts (tilt <= 30 deg about a random horizontal axis, random heading, |w_B| <= 1 rad/s) through PolyQP, AbortGuidance and every per-function check, on a separate RNG stream.\n', nAbort, nTilt);
        fprintf(fid, '- `AbortGuidance` is compared with `Guidance/AbortGuidance.m`, the abort block (moved out of `TrajectoryManager.m`) with the persistents passed explicitly and `ABORT` set to 1 (the file hard-codes 0, so the block can never run as shipped). Each case checks: fire at t = 2, no recompute on a second call (abortCalculated/tAbort persist), and no fire with the default trigger (C++ default ABORT = 0).\n');
        fprintf(fid, '- `flatnessInvert`, `allocateControls`, `evaluateConstraints`, `sampleOutputs`, `PolyEval` use `QPTraj` from MATLAB `PolyQP` and `Plan` from `resolveWaypoints` with `Plan.SegT(1:M) = QPTraj.T(1:M)`. `flatnessInvert` runs on the `evaluateConstraints` grid; `allocateControls` on MATLAB''s `Wrench`.\n');
        fprintf(fid, '- `solveLeg` inputs from runs are rebuilt exactly as `solveAllLegs` pass 1 builds `FixP`, at both the initial and the final segment times. For the abort, passes 2-3 do not change `FixP` (F04).\n');
        fprintf(fid, '- `resolveWaypoints`/`PolyQP` error IDs: 1 waypoint, r = 2, r = 7, zero first HDGRef, inconsistent HDGRef, zero MaxVel, negative MaxAccel, scaled-but-consistent HDGRef (no error). The MEX rejects more than 16 waypoints before the call (capacity), so the > MAX_WP branch is checked only through the count-of-1 case.\n\n');
        fprintf(fid, '## Abort summary (MATLAB results)\n\nFeasibility is reported separately for start velocities within the segment-1 caps [4; 4; 5] m/s and over them (FIX-F01 / D3).\n\n| Vehicle | Cases | Feasible | Within cap: cases / feasible | Over cap: cases / feasible | OuterIterations histogram (0,1,2,...) | T_total min / median / max (s) |\n|---|---|---|---|---|---|---|\n');
        for i = 1:numel(abortStats)
            h = accumarray(abortStats(i).outerHist(:) + 1, 1).';
            fprintf(fid, '| %s | %d | %d | %d / %d | %d / %d | %s | %.2f / %.2f / %.2f |\n', abortStats(i).veh, abortStats(i).n, abortStats(i).feasible, ...
                abortStats(i).n - abortStats(i).nOver, abortStats(i).feasible - abortStats(i).feasOver, abortStats(i).nOver, abortStats(i).feasOver, mat2str(h), ...
                min(abortStats(i).T_total), median(abortStats(i).T_total), max(abortStats(i).T_total));
        end
        fprintf(fid, '\nTilted set (%d; tilt <= 30 deg, |w| <= 1 rad/s): feasible %d / %d.\n', nTilt, tiltStats(1, 2), tiltStats(1, 1));
        fprintf(fid, '\n## Missions\n\n');
        fprintf(fid, '`BuildTrajectory.m` is a script (`clear; clc; close all`, CSV export, plots), so it is not run. Its one mission (Circle) defines its waypoints without Simulink; those lines are copied verbatim into `circleMissionWaypoints`. Two synthetic missions add an initial hold, a cruise split and free pass-through corners.\n\n');
        fprintf(fid, '| Mission / vehicle | Segs | Legs | Holds | Cruise | Feasible | OuterIter (MATLAB) | T_total (s) | OuterIter (C++) | Max end-to-end err |\n|---|---|---|---|---|---|---|---|---|---|\n');
        fprintf(fid, '%s\n', missionLog{:});
        fprintf(fid, '\n## MATLAB comment correspondence\n\n');
        fprintf(fid, 'The C++ comments of each file (C++-only tags such as `TRANSPILE-FLAG`, `ASTRA-ONLY`, `NO-NAN`, `STATIC`, `STREAM` excluded) must be the MATLAB comments in order. A MATLAB comment may be missing only where the code it annotates was deleted; every such comment is listed here.\n\n');
        fprintf(fid, '| MATLAB | C++ | MATLAB comments | Kept in C++ | Result |\n|---|---|---|---|---|\n');
        for i = 1:numel(selfReview)
            fprintf(fid, '| `%s` | `%s` | %d | %d | %s |\n', selfReview(i).m, selfReview(i).cpp, selfReview(i).nM, selfReview(i).nC, passStr(selfReview(i).ok));
        end
        for i = 1:numel(selfReview)
            if ~isempty(selfReview(i).dropped)
                fprintf(fid, '\nDropped from `%s` (code deleted):\n', selfReview(i).cpp);
                for j = 1:numel(selfReview(i).dropped)
                    fprintf(fid, '- `%s`\n', selfReview(i).dropped(j));
                end
            end
        end
        fprintf(fid, '\n# Abort-fix regressions (`fix_regressions.m`)\n\n');
        fprintf(fid, 'One section per FIX in `ABORT_FIXES.md`, checking its acceptance criteria on the fixed MATLAB and C++. Overall: %s.\n\n', passStr(all([fixR.pass])));
        for i = 1:numel(fixR)
            fprintf(fid, '%s\n', strrep(fixR(i).md, '## ', '## '));
        end
    end
end

% ================================================================ local helpers

function s = ternary(c, a, b)
    if c
        s = a;
    else
        s = b;
    end
end

function checkAbortCopyIsVerbatim(toadDir, here)
% Every waypoint line of the abort block (Guidance/AbortGuidance.m) must appear verbatim in abortWaypoints.m.
    tm = readlines(fullfile(toadDir, 'Guidance', 'AbortGuidance.m'));
    wp = fileread(fullfile(here, 'abortWaypoints.m'));
    i1 = find(contains(tm, 'if ABORT > 0 && t > ABORT && ~abortCalculated'), 1);
    i2 = find(contains(tm, 'tAbort = t;'), 1);
    assert(~isempty(i1) && ~isempty(i2) && i2 > i1, 'Abort block not found in AbortGuidance.m');
    for i = (i1 + 1):i2
        ln = strtrim(tm(i));
        if strlength(ln) == 0
            continue;
        end
        if contains(ln, 'TOADWaypoint(') || contains(ln, 'repmat(dummyWP')
            assert(contains(wp, ln), 'verify_polyqp:abortCopy', 'Abort waypoint line not copied verbatim into abortWaypoints.m: %s', ln);
        end
    end
end

function Waypoints = circleMissionWaypoints()
% Waypoint lines copied verbatim from BuildTrajectory.m (mission "Circle").
dummyWP   = TOADWaypoint([0; 0; 0]);
Waypoints = repmat(dummyWP, 1, 9);
Radius = 5;
Height = 7;

Waypoints(1) = TOADWaypoint([0; 0; 0], 'Label', 'FLT STRT');
Waypoints(2) = TOADWaypoint([0; 0; 2], 'Label', 'FLT 1', 'StraightLine', true);
Waypoints(3) = TOADWaypoint([Radius; 0; Height], 'Label', 'FLT 2', 'VelBC', [NaN, NaN, 0]);
Waypoints(4) = TOADWaypoint([0; Radius; Height], 'Label', 'FLT 3', 'VelBC', [NaN, NaN, 0]);
Waypoints(5) = TOADWaypoint([-Radius; 0; Height], 'Label', 'FLT 4', 'VelBC', [NaN, NaN, 0]);
Waypoints(6) = TOADWaypoint([0; -Radius; Height], 'Label', 'FLT 5', 'VelBC', [NaN, NaN, 0]);
Waypoints(7) = TOADWaypoint([Radius; 0; Height], 'Label', 'FLT 6');
Waypoints(8) = TOADWaypoint([0; 0; 2], 'Label', 'FLT 7', 'HoldTime', 1.0);
Waypoints(9) = TOADWaypoint([0; 0; 0], 'Label', 'FLT END', 'HoldTime', 1.0, 'MaxVel', [3; 3; 1], 'StraightLine', true);
end

function Waypoints = cruiseMissionWaypoints()
% Synthetic: initial hold, long straight rest-to-rest transitions (cruise split), StraightLine.
dummyWP   = TOADWaypoint([0; 0; 0]);
Waypoints = repmat(dummyWP, 1, 4);
Waypoints(1) = TOADWaypoint([0; 0; 0], 'Label', 'START', 'HoldTime', 1.0);
Waypoints(2) = TOADWaypoint([0; 0; 12], 'Label', 'UP', 'HoldTime', 0.5, 'StraightLine', true);
Waypoints(3) = TOADWaypoint([30; 0; 12], 'Label', 'CRUISE', 'HoldTime', 0.5, 'MaxVel', [3; 3; 2]);
Waypoints(4) = TOADWaypoint([30; 0; 0], 'Label', 'LAND', 'HoldTime', 1.0, 'MaxVel', [3; 3; 1], 'StraightLine', true);
end

function Waypoints = cornerMissionWaypoints()
% Synthetic: free pass-through corners (multi-segment leg, allocateTimes, P-space refinement).
dummyWP   = TOADWaypoint([0; 0; 0]);
Waypoints = repmat(dummyWP, 1, 5);
Waypoints(1) = TOADWaypoint([0; 0; 0], 'Label', 'START');
Waypoints(2) = TOADWaypoint([3; 1; 4], 'Label', 'C1');
Waypoints(3) = TOADWaypoint([6; -2; 6], 'Label', 'C2', 'AccelBC', [NaN; NaN; 0]);
Waypoints(4) = TOADWaypoint([2; -4; 3], 'Label', 'C3', 'VelBC', [NaN; 0; NaN]);
Waypoints(5) = TOADWaypoint([0; 0; 0], 'Label', 'END', 'HoldTime', 0.0, 'MaxVel', [3; 3; 1]);
end

function X = tiltedAbortState(rs, c)
% Tilt <= 30 deg about a random horizontal axis, random heading, |w_B| <= 1 rad/s, position/velocity as the verification box, TOAD tanks 0-100 %.
    tilt = deg2rad(30) * rand(rs);
    az = 2 * pi * rand(rs);
    psi = 2 * pi * rand(rs) - pi;
    ax = [cos(az); sin(az); 0];
    K = [0, -ax(3), ax(2); ax(3), 0, -ax(1); -ax(2), ax(1), 0];
    Rt = eye(3) + sin(tilt) * K + (1 - cos(tilt)) * (K * K);
    Rz = [cos(psi), -sin(psi), 0; sin(psi), cos(psi), 0; 0, 0, 1];
    q = DCM_Quat_Conversion(Rz * Rt);
    w = randn(rs, 3, 1);
    w = w / norm(w) * rand(rs);
    X = [q; -20 + 40 * rand(rs); -20 + 40 * rand(rs); 0.5 + 39.5 * rand(rs); -5 + 10 * rand(rs); -5 + 10 * rand(rs); -6 + 10 * rand(rs); w; ...
        rand(rs) * c.OxMass; rand(rs) * c.FuMass];
end

function FixLeg = abortLikeFixLeg(r)
% Shape of the abort leg: knot 1 pos+vel fixed, knots 2-3 lateral position free, higher derivatives 0.
    FixLeg = nan(3, r, 3);
    FixLeg(:, 1, 1) = [3; -2; 10.2];
    FixLeg(:, 2, 1) = [1; 0.5; -1];
    FixLeg(:, 1, 2) = [NaN; NaN; 2.2];
    FixLeg(:, 2:r, 2) = 0;
    FixLeg(:, 1, 3) = [NaN; NaN; 0.2];
    FixLeg(:, 2:r, 3) = 0;
end

function lc = logCondRPP(FixLeg, SegTimes, H0)
% max over axes of log10 cond(R_PP_eq), the matrix solveLeg.m factorises (FIX-F07 reduced system included; NaN if none).
    Req = rppEqAsSolved(FixLeg, SegTimes, H0);
    lc = NaN;
    for d = 1:numel(Req)
        if ~isempty(Req{d})
            lc = max(lc, log10(cond(Req{d})));
        end
    end
end

function fl = quadFormFloor(FixLeg, SegTimes, H0, Deriv, Cost)
% nv * eps * (|d|' |Rmat| |d|) / |Cost| with Rmat assembled as solveLeg does (all segments), d = Deriv per axis.
    r = size(FixLeg, 2);
    K = size(FixLeg, 3);
    R = zeros(r * K);
    for m = 1:(K - 1)
        s = [SegTimes(m) .^ (0:(r - 1)), SegTimes(m) .^ (0:(r - 1))]';
        ix = (m - 1) * r + (1:(2 * r));
        R(ix, ix) = R(ix, ix) + SegTimes(m) ^ (1 - 2 * r) * (s * s') .* H0;
    end
    t = 0;
    for d = 1:size(Deriv, 1)
        dv = reshape(Deriv(d, :, :), [], 1);
        t = t + abs(dv)' * abs(R) * abs(dv);
    end
    fl = r * K * eps * t / max(abs(Cost), realmin);
end

function rc = minRcondRPP(FixLeg, SegTimes, H0)
    Req = rppEqAsSolved(FixLeg, SegTimes, H0);
    rc = NaN;
    for d = 1:numel(Req)
        if ~isempty(Req{d})
            rc = min(rc, rcond(Req{d}));
        end
    end
end

function [ok, err] = cmpNum(a, b, mode, tol)
% mode: 'rel' | 'abs' | 'relabs' (pass if rel or abs within tol) | 'exact'
    err = 0;
    if islogical(a) ~= islogical(b) && strcmp(mode, 'exact')
        ok = false; err = Inf; return;
    end
    if ~isequal(size(a), size(b))
        ok = false; err = Inf; return;
    end
    a = double(a); b = double(b);
    if strcmp(mode, 'exact')
        ok = isequaln(a, b);
        err = double(~ok);
        return;
    end
    na = isnan(a); nb = isnan(b);
    ia = isinf(a); ib = isinf(b);
    if any(na(:) ~= nb(:)) || any(ia(:) ~= ib(:)) || any(a(ia) ~= b(ib))
        ok = false; err = Inf; return;
    end
    m = ~na & ~ia;
    if ~any(m(:))
        ok = true; return;
    end
    d = max(abs(a(m) - b(m)));
    s = max(max(abs(a(m))), 1e-300);
    rel = d / s;
    switch mode
        case 'rel'
            err = rel;
        case 'abs'
            err = d;
        case 'relabs'
            err = min(rel, d);
    end
    ok = err <= tol;
end

function [ok, err, msg] = cmpStruct(A, B, mode, tol, exactFields)
% Compares field names (and order), then each field: exactFields exactly, others numerically.
% mode 'exactonly' checks only exactFields.
    ok = true; err = 0; msg = '';
    fa = fieldnames(A); fb = fieldnames(B);
    if ~isequal(fa, fb)
        ok = false; err = Inf; msg = sprintf('field names differ: {%s} vs {%s}', strjoin(fa.', ','), strjoin(fb.', ','));
        return;
    end
    for i = 1:numel(fa)
        f = fa{i};
        a = A.(f); b = B.(f);
        isExact = any(strcmp(f, exactFields));
        if strcmp(mode, 'exactonly') && ~isExact
            continue;
        end
        if ~strcmp(mode, 'exactonly') && isExact && ~strcmp(mode, 'rel')
            continue;
        end
        if ischar(a) || ischar(b)
            okf = isequal(a, b); e = double(~okf);
        elseif isExact
            okf = isequal(class(a), class(b)) || (isnumeric(a) && isnumeric(b));
            okf = okf && isequal(islogical(a), islogical(b)) && isequal(size(a), size(b)) && isequaln(double(a), double(b));
            e = double(~okf);
            if strcmp(mode, 'rel')
                e = 0;
            end
        else
            [okf, e] = cmpNum(a, b, mode, tol);
        end
        if ~okf
            ok = false;
            msg = [msg, sprintf('%s(%s) ', f, sizeStr(a, b, e))]; %#ok<AGROW>
        end
        if ~isinf(e)
            err = max(err, e);
        elseif ~okf
            err = Inf;
        end
    end
    if ok
        msg = 'ok';
    end
end

function s = sizeStr(a, b, e)
    if ~isequal(size(a), size(b))
        s = sprintf('size %s vs %s', mat2str(size(a)), mat2str(size(b)));
    else
        s = sprintf('%.3g', e);
    end
end

function [ok, err, msg] = sameShapes(A, B)
    ok = isequal(fieldnames(A), fieldnames(B));
    err = NaN;
    msg = 'ok';
    if ~ok
        msg = 'field names differ';
        return;
    end
    f = fieldnames(A);
    for i = 1:numel(f)
        if ~isequal(size(A.(f{i})), size(B.(f{i})))
            ok = false;
            msg = sprintf('%s size %s vs %s', f{i}, mat2str(size(A.(f{i}))), mat2str(size(B.(f{i}))));
            return;
        end
    end
end

function ps = pstream()
% Separate RNG stream for roundoff-floor perturbations (keeps the rng(42) case stream independent of results).
    persistent s
    if isempty(s)
        s = RandStream('mt19937ar', 'Seed', 7);
    end
    ps = s;
end

function y = pert(x)
% Multiply every element by 1 + 2*eps*randn (NaN, Inf and 0 are unchanged).
    y = x .* (1 + 2 * eps * randn(pstream(), size(x)));
end

function fl = floorAllocateTimes(FixLeg, T0, A0, H0, MinSegTime, Tml, Jml)
    fl = 0;
    for p = 1:3
        try
            [T2, J2] = allocateTimes(pert(FixLeg), pert(T0), A0, H0, MinSegTime);
        catch
            continue;
        end
        [~, e1] = cmpNum(Tml, T2, 'rel', Inf);
        [~, e2] = cmpNum(Jml, J2, 'rel', Inf);
        fl = max([fl, e1, e2]);
    end
end

function fl = floorPolyEval(QT, tq, k, Out)
    fl = 0;
    for p = 1:3
        QT2 = QT;
        QT2.Coef = pert(QT.Coef);
        [~, e] = cmpNum(Out, PolyEval(QT2, tq, k), 'rel', Inf);
        fl = max(fl, e);
    end
end

function fl = floorPolyQP(Waypoints, c, r, IG, QT, RP, Opts)
% max of (a) MATLAB's own error on its final QP solve against exact arithmetic (abs, on QPTraj.Coef/Fix)
% and (b) MATLAB's self-sensitivity to ~2 ulp waypoint perturbations (same metric as the check).
    fl = exactFinalSolveFloor(Waypoints, c, r, QT);
    for p = 1:3
        W2 = Waypoints;
        for i = 1:numel(W2)
            W2(i).PosBC = pert(W2(i).PosBC);
            W2(i).VelBC = pert(W2(i).VelBC);
        end
        try
            if nargin < 7
                Opts = struct('R', r);
            end
            [IG2, QT2, RP2] = PolyQP(W2, c, Opts);
        catch
            continue;
        end
        if RP2.OuterIterations ~= RP.OuterIterations || ~isequal(size(IG2.X), size(IG.X))
            continue;
        end
        [~, e1] = cmpStruct(IG, IG2, 'abs', Inf, {});
        [~, e2] = cmpStruct(QT, QT2, 'abs', Inf, {'R', 'Count', 'IsHold', 'IsCruise'});
        [~, e3] = cmpStruct(RP, RP2, 'abs', Inf, {'Feasible', 'OuterIterations'});
        fl = max([fl, e1, e2, e3]);
    end
end

function fl = floorAbort(t, X_est, c, ABORT, S0, S1)
    fl = 0;
    for p = 1:3
        X2 = X_est;
        X2(5:10) = pert(X_est(5:10));
        try
            S2 = AbortGuidance(t, X2, c, ABORT, S0);
        catch
            continue;
        end
        if ~isequal(size(S2.States), size(S1.States))
            continue;
        end
        [~, e1] = cmpNum(S1.Time, S2.Time, 'abs', Inf);
        [~, e2] = cmpNum(S1.States, S2.States, 'abs', Inf);
        [~, e3] = cmpNum(S1.Inputs, S2.Inputs, 'abs', Inf);
        fl = max([fl, e1, e2, e3]);
    end
end

function fl = exactFinalSolveFloor(Waypoints, c, r, QT)
% Rebuild the last solveAllLegs pass: fixed entries of FixP are carried unchanged into QPTraj.Fix,
% the free pattern is that of Plan.FixCg; solve every leg exactly and measure MATLAB's abs error.
    Lim = buildVehicleLimits(c);
    Plan = resolveWaypoints(Waypoints, r, Lim);
    K = Plan.KnotCount;
    nc = 2 * r;
    FixP = QT.Fix(:, 1:r, 1:K);
    FixP(isnan(Plan.FixCg(:, 1:r, 1:K))) = NaN;
    [A0, ~, H0] = assembleQP(r);
    fl = 0;
    for leg = 1:Plan.LegCount
        m1 = Plan.LegStartSeg(leg);
        m2 = Plan.LegEndSeg(leg);
        [Cex, Dex] = exactSolveLeg(FixP(:, :, m1:(m2 + 1)), QT.T(m1:m2), A0, H0);
        fl = max(fl, max(abs(QT.Coef(1:nc, m1:m2, :) - Cex), [], 'all'));
        fl = max(fl, max(abs(QT.Fix(:, 1:r, m1:(m2 + 1)) - Dex), [], 'all'));
    end
end

function [CoefPhys, Deriv, Cost] = exactSolveLeg(FixLeg, SegTimes, A0, H0)
% solveLeg.m in 60-digit vpa from the identical double inputs (no equilibration: it is exact in exact arithmetic).
    oldDigits = digits(60);
    restoreDigits = onCleanup(@() digits(oldDigits));
    D = size(FixLeg, 1);
    r = size(FixLeg, 2);
    K = size(FixLeg, 3);
    M = K - 1;
    nc = 2 * r;
    nv = r * K;
    H0s = vpa(sym(H0, 'f'));
    A0s = vpa(sym(A0, 'f'));
    Rmat = vpa(zeros(nv, nv));
    Ts = vpa(sym(SegTimes, 'f'));
    for m = 1:M
        Tm = Ts(m);
        s = [Tm .^ (0:(r - 1)), Tm .^ (0:(r - 1))].';
        Hm = (Tm ^ (1 - 2 * r)) * (s * s.') .* H0s;
        ix = (m - 1) * r + (1:nc);
        Rmat(ix, ix) = Rmat(ix, ix) + Hm;
    end
    DerivS = vpa(zeros(D, nv));
    Cost = vpa(0);
    for d = 1:D
        fx = reshape(FixLeg(d, :, :), nv, 1);
        idxF = find(~isnan(fx));
        idxP = find(isnan(fx));
        dv = vpa(zeros(nv, 1));
        dv(idxF) = vpa(sym(fx(idxF), 'f'));
        if ~isempty(idxP)
            dv(idxP) = -(Rmat(idxP, idxP) \ (Rmat(idxP, idxF) * dv(idxF)));
        end
        DerivS(d, :) = dv.';
        Cost = Cost + dv.' * Rmat * dv;
    end
    D0_inv = diag(1 ./ diag(A0s(1:r, 1:r)));
    A10 = A0s((r + 1):nc, 1:r);
    A11 = A0s((r + 1):nc, (r + 1):nc);
    CoefPhys = zeros(nc, M, D);
    for m = 1:M
        Tm = Ts(m);
        s = [Tm .^ (0:(r - 1)), Tm .^ (0:(r - 1))].';
        scalePhys = (Tm .^ (0:(nc - 1))).';
        for d = 1:D
            d_seg = s .* DerivS(d, ((m - 1) * r + 1):((m + 1) * r)).';
            c_top = D0_inv * d_seg(1:r);
            c_bot = A11 \ (d_seg((r + 1):nc) - A10 * c_top);
            CoefPhys(:, m, d) = double([c_top; c_bot] ./ scalePhys);
        end
    end
    Deriv = reshape(double(DerivS), D, r, K);
    Cost = double(Cost);
end

function R = rotx(a)
    R = [1, 0, 0; 0, cos(a), -sin(a); 0, sin(a), cos(a)];
end
function R = roty(a)
    R = [cos(a), 0, sin(a); 0, 1, 0; -sin(a), 0, cos(a)];
end
function R = rotz(a)
    R = [cos(a), -sin(a), 0; sin(a), cos(a), 0; 0, 0, 1];
end

function R = quatToDCM(q)
    w = q(1); x = q(2); y = q(3); z = q(4);
    R = [1 - 2 * (y^2 + z^2), 2 * (x * y - w * z), 2 * (x * z + w * y);
         2 * (x * y + w * z), 1 - 2 * (x^2 + z^2), 2 * (y * z - w * x);
         2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x^2 + y^2)];
end

function p = prctile95(x)
    x = sort(x(:));
    p = x(max(1, ceil(0.95 * numel(x))));
end

function s = passStr(b)
    if b
        s = 'pass';
    else
        s = '**FAIL**';
    end
end

function R = selfReviewRows(toadDir, polyDir)
    pairs = {
        'Guidance/PolynomialQP/PolyQP.m',              'PolyQP.cpp',              [];
        'Guidance/PolynomialQP/resolveWaypoints.m',    'resolveWaypoints.cpp',    [];
        'Guidance/PolynomialQP/assembleQP.m',          'assembleQP.cpp',          [];
        'Guidance/PolynomialQP/solveLeg.m',            'solveLeg.cpp',            [];
        'Guidance/PolynomialQP/allocateTimes.m',       'allocateTimes.cpp',       [];
        'Guidance/PolynomialQP/evaluateConstraints.m', 'evaluateConstraints.cpp', [];
        'Guidance/PolynomialQP/flatnessInvert.m',      'flatnessInvert.cpp',      [];
        'Guidance/PolynomialQP/allocateControls.m',    'allocateControls.cpp',    [];
        'Guidance/PolynomialQP/sampleOutputs.m',       'sampleOutputs.cpp',       [];
        'Guidance/PolynomialQP/PolyEval.m',            'PolyEval.cpp',            [];
        'Guidance/PolynomialQP/buildVehicleLimits.m',  'buildVehicleLimits.cpp',  [];
        'Guidance/PolynomialQP/TOADWaypoint.m',        'TOADWaypoint.cpp',        [];
        'Helper/ComputeJtot.m',                        '../matlab_helpers.cpp',   'region';
        'Helper/DCM_Quat_Conversion.m',                '../matlab_funcs.h',       'region';
        'Guidance/AbortGuidance.m',                    'AbortGuidance.cpp',       'abort';
        'Guidance/PolynomialQP/abortFallback.m',       'abortFallback.cpp',       [];
        'Helper/quatRot.m',                            '../matlab_funcs.h',       'region';
        };
    R = struct('m', {}, 'cpp', {}, 'nM', {}, 'nC', {}, 'commentsOk', {}, 'scanOk', {}, 'ok', {}, 'msg', {}, 'dropped', {});
    % C++-only comment tags (not MATLAB comments): transpile flags, file headers, and the notes of the ASTRA integration
    % (ASTRA-ONLY: vehicle 1 / r ~= 5 removal, NO-NAN: NaN sentinel -> masks, STATIC: fixed-size storage, STREAM: streamed grids).
    cppOnly = {'TRANSPILE-FLAG', 'Transpiled from', 'ASTRA-ONLY', 'NO-NAN', 'STATIC', 'STREAM'};
    forbidden = '(mex\.h|\bmx[A-Z]\w*\(|iostream|fstream|std::cout|std::cerr|\bfopen\b|\bprintf\(|<chrono>|\bclock\()';
    for i = 1:size(pairs, 1)
        mLines = readlines(fullfile(toadDir, pairs{i, 1}));
        if strcmp(pairs{i, 3}, 'abort')
            i1 = find(contains(mLines, 'ABORT = 0'), 1);
            i2 = find(contains(mLines, 'tAbort = t;'), 1) + 1;
            mLines = mLines(i1:i2);
        end
        cText = fileread(fullfile(polyDir, pairs{i, 2}));
        cLines = splitlines(string(cText));
        if strcmp(pairs{i, 3}, 'region')
            % Helpers moved into cpp/matlab_helpers.cpp / matlab_funcs.h: the region from this helper's
            % "// Transpiled from TOAD/<file>.m" header to the next "// Transpiled from" header (or end of file).
            j1 = find(contains(cLines, ['// Transpiled from TOAD/', pairs{i, 1}]), 1);
            assert(~isempty(j1), 'selfReview:region', 'No "Transpiled from TOAD/%s" header in %s', pairs{i, 1}, pairs{i, 2});
            j2 = find(contains(cLines((j1 + 1):end), '// Transpiled from'), 1);
            if isempty(j2)
                j2 = numel(cLines);
            else
                j2 = j1 + j2 - 1;
            end
            cLines = cLines(j1:j2);
            cText = char(strjoin(cLines, newline));
        end
        mc = matlabComments(mLines);
        cc = strings(0, 1);
        for j = 1:numel(cLines)
            tok = regexp(cLines(j), '//\s?(.*)$', 'tokens', 'once');
            if isempty(tok) || mod(count(extractBefore(cLines(j), '//'), '"'), 2) == 1
                continue;
            end
            t = strip(tok{1}, 'right');
            if startsWith(t, cppOnly)
                continue;
            end
            cc(end + 1, 1) = t; %#ok<AGROW>
        end
        % The C++ comments must be the MATLAB comments in order; a MATLAB comment may be missing only where its code was
        % deleted (vehicle 1, r ~= 5, NaN/finiteness checks). Every dropped comment is listed in the report.
        [commentsOk, dropped] = orderedSubsequence(cc, mc);
        scanOk = isempty(regexp(cText, forbidden, 'once'));
        msg = 'ok';
        if ~commentsOk
            msg = sprintf('C++ comment not found (in order) among the MATLAB comments: "%s"', dropped);
        elseif ~scanOk
            msg = 'verification code in flight source';
        end
        if ~commentsOk
            dropped = strings(0, 1);
        end
        R(end + 1) = struct('m', pairs{i, 1}, 'cpp', pairs{i, 2}, 'nM', numel(mc), 'nC', numel(cc), ...
            'commentsOk', commentsOk, 'scanOk', scanOk, 'ok', commentsOk && scanOk, 'msg', msg, 'dropped', {dropped}); %#ok<AGROW>
    end
end

function [ok, out] = orderedSubsequence(cc, mc)
% ok: every element of cc appears in mc in the same order. out: the mc elements not matched (or, if not ok, the first unmatched cc element).
    j = 0;
    used = false(numel(mc), 1);
    for i = 1:numel(cc)
        k = find(mc((j + 1):end) == cc(i), 1);
        if isempty(k)
            ok = false;
            out = cc(i);
            return;
        end
        j = j + k;
        used(j) = true;
    end
    ok = true;
    out = mc(~used);
end

function c = matlabComments(lines)
% MATLAB comment text (after the % markers), skipping % inside char literals; ' after an identifier, ), ], }, . or ' is a transpose.
    c = strings(0, 1);
    for j = 1:numel(lines)
        l = char(lines(j));
        inStr = false;
        prev = ' ';
        k = 1;
        while k <= numel(l)
            ch = l(k);
            if inStr
                if ch == ''''
                    if k < numel(l) && l(k + 1) == ''''
                        k = k + 2;
                        continue;
                    end
                    inStr = false;
                    prev = '''';
                end
            else
                if ch == '%'
                    c(end + 1, 1) = strip(string(regexprep(l(k:end), '^%+ ?', '')), 'right'); %#ok<AGROW>
                    break;
                end
                if ch == '''' && ~(isstrprop(prev, 'alphanum') || any(prev == ')]}._'''))
                    inStr = true;
                end
                if ~isspace(ch)
                    prev = ch;
                else
                    prev = ' ';
                end
            end
            k = k + 1;
        end
    end
end

function s = okOr(id)
    if isempty(id)
        s = 'ok';
    else
        s = id;
    end
end
