function PlotMC_v5(filename)
% PlotMC_v5
% Monte Carlo plotting for the V5 framework (built around per-run
% pass/fail, the physical wind envelope and the V5 disturbances.)
%
%   PlotMC_v5                  % after MonteCarlo_v5: reads the base workspace
%   PlotMC_v5('MC5_...mat')    % reads a saved run
%
% Figures
%   1  Summary: failure rate (90% Clopper-Pearson), failure modes, max errors
%   2  Trajectories (passed grey, failed red) + projections
%   3  Position / velocity distributions
%   4  Error histories against the pass/fail thresholds
%   5  Tilt and angular rates
%   6  Wind envelope: failure rate vs mean wind, wind direction
%   7  Failure drivers: ranking + top-driver scatters
%   8  Sensor disturbance diagnostics
%   9  Vehicle / actuator disturbance diagnostics
%   10 Filter RMSE (control RMSE only when plot_control_rmse = true)
%   11 Landing angular rates
%   12-14 Spectral radius and LESO lead/lag diagnostics (from V4)
%
% Time-history bands are percentiles across runs (median, 5-95%, min-max),
% not mean +/- 3 sigma: a few diverging runs inflate sigma at every time
% step, and the runs are not Gaussian once any of them fail.

close all;
bkgColor = 'w';
alphaVal = 0.12;

% Control RMSE comes from SSE_Controls, which integrates the signed (not
% squared) error in the current model. Leave false until that is fixed.
plot_control_rmse = false;

% Convention (NWU): X north, Y west, Z up
oranX = [0.85 0.60 0.45];
blueY = [0.55 0.71 0.84];
yellZ = [0.93 0.81 0.54];
axColors  = {oranX, blueY, yellZ};
axis_lbls = {'X (North)', 'Y (West)', 'Z (Up)'};
passColor = [0.45 0.68 0.88];
failColor = [0.85 0.33 0.10];
bandColor = [0.30 0.50 0.80];
clr.pass = passColor; clr.fail = failColor; clr.band = bandColor;

%% 1. Load Data
varList = {'t_common', 'pos_all', 'vel_all', 'quat_all', 'target_state_all', 'ok', ...
           'seed', 'pad_hold_time', 'pos_err_fail', 'att_err_fail', ...
           'wind_sweep', 'dryden_W6', ...
           'Wind_Mean_all', 'Wind_Dir_all', 'Wind_Gain_all', 'Wind_Covar_all', ...
           'Lever_Radial', 'Lever_Axial', 'J_d_vals', ...
           'AccelBias_Mag', 'GyroBias_Mag', 'HardIron_Mag', 'SoftIron_Dev', ...
           'GimbalBias_Mag', 'ThrustBias_all', 'IMUDriftParams_vals', ...
           'RMSE_Controls_all', 'RMSE_Filter_all', 'MaxSpectralRad_all', ...
           'specrad_all', 'LESODist_all', 'out'};
S = struct();
if nargin < 1 || isempty(filename)
    for k = 1:numel(varList)
        try
            S.(varList{k}) = evalin('base', varList{k});
        catch
            % not present in this run; sections that need it are skipped
        end
    end
else
    w = warning('off', 'MATLAB:load:variableNotFound');
    S = load(filename, varList{:});
    warning(w);
end

reqVars = {'t_common', 'pos_all', 'quat_all', 'target_state_all'};
missing = reqVars(~isfield(S, reqVars));
if ~isempty(missing)
    error('PlotMC_v5:MissingData', ...
        'Missing required variables: %s. Run MonteCarlo_v5 first or pass its .mat file.', ...
        strjoin(missing, ', '));
end

t_common = S.t_common(:);
T        = numel(t_common);
pos_all  = S.pos_all;
num_sims = size(pos_all, 1);
quat_all = S.quat_all;
tgt      = S.target_state_all;
vel_all  = getField(S, 'vel_all', nan(num_sims, T, 3));
ok       = logical(reshape(getField(S, 'ok', true(num_sims, 1)), [], 1));

pad_hold_time = getField(S, 'pad_hold_time', 15);
pos_err_fail  = getField(S, 'pos_err_fail', 2.0);
att_err_fail  = getField(S, 'att_err_fail', 15);

%% 2. Derived Errors & Pass/Fail
eval_idx = t_common >= pad_hold_time;
if ~any(eval_idx)
    warning('PlotMC_v5:NoEvalWindow', 'No samples after the pad hold; evaluating the full run.');
    eval_idx(:) = true;
end

tgt_pos = tgt(:,:,5:7);
pos_err = pos_all - tgt_pos;                               % N x T x 3
pos_err_norm = sqrt(sum(pos_err.^2, 3));                   % N x T
horiz_err    = sqrt(sum(pos_err(:,:,1:2).^2, 3));          % N x T
vert_err     = pos_err(:,:,3);                             % N x T (signed)
has_tgt_vel  = size(tgt, 3) >= 10;

q_sim = normQuat(quat_all);
q_trg = normQuat(tgt(:,:,1:4));
q_dot = abs(sum(q_sim .* q_trg, 3));
q_dot(q_dot > 1) = 1;                                      % keeps NaN as NaN
att_err     = 2 * acosd(q_dot);                            % N x T, deg
tilt_actual = tiltFromQuat(q_sim);                         % N x T, deg
tilt_target = tiltFromQuat(q_trg);

Pos_Err_Max   = max(pos_err_norm(:, eval_idx), [], 2);
Horiz_Err_Max = max(horiz_err(:, eval_idx), [], 2);
Vert_Err_Max  = max(abs(vert_err(:, eval_idx)), [], 2);
Att_Err_Max   = max(att_err(:, eval_idx), [], 2);

Fail_Error    = ~ok;
Fail_NoData   = ~Fail_Error & (isnan(Pos_Err_Max) | isnan(Att_Err_Max));
Fail_Position = Pos_Err_Max > pos_err_fail;
Fail_Attitude = Att_Err_Max > att_err_fail;
MC_Fail = Fail_Error | Fail_NoData | Fail_Position | Fail_Attitude;
nFail   = nnz(MC_Fail);
[ciLo, ciHi] = clopperPearson(nFail, num_sims, 0.90);

% First threshold violation after the pad hold (when failures happen)
viol = (pos_err_norm > pos_err_fail | att_err > att_err_fail) & eval_idx(:)';
t_first_fail = nan(num_sims, 1);
for i = 1:num_sims
    k = find(viol(i,:), 1, 'first');
    if ~isempty(k), t_first_fail(i) = t_common(k); end
end

% Reference target (same for every run)
tgtFinite  = all(isfinite(reshape(tgt_pos, num_sims, [])), 2);
firstValid = find(ok & tgtFinite, 1, 'first');
if isempty(firstValid), firstValid = find(tgtFinite, 1, 'first'); end
if isempty(firstValid)
    error('PlotMC_v5:NoTarget', 'No run has a finite target trajectory.');
end
target_pos_ref = reshape(tgt_pos(firstValid,:,:), T, 3);

fprintf('\n=== PlotMC_v5: %d runs, %d failed (%.1f%%), 90%% CI on P(fail) [%.3f, %.3f] ===\n', ...
    num_sims, nFail, 100 * nFail / num_sims, ciLo, ciHi);
fprintf('   sim error %d | no data %d | position > %.2f m %d | attitude > %.1f deg %d\n', ...
    nnz(Fail_Error), nnz(Fail_NoData), pos_err_fail, nnz(Fail_Position), ...
    att_err_fail, nnz(Fail_Attitude));

%% 3. Summary
figure('Name', 'V5 Summary', 'Color', bkgColor, 'WindowStyle', 'docked');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); axis(ax, 'off');
txt = {sprintf('Runs:             %d', num_sims), ...
       sprintf('Failed:           %d (%.1f%%)', nFail, 100 * nFail / num_sims), ...
       sprintf('P(fail) 90%% CI:   [%.3f, %.3f]', ciLo, ciHi), ...
       sprintf('Position limit:   %.2f m', pos_err_fail), ...
       sprintf('Attitude limit:   %.1f deg', att_err_fail), ...
       sprintf('Evaluated after:  t = %.1f s', pad_hold_time)};
if isfield(S, 'seed')
    txt{end+1} = sprintf('Master seed:      %.0f', double(S.seed));
end
if isfield(S, 'wind_sweep')
    if S.wind_sweep
        txt{end+1} = 'Wind campaign:    envelope sweep';
    else
        txt{end+1} = 'Wind campaign:    operational (Weibull)';
    end
end
if nFail == 0
    txt{end+1} = sprintf('0 failures: P(fail) < %.2f%% (95%%, rule of 3)', 300 / num_sims);
end
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'FixedWidth', ...
    'Interpreter', 'none', 'FontSize', 10);
title(ax, 'Campaign');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
modeCounts = [nnz(Fail_Error), nnz(Fail_NoData), nnz(Fail_Position), nnz(Fail_Attitude)];
modeLbls   = {'Sim error', 'No state data', 'Position', 'Attitude'};
barh(ax, 1:4, modeCounts, 'FaceColor', failColor, 'FaceAlpha', 0.85);
for k = 1:4
    text(ax, modeCounts(k), k, sprintf('  %d', modeCounts(k)), 'VerticalAlignment', 'middle');
end
set(ax, 'YTick', 1:4, 'YTickLabel', modeLbls, 'YDir', 'reverse');
xlim(ax, [0, max(1, max(modeCounts)) * 1.25]);
xlabel(ax, 'Runs (a run can fail more than one way)');
title(ax, 'Failure Modes');

ax = nexttile(tl);
plotPassFail(ax, Pos_Err_Max, Att_Err_Max, MC_Fail, 'Max |pos err| [m]', ...
    'Max att err [deg]', 'Failure Space', att_err_fail, clr);
xline(ax, pos_err_fail, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
legend(ax, 'Location', 'best');

ax = nexttile(tl);
plotSmartHistogram(ax, Pos_Err_Max, axColors{1}, 'Max |pos err|');
xline(ax, pos_err_fail, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Limit');
title(ax, 'Max Position Error (after pad hold)'); xlabel(ax, 'm'); grid(ax, 'on'); legend(ax, 'show');

ax = nexttile(tl);
plotSmartHistogram(ax, Att_Err_Max, axColors{2}, 'Max att err');
xline(ax, att_err_fail, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Limit');
title(ax, 'Max Attitude Error (after pad hold)'); xlabel(ax, 'deg'); grid(ax, 'on'); legend(ax, 'show');

ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
tf = t_first_fail(isfinite(t_first_fail));
if isempty(tf)
    text(ax, 0.5, 0.5, 'No threshold violations', 'HorizontalAlignment', 'center', 'Units', 'normalized');
else
    histogram(ax, tf, 'FaceColor', failColor, 'BinLimits', [pad_hold_time, t_common(end)]);
end
xlabel(ax, 'Time [s]'); ylabel(ax, 'Runs');
title(ax, 'Time of First Threshold Violation');

%% 4. Kinematic Trajectory Overlay (3D + Projections)
figure('Name', 'MC 3D Trajectories', 'Color', bkgColor, 'WindowStyle', 'docked');
tl = tiledlayout(3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

axMain = nexttile(tl, 1, [3 3]);
hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
xlabel(axMain, 'X (North) [m]'); ylabel(axMain, 'Y (West) [m]'); zlabel(axMain, 'Z (Up) [m]');
title(axMain, sprintf('MC Trajectories vs Target (%d runs, %d failed)', num_sims, nFail));
for pass = [false true]                     % passed first, failed on top
    idx = find(MC_Fail == pass)';
    if pass, c = [failColor, min(1, alphaVal * 3)]; else, c = [0.25 0.25 0.25 alphaVal]; end
    for i = idx
        xyz = reshape(pos_all(i,:,:), T, 3);
        valid = all(isfinite(xyz), 2);
        plot3(axMain, xyz(valid,1), xyz(valid,2), xyz(valid,3), 'Color', c, ...
            'LineWidth', 0.5, 'HandleVisibility', 'off');
    end
end
plot3(axMain, target_pos_ref(:,1), target_pos_ref(:,2), target_pos_ref(:,3), 'b--', ...
    'LineWidth', 2.5, 'DisplayName', 'Target');
plot3(axMain, nan, nan, nan, '-', 'Color', [0.25 0.25 0.25], 'DisplayName', 'Passed');
plot3(axMain, nan, nan, nan, '-', 'Color', failColor, 'DisplayName', 'Failed');
legend(axMain, 'Location', 'best');

axTop = nexttile(tl, 4);
drawMCTraj_v5(axTop, pos_all, target_pos_ref, 1, 2, 'X (North) [m]', 'Y (West) [m]', alphaVal, MC_Fail, failColor);
title(axTop, 'Top View');
axSide = nexttile(tl, 8);
drawMCTraj_v5(axSide, pos_all, target_pos_ref, 1, 3, 'X (North) [m]', 'Z (Up) [m]', alphaVal, MC_Fail, failColor);
title(axSide, 'Side View');
axFront = nexttile(tl, 12);
drawMCTraj_v5(axFront, pos_all, target_pos_ref, 2, 3, 'Y (West) [m]', 'Z (Up) [m]', alphaVal, MC_Fail, failColor);
title(axFront, 'Front View');

%% 5. State Distributions (Percentile Bands)
figure('Name', 'State Distributions', 'Color', bkgColor, 'WindowStyle', 'docked');
tiledlayout(2, 3, 'TileSpacing', 'compact');
blocks = {pos_all, vel_all};
tgtIdx = {5:7, 8:10};
ylbls  = {'Pos', 'Vel'};
units  = {'m', 'm/s'};
for r = 1:2
    for c = 1:3
        ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
        X = blocks{r}(:,:,c);
        plotBands(ax, t_common, X, bandColor);
        plotFailedRuns(ax, t_common, X, MC_Fail, failColor);
        if r == 1 || has_tgt_vel
            plot(ax, t_common, reshape(tgt(firstValid,:,tgtIdx{r}(c)), T, 1), 'k', ...
                'LineWidth', 1.5, 'DisplayName', 'Target');
        end
        shadePadHold(ax, t_common(1), pad_hold_time);
        title(ax, sprintf('%s - %s', ylbls{r}, axis_lbls{c}));
        xlabel(ax, 'Time (s)'); ylabel(ax, sprintf('%s (%s)', ylbls{r}, units{r}));
        if r == 1 && c == 1, legend(ax, 'Location', 'best'); end
    end
end

%% 6. Error Histories vs Pass/Fail Thresholds
figure('Name', 'Error Histories', 'Color', bkgColor, 'WindowStyle', 'docked');
tiledlayout(2, 2, 'TileSpacing', 'compact');
errSets = {pos_err_norm, horiz_err, vert_err, att_err};
errTtl  = {'|Position Error|', 'Horizontal Error', 'Vertical Error (signed)', 'Attitude Error'};
errLbl  = {'m', 'm', 'm', 'deg'};
errThr  = [pos_err_fail, NaN, NaN, att_err_fail];
for k = 1:4
    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plotBands(ax, t_common, errSets{k}, bandColor);
    plotFailedRuns(ax, t_common, errSets{k}, MC_Fail, failColor);
    if isfinite(errThr(k))
        yline(ax, errThr(k), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Limit');
    end
    shadePadHold(ax, t_common(1), pad_hold_time);
    title(ax, errTtl{k}); xlabel(ax, 'Time (s)'); ylabel(ax, errLbl{k});
    if k == 1, legend(ax, 'Location', 'best'); end
end

%% 7. Rotational Kinematics: Tilt Angle & Angular Rates
% Body rates are not stored by MonteCarlo_v5; they are read from `out`
% (state_log columns 11:13) when it is available.
rates_all = [];
if isfield(S, 'out') && ~isempty(S.out)
    rates_all = nan(num_sims, T, 3);
    for i = 1:min(num_sims, numel(S.out))
        if ~ok(i), continue; end
        try
            ts = S.out(i).state_log;
            raw = extractLoggedState(ts);
            if size(raw, 2) >= 13
                for k = 1:3
                    rates_all(i,:,k) = interp1(ts.Time, raw(:,10+k), t_common, 'linear', NaN);
                end
            end
        catch
        end
    end
    if all(isnan(rates_all(:))), rates_all = []; end
end

figure('Name', 'Attitude & Angular Rates', 'Color', bkgColor, 'WindowStyle', 'docked');
if isempty(rates_all)
    tl_att = tiledlayout(1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
else
    tl_att = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
end
axTilt = nexttile(tl_att, 1, [1, tl_att.GridSize(2)]); hold(axTilt, 'on'); grid(axTilt, 'on');
plotBands(axTilt, t_common, tilt_actual, bandColor);
plotFailedRuns(axTilt, t_common, tilt_actual, MC_Fail, failColor);
plot(axTilt, t_common, tilt_target(firstValid,:), 'k', 'LineWidth', 1.5, 'DisplayName', 'Target');
shadePadHold(axTilt, t_common(1), pad_hold_time);
title(axTilt, 'Tilt Angle from Vertical'); xlabel(axTilt, 'Time (s)'); ylabel(axTilt, 'Tilt (deg)');
legend(axTilt, 'Location', 'best');

if ~isempty(rates_all)
    rate_titles = {'Rate X', 'Rate Y', 'Rate Z'};
    for c = 1:3
        ax = nexttile(tl_att); hold(ax, 'on'); grid(ax, 'on');
        X = rates_all(:,:,c) * (180/pi);
        plotBands(ax, t_common, X, bandColor);
        plotFailedRuns(ax, t_common, X, MC_Fail, failColor);
        if size(tgt, 3) >= 13
            plot(ax, t_common, reshape(tgt(firstValid,:,10+c), T, 1) * (180/pi), 'k', ...
                'LineWidth', 1.5, 'DisplayName', 'Target');
        end
        shadePadHold(ax, t_common(1), pad_hold_time);
        title(ax, rate_titles{c}); xlabel(ax, 'Time (s)'); ylabel(ax, 'Rate (deg/s)');
    end
else
    disp('Body rates unavailable (no `out`) -- skipping rate panels and landing rates.');
end

%% 8. Wind Envelope
if isfield(S, 'Wind_Mean_all')
    W = reshape(S.Wind_Mean_all, [], 1);
    figure('Name', 'Wind Envelope', 'Color', bkgColor, 'WindowStyle', 'docked');
    tl = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on');
    if max(W) > 0
        nb = max(2, min(10, floor(num_sims / 10)));
        edges = linspace(0, max(W) * (1 + 1e-9), nb + 1);
        ctr = (edges(1:end-1) + edges(2:end)) / 2;
        nBin = zeros(1, nb); kBin = zeros(1, nb);
        for b = 1:nb
            inb = W >= edges(b) & W < edges(b+1);
            nBin(b) = nnz(inb); kBin(b) = nnz(inb & MC_Fail);
        end
        has = nBin > 0;
        rate = kBin(has) ./ nBin(has);
        [lo, hi] = clopperPearson(kBin(has), nBin(has), 0.90);
        % raw outcomes as a jittered rug at 0 (pass) and 1 (fail)
        jit = 0.03 * (mod((1:num_sims)' * 0.6180339887, 1) - 0.5);   % deterministic, leaves rng alone
        scatter(ax, W(~MC_Fail), 0 + jit(~MC_Fail), 12, passColor, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'DisplayName', 'Run (pass)');
        scatter(ax, W(MC_Fail), 1 + jit(MC_Fail), 12, failColor, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'DisplayName', 'Run (fail)');
        errorbar(ax, ctr(has), rate, rate - lo, hi - rate, 'ko-', 'LineWidth', 1.3, ...
            'MarkerFaceColor', 'k', 'DisplayName', 'Binned rate (90% CI)');
        idxHas = find(has);
        for b = 1:numel(idxHas)
            text(ax, ctr(idxHas(b)), min(1.0, hi(b)) + 0.04, sprintf('n=%d', nBin(idxHas(b))), ...
                'HorizontalAlignment', 'center', 'FontSize', 7);
        end
        ylim(ax, [-0.05, 1.12]);
        legend(ax, 'Location', 'northwest');
    end
    xlabel(ax, 'Mean wind at 6 m [m/s]'); ylabel(ax, 'P(fail)');
    ttl = 'Failure Rate vs Mean Wind';
    if isfield(S, 'wind_sweep') && ~S.wind_sweep
        ttl = [ttl ' (operational draw: few runs at high wind)'];
    end
    title(ax, ttl);

    ax = nexttile(tl);
    plotPassFail(ax, W, Pos_Err_Max, MC_Fail, 'Mean wind [m/s]', 'Max |pos err| [m]', ...
        'Position Error vs Mean Wind', pos_err_fail, clr);
    ax = nexttile(tl);
    plotPassFail(ax, W, Att_Err_Max, MC_Fail, 'Mean wind [m/s]', 'Max att err [deg]', ...
        'Attitude Error vs Mean Wind', att_err_fail, clr);

    if isfield(S, 'Wind_Dir_all')
        th = reshape(S.Wind_Dir_all, [], 1);
        try
            pax = polaraxes(tl); pax.Layout.Tile = 4; hold(pax, 'on');
            if any(~MC_Fail)
                polarscatter(pax, th(~MC_Fail), W(~MC_Fail), 25, passColor, 'filled', 'MarkerFaceAlpha', 0.6);
            end
            if any(MC_Fail)
                polarscatter(pax, th(MC_Fail), W(MC_Fail), 25, failColor, 'filled', 'MarkerFaceAlpha', 0.8);
            end
            pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
            title(pax, 'Sampled Direction (0 = meridional+) vs Speed [m/s]');
        catch
            ax = nexttile(tl, 4);
            plotPassFail(ax, th * 180/pi, W, MC_Fail, 'Direction [deg]', 'Mean wind [m/s]', ...
                'Wind Direction vs Speed', NaN, clr);
        end
    end
else
    disp('Wind_Mean_all not found -- skipping the wind envelope (V4 file?).');
end

%% 9. Failure Drivers
% Every sampled disturbance, ranked two ways:
%   KS distance  : how differently the input is distributed in failed vs
%                  passed runs (0 = identical, 1 = fully separated)
%   Spearman rho : rank correlation of the input with the max errors
drv = struct('name', {}, 'x', {});
drv = addDriver(drv, 'Mean wind [m/s]',          getField(S, 'Wind_Mean_all', []), num_sims);
Wm = reshape(getField(S, 'Wind_Mean_all', nan(num_sims,1)), [], 1);
Wg = reshape(getField(S, 'Wind_Gain_all', nan(num_sims,1)), [], 1);
if isfield(S, 'dryden_W6') && numel(Wm) == num_sims && numel(Wg) == num_sims
    tm = Wg * S.dryden_W6 ./ Wm; tm(Wm <= 0) = NaN;
    drv = addDriver(drv, 'Turbulence multiplier [-]', tm, num_sims);
else
    drv = addDriver(drv, 'Wind gain [-]', Wg, num_sims);
end
drv = addDriver(drv, 'GMP covariance [-]',        getField(S, 'Wind_Covar_all', []), num_sims);
drv = addDriver(drv, 'Radial lever [m]',          getField(S, 'Lever_Radial', []), num_sims);
drv = addDriver(drv, 'Axial lever [m]',           getField(S, 'Lever_Axial', []), num_sims);
if isfield(S, 'J_d_vals') && iscell(S.J_d_vals) && numel(S.J_d_vals) == num_sims
    drv = addDriver(drv, 'dJxx [kg m^2]', cellfun(@(M) M(1,1), S.J_d_vals), num_sims);
    drv = addDriver(drv, 'dJyy [kg m^2]', cellfun(@(M) M(2,2), S.J_d_vals), num_sims);
    drv = addDriver(drv, 'dJzz [kg m^2]', cellfun(@(M) M(3,3), S.J_d_vals), num_sims);
end
drv = addDriver(drv, 'Accel bias |b| [m/s^2]',    getField(S, 'AccelBias_Mag', []), num_sims);
drv = addDriver(drv, 'Gyro bias |b| [deg/s]',     getField(S, 'GyroBias_Mag', []) * 180/pi, num_sims);
drv = addDriver(drv, 'Hard iron |b| [-]',         getField(S, 'HardIron_Mag', []), num_sims);
drv = addDriver(drv, 'Soft iron deviation [-]',   getField(S, 'SoftIron_Dev', []), num_sims);
drv = addDriver(drv, 'Gimbal bias [deg]',         getField(S, 'GimbalBias_Mag', []) * 180/pi, num_sims);
drv = addDriver(drv, 'Thrust bias [N]',           getField(S, 'ThrustBias_all', []), num_sims);
if isfield(S, 'IMUDriftParams_vals') && iscell(S.IMUDriftParams_vals) ...
        && numel(S.IMUDriftParams_vals) == num_sims ...
        && all(cellfun(@numel, S.IMUDriftParams_vals) >= 8)
    A = cell2mat(cellfun(@(v) reshape(v(1:8), 1, []), S.IMUDriftParams_vals(:), 'UniformOutput', false));
    drv = addDriver(drv, 'Accel drift |\mu| [m/s^2/s]',   sqrt(sum(A(:,1:3).^2, 2)), num_sims);
    drv = addDriver(drv, 'Accel drift \sigma',            A(:,4), num_sims);
    drv = addDriver(drv, 'Gyro drift |\mu| [deg/s/s]',    sqrt(sum(A(:,5:7).^2, 2)) * 180/pi, num_sims);
    drv = addDriver(drv, 'Gyro drift \sigma',             A(:,8) * 180/pi, num_sims);
end

nd = numel(drv);
ks   = nan(nd, 1);
rhoP = nan(nd, 1);
rhoA = nan(nd, 1);
haveBoth = any(MC_Fail) && any(~MC_Fail);
for d = 1:nd
    if haveBoth, ks(d) = ksDistance(drv(d).x(MC_Fail), drv(d).x(~MC_Fail)); end
    rhoP(d) = spearmanRho(drv(d).x, Pos_Err_Max);
    rhoA(d) = spearmanRho(drv(d).x, Att_Err_Max);
end

if nd > 0
    names = {drv.name};
    figure('Name', 'Failure Drivers', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    if haveBoth
        plotRanking(ax, ks, names, failColor, false);
        xlabel(ax, 'KS distance (failed vs passed)'); xlim(ax, [0 1]);
    else
        axis(ax, 'off');
        text(ax, 0.5, 0.5, sprintf('KS needs both failed and passed runs\n(%d failed of %d)', nFail, num_sims), ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
    title(ax, 'What Failed Runs Have in Common');

    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plotRanking(ax, rhoP, names, axColors{1}, true);
    xlabel(ax, 'Spearman \rho with max |pos err|'); xlim(ax, [-1 1]);
    title(ax, 'Drivers of Position Error');

    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    plotRanking(ax, rhoA, names, axColors{2}, true);
    xlabel(ax, 'Spearman \rho with max att err'); xlim(ax, [-1 1]);
    title(ax, 'Drivers of Attitude Error');

    % Top drivers in detail
    if haveBoth, score = ks; else, score = max(abs(rhoP), abs(rhoA)); end
    score(~isfinite(score)) = -inf;
    [~, ord] = sort(score, 'descend');
    nTop = min(6, nd);
    figure('Name', 'Top Driver Scatters', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    for k = 1:nTop
        d = ord(k);
        ax = nexttile;
        if abs(rhoA(d)) > abs(rhoP(d))
            plotPassFail(ax, drv(d).x, Att_Err_Max, MC_Fail, drv(d).name, 'Max att err [deg]', ...
                sprintf('#%d  (\\rho = %.2f)', k, rhoA(d)), att_err_fail, clr);
        else
            plotPassFail(ax, drv(d).x, Pos_Err_Max, MC_Fail, drv(d).name, 'Max |pos err| [m]', ...
                sprintf('#%d  (\\rho = %.2f)', k, rhoP(d)), pos_err_fail, clr);
        end
        if k == 1, legend(ax, 'Location', 'best'); end
    end
end

%% 10. Sensor Disturbance Diagnostics
% Accelerometer errors act through tilt -> horizontal drift (position);
% gyro and magnetometer errors act on the attitude estimate (attitude).
sensorPairs = {'Accel bias |b| [m/s^2]',      Pos_Err_Max, 'Max |pos err| [m]', pos_err_fail;
               'Accel drift |\mu| [m/s^2/s]', Pos_Err_Max, 'Max |pos err| [m]', pos_err_fail;
               'Gyro bias |b| [deg/s]',       Att_Err_Max, 'Max att err [deg]', att_err_fail;
               'Gyro drift |\mu| [deg/s/s]',  Att_Err_Max, 'Max att err [deg]', att_err_fail;
               'Hard iron |b| [-]',           Att_Err_Max, 'Max att err [deg]', att_err_fail;
               'Soft iron deviation [-]',     Att_Err_Max, 'Max att err [deg]', att_err_fail};
plotDiagGrid('Sensor Disturbance Diagnostics', sensorPairs, drv, MC_Fail, clr, bkgColor);

%% 11. Vehicle / Actuator Disturbance Diagnostics
vehPairs = {'Gimbal bias [deg]',  Pos_Err_Max,  'Max |pos err| [m]',   pos_err_fail;
            'Gimbal bias [deg]',  Att_Err_Max,  'Max att err [deg]',   att_err_fail;
            'Thrust bias [N]',    Vert_Err_Max, 'Max |vert err| [m]',  NaN;
            'Radial lever [m]',   Att_Err_Max,  'Max att err [deg]',   att_err_fail;
            'dJxx [kg m^2]',      Att_Err_Max,  'Max att err [deg]',   att_err_fail;
            'dJzz [kg m^2]',      Att_Err_Max,  'Max att err [deg]',   att_err_fail};
plotDiagGrid('Vehicle & Actuator Diagnostics', vehPairs, drv, MC_Fail, clr, bkgColor);

%% 12. Filter (and Control) RMSE Distributions
if isfield(S, 'RMSE_Filter_all') && size(S.RMSE_Filter_all, 1) >= 3
    RF = S.RMSE_Filter_all;
    if any(imag(RF(:)) ~= 0)
        warning('PlotMC_v5:ComplexRMSE', ['RMSE_Filter_all is complex: SSE_Filter went negative, ' ...
            'so it is not integrating a squared error. Plotting the real part only.']);
    end
    figure('Name', 'Filter Attitude RMSE', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    filter_lbls = {'Filter RMSE 1', 'Filter RMSE 2', 'Filter RMSE 3'};
    for c = 1:3
        ax = nexttile;
        plotSmartHistogram(ax, real(RF(c,:)), axColors{c}, filter_lbls{c});
        title(ax, filter_lbls{c}); xlabel(ax, 'RMSE (as logged)'); grid(ax, 'on'); legend(ax, 'show');
    end
end

if plot_control_rmse && isfield(S, 'RMSE_Controls_all') && size(S.RMSE_Controls_all, 1) == 12
    labels = {'Attitude', 'Angular Rate', 'Position', 'Velocity'};
    cunits = {'deg', 'deg/s', 'm', 'm/s'};
    multipliers = [(180/pi), (180/pi), 1, 1];
    for g = 1:4
        figure('Name', sprintf('%s Control RMSE', labels{g}), 'Color', bkgColor, 'WindowStyle', 'docked');
        tiledlayout(3, 1, 'TileSpacing', 'compact');
        for c = 1:3
            ax = nexttile;
            plotSmartHistogram(ax, real(S.RMSE_Controls_all((g-1)*3 + c, :)) * multipliers(g), ...
                axColors{c}, axis_lbls{c});
            title(ax, sprintf('%s - %s', labels{g}, axis_lbls{c}));
            xlabel(ax, sprintf('RMSE (%s)', cunits{g})); grid(ax, 'on'); legend(ax, 'show');
        end
    end
end

%% 13. Landing Angular-Rate Distribution
if ~isempty(rates_all)
    landing_rate = nan(num_sims, 3);
    for i = 1:num_sims
        rates = reshape(rates_all(i,:,:), T, 3) * (180/pi);
        z = pos_all(i,:,3)';
        validRows = all(isfinite(rates), 2) & isfinite(z);
        if ~any(validRows), continue; end
        z0 = z(find(validRows, 1, 'first'));            % pad height
        above = (z - z0) > 0.05 & validRows;
        crossings = find(diff(above) == -1);
        if ~isempty(crossings)
            land_idx = crossings(end) + 1;
        else
            land_idx = find(validRows, 1, 'last');
        end
        landing_rate(i,:) = rates(land_idx,:);
    end

    figure('Name', 'Landing Rates', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(3, 1, 'TileSpacing', 'compact');
    rate_lbls = {'Rate X', 'Rate Y', 'Rate Z'};
    for c = 1:3
        ax = nexttile;
        plotSmartHistogram(ax, landing_rate(:,c), axColors{c}, rate_lbls{c});
        title(ax, sprintf('Landing %s', rate_lbls{c})); xlabel(ax, 'Rate (deg/s)');
        grid(ax, 'on'); legend(ax, 'show');
    end
end

%% 14. Spectral Radius Sensitivities
if isfield(S, 'MaxSpectralRad_all') && size(S.MaxSpectralRad_all, 1) == 2
    figure('Name', 'Spectral Radius Sensitivities', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(1, 2, 'TileSpacing', 'compact');
    for c = 1:2
        ax = nexttile;
        plotPassFail(ax, S.MaxSpectralRad_all(c,:), Pos_Err_Max, MC_Fail, sprintf('Max MSR(%d)', c), ...
            'Max |pos err| [m]', sprintf('Position Error vs MSR(%d)', c), pos_err_fail, clr);
        xline(ax, 1, 'k--', 'HandleVisibility', 'off');
        if c == 1, legend(ax, 'Location', 'best'); end
    end
end

%% 15. Spectral Radius Timeseries Overlay & Causality Check
baseWin = [pad_hold_time + 2, pad_hold_time + 8];      % settled, after the pad hold
if isfield(S, 'specrad_all') && ~isempty(S.specrad_all) && size(S.specrad_all, 3) >= 2
    specrad_all = S.specrad_all;
    [worst_idx, clean_idx] = pickWorstClean(Pos_Err_Max, 8, 4);
    plot_set   = [worst_idx(:); clean_idx(:)];
    run_labels = [repmat({'WORST'}, numel(worst_idx), 1); repmat({'clean'}, numel(clean_idx), 1)];

    if ~isempty(plot_set)
        figure('Name', 'MSR(2) vs Position Error Overlay', 'Color', bkgColor, 'WindowStyle', 'docked');
        nCols = 4; nRows = ceil(numel(plot_set) / nCols);
        tlOverlay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:numel(plot_set)
            i = plot_set(k);
            ax1 = nexttile(tlOverlay); hold(ax1, 'on'); grid(ax1, 'on');
            yyaxis(ax1, 'left');
            plot(ax1, t_common, reshape(specrad_all(i,:,2), T, 1), 'b-', 'LineWidth', 1.3);
            yline(ax1, 1, 'b--');
            ylabel(ax1, 'MSR(2)');
            yyaxis(ax1, 'right');
            plot(ax1, t_common, pos_err_norm(i,:), 'r-', 'LineWidth', 1.3);
            ylabel(ax1, '||pos err|| [m]');
            title(ax1, sprintf('Run %d (%s), max err %.2f m', i, run_labels{k}, Pos_Err_Max(i)));
            xlabel(ax1, 'Time [s]');
        end
        title(tlOverlay, 'Blue = MSR(2) (left, ref line at 1)   Red = position error norm (right)');
    end

    msr2 = specrad_all(:,:,2);
    plotLagAnalysis(t_common, msr2, pos_err_norm, baseWin, @(sig, ~) 1, pos_err_fail, ...
        'MSR(2)', 'crosses 1', bkgColor);
else
    disp('specrad_all not found -- skipping the spectral radius timeseries overlay.');
end

%% 16. LESO Disturbance Causality: Accel & Torque vs Divergence
if isfield(S, 'LESODist_all') && ~isempty(S.LESODist_all) && size(S.LESODist_all, 3) == 6
    L = S.LESODist_all;
    accel_dist_norm  = sqrt(sum(L(:,:,1:3).^2, 3));    % N x T, m/s^2
    torque_dist_norm = sqrt(sum(L(:,:,4:6).^2, 3));    % N x T, N*m

    plotLESOCausality(t_common, accel_dist_norm, pos_err_norm, Pos_Err_Max, baseWin, ...
        bkgColor, 'LESO Accel Disturbance', '||pos err|| [m]', 'Accel Dist [m/s^2]', pos_err_fail);
    plotLESOCausality(t_common, torque_dist_norm, att_err, Att_Err_Max, baseWin, ...
        bkgColor, 'LESO Torque Disturbance', 'att err [deg]', 'Torque Dist [N\cdotm]', att_err_fail);
else
    disp('LESODist_all not found or malformed -- skipping LESO causality diagnostics.');
end

end

%% ======================= Local Helpers =======================
function v = getField(S, name, default)
    if isfield(S, name) && ~isempty(S.(name))
        v = S.(name);
    else
        v = default;
    end
end

function q = normQuat(q)
    q = q ./ sqrt(sum(q.^2, 3));
end

function tilt = tiltFromQuat(q)
    % Angle between body z and inertial z: C33 = q0^2 - q1^2 - q2^2 + q3^2
    c = q(:,:,1).^2 - q(:,:,2).^2 - q(:,:,3).^2 + q(:,:,4).^2;
    c(c > 1) = 1; c(c < -1) = -1;                  % keeps NaN as NaN
    tilt = acosd(c);
end

function data = extractLoggedState(ts)
    if isempty(ts), data = []; return; end
    data = squeeze(ts.Data);
    if size(data, 1) ~= numel(ts.Time), data = data'; end
end

function [lo, hi] = clopperPearson(k, n, conf)
    % Exact binomial interval for a failure rate (base MATLAB only).
    a = (1 - conf) / 2;
    lo = zeros(size(k)); hi = ones(size(k));
    for j = 1:numel(k)
        if n(j) == 0, lo(j) = NaN; hi(j) = NaN; continue; end
        if k(j) > 0,    lo(j) = betaincinv(a,     k(j),     n(j) - k(j) + 1); end
        if k(j) < n(j), hi(j) = betaincinv(1 - a, k(j) + 1, n(j) - k(j));     end
    end
end

function P = pctlCols(X, p)
    % Percentiles of each column of X (N x T), ignoring NaN.
    T = size(X, 2);
    P = nan(numel(p), T);
    Xs = sort(X, 1);                               % NaN sorts to the end
    n = sum(~isnan(X), 1);
    for c = 1:T
        if n(c) > 0, P(:, c) = pctlSorted(Xs(1:n(c), c), p); end
    end
end

function v = pctl1(x, p)
    x = sort(x(~isnan(x)));
    if isempty(x), v = nan(numel(p), 1); else, v = pctlSorted(x(:), p); end
end

function v = pctlSorted(x, p)
    % Same interpolation as prctile's default method.
    n = numel(x);
    if n == 1, v = repmat(x, numel(p), 1); return; end
    pos = n * p(:) / 100 + 0.5;
    pos = min(max(pos, 1), n);
    lo = floor(pos); hi = ceil(pos); f = pos - lo;
    v = (1 - f) .* x(lo) + f .* x(hi);
end

function plotBands(ax, t, X, clr)
    % Median, 5-95% and min-max across runs (rows of X).
    P = pctlCols(X, [0 5 50 95 100]);
    good = all(isfinite(P), 1);
    if ~any(good), return; end
    tg = t(good)'; P = P(:, good);
    fill(ax, [tg, fliplr(tg)], [P(1,:), fliplr(P(5,:))], clr, 'FaceAlpha', 0.10, ...
        'EdgeColor', 'none', 'DisplayName', 'Min-max');
    fill(ax, [tg, fliplr(tg)], [P(2,:), fliplr(P(4,:))], clr, 'FaceAlpha', 0.30, ...
        'EdgeColor', 'none', 'DisplayName', '5-95%');
    plot(ax, tg, P(3,:), '-', 'Color', clr, 'LineWidth', 1.5, 'DisplayName', 'Median');
end

function plotFailedRuns(ax, t, X, fail, failColor)
    idx = find(fail(:) & any(isfinite(X), 2));
    for j = 1:numel(idx)
        h = plot(ax, t, X(idx(j), :), '-', 'Color', [failColor 0.45], 'LineWidth', 0.7);
        if j == 1
            set(h, 'DisplayName', 'Failed runs');
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

function shadePadHold(ax, t0, t1)
    if t1 <= t0, return; end
    yl = ylim(ax);
    h = patch(ax, [t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [0.85 0.85 0.85], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
    uistack(h, 'bottom');
    ylim(ax, yl);
end

function plotPassFail(ax, x, y, fail, xlbl, ylbl, ttl, thr, clr)
    x = x(:); y = y(:); fail = fail(:);
    hold(ax, 'on'); grid(ax, 'on');
    good = isfinite(x) & isfinite(y);
    if any(good & ~fail)
        scatter(ax, x(good & ~fail), y(good & ~fail), 28, clr.pass, 'filled', ...
            'MarkerFaceAlpha', 0.55, 'MarkerEdgeColor', 'none', 'DisplayName', 'Pass');
    end
    if any(good & fail)
        scatter(ax, x(good & fail), y(good & fail), 28, clr.fail, 'filled', ...
            'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none', 'DisplayName', 'Fail');
    end
    if isfinite(thr)
        yline(ax, thr, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Limit');
    end
    xlabel(ax, xlbl); ylabel(ax, ylbl); title(ax, ttl);
end

function drv = addDriver(drv, name, x, N)
    % Keep an input only if it has one value per run and actually varies.
    x = double(x(:));
    if numel(x) ~= N, return; end
    xf = x(isfinite(x));
    if numel(xf) < 2 || max(xf) - min(xf) <= 0, return; end
    drv(end+1) = struct('name', name, 'x', x);
end

function x = findDriver(drv, name)
    x = [];
    k = find(strcmp({drv.name}, name), 1);
    if ~isempty(k), x = drv(k).x; end
end

function plotRanking(ax, vals, names, clr, signed)
    v = vals(:);
    keep = isfinite(v);
    v = v(keep); nm = names(keep);
    if isempty(v)
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'Not enough data', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        return;
    end
    if signed, [~, o] = sort(abs(v), 'ascend'); else, [~, o] = sort(v, 'ascend'); end
    barh(ax, 1:numel(v), v(o), 'FaceColor', clr, 'FaceAlpha', 0.85);   % largest at the top
    set(ax, 'YTick', 1:numel(v), 'YTickLabel', nm(o));
    ylim(ax, [0.4, numel(v) + 0.6]);
end

function plotDiagGrid(figName, pairs, drv, fail, clr, bkgColor)
    have = false(size(pairs, 1), 1);
    for r = 1:size(pairs, 1), have(r) = ~isempty(findDriver(drv, pairs{r,1})); end
    if ~any(have)
        fprintf('%s: none of its inputs were dispersed -- skipped.\n', figName);
        return;
    end
    rows = find(have);
    figure('Name', figName, 'Color', bkgColor, 'WindowStyle', 'docked');
    nC = min(3, numel(rows)); nR = ceil(numel(rows) / nC);
    tiledlayout(nR, nC, 'TileSpacing', 'compact', 'Padding', 'compact');
    for j = 1:numel(rows)
        r = rows(j);
        x = findDriver(drv, pairs{r,1});
        y = pairs{r,2};
        ax = nexttile;
        plotPassFail(ax, x, y, fail, pairs{r,1}, pairs{r,3}, ...
            sprintf('%s  (\\rho = %.2f)', pairs{r,3}, spearmanRho(x, y)), pairs{r,4}, clr);
        if j == 1, legend(ax, 'Location', 'best'); end
    end
end

function r = tiedRank(v)
    v = v(:); n = numel(v);
    [s, idx] = sort(v);
    r = zeros(n, 1);
    k = 1;
    while k <= n
        m = k;
        while m < n && s(m+1) == s(k), m = m + 1; end
        r(idx(k:m)) = (k + m) / 2;
        k = m + 1;
    end
end

function rho = spearmanRho(x, y)
    x = x(:); y = y(:);
    good = isfinite(x) & isfinite(y);
    rho = NaN;
    if nnz(good) < 3, return; end
    rx = tiedRank(x(good)); ry = tiedRank(y(good));
    rx = rx - mean(rx); ry = ry - mean(ry);
    den = sqrt(sum(rx.^2) * sum(ry.^2));
    if den > 0, rho = sum(rx .* ry) / den; end
end

function d = ksDistance(a, b)
    % Two-sample Kolmogorov-Smirnov statistic: max gap between the CDFs.
    a = a(~isnan(a)); b = b(~isnan(b));
    d = NaN;
    if isempty(a) || isempty(b), return; end
    v = sort([a(:); b(:)]);
    Fa = sum(a(:)' <= v, 2) / numel(a);
    Fb = sum(b(:)' <= v, 2) / numel(b);
    d = max(abs(Fa - Fb));
end

function [worst, clean] = pickWorstClean(metric, nWorst, nClean)
    % NaN-safe: sort('descend') would otherwise put NaN (crashed) runs first.
    m = metric(:);
    fin = find(isfinite(m));
    [~, o] = sort(m(fin), 'descend');
    worst = fin(o(1:min(nWorst, numel(o))));
    rest = setdiff(fin, worst, 'stable');
    [~, o2] = sort(m(rest), 'ascend');
    clean = rest(o2(1:min(nClean, numel(o2))));
end

function plotSmartHistogram(ax, data, clr, name)
    data = data(isfinite(data));
    if isempty(data), return; end
    hold(ax, 'on');

    % Bin the 2nd-98th percentile core at full resolution; collapse the
    % tails into single labeled overflow/underflow bins instead of
    % stretching the axis to fit a handful of outlier runs.
    b = pctl1(data, [2 98]);
    loB = b(1); hiB = b(2);
    if hiB - loB < eps
        hiB = loB + max(abs(loB) * 0.05, eps);
    end

    nBins = max(10, min(40, round(sqrt(numel(data)))));
    edges = linspace(loB, hiB, nBins + 1);
    binW  = edges(2) - edges(1);
    centers = edges(1:end-1) + binW / 2;

    counts    = histcounts(data, edges);
    underflow = sum(data < loB);
    overflow  = sum(data > hiB);

    bar(ax, centers, counts, 1, 'FaceColor', clr, 'FaceAlpha', 0.85, 'EdgeColor', 'k', 'DisplayName', name);

    if underflow > 0
        bar(ax, loB - binW, underflow, binW, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, loB - binW, underflow, num2str(underflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    if overflow > 0
        bar(ax, hiB + binW, overflow, binW, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, hiB + binW, overflow, num2str(overflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end

    % Log-scale the y-axis when one bin dwarfs the rest, so smaller bars
    % stay visible instead of flattening to nothing next to the spike.
    allCounts = [counts, underflow, overflow];
    peak = max(allCounts);
    rest = allCounts(allCounts < peak & allCounts > 0);
    if ~isempty(rest) && peak > 4 * max(rest)
        set(ax, 'YScale', 'log');
        ylim(ax, [0.8, peak * 1.3]);
    end
end

function plotLagAnalysis(t, sig_all, err_all, baseWin, sigThreshFn, err_abs_floor, sigName, eventWord, bkgColor)
    % For every run: time from the indicator signal's event (sigThreshFn) to
    % the error crossing max(err_abs_floor, 5 x its settled-window std).
    % Positive lag = the indicator leads divergence.
    base = t >= baseWin(1) & t <= baseWin(2);
    lags = [];
    crossed_never_diverged = 0;
    diverged_never_crossed = 0;
    for i = 1:size(sig_all, 1)
        sig  = sig_all(i,:);
        perr = err_all(i,:);
        if all(isnan(sig)) || all(isnan(perr)), continue; end

        sig_thresh = sigThreshFn(sig, base);
        err_thresh = max(err_abs_floor, 5 * std(perr(base), 'omitnan'));

        cross_i = find(sig > sig_thresh, 1, 'first');
        div_i   = find(perr > err_thresh, 1, 'first');

        if ~isempty(cross_i) && ~isempty(div_i)
            lags(end+1) = t(div_i) - t(cross_i); %#ok<AGROW>
        elseif ~isempty(cross_i)
            crossed_never_diverged = crossed_never_diverged + 1;
        elseif ~isempty(div_i)
            diverged_never_crossed = diverged_never_crossed + 1;
        end
    end

    figure('Name', sprintf('%s -> Divergence Lag', sigName), 'Color', bkgColor, 'WindowStyle', 'docked');
    hold on; grid on;
    if ~isempty(lags)
        histogram(lags, 'FaceColor', [0.3 0.5 0.8]);
    end
    xline(0, 'k--', 'LineWidth', 1.5);
    xlabel(sprintf('t_{divergence} - t_{%s %s}  [s]', sigName, eventWord));
    ylabel('Count');
    title(sprintf(['Lag distribution (n=%d runs with both events)\n' ...
        '%d %s but never diverged | %d diverged without it'], ...
        numel(lags), crossed_never_diverged, eventWord, diverged_never_crossed));

    fprintf('\n--- %s %s -> divergence timing ---\n', sigName, eventWord);
    fprintf('Runs with both events: %d\n', numel(lags));
    if ~isempty(lags)
        fprintf('  median lag = %.2f s, mean lag = %.2f s\n', median(lags), mean(lags));
        fprintf('  fraction with POSITIVE lag (%s leads divergence): %.0f%%\n', sigName, 100 * mean(lags > 0));
    end
    fprintf('%s but NEVER diverged: %d (harmless bumps)\n', eventWord, crossed_never_diverged);
    fprintf('Diverged WITHOUT it: %d (another cause)\n', diverged_never_crossed);
end

function plotLESOCausality(t, sig_all, err_all, sortMetric, baseWin, bkgColor, sigName, errYLabel, sigYLabel, err_abs_floor)
    [worst_idx, clean_idx] = pickWorstClean(sortMetric, 8, 4);
    plot_set   = [worst_idx(:); clean_idx(:)];
    run_labels = [repmat({'WORST'}, numel(worst_idx), 1); repmat({'clean'}, numel(clean_idx), 1)];

    if ~isempty(plot_set)
        figure('Name', sprintf('%s vs %s Overlay', sigName, errYLabel), 'Color', bkgColor, 'WindowStyle', 'docked');
        nCols = 4; nRows = ceil(numel(plot_set) / nCols);
        tlOverlay = tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');
        for k = 1:numel(plot_set)
            i = plot_set(k);
            ax1 = nexttile(tlOverlay); hold(ax1, 'on'); grid(ax1, 'on');
            yyaxis(ax1, 'left');
            plot(ax1, t, sig_all(i,:), 'b-', 'LineWidth', 1.3);
            ylabel(ax1, sigYLabel);
            yyaxis(ax1, 'right');
            plot(ax1, t, err_all(i,:), 'r-', 'LineWidth', 1.3);
            ylabel(ax1, errYLabel);
            title(ax1, sprintf('Run %d (%s), max %.2f', i, run_labels{k}, sortMetric(i)));
            xlabel(ax1, 'Time [s]');
        end
        title(tlOverlay, sprintf('Blue = %s (left)   Red = %s (right)', sigName, errYLabel));
    end

    % No fixed stability boundary for the LESO estimate (unlike MSR = 1),
    % so its threshold is adaptive per run: baseline mean + 5 sigma over the
    % settled window after the pad hold.
    threshFn = @(sig, base) mean(sig(base), 'omitnan') + 5 * std(sig(base), 'omitnan');
    plotLagAnalysis(t, sig_all, err_all, baseWin, threshFn, err_abs_floor, sigName, 'elevated', bkgColor);
end

function drawMCTraj_v5(ax, pos_all, target_pos, x_idx, y_idx, x_lbl, y_lbl, alphaVal, fail, failColor)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    T = size(pos_all, 2);
    for pass = [false true]                    % passed first, failed on top
        idx = find(fail(:) == pass)';
        if pass, c = [failColor, min(1, alphaVal * 3)]; else, c = [0.25 0.25 0.25 alphaVal]; end
        for i = idx
            xy = reshape(pos_all(i, :, [x_idx, y_idx]), T, 2);
            valid = all(isfinite(xy), 2);
            plot(ax, xy(valid,1), xy(valid,2), 'Color', c, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        end
    end
    plot(ax, target_pos(:, x_idx), target_pos(:, y_idx), 'b--', 'LineWidth', 2);
    xlabel(ax, x_lbl); ylabel(ax, y_lbl);
end