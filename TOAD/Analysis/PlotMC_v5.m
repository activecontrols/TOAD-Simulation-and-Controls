function PlotMC_v5(filename, withinValidWind, windLimit)
% PlotMC_v5
% Monte Carlo plotting for the V5 framework (built around per-run
% pass/fail, physical wind envelope and verified plant disturbances).
% Compatible with TOAD_Simulation_CPP To Workspace log blocks.
%
%   PlotMC_v5                        % after MonteCarlo_v5: reads the base workspace
%   PlotMC_v5('MC5_...mat')          % reads a saved run
%   PlotMC_v5('MC5_...mat', true, 6) % optionally restrict to mean wind < 6 m/s
%
% Figures
%   1  Summary: failure rate (90% Clopper-Pearson), failure modes, max errors
%   2  Trajectories (passed grey, failed red) + projections
%   3  Position / velocity distributions
%   4  Error histories against the pass/fail thresholds
%   5  Tilt and angular rates
%   6  Wind envelope: failure rate vs mean wind, wind direction
%   7  Failure drivers: ranking + top-driver scatters
%   8  Sensor disturbance diagnostics (Accel bias, Gyro bias, Soft iron)
%   9  Wind analysis: 3-axis wind timeseries + wind-speed error/gimbal correlations
%   10 Controls: gimbal command timeseries + frequency PSD (oscillatory modes)
%   11 Landing angular rates

close all;
if nargin < 2 || isempty(withinValidWind)
    withinValidWind = false;   % true: plot & evaluate only runs with mean wind < windLimit
end
if nargin < 3 || isempty(windLimit)
    windLimit       = 6;       % m/s, upper mean-wind cutoff when withinValidWind is true
end

bkgColor = 'w';
alphaVal = 0.12;

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
varList = {'t_common', 't_u', 'pos_all', 'vel_all', 'quat_all', 'rates_all', 'wind_all', 'u_cmd_all', ...
           'target_state_all', 'ok', ...
           'seed', 'pad_hold_time', 'pos_err_fail', 'att_err_fail', ...
           'wind_sweep', 'dryden_W6', ...
           'Wind_Mean_all', 'Wind_Dir_all', 'Wind_Gain_all', 'Wind_Covar_all', ...
           'Lever_Radial', 'Lever_Axial', 'J_d_vals', ...
           'AccelBias_Mag', 'GyroBias_Mag', 'SoftIron_Dev', ...
           'RMSE_Filter_all', 'RMSE_Wind_all', 'out'};
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
    avail = {whos('-file', filename).name};
    toLoad = varList(ismember(varList, avail));
    S = load(filename, toLoad{:});
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
num_sims_total = size(S.pos_all, 1);

if withinValidWind
    if ~isfield(S, 'Wind_Mean_all') || numel(S.Wind_Mean_all) ~= num_sims_total
        error('PlotMC_v5:MissingWind', ...
            'Wind_Mean_all is required when withinValidWind is true.');
    end
    W_raw = reshape(S.Wind_Mean_all, [], 1);
    keep  = isfinite(W_raw) & (W_raw < windLimit);
    if ~any(keep)
        error('PlotMC_v5:NoRunsWithinWind', ...
            'No simulations have mean wind below %.2f m/s (min sampled: %.2f m/s).', ...
            windLimit, min(W_raw));
    end
    S = filterStructRuns(S, keep, num_sims_total);
end

pos_all  = S.pos_all;
num_sims = size(pos_all, 1);
quat_all = S.quat_all;
tgt      = S.target_state_all;
vel_all  = getField(S, 'vel_all', nan(num_sims, T, 3));
ok       = logical(reshape(getField(S, 'ok', true(num_sims, 1)), [], 1));

pad_hold_time = getField(S, 'pad_hold_time', 15);
pos_err_fallback = getField(S, 'pos_err_fail', 2.0);
att_err_fallback = getField(S, 'att_err_fail', 15);

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
Horiz_Err_Max = max(horiz_err(:, eval_idx), [], 2); %#ok<NASGU>
Vert_Err_Max  = max(abs(vert_err(:, eval_idx)), [], 2); %#ok<NASGU>
Att_Err_Max   = max(att_err(:, eval_idx), [], 2);

Fail_Error    = ~ok;
Fail_NoData   = ~Fail_Error & (isnan(Pos_Err_Max) | isnan(Att_Err_Max));

validStats = ~Fail_Error & ~Fail_NoData;
[pos_err_fail, pos_mu, pos_sig] = twoSigmaLimit(Pos_Err_Max(validStats), pos_err_fallback);
[att_err_fail, att_mu, att_sig] = twoSigmaLimit(Att_Err_Max(validStats), att_err_fallback);

Fail_Position = Pos_Err_Max > pos_err_fail;
Fail_Attitude = Att_Err_Max > att_err_fail;
MC_Fail = Fail_Error | Fail_NoData | Fail_Position | Fail_Attitude;
nFail   = nnz(MC_Fail);
[ciLo, ciHi] = clopperPearson(nFail, num_sims, 0.90);

% First threshold violation after the pad hold
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

if withinValidWind
    windNote = sprintf(' (of %d, wind < %.1f m/s)', num_sims_total, windLimit);
else
    windNote = '';
end
fprintf('\n=== PlotMC_v5: %d runs%s, %d failed (%.1f%%), 90%% CI on P(fail) [%.3f, %.3f] ===\n', ...
    num_sims, windNote, nFail, 100 * nFail / num_sims, ciLo, ciHi);
fprintf('   sim error %d | no data %d | pos > %.2f m (mu=%.2f, 2sig=%.2f) %d | att > %.1f deg (mu=%.1f, 2sig=%.1f) %d\n', ...
    nnz(Fail_Error), nnz(Fail_NoData), pos_err_fail, pos_mu, 2 * pos_sig, nnz(Fail_Position), ...
    att_err_fail, att_mu, 2 * att_sig, nnz(Fail_Attitude));

%% 3. Summary (Figure 1)
figure('Name', 'V5 Summary', 'Color', bkgColor, 'WindowStyle', 'docked');
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); axis(ax, 'off');
if withinValidWind
    runsTxt = sprintf('Runs:          %d/%d (<%.1f m/s)', num_sims, num_sims_total, windLimit);
else
    runsTxt = sprintf('Runs:          %d', num_sims);
end
txt = {runsTxt, ...
       sprintf('Failed:        %d (%.1f%%)', nFail, 100 * nFail / num_sims), ...
       sprintf('P(fail) 90CI:  [%.3f, %.3f]', ciLo, ciHi), ...
       sprintf('Pos lim 2sig:  %.2f m', pos_err_fail), ...
       sprintf('Att lim 2sig:  %.1f deg', att_err_fail), ...
       sprintf('Eval after:    t = %.1f s', pad_hold_time)};
if isfield(S, 'seed')
    txt{end+1} = sprintf('Master seed:   %.0f', double(S.seed));
end
if isfield(S, 'wind_sweep')
    if S.wind_sweep
        txt{end+1} = 'Wind campaign: envelope';
    else
        txt{end+1} = 'Wind campaign: Weibull';
    end
end
if nFail == 0
    txt{end+1} = sprintf('0 fails: P < %.2f%% (95%%, rule3)', 300 / num_sims);
end
text(ax, 0, 1, txt, 'VerticalAlignment', 'top', 'FontName', 'FixedWidth', ...
    'Interpreter', 'none', 'FontSize', 8.5);
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
xline(ax, pos_err_fail, 'r--', 'LineWidth', 1.5, 'DisplayName', '2\sigma Limit');
title(ax, 'Max Position Error (after pad hold)'); xlabel(ax, 'm'); grid(ax, 'on'); legend(ax, 'show');

ax = nexttile(tl);
plotSmartHistogram(ax, Att_Err_Max, axColors{2}, 'Max att err');
xline(ax, att_err_fail, 'r--', 'LineWidth', 1.5, 'DisplayName', '2\sigma Limit');
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

%% 4. Kinematic Trajectory Overlay (3D + Projections - Figure 2)
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

%% 5. State Distributions (Percentile Bands - Figure 3)
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

%% 6. Error Histories vs Pass/Fail Thresholds (Figure 4)
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
        yline(ax, errThr(k), 'r--', 'LineWidth', 1.5, 'DisplayName', '2\sigma Limit');
    end
    shadePadHold(ax, t_common(1), pad_hold_time);
    title(ax, errTtl{k}); xlabel(ax, 'Time (s)'); ylabel(ax, errLbl{k});
    if k == 1, legend(ax, 'Location', 'best'); end
end

%% 7. Rotational Kinematics: Tilt Angle & Angular Rates (Figure 5)
rates_all = getField(S, 'rates_all', []);
if isempty(rates_all) && isfield(S, 'out') && ~isempty(S.out)
    rates_all = nan(num_sims, T, 3);
    for i = 1:min(num_sims, numel(S.out))
        if ~ok(i), continue; end
        ts = S.out(i).state_log;
        raw = extractLoggedState(ts);
        if size(raw, 2) >= 13
            for k = 1:3
                rates_all(i,:,k) = interp1(ts.Time, raw(:,10+k), t_common, 'linear', NaN);
            end
        end
    end
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
end

%% 8. Wind Envelope (Figure 6)
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
        jit = 0.03 * (mod((1:num_sims)' * 0.6180339887, 1) - 0.5);
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
        pax = polaraxes(tl); pax.Layout.Tile = 4; hold(pax, 'on');
        if any(~MC_Fail)
            polarscatter(pax, th(~MC_Fail), W(~MC_Fail), 25, passColor, 'filled', 'MarkerFaceAlpha', 0.6);
        end
        if any(MC_Fail)
            polarscatter(pax, th(MC_Fail), W(MC_Fail), 25, failColor, 'filled', 'MarkerFaceAlpha', 0.8);
        end
        pax.ThetaZeroLocation = 'top'; pax.ThetaDir = 'clockwise';
        title(pax, 'Sampled Direction (0 = meridional+) vs Speed [m/s]');
    end
end

%% 9. Failure Drivers (Figure 7 & Figure 8)
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
drv = addDriver(drv, 'Soft iron deviation [-]',   getField(S, 'SoftIron_Dev', []), num_sims);

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

%% 10. Sensor Disturbance Diagnostics (Figure 8 / Grid)
sensorPairs = {'Accel bias |b| [m/s^2]',    Pos_Err_Max, 'Max |pos err| [m]', pos_err_fail;
               'Gyro bias |b| [deg/s]',     Att_Err_Max, 'Max att err [deg]', att_err_fail;
               'Soft iron deviation [-]',   Att_Err_Max, 'Max att err [deg]', att_err_fail};
plotDiagGrid('Sensor Disturbance Diagnostics', sensorPairs, drv, MC_Fail, clr, bkgColor);

%% 11. Wind Time Series & Correlations (Figure 9)
wind_all = extractWindAll(S, t_common, num_sims, ok);
[u_cmd_all, t_u] = extractCmdAll(S, t_common, num_sims, ok);
if ~isempty(u_cmd_all)
    gimbal_deg = u_cmd_all(:, :, 1:2) * (180 / pi);
    u_eval_idx = t_u >= pad_hold_time;
    if ~any(u_eval_idx), u_eval_idx(:) = true; end
else
    gimbal_deg = [];
    u_eval_idx = [];
end

if ~isempty(wind_all) || isfield(S, 'Wind_Mean_all')
    if ~isempty(wind_all)
        w_mag = sqrt(sum(wind_all(:, eval_idx, :).^2, 3));
        W_spd = mean(w_mag, 2, 'omitnan');
        wSpdLbl = 'Mean Flight Wind Speed [m/s]';
    else
        W_spd = reshape(S.Wind_Mean_all, [], 1);
        wSpdLbl = 'Mean Wind [m/s]';
    end

    figure('Name', 'Wind Analysis', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
    if ~isempty(wind_all)
        plotOverlayBands(ax, t_common, wind_all, axColors, axis_lbls);
        shadePadHold(ax, t_common(1), pad_hold_time);
        xlabel(ax, 'Time (s)'); ylabel(ax, 'Wind Velocity (m/s)');
        title(ax, '3-Axis Wind Time Series (Median & 5-95%)');
        legend(ax, 'Location', 'best');
    else
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'Wind time series (SSE_Wind / wind_all) unavailable', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'Interpreter', 'none');
        title(ax, '3-Axis Wind Time Series');
    end

    ax = nexttile;
    plotPassFail(ax, W_spd, Pos_Err_Max, MC_Fail, wSpdLbl, 'Max |pos err| [m]', ...
        sprintf('Pos Error vs Wind Speed  (\\rho = %.2f)', spearmanRho(W_spd, Pos_Err_Max)), ...
        pos_err_fail, clr);
    legend(ax, 'Location', 'best');

    ax = nexttile;
    plotPassFail(ax, W_spd, Att_Err_Max, MC_Fail, wSpdLbl, 'Max att err [deg]', ...
        sprintf('Att Error vs Wind Speed  (\\rho = %.2f)', spearmanRho(W_spd, Att_Err_Max)), ...
        att_err_fail, clr);

    ax = nexttile;
    if ~isempty(gimbal_deg)
        Gimbal_Max = max(max(abs(gimbal_deg(:, u_eval_idx, :)), [], 3), [], 2);
        plotPassFail(ax, W_spd, Gimbal_Max, MC_Fail, wSpdLbl, 'Max |gimbal cmd| [deg]', ...
            sprintf('Gimbal Cmd vs Wind Speed  (\\rho = %.2f)', spearmanRho(W_spd, Gimbal_Max)), ...
            NaN, clr);
        yline(ax, 10, 'r--', 'LineWidth', 1.2, 'DisplayName', '10\circ Limit');
        legend(ax, 'Location', 'best');
    else
        axis(ax, 'off');
        text(ax, 0.5, 0.5, {'Gimbal commands not logged in this file.', ...
            'Enable input_log in MonteCarlo_v5 for correlation.'}, ...
            'HorizontalAlignment', 'center', 'Units', 'normalized', 'Interpreter', 'none');
        title(ax, 'Gimbal Cmd vs Wind Speed');
    end
end

%% 12. Controls: Gimbal Commands & Frequency Analysis (Figure 10)
if ~isempty(gimbal_deg) || ~isempty(rates_all)
    figure('Name', 'Gimbal & Frequency Analysis', 'Color', bkgColor, 'WindowStyle', 'docked');
    tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    if ~isempty(gimbal_deg)
        gim_lbls = {'Gimbal \theta (Ch 1)', 'Gimbal \phi (Ch 2)'};
        for c = 1:2
            ax = nexttile; hold(ax, 'on'); grid(ax, 'on');
            X = gimbal_deg(:, :, c);
            plotBands(ax, t_u, X, bandColor);
            plotFailedRuns(ax, t_u, X, MC_Fail, failColor);
            yline(ax, 10, 'r--', 'LineWidth', 1.2, 'DisplayName', '\pm10\circ Limit');
            yline(ax, -10, 'r--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            shadePadHold(ax, t_u(1), pad_hold_time);
            nSat = nnz(max(abs(X(:, u_eval_idx)), [], 2) >= 9.95);
            title(ax, sprintf('%s (%d/%d hit \\pm10\\circ)', gim_lbls{c}, nSat, num_sims));
            xlabel(ax, 'Time (s)'); ylabel(ax, 'Command (deg)');
            if c == 1, legend(ax, 'Location', 'best'); end
        end

        ax = nexttile;
        plotPSDGroup(ax, t_u, gimbal_deg, u_eval_idx, [1 2], {'\theta cmd', '\phi cmd'}, ...
            {oranX, blueY}, MC_Fail, failColor, 'deg', 'Gimbal Command PSD (Active Flight)');
    else
        for k = 1:3
            ax = nexttile; axis(ax, 'off');
            if k <= 2
                ttl = sprintf('Gimbal Ch %d Command', k);
            else
                ttl = 'Gimbal Command PSD';
            end
            text(ax, 0.5, 0.5, {'Gimbal commands (input_log / u_cmd_all) unavailable.', ...
                'Re-run MonteCarlo_v5 with input_log enabled.'}, ...
                'HorizontalAlignment', 'center', 'Units', 'normalized', 'Interpreter', 'none');
            title(ax, ttl);
        end
    end

    ax = nexttile;
    if ~isempty(rates_all)
        plotPSDGroup(ax, t_common, rates_all * (180/pi), eval_idx, [1 2 3], ...
            {'Rate X', 'Rate Y', 'Rate Z'}, axColors, MC_Fail, failColor, 'deg/s', ...
            'Body Angular Rate PSD (Active Flight)');
    else
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'Angular rate logs unavailable', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(ax, 'Body Angular Rate PSD');
    end
end

%% 13. Landing Angular-Rate Distribution (Figure 11)
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

end

%% ======================= Local Helpers =======================
function S = filterStructRuns(S, keep, N)
    % Subset all per-run arrays by logical mask keep without touching scalars/time.
    runFields = {'pos_all', 'vel_all', 'quat_all', 'rates_all', 'wind_all', 'u_cmd_all', ...
                 'target_state_all', 'ok', ...
                 'Wind_Mean_all', 'Wind_Dir_all', 'Wind_Gain_all', 'Wind_Covar_all', ...
                 'Lever_Radial', 'Lever_Axial', 'J_d_vals', ...
                 'AccelBias_Mag', 'GyroBias_Mag', 'SoftIron_Dev', 'RMSE_Wind_all', 'out'};
    for k = 1:numel(runFields)
        fn = runFields{k};
        if ~isfield(S, fn) || isempty(S.(fn)), continue; end
        v = S.(fn);
        if isvector(v) && numel(v) == N
            S.(fn) = v(keep);
        elseif size(v, 1) == N
            idx = repmat({':'}, 1, ndims(v));
            idx{1} = keep;
            S.(fn) = v(idx{:});
        end
    end
    if isfield(S, 'RMSE_Filter_all') && size(S.RMSE_Filter_all, 2) == N
        S.RMSE_Filter_all = S.RMSE_Filter_all(:, keep);
    end
end

function [lim, mu, sig] = twoSigmaLimit(vals, fallback)
    % Upper 2-sigma bound (mean + 2*std) on valid finite run-level peak errors.
    vals = vals(isfinite(vals));
    if numel(vals) >= 2
        mu  = mean(vals);
        sig = std(vals, 0);
        lim = mu + 2 * sig;
    elseif isscalar(vals)
        mu  = vals;
        sig = 0;
        lim = fallback;
    else
        mu  = NaN;
        sig = NaN;
        lim = fallback;
    end
end

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
    c(c > 1) = 1; c(c < -1) = -1;
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
    Xs = sort(X, 1);
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
        yline(ax, thr, 'r--', 'LineWidth', 1.2, 'DisplayName', '2\sigma Limit');
    end
    xlabel(ax, xlbl); ylabel(ax, ylbl); title(ax, ttl);
end

function drv = addDriver(drv, name, x, N)
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
    barh(ax, 1:numel(v), v(o), 'FaceColor', clr, 'FaceAlpha', 0.85);
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
    a = a(~isnan(a)); b = b(~isnan(b));
    d = NaN;
    if isempty(a) || isempty(b), return; end
    v = sort([a(:); b(:)]);
    Fa = sum(a(:)' <= v, 2) / numel(a);
    Fb = sum(b(:)' <= v, 2) / numel(b);
    d = max(abs(Fa - Fb));
end

function plotSmartHistogram(ax, data, clr, name)
    data = data(isfinite(data));
    if isempty(data), return; end
    hold(ax, 'on');

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
        bar(ax, [loB - binW, loB], [underflow, 0], 1, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, loB - binW, underflow, num2str(underflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
    if overflow > 0
        bar(ax, [hiB, hiB + binW], [0, overflow], 1, 'FaceColor', [0.4 0.4 0.4], 'FaceAlpha', 0.6, 'EdgeColor', 'k', 'HandleVisibility', 'off');
        text(ax, hiB + binW, overflow, num2str(overflow), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end

    allCounts = [counts, underflow, overflow];
    peak = max(allCounts);
    rest = allCounts(allCounts < peak & allCounts > 0);
    if ~isempty(rest) && peak > 4 * max(rest)
        set(ax, 'YScale', 'log');
        ylim(ax, [0.8, peak * 1.3]);
    end
end

function drawMCTraj_v5(ax, pos_all, target_pos, x_idx, y_idx, x_lbl, y_lbl, alphaVal, fail, failColor)
    hold(ax, 'on'); grid(ax, 'on'); axis(ax, 'equal');
    T = size(pos_all, 2);
    for pass = [false true]
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

function tf = hasOutVar(obj, name)
    if isstruct(obj)
        tf = isfield(obj, name);
    else
        tf = any(strcmp(who(obj), name));
    end
end

function wind_all = extractWindAll(S, t_common, num_sims, ok)
    wind_all = getField(S, 'wind_all', []);
    if ~isempty(wind_all) && size(wind_all, 1) == num_sims
        return;
    end
    wind_all = [];
    if ~isfield(S, 'out') || isempty(S.out), return; end
    buf = nan(num_sims, numel(t_common), 3);
    found = false;
    for i = 1:min(num_sims, numel(S.out))
        if ~ok(i) || ~hasOutVar(S.out(i), 'SSE_Wind'), continue; end
        w_raw = S.out(i).SSE_Wind;
        if isa(w_raw, 'timeseries')
            w_t = w_raw.Time(:);
            w_vec = extractLoggedState(w_raw);
        else
            w_vec = squeeze(w_raw);
            if size(w_vec, 1) == 3 && size(w_vec, 2) ~= 3, w_vec = w_vec.'; end
            w_t = linspace(t_common(1), t_common(end), size(w_vec, 1))';
        end
        if size(w_vec, 2) >= 3 && numel(w_t) >= 2
            buf(i, :, :) = interp1(w_t, w_vec(:, 1:3), t_common, 'linear', NaN);
            found = true;
        end
    end
    if found, wind_all = buf; end
end

function [u_cmd_all, t_u] = extractCmdAll(S, t_common, num_sims, ok)
    u_cmd_all = getField(S, 'u_cmd_all', []);
    if ~isempty(u_cmd_all) && size(u_cmd_all, 1) == num_sims
        t_u = reshape(getField(S, 't_u', linspace(t_common(1), t_common(end), size(u_cmd_all, 2))'), [], 1);
        return;
    end
    u_cmd_all = []; t_u = [];
    if ~isfield(S, 'out') || isempty(S.out), return; end
    for i = 1:min(num_sims, numel(S.out))
        if ok(i) && hasOutVar(S.out(i), 'input_log') && ~isempty(S.out(i).input_log)
            ts = S.out(i).input_log;
            if numel(ts.Time) >= 2
                t_u = ts.Time(:);
                break;
            end
        end
    end
    if isempty(t_u), return; end
    buf = nan(num_sims, numel(t_u), 4);
    found = false;
    for i = 1:min(num_sims, numel(S.out))
        if ~ok(i) || ~hasOutVar(S.out(i), 'input_log'), continue; end
        ts = S.out(i).input_log;
        u_raw = extractLoggedState(ts);
        if size(u_raw, 2) >= 2 && numel(ts.Time) >= 2
            nCh = min(4, size(u_raw, 2));
            buf(i, :, 1:nCh) = interp1(ts.Time(:), u_raw(:, 1:nCh), t_u, 'linear', NaN);
            found = true;
        end
    end
    if found, u_cmd_all = buf; else, t_u = []; end
end

function plotOverlayBands(ax, t, data3d, axColors, axis_lbls)
    % Draw 5-95% shaded bands first, then medians on top using axis colors.
    C = size(data3d, 3);
    P_all = cell(1, C);
    g_all = cell(1, C);
    for c = 1:C
        P = pctlCols(data3d(:, :, c), [5 50 95]);
        good = all(isfinite(P), 1);
        P_all{c} = P(:, good);
        g_all{c} = t(good)';
        if any(good)
            fill(ax, [g_all{c}, fliplr(g_all{c})], [P_all{c}(1,:), fliplr(P_all{c}(3,:))], ...
                axColors{c}, 'FaceAlpha', 0.22, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
    for c = 1:C
        if ~isempty(g_all{c})
            plot(ax, g_all{c}, P_all{c}(2,:), '-', 'Color', axColors{c}, ...
                'LineWidth', 1.6, 'DisplayName', axis_lbls{c});
        end
    end
end

function plotPSDGroup(ax, t_sig, data3d, evalMask, chIdx, chNames, chColors, failMask, failColor, unitLbl, ttl)
    % Hann-windowed one-sided PSD (median + 5-95% band) in pure base MATLAB.
    hold(ax, 'on'); grid(ax, 'on');
    dt_sig = median(diff(t_sig));
    M = nnz(evalMask);
    if ~(dt_sig > 0) || M < 8
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'Insufficient samples for PSD', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(ax, ttl); return;
    end

    fs = 1 / dt_sig;
    win = 0.5 * (1 - cos(2 * pi * (0:M-1)' / max(1, M - 1)));
    nfft = max(2^nextpow2(M), 1024);
    f = (0:nfft/2)' * (fs / nfft);
    fKeep = f >= max(0.05, f(2));
    fp = f(fKeep)';

    P_med_all = zeros(numel(chIdx), numel(fp));
    Pf_overlay = [];
    anyPlotted = false;
    for k = 1:numel(chIdx)
        X = data3d(:, evalMask, chIdx(k));
        good = all(isfinite(X), 2);
        if ~any(good), continue; end
        Xg = X(good, :) - mean(X(good, :), 2);
        Xf = fft(Xg .* win', nfft, 2);
        Pxx = (abs(Xf(:, 1:nfft/2+1)).^2) / (fs * sum(win.^2));
        Pxx(:, 2:end-1) = 2 * Pxx(:, 2:end-1);
        Ps = max(pctlCols(Pxx(:, fKeep), [5 50 95]), 1e-12);
        fill(ax, [fp, fliplr(fp)], [Ps(1,:), fliplr(Ps(3,:))], chColors{k}, ...
            'FaceAlpha', 0.16, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        plot(ax, fp, Ps(2,:), '-', 'Color', chColors{k}, 'LineWidth', 1.6, 'DisplayName', chNames{k});
        P_med_all(k, :) = Ps(2, :);
        anyPlotted = true;

        if k == 1 && any(failMask(:) & good) && any(~failMask(:) & good)
            Pf_overlay = max(pctlCols(Pxx(failMask(:) & good, fKeep), 50), 1e-12);
        end
    end

    if ~anyPlotted
        axis(ax, 'off');
        text(ax, 0.5, 0.5, 'No finite data for PSD', 'HorizontalAlignment', 'center', 'Units', 'normalized');
        title(ax, ttl); return;
    end
    if ~isempty(Pf_overlay)
        plot(ax, fp, Pf_overlay(1,:), '--', 'Color', failColor, 'LineWidth', 1.3, ...
            'DisplayName', sprintf('%s (failed)', chNames{1}));
    end

    % Highlight top 2 distinct oscillatory peaks above 0.25 Hz (skipping DC/circle fundamental)
    env = max(P_med_all, [], 1);
    oscMask = fp >= 0.25;
    if nnz(oscMask) >= 3
        fo = fp(oscMask); eo = env(oscMask);
        pkIdx = find(eo(2:end-1) > eo(1:end-2) & eo(2:end-1) >= eo(3:end)) + 1;
        if ~isempty(pkIdx)
            [~, ord] = sort(eo(pkIdx), 'descend');
            sel = [];
            for m = ord
                if isempty(sel) || all(abs(log10(fo(pkIdx(m)) ./ fo(sel))) >= 0.15)
                    sel(end+1) = pkIdx(m); %#ok<AGROW>
                end
                if numel(sel) >= 2, break; end
            end
            for j = 1:numel(sel)
                fx = fo(sel(j)); fy = eo(sel(j));
                plot(ax, fx, fy, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'w', 'HandleVisibility', 'off');
                text(ax, fx, fy * (1.45 + 0.45 * mod(j - 1, 2)), sprintf('%.2f Hz', fx), ...
                    'HorizontalAlignment', 'center', 'FontSize', 7.5);
            end
        end
    end

    set(ax, 'XScale', 'log', 'YScale', 'log');
    xlim(ax, [fp(1), fs / 2]);
    xlabel(ax, 'Frequency (Hz)'); ylabel(ax, sprintf('PSD ((%s)^2/Hz)', unitLbl));
    title(ax, ttl);
    legend(ax, 'Location', 'southwest');
end
