function abortHarness()
% abortHarness  Monte Carlo test of ControllerAndEstimator::compute_abort (via compute_abort_mex).
%
% Build once (adjust paths; needs Eigen, matlab_funcs.*, estimator/controller sources that
% controller_and_estimator.cpp links against, and the Simulink headers for simstruc.h):
%   mex compute_abort_mex.cpp controller_and_estimator.cpp computeAbort.cpp flatnessInvert.cpp ...
%       allocateControls.cpp matlab_funcs.cpp -I<eigen> -I(fullfile(matlabroot,'simulink','include'))
%
% Frame: NWU (X north, Y west, Z up), pad at origin. Quaternion scalar-first [w x y z].

%% ---- Pass/fail limits (evaluated on the generated trajectory, CG states) ----
MAX_VEL_LIMIT  = 30;       % [m/s]   fail if max ||v||  over trajectory exceeds this
MAX_RATE_LIMIT = 5*pi;     % [rad/s] fail if max ||w||  over trajectory exceeds this

%% ---- Dispersions (initial abort state) ----
N         = 500;           % runs
SEED      = 1;             % rng seed
TILT_MAX  = 30;            % [deg]   body-z tilt from vertical
ALT_MAX   = 15;            % [m]     altitude drawn U(0, ALT_MAX)
XY_MAX    = 5;             % [m]     x,y drawn U(-XY_MAX, XY_MAX) about the pad
V0_MAX    = 5;             % [m/s]   ||v0|| drawn U(0, V0_MAX), random direction
W0_MAX    = pi/4;          % [rad/s] ||w0|| drawn U(0, W0_MAX), random direction
LAST_THRUST = 1.169*9.8015; % [N]     hover thrust (m*g from constantsASTRA); computeAbort clamps it
T_NOW     = 0;             % [s]     abort start time
MAXPTS    = 100;           % = ABORT_MAX_POINTS in AbortCode.h (row padding for NaN-filled logs)

if exist('compute_abort_mex', 'file') ~= 3
    error('abortHarness:NoMex', 'compute_abort_mex not built. See build line in header.');
end

%% ---- Sample initial states ----
rng(SEED);
tilt = deg2rad(TILT_MAX) * rand(N,1);              % tilt magnitude (uniform in angle)
az   = 2*pi * rand(N,1);                           % direction of tilt
yaw  = 2*pi * rand(N,1);                           % heading about vertical
% q = q_yaw(z) (x) q_tilt(axis in xy-plane). C33 = 1-2(qx^2+qy^2) = cos(tilt) by construction.
qt = [cos(tilt/2), sin(tilt/2).*cos(az), sin(tilt/2).*sin(az), zeros(N,1)];
a = cos(yaw/2); d = sin(yaw/2);
q0   = [a.*qt(:,1), a.*qt(:,2)-d.*qt(:,3), a.*qt(:,3)+d.*qt(:,2), d.*qt(:,1)];
pos0 = [XY_MAX*(2*rand(N,2)-1), ALT_MAX*rand(N,1)];
dirV = randn(N,3); v0 = V0_MAX * rand(N,1) .* (dirV ./ vecnorm(dirV,2,2));
dirW = randn(N,3); w0 = W0_MAX * rand(N,1) .* (dirW ./ vecnorm(dirW,2,2));

%% ---- Run aborts ----
tr_t = nan(N,MAXPTS); tr_q = nan(N,MAXPTS,4); tr_pos = nan(N,MAXPTS,3);
tr_vel = nan(N,MAXPTS,3); tr_w = nan(N,MAXPTS,3);
bad = false(N,1);                                  % mex error or non-finite output
for i = 1:N
    X0 = [q0(i,:), pos0(i,:), v0(i,:), w0(i,:), 0, 0];   % layout of last_X_ctrl
    try
        [t_i, X_i, ~] = compute_abort_mex(T_NOW, X0, LAST_THRUST);
        if any(~isfinite(X_i(:))), bad(i) = true; end
        n = numel(t_i);
        tr_t(i,1:n)     = t_i;
        tr_q(i,1:n,:)   = reshape(X_i(:,1:4),   1, n, 4);
        tr_pos(i,1:n,:) = reshape(X_i(:,5:7),   1, n, 3);
        tr_vel(i,1:n,:) = reshape(X_i(:,8:10),  1, n, 3);
        tr_w(i,1:n,:)   = reshape(X_i(:,11:13), 1, n, 3);
    catch ME
        bad(i) = true; fprintf('run %d: %s\n', i, ME.message);
    end
end

%% ---- Metrics and pass/fail ----
vel_norm  = vecnorm(tr_vel, 2, 3);                 % N x MAXPTS [m/s]
rate_norm = vecnorm(tr_w,   2, 3);                 % N x MAXPTS [rad/s]
pos_norm  = vecnorm(tr_pos, 2, 3);
alt_sign  = tr_pos(:, :, 3);
tilt_deg  = acosd(max(-1, min(1, 1 - 2*(tr_q(:,:,2).^2 + tr_q(:,:,3).^2))));  % body-z vs vertical
maxVel    = max(vel_norm,  [], 2);
maxRate   = max(rate_norm, [], 2);
maxPos   = max(pos_norm,  [], 2);
minAlt   = min(alt_sign, [], 2);
failVel   = maxVel  > MAX_VEL_LIMIT;
failRate  = maxRate > MAX_RATE_LIMIT;
failPos   = maxPos  > 30;
failAlt   = minAlt < 0;

fail      = bad | failVel | failRate | failPos | failAlt;
fprintf('\n=== Abort MC: %d runs, seed %d, %d failed (%.1f%%) ===\n', N, SEED, nnz(fail), 100*nnz(fail)/N);
fprintf('   mex error/non-finite %d | max|v| > %.1f m/s: %d | max|w| > %.2f rad/s: %d\n', ...
    nnz(bad), MAX_VEL_LIMIT, nnz(failVel), MAX_RATE_LIMIT, nnz(failRate));

%% ---- Plots (colors/layout follow PlotMC_v5) ----
passColor = [0.15 0.68 0.20]; failColor = [0.85 0.33 0.10];
lineClr = {[0.25 0.65 0.25 0.82], [failColor 0.36]};   % {passed, failed}; failed drawn on top

% Figure 1: summary
figure('Name','Abort Summary','Color','w','WindowStyle','docked');
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); axis(ax,'off');
text(ax,0,1,{sprintf('Runs:         %d  (seed %d)',N,SEED), ...
    sprintf('Failed:       %d (%.1f%%)',nnz(fail),100*nnz(fail)/N), ...
    sprintf('Mex/NaN:      %d',nnz(bad)), ...
    sprintf('Vel limit:    %.1f m/s  (fails %d)',MAX_VEL_LIMIT,nnz(failVel)), ...
    sprintf('Rate limit:   %.2f rad/s (fails %d)',MAX_RATE_LIMIT,nnz(failRate))}, ...
    'VerticalAlignment','top','FontName','FixedWidth','FontSize',8.5);
title(ax,'Campaign');

pairs = {maxVel,        maxRate,        'Max ||v|| [m/s]',      'Max ||w|| [rad/s]',  'Failure Space';
         pos0(:,3),     vecnorm(v0,2,2),'Initial altitude [m]', 'Initial ||v|| [m/s]', 'Initial Altitude vs Speed';
         tilt*180/pi,   vecnorm(w0,2,2),'Initial tilt [deg]',   'Initial ||w|| [rad/s]','Initial Tilt vs Rate'};
for k = 1:3
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    scatter(ax,pairs{k,1}(~fail),pairs{k,2}(~fail),28,passColor,'filled','MarkerFaceAlpha',0.55,'DisplayName','Pass');
    scatter(ax,pairs{k,1}(fail), pairs{k,2}(fail), 28,failColor,'filled','MarkerFaceAlpha',0.85,'DisplayName','Fail');
    xlabel(ax,pairs{k,3}); ylabel(ax,pairs{k,4}); title(ax,pairs{k,5});
    if k == 1
        xline(ax,MAX_VEL_LIMIT,'r--','LineWidth',1.2,'DisplayName','Limit');
        yline(ax,MAX_RATE_LIMIT,'r--','LineWidth',1.2,'HandleVisibility','off');
        legend(ax,'Location','best');
    end
end

% Figure 2: 3D trajectories + projections
figure('Name','Abort 3D Trajectories','Color','w','WindowStyle','docked');
tl = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
dims = {[1 2 3],[1 2],[1 3],[2 3]};  tiles = {[],4,8,12};  % main (3x3 span), top, side, front
lbl = {'X (North) [m]','Y (West) [m]','Z (Up) [m]'};
ttl = {sprintf('Abort Trajectories (%d runs, %d failed)',N,nnz(fail)),'Top View','Side View','Front View'};
for a = 1:4
    if a == 1, ax = nexttile(tl,1,[3 3]); else, ax = nexttile(tl,tiles{a}); end
    hold(ax,'on'); grid(ax,'on'); axis(ax,'equal');
    for isFail = [false true]
        idx = find(fail == isFail)';
        if isempty(idx), continue; end
        c = lineClr{isFail+1};
        if a == 1
            plot3(ax,tr_pos(idx,:,1)',tr_pos(idx,:,2)',tr_pos(idx,:,3)','Color',c,'LineWidth',0.5,'HandleVisibility','off');
        else
            plot(ax,tr_pos(idx,:,dims{a}(1))',tr_pos(idx,:,dims{a}(2))','Color',c,'LineWidth',0.5,'HandleVisibility','off');
        end
    end
    if a == 1
        plot3(ax,0,0,0,'b^','MarkerFaceColor','b','MarkerSize',8,'DisplayName','Pad'); view(ax,3);
        zlabel(ax,lbl{3});
        plot3(ax,nan,nan,nan,'-','Color',[0.25 0.25 0.25],'DisplayName','Passed');
        plot3(ax,nan,nan,nan,'-','Color',failColor,'DisplayName','Failed');
        legend(ax,'Location','best');
    else
        plot(ax,0,0,'b^','MarkerFaceColor','b','MarkerSize',6);
    end
    xlabel(ax,lbl{dims{a}(1)}); ylabel(ax,lbl{dims{a}(2)}); title(ax,ttl{a});
end

% Figure 3: time histories vs pass/fail limits
figure('Name','Abort Time Histories','Color','w','WindowStyle','docked');
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
Y     = {vel_norm, rate_norm, tr_pos(:,:,3), tilt_deg};
yLbl  = {'||v|| [m/s]','||w|| [rad/s]','Altitude [m]','Tilt [deg]'};
yLim  = [MAX_VEL_LIMIT, MAX_RATE_LIMIT, NaN, NaN];
for k = 1:4
    ax = nexttile(tl); hold(ax,'on'); grid(ax,'on');
    for isFail = [false true]
        idx = find(fail == isFail)';
        if ~isempty(idx), plot(ax,tr_t(idx,:)'-T_NOW,Y{k}(idx,:)','Color',lineClr{isFail+1},'LineWidth',0.5); end
    end
    if isfinite(yLim(k)), yline(ax,yLim(k),'r--','LineWidth',1.5); end
    xlabel(ax,'Time since abort [s]'); ylabel(ax,yLbl{k}); title(ax,yLbl{k});
end
end
