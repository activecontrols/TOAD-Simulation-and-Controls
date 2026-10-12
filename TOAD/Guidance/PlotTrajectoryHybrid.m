function figHandles = PlotTrajectoryHybrid(opt, varargin)
%% PLOTTRAJECTORYHYBRID  Publication-quality 3D trajectory and control dashboard.
% Standalone visualization module for TrajectoryOptimizerHybrid.
% Supports both full 6-DoF solution dashboards and Stage-1 analytical initial guess figures.
%
% Usage:
%   fig = PlotTrajectoryHybrid(opt);                  % Full 4-subplot dashboard
%   fig = PlotTrajectoryHybrid(opt, 'InitialGuess', true); % Stage-1 initial guess
%
% Authors: PSP Active Controls & Antigravity

    p = inputParser;
    addParameter(p, 'InitialGuess', false, @islogical);
    parse(p, varargin{:});

    if p.Results.InitialGuess
        figHandles = plotInitialGuessDashboard(opt);
    else
        figHandles = plotSolutionDashboard(opt);
    end
end

function f1 = plotSolutionDashboard(opt)
    if isempty(opt.Solution) || ~isfield(opt.Solution, 'X')
        error('PlotTrajectoryHybrid:NoSolution', 'No solution available to plot. Call solve() first.');
    end

    tv = opt.Solution.Time;
    Xs = opt.Solution.X;
    Us = opt.Solution.U;

    f1 = figure('Name', sprintf('%s 6-DoF Trajectory (%s v%03d)', ...
                opt.Maneuver, opt.getVehicleName(), opt.Version), ...
                'Position', [100, 100, 1300, 800], 'Visible', 'on', 'Color', 'w');

    % -----------------------------------------------------------------
    % 1. 3D Flight Profile with Glideslope Cone & Landing Funnel
    % -----------------------------------------------------------------
    ax3d = subplot(2, 2, [1, 3]);
    set(ax3d, 'Tag', 'FlightProfile3D');
    hold(ax3d, 'on'); grid(ax3d, 'on'); box(ax3d, 'on');

    % Ground shadow projection on z = 0
    plot3(ax3d, Xs(5,:), Xs(6,:), zeros(size(Xs(7,:))), ':', ...
          'Color', [0.72 0.72 0.72], 'LineWidth', 1.2);

    % 3D trajectory colored continuously by mission time
    h_traj = surface(ax3d, [Xs(5,:); Xs(5,:)], [Xs(6,:); Xs(6,:)], [Xs(7,:); Xs(7,:)], ...
                     [tv; tv], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.8);

    % Launch pad & Landing target
    p_launch = plot3(ax3d, opt.r0(1), opt.r0(2), opt.r0(3), 'o', 'MarkerSize', 9, ...
                     'MarkerFaceColor', [0.15 0.75 0.20], 'MarkerEdgeColor', [0.05 0.40 0.10], 'LineWidth', 1.8);
    p_land   = plot3(ax3d, opt.r_f(1), opt.r_f(2), opt.r_f(3), 's', 'MarkerSize', 9, ...
                     'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', [0.50 0.05 0.05], 'LineWidth', 1.8);

    leg_handles = [h_traj, p_launch, p_land];
    leg_labels  = {'6-DoF Trajectory', 'Launch Pad', 'Landing Target'};

    % Takeoff Glideslope Cone wireframe (GlideslopeAngle half-angle)
    z_max_traj = max(Xs(7,:));
    z_to_max = opt.getTakeoffClearanceAlt();
    z_cone = linspace(0, z_to_max, 22);
    th_c = linspace(0, 2*pi, 36);
    [TH_to, ZC_to] = meshgrid(th_c, z_cone);
    R_to = ZC_to * tand(opt.GlideslopeAngle) + 0.15;
    XC_to = opt.r0(1) + R_to .* cos(TH_to);
    YC_to = opt.r0(2) + R_to .* sin(TH_to);
    p_to = mesh(ax3d, XC_to, YC_to, ZC_to, 'FaceColor', [0.15 0.75 0.20], 'FaceAlpha', 0.05, ...
                'EdgeColor', [0.15 0.70 0.20], 'EdgeAlpha', 0.35, 'LineStyle', ':');
    leg_handles(end+1) = p_to;
    leg_labels{end+1}  = sprintf('Takeoff Cone (%.0f°)', opt.GlideslopeAngle);

    % Flaring Landing Funnel wireframe (quadratic flare)
    z_ld_max = opt.getLandingFlareAlt();
    cf_plot = opt.FunnelCurvature;
    z_fun = linspace(0, z_ld_max, 25);
    [TH_ld, ZC_ld] = meshgrid(th_c, z_fun);
    gsa = opt.GlideslopeAngle;
    R_ld = ZC_ld * tand(gsa) + cf_plot * (ZC_ld.^2) + 0.15;
    XC_ld = opt.r_f(1) + R_ld .* cos(TH_ld);
    YC_ld = opt.r_f(2) + R_ld .* sin(TH_ld);
    p_ld = mesh(ax3d, XC_ld, YC_ld, ZC_ld, 'FaceColor', [0.85 0.15 0.15], 'FaceAlpha', 0.05, ...
                'EdgeColor', [0.85 0.15 0.15], 'EdgeAlpha', 0.35, 'LineStyle', ':');
    leg_handles(end+1) = p_ld;
    leg_labels{end+1}  = sprintf('Landing Funnel (%.0f°)', gsa);

    % Maneuver overlays: Circle orbit or Waypoints
    if opt.Maneuver == "Circle"
        thc = linspace(0, 2*pi, 100);
        mp = opt.ManeuverParams;
        p_orb = plot3(ax3d, mp.circle_center(1) + mp.circle_radius * cos(thc), ...
                      mp.circle_center(2) + mp.circle_radius * sin(thc), ...
                      mp.circle_alt * ones(size(thc)), 'k--', 'LineWidth', 1.5);
        leg_handles(end+1) = p_orb;
        leg_labels{end+1}  = 'Nominal Orbit';
    elseif (opt.Maneuver == "Waypoint" || opt.Maneuver == "Custom") && ~isempty(opt.Waypoints)
        p_wp = plot3(ax3d, opt.Waypoints(1,:), opt.Waypoints(2,:), opt.Waypoints(3,:), 'd--', ...
                     'Color', [0.90 0.50 0.10], 'LineWidth', 1.6, 'MarkerSize', 7, ...
                     'MarkerFaceColor', [1.0 0.75 0.20], 'MarkerEdgeColor', [0.65 0.30 0.05]);
        for kw = 1:size(opt.Waypoints, 2)
            text(ax3d, opt.Waypoints(1,kw), opt.Waypoints(2,kw), opt.Waypoints(3,kw) + 0.04*z_max_traj, ...
                 sprintf(' WP%d', kw), 'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.65 0.30 0.05]);
        end
        leg_handles(end+1) = p_wp;
        leg_labels{end+1}  = 'Survey Waypoints';
    end

    % Apex Marker
    [z_peak, k_peak] = max(Xs(7,:));
    plot3(ax3d, Xs(5,k_peak), Xs(6,k_peak), z_peak, 'k^', 'MarkerSize', 6, ...
          'MarkerFaceColor', [0.1 0.1 0.1]);
    text(ax3d, Xs(5,k_peak), Xs(6,k_peak), z_peak + 0.03*z_max_traj, ...
         sprintf('Apex: %.1fm', z_peak), 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

    % Time colorbar
    colormap(ax3d, jet(256));
    cb = colorbar(ax3d, 'Location', 'southoutside');
    cb.Label.String = 'Mission Time [s]';
    cb.FontSize = 8.5;

    legend(ax3d, leg_handles, leg_labels, 'Location', 'best', 'FontSize', 8);
    axis(ax3d, 'equal');

    % Square ground footprint bounds
    x_pts = [Xs(5,:), opt.r0(1), opt.r_f(1), XC_to(:)', XC_ld(:)'];
    y_pts = [Xs(6,:), opt.r0(2), opt.r_f(2), YC_to(:)', YC_ld(:)'];
    if ~isempty(opt.Waypoints)
        x_pts = [x_pts, opt.Waypoints(1,:)];
        y_pts = [y_pts, opt.Waypoints(2,:)];
    end
    if opt.Maneuver == "Circle" && isfield(opt.ManeuverParams, 'circle_center')
        mp = opt.ManeuverParams;
        x_pts = [x_pts, mp.circle_center(1) - mp.circle_radius, mp.circle_center(1) + mp.circle_radius];
        y_pts = [y_pts, mp.circle_center(2) - mp.circle_radius, mp.circle_center(2) + mp.circle_radius];
    end

    x_min = min(x_pts); x_max = max(x_pts);
    y_min = min(y_pts); y_max = max(y_pts);
    base_span = max([x_max - x_min, y_max - y_min, 4.0]);
    half_w = 0.58 * base_span;
    x_mid = 0.5 * (x_min + x_max);
    y_mid = 0.5 * (y_min + y_max);
    xlim(ax3d, [x_mid - half_w, x_mid + half_w]);
    ylim(ax3d, [y_mid - half_w, y_mid + half_w]);

    z_top = max([z_max_traj * 1.08, z_ld_max * 1.15, 5.0]);
    zlim(ax3d, [0, z_top]);
    xlabel(ax3d, 'X (North) [m]', 'FontWeight', 'bold');
    ylabel(ax3d, 'Y (East) [m]', 'FontWeight', 'bold');
    zlabel(ax3d, 'Z (Altitude) [m]', 'FontWeight', 'bold');
    title(ax3d, sprintf('%s Flight Profile (%s, Duration: %.2f s, Path: %.1f m)', ...
        opt.Maneuver, opt.getVehicleName(), opt.Solution.T_total, opt.Solution.L_path), ...
        'FontWeight', 'bold');

    switch opt.Maneuver
        case "Circle",    view(ax3d, -35, 26);
        case "Hop",       view(ax3d, -25, 22);
        case "Backflip",  view(ax3d, -30, 20);
        case "Waypoint",  view(ax3d, -32, 28);
        otherwise,        view(ax3d, -35, 25);
    end

    % -----------------------------------------------------------------
    % 2. TVC Gimbal Deflections
    % -----------------------------------------------------------------
    subplot(2, 2, 2);
    plot(tv(1:end-1), rad2deg(Us(1,:)), 'r-', 'LineWidth', 1.6); hold on;
    plot(tv(1:end-1), rad2deg(Us(2,:)), 'b-', 'LineWidth', 1.6);
    yline(rad2deg(opt.max_gimbal_angle), 'k:', 'Max Gimbal', 'LabelHorizontalAlignment', 'left');
    yline(-rad2deg(opt.max_gimbal_angle), 'k:', '-Max Gimbal', 'LabelHorizontalAlignment', 'left');
    grid on; xlabel('Time [s]', 'FontWeight', 'bold');
    ylabel('Gimbal Angle [deg]', 'FontWeight', 'bold');
    legend({'\theta (Pitch)', '\phi (Yaw)'}, 'Location', 'best');
    title('TVC Gimbal Deflections', 'FontWeight', 'bold');

    % -----------------------------------------------------------------
    % 3. Commanded Thrust vs Time
    % -----------------------------------------------------------------
    subplot(2, 2, 4);
    plot(tv(1:end-1), Us(3,:), 'k-', 'LineWidth', 1.8); hold on;
    hover_thrust = opt.constants.m_dry * opt.constants.g;
    yline(hover_thrust, 'r--', sprintf('Dry Hover (%.1f N)', hover_thrust), 'LabelHorizontalAlignment', 'left');
    yline(opt.constants.MaxThrust, 'b:', sprintf('Max Thrust (%.1f N)', opt.constants.MaxThrust), 'LabelHorizontalAlignment', 'left');
    grid on; xlabel('Time [s]', 'FontWeight', 'bold');
    ylabel('Thrust [N]', 'FontWeight', 'bold');
    title(sprintf('Commanded Thrust Profile (Status: %s)', opt.Solution.Status), 'FontWeight', 'bold');
    legend({'Thrust Cmd', 'Hover Ref', 'Max Limit'}, 'Location', 'best');
end

function figHandle = plotInitialGuessDashboard(opt)
    if isempty(opt.InitialGuess) || ~isfield(opt.InitialGuess, 'X')
        opt.buildInitialGuess();
    end
    tv = opt.InitialGuess.Time;
    Xg = opt.InitialGuess.X;
    Ug = opt.InitialGuess.U;

    figHandle = figure('Name', sprintf('%s Initial Guess (%s)', opt.Maneuver, opt.getVehicleName()), ...
                       'Position', [150, 150, 950, 650], 'Visible', 'on', 'Color', 'w');
    ax_g3d = subplot(2, 1, 1);
    set(ax_g3d, 'Tag', 'InitialGuess3D');
    hold(ax_g3d, 'on'); grid(ax_g3d, 'on'); box(ax_g3d, 'on');
    plot3(ax_g3d, Xg(5,:), Xg(6,:), Xg(7,:), 'm--', 'LineWidth', 2);
    plot3(ax_g3d, opt.r0(1), opt.r0(2), opt.r0(3), 'go', 'MarkerSize', 8, 'LineWidth', 2);
    plot3(ax_g3d, opt.r_f(1), opt.r_f(2), opt.r_f(3), 'rs', 'MarkerSize', 8, 'LineWidth', 2);

    axis(ax_g3d, 'equal');
    x_pts_g = [Xg(5,:), opt.r0(1), opt.r_f(1)];
    y_pts_g = [Xg(6,:), opt.r0(2), opt.r_f(2)];
    base_g = max([max(x_pts_g) - min(x_pts_g), max(y_pts_g) - min(y_pts_g), 4.0]);
    half_g = 0.58 * base_g;
    x_mid_g = 0.5 * (min(x_pts_g) + max(x_pts_g));
    y_mid_g = 0.5 * (min(y_pts_g) + max(y_pts_g));
    xlim(ax_g3d, [x_mid_g - half_g, x_mid_g + half_g]);
    ylim(ax_g3d, [y_mid_g - half_g, y_mid_g + half_g]);
    zlim(ax_g3d, [0, max(max(Xg(7,:)) * 1.08, 5.0)]);
    xlabel('X [m]', 'FontWeight', 'bold');
    ylabel('Y [m]', 'FontWeight', 'bold');
    zlabel('Z [m]', 'FontWeight', 'bold');
    title(sprintf('%s Stage-1 QP Seed 3D Trajectory (%s)', opt.Maneuver, opt.getVehicleName()), 'FontWeight', 'bold');
    legend({'QP Seed', 'Launch', 'Landing'}, 'Location', 'best');

    switch opt.Maneuver
        case "Circle",    view(-35, 26);
        case "Hop",       view(-25, 22);
        case "Backflip",  view(-30, 20);
        case "Waypoint",  view(-32, 28);
        otherwise,        view(-35, 25);
    end

    subplot(2, 1, 2);
    plot(tv(1:end-1), Ug(3,:), 'k-', 'LineWidth', 1.5);
    grid on; xlabel('Time [s]', 'FontWeight', 'bold');
    ylabel('Thrust [N]', 'FontWeight', 'bold');
    title('Stage-1 Seed Commanded Thrust', 'FontWeight', 'bold');
end
