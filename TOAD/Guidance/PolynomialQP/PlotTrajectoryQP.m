function figHandles = PlotTrajectoryQP(InitialGuess, QPTraj, Report, Waypoints, Constants, varargin)
%% 4-Figure Docked Engineering Suite for PolyQP v2 Trajectory Generation
    p = inputParser;
    addParameter(p, 'WindowStyle', 'docked');
    addParameter(p, 'Visible', 'on');
    parse(p, varargin{:});
    winStyle = p.Results.WindowStyle;
    figVis   = p.Results.Visible;

    Limits = buildVehicleLimits(Constants);
    isTOAD = isfield(Constants, 'Vehicle') && Constants.Vehicle == 1;
    if isTOAD
        vehName = 'TOAD (Vehicle 1)';
    else
        vehName = 'ASTRA (Vehicle 0)';
    end

    % Unpack sampled trajectory (DtSample = 0.05 s)
    t   = InitialGuess.Time(:)';
    X   = InitialGuess.X;
    U   = InitialGuess.U;
    Nn  = numel(t);
    Nc  = size(U, 2);
    tU  = t(1:Nc);
    dt  = InitialGuess.dt_row(1);

    q_cg = X(1:4, :);
    r_cg = X(5:7, :);
    v_cg = X(8:10, :);
    w_B  = X(11:13, :);
    mLox = X(14, :);
    mIpa = X(15, :);
    mTot = Constants.m_dry + mLox + mIpa;

    % Evaluate exact polynomial P-space derivatives across sample grid
    M  = QPTraj.Count;
    tK = QPTraj.KeyTimes(1:(M + 1));
    P_pos  = zeros(3, Nn);
    P_vel  = zeros(3, Nn);
    P_acc  = zeros(3, Nn);
    P_jerk = zeros(3, Nn);
    P_snap = zeros(3, Nn);
    CGz_t  = zeros(1, Nn);
    tiltDeg = zeros(1, Nn);
    hdgDeg  = zeros(1, Nn);
    b3_hist = zeros(3, Nn);

    for k = 1:Nn
        P_pos(:, k)  = PolyEval(QPTraj, t(k), 0);
        P_vel(:, k)  = PolyEval(QPTraj, t(k), 1);
        P_acc(:, k)  = PolyEval(QPTraj, t(k), 2);
        P_jerk(:, k) = PolyEval(QPTraj, t(k), 3);
        P_snap(:, k) = PolyEval(QPTraj, t(k), 4);

        C_BI = quatRot(q_cg(:, k));
        C_IB = C_BI.';
        b1 = C_IB(:, 1);
        b3 = C_IB(:, 3);
        b3_hist(:, k) = b3;
        tiltDeg(k) = rad2deg(acos(max(-1.0, min(1.0, b3(3)))));
        hdgDeg(k)  = rad2deg(atan2(b1(2), b1(1)));

        if isTOAD
            [~, CGz_k] = ComputeJtot(mLox(k), mIpa(k), Constants);
            CGz_t(k) = CGz_k;
        else
            CGz_t(k) = Constants.rTB;
        end
    end

    % Compute CG acceleration via central differences for comparison with P_acc
    a_cg = zeros(3, Nn);
    if Nn >= 3
        a_cg(:, 2:(Nn - 1)) = (v_cg(:, 3:Nn) - v_cg(:, 1:(Nn - 2))) / (2 * dt);
        a_cg(:, 1)  = (v_cg(:, 2) - v_cg(:, 1)) / dt;
        a_cg(:, Nn) = (v_cg(:, Nn) - v_cg(:, Nn - 1)) / dt;
    end

    % Palette
    cX = [0.85, 0.22, 0.18];
    cY = [0.12, 0.58, 0.28];
    cZ = [0.12, 0.38, 0.82];
    cNorm = [0.15, 0.15, 0.18];
    cLim  = [0.80, 0.15, 0.15];

    figHandles = gobjects(1, 4);

    %% =====================================================================
    %% FIGURE 1: 3D Mission Profile & Projections
    %% =====================================================================
    figHandles(1) = createDockedFig('1: 3D Mission Profile & Projections', winStyle, figVis);
    tl1 = tiledlayout(figHandles(1), 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl1, sprintf('PolyQP v2 Mission Profile — %s  (T = %.2f s, %d Segments, Feasible = %d)', ...
        vehName, QPTraj.T_total, M, Report.Feasible), 'FontWeight', 'bold', 'FontSize', 12);

    % 3D Isometric View (spans 3x2 on left)
    ax3d = nexttile(tl1, 1, [3, 2]);
    hold(ax3d, 'on'); grid(ax3d, 'on'); box(ax3d, 'on');
    zMin = min(0, min(r_cg(3, :)));

    % Ground shadow
    plot3(ax3d, r_cg(1, :), r_cg(2, :), zMin * ones(1, Nn), '--', ...
        'Color', [0.70, 0.70, 0.74], 'LineWidth', 1.2, 'DisplayName', 'Ground Track');

    % CoP trajectory P(t)
    plot3(ax3d, P_pos(1, :), P_pos(2, :), P_pos(3, :), ':', ...
        'Color', [0.55, 0.35, 0.75], 'LineWidth', 1.3, 'DisplayName', 'CoP Flat Output P(t)');

    % Color-coded CG segments by mode (Hold / Cruise / Flight)
    plotSegmented3D(ax3d, t, r_cg, QPTraj);

    % Attitude b3 thrust vectors along trajectory
    nTriads = min(22, Nn);
    idxTriad = unique(round(linspace(1, Nn, nTriads)));
    span3D = max([max(r_cg(1, :)) - min(r_cg(1, :)), ...
                  max(r_cg(2, :)) - min(r_cg(2, :)), ...
                  max(r_cg(3, :)) - min(r_cg(3, :)), 5.0]);
    vecScale = 0.065 * span3D;
    quiver3(ax3d, r_cg(1, idxTriad), r_cg(2, idxTriad), r_cg(3, idxTriad), ...
        b3_hist(1, idxTriad) * vecScale, b3_hist(2, idxTriad) * vecScale, b3_hist(3, idxTriad) * vecScale, ...
        0, 'Color', [0.92, 0.45, 0.10], 'LineWidth', 1.2, 'MaxHeadSize', 0.4, 'DisplayName', 'Body b_3 Axis');

    % Plot Waypoints & Labels
    plotWaypoints3D(ax3d, Waypoints, span3D);

    xlabel(ax3d, 'North X [m]', 'FontWeight', 'bold');
    ylabel(ax3d, 'West Y [m]', 'FontWeight', 'bold');
    zlabel(ax3d, 'Up Z [m]', 'FontWeight', 'bold');
    title(ax3d, '3D NWU Flight Path & Attitude Vectors');
    axis(ax3d, 'equal');
    zlim(ax3d, [zMin - 0.05 * span3D, max(r_cg(3, :)) + 0.15 * span3D]);
    view(ax3d, -38, 24);
    legend(ax3d, 'Location', 'northwest', 'FontSize', 8);

    % Right stack: 3 Orthographic Projections (XY, XZ, YZ)
    axXY = nexttile(tl1, 3);
    plotOrtho2D(axXY, r_cg(1, :), r_cg(2, :), P_pos(1, :), P_pos(2, :), Waypoints, 1, 2, ...
        'North X [m]', 'West Y [m]', 'Top-Down (X-Y)');

    axXZ = nexttile(tl1, 6);
    plotOrtho2D(axXZ, r_cg(1, :), r_cg(3, :), P_pos(1, :), P_pos(3, :), Waypoints, 1, 3, ...
        'North X [m]', 'Up Z [m]', 'North-Up Profile (X-Z)');

    axYZ = nexttile(tl1, 9);
    plotOrtho2D(axYZ, r_cg(2, :), r_cg(3, :), P_pos(2, :), P_pos(3, :), Waypoints, 2, 3, ...
        'West Y [m]', 'Up Z [m]', 'West-Up Profile (Y-Z)');

    %% =====================================================================
    %% FIGURE 2: Translational Kinematics & High-Order Derivatives
    %% =====================================================================
    figHandles(2) = createDockedFig('2: Translational Kinematics & Derivatives', winStyle, figVis);
    tl2 = tiledlayout(figHandles(2), 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl2, sprintf('Translational Kinematics & Flat-Output Derivatives — %s', vehName), ...
        'FontWeight', 'bold', 'FontSize', 12);

    % Tile 1: Position r_cg vs P
    axPos = nexttile(tl2, 1);
    hold(axPos, 'on'); grid(axPos, 'on'); box(axPos, 'on');
    shadeSegments(axPos, QPTraj);
    plot(axPos, t, r_cg(1, :), '-', 'Color', cX, 'LineWidth', 1.6, 'DisplayName', 'r_{cg,X}');
    plot(axPos, t, r_cg(2, :), '-', 'Color', cY, 'LineWidth', 1.6, 'DisplayName', 'r_{cg,Y}');
    plot(axPos, t, r_cg(3, :), '-', 'Color', cZ, 'LineWidth', 1.6, 'DisplayName', 'r_{cg,Z}');
    plot(axPos, t, P_pos(3, :), '--', 'Color', cZ, 'LineWidth', 1.0, 'DisplayName', 'P_Z (CoP)');
    drawKnotLines(axPos, tK);
    ylabel(axPos, 'Position [m]', 'FontWeight', 'bold');
    title(axPos, 'CG Position r_{cg}(t) vs CoP Flat Output P(t)');
    legend(axPos, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 2: Velocity & Speed
    axVel = nexttile(tl2, 2);
    hold(axVel, 'on'); grid(axVel, 'on'); box(axVel, 'on');
    shadeSegments(axVel, QPTraj);
    vMag = vecnorm(v_cg, 2, 1);
    plot(axVel, t, v_cg(1, :), '-', 'Color', cX, 'LineWidth', 1.5, 'DisplayName', 'v_X');
    plot(axVel, t, v_cg(2, :), '-', 'Color', cY, 'LineWidth', 1.5, 'DisplayName', 'v_Y');
    plot(axVel, t, v_cg(3, :), '-', 'Color', cZ, 'LineWidth', 1.5, 'DisplayName', 'v_Z');
    plot(axVel, t, vMag, '-.', 'Color', cNorm, 'LineWidth', 1.4, 'DisplayName', '||v||');
    drawKnotLines(axVel, tK);
    ylabel(axVel, 'Velocity [m/s]', 'FontWeight', 'bold');
    title(axVel, sprintf('CG Velocity (Peak UtilVel = %.1f%%)', 100 * Report.UtilVel));
    legend(axVel, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 3: Acceleration
    axAcc = nexttile(tl2, 3);
    hold(axAcc, 'on'); grid(axAcc, 'on'); box(axAcc, 'on');
    shadeSegments(axAcc, QPTraj);
    aMag = vecnorm(a_cg, 2, 1);
    plot(axAcc, t, a_cg(1, :), '-', 'Color', cX, 'LineWidth', 1.5, 'DisplayName', 'a_{cg,X}');
    plot(axAcc, t, a_cg(2, :), '-', 'Color', cY, 'LineWidth', 1.5, 'DisplayName', 'a_{cg,Y}');
    plot(axAcc, t, a_cg(3, :), '-', 'Color', cZ, 'LineWidth', 1.5, 'DisplayName', 'a_{cg,Z}');
    plot(axAcc, t, aMag, '-.', 'Color', cNorm, 'LineWidth', 1.3, 'DisplayName', '||a_{cg}||');
    drawKnotLines(axAcc, tK);
    xlabel(axAcc, 'Time [s]', 'FontWeight', 'bold');
    ylabel(axAcc, 'Acceleration [m/s^2]', 'FontWeight', 'bold');
    title(axAcc, sprintf('CG Acceleration (Peak UtilAccel = %.1f%%)', 100 * Report.UtilAccel));
    legend(axAcc, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 4: Jerk & Snap Norm
    axJS = nexttile(tl2, 4);
    hold(axJS, 'on'); grid(axJS, 'on'); box(axJS, 'on');
    shadeSegments(axJS, QPTraj);
    yyaxis(axJS, 'left');
    plot(axJS, t, P_jerk(1, :), '-', 'Color', cX, 'LineWidth', 1.3, 'DisplayName', 'j_X');
    plot(axJS, t, P_jerk(2, :), '-', 'Color', cY, 'LineWidth', 1.3, 'DisplayName', 'j_Y');
    plot(axJS, t, P_jerk(3, :), '-', 'Color', cZ, 'LineWidth', 1.3, 'DisplayName', 'j_Z');
    ylabel(axJS, 'Jerk [m/s^3]', 'FontWeight', 'bold');
    axJS.YColor = cNorm;

    yyaxis(axJS, 'right');
    sMag = vecnorm(P_snap, 2, 1);
    plot(axJS, t, sMag, '--', 'Color', [0.55, 0.20, 0.70], 'LineWidth', 1.4, 'DisplayName', '||Snap||');
    ylabel(axJS, 'Snap Norm [m/s^4]', 'FontWeight', 'bold');
    axJS.YColor = [0.55, 0.20, 0.70];
    drawKnotLines(axJS, tK);
    xlabel(axJS, 'Time [s]', 'FontWeight', 'bold');
    title(axJS, 'CoP Jerk Components & Snap Norm ||d^4P/dt^4||');
    legend(axJS, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    linkaxes([axPos, axVel, axAcc, axJS], 'x');
    xlim(axPos, [0, QPTraj.T_total]);

    %% =====================================================================
    %% FIGURE 3: Attitude, Body Rates & Mass Depletion
    %% =====================================================================
    figHandles(3) = createDockedFig('3: Attitude, Body Rates & Mass Depletion', winStyle, figVis);
    tl3 = tiledlayout(figHandles(3), 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl3, sprintf('Attitude Kinematics, Angular Rates & Mass Depletion — %s', vehName), ...
        'FontWeight', 'bold', 'FontSize', 12);

    % Tile 1: Tilt Angle & Heading
    axTilt = nexttile(tl3, 1);
    hold(axTilt, 'on'); grid(axTilt, 'on'); box(axTilt, 'on');
    shadeSegments(axTilt, QPTraj);
    maxTiltDeg = rad2deg(Limits.MaxTilt);
    plot(axTilt, t, tiltDeg, '-', 'Color', cZ, 'LineWidth', 1.7, 'DisplayName', 'Tilt \theta_{tilt} [deg]');
    plot(axTilt, t, hdgDeg, '--', 'Color', cY, 'LineWidth', 1.3, 'DisplayName', 'Heading \psi [deg]');
    yline(axTilt, maxTiltDeg, '--', 'MaxTilt', 'Color', cLim, 'LineWidth', 1.3, 'DisplayName', 'MaxTilt Limit');
    drawKnotLines(axTilt, tK);
    ylabel(axTilt, 'Angle [deg]', 'FontWeight', 'bold');
    title(axTilt, sprintf('Vehicle Tilt & Heading (Peak UtilTilt = %.1f%%)', 100 * Report.UtilTilt));
    legend(axTilt, 'Location', 'best', 'FontSize', 8);

    % Tile 2: Quaternion Components
    axQuat = nexttile(tl3, 2);
    hold(axQuat, 'on'); grid(axQuat, 'on'); box(axQuat, 'on');
    shadeSegments(axQuat, QPTraj);
    plot(axQuat, t, q_cg(1, :), '-', 'Color', cNorm, 'LineWidth', 1.6, 'DisplayName', 'q_w');
    plot(axQuat, t, q_cg(2, :), '-', 'Color', cX, 'LineWidth', 1.4, 'DisplayName', 'q_x');
    plot(axQuat, t, q_cg(3, :), '-', 'Color', cY, 'LineWidth', 1.4, 'DisplayName', 'q_y');
    plot(axQuat, t, q_cg(4, :), '-', 'Color', cZ, 'LineWidth', 1.4, 'DisplayName', 'q_z');
    drawKnotLines(axQuat, tK);
    ylabel(axQuat, 'Quaternion [-]', 'FontWeight', 'bold');
    ylim(axQuat, [-1.05, 1.05]);
    title(axQuat, 'Attitude Unit Quaternion q_{IB}(t)');
    legend(axQuat, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 3: Body Angular Velocity omega_B
    axOmega = nexttile(tl3, 3);
    hold(axOmega, 'on'); grid(axOmega, 'on'); box(axOmega, 'on');
    shadeSegments(axOmega, QPTraj);
    wDeg = rad2deg(w_B);
    plot(axOmega, t, wDeg(1, :), '-', 'Color', cX, 'LineWidth', 1.5, 'DisplayName', '\omega_x (Roll)');
    plot(axOmega, t, wDeg(2, :), '-', 'Color', cY, 'LineWidth', 1.5, 'DisplayName', '\omega_y (Pitch)');
    plot(axOmega, t, wDeg(3, :), '-', 'Color', cZ, 'LineWidth', 1.5, 'DisplayName', '\omega_z (Yaw)');
    drawKnotLines(axOmega, tK);
    xlabel(axOmega, 'Time [s]', 'FontWeight', 'bold');
    ylabel(axOmega, 'Body Rate [deg/s]', 'FontWeight', 'bold');
    title(axOmega, 'Body Angular Rates \omega_B(t)');
    legend(axOmega, 'Location', 'best', 'FontSize', 8);

    % Tile 4: Mass Depletion & CGz Shift
    axMass = nexttile(tl3, 4);
    hold(axMass, 'on'); grid(axMass, 'on'); box(axMass, 'on');
    shadeSegments(axMass, QPTraj);
    yyaxis(axMass, 'left');
    plot(axMass, t, mTot, '-', 'Color', cNorm, 'LineWidth', 1.8, 'DisplayName', 'm_{total} [kg]');
    if isTOAD
        plot(axMass, t, Constants.m_dry + mLox, '--', 'Color', cZ, 'LineWidth', 1.2, 'DisplayName', 'm_{dry} + m_{LOX}');
        yline(axMass, Constants.m_dry, ':', 'm_{dry}', 'Color', cLim, 'LineWidth', 1.2, 'DisplayName', 'Dry Mass');
    end
    ylabel(axMass, 'Mass [kg]', 'FontWeight', 'bold');
    axMass.YColor = cNorm;

    yyaxis(axMass, 'right');
    plot(axMass, t, CGz_t * 100, '-.', 'Color', [0.85, 0.40, 0.10], 'LineWidth', 1.5, 'DisplayName', 'CG_z [cm]');
    ylabel(axMass, 'Moment Arm CG_z [cm]', 'FontWeight', 'bold');
    axMass.YColor = [0.85, 0.40, 0.10];
    drawKnotLines(axMass, tK);
    xlabel(axMass, 'Time [s]', 'FontWeight', 'bold');
    dmBurn = mTot(1) - mTot(end);
    title(axMass, sprintf('Mass Depletion (\\Deltam = %.2f kg) & Moment Arm CG_z', dmBurn));
    legend(axMass, 'Location', 'best', 'FontSize', 8);

    linkaxes([axTilt, axQuat, axOmega, axMass], 'x');
    xlim(axTilt, [0, QPTraj.T_total]);

    %% =====================================================================
    %% FIGURE 4: Actuators, Slew Rates & Constraint Utilization
    %% =====================================================================
    figHandles(4) = createDockedFig('4: Actuators, Slew Rates & Utilization', winStyle, figVis);
    tl4 = tiledlayout(figHandles(4), 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl4, sprintf('Actuator Allocations, Slew Rates & Constraint Utilization — %s', vehName), ...
        'FontWeight', 'bold', 'FontSize', 12);

    % Tile 1: Thrust Magnitude vs Limits & Hover Reference
    axThr = nexttile(tl4, 1);
    hold(axThr, 'on'); grid(axThr, 'on'); box(axThr, 'on');
    shadeSegments(axThr, QPTraj);
    if isTOAD
        minThrustHw = 0.50 * Limits.MaxThrust;
    else
        minThrustHw = 0.25 * Limits.MaxThrust;
    end
    plot(axThr, tU, U(3, :), '-', 'Color', cZ, 'LineWidth', 1.8, 'DisplayName', 'Thrust T(t)');
    plot(axThr, t, mTot * Limits.g, '-.', 'Color', [0.35, 0.35, 0.38], 'LineWidth', 1.3, 'DisplayName', 'Hover m(t)g');
    yline(axThr, Limits.T_hi, '--', 'T_{hi} (Eff)', 'Color', cLim, 'LineWidth', 1.3, 'DisplayName', 'Effective Bounds');
    yline(axThr, Limits.T_lo, '--', 'T_{lo} (Eff)', 'Color', cLim, 'LineWidth', 1.3, 'HandleVisibility', 'off');
    yline(axThr, Limits.MaxThrust, ':', 'MaxThrust', 'Color', [0.55, 0.15, 0.15], 'LineWidth', 1.0, 'DisplayName', 'Hardware Bounds');
    yline(axThr, minThrustHw, ':', 'MinThrust', 'Color', [0.55, 0.15, 0.15], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    drawKnotLines(axThr, tK);
    ylabel(axThr, 'Thrust [N]', 'FontWeight', 'bold');
    title(axThr, sprintf('Thrust Command (UtilHi = %.1f%%, UtilLo = %.1f%%)', ...
        100 * Report.UtilThrustHi, 100 * Report.UtilThrustLo));
    legend(axThr, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 2: Gimbal Deflections (theta, phi) & Roll Torque
    axGim = nexttile(tl4, 2);
    hold(axGim, 'on'); grid(axGim, 'on'); box(axGim, 'on');
    shadeSegments(axGim, QPTraj);
    gimLimEff = rad2deg(Limits.GimbalLimit);
    gimLimHw  = 10.0;
    plot(axGim, tU, rad2deg(U(1, :)), '-', 'Color', cX, 'LineWidth', 1.6, 'DisplayName', 'Gimbal \theta [deg]');
    plot(axGim, tU, rad2deg(U(2, :)), '-', 'Color', cY, 'LineWidth', 1.6, 'DisplayName', 'Gimbal \phi [deg]');
    yline(axGim,  gimLimEff, '--', '+\delta_{eff}', 'Color', cLim, 'LineWidth', 1.2, 'DisplayName', '\pm\delta_{eff} (8.5^\circ)');
    yline(axGim, -gimLimEff, '--', '-\delta_{eff}', 'Color', cLim, 'LineWidth', 1.2, 'HandleVisibility', 'off');
    yline(axGim,  gimLimHw,  ':',  '+\delta_{max}', 'Color', [0.55, 0.15, 0.15], 'LineWidth', 1.0, 'DisplayName', '\pm\delta_{hw} (10^\circ)');
    yline(axGim, -gimLimHw,  ':',  '-\delta_{max}', 'Color', [0.55, 0.15, 0.15], 'LineWidth', 1.0, 'HandleVisibility', 'off');
    drawKnotLines(axGim, tK);
    ylabel(axGim, 'Gimbal Angle [deg]', 'FontWeight', 'bold');
    title(axGim, sprintf('TVC Gimbal Angles (Peak UtilGimbal = %.1f%%)', 100 * Report.UtilGimbal));
    legend(axGim, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    % Tile 3: Actuator Slew Rates (Gimbal Rate & Thrust Rate)
    axSlew = nexttile(tl4, 3);
    hold(axSlew, 'on'); grid(axSlew, 'on'); box(axSlew, 'on');
    shadeSegments(axSlew, QPTraj);
    if Nc >= 2
        dTheta = diff(rad2deg(U(1, :))) / dt;
        dPhi   = diff(rad2deg(U(2, :))) / dt;
        dThr   = diff(U(3, :)) / dt;
        tSlew  = tU(1:(Nc - 1));
    else
        dTheta = 0; dPhi = 0; dThr = 0; tSlew = 0;
    end
    yyaxis(axSlew, 'left');
    plot(axSlew, tSlew, dTheta, '-', 'Color', cX, 'LineWidth', 1.4, 'DisplayName', 'd\theta/dt [deg/s]');
    plot(axSlew, tSlew, dPhi,   '-', 'Color', cY, 'LineWidth', 1.4, 'DisplayName', 'd\phi/dt [deg/s]');
    yline(axSlew,  rad2deg(Limits.GimbalRateLim), '--', 'Color', cLim, 'LineWidth', 1.1, 'DisplayName', 'MaxGimbalRate');
    yline(axSlew, -rad2deg(Limits.GimbalRateLim), '--', 'Color', cLim, 'LineWidth', 1.1, 'HandleVisibility', 'off');
    ylabel(axSlew, 'Gimbal Slew [deg/s]', 'FontWeight', 'bold');
    axSlew.YColor = cNorm;

    yyaxis(axSlew, 'right');
    plot(axSlew, tSlew, dThr, '-.', 'Color', cZ, 'LineWidth', 1.3, 'DisplayName', 'dT/dt [N/s]');
    ylabel(axSlew, 'Thrust Slew [N/s]', 'FontWeight', 'bold');
    axSlew.YColor = cZ;
    drawKnotLines(axSlew, tK);
    xlabel(axSlew, 'Time [s]', 'FontWeight', 'bold');
    title(axSlew, sprintf('Actuator Slew Rates (UtilGimRate = %.1f%%, UtilThrRate = %.1f%%)', ...
        100 * Report.UtilGimbalRate, 100 * Report.UtilThrustRate));
    legend(axSlew, 'Location', 'best', 'NumColumns', 2, 'FontSize', 8);

    linkaxes([axThr, axGim, axSlew], 'x');
    xlim(axThr, [0, QPTraj.T_total]);

    % Tile 4: Normalized Constraint Utilization Summary Bar Chart
    axUtil = nexttile(tl4, 4);
    hold(axUtil, 'on'); grid(axUtil, 'on'); box(axUtil, 'on');
    utilNames = {'Velocity', 'Acceleration', 'Thrust High', 'Thrust Low', ...
                 'Gimbal Angle', 'Gimbal Rate', 'Thrust Rate', 'Roll Torque', 'Tilt Angle'};
    utilVals = 100 * [Report.UtilVel, Report.UtilAccel, Report.UtilThrustHi, Report.UtilThrustLo, ...
                      Report.UtilGimbal, Report.UtilGimbalRate, Report.UtilThrustRate, ...
                      Report.UtilRoll, Report.UtilTilt];
    b = barh(axUtil, 1:9, utilVals, 0.62, 'FaceColor', 'flat');
    for i = 1:9
        if utilVals(i) > 100.1
            b.CData(i, :) = [0.85, 0.20, 0.20];
        elseif utilVals(i) > 85.0
            b.CData(i, :) = [0.92, 0.58, 0.12];
        else
            b.CData(i, :) = [0.18, 0.52, 0.82];
        end
        text(axUtil, utilVals(i) + 1.8, i, sprintf('%.1f%%', utilVals(i)), ...
            'VerticalAlignment', 'middle', 'FontSize', 8, 'FontWeight', 'bold');
    end
    xline(axUtil, 100.0, '--', '100% Limit', 'Color', cLim, 'LineWidth', 1.5, ...
        'LabelVerticalAlignment', 'bottom', 'FontWeight', 'bold');
    yticks(axUtil, 1:9);
    yticklabels(axUtil, utilNames);
    set(axUtil, 'YDir', 'reverse');
    xlim(axUtil, [0, max(118, max(utilVals) + 15)]);
    xlabel(axUtil, 'Peak Constraint Utilization [%]', 'FontWeight', 'bold');
    title(axUtil, sprintf('Constraint Utilization Summary (Outer Dilation Iters = %d)', Report.OuterIterations));
end

function fig = createDockedFig(figName, winStyle, figVis)
%% Create Clean White Docked Figure Tab
    fig = figure('Name', figName, 'NumberTitle', 'off', ...
        'Color', 'w', 'WindowStyle', winStyle, 'Visible', figVis);
    if ~strcmpi(winStyle, 'docked')
        fig.Position = [80, 60, 1420, 840];
    end
end

function plotSegmented3D(ax, t, r_cg, QPTraj)
%% Plot 3D Trajectory Colored by Segment Mode (Hold, Cruise, Flight)
    M  = QPTraj.Count;
    tK = QPTraj.KeyTimes(1:(M + 1));
    drewFlight = false;
    drewCruise = false;
    drewHold   = false;

    for m = 1:M
        mask = (t >= tK(m) - 1e-9) & (t <= tK(m + 1) + 1e-9);
        if QPTraj.IsHold(m)
            lbl = 'Hold Segment';
            vis = ~drewHold; drewHold = true;
            plot3(ax, r_cg(1, mask), r_cg(2, mask), r_cg(3, mask), 'o', ...
                'Color', [0.85, 0.20, 0.20], 'MarkerFaceColor', [0.85, 0.20, 0.20], ...
                'MarkerSize', 6, 'DisplayName', lbl, 'HandleVisibility', onOff(vis));
        elseif QPTraj.IsCruise(m)
            lbl = 'Constant-Vel Cruise';
            vis = ~drewCruise; drewCruise = true;
            plot3(ax, r_cg(1, mask), r_cg(2, mask), r_cg(3, mask), '-', ...
                'Color', [0.12, 0.65, 0.30], 'LineWidth', 3.0, ...
                'DisplayName', lbl, 'HandleVisibility', onOff(vis));
        else
            lbl = 'PolyQP Flight Segment';
            vis = ~drewFlight; drewFlight = true;
            plot3(ax, r_cg(1, mask), r_cg(2, mask), r_cg(3, mask), '-', ...
                'Color', [0.12, 0.38, 0.82], 'LineWidth', 2.4, ...
                'DisplayName', lbl, 'HandleVisibility', onOff(vis));
        end
        % Segment knot marker
        idxKnot = find(mask, 1, 'first');
        if ~isempty(idxKnot)
            plot3(ax, r_cg(1, idxKnot), r_cg(2, idxKnot), r_cg(3, idxKnot), 's', ...
                'Color', [0.20, 0.20, 0.25], 'MarkerFaceColor', 'w', 'MarkerSize', 5, ...
                'HandleVisibility', 'off');
        end
    end
end

function plotWaypoints3D(ax, Waypoints, span3D)
%% Overlay Waypoint Pins and Labels in 3D
    nW = numel(Waypoints);
    dx = 0.025 * span3D;
    for w = 1:nW
        pos = Waypoints(w).PosBC;
        if all(~isnan(pos))
            plot3(ax, pos(1), pos(2), pos(3), 'p', 'Color', [0.80, 0.15, 0.15], ...
                'MarkerFaceColor', [0.95, 0.30, 0.25], 'MarkerSize', 10, ...
                'DisplayName', 'Waypoints', 'HandleVisibility', onOff(w == 1));
            lbl = char(Waypoints(w).Label);
            if isempty(lbl)
                lbl = sprintf('WP %d', w);
            end
            if Waypoints(w).HoldTime > 0
                lbl = sprintf('%s (%.0fs)', lbl, Waypoints(w).HoldTime);
            end
            text(ax, pos(1) + dx, pos(2) + dx, pos(3) + dx, lbl, ...
                'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.15, 0.15, 0.20], ...
                'BackgroundColor', [1, 1, 1, 0.75], 'Margin', 1);
        end
    end
end

function plotOrtho2D(ax, x_cg, y_cg, x_p, y_p, Waypoints, idx1, idx2, xLbl, yLbl, ttl)
%% Plot 2D Orthographic Projection with Waypoint Pins
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    plot(ax, x_p, y_p, ':', 'Color', [0.55, 0.35, 0.75], 'LineWidth', 1.1);
    plot(ax, x_cg, y_cg, '-', 'Color', [0.12, 0.38, 0.82], 'LineWidth', 1.8);
    for w = 1:numel(Waypoints)
        pos = Waypoints(w).PosBC;
        if ~isnan(pos(idx1)) && ~isnan(pos(idx2))
            plot(ax, pos(idx1), pos(idx2), 'o', 'Color', [0.80, 0.15, 0.15], ...
                'MarkerFaceColor', [0.95, 0.30, 0.25], 'MarkerSize', 5);
        end
    end
    xlabel(ax, xLbl, 'FontSize', 8, 'FontWeight', 'bold');
    ylabel(ax, yLbl, 'FontSize', 8, 'FontWeight', 'bold');
    title(ax, ttl, 'FontSize', 9, 'FontWeight', 'bold');
    axis(ax, 'equal');
end

function shadeSegments(ax, QPTraj)
%% Light Background Shading for Hold and Constant-Velocity Cruise Segments
    M  = QPTraj.Count;
    tK = QPTraj.KeyTimes(1:(M + 1));
    for m = 1:M
        if QPTraj.IsHold(m)
            xregion(ax, tK(m), tK(m + 1), 'FaceColor', [0.92, 0.82, 0.82], ...
                'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        elseif QPTraj.IsCruise(m)
            xregion(ax, tK(m), tK(m + 1), 'FaceColor', [0.80, 0.94, 0.84], ...
                'FaceAlpha', 0.28, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
    end
end

function drawKnotLines(ax, tK)
%% Draw Vertical Dashed Knot Lines at Segment Transitions
    for k = 2:(numel(tK) - 1)
        xline(ax, tK(k), ':', 'Color', [0.65, 0.65, 0.70], ...
            'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

function s = onOff(flag)
    if flag
        s = 'on';
    else
        s = 'off';
    end
end
