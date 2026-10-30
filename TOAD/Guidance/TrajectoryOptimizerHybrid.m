classdef TrajectoryOptimizerHybrid < handle
    %% TrajectoryOptimizerHybrid - Distilled Two-Stage 6-DoF Optimal Trajectory Generator.
    % High-performance, mathematically rigorous hybrid guidance engine combining:
    %   Stage 1: Minimum-snap polynomial trajectory optimization (PolyQP)
    %            with corridor enforcement & differential-flatness state inversion.
    %   Stage 2: 6-DoF direct collocation refinement (IPOPT / CasADi RK4) with
    %            actuator chatter suppression and takeoff/landing flight corridors.
    %
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    %% Properties
    properties
        constants; Vehicle double = 0; N double = 60
        T_bounds = [8, 35]; T_initial double = 16
        CircleTightness double = 0.5; GlideslopeAngle double = 5.0; FunnelCurvature double = 0.01
        TakeoffClearanceAlt double = 0; LandingFlareAlt double = 0
        SigPos double = 0.30; SigVel double = 0.40; SigAtt double = deg2rad(5)
        SigT double = 0.03; TightFactor double = 0.25

        %% Costs
        WeightTrackPos double = 1.00; WeightTrackVel double = 0.50; WeightTrackAtt double = 0.1
        WeightTrackTime double = 0; WeightCorridor double = 1.00
        WeightSlew double = 0.1; WeightCurv double = 0.1

        Maneuver string = "Circle"; ManeuverParams struct = struct()
        Waypoints double = []; WaypointTolerances double = []
        q0 = [1; 0; 0; 0]; r0 = [0; 0; 0]; v0 = [0; 0; 0]; w0 = [0; 0; 0]; r_f = [0; 0; 0]
        LandingVelTol double = 0.35; v_f_tol double = 0.35
        ThrustMargin double = 0.05; GimbalMargin double = 0.15
        MaxGimbalRate double = deg2rad(35); MaxThrustRate double = 1500
        MaxRollRate double = 5; MaxGimbalAngle double = pi/15
        Sx double; Su double; L_c double = 50
        MaxIter double = 120; MaxCpuTime double = 45.0; Tol double = 1.5e-2
        ConstrViolTol double = 5e-3; PrintLevel double = 0
        SaveDir string = ""; Version double = 1; Filename string = ""
        Stage1Diag struct = struct(); QPTraj; InitialGuess struct = struct()
        Solution struct = struct(); OptiVars struct = struct()
        RecordProgression logical = false; ProgressionLogger; Progression struct = struct()
    end

    %% Methods
    methods
        %% Object Lifecycle & Configuration
        function obj = TrajectoryOptimizerHybrid(constants6DoF, varargin)
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizerHybrid:MissingConstants', 'constants6DoF struct required.');
            end
            obj.constants = constants6DoF;
            if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir') && isempty(which('casadi.Opti'))
                addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                if isprop(obj, name), obj.(name) = varargin{k+1}; end
            end
            obj.v_f_tol = obj.LandingVelTol;
            obj.Maneuver = string(obj.Maneuver);

            % State and control scaling normalization vectors
            V_c = 15; Omega_c = 2.0;
            if obj.Vehicle == 0
                obj.Sx = [1; 1; 1; 1; obj.L_c; obj.L_c; obj.L_c; ...
                          V_c; V_c; V_c; Omega_c; Omega_c; Omega_c; 1; 1];
            else
                obj.Sx = [1; 1; 1; 1; obj.L_c; obj.L_c; obj.L_c; ...
                          V_c; V_c; V_c; Omega_c; Omega_c; Omega_c; ...
                          obj.constants.OxMass; obj.constants.FuMass];
            end
            obj.Su = [obj.MaxGimbalAngle; obj.MaxGimbalAngle; ...
                      obj.constants.MaxThrust; obj.MaxRollRate];

            obj.setManeuver(obj.Maneuver);
        end

        function d = getVehicleDefaults(obj)
            %% GETVEHICLEDEFAULTS  Centralized vehicle-specific physical heuristics.
            d = struct('MinThrottle', 0.25, 'v_desc_lim', 3.5, 'h_to_clear', 1.0, ...
                'h_ld_clear', 1.2, 'dR', 5.0, 'dA', 7.0, 'dHA', 15.0, 'dBA', 20.0, ...
                'accel_est', 3.5, 't_base_pad', 3.0);
            if obj.Vehicle == 1
                d = struct('MinThrottle', 0.50, 'v_desc_lim', 6.0, 'h_to_clear', 2.5, ...
                    'h_ld_clear', 3.0, 'dR', 15.0, 'dA', 25.0, 'dHA', 50.0, 'dBA', 45.0, ...
                    'accel_est', 2.2, 't_base_pad', 4.0);
            end
        end

        function h = getTakeoffClearanceAlt(obj)
            if obj.TakeoffClearanceAlt > 0, h = obj.TakeoffClearanceAlt;
            else, d = obj.getVehicleDefaults(); h = d.h_to_clear; end
        end

        function h = getLandingFlareAlt(obj)
            if obj.LandingFlareAlt > 0, h = obj.LandingFlareAlt;
            else, d = obj.getVehicleDefaults(); h = d.h_ld_clear; end
        end

        function setBoundaries(obj, r0_in, r_f_in, varargin)
            %% SETBOUNDARIES  Configure launch pad and landing target coordinates.
            obj.r0 = r0_in(:); obj.r_f = r_f_in(:);
            p = inputParser;
            addParameter(p, 'v0', [0; 0; 0]); addParameter(p, 'q0', [1; 0; 0; 0]);
            addParameter(p, 'LandingVelTol', obj.LandingVelTol);
            addParameter(p, 'v_f_tol', obj.LandingVelTol);
            parse(p, varargin{:});
            obj.v0 = p.Results.v0(:); obj.q0 = p.Results.q0(:);
            if ~ismember('v_f_tol', p.UsingDefaults), obj.LandingVelTol = double(p.Results.v_f_tol);
            else, obj.LandingVelTol = double(p.Results.LandingVelTol); end
            obj.v_f_tol = obj.LandingVelTol;
            obj.updateDurationEstimate();
            obj.QPTraj = []; obj.InitialGuess = struct();
        end

        function updateDurationEstimate(obj)
            %% UPDATEDURATIONESTIMATE  Kinematics-based baseline duration and bounds.
            if obj.Maneuver == "Hop" && isstruct(obj.ManeuverParams) && isfield(obj.ManeuverParams, 'apex_alt')
                dh = max(1.5, obj.ManeuverParams.apex_alt - max(obj.r0(3), obj.r_f(3)));
                d_h = norm(obj.r_f(1:2) - obj.r0(1:2));
                vd = obj.getVehicleDefaults();
                t_phys = 2.0 * sqrt(2.0 * dh / vd.accel_est) + max(vd.t_base_pad, sqrt(d_h));
                if obj.Vehicle == 1
                    obj.T_initial = max(20.0, min(36.0, round(t_phys, 1)));
                    obj.T_bounds  = [max(14.0, round(0.65 * obj.T_initial, 1)), ...
                                     min(45.0, round(1.60 * obj.T_initial, 1))];
                else
                    obj.T_initial = max(10.0, min(24.0, round(t_phys, 1)));
                    obj.T_bounds  = [max(8.0, round(0.65 * obj.T_initial, 1)), ...
                                     min(28.0, round(1.50 * obj.T_initial, 1))];
                end
            end
        end

        function setManeuver(obj, name, varargin)
            %% SETMANEUVER  Configure standard flight maneuvers and geometry presets.
            obj.Maneuver = string(name);
            p = inputParser; vd = obj.getVehicleDefaults();
            switch obj.Maneuver
                case "Circle"
                    addParameter(p, 'circle_radius', vd.dR);
                    addParameter(p, 'circle_alt', vd.dA);
                    addParameter(p, 'circle_center', [0; 0]);
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('circle_radius', double(p.Results.circle_radius), ...
                        'circle_alt', double(p.Results.circle_alt), ...
                        'circle_center', double(p.Results.circle_center(:)), ...
                        'f_orbit_start', 0.25, 'f_orbit_end', 0.75);
                    if obj.Vehicle == 0, obj.T_bounds = [12, 50]; obj.T_initial = 20;
                    else,                obj.T_bounds = [15, 60]; obj.T_initial = 30; end
                case "Backflip"
                    addParameter(p, 'apex_alt', vd.dBA);
                    addParameter(p, 'flip_start_frac', 0.30);
                    addParameter(p, 'flip_end_frac', 0.70);
                    addParameter(p, 'theta_tol', deg2rad(25));
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt), ...
                        'flip_start_frac', double(p.Results.flip_start_frac), ...
                        'flip_end_frac', double(p.Results.flip_end_frac), ...
                        'theta_tol', double(p.Results.theta_tol), 'q_inverted', [0; 0; -1; 0]);
                    if obj.Vehicle == 0, obj.T_bounds = [12, 50]; obj.T_initial = 20;
                    else,                obj.T_bounds = [15, 60]; obj.T_initial = 30; end
                case "Hop"
                    addParameter(p, 'apex_alt', vd.dHA);
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt));
                    obj.updateDurationEstimate();
                case {"Waypoint", "Custom"}
                    addParameter(p, 'Waypoints', []); addParameter(p, 'Tolerances', []);
                    parse(p, varargin{:});
                    if ~isempty(p.Results.Waypoints)
                        obj.setWaypoints(p.Results.Waypoints, 'Tolerances', p.Results.Tolerances);
                    else
                        obj.ManeuverParams = struct();
                    end
                otherwise
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct();
            end
            obj.QPTraj = []; obj.InitialGuess = struct();
        end

        function setWaypoints(obj, waypoints, varargin)
            %% SETWAYPOINTS  Configure arbitrary 3D multi-waypoint trajectory.
            if size(waypoints, 1) ~= 3 || size(waypoints, 2) < 2
                error('TrajectoryOptimizerHybrid:InvalidWaypoints', 'Waypoints must be a 3 x K matrix with K >= 2.');
            end
            obj.Waypoints = double(waypoints);
            K = size(obj.Waypoints, 2);
            p = inputParser;
            addParameter(p, 'Tolerances', 0.25 * ones(1, K));
            addParameter(p, 'T_total', obj.T_initial);
            parse(p, varargin{:});
            obj.WaypointTolerances = double(p.Results.Tolerances);
            if isscalar(obj.WaypointTolerances), obj.WaypointTolerances = repmat(obj.WaypointTolerances, 1, K); end
            obj.T_initial = double(p.Results.T_total);
            obj.T_bounds = [max(4.0, round(0.50 * obj.T_initial, 1)), ...
                            min(70.0, round(2.0 * obj.T_initial, 1))];
            obj.Maneuver = "Waypoint"; obj.ManeuverParams = struct();
            obj.r0 = obj.Waypoints(:, 1); obj.r_f = obj.Waypoints(:, end);
            obj.QPTraj = []; obj.InitialGuess = struct();
        end

        %% Stage 1: Polynomial QP Engine
        function [viol, tgt] = evalTakeoffCone(obj, ts, rs, tg)
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts); t_end = 0.15 * t_max;
            if isfield(tg, 'orb_entry'), t_end = min(0.20 * t_max, tg.orb_entry);
            elseif isfield(tg, 'wp_2'),  t_end = min(0.20 * t_max, tg.wp_2); end

            r0p = obj.r0; gsa = tand(obj.GlideslopeAngle);
            h_liftoff = obj.getTakeoffClearanceAlt();
            for i = 1:n
                if ts(i) > t_end, continue; end
                dz = rs(3,i) - r0p(3);
                if dz < 0.01 || dz > h_liftoff, continue; end
                rxy = sqrt((rs(1,i)-r0p(1))^2 + (rs(2,i)-r0p(2))^2 + 1e-4);
                lim = dz * gsa + 0.15; v = rxy - lim;
                if v > 0
                    viol(i) = v; sc = lim / max(1e-6, rxy);
                    tgt(1:2, i) = r0p(1:2) + (rs(1:2, i) - r0p(1:2)) * sc;
                end
            end
        end

        function [viol, tgt] = evalLandingFunnel(obj, ts, rs, tg)
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts); t_start = 0.70 * t_max;
            if isfield(tg, 'orb_exit'),   t_start = max(0.65 * t_max, tg.orb_exit);
            elseif isfield(tg, 'apex'),   t_start = max(0.60 * t_max, tg.apex);
            elseif ~isempty(obj.Waypoints) && size(obj.Waypoints, 2) >= 2
                tag_penult = sprintf('wp_%d', size(obj.Waypoints, 2) - 1);
                if isfield(tg, tag_penult), t_start = tg.(tag_penult); end
            end

            gsa = tand(obj.GlideslopeAngle); cf = obj.FunnelCurvature;
            rfp = obj.r_f; h_flare = obj.getLandingFlareAlt();
            for i = 1:n
                if ts(i) < t_start, continue; end
                dz = rs(3,i) - rfp(3);
                if dz < 0.01 || dz > h_flare, continue; end
                rxy = sqrt((rs(1,i)-rfp(1))^2 + (rs(2,i)-rfp(2))^2 + 1e-4);
                lim = dz * gsa + cf * dz^2 + 0.15; v = rxy - lim;
                if v > 0
                    viol(i) = v; sc = lim / max(1e-6, rxy);
                    tgt(1:2, i) = rfp(1:2) + (rs(1:2, i) - rfp(1:2)) * sc;
                end
            end
        end

        function guess = buildInitialGuess(obj)
            %% BUILDINITIALGUESS  Sample QP trajectory and invert via differential flatness.
            t0 = tic;

            % 1. Convert maneuver geometry to keyframe arrays
            R = 4; vL = [0; 0; -0.10];
            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    Rad = mp.circle_radius; h = mp.circle_alt;
                    N_orb = 16; angles = linspace(0, 2*pi, N_orb + 1);
                    K = 1 + (N_orb + 1) + 1;
                    Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    for j = 1:(N_orb + 1)
                        Fix(:,1,1 + j) = [cx + Rad*cos(angles(j)); cy + Rad*sin(angles(j)); h];
                    end
                    Tag{2} = 'orb_entry'; Tag{N_orb + 2} = 'orb_exit';
                    Fix(:,1,K) = obj.r_f; Fix(:,2,K) = vL; Tag{K} = 'land';
                    M = K - 1; Ts = obj.T_initial;
                    T_seg = zeros(1, M); T_seg(1) = 0.25 * Ts;
                    for j = 1:N_orb, T_seg(1 + j) = 0.50 * Ts / N_orb; end
                    T_seg(M) = 0.25 * Ts;

                case {"Hop", "Backflip"}
                    mp = obj.ManeuverParams;
                    apex = [0.5*(obj.r0(1)+obj.r_f(1)); 0.5*(obj.r0(2)+obj.r_f(2)); mp.apex_alt];
                    K = 3; Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    Fix(:,1,2) = apex;  Fix(3,2,2) = 0;      Tag{2} = 'apex';
                    Fix(:,1,3) = obj.r_f; Fix(:,2,3) = vL;   Tag{3} = 'land';
                    T_seg = [0.50, 0.50] * obj.T_initial;

                case {"Waypoint", "Custom"}
                    if isempty(obj.Waypoints)
                        obj.Waypoints = [obj.r0, 0.5*(obj.r0 + obj.r_f) + [0; 0; 10], obj.r_f];
                    end
                    K = size(obj.Waypoints, 2);
                    Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.Waypoints(:,1); Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    for k = 2:K-1
                        Fix(:,1,k) = obj.Waypoints(:,k);
                        Tag{k} = sprintf('wp_%d', k);
                    end
                    Fix(:,1,K) = obj.Waypoints(:,K); Fix(:,2,K) = vL; Tag{K} = 'land';
                    M = K - 1; dists = zeros(1, M);
                    for m = 1:M, dists(m) = norm(obj.Waypoints(:,m+1) - obj.Waypoints(:,m)); end
                    td = sum(dists);
                    if td < 1e-3, T_seg = (obj.T_initial / M) * ones(1, M);
                    else,         T_seg = max(1.5, obj.T_initial * dists / td); T_seg = obj.T_initial * T_seg / sum(T_seg); end

                otherwise
                    K = 2; Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    Fix(:,1,2) = obj.r_f; Fix(:,2,2) = vL; Tag{2} = 'land';
                    T_seg = obj.T_initial;
            end
            corridors = { ...
                struct('fn', @(ts,rs,tg) obj.evalTakeoffCone(ts,rs,tg)), ...
                struct('fn', @(ts,rs,tg) obj.evalLandingFunnel(ts,rs,tg)) ...
            };
            do_opt_times = (obj.Maneuver ~= "Circle" && numel(T_seg) >= 2);

            opts = struct('Tag', {Tag}, 'Corridors', {corridors}, 'N', obj.N, ...
                'Vehicle', obj.Vehicle, 'Maneuver', obj.Maneuver, ...
                'ManeuverParams', obj.ManeuverParams, ...
                'q0', obj.q0, 'w0', obj.w0, 'r0', obj.r0, 'rf', obj.r_f, ...
                'ThrustMargin', obj.ThrustMargin, 'GimbalMargin', obj.GimbalMargin, ...
                'MaxGimbalAngle', obj.MaxGimbalAngle, 'MaxGimbalRate', obj.MaxGimbalRate, ...
                'MaxThrustRate', obj.MaxThrustRate, 'MaxRollRate', obj.MaxRollRate, ...
                'OptimizeTimes', do_opt_times);

            %% Call polynomial optimizer
            [obj.InitialGuess, obj.QPTraj] = PolyQP(Fix, T_seg, obj.constants, opts);

            tg = obj.QPTraj.TagTimes; Ttot = obj.InitialGuess.T_total;
            if isfield(tg, 'orb_entry'), obj.ManeuverParams.f_orbit_start = tg.orb_entry / Ttot; end
            if isfield(tg, 'orb_exit'),  obj.ManeuverParams.f_orbit_end   = tg.orb_exit / Ttot;  end
            if isfield(tg, 'apex'),      obj.ManeuverParams.f_apex        = tg.apex / Ttot;      end

            obj.T_initial = Ttot;
            obj.InitialGuess.Xhat = obj.InitialGuess.X ./ obj.Sx;
            obj.InitialGuess.Uhat = obj.InitialGuess.U ./ obj.Su;

            obj.Stage1Diag.t_qp_ms = toc(t0) * 1000;
            obj.Stage1Diag.t_stage1_ms = toc(t0) * 1000;
            obj.Stage1Diag.T_total = Ttot;
            obj.Stage1Diag.max_thrust = max(obj.InitialGuess.U(3,:));
            obj.Stage1Diag.max_gimbal_deg = rad2deg(max(sqrt(obj.InitialGuess.U(1,:).^2 + obj.InitialGuess.U(2,:).^2)));

            guess = obj.InitialGuess;

            if obj.RecordProgression
                obj.ensureProgressionLogger();
                obj.ProgressionLogger.recordQP('Stage 1A: QP Baseline Spline', obj.QPTraj, Fix);
                obj.ProgressionLogger.recordGuess('Stage 1D: Flatness Inversion (Initial Guess)', guess);
                obj.Progression = obj.ProgressionLogger.Progression;
            end
        end

        %% Stage 2: Direct Collocation
        function sol = solve(obj)
            %% SOLVE  Two-Stage Pipeline: QP Polynomial Seed -> 6-DoF Direct Collocation
            import casadi.*
            t_total_start = tic;

            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'Xhat')
                obj.buildInitialGuess();
            end

            N_nodes = obj.N; opti = Opti(); dyn_fnc = obj.getCasADiDynamics();
            obj.OptiVars = struct();

            % Decision variables
            T_total = opti.variable();
            opti.subject_to(obj.T_bounds(1) <= T_total);
            opti.subject_to(T_total <= obj.T_bounds(2));
            opti.set_initial(T_total, obj.InitialGuess.T_total);
            dt_row = repmat(T_total / N_nodes, 1, N_nodes);

            Xhat = opti.variable(15, N_nodes + 1); Uhat = opti.variable(4, N_nodes);
            X = obj.Sx .* Xhat; U = obj.Su .* Uhat;

            % RK4 Discretization Map
            params_val = [obj.constants.m_dry; obj.constants.g; obj.constants.rTB; ...
                obj.constants.Ox_Z; obj.constants.OxMass; obj.constants.OxHeight; ...
                obj.constants.Fu_Z; obj.constants.FuMass; obj.constants.FuHeight; ...
                obj.constants.J(:); obj.constants.OxRadius; obj.constants.FuRadius; ...
                obj.constants.MaxThrust; obj.constants.OF; obj.constants.MaxMdot];

            xs = MX.sym('xh', 15); us = MX.sym('uh', 4); ds = MX.sym('dt');
            xp = obj.Sx .* xs; up = obj.Su .* us;
            k1 = dyn_fnc(xp,               up, params_val);
            k2 = dyn_fnc(xp + ds/2 * k1,   up, params_val);
            k3 = dyn_fnc(xp + ds/2 * k2,   up, params_val);
            k4 = dyn_fnc(xp + ds * k3,     up, params_val);
            xn = xp + ds / 6 * (k1 + 2*k2 + 2*k3 + k4);
            xn = [xn(1:4) / sqrt(sum(xn(1:4).^2) + 1e-12); xn(5:end)];
            Fstep = Function('F_step', {xs, us, ds}, {xn ./ obj.Sx});
            Fmap  = Fstep.map(N_nodes);
            opti.subject_to(Xhat(:, 2:end) == Fmap(Xhat(:,1:N_nodes), Uhat, dt_row));

            % Horizontal and altitude bounds
            max_h = max([norm(obj.r0(1:2)), norm(obj.r_f(1:2)), norm(obj.r_f(1:2) - obj.r0(1:2)), 30.0]);
            if ~isempty(obj.Waypoints), max_h = max(max_h, max(sqrt(sum(obj.Waypoints(1:2,:).^2, 1)))); end
            Rh = max(150.0, 3.0 * max_h);
            opti.subject_to(-Rh <= X(5,:)); opti.subject_to(X(5,:) <= Rh);
            opti.subject_to(-Rh <= X(6,:)); opti.subject_to(X(6,:) <= Rh);

            alts = [obj.r0(3), obj.r_f(3), 20.0];
            if isfield(obj.ManeuverParams, 'apex_alt'),   alts = [alts, obj.ManeuverParams.apex_alt];   end
            if isfield(obj.ManeuverParams, 'circle_alt'), alts = [alts, obj.ManeuverParams.circle_alt]; end
            if ~isempty(obj.Waypoints),                   alts = [alts, obj.Waypoints(3,:)];            end
            opti.subject_to(0 <= X(7,:)); opti.subject_to(X(7,:) <= max(180.0, 2.5 * max(alts)));

            % Boundary constraints
            m_lox0 = obj.constants.OxMass; m_ipa0 = obj.constants.FuMass;
            opti.subject_to(X(:,1) == [obj.q0; obj.r0; obj.v0; obj.w0; m_lox0; m_ipa0]);
            if obj.Maneuver == "Backflip"
                opti.subject_to(X(2:4,end) == -obj.q0(2:4));
                opti.subject_to(X(1,end) <= -0.5);
            else
                opti.subject_to(X(2:4,end) == obj.q0(2:4));
                opti.subject_to(X(1,end) >= 0.5);
            end
            opti.subject_to(X(5:7,end) == obj.r_f);
            opti.subject_to(sum(X(8:10,end).^2) <= obj.LandingVelTol^2);
            if obj.Vehicle == 1
                opti.subject_to(X(14,end) >= 0.02*m_lox0);
                opti.subject_to(X(15,end) >= 0.02*m_ipa0);
            end

            % Control bounds & slew rates
            MT = obj.constants.MaxThrust; tm = obj.ThrustMargin;
            gm = obj.GimbalMargin; mg = obj.MaxGimbalAngle;
            vd = obj.getVehicleDefaults();

            opti.subject_to((vd.MinThrottle + tm)*MT <= U(3,:));
            opti.subject_to(U(3,:) <= (1 - tm)*MT);
            opti.subject_to(-(1 - gm)*mg <= U(1,:)); opti.subject_to(U(1,:) <= (1 - gm)*mg);
            opti.subject_to(-(1 - gm)*mg <= U(2,:)); opti.subject_to(U(2,:) <= (1 - gm)*mg);
            opti.subject_to(-(1 - tm)*obj.MaxRollRate <= U(4,:));
            opti.subject_to(U(4,:) <= (1 - tm)*obj.MaxRollRate);

            dU = U(:,2:end) - U(:,1:end-1); dts = dt_row(1:end-1);
            for ch = [1, 2]
                opti.subject_to(-obj.MaxGimbalRate * dts <= dU(ch,:));
                opti.subject_to(dU(ch,:) <= obj.MaxGimbalRate * dts);
            end
            opti.subject_to(-obj.MaxThrustRate * dts <= dU(3,:));
            opti.subject_to(dU(3,:) <= obj.MaxThrustRate * dts);
            opti.subject_to(-obj.MaxRollRate * dts <= dU(4,:));
            opti.subject_to(dU(4,:) <= obj.MaxRollRate * dts);

            % Decoupled corridor metrics
            [viol_to, viol_ld, tight_idx] = obj.computeCorridorMetrics(X);

            % Apply physics constraints & projection cost
            obj.applyManeuverConstraints(opti, X, tight_idx);
            obj.applyCostFunction(opti, X, U, T_total, viol_to, viol_ld, tight_idx);

            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);

            % Solver configuration
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, 'tol', obj.Tol, ...
                'constr_viol_tol', obj.ConstrViolTol, 'dual_inf_tol', 1.0, ...
                'acceptable_tol', max(obj.Tol, 2e-2), 'acceptable_constr_viol_tol', 1e-3, ...
                'acceptable_dual_inf_tol', 1e5, 'acceptable_iter', 2, ...
                'max_cpu_time', obj.MaxCpuTime, 'mu_strategy', 'adaptive', 'mu_init', 1e-3, ...
                'bound_push', 1e-6, 'bound_frac', 1e-6, 'print_level', obj.PrintLevel);
            if obj.Vehicle == 0 && (obj.Maneuver == "Waypoint" || obj.Maneuver == "Custom")
                s_opts.hessian_approximation = 'limited-memory';
            end
            if obj.RecordProgression
                obj.ensureProgressionLogger();
                opti.callback(@(iter) obj.ProgressionLogger.recordCollocationIter(opti, Xhat, Uhat, T_total, iter));
            end
            opti.solver('ipopt', p_opts, s_opts);

            t_stg2_start = tic;
            try
                sc = opti.solve(); status = 'Success';
                X_r = obj.Sx .* sc.value(Xhat); U_r = obj.Su .* sc.value(Uhat);
                T_r = sc.value(T_total); dt_r = sc.value(dt_row); so = sc;
                stats = sc.stats(); stats.t_wall_total = toc(t_stg2_start);
                if obj.RecordProgression
                    obj.ProgressionLogger.recordCollocationFrame( ...
                        sprintf('Stage 2: Optimal Solution Converged (Iter %d)', stats.iter_count), ...
                        X_r, U_r, T_r, stats.iter_count);
                end
            catch ME
                status = 'Infeasible';
                X_r = obj.Sx .* opti.debug.value(Xhat); U_r = obj.Su .* opti.debug.value(Uhat);
                T_r = opti.debug.value(T_total); dt_r = opti.debug.value(dt_row); so = opti;
                it_cnt = obj.MaxIter; st = opti.stats();
                if isfield(st, 'iter_count'), it_cnt = st.iter_count; end
                stats = struct('t_wall_total', toc(t_stg2_start), 'iter_count', it_cnt, ...
                    'return_status', 'Infeasible', 'error_msg', ME.message);
                if obj.RecordProgression
                    obj.ProgressionLogger.recordCollocationFrame( ...
                        sprintf('Stage 2: Non-Converged (Iter %d)', it_cnt), ...
                        X_r, U_r, T_r, it_cnt);
                end
            end

            stats.t_stage1_ms = obj.Stage1Diag.t_stage1_ms;
            stats.t_stage2_s  = stats.t_wall_total;
            stats.t_total_s   = toc(t_total_start);
            stats.stage1_diag = obj.Stage1Diag;

            dr_r = diff(X_r(5:7,:), 1, 2);
            t_vec = [0, cumsum(dt_r)];
            sol = struct('Status', status, 'T_total', T_r, 'L_path', sum(sqrt(sum(dr_r.^2, 1))), ...
                         'X', X_r, 'x', X_r, 'U', U_r, 'u', U_r, ...
                         'Time', t_vec, 't', t_vec, 'stats', stats, 'opti_obj', so, ...
                         'Progression', obj.Progression);
            obj.Solution = sol;
        end

        function [viol_to, viol_ld, tight_idx] = computeCorridorMetrics(obj, X)
            %% COMPUTECORRIDORMETRICS  Decoupled geometric evaluation for cones & funnels.
            import casadi.*
            N_n = size(X, 2);
            h_to_clear = obj.getTakeoffClearanceAlt();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                idx_to = find(obj.InitialGuess.X(7,:) - obj.r0(3) <= h_to_clear & (1:N_n) <= round(0.12*N_n));
                if isempty(idx_to), N_to = 2; else, N_to = max(2, max(idx_to)); end
            else
                N_to = max(2, round(0.06 * N_n));
            end
            dz_to = fmax(0, X(7, 1:N_to) - obj.r0(3));
            R_to = dz_to * tand(obj.GlideslopeAngle) + 0.15;
            rxy_to_sq = (X(5, 1:N_to) - obj.r0(1)).^2 + (X(6, 1:N_to) - obj.r0(2)).^2;
            viol_to = fmax(0, rxy_to_sq - R_to.^2);

            h_ld_clear = obj.getLandingFlareAlt();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                idx_ld = find(obj.InitialGuess.X(7,:) - obj.r_f(3) <= h_ld_clear & (1:N_n) >= round(0.65*N_n));
                if isempty(idx_ld), N_land = N_n - 2; else, N_land = min(N_n - 2, min(idx_ld)); end
            else
                N_land = min(N_n - 2, round(0.85 * N_n));
            end
            tight_idx = unique([1:N_to, N_land:N_n]);

            dz_ld = fmax(0, X(7, N_land:end) - obj.r_f(3));
            R_ld = dz_ld * tand(obj.GlideslopeAngle) + obj.FunnelCurvature * (dz_ld.^2) + 0.15;
            rxy_ld_sq = (X(5, N_land:end) - obj.r_f(1)).^2 + (X(6, N_land:end) - obj.r_f(2)).^2;
            viol_ld = fmax(0, rxy_ld_sq - R_ld.^2);
        end

        function applyManeuverConstraints(obj, opti, X, ~)
            %% APPLYMANEUVERCONSTRAINTS  Universal flight constraints and maneuver geometries.
            N_nodes = obj.N; N_n = size(X, 2);
            R33 = @(q) q(1,:).^2 - q(2,:).^2 - q(3,:).^2 + q(4,:).^2;

            % Liftoff alignment
            N_lo = max(2, round(0.04 * N_nodes));
            v_lo_max = max(0.60, 1.5 * norm(obj.r_f(1:2) - obj.r0(1:2)) / max(1.0, obj.T_initial));
            for ch = 8:9
                opti.subject_to(-v_lo_max <= X(ch, 1:N_lo));
                opti.subject_to(X(ch, 1:N_lo) <= v_lo_max);
                opti.subject_to(-0.35 <= X(ch, end-2:end));
                opti.subject_to(X(ch, end-2:end) <= 0.35);
            end
            opti.subject_to(R33(X(1:4, 1:N_lo)) >= cosd(28));
            opti.subject_to(X(10, end-2:end) >= -2.5);
            opti.subject_to(R33(X(1:4, end-2:end)) >= cosd(23));

            vd = obj.getVehicleDefaults();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                idx_ld = find(obj.InitialGuess.X(7,:) - obj.r_f(3) <= vd.h_ld_clear & (1:N_n) >= round(0.65*N_n));
                if isempty(idx_ld), N_land = N_n - 2; else, N_land = min(N_n - 2, min(idx_ld)); end
            else
                N_land = min(N_n - 2, round(0.85 * N_n));
            end
            opti.subject_to(X(10, N_land:end) >= -vd.v_desc_lim);
            opti.subject_to(R33(X(1:4, N_land:end)) >= cosd(38));

            if obj.Maneuver ~= "Backflip"
                opti.subject_to(R33(X(1:4,:)) >= cosd(55));
            end

            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    Nos = max(3, round(mp.f_orbit_start * N_nodes));
                    Noe = min(N_nodes - 2, round(mp.f_orbit_end * N_nodes));
                    Norb = max(2, Noe - Nos);
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    R = mp.circle_radius; h = mp.circle_alt; tr = max(0.35, obj.CircleTightness);

                    opti.subject_to(X(10, 1:Nos) >= 0.0); opti.subject_to(X(10, Noe:end) <= 0.0);
                    rsq = (X(5,Nos:Noe)-cx).^2 + (X(6,Nos:Noe)-cy).^2;
                    opti.subject_to((R-tr)^2 <= rsq); opti.subject_to(rsq <= (R+tr)^2);
                    opti.subject_to(abs(X(7,Nos:Noe) - h) <= 0.35);

                    rx = X(5,Nos:Noe)-cx; ry = X(6,Nos:Noe)-cy;
                    opti.subject_to(rx(1:end-1) .* ry(2:end) - ry(1:end-1) .* rx(2:end) >= 0.0);

                    kq1 = Nos + round(0.25*Norb); kq2 = Nos + round(0.50*Norb); kq3 = Nos + round(0.75*Norb);
                    opti.subject_to(X(6, kq1) - cy >= 0);
                    opti.subject_to(X(5, kq2) - cx <= 0);
                    opti.subject_to(X(6, kq3) - cy <= 0);
                    opti.subject_to((X(5, Noe) - (cx + R))^2 + (X(6, Noe) - cy)^2 <= tr^2);

                case "Backflip"
                    mp = obj.ManeuverParams;
                    Na  = max(3, round(mp.flip_start_frac * N_nodes)) + 1;
                    Nf  = round(0.5*(mp.flip_start_frac + mp.flip_end_frac) * N_nodes) + 1;
                    Nap = min(N_nodes, round(mp.flip_end_frac * N_nodes)) + 1;
                    opti.subject_to(X(10, 1:Na) >= 0.0); opti.subject_to(X(10, Nap:end) <= 0.0);
                    opti.subject_to(X(7,:) <= mp.apex_alt + 3.5);
                    opti.subject_to(X(7,Nf) >= mp.apex_alt - 3.5);
                    opti.subject_to(mp.q_inverted' * X(1:4,Nf) >= cos(mp.theta_tol / 2));
                    opti.subject_to(X(12, Na:Nap) <= -0.05);
                    for ch = [6, 9, 11, 13]
                        opti.subject_to(-0.75 <= X(ch,:)); opti.subject_to(X(ch,:) <= 0.75);
                    end

                case "Hop"
                    mp = obj.ManeuverParams;
                    if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                        [~, Nap_h] = max(obj.InitialGuess.X(7,:));
                    else
                        Nap_h = round(0.50 * N_nodes);
                    end
                    Nap_h = max(round(0.20*N_nodes), min(round(0.80*N_nodes), Nap_h));

                    opti.subject_to(X(7,Nap_h) >= mp.apex_alt - 2.5);
                    opti.subject_to(X(7,:) <= mp.apex_alt + 3.5);
                    opti.subject_to(X(10, 1:round(0.15*N_nodes)) >= 0.0);
                    opti.subject_to(X(10, round(0.85*N_nodes):end) <= 0.0);

                case {"Waypoint", "Custom"}
                    opti.subject_to(-0.25 <= X(4,:)); opti.subject_to(X(4,:) <= 0.25);
                    for ch = 11:12
                        opti.subject_to(-2.5 <= X(ch,:)); opti.subject_to(X(ch,:) <= 2.5);
                    end
                    opti.subject_to(-1.0 <= X(13,:)); opti.subject_to(X(13,:) <= 1.0);

                    if ~isempty(obj.Waypoints) && size(obj.Waypoints, 2) >= 3
                        K = size(obj.Waypoints, 2); Ttot = sum(obj.QPTraj.T); tg = obj.QPTraj.TagTimes;
                        for m = 2:K-1
                            tag_m = sprintf('wp_%d', m);
                            if isfield(tg, tag_m), t_wp = tg.(tag_m);
                            else,                  t_wp = (m - 1) * Ttot / (K - 1); end
                            k_wp = max(2, min(N_nodes, round(N_nodes * t_wp / Ttot) + 1));
                            tol_m = 0.80;
                            if ~isempty(obj.WaypointTolerances) && numel(obj.WaypointTolerances) >= m
                                tol_m = max(0.80, obj.WaypointTolerances(m));
                            end
                            opti.subject_to(sum((X(5:7, k_wp) - obj.Waypoints(:, m)).^2) <= tol_m^2);
                        end
                    end
            end
        end

        function applyCostFunction(obj, opti, X, U, T_total, viol_to, viol_ld, tight_idx)
            %% APPLYCOSTFUNCTION  Weighted projection onto Stage-1 spline + chatter suppression.
            import casadi.*
            N_nodes = obj.N; G = obj.InitialGuess; Xg = G.X; Tg = G.T_total;

            sp = obj.SigPos * ones(1, N_nodes + 1);
            if ~isempty(tight_idx), sp(tight_idx) = obj.SigPos * obj.TightFactor; end
            J_pos = sumsqr((X(5:7,:) - Xg(5:7,:)) ./ repmat(sp, 3, 1)) / (N_nodes + 1);
            J_vel = sumsqr((X(8:10,:) - Xg(8:10,:)) / obj.SigVel) / (N_nodes + 1);

            q_dots = sum1(X(1:4,:) .* Xg(1:4,:));
            J_att  = sum2(4 * (1 - q_dots.^2)) / (N_nodes + 1) / obj.SigAtt^2;
            J_T    = ((T_total - Tg) / (obj.SigT * Tg))^2;

            rmax = [obj.MaxGimbalRate; obj.MaxGimbalRate; obj.MaxThrustRate; obj.MaxRollRate];
            dt   = T_total / N_nodes;
            dU   = U(:,2:end) - U(:,1:end-1);
            d2U  = U(:,3:end) - 2*U(:,2:end-1) + U(:,1:end-2);
            J_du  = sumsqr(dU  ./ (repmat(rmax, 1, N_nodes-1) * dt)) / (N_nodes - 1);
            J_ddu = sumsqr(d2U ./ (repmat(rmax, 1, N_nodes-2) * dt)) / max(1, N_nodes - 2);

            v_ref = 0.15^2;
            J_corr = sumsqr(viol_to / v_ref) / size(viol_to, 2) + ...
                     sumsqr(viol_ld / v_ref) / size(viol_ld, 2);

            J = obj.WeightTrackPos  * J_pos  + ...
                obj.WeightTrackVel  * J_vel  + ...
                obj.WeightTrackAtt  * J_att  + ...
                obj.WeightTrackTime * J_T    + ...
                obj.WeightCorridor  * J_corr + ...
                obj.WeightSlew      * J_du   + ...
                obj.WeightCurv      * J_ddu;
            opti.minimize(J);
        end

        %% Dynamics & Kinematics
        function dyn_fnc = getCasADiDynamics(obj)
            import casadi.*
            q = MX.sym('q', 4); r = MX.sym('r', 3); v = MX.sym('v', 3);
            omegaB = MX.sym('omegaB', 3); m_lox = MX.sym('m_lox'); m_ipa = MX.sym('m_ipa');
            m_dry = MX.sym('m_dry'); gSym = MX.sym('g'); rTB = MX.sym('rTB');
            Ox_Z = MX.sym('Ox_Z'); OxMassI = MX.sym('OxMassI'); OxHeight = MX.sym('OxHeight');
            Fu_Z = MX.sym('Fu_Z'); FuMassI = MX.sym('FuMassI'); FuHeight = MX.sym('FuHeight');
            J = MX.sym('J', 3, 3); OxRadius = MX.sym('OxR'); FuRadius = MX.sym('FuR');
            MaxThrust = MX.sym('MT'); OF = MX.sym('OF'); MaxMdot = MX.sym('Mdot');
            theta = MX.sym('th'); phi = MX.sym('ph'); thrust = MX.sym('T'); roll = MX.sym('rl');

            m = m_dry + m_lox + m_ipa;
            TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            FI = quatRot(q).' * TB + [0; 0; -m * gSym];

            if obj.Vehicle == 0
                mdl = 0; mdi = 0; OFH = 0; FFH = 0;
            else
                mdl = -thrust/MaxThrust * OF/(1+OF) * MaxMdot;
                mdi = -thrust/MaxThrust * 1/(1+OF)  * MaxMdot;
                OFH = (m_lox / OxMassI) * OxHeight * 0.9;
                FFH = (m_ipa / FuMassI) * FuHeight * 0.9;
            end

            Jlox = diag([1/12*m_lox*(3*OxRadius^2+OFH^2), 1/12*m_lox*(3*OxRadius^2+OFH^2), 1/2*m_lox*OxRadius^2]);
            Jipa = diag([1/12*m_ipa*(3*FuRadius^2+FFH^2), 1/12*m_ipa*(3*FuRadius^2+FFH^2), 1/2*m_ipa*FuRadius^2]);
            OFL = Ox_Z + OFH/2; FFL = Fu_Z + FFH/2;
            CGz = (m_dry*rTB + m_lox*OFL + m_ipa*FFL) / m;
            dd = rTB - CGz; dl = OFL - CGz; di = FFL - CGz;
            J_tot = J + m_dry*diag([dd^2,dd^2,0]) + Jlox + m_lox*diag([dl^2,dl^2,0]) + Jipa + m_ipa*diag([di^2,di^2,0]);

            tDir = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            if obj.Vehicle == 0
                MB = zetaCross([0; 0; -CGz]) * TB + roll * tDir;
            else
                MB = zetaCross([0; 0; -CGz]) * TB + [0; 0; roll];
            end

            qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
            wdot = (MB - zetaCross(omegaB) * J_tot * omegaB) ./ diag(J_tot);

            params = [m_dry; gSym; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight; ...
                      J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot];
            dyn_fnc = Function('dyn_fnc', {[q; r; v; omegaB; m_lox; m_ipa], [theta; phi; thrust; roll], params}, ...
                               {[qdot; v; FI/m; wdot; mdl; mdi]});
        end

        function name = getVehicleName(obj)
            if obj.Vehicle == 1, name = "TOAD"; else, name = "ASTRA"; end
        end

        %% Public API & Delegations
        function figHandles = plot(obj, varargin)
            %% PLOT  Delegate dashboard plotting to PlotTrajectoryHybrid.
            figHandles = PlotTrajectoryHybrid(obj, varargin{:});
        end

        function ensureProgressionLogger(obj)
            if isempty(obj.ProgressionLogger)
                obj.ProgressionLogger = TrajectoryProgressionLogger(obj);
            end
        end

        function animOut = animateBuildProcess(obj, varargin)
            %% ANIMATEBUILDPROCESS  Forward progression animation to ProgressionLogger.
            obj.ensureProgressionLogger();
            animOut = obj.ProgressionLogger.animate(varargin{:});
        end

        function [tbl, filepath] = exportCSV(obj, filepath)
            %% EXPORTCSV  Export trajectory table matching TOAD simulation schema.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                error('TrajectoryOptimizerHybrid:NoSolution', 'Call solve() first.');
            end
            if nargin < 2 || isempty(filepath)
                tdir = fullfile(pwd, 'Guidance', 'Trajectories');
                if obj.SaveDir ~= "", tdir = obj.SaveDir; end
                if ~exist(tdir, 'dir'), mkdir(tdir); end
                fn = sprintf('%s_%s_v%03d', obj.getVehicleName(), obj.Maneuver, obj.Version);
                if obj.Filename ~= "", fn = obj.Filename; end
                filepath = fullfile(tdir, fn + ".csv");
            end
            ts = obj.Solution.Time(:); Xr = obj.Solution.X; Ur = obj.Solution.U;
            tbl = table(ts, Xr(1,:)', Xr(2,:)', Xr(3,:)', Xr(4,:)', ...
                Xr(5,:)', Xr(6,:)', Xr(7,:)', Xr(8,:)', Xr(9,:)', Xr(10,:)', ...
                Xr(11,:)', Xr(12,:)', Xr(13,:)', Xr(14,:)', Xr(15,:)', ...
                [Ur(1,:), Ur(1,end)]', [Ur(2,:), Ur(2,end)]', ...
                [Ur(3,:), Ur(3,end)]', [Ur(4,:), Ur(4,end)]', ...
                'VariableNames', {'Time', 'QuatW', 'QuatX', 'QuatY', 'QuatZ', ...
                                  'PosX', 'PosY', 'PosZ', 'VelX', 'VelY', 'VelZ', ...
                                  'AngRateX', 'AngRateY', 'AngRateZ', 'MassLox', 'MassFuel', ...
                                  'GimbalTheta', 'GimbalPhi', 'ThrustMag', 'RollCmd'});
            writetable(tbl, filepath);
            fprintf('  [TrajectoryOptimizerHybrid] Exported trajectory to: %s\n', filepath);
        end
    end
end
