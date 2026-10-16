classdef TrajectoryOptimizerHybrid < handle
    %% TrajectoryOptimizerHybrid - Distilled Two-Stage 6-DoF Optimal Trajectory Generator.
    % High-performance, mathematically rigorous hybrid guidance engine combining:
    %   Stage 1: Minimum-snap polynomial trajectory optimization (PolyTrajectoryQP)
    %            with corridor enforcement & differential-flatness state inversion.
    %   Stage 2: 6-DoF direct collocation refinement (IPOPT / CasADi RK4) with
    %            actuator chatter suppression and takeoff/landing flight corridors.
    %
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    properties
        % System & Vehicle
        constants                       % Vehicle constants struct from LoadTOADParams
        Vehicle double = 0              % 0 = ASTRA (electric), 1 = TOAD (liquid biprop)

        % Discretization & Timing
        N double = 60                   % Number of control intervals
        T_bounds = [8, 35]              % Total duration bounds [s]
        T_initial double = 16           % Initial duration guess [s]
        FixDuration logical = false     % When false, T is free within T_bounds guided by smooth J_T penalty

        % Flight Corridor Parameters
        CircleTightness double = 0.5    % Radial tolerance band [m]
        GlideslopeAngle double = 3.0    % Takeoff and landing cone half-angle [deg]
        FunnelCurvature double = 0.01  % Quadratic flaring parameter
        TakeoffClearanceAlt double = 0  % Liftoff clearance altitude [m] (0 = auto: 1.0m ASTRA, 2.5m TOAD)
        LandingFlareAlt double = 0      % Terminal descent flare altitude [m] (0 = auto: 1.2m ASTRA, 3.0m TOAD)

        % Stage-2 Projection Tolerances (1-sigma bounds)
        sig_pos  double = 0.30          % Position molding tolerance [m]
        sig_vel  double = 0.40          % Velocity molding tolerance [m/s]
        sig_att  double = deg2rad(5)    % Attitude molding tolerance [rad]
        sig_T    double = 0.03          % Duration tolerance [fraction of Stage-1 T] (fallback mode)
        tight_factor double = 0.25      % sig_pos multiplier at takeoff/landing corridor nodes

        % Cost Weights (Stage 2 = Projection onto Stage 1 + Actuator Smoothing)
        w_track_pos  double = 1.00      % Stage-1 position tracking
        w_track_vel  double = 0.50      % Stage-1 velocity tracking
        w_track_att  double = 0.25      % Stage-1 attitude tracking (invariant to q == -q)
        w_track_time double = 0.01      % Duration deviation (only active if T is unpinned)
        w_corridor   double = 1.00      % Soft takeoff cone / landing funnel violation penalty
        w_slew       double = 0.2      % Actuator first difference / (rate limit * dt)
        w_curv       double = 0.2      % Actuator second difference / (rate limit * dt)

        % Maneuver Definition & Boundary Conditions
        Maneuver string = "Circle"      % "Circle", "Backflip", "Hop", "Waypoint", "Custom"
        ManeuverParams struct = struct()
        Waypoints double = []           % 3 x K matrix of waypoints [m]
        WaypointTolerances double = []  % 1 x K tolerance sphere radii [m]
        q0 = [1; 0; 0; 0]
        r0 = [0; 0; 0]
        v0 = [0; 0; 0]
        w0 = [0; 0; 0]
        r_f = [0; 0; 0]
        v_f_tol = 0.35                  % Max allowable touchdown speed [m/s]

        % Actuator Physical Limits & Margins
        thrust_margin = 0.05
        gimbal_margin = 0.15
        max_gimbal_rate = deg2rad(35)   % Max gimbal slew rate [rad/s] (slightly relaxed to 35 deg/s)
        max_thrust_rate = 1500          % Max thrust rate [N/s] (slightly relaxed to 1500 N/s)
        max_roll_rate = 5               % Max roll torque rate [N*m/s] (slightly relaxed to 5 N*m/s)
        max_gimbal_angle = pi/15        % Max gimbal angle [rad] (12 deg)

        % Scaling & Solver Settings
        Sx double
        Su double
        L_c double = 50                 % Characteristic position length [m]
        MaxIter double = 120
        MaxCpuTime double = 45.0
        Tol double = 1.5e-2
        ConstrViolTol double = 5e-3
        PrintLevel double = 0

        % File Export & Diagnostics
        SaveDir string = ""
        Version double = 1
        Filename string = ""
        Stage1Diag struct = struct()
        QPTraj                          % PolyTrajectoryQP object
        InitialGuess struct = struct()
        Solution struct = struct()
        OptiVars struct = struct()

        % Progression Recording & Visualization
        RecordProgression logical = false
        ProgressionLogger               % TrajectoryProgressionLogger instance
    end

    properties (Dependent)
        Progression struct
    end

    methods
        function val = get.Progression(obj)
            if ~isempty(obj.ProgressionLogger)
                val = obj.ProgressionLogger.Progression;
            else
                val = struct();
            end
        end

        %% ======================== OBJECT LIFECYCLE & CONFIG ========================

        function obj = TrajectoryOptimizerHybrid(constants6DoF, varargin)
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizerHybrid:MissingConstants', ...
                    'TrajectoryOptimizerHybrid requires constants6DoF struct from LoadTOADParams.');
            end
            obj.constants = constants6DoF;

            % Ensure CasADi toolbox is on the search path
            if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir') && isempty(which('casadi.Opti'))
                addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
            end

            p = inputParser;
            p.KeepUnmatched = true;
            props = {'Vehicle','Maneuver','N','CircleTightness','GlideslopeAngle','FunnelCurvature',...
                     'TakeoffClearanceAlt','LandingFlareAlt','sig_pos','sig_vel','sig_att','sig_T',...
                     'tight_factor','FixDuration','w_track_pos','w_track_vel','w_track_att',...
                     'w_track_time','w_corridor','w_slew','w_curv',...
                     'PrintLevel','MaxIter','Tol','ConstrViolTol',...
                     'SaveDir','Version','Filename','RecordProgression'};
            for i = 1:numel(props)
                addParameter(p, props{i}, obj.(props{i}));
            end
            addParameter(p, 'T_initial', []);
            addParameter(p, 'T_bounds', []);
            parse(p, varargin{:});
            for i = 1:numel(props)
                obj.(props{i}) = p.Results.(props{i});
            end
            obj.Maneuver = string(obj.Maneuver);

            obj.updateScaling();
            obj.setManeuver(obj.Maneuver);

            if ~isempty(p.Results.T_bounds), obj.T_bounds = p.Results.T_bounds; end
            if ~isempty(p.Results.T_initial)
                obj.T_initial = p.Results.T_initial;
                if isempty(p.Results.T_bounds)
                    obj.T_bounds = [max(4.0, round(0.50 * obj.T_initial, 1)), ...
                                    min(60.0, round(2.0 * obj.T_initial, 1))];
                end
            end
        end

        function d = getVehicleDefaults(obj)
            %% GETVEHICLEDEFAULTS  Centralized vehicle-specific physical heuristics.
            d = struct();
            if obj.Vehicle == 0
                d.MinThrottle = 0.25;
                d.v_desc_lim  = 3.5;
                d.h_to_clear  = 1.0;
                d.h_ld_clear  = 1.2;
                d.dR = 5.0;  d.dA = 7.0;  d.dHA = 15.0; d.dBA = 20.0;
                d.accel_est   = 3.5;
                d.t_base_pad  = 3.0;
            else
                d.MinThrottle = 0.50;  % Unchanged per user directive
                d.v_desc_lim  = 6.0;
                d.h_to_clear  = 2.5;
                d.h_ld_clear  = 3.0;
                d.dR = 15.0; d.dA = 25.0; d.dHA = 50.0; d.dBA = 45.0;
                d.accel_est   = 2.2;   % Adjusted for 55% throttle floor descent dynamics
                d.t_base_pad  = 4.0;
            end
        end

        function h = getTakeoffClearanceAlt(obj)
            if obj.TakeoffClearanceAlt > 0
                h = obj.TakeoffClearanceAlt;
            else
                d = obj.getVehicleDefaults();
                h = d.h_to_clear;
            end
        end

        function h = getLandingFlareAlt(obj)
            if obj.LandingFlareAlt > 0
                h = obj.LandingFlareAlt;
            else
                d = obj.getVehicleDefaults();
                h = d.h_ld_clear;
            end
        end

        function updateScaling(obj)
            %% UPDATESCALING  O(1) state and control normalization vectors.
            V_c = 15; Omega_c = 2.0;
            if obj.Vehicle == 0
                obj.Sx = [1; 1; 1; 1; obj.L_c; obj.L_c; obj.L_c; ...
                          V_c; V_c; V_c; Omega_c; Omega_c; Omega_c; 1; 1];
            else
                obj.Sx = [1; 1; 1; 1; obj.L_c; obj.L_c; obj.L_c; ...
                          V_c; V_c; V_c; Omega_c; Omega_c; Omega_c; ...
                          obj.constants.OxMass; obj.constants.FuMass];
            end
            obj.Su = [obj.max_gimbal_angle; obj.max_gimbal_angle; ...
                      obj.constants.MaxThrust; obj.max_roll_rate];
        end

        function setBoundaries(obj, r0_in, r_f_in, varargin)
            %% SETBOUNDARIES  Configure launch pad and landing target coordinates.
            obj.r0 = r0_in(:);
            obj.r_f = r_f_in(:);
            p = inputParser;
            addParameter(p, 'v0', [0; 0; 0]);
            addParameter(p, 'q0', [1; 0; 0; 0]);
            addParameter(p, 'v_f_tol', obj.v_f_tol);
            parse(p, varargin{:});
            obj.v0 = p.Results.v0(:);
            obj.q0 = p.Results.q0(:);
            obj.v_f_tol = double(p.Results.v_f_tol);

            obj.updateDurationEstimate();
            obj.QPTraj = [];
            obj.InitialGuess = struct();
        end

        function updateDurationEstimate(obj)
            %% UPDATEDURATIONESTIMATE  Centralized kinematics-based baseline duration and bounds.
            if obj.Maneuver == "Hop" && isstruct(obj.ManeuverParams) && isfield(obj.ManeuverParams, 'apex_alt')
                mp = obj.ManeuverParams;
                dh = max(1.5, mp.apex_alt - max(obj.r0(3), obj.r_f(3)));
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
            p = inputParser;
            vd = obj.getVehicleDefaults();

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
                    if obj.Vehicle == 0
                        obj.T_bounds = [16, 45]; obj.T_initial = 28;
                    else
                        obj.T_bounds = [20, 65]; obj.T_initial = 42;
                    end

                case "Backflip"
                    addParameter(p, 'apex_alt', vd.dBA);
                    addParameter(p, 'flip_start_frac', 0.30);
                    addParameter(p, 'flip_end_frac', 0.70);
                    addParameter(p, 'theta_tol', deg2rad(25));
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt), ...
                        'flip_start_frac', double(p.Results.flip_start_frac), ...
                        'flip_end_frac', double(p.Results.flip_end_frac), ...
                        'theta_tol', double(p.Results.theta_tol), ...
                        'q_inverted', [0; 0; -1; 0]);
                    if obj.Vehicle == 0
                        obj.T_bounds = [10, 24]; obj.T_initial = 16;
                    else
                        obj.T_bounds = [16, 36]; obj.T_initial = 24;
                    end

                case "Hop"
                    addParameter(p, 'apex_alt', vd.dHA);
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt));
                    obj.updateDurationEstimate();

                case {"Waypoint", "Custom"}
                    addParameter(p, 'Waypoints', []);
                    addParameter(p, 'Tolerances', []);
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
            obj.QPTraj = [];
            obj.InitialGuess = struct();
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
            if isscalar(obj.WaypointTolerances)
                obj.WaypointTolerances = repmat(obj.WaypointTolerances, 1, K);
            end
            obj.T_initial = double(p.Results.T_total);
            obj.T_bounds = [max(4.0, round(0.50 * obj.T_initial, 1)), ...
                            min(70.0, round(2.0 * obj.T_initial, 1))];
            obj.Maneuver = "Waypoint";
            obj.ManeuverParams = struct();
            obj.r0 = obj.Waypoints(:, 1);
            obj.r_f = obj.Waypoints(:, end);
            obj.QPTraj = [];
            obj.InitialGuess = struct();
        end

        %% ======================== STAGE 1: POLYNOMIAL QP ENGINE ========================

        function [Fix, T_seg, Tag] = buildKeyframes(obj)
            %% BUILDKEYFRAMES  Convert maneuver geometry to PolyTrajectoryQP keyframe arrays.
            R = 4;               % Minimum-snap (7th-order septic polynomials)
            vL = [0; 0; -0.10];  % Gentle terminal landing velocity

            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    Rad = mp.circle_radius; h = mp.circle_alt;
                    N_orb = 16;
                    angles = linspace(0, 2*pi, N_orb + 1);
                    K = 1 + (N_orb + 1) + 1;
                    Fix = nan(3, R, K); Tag = repmat({''}, 1, K);

                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    for j = 1:(N_orb + 1)
                        Fix(:,1,1 + j) = [cx + Rad*cos(angles(j)); cy + Rad*sin(angles(j)); h];
                    end
                    Tag{2} = 'orb_entry'; Tag{N_orb + 2} = 'orb_exit';
                    Fix(:,1,K) = obj.r_f; Fix(:,2,K) = vL; Tag{K} = 'land';

                    M = K - 1; Ts = obj.T_initial;
                    T_seg = zeros(1, M);
                    T_seg(1) = 0.25 * Ts;
                    for j = 1:N_orb, T_seg(1 + j) = 0.50 * Ts / N_orb; end
                    T_seg(M) = 0.25 * Ts;

                case {"Hop", "Backflip"}
                    % Shared parabolic waypoint topology
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
                    if td < 1e-3
                        T_seg = (obj.T_initial / M) * ones(1, M);
                    else
                        T_seg = max(1.5, obj.T_initial * dists / td);
                        T_seg = obj.T_initial * T_seg / sum(T_seg);
                    end

                otherwise
                    K = 2; Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0;
                    Fix(:,3,1) = [0; 0; 0]; Fix(:,4,1) = [0; 0; 0]; Tag{1} = 'launch';
                    Fix(:,1,2) = obj.r_f; Fix(:,2,2) = vL;
                    Fix(:,3,2) = [0; 0; 0]; Fix(:,4,2) = [0; 0; 0]; Tag{2} = 'land';
                    T_seg = obj.T_initial;
            end
        end

        function buildQPTrajectory(obj)
            %% BUILDQPTRAJECTORY  Construct PolyTrajectoryQP, optimize times, enforce corridors.
            t0 = tic;
            [Fix, T_seg, Tag] = obj.buildKeyframes();
            obj.QPTraj = PolyTrajectoryQP(Fix, T_seg, 'Tag', Tag);

            if obj.RecordProgression
                obj.ensureProgressionLogger();
                obj.ProgressionLogger.recordQP('Stage 1A: QP Baseline Spline', obj.QPTraj, Fix);
            end

            % Optimize segment durations for multi-segment trajectories (Richter Sec. V)
            if obj.Maneuver ~= "Circle" && numel(obj.QPTraj.T) >= 2
                obj.QPTraj.optimizeTimes(10);
                if obj.RecordProgression
                    obj.ProgressionLogger.recordQP('Stage 1B: Time-Optimized Spline', obj.QPTraj, Fix);
                end
            end

            % Enforce geometric flight corridors via iterative keyframe insertion
            corridors = { ...
                struct('fn', @(ts,rs,tg) obj.evalTakeoffCone(ts,rs,tg)), ...
                struct('fn', @(ts,rs,tg) obj.evalLandingFunnel(ts,rs,tg)) ...
            };
            obj.QPTraj.enforceCorridors(corridors, 20, 0.05);
            if obj.RecordProgression
                obj.ProgressionLogger.recordQP('Stage 1C: Corridor-Enforced Spline', obj.QPTraj, Fix);
            end

            % Synchronize keyframe tags to ManeuverParams
            tg = obj.QPTraj.tagTimes();
            Ttot = sum(obj.QPTraj.T);
            if isfield(tg, 'orb_entry'), obj.ManeuverParams.f_orbit_start = tg.orb_entry / Ttot; end
            if isfield(tg, 'orb_exit'),  obj.ManeuverParams.f_orbit_end   = tg.orb_exit / Ttot;  end
            if isfield(tg, 'apex'),      obj.ManeuverParams.f_apex        = tg.apex / Ttot;      end

            obj.Stage1Diag = struct('t_qp_ms', toc(t0) * 1000);
        end

        function [viol, tgt] = evalTakeoffCone(obj, ts, rs, tg)
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts);
            t_end = 0.15 * t_max;
            if isfield(tg, 'orb_entry')
                t_end = min(0.20 * t_max, tg.orb_entry);
            elseif isfield(tg, 'wp_2')
                t_end = min(0.20 * t_max, tg.wp_2);
            end

            r0p = obj.r0; gsa = tand(obj.GlideslopeAngle);
            h_liftoff = obj.getTakeoffClearanceAlt();

            for i = 1:n
                if ts(i) > t_end, continue; end
                dz = rs(3,i) - r0p(3);
                if dz < 0.01 || dz > h_liftoff, continue; end
                rxy = sqrt((rs(1,i)-r0p(1))^2 + (rs(2,i)-r0p(2))^2 + 1e-4);
                lim = dz * gsa + 0.15;
                v = rxy - lim;
                if v > 0
                    viol(i) = v;
                    sc = lim / max(1e-6, rxy);
                    tgt(1,i) = r0p(1) + (rs(1,i)-r0p(1)) * sc;
                    tgt(2,i) = r0p(2) + (rs(2,i)-r0p(2)) * sc;
                end
            end
        end

        function [viol, tgt] = evalLandingFunnel(obj, ts, rs, tg)
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts);
            t_start = 0.70 * t_max;
            if isfield(tg, 'orb_exit')
                t_start = max(0.65 * t_max, tg.orb_exit);
            elseif isfield(tg, 'apex')
                t_start = max(0.60 * t_max, tg.apex);
            elseif ~isempty(obj.Waypoints) && size(obj.Waypoints, 2) >= 2
                K = size(obj.Waypoints, 2);
                tag_penult = sprintf('wp_%d', K-1);
                if isfield(tg, tag_penult), t_start = tg.(tag_penult); end
            end

            gsa = tand(obj.GlideslopeAngle);
            cf = obj.FunnelCurvature;
            rfp = obj.r_f;
            h_flare = obj.getLandingFlareAlt();

            for i = 1:n
                if ts(i) < t_start, continue; end
                dz = rs(3,i) - rfp(3);
                if dz < 0.01 || dz > h_flare, continue; end
                rxy = sqrt((rs(1,i)-rfp(1))^2 + (rs(2,i)-rfp(2))^2 + 1e-4);
                lim = dz * gsa + cf * dz^2 + 0.15;
                v = rxy - lim;
                if v > 0
                    viol(i) = v;
                    sc = lim / max(1e-6, rxy);
                    tgt(1,i) = rfp(1) + (rs(1,i)-rfp(1)) * sc;
                    tgt(2,i) = rfp(2) + (rs(2,i)-rfp(2)) * sc;
                end
            end
        end

        %% ======================== STAGE 1: FLATNESS INVERSION ========================

        function guess = buildInitialGuess(obj)
            %% BUILDINITIALGUESS  Sample QP trajectory and invert via differential flatness.
            if isempty(obj.QPTraj), obj.buildQPTrajectory(); end
            t0 = tic;
            Nn = obj.N + 1; Nc = obj.N;
            g = obj.constants.g; MT = obj.constants.MaxThrust;
            thr_allow = (1 - obj.thrust_margin) * MT;

            % Dynamic thrust sizing: scale polynomial time if thrust exceeds vehicle capability
            if obj.Vehicle == 0
                m_chk = obj.constants.m_dry;
            else
                m_chk = obj.constants.m_dry + obj.constants.OxMass + obj.constants.FuMass;
            end
            Tq = sum(obj.QPTraj.T);
            tv_c = linspace(0, Tq, 40);
            a_c = obj.QPTraj.evalAt(tv_c, 2);
            Fmag = m_chk * sqrt(sum((a_c + [0; 0; g]).^2, 1));
            kap = sqrt(max(Fmag) / thr_allow);
            if kap > 1.0
                obj.QPTraj.scaleTime(kap);
            end

            % Sample QP trajectory at NLP grid points
            T_total = sum(obj.QPTraj.T);
            tv = linspace(0, T_total, Nn);
            dt = T_total / Nc; dt_row = dt * ones(1, Nc);
            r_s = obj.QPTraj.evalAt(tv, 0);
            v_s = obj.QPTraj.evalAt(tv, 1);
            a_s = obj.QPTraj.evalAt(tv, 2);
            j_s = obj.QPTraj.evalAt(tv, 3);

            % Mass depletion profile
            m_dry = obj.constants.m_dry;
            ml = zeros(1, Nn); mi = zeros(1, Nn);
            if obj.Vehicle == 1
                ml(1) = obj.constants.OxMass;
                mi(1) = obj.constants.FuMass;
            end
            m_vec = (m_dry + ml(1) + mi(1)) * ones(1, Nn);

            q_s = zeros(4, Nn); w_s = zeros(3, Nn); U_s = zeros(4, Nc);
            Jyy = obj.constants.J(2,2);
            cgz = obj.constants.rTB;
            plim = sin((1 - obj.gimbal_margin) * obj.max_gimbal_angle);
            OF = obj.constants.OF;
            Mdot = obj.constants.MaxMdot;
            t_lo = (0.25 + obj.thrust_margin) * MT;
            t_hi = (1.00 - obj.thrust_margin) * MT;

            if obj.Maneuver == "Backflip"
                % Backflip cycloidal attitude & pitch rate ramp
                mp = obj.ManeuverParams;
                ks = max(3, round(mp.flip_start_frac * Nc)) + 1;
                ke = min(Nc - 2, round(mp.flip_end_frac * Nc)) + 1;
                Lf = max(2, ke - ks); Tf = tv(ke) - tv(ks);

                for k = 1:Nn
                    if k <= ks
                        q_s(:,k) = obj.q0; w_s(:,k) = [0; 0; 0];
                    elseif k <= ke
                        tau = (k - ks) / Lf;
                        th  = -2*pi * (10*tau^3 - 15*tau^4 + 6*tau^5);
                        dth = -(2*pi/Tf) * (30*tau^2 - 60*tau^3 + 30*tau^4);
                        q_s(:,k) = [cos(th/2); 0; sin(th/2); 0];
                        w_s(:,k) = [0; dth; 0];
                    else
                        q_s(:,k) = -obj.q0; w_s(:,k) = [0; 0; 0];
                    end
                end

                for k = 1:Nc
                    F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                    U_s(3,k) = max(t_lo, min(t_hi, norm(F_req)));
                    if k > ks && k <= ke
                        tau = (k - ks) / Lf;
                        th  = -2*pi * (10*tau^3 - 15*tau^4 + 6*tau^5);
                        ddth = -(2*pi/Tf^2) * (60*tau - 180*tau^2 + 120*tau^3);
                        cos_th = cos(th);
                        if cos_th < 0
                            throttle_mult = 1.0 + 0.20 * cos_th;
                            U_s(3,k) = max(t_lo, U_s(3,k) * throttle_mult);
                        end
                        sp = -(Jyy * ddth) / max(1.0, cgz * U_s(3,k));
                        U_s(2,k) = asin(max(-plim, min(plim, sp)));
                    end
                    if obj.Vehicle == 1
                        mdl_k = -(U_s(3,k)/MT) * (OF/(1+OF)) * Mdot;
                        mdi_k = -(U_s(3,k)/MT) * (1/(1+OF))  * Mdot;
                        ml(k+1) = max(0.02*obj.constants.OxMass, ml(k) + dt*mdl_k);
                        mi(k+1) = max(0.02*obj.constants.FuMass, mi(k) + dt*mdi_k);
                        m_vec(k+1) = m_dry + ml(k+1) + mi(k+1);
                    end
                end
            else
                % Flatness inversion: thrust vector -> quaternion
                for k = 1:Nc
                    F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                    U_s(3,k) = max(t_lo, min(t_hi, norm(F_req)));
                    if obj.Vehicle == 1
                        mdl_k = -(U_s(3,k)/MT) * (OF/(1+OF)) * Mdot;
                        mdi_k = -(U_s(3,k)/MT) * (1/(1+OF))  * Mdot;
                        ml(k+1) = max(0.02*obj.constants.OxMass, ml(k) + dt*mdl_k);
                        mi(k+1) = max(0.02*obj.constants.FuMass, mi(k) + dt*mdi_k);
                        m_vec(k+1) = m_dry + ml(k+1) + mi(k+1);
                    end
                end

                for k = 1:Nn
                    F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                    T_mag = norm(F_req);
                    b3 = F_req / max(1e-3, T_mag);
                    q_flat = TrajectoryOptimizerHybrid.vectorToQuat(b3);
                    if k <= 4
                        alpha = (k - 1) / 4;
                        if dot(obj.q0, q_flat) < 0, q_flat = -q_flat; end
                        qb = (1-alpha)*obj.q0 + alpha*q_flat;
                        q_s(:,k) = qb / norm(qb);
                    elseif k >= Nn - 3
                        alpha = (Nn - k) / 3;
                        if dot(obj.q0, q_flat) < 0, q_flat = -q_flat; end
                        qb = (1-alpha)*obj.q0 + alpha*q_flat;
                        q_s(:,k) = qb / norm(qb);
                    else
                        q_s(:,k) = q_flat;
                    end
                    b3_dot = (j_s(:,k) - b3*(b3'*j_s(:,k))) * m_vec(k) / max(1e-3, T_mag);
                    R_IB = TrajectoryOptimizerHybrid.quatToRot(q_s(:,k));
                    b3d_B = R_IB' * b3_dot;
                    w_s(:,k) = [-b3d_B(2); b3d_B(1); 0];
                end
                q_s(:,1) = obj.q0; q_s(:,end) = obj.q0;
                w_s(:,1) = obj.w0; w_s(:,end) = [0; 0; 0];

                % TVC gimbal synthesis
                J_mat = obj.constants.J; rTB = obj.constants.rTB;
                for k = 1:Nc
                    w_dot = (w_s(:,k+1) - w_s(:,k)) / dt;
                    w_mid = 0.5 * (w_s(:,k) + w_s(:,k+1));
                    tau_req = J_mat * w_dot + cross(w_mid, J_mat * w_mid);
                    denom = max(1.0, U_s(3,k) * rTB);
                    U_s(1,k) = asin(max(-plim, min(plim, -tau_req(2)/denom)));
                    U_s(2,k) = asin(max(-plim, min(plim,  tau_req(1)/denom)));
                    U_s(4,k) = 0;
                end
            end

            % Vectorized forward-backward rate smoothing across all 4 control channels
            rate_limits = [obj.max_gimbal_rate; obj.max_gimbal_rate; obj.max_thrust_rate; obj.max_roll_rate];
            max_du = 0.90 * rate_limits * dt;
            for k = 2:Nc
                U_s(:,k) = max(U_s(:,k-1) - max_du, min(U_s(:,k-1) + max_du, U_s(:,k)));
            end
            for k = Nc-1:-1:1
                U_s(:,k) = max(U_s(:,k+1) - max_du, min(U_s(:,k+1) + max_du, U_s(:,k)));
            end

            % Fix boundary points
            r_s(:,1) = obj.r0; r_s(:,end) = obj.r_f; v_s(:,end) = [0; 0; -0.05];

            % Package initial guess
            Xp = [q_s; r_s; v_s; w_s; ml; mi];
            obj.T_initial = T_total;
            if isempty(obj.T_bounds)
                obj.T_bounds = [max(4.0, round(0.70 * T_total, 1)), ...
                                min(70.0, round(1.40 * T_total, 1))];
            end

            obj.Stage1Diag.T_total = T_total;
            obj.Stage1Diag.max_thrust = max(U_s(3,:));
            obj.Stage1Diag.max_gimbal_deg = rad2deg(max(sqrt(U_s(1,:).^2 + U_s(2,:).^2)));
            obj.Stage1Diag.t_stage1_ms = toc(t0) * 1000;

            guess = struct('Time', tv, 'X', Xp, 'U', U_s, 'T_total', T_total, ...
                           'dt_row', dt_row, 'Xhat', Xp ./ obj.Sx, 'Uhat', U_s ./ obj.Su);
            obj.InitialGuess = guess;

            if obj.RecordProgression
                obj.ensureProgressionLogger();
                obj.ProgressionLogger.recordGuess('Stage 1D: Flatness Inversion (Initial Guess)', guess);
            end
        end

        %% ======================== STAGE 2: DIRECT COLLOCATION ========================

        function sol = solve(obj)
            %% SOLVE  Two-Stage Pipeline: QP Polynomial Seed -> 6-DoF Direct Collocation
            import casadi.*
            t_total_start = tic;

            if isempty(obj.QPTraj), obj.buildQPTrajectory(); end
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
            if obj.FixDuration
                opti.subject_to(T_total == obj.InitialGuess.T_total);
            end
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

            % Universal dynamic domain bounds
            max_h_extent = max([norm(obj.r0(1:2)), norm(obj.r_f(1:2)), norm(obj.r_f(1:2) - obj.r0(1:2)), 30.0]);
            if ~isempty(obj.Waypoints)
                max_h_extent = max(max_h_extent, max(sqrt(sum(obj.Waypoints(1:2,:).^2, 1))));
            end
            R_horiz = max(150.0, 3.0 * max_h_extent);
            opti.subject_to(-R_horiz <= X(5,:)); opti.subject_to(X(5,:) <= R_horiz);
            opti.subject_to(-R_horiz <= X(6,:)); opti.subject_to(X(6,:) <= R_horiz);

            alts = [obj.r0(3), obj.r_f(3), 20.0];
            if isfield(obj.ManeuverParams, 'apex_alt'),   alts = [alts, obj.ManeuverParams.apex_alt];   end
            if isfield(obj.ManeuverParams, 'circle_alt'), alts = [alts, obj.ManeuverParams.circle_alt]; end
            if ~isempty(obj.Waypoints),                   alts = [alts, obj.Waypoints(3,:)];            end
            alt_ceil = max(180.0, 2.5 * max(alts));
            opti.subject_to(0 <= X(7,:)); opti.subject_to(X(7,:) <= alt_ceil);

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
            opti.subject_to(sum(X(8:10,end).^2) <= obj.v_f_tol^2);
            if obj.Vehicle == 1
                opti.subject_to(X(14,end) >= 0.02*m_lox0);
                opti.subject_to(X(15,end) >= 0.02*m_ipa0);
            end

            % Control bounds & slew rates
            MT = obj.constants.MaxThrust; tm = obj.thrust_margin;
            gm = obj.gimbal_margin; mg = obj.max_gimbal_angle;
            vd = obj.getVehicleDefaults();
            opti.subject_to((vd.MinThrottle + tm)*MT <= U(3,:));
            opti.subject_to(U(3,:) <= (1 - tm)*MT);
            opti.subject_to(-(1 - gm)*mg <= U(1,:)); opti.subject_to(U(1,:) <= (1 - gm)*mg);
            opti.subject_to(-(1 - gm)*mg <= U(2,:)); opti.subject_to(U(2,:) <= (1 - gm)*mg);
            opti.subject_to(-(1 - tm)*obj.max_roll_rate <= U(4,:));
            opti.subject_to(U(4,:) <= (1 - tm)*obj.max_roll_rate);

            % Control slew rate limits
            dU = U(:,2:end) - U(:,1:end-1); dts = dt_row(1:end-1);
            for ch = [1, 2]
                opti.subject_to(-obj.max_gimbal_rate * dts <= dU(ch,:));
                opti.subject_to(dU(ch,:) <= obj.max_gimbal_rate * dts);
            end
            opti.subject_to(-obj.max_thrust_rate * dts <= dU(3,:));
            opti.subject_to(dU(3,:) <= obj.max_thrust_rate * dts);
            opti.subject_to(-obj.max_roll_rate * dts <= dU(4,:));
            opti.subject_to(dU(4,:) <= obj.max_roll_rate * dts);

            % Decoupled corridor metrics
            [viol_to, viol_ld, tight_idx] = obj.computeCorridorMetrics(X);

            % Apply physics constraints & projection cost
            obj.applyManeuverConstraints(opti, X, tight_idx);
            obj.applyCostFunction(opti, X, U, T_total, viol_to, viol_ld, tight_idx);

            % Seed IPOPT with analytical flatness guess
            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);

            % Solver configuration
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, 'tol', obj.Tol, ...
                'constr_viol_tol', obj.ConstrViolTol, ...
                'dual_inf_tol', 1.0, ...
                'acceptable_tol', max(obj.Tol, 2e-2), ...
                'acceptable_constr_viol_tol', 1e-3, ...
                'acceptable_dual_inf_tol', 1e5, ...
                'acceptable_iter', 2, ...
                'max_cpu_time', obj.MaxCpuTime, ...
                'mu_strategy', 'adaptive', 'mu_init', 1e-3, ...
                'bound_push', 1e-6, 'bound_frac', 1e-6, ...
                'print_level', obj.PrintLevel);
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
                it_cnt = obj.MaxIter;
                try
                    st = opti.stats();
                    if isfield(st, 'iter_count'), it_cnt = st.iter_count; end
                catch
                end
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
            LP = sum(sqrt(sum(dr_r.^2, 1)));
            t_vec = [0, cumsum(dt_r)];
            sol = struct('Status', status, 'T_total', T_r, 'L_path', LP, ...
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
                z_g = obj.InitialGuess.X(7, :);
                idx_to = find(z_g - obj.r0(3) <= h_to_clear & (1:N_n) <= round(0.12*N_n));
                if isempty(idx_to), N_to = 2; else, N_to = max(2, max(idx_to)); end
            else
                N_to = max(2, round(0.06 * N_n));
            end
            r0_p = obj.r0;
            dz_to = fmax(0, X(7, 1:N_to) - r0_p(3));
            R_to = dz_to * tand(obj.GlideslopeAngle) + 0.15;
            rxy_to_sq = (X(5, 1:N_to) - r0_p(1)).^2 + (X(6, 1:N_to) - r0_p(2)).^2;
            viol_to = fmax(0, rxy_to_sq - R_to.^2);

            h_ld_clear = obj.getLandingFlareAlt();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                z_g = obj.InitialGuess.X(7, :);
                idx_ld = find(z_g - obj.r_f(3) <= h_ld_clear & (1:N_n) >= round(0.65*N_n));
                if isempty(idx_ld), N_land = N_n - 2; else, N_land = min(N_n - 2, min(idx_ld)); end
            else
                N_land = min(N_n - 2, round(0.85 * N_n));
            end

            tight_idx = unique([1:N_to, N_land:N_n]);

            rf_p = obj.r_f;
            dz_ld = fmax(0, X(7, N_land:end) - rf_p(3));
            R_ld = dz_ld * tand(obj.GlideslopeAngle) + obj.FunnelCurvature * (dz_ld.^2) + 0.15;
            rxy_ld_sq = (X(5, N_land:end) - rf_p(1)).^2 + (X(6, N_land:end) - rf_p(2)).^2;
            viol_ld = fmax(0, rxy_ld_sq - R_ld.^2);
        end

        function applyManeuverConstraints(obj, opti, X, ~)
            %% APPLYMANEUVERCONSTRAINTS  Universal flight constraints and maneuver geometries.
            N_nodes = obj.N;
            N_n = size(X, 2);
            R33 = @(q) q(1,:).^2 - q(2,:).^2 - q(3,:).^2 + q(4,:).^2;

            % Liftoff alignment
            N_lo = max(2, round(0.04 * N_nodes));
            v_lo_max = max(0.60, 1.5 * norm(obj.r_f(1:2) - obj.r0(1:2)) / max(1.0, obj.T_initial));
            for ch = 8:9
                opti.subject_to(-v_lo_max <= X(ch, 1:N_lo));
                opti.subject_to(X(ch, 1:N_lo) <= v_lo_max);
            end
            opti.subject_to(R33(X(1:4, 1:N_lo)) >= cosd(28));

            % Touchdown alignment (terminal 3 nodes)
            for ch = 8:9
                opti.subject_to(-0.35 <= X(ch, end-2:end));
                opti.subject_to(X(ch, end-2:end) <= 0.35);
            end
            opti.subject_to(X(10, end-2:end) >= -2.5);
            opti.subject_to(R33(X(1:4, end-2:end)) >= cosd(23));

            % Flared descent vertical velocity limit & funnel upright tilt
            vd = obj.getVehicleDefaults();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                z_g = obj.InitialGuess.X(7, :);
                idx_ld = find(z_g - obj.r_f(3) <= vd.h_ld_clear & (1:N_n) >= round(0.65*N_n));
                if isempty(idx_ld), N_land = N_n - 2; else, N_land = min(N_n - 2, min(idx_ld)); end
            else
                N_land = min(N_n - 2, round(0.85 * N_n));
            end
            opti.subject_to(X(10, N_land:end) >= -vd.v_desc_lim);
            opti.subject_to(R33(X(1:4, N_land:end)) >= cosd(38));

            % Universal upright tilt envelope for non-Backflip maneuvers (slightly relaxed to 55 deg)
            if obj.Maneuver ~= "Backflip"
                opti.subject_to(R33(X(1:4,:)) >= cosd(55));
            end

            % Maneuver-specific geometric constraints
            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    Nos = max(3, round(mp.f_orbit_start * N_nodes));
                    Noe = min(N_nodes - 2, round(mp.f_orbit_end * N_nodes));
                    Norb = max(2, Noe - Nos);
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    R = mp.circle_radius; h = mp.circle_alt;
                    tr = max(0.35, obj.CircleTightness);

                    opti.subject_to(X(10, 1:Nos) >= 0.0);
                    opti.subject_to(X(10, Noe:end) <= 0.0);

                    rsq = (X(5,Nos:Noe)-cx).^2 + (X(6,Nos:Noe)-cy).^2;
                    opti.subject_to((R-tr)^2 <= rsq); opti.subject_to(rsq <= (R+tr)^2);
                    opti.subject_to(abs(X(7,Nos:Noe) - h) <= 0.35);

                    rx = X(5,Nos:Noe)-cx; ry = X(6,Nos:Noe)-cy;
                    cprog = rx(1:end-1) .* ry(2:end) - ry(1:end-1) .* rx(2:end);
                    opti.subject_to(cprog >= 0.0);

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
                    att_tol = cos(mp.theta_tol / 2);

                    opti.subject_to(X(10, 1:Na) >= 0.0);
                    opti.subject_to(X(10, Nap:end) <= 0.0);
                    opti.subject_to(X(7,:) <= mp.apex_alt + 3.5);
                    opti.subject_to(X(7,Nf) >= mp.apex_alt - 3.5);
                    opti.subject_to(mp.q_inverted' * X(1:4,Nf) >= att_tol);
                    opti.subject_to(X(12, Na:Nap) <= -0.05);

                    % Planar symmetry constraints
                    opti.subject_to(-0.75 <= X(6,:));  opti.subject_to(X(6,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(9,:));  opti.subject_to(X(9,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(11,:)); opti.subject_to(X(11,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(13,:)); opti.subject_to(X(13,:) <= 0.75);

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
                        opti.subject_to(-2.5 <= X(ch,:));
                        opti.subject_to(X(ch,:) <= 2.5);
                    end
                    opti.subject_to(-1.0 <= X(13,:)); opti.subject_to(X(13,:) <= 1.0);

                    if ~isempty(obj.Waypoints)
                        K = size(obj.Waypoints, 2);
                        Ttot = sum(obj.QPTraj.T);
                        tg = obj.QPTraj.tagTimes();
                        for m = 2:(K-1)
                            tag_name = sprintf('wp_%d', m);
                            if isfield(tg, tag_name), t_wp = tg.(tag_name);
                            else,                     t_wp = (m - 1) * Ttot / (K - 1); end
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
            N_nodes = obj.N;
            G  = obj.InitialGuess;
            Xg = G.X;
            Tg = G.T_total;

            % 1. Position tracking (tightened at corridor nodes)
            sp = obj.sig_pos * ones(1, N_nodes + 1);
            if ~isempty(tight_idx)
                sp(tight_idx) = obj.sig_pos * obj.tight_factor;
            end
            J_pos = sumsqr((X(5:7,:) - Xg(5:7,:)) ./ repmat(sp, 3, 1)) / (N_nodes + 1);

            % 2. Velocity tracking
            J_vel = sumsqr((X(8:10,:) - Xg(8:10,:)) / obj.sig_vel) / (N_nodes + 1);

            % 3. Attitude tracking (invariant to quaternion double-cover)
            q_dots = sum1(X(1:4,:) .* Xg(1:4,:));
            J_att  = sum2(4 * (1 - q_dots.^2)) / (N_nodes + 1) / obj.sig_att^2;

            % 4. Duration penalty (only active when duration is unpinned)
            J_T = ((T_total - Tg) / (obj.sig_T * Tg))^2;

            % 5. Actuator smoothness (normalized first and second differences)
            rmax = [obj.max_gimbal_rate; obj.max_gimbal_rate; obj.max_thrust_rate; obj.max_roll_rate];
            dt   = T_total / N_nodes;
            dU   = U(:,2:end) - U(:,1:end-1);
            d2U  = U(:,3:end) - 2*U(:,2:end-1) + U(:,1:end-2);
            J_du  = sumsqr(dU  ./ (repmat(rmax, 1, N_nodes-1) * dt)) / (N_nodes - 1);
            J_ddu = sumsqr(d2U ./ (repmat(rmax, 1, N_nodes-2) * dt)) / max(1, N_nodes - 2);

            % 6. Soft flight corridor penalty
            v_ref = 0.15^2;
            J_corr = sumsqr(viol_to / v_ref) / size(viol_to, 2) + ...
                     sumsqr(viol_ld / v_ref) / size(viol_ld, 2);

            % Total Stage-2 Objective
            J = obj.w_track_pos  * J_pos  + ...
                obj.w_track_vel  * J_vel  + ...
                obj.w_track_att  * J_att  + ...
                obj.w_track_time * J_T    + ...
                obj.w_corridor   * J_corr + ...
                obj.w_slew       * J_du   + ...
                obj.w_curv       * J_ddu;
            opti.minimize(J);
        end

        %% ======================== DYNAMICS & KINEMATICS ========================

        function dyn_fnc = getCasADiDynamics(obj)
            %% GETCASADIDYNAMICS  6-DoF rigid body dynamics with LOX/IPA mass depletion.
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
            C_IB = quatRot(q).';
            TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            FI = C_IB * TB + [0; 0; -m * gSym];

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

            x_sym = [q; r; v; omegaB; m_lox; m_ipa];
            u_sym = [theta; phi; thrust; roll];
            params = [m_dry; gSym; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight; ...
                      J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot];
            xdot = [qdot; v; FI/m; wdot; mdl; mdi];
            dyn_fnc = Function('dyn_fnc', {x_sym, u_sym, params}, {xdot});
        end

        function name = getVehicleName(obj)
            if obj.Vehicle == 1, name = "TOAD"; else, name = "ASTRA"; end
        end

        %% ======================== PUBLIC API & DELEGATIONS ========================

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
                error('TrajectoryOptimizerHybrid:NoSolution', 'No solution available for export. Call solve() first.');
            end

            if nargin < 2 || isempty(filepath)
                if obj.SaveDir ~= ""
                    target_dir = obj.SaveDir;
                else
                    target_dir = fullfile(pwd, 'Guidance', 'Trajectories');
                end
                if ~exist(target_dir, 'dir'), mkdir(target_dir); end
                if obj.Filename ~= ""
                    fn = obj.Filename;
                else
                    fn = sprintf('%s_%s_v%03d', obj.getVehicleName(), obj.Maneuver, obj.Version);
                end
                filepath = fullfile(target_dir, fn + ".csv");
            end

            ts = obj.Solution.Time(:);
            Xr = obj.Solution.X;
            Ur = obj.Solution.U;

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

    methods (Static)
        function qout = vectorToQuat(v_in)
            %% VECTORTOQUAT  Compute unit quaternion rotating z_body [0;0;1] to v_in.
            v_in = v_in(:) / max(1e-6, norm(v_in));
            z_body = [0; 0; 1]; c = cross(z_body, v_in); d = dot(z_body, v_in);
            if d < -0.9999
                qout = [0; 0; 1; 0];
            else
                s = sqrt(2*(1+d));
                qout = [0.5*s; c/s];
                qout = qout / norm(qout);
            end
        end

        function R = quatToRot(q_in)
            %% QUATTOROT  Compute 3x3 direction cosine matrix from quaternion.
            w = q_in(1); x = q_in(2); y = q_in(3); z = q_in(4);
            R = [1-2*(y^2+z^2), 2*(x*y-w*z), 2*(x*z+w*y); ...
                 2*(x*y+w*z), 1-2*(x^2+z^2), 2*(y*z-w*x); ...
                 2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x^2+y^2)];
        end
    end
end
