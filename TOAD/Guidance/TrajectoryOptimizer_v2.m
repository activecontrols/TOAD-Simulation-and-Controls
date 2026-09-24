classdef TrajectoryOptimizer_v2 < handle
    %% TRAJECTORYOPTIMIZER_V2  Unified Two-Stage 6-DoF Optimal Trajectory Generator.
    % Stage 1: Closed-form differential flatness with dynamic actuator sizing (<10 ms).
    % Stage 2: Full 6-DoF direct collocation (IPOPT) with chatter-free regularizations.
    %
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    properties
        % System & Vehicle Properties
        constants           % Vehicle constants struct from LoadTOADParams
        Vehicle double = 0  % 0 for ASTRA (electric), 1 for TOAD (liquid biprop)
    end

    properties (Dependent)
        isElectric logical
    end

    properties
        % Discretization & Mesh
        N double = 80               % Number of control intervals (balanced resolution & speed)
        T_bounds = [8, 35]          % Total duration bounds [s]
        T_initial double = 16       % Initial duration guess [s]
        TimeParamMode string = "Single"

        % Formulation Settings
        CircleTightness double = 0.25   % Radial tolerance band [m]
        GlideslopeAngle double = 10     % Takeoff and landing cone half-angle [deg]
        FunnelCurvature double = 0.015  % Flaring curvature for landing funnel [1/m]

        % Multi-Objective Cost Function Weights
        w_time double   = 0.25      % Mission duration weight
        w_length double = 0.45      % 3D spatial path length weight
        w_effort double = 0.15      % Control effort weight (hover deviation + gimbal tilt)
        w_slew double   = 0.15      % Thrust/roll actuator slew rate regularization
        w_slew_gim double = 0.80    % Dedicated TVC gimbal slew penalty (kills wiggles)
        w_curv double   = 0.25      % Thrust/roll actuator curvature regularization
        w_curv_gim double = 0.50    % Dedicated TVC gimbal curvature penalty (kills chatter)
        w_smooth double = 0.02      % Velocity step smoothness
        w_rate double   = 0.025     % Body angular rate penalty
        w_qz double     = 0.025     % Yaw deflection penalty

        % Maneuver Definition
        Maneuver string = "Circle"  % "Circle", "Backflip", "Hop", "Waypoint", "Custom"
        ManeuverParams struct
        CustomWaypoints = []
        Waypoints double = []           % 3 x K matrix of waypoints [m]
        WaypointTolerances double = []  % 1 x K tolerance sphere radii [m]
        T_segments double = []          % 1 x M segment durations [s]
        UseMinimumSnap logical = true   % Enable closed-form minimum-snap spline
        TimeOptimization logical = true % Refine segment times via minimum snap

        % Boundary Conditions
        q0 = [1; 0; 0; 0]           % Initial attitude quaternion [w; x; y; z]
        r0 = [0; 0; 0]              % Launch pad position [m]
        v0 = [0; 0; 0]              % Launch pad velocity [m/s]
        w0 = [0; 0; 0]              % Initial angular rate [rad/s]
        r_f = [0; 0; 0]             % Landing zone position [m]
        v_f_tol = 0.25              % Max allowable touchdown speed [m/s]

        % Control Limits
        thrust_margin = 0.05
        gimbal_margin = 0.15
        max_gimbal_rate = deg2rad(30)   % Gimbal slew limit [rad/s]
        max_thrust_rate = 1000          % Thrust slew limit [N/s]
        max_roll_rate = 4               % Roll torque limit [N*m]
        max_gimbal_angle = pi/15        % Max physical gimbal angle [rad] (12 deg)

        % Scaling & Solver Options
        Sx double
        Su double
        L_c double = 50                 % Characteristic position length [m]
        MaxIter double = 120
        Tol double = 1e-2
        ConstrViolTol double = 5e-3
        PrintLevel double = 0

        % Storage
        Stage1Diag struct
        InitialGuess struct
        Solution struct
        OptiVars struct
    end

    methods
        function val = get.isElectric(obj), val = (obj.Vehicle == 0); end

        function obj = TrajectoryOptimizer_v2(constants6DoF, varargin)
            %% CONSTRUCTOR
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizer_v2 requires constants6DoF struct from LoadTOADParams.');
            end
            obj.constants = constants6DoF;

            p = inputParser;
            p.KeepUnmatched = true;
            addParameter(p, 'Vehicle', obj.Vehicle, @isnumeric);
            addParameter(p, 'Maneuver', obj.Maneuver, @(x) ischar(x)||isstring(x));
            addParameter(p, 'N', obj.N, @isnumeric);
            addParameter(p, 'T_initial', [], @isnumeric);
            addParameter(p, 'T_bounds', [], @isnumeric);
            addParameter(p, 'TimeParamMode', obj.TimeParamMode, @(x) ischar(x)||isstring(x));
            addParameter(p, 'CircleTightness', obj.CircleTightness, @isnumeric);
            addParameter(p, 'GlideslopeAngle', obj.GlideslopeAngle, @isnumeric);
            addParameter(p, 'FunnelCurvature', obj.FunnelCurvature, @isnumeric);
            addParameter(p, 'w_time', obj.w_time, @isnumeric);
            addParameter(p, 'w_length', obj.w_length, @isnumeric);
            addParameter(p, 'w_effort', obj.w_effort, @isnumeric);
            addParameter(p, 'w_slew', obj.w_slew, @isnumeric);
            addParameter(p, 'w_slew_gim', obj.w_slew_gim, @isnumeric);
            addParameter(p, 'w_curv', obj.w_curv, @isnumeric);
            addParameter(p, 'w_curv_gim', obj.w_curv_gim, @isnumeric);
            addParameter(p, 'w_smooth', obj.w_smooth, @isnumeric);
            addParameter(p, 'w_rate', obj.w_rate, @isnumeric);
            addParameter(p, 'w_qz', obj.w_qz, @isnumeric);
            addParameter(p, 'PrintLevel', obj.PrintLevel, @isnumeric);
            addParameter(p, 'MaxIter', obj.MaxIter, @isnumeric);
            addParameter(p, 'Tol', obj.Tol, @isnumeric);
            addParameter(p, 'ConstrViolTol', obj.ConstrViolTol, @isnumeric);
            parse(p, varargin{:});

            flds = {'Vehicle','Maneuver','N','TimeParamMode','CircleTightness','GlideslopeAngle', ...
                    'FunnelCurvature','w_time','w_length','w_effort','w_slew','w_slew_gim', ...
                    'w_curv','w_curv_gim','w_smooth','w_rate','w_qz','PrintLevel','MaxIter','Tol','ConstrViolTol'};
            for i = 1:numel(flds)
                f = flds{i}; obj.(f) = p.Results.(f);
            end

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

        function updateScaling(obj)
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

        function setBoundaries(obj, r0, r_f, varargin)
            obj.r0 = r0(:); obj.r_f = r_f(:);
            p = inputParser;
            addParameter(p, 'v0', [0; 0; 0], @isnumeric);
            addParameter(p, 'q0', [1; 0; 0; 0], @isnumeric);
            addParameter(p, 'v_f_tol', obj.v_f_tol, @isnumeric);
            parse(p, varargin{:});
            obj.v0 = p.Results.v0(:); obj.q0 = p.Results.q0(:);
            obj.v_f_tol = double(p.Results.v_f_tol);
            if obj.Maneuver == "Hop"
                obj.adaptHopTiming();
            end
            obj.InitialGuess = struct();
        end

        function setManeuver(obj, name, varargin)
            obj.Maneuver = string(name);
            p = inputParser;
            if obj.Vehicle == 0
                dR = 5.0;  dA = 7.0;  dDR = 2.5; dCB = [14.0, 30.0]; dCI = 20.0;
                dBA = 20.0; dBB = [10.0, 24.0]; dBI = 16.0;
                dHA = 15.0; % ASTRA Hop apex requirement: ~15 m
            else
                dR = 15.0; dA = 25.0; dDR = 5.0; dCB = [18.0, 42.0]; dCI = 26.0;
                dBA = 45.0; dBB = [16.0, 38.0]; dBI = 24.0;
                dHA = 50.0; % TOAD Hop apex requirement: ~50 m
            end

            switch obj.Maneuver
                case "Circle"
                    addParameter(p, 'circle_radius', dR, @isnumeric);
                    addParameter(p, 'circle_alt', dA, @isnumeric);
                    addParameter(p, 'circle_center', [0; 0], @isnumeric);
                    addParameter(p, 'max_descent_rate', dDR, @isnumeric);
                    addParameter(p, 'f_orbit_start', 0.25, @isnumeric);
                    addParameter(p, 'f_orbit_end', 0.75, @isnumeric);
                    parse(p, varargin{:});
                    obj.T_bounds = dCB; obj.T_initial = dCI;
                    fs = p.Results.f_orbit_start; fe = p.Results.f_orbit_end;
                    Ns = round(fs * obj.N); Ne = round(fe * obj.N);
                    obj.ManeuverParams = struct('circle_radius', double(p.Results.circle_radius), ...
                        'circle_alt', double(p.Results.circle_alt), 'circle_center', double(p.Results.circle_center(:)), ...
                        'max_descent_rate', double(p.Results.max_descent_rate), 'f_orbit_start', fs, ...
                        'f_orbit_end', fe, 'N_orbit_start', Ns, 'N_orbit_end', Ne);

                case "Backflip"
                    addParameter(p, 'apex_alt', dBA, @isnumeric);
                    addParameter(p, 'flip_start_frac', 0.35, @isnumeric);
                    addParameter(p, 'flip_end_frac', 0.65, @isnumeric);
                    addParameter(p, 'theta_tol', deg2rad(22), @isnumeric);
                    parse(p, varargin{:});
                    obj.T_bounds = dBB; obj.T_initial = dBI;
                    fa = p.Results.flip_start_frac; fb = p.Results.flip_end_frac;
                    Na = round(fa * obj.N); Nap = round(fb * obj.N); Nf = round(0.5 * (Na + Nap));
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt), ...
                        'flip_start_frac', fa, 'flip_end_frac', fb, 'theta_tol', double(p.Results.theta_tol), ...
                        'N_ascent', Na, 'N_flip', Nf, 'N_approach', Nap, 'q_inverted', [0; 0; -1; 0]);

                case "Hop"
                    addParameter(p, 'apex_alt', dHA, @isnumeric);
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('apex_alt', double(p.Results.apex_alt), ...
                        'N_ascent', round(0.5 * obj.N), 'N_descent', obj.N);
                    obj.adaptHopTiming();

                otherwise
                    obj.ManeuverParams = struct();
            end
            obj.InitialGuess = struct();
        end

        function adaptHopTiming(obj)
            %% ADAPTHOPTIMING  Physically scale duration bounds and initial guess from geometry.
            if obj.Maneuver ~= "Hop", return; end
            mp = obj.ManeuverParams;
            z_base = max(obj.r0(3), obj.r_f(3));
            d_horiz = norm(obj.r_f(1:2) - obj.r0(1:2));

            dh = max(1.5, mp.apex_alt - z_base);
            % Kinematic estimate: climb with vehicle-specific thrust-to-weight margin
            if obj.Vehicle == 1
                % TOAD (biprop): wet mass ~182 kg, MaxThrust ~2446 N
                % Realistic gentle climb ~0.8 m/s^2, descent ~0.8 m/s^2, 8s liftoff/landing flare
                a_climb = 0.8;
                a_desc  = 0.8;
                pad_time = 8.0;
            else
                % ASTRA (electric quad): climb ~1.8 m/s^2, descent ~1.5 m/s^2
                a_climb = max(1.5, 0.20 * obj.constants.g);
                a_desc  = max(1.5, 0.20 * obj.constants.g);
                pad_time = 4.0;
            end
            t_vert = sqrt(2 * dh / a_climb) + sqrt(2 * dh / a_desc);

            % Horizontal translation with 10-12 deg tilt: a_h ~ g * sin(10 deg) ~ 1.7 m/s^2
            a_h = max(1.0, obj.constants.g * sind(10));
            t_horiz = 2 * sqrt(max(0.5, d_horiz) / a_h);

            t_phys = max(t_vert, t_horiz) + pad_time;
            t_phys = max(10.0, min(55.0, t_phys));

            obj.T_initial = round(t_phys, 1);
            obj.T_bounds = [max(8.0, round(0.40 * t_phys, 1)), min(65.0, round(1.8 * t_phys, 1))];
        end

        function setWaypoints(obj, waypoints, varargin)
            %% SETWAYPOINTS  Configure arbitrary 3D multi-waypoint trajectory.
            %   waypoints: 3 x K matrix of 3D positions [m] (K >= 2).
            if size(waypoints, 1) ~= 3 || size(waypoints, 2) < 2
                error('Waypoints must be a 3 x K matrix with K >= 2.');
            end
            obj.Waypoints = double(waypoints);
            K = size(obj.Waypoints, 2);
            M = K - 1;

            p = inputParser;
            addParameter(p, 'T_segments', [], @isnumeric);
            addParameter(p, 'Tolerances', 0.25 * ones(1, K), @isnumeric);
            addParameter(p, 'T_total', obj.T_initial, @isnumeric);
            parse(p, varargin{:});

            obj.WaypointTolerances = double(p.Results.Tolerances);
            if isscalar(obj.WaypointTolerances)
                obj.WaypointTolerances = repmat(obj.WaypointTolerances, 1, K);
            end

            if ~isempty(p.Results.T_segments)
                obj.T_segments = double(p.Results.T_segments(:))';
                if numel(obj.T_segments) ~= M
                    error('T_segments must have length K - 1 = %d.', M);
                end
                obj.T_initial = sum(obj.T_segments);
            else
                T_tot = double(p.Results.T_total);
                obj.T_segments = obj.allocateSegmentTimes(obj.Waypoints, T_tot);
                obj.T_initial = sum(obj.T_segments);
            end

            obj.T_bounds = [max(4.0, round(0.50 * obj.T_initial, 1)), ...
                            min(70.0, round(2.0 * obj.T_initial, 1))];
            obj.Maneuver = "Waypoint";
            obj.r0 = obj.Waypoints(:, 1);
            obj.r_f = obj.Waypoints(:, end);
            obj.InitialGuess = struct();
        end

        function T_seg = allocateSegmentTimes(~, waypoints, T_total)
            %% ALLOCATESEGMENTTIMES  Heuristic segment time allocation based on Euclidean distance.
            % Based on Richter, Bry, and Roy (2016).
            M = size(waypoints, 2) - 1;
            dists = zeros(1, M);
            for m = 1:M
                dists(m) = norm(waypoints(:, m+1) - waypoints(:, m));
            end
            total_dist = sum(dists);
            if total_dist < 1e-3
                T_seg = (T_total / M) * ones(1, M);
                return;
            end

            % Initial allocation proportional to distance
            T_seg = T_total * (dists / total_dist);

            % Enforce minimum duration per segment
            min_T = 1.5;
            T_seg = max(min_T, T_seg);
            % Renormalize to T_total
            T_seg = T_total * (T_seg / sum(T_seg));
        end

        function T_opt = refineSegmentTimes(obj, waypoints, T_seg, v0, a0, vf, af)
            %% REFINESEGMENTTIMES  Refine segment times to minimize integrated snap.
            % Based on Richter, Bry, and Roy (2016) and Mellinger & Kumar (2011).
            M = numel(T_seg);
            if M <= 1, T_opt = T_seg; return; end
            if nargin < 4, v0 = obj.v0; end
            if nargin < 5, a0 = [0; 0; 0]; end
            if nargin < 6, vf = [0; 0; -0.15]; end
            if nargin < 7, af = [0; 0; 0]; end

            T_total = sum(T_seg);
            % Softmax mapping to maintain sum = T_total and T_m > 0
            alpha0 = log(max(0.1, T_seg(1:M-1)) / max(0.1, T_seg(M)));

            cost_fun = @(a) obj.evalSnapCost(waypoints, a, T_total, v0, a0, vf, af);
            opts = optimset('Display', 'off', 'MaxIter', 25, 'MaxFunEvals', 50, 'TolX', 1e-2);
            try
                alpha_opt = fminsearch(cost_fun, alpha0, opts);
                exp_a = [exp(alpha_opt(:))', 1.0];
                T_opt = T_total * (exp_a / sum(exp_a));
                min_T = 1.0;
                T_opt = max(min_T, T_opt);
                T_opt = T_total * (T_opt / sum(T_opt));
            catch
                T_opt = T_seg;
            end
        end

        function cost = evalSnapCost(obj, waypoints, alpha, T_total, v0, a0, vf, af)
            exp_a = [exp(alpha(:))', 1.0];
            T_s = T_total * (exp_a / sum(exp_a));
            try
                [~, ~, ~, ~, ~, ~, s_s] = obj.solveMinimumSnapSpline(waypoints, T_s, v0, a0, vf, af, 40);
                cost = sum(sum(s_s.^2, 1)) * (T_total / 40);
            catch
                cost = 1e9;
            end
        end

        %% ===================== STAGE 1: FLATNESS & DYNAMIC SIZING =====================

        function guess = generateInitialGuess(obj)
            %% GENERATEINITIALGUESS  Stage-1 closed-form flatness with dynamic actuator sizing.
            t0 = tic;
            Nn = obj.N + 1; Nc = obj.N;
            m_dry = obj.constants.m_dry; MT = obj.constants.MaxThrust;
            thr_max_allow = (1 - obj.thrust_margin) * MT;
            gim_max_allow = (1 - obj.gimbal_margin) * obj.max_gimbal_angle;

            T_cur = obj.T_initial;
            if obj.Vehicle == 0
                m_vec = m_dry * ones(1, Nn); ml = zeros(1, Nn); mi = zeros(1, Nn);
            else
                ml = linspace(obj.constants.OxMass, 0.25 * obj.constants.OxMass, Nn);
                mi = linspace(obj.constants.FuMass, 0.25 * obj.constants.FuMass, Nn);
                m_vec = m_dry + ml + mi;
            end

            for iter = 1:5
                tv = linspace(0, T_cur, Nn);
                dt = T_cur / Nc;
                dt_row = dt * ones(1, Nc);

                [r_s, v_s, ~, q_s, w_s, U_s] = obj.evalFlatOutputs(tv, dt_row, m_vec);

                max_T = max(U_s(3, :));
                max_gim = max(sqrt(U_s(1,:).^2 + U_s(2,:).^2));

                scale_factor = 1.0;
                if max_T > thr_max_allow, scale_factor = max(scale_factor, sqrt(max_T / thr_max_allow)); end
                if obj.Maneuver == "Backflip" && max_gim > gim_max_allow
                    scale_factor = max(scale_factor, sqrt(max_gim / gim_max_allow));
                end
                if scale_factor <= 1.03, break; end
                T_cur = T_cur * min(scale_factor, 1.30);
            end

            t_stg1 = toc(t0) * 1000;
            obj.Stage1Diag = struct('t_stage1_ms', t_stg1, 'T_total', T_cur, ...
                                    'max_thrust', max_T, 'max_gimbal_deg', rad2deg(max_gim));

            obj.T_initial = T_cur;
            obj.T_bounds = [min(obj.T_bounds(1), max(4.0, 0.50 * T_cur)), ...
                            max(obj.T_bounds(2), min(60.0, 1.80 * T_cur))];

            Xp = [q_s; r_s; v_s; w_s; ml; mi];
            guess = struct('Time', tv, 'X', Xp, 'U', U_s, 'T_total', T_cur, 'dt_row', dt_row, ...
                           'Xhat', Xp ./ obj.Sx, 'Uhat', U_s ./ obj.Su);
            obj.InitialGuess = guess;
        end

        function [r_s, v_s, a_s, q_s, w_s, U_s] = evalFlatOutputs(obj, tv, dt_row, m_vec)
            Nn = obj.N + 1; Nc = obj.N; g = obj.constants.g; Ts = tv(end);
            r_s = zeros(3, Nn); v_s = zeros(3, Nn); a_s = zeros(3, Nn);
            j_s = zeros(3, Nn); s_s = zeros(3, Nn);
            q_s = zeros(4, Nn); w_s = zeros(3, Nn); U_s = zeros(4, Nc);

            switch obj.Maneuver
                case "Waypoint"
                    K = size(obj.Waypoints, 2);
                    if isempty(obj.T_segments) || numel(obj.T_segments) ~= (K - 1)
                        obj.T_segments = obj.allocateSegmentTimes(obj.Waypoints, Ts);
                    else
                        obj.T_segments = Ts * (obj.T_segments / sum(obj.T_segments));
                    end
                    if obj.TimeOptimization && (K - 1) > 1
                        T_use = obj.refineSegmentTimes(obj.Waypoints, obj.T_segments);
                    else
                        T_use = obj.T_segments;
                    end
                    [~, ~, r_s, v_s, a_s, j_s, ~] = obj.solveMinimumSnapSpline( ...
                        obj.Waypoints, T_use, obj.v0, [0; 0; 0], [0; 0; -0.15], [0; 0; 0], Nn);

                case "Hop"
                    p = obj.ManeuverParams;
                    apex = [0.5 * (obj.r0(1) + obj.r_f(1)); 0.5 * (obj.r0(2) + obj.r_f(2)); p.apex_alt];
                    hop_wps = [obj.r0, apex, obj.r_f];
                    if obj.Vehicle == 1
                        t_ratio = 0.53;
                    else
                        t_ratio = 0.50;
                    end
                    T_hop = [t_ratio * Ts, (1 - t_ratio) * Ts];
                    if obj.TimeOptimization
                        T_hop = obj.refineSegmentTimes(hop_wps, T_hop);
                    end
                    [~, ~, r_s, v_s, a_s, j_s, ~] = obj.solveMinimumSnapSpline( ...
                        hop_wps, T_hop, obj.v0, [0; 0; 0], [0; 0; -0.15], [0; 0; 0], Nn);

                case "Circle"
                    p = obj.ManeuverParams;
                    cx = p.circle_center(1); cy = p.circle_center(2);
                    R = p.circle_radius;     h = p.circle_alt;
                    ks = max(3, round(p.f_orbit_start * obj.N)) + 1;
                    ke = min(obj.N - 2, round(p.f_orbit_end * obj.N)) + 1;
                    t1 = tv(ks); t2 = tv(ke); Torb = max(1e-2, t2 - t1); wc = 2 * pi / Torb;
                    re = [cx + R; cy; h]; ve = [0; R * wc; 0]; ae = [-R * wc^2; 0; 0];
                    je = [0; -R * wc^3; 0];

                    for k = 1:ks
                        tau = tv(k) / max(1e-2, t1);
                        [r_s(:,k), v_s(:,k), a_s(:,k), j_s(:,k), ~] = obj.evalSepticSpline( ...
                            obj.r0, [0; 0; 0.6], [0; 0; 0], [0; 0; 0], re, ve, ae, je, tau, t1);
                    end
                    for k = (ks + 1):ke
                        th = 2 * pi * (tv(k) - t1) / Torb;
                        r_s(:,k) = [cx + R * cos(th); cy + R * sin(th); h];
                        v_s(:,k) = [-R * wc * sin(th); R * wc * cos(th); 0];
                        a_s(:,k) = [-R * wc^2 * cos(th); -R * wc^2 * sin(th); 0];
                        j_s(:,k) = [R * wc^3 * sin(th); -R * wc^3 * cos(th); 0];
                    end
                    td = max(1e-2, tv(end) - t2);
                    for k = (ke + 1):Nn
                        tau = (tv(k) - t2) / td;
                        [r_s(:,k), v_s(:,k), a_s(:,k), j_s(:,k), ~] = obj.evalSepticSpline( ...
                            re, ve, ae, je, obj.r_f, [0; 0; -0.15], [0; 0; 0], [0; 0; 0], tau, td);
                    end

                case "Backflip"
                    p = obj.ManeuverParams;
                    apex = [0.5 * (obj.r0(1) + obj.r_f(1)); 0.5 * (obj.r0(2) + obj.r_f(2)); p.apex_alt];
                    ta = 0.5 * Ts; vh = (obj.r_f - obj.r0) / max(1e-2, Ts); va = [vh(1); vh(2); 0];

                    for k = 1:Nn
                        t = tv(k);
                        if t <= ta
                            tau = t / ta;
                            [r_s(:,k), v_s(:,k), a_s(:,k), j_s(:,k), s_s(:,k)] = obj.evalSepticSpline( ...
                                obj.r0, [0; 0; 0.5], [0; 0; 0], [0; 0; 0], apex, va, [0; 0; 0], [0; 0; 0], tau, ta);
                        else
                            tau = (t - ta) / (Ts - ta);
                            [r_s(:,k), v_s(:,k), a_s(:,k), j_s(:,k), s_s(:,k)] = obj.evalSepticSpline( ...
                                apex, va, [0; 0; 0], [0; 0; 0], obj.r_f, [0; 0; -0.15], [0; 0; 0], [0; 0; 0], tau, Ts - ta);
                        end
                    end

                otherwise
                    for k = 1:Nn
                        tau = tv(k) / Ts;
                        [r_s(:,k), v_s(:,k), a_s(:,k), j_s(:,k), s_s(:,k)] = obj.evalSepticSpline( ...
                            obj.r0, obj.v0, [0;0;0], [0;0;0], obj.r_f, [0;0;-0.15], [0;0;0], [0;0;0], tau, Ts);
                    end
            end

            r_s(:,1) = obj.r0; v_s(:,1) = obj.v0; r_s(:,end) = obj.r_f; v_s(:,end) = [0; 0; -0.05];

            for k = 1:Nc
                F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                U_s(3, k) = max(0.1, norm(F_req));
            end

            if obj.Maneuver == "Backflip"
                p = obj.ManeuverParams;
                ks = max(3, round(p.flip_start_frac * obj.N)) + 1;
                ke = min(obj.N - 2, round(p.flip_end_frac * obj.N)) + 1;
                Lf = max(2, ke - ks); Tf = tv(ke) - tv(ks);
                Jyy = obj.constants.J(2, 2); cgz = max(0.15, obj.constants.rTB);

                for k = 1:Nn
                    if k <= ks
                        q_s(:,k) = obj.q0; w_s(:,k) = [0; 0; 0];
                    elseif k <= ke
                        tau = (k - ks) / Lf;
                        th   = -2 * pi * (10 * tau^3 - 15 * tau^4 + 6 * tau^5);
                        dth  = -(2 * pi / Tf) * (30 * tau^2 - 60 * tau^3 + 30 * tau^4);
                        ddth = -(2 * pi / Tf^2) * (60 * tau - 180 * tau^2 + 120 * tau^3);
                        q_s(:,k) = [cos(th/2); 0; sin(th/2); 0]; w_s(:,k) = [0; dth; 0];
                        if k <= Nc
                            cos_th = cos(th);
                            if cos_th < 0
                                throttle_mult = 1.0 + 0.20 * cos_th;
                                U_s(3, k) = max(0.28 * obj.constants.MaxThrust, U_s(3, k) * throttle_mult);
                            end
                            sin_phi = (Jyy * ddth) / max(1.0, cgz * U_s(3, k));
                            phi_lim = sin((1 - obj.gimbal_margin) * obj.max_gimbal_angle);
                            U_s(2, k) = asin(max(-phi_lim, min(phi_lim, sin_phi)));
                        end
                    else
                        q_s(:,k) = -obj.q0; w_s(:,k) = [0; 0; 0];
                    end
                end
            else
                for k = 1:Nn
                    F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                    T_mag = norm(F_req);
                    b3 = F_req / max(1e-3, T_mag);
                    q_flat = obj.vectorToQuat(b3);
                    if k <= 4
                        alpha = (k - 1) / 4;
                        if dot(obj.q0, q_flat) < 0, q_flat = -q_flat; end
                        q_blend = (1 - alpha) * obj.q0 + alpha * q_flat;
                        q_s(:,k) = q_blend / norm(q_blend);
                    elseif k >= Nn - 3
                        alpha = (Nn - k) / 3;
                        if dot(obj.q0, q_flat) < 0, q_flat = -q_flat; end
                        q_blend = (1 - alpha) * obj.q0 + alpha * q_flat;
                        q_s(:,k) = q_blend / norm(q_blend);
                    else
                        q_s(:,k) = q_flat;
                    end

                    % Body angular velocity w from jerk j_s
                    b3_dot = (j_s(:,k) - b3 * (b3' * j_s(:,k))) * m_vec(k) / max(1e-3, T_mag);
                    R_IB = obj.quatToRot(q_s(:,k));
                    b3_dot_B = R_IB' * b3_dot;
                    w_s(:,k) = [-b3_dot_B(2); b3_dot_B(1); 0];
                end
                q_s(:,1) = obj.q0; q_s(:,end) = obj.q0;
                w_s(:,1) = obj.w0; w_s(:,end) = [0; 0; 0];

                % Synthesize TVC gimbal angles from torque balance
                J_mat = obj.constants.J;
                rTB = max(0.15, obj.constants.rTB);
                for k = 1:Nc
                    F_req = m_vec(k) * (a_s(:,k) + [0; 0; g]);
                    T_mag = norm(F_req);
                    U_s(3, k) = max(0.28 * obj.constants.MaxThrust, min(0.95 * obj.constants.MaxThrust, T_mag));

                    w_dot = (w_s(:,k+1) - w_s(:,k)) / dt_row(k);
                    w_mid = 0.5 * (w_s(:,k) + w_s(:,k+1));
                    tau_req = J_mat * w_dot + cross(w_mid, J_mat * w_mid);

                    denom = max(1.0, U_s(3, k) * rTB);
                    phi_lim = sin((1 - obj.gimbal_margin) * obj.max_gimbal_angle);
                    phi_cmd   = asin(max(-phi_lim, min(phi_lim,  tau_req(1) / denom)));
                    theta_cmd = asin(max(-phi_lim, min(phi_lim, -tau_req(2) / denom)));

                    U_s(1, k) = theta_cmd;
                    U_s(2, k) = phi_cmd;
                    U_s(4, k) = 0;
                end
            end
        end

        function [poly_coeffs, tv, r_s, v_s, a_s, j_s, s_s] = solveMinimumSnapSpline(~, waypoints, T_segments, v0, a0, vf, af, N_eval)
            %% SOLVEMINIMUMSNAPSPLINE  Closed-form multi-segment septic minimum-snap spline.
            % Richter, C., Bry, A., and Roy, N. (2016) and Mellinger & Kumar (2011).
            M = size(waypoints, 2) - 1;
            N_deg = 7; N_c = 8;
            Q_blocks = cell(M, 1); A_blocks = cell(M, 1);
            for m = 1:M
                Tm = T_segments(m);
                Qm = zeros(N_c, N_c);
                for i = 4:N_deg
                    for j = 4:N_deg
                        Qm(i+1, j+1) = (i*(i-1)*(i-2)*(i-3)*j*(j-1)*(j-2)*(j-3)) / (i+j-7) * Tm^(i+j-7);
                    end
                end
                Q_blocks{m} = Qm;
                Am = zeros(N_c, N_c);
                Am(1, 1) = 1; Am(2, 2) = 1; Am(3, 3) = 2; Am(4, 4) = 6;
                for k = 0:N_deg
                    Am(5, k+1) = Tm^k;
                    if k >= 1, Am(6, k+1) = k * Tm^(k-1); end
                    if k >= 2, Am(7, k+1) = k*(k-1) * Tm^(k-2); end
                    if k >= 3, Am(8, k+1) = k*(k-1)*(k-2) * Tm^(k-3); end
                end
                A_blocks{m} = Am;
            end
            Q = blkdiag(Q_blocks{:}); A = blkdiag(A_blocks{:});
            A_inv = inv(A); R = A_inv' * Q * A_inv;
            n_fixed = M + 7; n_free = 3 * (M - 1);
            C = zeros(8 * M, n_fixed + n_free);
            C(1:4, 1:4) = eye(4);
            for m = 1:(M - 1)
                r_idx = 4 + m; v_idx = n_fixed + (m-1)*3 + 1;
                a_idx = n_fixed + (m-1)*3 + 2; j_idx = n_fixed + (m-1)*3 + 3;
                row_end = (m - 1)*8 + 5;
                C(row_end:row_end+3, [r_idx, v_idx, a_idx, j_idx]) = eye(4);
                row_start = m*8 + 1;
                C(row_start:row_start+3, [r_idx, v_idx, a_idx, j_idx]) = eye(4);
            end
            row_last = (M - 1)*8 + 5;
            C(row_last:row_last+3, (4+M):(4+M+3)) = eye(4);

            R_tilde = C' * R * C;
            R_PF = R_tilde(n_fixed+1:end, 1:n_fixed);
            R_PP = R_tilde(n_fixed+1:end, n_fixed+1:end);

            poly_coeffs = zeros(8 * M, 3);
            j0 = [0; 0; 0]; jf = [0; 0; 0];
            for dim = 1:3
                d_F = zeros(n_fixed, 1);
                d_F(1:4) = [waypoints(dim, 1); v0(dim); a0(dim); j0(dim)];
                for m = 1:(M - 1), d_F(4 + m) = waypoints(dim, m + 1); end
                d_F((4+M):(4+M+3)) = [waypoints(dim, end); vf(dim); af(dim); jf(dim)];
                if n_free > 0
                    d_P = -R_PP \ (R_PF * d_F);
                    d_opt = C * [d_F; d_P];
                else
                    d_opt = C * d_F;
                end
                poly_coeffs(:, dim) = A \ d_opt;
            end

            T_total = sum(T_segments);
            tv = linspace(0, T_total, N_eval);
            r_s = zeros(3, N_eval); v_s = zeros(3, N_eval);
            a_s = zeros(3, N_eval); j_s = zeros(3, N_eval); s_s = zeros(3, N_eval);
            cum_T = [0, cumsum(T_segments)];
            for k = 1:N_eval
                t = tv(k);
                m = find(t >= cum_T(1:end-1) & t <= cum_T(2:end), 1, 'last');
                if isempty(m), m = M; end
                m = min(m, M);
                tau = t - cum_T(m);
                c_m = poly_coeffs((m - 1)*8 + (1:8), :);
                p = [1; tau; tau^2; tau^3; tau^4; tau^5; tau^6; tau^7];
                p1 = [0; 1; 2*tau; 3*tau^2; 4*tau^3; 5*tau^4; 6*tau^5; 7*tau^6];
                p2 = [0; 0; 2; 6*tau; 12*tau^2; 20*tau^3; 30*tau^4; 42*tau^5];
                p3 = [0; 0; 0; 6; 24*tau; 60*tau^2; 120*tau^3; 210*tau^4];
                p4 = [0; 0; 0; 0; 24; 120*tau; 360*tau^2; 840*tau^3];
                r_s(:, k) = c_m' * p;
                v_s(:, k) = c_m' * p1;
                a_s(:, k) = c_m' * p2;
                j_s(:, k) = c_m' * p3;
                s_s(:, k) = c_m' * p4;
            end
        end

        function [r, v, a, j, s] = evalSepticSpline(~, r0, v0, a0, j0, r1, v1, a1, j1, tau, T)
            tau = max(0, min(1, tau)); h = r1 - r0;
            T2 = T^2; T3 = T^3; T4 = T^4;
            v0_s = v0 * T;   v1_s = v1 * T;
            a0_s = a0 * T2;  a1_s = a1 * T2;
            j0_s = j0 * T3;  j1_s = j1 * T3;

            dr = h - (v0_s + 0.5 * a0_s + (1/6) * j0_s);
            dv = v1_s - (v0_s + a0_s + 0.5 * j0_s);
            da = a1_s - (a0_s + j0_s);
            dj = j1_s - j0_s;

            c4 =  35 * dr - 15 * dv + 2.5 * da - (1/6) * dj;
            c5 = -84 * dr + 39 * dv - 7.0 * da + 0.5   * dj;
            c6 =  70 * dr - 34 * dv + 6.5 * da - 0.5   * dj;
            c7 = -20 * dr + 10 * dv - 2.0 * da + (1/6) * dj;

            tau2 = tau * tau; tau3 = tau2 * tau; tau4 = tau3 * tau;
            tau5 = tau4 * tau; tau6 = tau5 * tau; tau7 = tau6 * tau;

            r = r0 + v0_s * tau + 0.5 * a0_s * tau2 + (1/6) * j0_s * tau3 + ...
                c4 * tau4 + c5 * tau5 + c6 * tau6 + c7 * tau7;
            v = (v0_s + a0_s * tau + 0.5 * j0_s * tau2 + ...
                 4 * c4 * tau3 + 5 * c5 * tau4 + 6 * c6 * tau5 + 7 * c7 * tau6) / T;
            a = (a0_s + j0_s * tau + ...
                 12 * c4 * tau2 + 20 * c5 * tau3 + 30 * c6 * tau4 + 42 * c7 * tau5) / T2;
            j = (j0_s + 24 * c4 * tau + 60 * c5 * tau2 + 120 * c6 * tau3 + 210 * c7 * tau4) / T3;
            s = (24 * c4 + 120 * c5 * tau + 360 * c6 * tau2 + 840 * c7 * tau3) / T4;
        end

        function q = vectorToQuat(~, v)
            v = v(:) / max(1e-6, norm(v));
            z_body = [0; 0; 1]; c = cross(z_body, v); d = dot(z_body, v);
            if d < -0.9999, q = [0; 0; 1; 0];
            else, s = sqrt(2 * (1 + d)); q = [0.5 * s; c / s]; q = q / norm(q); end
        end

        function R = quatToRot(~, q)
            w = q(1); x = q(2); y = q(3); z = q(4);
            R = [1 - 2*(y^2 + z^2), 2*(x*y - w*z), 2*(x*z + w*y); ...
                 2*(x*y + w*z), 1 - 2*(x^2 + z^2), 2*(y*z - w*x); ...
                 2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x^2 + y^2)];
        end

        %% ===================== STAGE 2: 6-DOF DIRECT COLLOCATION NLP =====================

        function sol = optimize(obj)
            %% OPTIMIZE  Alias for solve()
            sol = obj.solve();
        end

        function sol = solve(obj)
            %% SOLVE  Two-Stage Solver: Stage-1 Flatness Seed -> Stage-2 6-DoF Collocation
            import casadi.*
            t_total_start = tic;

            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'Xhat')
                obj.generateInitialGuess();
            end

            N = obj.N; opti = Opti(); dyn_fnc = obj.getCasADiDynamics();
            obj.OptiVars = struct();

            % 1. Uniform Timing Parameterization (Direct Total Duration)
            T_total = opti.variable();
            opti.subject_to(obj.T_bounds(1) <= T_total);
            opti.subject_to(T_total <= obj.T_bounds(2));
            opti.set_initial(T_total, obj.InitialGuess.T_total);
            dt_row = repmat(T_total / N, 1, N);
            obj.OptiVars.T_total = T_total; obj.OptiVars.dt_row = dt_row;

            Xhat = opti.variable(15, N + 1); Uhat = opti.variable(4, N);
            X = obj.Sx .* Xhat; U = obj.Su .* Uhat;

            % 2. RK4 Discretization
            params_val = [obj.constants.m_dry; obj.constants.g; obj.constants.rTB; ...
                obj.constants.Ox_Z; obj.constants.OxMass; obj.constants.OxHeight; ...
                obj.constants.Fu_Z; obj.constants.FuMass; obj.constants.FuHeight; ...
                obj.constants.J(:); obj.constants.OxRadius; obj.constants.FuRadius; ...
                obj.constants.MaxThrust; obj.constants.OF; obj.constants.MaxMdot; ...
                0; zeros(9, 1); zeros(3, 1)];

            xs = MX.sym('xh', 15); us = MX.sym('uh', 4); ds = MX.sym('dt');
            xp = obj.Sx .* xs;     up = obj.Su .* us;
            k1 = dyn_fnc(xp,                 up, params_val);
            k2 = dyn_fnc(xp + ds / 2 * k1,   up, params_val);
            k3 = dyn_fnc(xp + ds / 2 * k2,   up, params_val);
            k4 = dyn_fnc(xp + ds * k3,       up, params_val);
            xn = xp + ds / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
            xn = [xn(1:4) / sqrt(sum(xn(1:4).^2) + 1e-12); xn(5:end)];
            Fstep = Function('F_step', {xs, us, ds}, {xn ./ obj.Sx});
            Fmap  = Fstep.map(N);
            opti.subject_to(Xhat(:, 2:end) == Fmap(Xhat(:, 1:N), Uhat, dt_row));

            % 3. State & Boundary Constraints
            opti.subject_to(-40 <= X(5, :)); opti.subject_to(X(5, :) <= 40);
            opti.subject_to(-40 <= X(6, :)); opti.subject_to(X(6, :) <= 40);
            alt_ceil = 85 * (obj.Vehicle == 0) + 160 * (obj.Vehicle == 1);
            opti.subject_to(-1 <= X(7, :)); opti.subject_to(X(7, :) <= alt_ceil);

            m_lox0 = obj.constants.OxMass; m_ipa0 = obj.constants.FuMass;
            opti.subject_to(X(:, 1) == [obj.q0; obj.r0; obj.v0; obj.w0; m_lox0; m_ipa0]);
            if obj.Maneuver == "Backflip", opti.subject_to(X(1:4, end) == -obj.q0);
            else, opti.subject_to(X(1:4, end) == obj.q0); end
            opti.subject_to(X(5:7, end) == obj.r_f);
            opti.subject_to(sum(X(8:10, end).^2) <= obj.v_f_tol^2);

            if obj.Vehicle == 1
                opti.subject_to(X(14, end) >= 0.05 * m_lox0);
                opti.subject_to(X(15, end) >= 0.05 * m_ipa0);
            else
                opti.subject_to(X(14:15, :) == 0);
            end

            % 4. Control Bounds & Rates
            MT = obj.constants.MaxThrust; tm = obj.thrust_margin;
            gm = obj.gimbal_margin;      mg = obj.max_gimbal_angle;
            opti.subject_to((0.25 + tm) * MT <= U(3, :)); opti.subject_to(U(3, :) <= (1 - tm) * MT);
            opti.subject_to(-(1 - gm) * mg <= U(1, :));   opti.subject_to(U(1, :) <= (1 - gm) * mg);
            opti.subject_to(-(1 - gm) * mg <= U(2, :));   opti.subject_to(U(2, :) <= (1 - gm) * mg);
            opti.subject_to(-(1 - tm) * obj.max_roll_rate <= U(4, :));
            opti.subject_to(U(4, :) <= (1 - tm) * obj.max_roll_rate);

            dU = U(:, 2:end) - U(:, 1:end-1); dts = dt_row(1:end-1);
            for ch = [1, 2]
                opti.subject_to(-obj.max_gimbal_rate * dts <= dU(ch, :));
                opti.subject_to(dU(ch, :) <= obj.max_gimbal_rate * dts);
            end
            opti.subject_to(-obj.max_thrust_rate * dts <= dU(3, :)); opti.subject_to(dU(3, :) <= obj.max_thrust_rate * dts);
            opti.subject_to(-obj.max_roll_rate * dts <= dU(4, :));   opti.subject_to(dU(4, :) <= obj.max_roll_rate * dts);

            % 5. Maneuvers & Multi-Objective Cost
            obj.applyManeuverConstraints(opti, X, U, dt_row);
            obj.applyCostFunction(opti, X, Xhat, U, Uhat, T_total, dt_row);

            % Seed Initial Guess
            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);

            % 6. Refinement-Oriented IPOPT Solver Setup
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, 'tol', obj.Tol, ...
                'constr_viol_tol', obj.ConstrViolTol, 'acceptable_tol', 1e-2, ...
                'acceptable_constr_viol_tol', 1e-2, 'acceptable_iter', 5, ...
                'max_cpu_time', 40.0, ...
                'mu_strategy', 'adaptive', 'print_level', obj.PrintLevel);
            opti.solver('ipopt', p_opts, s_opts);

            t_stg2_start = tic;
            try
                sc = opti.solve(); status = 'Success';
                X_r = obj.Sx .* sc.value(Xhat); U_r = obj.Su .* sc.value(Uhat);
                T_r = sc.value(T_total); dt_r = sc.value(dt_row); so = sc;
                stats = sc.stats(); stats.t_wall_total = toc(t_stg2_start);
            catch ME
                status = 'Infeasible';
                X_r = obj.Sx .* opti.debug.value(Xhat); U_r = obj.Su .* opti.debug.value(Uhat);
                T_r = opti.debug.value(T_total); dt_r = opti.debug.value(dt_row); so = opti;
                stats = struct('t_wall_total', toc(t_stg2_start), 'iter_count', obj.MaxIter, 'return_status', 'Infeasible', 'error_msg', ME.message);
            end

            t_total_all = toc(t_total_start);
            stats.t_stage1_ms = obj.Stage1Diag.t_stage1_ms;
            stats.t_stage2_s  = stats.t_wall_total;
            stats.t_total_s   = t_total_all;
            stats.stage1_diag = obj.Stage1Diag;

            dr_r = diff(X_r(5:7, :), 1, 2); LP = sum(sqrt(sum(dr_r.^2, 1)));
            t_vec = [0, cumsum(dt_r)];
            sol = struct('Status', status, 'T_total', T_r, 'L_path', LP, ...
                         'X', X_r, 'x', X_r, 'U', U_r, 'u', U_r, ...
                         'Time', t_vec, 't', t_vec, 'stats', stats, 'opti_obj', so);
            obj.Solution = sol;
        end

        function applyManeuverConstraints(obj, opti, X, ~, ~)
            %% APPLYMANEUVERCONSTRAINTS  Universal corridors & strict vertical monotonicity
            p = obj.ManeuverParams; N = obj.N;
            GSA = tand(obj.GlideslopeAngle); cf = obj.FunnelCurvature;
            r0_p = obj.r0; rf_p = obj.r_f;

            % 1. Universal Takeoff & Landing Subsets
            N_to   = max(3, round(0.15 * N));
            if obj.Maneuver == "Hop"
                N_land = min(N - 2, round(0.92 * N));
            else
                N_land = min(N - 2, round(0.85 * N));
            end

            % 2. Universal Takeoff Glideslope Cone
            dz_to = X(7, 1:N_to) - r0_p(3);
            rxy_to = sqrt((X(5, 1:N_to) - r0_p(1)).^2 + (X(6, 1:N_to) - r0_p(2)).^2 + 1e-4);
            opti.subject_to(rxy_to <= dz_to * GSA + 0.15);

            % Dedicated Liftoff Alignment Column
            N_liftoff = max(2, round(0.04 * N));
            opti.subject_to(rxy_to(1:N_liftoff) <= 0.15);
            opti.subject_to(abs(X(8, 1:N_liftoff)) <= 0.25);
            opti.subject_to(abs(X(9, 1:N_liftoff)) <= 0.25);
            R33_to = X(1, 1:N_liftoff).^2 - X(2, 1:N_liftoff).^2 - X(3, 1:N_liftoff).^2 + X(4, 1:N_liftoff).^2;
            opti.subject_to(R33_to >= cosd(15));

            % 3. Universal Landing Funnel
            dz_land = X(7, N_land:end) - rf_p(3);
            rxy_land = sqrt((X(5, N_land:end) - rf_p(1)).^2 + (X(6, N_land:end) - rf_p(2)).^2 + 1e-4);
            if obj.Maneuver == "Hop"
                GSA_land = tand(max(20, obj.GlideslopeAngle));
                opti.subject_to(rxy_land <= dz_land * GSA_land + cf * dz_land.^2 + 0.25);
            else
                opti.subject_to(rxy_land <= dz_land * GSA + cf * dz_land.^2 + 0.15);
            end
            v_desc_lim = 3.0 * (obj.Vehicle == 0) + 6.0 * (obj.Vehicle == 1);
            opti.subject_to(X(10, N_land:end) >= -v_desc_lim);
            R33_land = X(1, N_land:end).^2 - X(2, N_land:end).^2 - X(3, N_land:end).^2 + X(4, N_land:end).^2;
            opti.subject_to(R33_land >= cosd(22));

            % Dedicated Terminal Landing Alignment
            opti.subject_to(abs(X(8, end-2:end)) <= 0.25);
            opti.subject_to(abs(X(9, end-2:end)) <= 0.25);
            R33_term = X(1, end-2:end).^2 - X(2, end-2:end).^2 - X(3, end-2:end).^2 + X(4, end-2:end).^2;
            opti.subject_to(R33_term >= cosd(12));

            % 4. Strict Vertical Velocity Monotonicity (Disjoint partitions - zero collinear rows)
            switch obj.Maneuver
                case "Circle"
                    Nos = p.N_orbit_start; Noe = p.N_orbit_end; Norb = Noe - Nos;
                    cx = p.circle_center(1); cy = p.circle_center(2);
                    R  = p.circle_radius;     h  = p.circle_alt;

                    if Nos > N_to + 1
                        opti.subject_to(X(10, 1:N_to-1) >= 0.0);
                        opti.subject_to(X(10, N_to:Nos-2) >= 0.15);
                        opti.subject_to(X(10, Nos-1:Nos) >= 0.0);
                    else
                        opti.subject_to(X(10, 1:Nos) >= 0.0);
                    end

                    if N_land > Noe + 1
                        opti.subject_to(X(10, Noe:Noe+1) <= 0.0);
                        opti.subject_to(X(10, Noe+2:N_land) <= -0.15);
                        opti.subject_to(X(10, N_land+1:end) <= 0.0);
                    else
                        opti.subject_to(X(10, Noe:end) <= 0.0);
                    end

                    rsq = (X(5, Nos:Noe) - cx).^2 + (X(6, Nos:Noe) - cy).^2;
                    tr = obj.CircleTightness;
                    opti.subject_to((R - tr)^2 <= rsq); opti.subject_to(rsq <= (R + tr)^2);
                    opti.subject_to(abs(X(7, Nos:Noe) - h) <= 0.25);

                    rx = X(5, Nos:Noe) - cx; ry = X(6, Nos:Noe) - cy;
                    cprog = rx(1:end-1) .* ry(2:end) - ry(1:end-1) .* rx(2:end);
                    opti.subject_to(cprog >= (R^2) * sin(2 * pi / Norb * 0.70));

                    kq1 = Nos + round(0.25*Norb); kq2 = Nos + round(0.50*Norb); kq3 = Nos + round(0.75*Norb);
                    opti.subject_to(X(6, kq1) - cy >= 0);
                    opti.subject_to(X(5, kq2) - cx <= 0);
                    opti.subject_to(X(6, kq3) - cy <= 0);
                    opti.subject_to((X(5, Noe) - (cx + R))^2 + (X(6, Noe) - cy)^2 <= tr^2);

                case "Backflip"
                    Na = p.N_ascent; Nf = p.N_flip; Nap = p.N_approach;


                    if Na > N_to + 1
                        opti.subject_to(X(10, 1:N_to-1) >= 0.0);
                        opti.subject_to(X(10, N_to:Na-2) >= 0.15);
                        opti.subject_to(X(10, Na-1:Na) >= 0.0);
                    else
                        opti.subject_to(X(10, 1:Na) >= 0.0);
                    end

                    opti.subject_to(X(7, :) <= p.apex_alt + 3.0);
                    opti.subject_to(X(7, Nf) >= p.apex_alt - 2.5);
                    opti.subject_to(X(7, Na) >= 0.50 * p.apex_alt);
                    att_tol = cos(p.theta_tol / 2);
                    opti.subject_to(p.q_inverted' * X(1:4, Nf) >= att_tol);
                    opti.subject_to([-1; 0; 0; 0]' * X(1:4, Nap) >= att_tol);
                    opti.subject_to(X(12, Na:Nf+2) <= 0.05);

                    if N_land > Nap + 1
                        opti.subject_to(X(10, Nap:Nap+1) <= 0.0);
                        opti.subject_to(X(10, Nap+2:N_land) <= -0.15);
                        opti.subject_to(X(10, N_land+1:end) <= 0.0);
                    else
                        opti.subject_to(X(10, Nap:end) <= 0.0);
                    end

                case "Hop"
                    % Determine apex node from Stage 1 initial guess
                    if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                        [~, Nap_h] = max(obj.InitialGuess.X(7, :));
                    else
                        Nap_h = round(0.50 * N);
                    end
                    Nap_h = max(N_to + 3, min(N - 8, Nap_h));

                    % Apex altitude target and overall ceiling
                    opti.subject_to(X(7, Nap_h) >= p.apex_alt - 1.5);
                    opti.subject_to(X(7, :) <= p.apex_alt + 3.0);

                    % Upright attitude margin
                    R33 = X(1, :).^2 - X(2, :).^2 - X(3, :).^2 + X(4, :).^2;
                    opti.subject_to(R33 >= cosd(42));

                case "Waypoint"
                    % Waypoint passage spherical tolerance constraints
                    if ~isempty(obj.Waypoints)
                        K = size(obj.Waypoints, 2);
                        cum_T = [0, cumsum(obj.T_segments)];
                        T_tot_guess = cum_T(end);
                        for m = 2:(K - 1)
                            t_wp = cum_T(m);
                            k_wp = max(1, min(N, round(N * t_wp / T_tot_guess) + 1));
                            tol_m = obj.WaypointTolerances(m);
                            wp_pos = obj.Waypoints(:, m);
                            opti.subject_to(sum((X(5:7, k_wp) - wp_pos).^2) <= tol_m^2);
                        end
                    end
                    % Upright attitude margin
                    R33_wp = X(1, :).^2 - X(2, :).^2 - X(3, :).^2 + X(4, :).^2;
                    opti.subject_to(R33_wp >= cosd(45));
            end
        end

        function applyCostFunction(obj, opti, X, Xhat, U, Uhat, T_total, dt_row)
            %% APPLYCOSTFUNCTION  Unified multi-objective cost with dedicated gimbal curvature
            N = obj.N; g = obj.constants.g; MT = obj.constants.MaxThrust;

            % 1. Time Optimality
            J_time = T_total / obj.T_initial;

            % 2. 3D Spatial Path Length (Promotes natural planarity)
            dr = X(5:7, 2:end) - X(5:7, 1:end-1);
            L_path = sum(sqrt(sum(dr.^2, 1) + 1e-6));
            if obj.Maneuver == "Circle"
                L_ref = norm(obj.r_f - obj.r0) + 2 * pi * obj.ManeuverParams.circle_radius + 2 * obj.ManeuverParams.circle_alt;
            elseif obj.Maneuver == "Waypoint" && ~isempty(obj.Waypoints)
                L_ref = sum(sqrt(sum(diff(obj.Waypoints, 1, 2).^2, 1))) + 5.0;
            else
                L_ref = max(1.0, norm(obj.r_f - obj.r0) + 20.0);
            end
            J_length = L_path / max(1.0, L_ref);

            % 3. Control Effort (Hover baseline + gimbal tilt)
            if obj.Vehicle == 0, mv = obj.constants.m_dry * ones(1, N);
            else, mv = obj.constants.m_dry + X(14, 1:N) + X(15, 1:N); end
            uh = mv * g / MT; gb = max(1e-3, (1 - obj.gimbal_margin) * obj.max_gimbal_angle);
            esq = (U(3, :) / MT - uh).^2 + (U(1, :) / gb).^2 + (U(2, :) / gb).^2 + ...
                  (U(4, :) / max(1e-3, (1 - obj.thrust_margin) * obj.max_roll_rate)).^2;
            J_effort = sum(dt_row .* esq) / (0.50 * obj.T_initial);

            % 4. Actuator Slew Regularization (Dedicated weighting on TVC gimbal vs thrust/roll)
            dUh = Uhat(:, 2:end) - Uhat(:, 1:end-1);
            J_slew_gim = sum(sum(dUh(1:2, :).^2, 1)) / (N - 1);
            J_slew_thr = sum(dUh(3, :).^2) / (N - 1);
            J_slew_roll = sum(dUh(4, :).^2) / (N - 1);

            % 5. Actuator Curvature Regularization (Dedicated weighting on gimbal channels kills chatter!)
            d2Uh = Uhat(:, 3:end) - 2 * Uhat(:, 2:end-1) + Uhat(:, 1:end-2);
            J_curv_gim = sum(sum(d2Uh(1:2, :).^2, 1)) / max(1, N - 2);
            J_curv_thr = sum(d2Uh(3, :).^2) / max(1, N - 2);
            J_curv_roll = sum(d2Uh(4, :).^2) / max(1, N - 2);

            % 6. Kinematic Velocity Step Regularization
            dVh = Xhat(8:10, 2:end) - Xhat(8:10, 1:end-1);
            J_smooth = sum(sum(dVh.^2, 1)) / (N - 1);

            % 7. Angular Rates & Yaw Deflection
            J_rate = sum(sum(Xhat(11:13, 1:end-1).^2, 1)) / N;
            J_qz   = sum(Xhat(4, :).^2) / (N + 1);

            J = obj.w_time     * J_time     + ...
                obj.w_length   * J_length   + ...
                obj.w_effort   * J_effort   + ...
                obj.w_slew_gim * J_slew_gim + ...
                obj.w_slew     * (J_slew_thr + J_slew_roll) + ...
                obj.w_curv_gim * J_curv_gim + ...
                obj.w_curv     * (J_curv_thr + J_curv_roll) + ...
                obj.w_smooth   * J_smooth   + ...
                obj.w_rate     * J_rate     + ...
                obj.w_qz       * J_qz;
            opti.minimize(J);
        end

        %% ===================== DYNAMICS & KINEMATICS =====================

        function dyn_fnc = getCasADiDynamics(obj)
            import casadi.*
            q = MX.sym('q', 4); r = MX.sym('r', 3); v = MX.sym('v', 3);
            omegaB = MX.sym('omegaB', 3); m_lox = MX.sym('m_lox'); m_ipa = MX.sym('m_ipa');
            m_dry = MX.sym('m_dry'); g = MX.sym('g'); rTB = MX.sym('rTB');
            Ox_Z = MX.sym('Ox_Z'); OxMassI = MX.sym('OxMassI'); OxHeight = MX.sym('OxHeight');
            Fu_Z = MX.sym('Fu_Z'); FuMassI = MX.sym('FuMassI'); FuHeight = MX.sym('FuHeight');
            J = MX.sym('J', 3, 3); OxRadius = MX.sym('OxRadius'); FuRadius = MX.sym('FuRadius');
            MaxThrust = MX.sym('MaxThrust'); OF = MX.sym('OF'); MaxMdot = MX.sym('MaxMdot');
            MaxMdot_d = MX.sym('MaxMdot_d'); J_d = MX.sym('J_d', 3, 3); TB_d = MX.sym('TB_d', 3);
            theta = MX.sym('theta'); phi = MX.sym('phi'); thrust = MX.sym('thrust'); roll = MX.sym('roll');

            m = m_dry + m_lox + m_ipa;
            C_IB = quatRot(q).';
            TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            FI = C_IB * TB + [0; 0; -m * g];

            if obj.Vehicle == 0, mdl = 0; mdi = 0; OFH = 0; FFH = 0;
            else
                mdl = -thrust / MaxThrust * OF / (1 + OF) * (MaxMdot + MaxMdot_d);
                mdi = -thrust / MaxThrust * 1 / (1 + OF) * (MaxMdot + MaxMdot_d);
                OFH = (m_lox / OxMassI) * OxHeight * 0.9;
                FFH = (m_ipa / FuMassI) * FuHeight * 0.9;
            end

            Jlox = diag([1/12*m_lox*(3*OxRadius^2 + OFH^2), 1/12*m_lox*(3*OxRadius^2 + OFH^2), 1/2*m_lox*OxRadius^2]);
            Jipa = diag([1/12*m_ipa*(3*FuRadius^2 + FFH^2), 1/12*m_ipa*(3*FuRadius^2 + FFH^2), 1/2*m_ipa*FuRadius^2]);
            OFL = Ox_Z + OFH / 2; FFL = Fu_Z + FFH / 2;
            CGz = (m_dry * rTB + m_lox * OFL + m_ipa * FFL) / m;
            dd = rTB - CGz + TB_d(3); dl = OFL - CGz; di = FFL - CGz;
            J_tot = J + m_dry * diag([dd^2, dd^2, 0]) + Jlox + m_lox * diag([dl^2, dl^2, 0]) + Jipa + m_ipa * diag([di^2, di^2, 0]) + J_d;

            tDir = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            if obj.Vehicle == 0, MB = zetaCross([0; 0; -CGz] + TB_d) * TB + roll * tDir;
            else, MB = zetaCross([0; 0; -CGz] + TB_d) * TB + [0; 0; roll]; end

            qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
            wdot = (MB - zetaCross(omegaB) * J_tot * omegaB) ./ diag(J_tot);

            x = [q; r; v; omegaB; m_lox; m_ipa]; u = [theta; phi; thrust; roll];
            params = [m_dry; g; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight; ...
                      J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot; MaxMdot_d; J_d(:); TB_d];
            xdot = [qdot; v; FI / m; wdot; mdl; mdi];
            dyn_fnc = Function('dyn_fnc', {x, u, params}, {xdot});
        end
    end
end
