classdef TrajectoryOptimizerHybrid < handle
    %% TrajectoryOptimizerHybrid - Two-Stage 6-DoF Optimal Trajectory Generator.
    % The following file contains all the methods and algorithms required
    % for a two stage trajectory generator, which combines a
    % polynomial-based trajectory generator via solving QP programs for an
    % initial guess, and then passes it through a Direct Collocation stage
    % using CasADi to refine trajectories and map commands to our
    % actuators. Meant to be a Universal system for both TOAD and ASTRA.
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    properties
        % System & Vehicle
        constants                       % Vehicle constants struct from LoadTOADParams
        Vehicle double = 0              % 0 = ASTRA (electric), 1 = TOAD (liquid biprop)

        % Discretization
        N double = 60                   % Number of control intervals
        T_bounds = [8, 35]              % Total duration bounds [s]
        T_initial double = 16           % Initial duration guess [s]

        % Corridor Parameters
        CircleTightness double = 0.5   % Radial tolerance band [m]
        GlideslopeAngle double = 8.0    % Takeoff and landing cone half-angle [deg]
        FunnelCurvature double = 0.005  % Larger quadratic flaring
        TakeoffClearanceAlt double = 0  % Liftoff clearance altitude [m] (0 = auto: 1.0m ASTRA, 2.5m TOAD)
        LandingFlareAlt double = 0      % Terminal descent flare altitude [m] (0 = auto: 1.2m ASTRA, 3.0m TOAD)

        % Multi-Objective Cost Weights
        w_time double   = 0.5          % Mission duration
        w_length double = 0.45          % 3D spatial path length
        w_effort double = 0.15          % Control effort (hover deviation + gimbal tilt)
        w_slew double   = 0.10          % Thrust/roll actuator slew rate
        w_slew_gim double = 0.20        % TVC gimbal slew penalty
        w_curv double   = 0.05          % Thrust/roll actuator curvature
        w_curv_gim double = 0.10        % TVC gimbal curvature penalty
        w_smooth double = 0.02          % Velocity step smoothness
        w_rate double   = 0.025         % Body angular rate penalty
        w_qz double     = 0.025         % Yaw deflection penalty

        % Maneuver Definition
        Maneuver string = "Circle"      % "Circle", "Backflip", "Hop", "Waypoint", "Custom"
        ManeuverParams struct
        Waypoints double = []           % 3 x K matrix of waypoints [m]
        WaypointTolerances double = []  % 1 x K tolerance sphere radii [m]

        % Boundary Conditions
        q0 = [1; 0; 0; 0]
        r0 = [0; 0; 0]
        v0 = [0; 0; 0]
        w0 = [0; 0; 0]
        r_f = [0; 0; 0]
        v_f_tol = 0.35                  % Max allowable touchdown speed [m/s]

        % Control Limits
        thrust_margin = 0.05
        gimbal_margin = 0.15
        max_gimbal_rate = deg2rad(30)   % Max gimbal slew rate [rad/s]
        max_thrust_rate = 1000          % Max thrust rate [N/s]
        max_roll_rate = 4               % Max roll torque [N*m]
        max_gimbal_angle = pi/15        % Max gimbal angle [rad] (12 deg)

        % Scaling & Solver
        Sx double
        Su double
        L_c double = 50                 % Characteristic position length [m]
        MaxIter double = 120
        MaxCpuTime double = 45.0
        Tol double = 1.5e-2
        ConstrViolTol double = 5e-3
        PrintLevel double = 0

        % Storage & File Export
        SaveDir string = ""
        Version double = 1
        Filename string = ""
        Stage1Diag struct = struct()
        QPTraj                          % PolyTrajectoryQP object
        InitialGuess struct = struct()
        Solution struct = struct()
        OptiVars struct = struct()
    end

    properties (Dependent)
        isElectric logical
    end

    methods
        function val = get.isElectric(obj)
            val = (obj.Vehicle == 0);
        end

        %% Object Configuration

        function obj = TrajectoryOptimizerHybrid(constants6DoF, varargin)
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizerHybrid requires constants6DoF struct from LoadTOADParams.');
            end
            obj.constants = constants6DoF;

            % Ensure CasADi toolbox is on the search path
            if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir') && isempty(which('casadi.Opti'))
                addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
            end

            p = inputParser;
            p.KeepUnmatched = true;
            props = {'Vehicle','Maneuver','N','CircleTightness','GlideslopeAngle','FunnelCurvature',...
                     'TakeoffClearanceAlt','LandingFlareAlt',...
                     'w_time','w_length','w_effort','w_slew','w_slew_gim','w_curv','w_curv_gim',...
                     'w_smooth','w_rate','w_qz','PrintLevel','MaxIter','Tol','ConstrViolTol',...
                     'SaveDir','Version','Filename'};
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

        function h = getTakeoffClearanceAlt(obj)
            %% GETTAKEOFFCLEARANCEALT  Liftoff altitude clearance threshold.
            if obj.TakeoffClearanceAlt > 0
                h = obj.TakeoffClearanceAlt;
            else
                h = 1.0 * (obj.Vehicle == 0) + 2.5 * (obj.Vehicle == 1);
            end
        end

        function h = getLandingFlareAlt(obj)
            %% GETLANDINGFLAREALT  Terminal descent landing flare altitude threshold.
            if obj.LandingFlareAlt > 0
                h = obj.LandingFlareAlt;
            else
                h = 1.2 * (obj.Vehicle == 0) + 3.0 * (obj.Vehicle == 1);
            end
        end

        function updateScaling(obj)
            % Updates scaling to O(1) depending on the vehicle selected
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
                if obj.Vehicle == 1
                    t_phys = 2.0 * sqrt(2.0 * dh / 3.0) + max(4.0, sqrt(d_h));
                    obj.T_initial = max(16.0, min(36.0, round(t_phys, 1)));
                    obj.T_bounds = [max(12.0, round(0.65 * obj.T_initial, 1)), ...
                                    min(42.0, round(1.50 * obj.T_initial, 1))];
                else
                    t_phys = 2.0 * sqrt(2.0 * dh / 3.5) + max(3.0, sqrt(d_h));
                    obj.T_initial = max(10.0, min(24.0, round(t_phys, 1)));
                    obj.T_bounds = [max(8.0, round(0.65 * obj.T_initial, 1)), ...
                                    min(28.0, round(1.50 * obj.T_initial, 1))];
                end
            end
        end

        function setManeuver(obj, name, varargin)
            obj.Maneuver = string(name);
            p = inputParser;
            % Parameters for manouvers
            if obj.Vehicle == 0
                dR = 5; dA = 7; dHA = 15; dBA = 20;
            else
                dR = 15; dA = 25; dHA = 50; dBA = 45;
            end

            switch obj.Maneuver
                case "Circle"
                    addParameter(p, 'circle_radius', dR);
                    addParameter(p, 'circle_alt', dA);
                    addParameter(p, 'circle_center', [0; 0]);
                    parse(p, varargin{:});
                    obj.ManeuverParams = struct('circle_radius', double(p.Results.circle_radius), ...
                        'circle_alt', double(p.Results.circle_alt), ...
                        'circle_center', double(p.Results.circle_center(:)), ...
                        'f_orbit_start', 0.25, 'f_orbit_end', 0.75);
                    if obj.Vehicle == 0
                        obj.T_bounds = [12, 26]; obj.T_initial = 18;
                    else
                        obj.T_bounds = [18, 38]; obj.T_initial = 26;
                    end

                case "Backflip"
                    addParameter(p, 'apex_alt', dBA);
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
                    addParameter(p, 'apex_alt', dHA);
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
                error('Waypoints must be a 3 x K matrix with K >= 2.');
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

        %% Stage 1, Polynomial Optimizer

        function [Fix, T_seg, Tag] = buildKeyframes(obj)
            %% BUILDKEYFRAMES  Convert maneuver geometry to PolyTrajectoryQP keyframe arrays.
            % Initial segment times are seeded proportionally from T_initial;
            % QPTraj.optimizeTimes() subsequently computes the optimal segment allocation.
            R = 4;  % minimum-snap (septic polynomials)
            vL = [0; 0; -0.10];  % gentle landing velocity

            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    Rad = mp.circle_radius; h = mp.circle_alt;
                    N_orb = 16;  % polygon vertices around orbit.
                    angles = linspace(0, 2*pi, N_orb + 1);
                    K = 1 + (N_orb + 1) + 1;
                    Fix = nan(3, R, K); Tag = repmat({''}, 1, K);

                    % Launch
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    % Orbit polygon
                    for j = 1:(N_orb + 1)
                        k = 1 + j;
                        Fix(:,1,k) = [cx + Rad*cos(angles(j)); cy + Rad*sin(angles(j)); h];
                    end
                    Tag{2} = 'orb_entry'; Tag{N_orb + 2} = 'orb_exit';
                    % Landing
                    Fix(:,1,K) = obj.r_f; Fix(:,2,K) = vL; Tag{K} = 'land';

                    % Segment times: ascent / orbit / descent
                    M = K - 1; Ts = obj.T_initial;
                    T_seg = zeros(1, M);
                    T_seg(1) = 0.25 * Ts;
                    for j = 1:N_orb, T_seg(1 + j) = 0.50 * Ts / N_orb; end
                    T_seg(M) = 0.25 * Ts;

                case "Backflip"
                    mp = obj.ManeuverParams;
                    apex = [0.5*(obj.r0(1)+obj.r_f(1)); 0.5*(obj.r0(2)+obj.r_f(2)); mp.apex_alt];
                    K = 3; Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    Fix(:,1,2) = apex;  Fix(3,2,2) = 0;  Tag{2} = 'apex';
                    Fix(:,1,3) = obj.r_f; Fix(:,2,3) = vL; Tag{3} = 'land';
                    T_seg = [0.50, 0.50] * obj.T_initial;

                case "Hop"
                    mp = obj.ManeuverParams;
                    apex = [0.5*(obj.r0(1)+obj.r_f(1)); 0.5*(obj.r0(2)+obj.r_f(2)); mp.apex_alt];
                    K = 3; Fix = nan(3, R, K); Tag = repmat({''}, 1, K);
                    Fix(:,1,1) = obj.r0; Fix(:,2,1) = obj.v0; Tag{1} = 'launch';
                    Fix(:,1,2) = apex;  Fix(3,2,2) = 0;  Tag{2} = 'apex';
                    Fix(:,1,3) = obj.r_f; Fix(:,2,3) = vL; Tag{3} = 'land';
                    % Initial symmetric 50/50 seed; QPTraj.optimizeTimes refines segment durations
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
                    
                    % Distance-proportional time allocation
                    M = K - 1; dists = zeros(1, M);
                    for m = 1:M, dists(m) = norm(obj.Waypoints(:,m+1) - obj.Waypoints(:,m)); end
                    td = sum(dists);
                    if td < 1e-3
                        T_seg = (obj.T_initial / M) * ones(1, M);
                    else
                        T_seg = max(1.5, obj.T_initial * dists / td);
                        T_seg = obj.T_initial * T_seg / sum(T_seg);
                    end

                otherwise  % Straight-line default
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

            % Optimize segment durations for multi-segment trajectories (Richter Sec. V)
            if obj.Maneuver ~= "Circle" && numel(obj.QPTraj.T) >= 2
                obj.QPTraj.optimizeTimes(10);
            end

            % Enforce geometric flight corridors via iterative keyframe insertion
            corridors = obj.makeCorridors();
            if ~isempty(corridors)
                obj.QPTraj.enforceCorridors(corridors, 20, 0.05);
            end

            % Universal keyframe tag synchronization to ManeuverParams
            tg = obj.QPTraj.tagTimes();
            Ttot = sum(obj.QPTraj.T);
            if isfield(tg, 'orb_entry')
                obj.ManeuverParams.f_orbit_start = tg.orb_entry / Ttot;
            end
            if isfield(tg, 'orb_exit')
                obj.ManeuverParams.f_orbit_end = tg.orb_exit / Ttot;
            end
            if isfield(tg, 'apex')
                obj.ManeuverParams.f_apex = tg.apex / Ttot;
            end

            obj.Stage1Diag = struct('t_qp_ms', toc(t0) * 1000);
        end

        function corridors = makeCorridors(obj)
            %% MAKECORRIDORS  Universal glideslope and landing-funnel corridor function handles.
            corridors = { ...
                struct('fn', @(ts,rs,tg) obj.evalTakeoffCone(ts,rs,tg)), ...
                struct('fn', @(ts,rs,tg) obj.evalLandingFunnel(ts,rs,tg)) ...
            };
        end

        function [viol, tgt] = evalTakeoffCone(obj, ts, rs, tg)
            %% EVALTAKEOFFCONE  Universal violation function for takeoff glideslope cone.
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts);

            % Liftoff transition time: first waypoint/orbit tag or default 15% duration
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
            %% EVALLANDINGFUNNEL  Universal post-last-checkpoint altitude-triggered landing funnel.
            n = numel(ts); viol = -inf(1, n); tgt = rs;
            t_max = max(ts);

            % Terminal descent begins post last checkpoint
            t_start = 0.70 * t_max;
            if isfield(tg, 'orb_exit')
                t_start = max(0.65 * t_max, tg.orb_exit);
            elseif isfield(tg, 'apex')
                t_start = max(0.60 * t_max, tg.apex);
            elseif ~isempty(obj.Waypoints) && size(obj.Waypoints, 2) >= 2
                K = size(obj.Waypoints, 2);
                tag_penult = sprintf('wp_%d', K-1);
                if isfield(tg, tag_penult)
                    t_start = tg.(tag_penult);
                end
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

        %% Inversion to get inputs

        function guess = buildInitialGuess(obj)
            %% Sample QP trajectory and invert via differential flatness.
            if isempty(obj.QPTraj)
                obj.buildQPTrajectory();
            end
            t0 = tic;
            Nn = obj.N + 1; Nc = obj.N;
            g = obj.constants.g; MT = obj.constants.MaxThrust;
            thr_allow = (1 - obj.thrust_margin) * MT;

            % Dynamic time sizing: coarsely check thrust feasibility and rescale
            if obj.Vehicle == 0
                m_chk = obj.constants.m_dry;
            else
                m_chk = obj.constants.m_dry + obj.constants.OxMass + obj.constants.FuMass;
            end
            for sizing = 1:3
                Tq = sum(obj.QPTraj.T);
                tv_c = linspace(0, Tq, 20);
                a_c = obj.QPTraj.evalAt(tv_c, 2);
                Fmag = m_chk * sqrt(sum((a_c + [0; 0; g]).^2, 1));
                kap = sqrt(max(Fmag) / thr_allow);
                if kap <= 1.03, break; end
                obj.QPTraj.scaleTime(min(kap, 1.30));
            end

            % Sample QP trajectory at NLP resolution
            T_total = sum(obj.QPTraj.T);
            tv = linspace(0, T_total, Nn);
            dt = T_total / Nc; dt_row = dt * ones(1, Nc);
            r_s = obj.QPTraj.evalAt(tv, 0);
            v_s = obj.QPTraj.evalAt(tv, 1);
            a_s = obj.QPTraj.evalAt(tv, 2);
            j_s = obj.QPTraj.evalAt(tv, 3);

            % Exact zero-defect mass profile & control initialization
            m_dry = obj.constants.m_dry;
            m_vec = m_dry * ones(1, Nn);
            ml = zeros(1, Nn); mi = zeros(1, Nn);
            if obj.Vehicle == 1
                ml(1) = obj.constants.OxMass;
                mi(1) = obj.constants.FuMass;
                m_vec(1) = m_dry + ml(1) + mi(1);
            end

            q_s = zeros(4, Nn); w_s = zeros(3, Nn); U_s = zeros(4, Nc);
            Jyy = obj.constants.J(2,2);
            cgz = max(0.15, obj.constants.rTB);
            plim = sin((1 - obj.gimbal_margin) * obj.max_gimbal_angle);
            OF = obj.constants.OF;
            Mdot = obj.constants.MaxMdot;
            t_lo = (0.25 + obj.thrust_margin) * MT;
            t_hi = (1.00 - obj.thrust_margin) * MT;

            if obj.Maneuver == "Backflip"
                % Backflip: position from QP, attitude from smooth cycloidal ramp
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
                end

                if obj.Vehicle == 1
                    for k = 1:Nc
                        mdl_k = -(U_s(3,k)/MT) * (OF/(1+OF)) * Mdot;
                        mdi_k = -(U_s(3,k)/MT) * (1/(1+OF))  * Mdot;
                        ml(k+1) = max(0.02*obj.constants.OxMass, ml(k) + dt*mdl_k);
                        mi(k+1) = max(0.02*obj.constants.FuMass, mi(k) + dt*mdi_k);
                        m_vec(k+1) = m_dry + ml(k+1) + mi(k+1);
                    end
                end
            else
                % Normal flatness inversion: thrust-vector -> quaternion
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
                    q_flat = obj.vectorToQuat(b3);
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
                    R_IB = obj.quatToRot(q_s(:,k));
                    b3d_B = R_IB' * b3_dot;
                    w_s(:,k) = [-b3d_B(2); b3d_B(1); 0];
                end
                q_s(:,1) = obj.q0; q_s(:,end) = obj.q0;
                w_s(:,1) = obj.w0; w_s(:,end) = [0; 0; 0];

                % TVC gimbal synthesis
                J_mat = obj.constants.J; rTB = max(0.15, obj.constants.rTB);
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

            % Universal forward-backward rate smoothing across all control channels
            for ch = 1:2
                for k = 2:Nc
                    max_dg = 0.90 * obj.max_gimbal_rate * dt;
                    U_s(ch,k) = max(U_s(ch,k-1) - max_dg, min(U_s(ch,k-1) + max_dg, U_s(ch,k)));
                end
                for k = Nc-1:-1:1
                    max_dg = 0.90 * obj.max_gimbal_rate * dt;
                    U_s(ch,k) = max(U_s(ch,k+1) - max_dg, min(U_s(ch,k+1) + max_dg, U_s(ch,k)));
                end
            end
            for k = 2:Nc
                max_d = 0.90 * obj.max_thrust_rate * dt;
                U_s(3,k) = max(U_s(3,k-1) - max_d, min(U_s(3,k-1) + max_d, U_s(3,k)));
            end
            for k = Nc-1:-1:1
                max_d = 0.90 * obj.max_thrust_rate * dt;
                U_s(3,k) = max(U_s(3,k+1) - max_d, min(U_s(3,k+1) + max_d, U_s(3,k)));
            end
            for k = 2:Nc
                max_dr = 0.90 * obj.max_roll_rate * dt;
                U_s(4,k) = max(U_s(4,k-1) - max_dr, min(U_s(4,k-1) + max_dr, U_s(4,k)));
            end
            for k = Nc-1:-1:1
                max_dr = 0.90 * obj.max_roll_rate * dt;
                U_s(4,k) = max(U_s(4,k+1) - max_dr, min(U_s(4,k+1) + max_dr, U_s(4,k)));
            end

            % Fix boundary values
            r_s(:,1) = obj.r0; r_s(:,end) = obj.r_f; v_s(:,end) = [0; 0; -0.05];

            % Package initial guess
            Xp = [q_s; r_s; v_s; w_s; ml; mi];
            obj.T_initial = T_total;
            if obj.Maneuver == "Waypoint" || obj.Maneuver == "Custom"
                obj.T_bounds = [round(0.85 * T_total, 1), round(1.20 * T_total, 1)];
            else
                obj.T_bounds = [min(obj.T_bounds(1), max(4.0, 0.50*T_total)), ...
                                max(obj.T_bounds(2), min(60.0, 1.80*T_total))];
            end

            obj.Stage1Diag.T_total = T_total;
            obj.Stage1Diag.max_thrust = max(U_s(3,:));
            obj.Stage1Diag.max_gimbal_deg = rad2deg(max(sqrt(U_s(1,:).^2 + U_s(2,:).^2)));
            obj.Stage1Diag.t_stage1_ms = toc(t0) * 1000;

            guess = struct('Time', tv, 'X', Xp, 'U', U_s, 'T_total', T_total, ...
                           'dt_row', dt_row, 'Xhat', Xp ./ obj.Sx, 'Uhat', U_s ./ obj.Su);
            obj.InitialGuess = guess;
        end

        %% Direct Collocation Stage

        function sol = optimize(obj)
            %% OPTIMIZE  Alias for solve()
            sol = obj.solve();
        end

        function sol = solve(obj)
            %% SOLVE  Two-Stage Pipeline: QP Polynomial Seed -> 6-DoF Collocation Refinement
            import casadi.*
            t_total_start = tic;

            if isempty(obj.QPTraj), obj.buildQPTrajectory(); end
            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'Xhat')
                obj.buildInitialGuess();
            end

            N = obj.N; opti = Opti(); dyn_fnc = obj.getCasADiDynamics();
            obj.OptiVars = struct();

            % Decision variables
            T_total = opti.variable();
            opti.subject_to(obj.T_bounds(1) <= T_total);
            opti.subject_to(T_total <= obj.T_bounds(2));
            opti.set_initial(T_total, obj.InitialGuess.T_total);
            dt_row = repmat(T_total / N, 1, N);
            obj.OptiVars.T_total = T_total; obj.OptiVars.dt_row = dt_row;

            Xhat = opti.variable(15, N + 1); Uhat = opti.variable(4, N);
            X = obj.Sx .* Xhat; U = obj.Su .* Uhat;

            % RK4 Discretization
            params_val = [obj.constants.m_dry; obj.constants.g; obj.constants.rTB; ...
                obj.constants.Ox_Z; obj.constants.OxMass; obj.constants.OxHeight; ...
                obj.constants.Fu_Z; obj.constants.FuMass; obj.constants.FuHeight; ...
                obj.constants.J(:); obj.constants.OxRadius; obj.constants.FuRadius; ...
                obj.constants.MaxThrust; obj.constants.OF; obj.constants.MaxMdot; ...
                0; zeros(9, 1); zeros(3, 1)];

            xs = MX.sym('xh', 15); us = MX.sym('uh', 4); ds = MX.sym('dt');
            xp = obj.Sx .* xs; up = obj.Su .* us;
            k1 = dyn_fnc(xp,               up, params_val);
            k2 = dyn_fnc(xp + ds/2 * k1,   up, params_val);
            k3 = dyn_fnc(xp + ds/2 * k2,   up, params_val);
            k4 = dyn_fnc(xp + ds * k3,     up, params_val);
            xn = xp + ds / 6 * (k1 + 2*k2 + 2*k3 + k4);
            xn = [xn(1:4) / sqrt(sum(xn(1:4).^2) + 1e-12); xn(5:end)];
            Fstep = Function('F_step', {xs, us, ds}, {xn ./ obj.Sx});
            Fmap  = Fstep.map(N);
            opti.subject_to(Xhat(:, 2:end) == Fmap(Xhat(:,1:N), Uhat, dt_row));

            % Universal dynamic domain bounds
            max_h_extent = max([norm(obj.r0(1:2)), norm(obj.r_f(1:2)), norm(obj.r_f(1:2) - obj.r0(1:2)), 30.0]);
            if ~isempty(obj.Waypoints)
                max_h_extent = max(max_h_extent, max(sqrt(sum(obj.Waypoints(1:2,:).^2, 1))));
            end
            R_horiz = max(150.0, 3.0 * max_h_extent);
            opti.subject_to(-R_horiz <= X(5,:)); opti.subject_to(X(5,:) <= R_horiz);
            opti.subject_to(-R_horiz <= X(6,:)); opti.subject_to(X(6,:) <= R_horiz);
            
            % Determine domain altitude ceiling from all active targets
            alts = [obj.r0(3), obj.r_f(3), 20.0];
            if isfield(obj.ManeuverParams, 'apex_alt'),   alts(end+1) = obj.ManeuverParams.apex_alt;   end
            if isfield(obj.ManeuverParams, 'circle_alt'), alts(end+1) = obj.ManeuverParams.circle_alt; end
            if ~isempty(obj.Waypoints),                   alts = [alts, obj.Waypoints(3,:)];          end
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
            MinThrottle = obj.Vehicle * 0.5 + (1 - obj.Vehicle) * 0.25;
            opti.subject_to((MinThrottle+tm)*MT <= U(3,:)); opti.subject_to(U(3,:) <= (1-tm)*MT);
            opti.subject_to(-(1-gm)*mg <= U(1,:));   opti.subject_to(U(1,:) <= (1-gm)*mg);
            opti.subject_to(-(1-gm)*mg <= U(2,:));   opti.subject_to(U(2,:) <= (1-gm)*mg);
            opti.subject_to(-(1-tm)*obj.max_roll_rate <= U(4,:));
            opti.subject_to(U(4,:) <= (1-tm)*obj.max_roll_rate);
            dU = U(:,2:end) - U(:,1:end-1); dts = dt_row(1:end-1);
            for ch = [1, 2]
                opti.subject_to(-obj.max_gimbal_rate * dts <= dU(ch,:));
                opti.subject_to(dU(ch,:) <= obj.max_gimbal_rate * dts);
            end
            opti.subject_to(-obj.max_thrust_rate*dts <= dU(3,:));
            opti.subject_to(dU(3,:) <= obj.max_thrust_rate*dts);
            opti.subject_to(-obj.max_roll_rate*dts <= dU(4,:));
            opti.subject_to(dU(4,:) <= obj.max_roll_rate*dts);

            % Maneuver constraints & cost
            obj.applyManeuverConstraints(opti, X);
            obj.applyCostFunction(opti, X, Xhat, U, Uhat, T_total, dt_row);

            % Initial guess
            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);

            % IPOPT solver
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, 'tol', obj.Tol, ...
                'constr_viol_tol', obj.ConstrViolTol, ...
                'dual_inf_tol', 1.0, ...
                'acceptable_tol', max(obj.Tol, 2e-2), ...
                'acceptable_constr_viol_tol', 1e-3, ...
                'acceptable_dual_inf_tol', 1e5, ...
                'acceptable_iter', 2, ...
                'max_cpu_time', obj.MaxCpuTime, ...
                'mu_strategy', 'adaptive', 'print_level', obj.PrintLevel);
            if obj.Vehicle == 0 && (obj.Maneuver == "Waypoint" || obj.Maneuver == "Custom")
                s_opts.hessian_approximation = 'limited-memory';
            end
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
                it_cnt = obj.MaxIter;
                try
                    st = opti.stats();
                    if isfield(st, 'iter_count'), it_cnt = st.iter_count; end
                catch
                end
                stats = struct('t_wall_total', toc(t_stg2_start), 'iter_count', it_cnt, ...
                    'return_status', 'Infeasible', 'error_msg', ME.message);
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
                         'Time', t_vec, 't', t_vec, 'stats', stats, 'opti_obj', so);
            obj.Solution = sol;
        end

        function applyManeuverConstraints(obj, opti, X)
            %% APPLYMANEUVERCONSTRAINTS  Streamlined, minimal, universal flight constraints.
            N = obj.N;

            % Universal Liftoff & Touchdown Alignments
            N_lo = max(2, round(0.04 * N));
            for ch = 8:9, opti.subject_to(abs(X(ch, 1:N_lo)) <= 0.60); end
            R33_lo = X(1,1:N_lo).^2 - X(2,1:N_lo).^2 - X(3,1:N_lo).^2 + X(4,1:N_lo).^2;
            opti.subject_to(R33_lo >= cosd(25));

            % Dedicated terminal touchdown alignment (last 3 nodes)
            for ch = 8:9, opti.subject_to(abs(X(ch, end-2:end)) <= 0.35); end
            opti.subject_to(X(10, end-2:end) >= -2.5);
            R33_t = X(1,end-2:end).^2 - X(2,end-2:end).^2 - X(3,end-2:end).^2 + X(4,end-2:end).^2;
            opti.subject_to(R33_t >= cosd(20));

            % Universal Altitude-Triggered Takeoff Glideslope Cone
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
            opti.subject_to(rxy_to_sq <= R_to.^2);

            % Universal Altitude-Triggered Flaring Landing Funnel
            h_ld_clear = obj.getLandingFlareAlt();
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                z_g = obj.InitialGuess.X(7, :);
                idx_ld = find(z_g - obj.r_f(3) <= h_ld_clear & (1:N_n) >= round(0.65*N_n));
                if isempty(idx_ld), N_land = N_n - 2; else, N_land = min(N_n - 2, min(idx_ld)); end
            else
                N_land = min(N_n - 2, round(0.85 * N_n));
            end

            rf_p = obj.r_f;
            dz_ld = fmax(0, X(7, N_land:end) - rf_p(3));
            cf_ld = obj.FunnelCurvature;
            R_ld = dz_ld * tand(obj.GlideslopeAngle) + cf_ld * (dz_ld.^2) + 0.15;
            rxy_ld_sq = (X(5, N_land:end) - rf_p(1)).^2 + (X(6, N_land:end) - rf_p(2)).^2;
            opti.subject_to(rxy_ld_sq <= R_ld.^2);

            % Flared descent vertical velocity limit
            v_desc_lim = 3.5 * (obj.Vehicle == 0) + 6.0 * (obj.Vehicle == 1);
            opti.subject_to(X(10, N_land:end) >= -v_desc_lim);

            % Upright attitude alignment in funnel
            R33_ld = X(1, N_land:end).^2 - X(2, N_land:end).^2 - X(3, N_land:end).^2 + X(4, N_land:end).^2;
            opti.subject_to(R33_ld >= cosd(25));

            % Universal upright tilt envelope for non-Backflip maneuvers
            if obj.Maneuver ~= "Backflip"
                R33_all = X(1,:).^2 - X(2,:).^2 - X(3,:).^2 + X(4,:).^2;
                opti.subject_to(R33_all >= cosd(50));
            end

            %%  Maneuver-Specific Formulations
            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    Nos = max(3, round(mp.f_orbit_start * N));
                    Noe = min(N-2, round(mp.f_orbit_end * N));
                    Norb = max(2, Noe - Nos);
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    R = mp.circle_radius; h = mp.circle_alt;
                    tr = max(0.35, obj.CircleTightness);

                    % Monotonic climb to orbit, monotonic descent from orbit
                    opti.subject_to(X(10, 1:Nos) >= 0.0);
                    opti.subject_to(X(10, Noe:end) <= 0.0);

                    % Radial annulus + altitude band
                    rsq = (X(5,Nos:Noe)-cx).^2 + (X(6,Nos:Noe)-cy).^2;
                    opti.subject_to((R-tr)^2 <= rsq); opti.subject_to(rsq <= (R+tr)^2);
                    opti.subject_to(abs(X(7,Nos:Noe) - h) <= 0.35);

                    % Continuous CCW angular progression
                    rx = X(5,Nos:Noe)-cx; ry = X(6,Nos:Noe)-cy;
                    cprog = rx(1:end-1) .* ry(2:end) - ry(1:end-1) .* rx(2:end);
                    opti.subject_to(cprog >= 0.0);

                    % Quadrant progression checkpoints
                    kq1 = Nos + round(0.25*Norb); kq2 = Nos + round(0.50*Norb); kq3 = Nos + round(0.75*Norb);
                    opti.subject_to(X(6, kq1) - cy >= 0);
                    opti.subject_to(X(5, kq2) - cx <= 0);
                    opti.subject_to(X(6, kq3) - cy <= 0);
                    opti.subject_to((X(5, Noe) - (cx + R))^2 + (X(6, Noe) - cy)^2 <= tr^2);

                case "Backflip"
                    mp = obj.ManeuverParams;
                    Na  = max(3, round(mp.flip_start_frac * N)) + 1;
                    Nf  = round(0.5*(mp.flip_start_frac + mp.flip_end_frac) * N) + 1;
                    Nap = min(N, round(mp.flip_end_frac * N)) + 1;
                    att_tol = cos(mp.theta_tol / 2);

                    opti.subject_to(X(10, 1:Na) >= 0.0);
                    opti.subject_to(X(10, Nap:end) <= 0.0);
                    opti.subject_to(X(7,:) <= mp.apex_alt + 3.5);
                    opti.subject_to(X(7,Nf) >= mp.apex_alt - 3.5);
                    opti.subject_to(mp.q_inverted' * X(1:4,Nf) >= att_tol);

                    % Unidirectional pitch rate during flip: guarantees full 360-degree rotation (no pendulum rocking)
                    opti.subject_to(X(12, Na:Nap) <= -0.05);

                    % Universal planar constraints for both vehicles (suppress lateral sway and yaw/roll cross-coupling)
                    opti.subject_to(-0.75 <= X(6,:));  opti.subject_to(X(6,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(9,:));  opti.subject_to(X(9,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(11,:)); opti.subject_to(X(11,:) <= 0.75);
                    opti.subject_to(-0.75 <= X(13,:)); opti.subject_to(X(13,:) <= 0.75);

                case "Hop"
                    mp = obj.ManeuverParams;
                    if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                        [~, Nap_h] = max(obj.InitialGuess.X(7,:));
                    else
                        Nap_h = round(0.50 * N);
                    end
                    Nap_h = max(round(0.20*N), min(round(0.80*N), Nap_h));

                    opti.subject_to(X(7,Nap_h) >= mp.apex_alt - 2.5);
                    opti.subject_to(X(7,:) <= mp.apex_alt + 3.5);
                    opti.subject_to(X(10, 1:round(0.15*N)) >= 0.0);
                    opti.subject_to(X(10, round(0.85*N):end) <= 0.0);

                case {"Waypoint", "Custom"}
                    % Heading deflection clamping: keep nominal survey heading (+/- 28 deg)
                    % Eliminates TVC horizontal spin / corkscrews in low-inertia vehicles
                    opti.subject_to(-0.25 <= X(4,:)); opti.subject_to(X(4,:) <= 0.25);
                    opti.subject_to(-2.5 <= X(11,:)); opti.subject_to(X(11,:) <= 2.5);
                    opti.subject_to(-2.5 <= X(12,:)); opti.subject_to(X(12,:) <= 2.5);
                    opti.subject_to(-1.0 <= X(13,:)); opti.subject_to(X(13,:) <= 1.0);

                    if ~isempty(obj.Waypoints)
                        K = size(obj.Waypoints, 2);
                        Ttot = sum(obj.QPTraj.T);
                        tg = obj.QPTraj.tagTimes();
                        for m = 2:(K-1)
                            tag_name = sprintf('wp_%d', m);
                            if isfield(tg, tag_name)
                                t_wp = tg.(tag_name);
                            else
                                t_wp = (m - 1) * Ttot / (K - 1);
                            end
                            k_wp = max(2, min(N, round(N * t_wp / Ttot) + 1));
                            tol_m = 0.80;
                            if ~isempty(obj.WaypointTolerances) && numel(obj.WaypointTolerances) >= m
                                tol_m = max(0.80, obj.WaypointTolerances(m));
                            end
                            opti.subject_to(sum((X(5:7, k_wp) - obj.Waypoints(:, m)).^2) <= tol_m^2);
                        end
                    end
            end
        end

        function applyCostFunction(obj, opti, X, Xhat, U, Uhat, T_total, dt_row)
            %% APPLYCOSTFUNCTION  Unified multi-objective cost with chatter suppression.
            N = obj.N; g = obj.constants.g; MT = obj.constants.MaxThrust;

            % Time optimality
            J_time = T_total / obj.T_initial;

            % 3D spatial path length
            dr = X(5:7, 2:end) - X(5:7, 1:end-1);
            L_path = sum(sqrt(sum(dr.^2, 1) + 1e-6));
            % Universal path length normalization using Stage-1 QP trajectory length
            if ~isempty(obj.InitialGuess) && isfield(obj.InitialGuess, 'X')
                dr_guess = obj.InitialGuess.X(5:7, 2:end) - obj.InitialGuess.X(5:7, 1:end-1);
                L_ref = max(1.0, sum(sqrt(sum(dr_guess.^2, 1))));
            else
                L_ref = max(1.0, norm(obj.r_f - obj.r0) + 20.0);
            end
            J_length = L_path / L_ref;

            % Control effort
            if obj.Vehicle == 0
                mv = obj.constants.m_dry * ones(1, N);
            else
                mv = obj.constants.m_dry + X(14,1:N) + X(15,1:N);
            end
            uh = mv * g / MT;
            gb = max(1e-3, (1-obj.gimbal_margin) * obj.max_gimbal_angle);
            esq = (U(3,:)/MT - uh).^2 + (U(1,:)/gb).^2 + (U(2,:)/gb).^2 + ...
                  (U(4,:) / max(1e-3, (1-obj.thrust_margin)*obj.max_roll_rate)).^2;
            J_effort = sum(dt_row .* esq) / (0.50 * obj.T_initial);

            % Actuator slew (dedicated gimbal vs thrust/roll)
            dUh = Uhat(:,2:end) - Uhat(:,1:end-1);
            J_slew_gim  = sum(sum(dUh(1:2,:).^2, 1)) / (N-1);
            J_slew_thr  = sum(dUh(3,:).^2) / (N-1);
            J_slew_roll = sum(dUh(4,:).^2) / (N-1);

            % Actuator curvature (kills chatter)
            d2Uh = Uhat(:,3:end) - 2*Uhat(:,2:end-1) + Uhat(:,1:end-2);
            J_curv_gim  = sum(sum(d2Uh(1:2,:).^2, 1)) / max(1, N-2);
            J_curv_thr  = sum(d2Uh(3,:).^2) / max(1, N-2);
            J_curv_roll = sum(d2Uh(4,:).^2) / max(1, N-2);

            % Velocity step smoothness
            dVh = Xhat(8:10, 2:end) - Xhat(8:10, 1:end-1);
            J_smooth = sum(sum(dVh.^2, 1)) / (N-1);

            % Angular rates & yaw deflection
            J_rate = sum(sum(Xhat(11:13, 1:end-1).^2, 1)) / N;
            J_qz   = sum(Xhat(4,:).^2) / (N+1);

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

        %% ======================== DYNAMICS & KINEMATICS ========================

        function dyn_fnc = getCasADiDynamics(obj)
            import casadi.*
            q = MX.sym('q', 4); r = MX.sym('r', 3); v = MX.sym('v', 3);
            omegaB = MX.sym('omegaB', 3); m_lox = MX.sym('m_lox'); m_ipa = MX.sym('m_ipa');
            m_dry = MX.sym('m_dry'); gSym = MX.sym('g'); rTB = MX.sym('rTB');
            Ox_Z = MX.sym('Ox_Z'); OxMassI = MX.sym('OxMassI'); OxHeight = MX.sym('OxHeight');
            Fu_Z = MX.sym('Fu_Z'); FuMassI = MX.sym('FuMassI'); FuHeight = MX.sym('FuHeight');
            J = MX.sym('J', 3, 3); OxRadius = MX.sym('OxR'); FuRadius = MX.sym('FuR');
            MaxThrust = MX.sym('MT'); OF = MX.sym('OF'); MaxMdot = MX.sym('Mdot');
            MaxMdot_d = MX.sym('Mdot_d'); J_d = MX.sym('Jd', 3, 3); TB_d = MX.sym('TBd', 3);
            theta = MX.sym('th'); phi = MX.sym('ph'); thrust = MX.sym('T'); roll = MX.sym('rl');

            m = m_dry + m_lox + m_ipa;
            C_IB = quatRot(q).';
            TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            FI = C_IB * TB + [0; 0; -m * gSym];

            if obj.Vehicle == 0
                mdl = 0; mdi = 0; OFH = 0; FFH = 0;
            else
                mdl = -thrust/MaxThrust * OF/(1+OF) * (MaxMdot + MaxMdot_d);
                mdi = -thrust/MaxThrust * 1/(1+OF)  * (MaxMdot + MaxMdot_d);
                OFH = (m_lox / OxMassI) * OxHeight * 0.9;
                FFH = (m_ipa / FuMassI) * FuHeight * 0.9;
            end

            Jlox = diag([1/12*m_lox*(3*OxRadius^2+OFH^2), ...
                         1/12*m_lox*(3*OxRadius^2+OFH^2), 1/2*m_lox*OxRadius^2]);
            Jipa = diag([1/12*m_ipa*(3*FuRadius^2+FFH^2), ...
                         1/12*m_ipa*(3*FuRadius^2+FFH^2), 1/2*m_ipa*FuRadius^2]);
            OFL = Ox_Z + OFH/2; FFL = Fu_Z + FFH/2;
            CGz = (m_dry*rTB + m_lox*OFL + m_ipa*FFL) / m;
            dd = rTB - CGz + TB_d(3); dl = OFL - CGz; di = FFL - CGz;
            J_tot = J + m_dry*diag([dd^2,dd^2,0]) + Jlox + m_lox*diag([dl^2,dl^2,0]) ...
                    + Jipa + m_ipa*diag([di^2,di^2,0]) + J_d;

            tDir = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            if obj.Vehicle == 0
                MB = zetaCross([0; 0; -CGz] + TB_d) * TB + roll * tDir;
            else
                MB = zetaCross([0; 0; -CGz] + TB_d) * TB + [0; 0; roll];
            end

            qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
            wdot = (MB - zetaCross(omegaB) * J_tot * omegaB) ./ diag(J_tot);

            x_sym = [q; r; v; omegaB; m_lox; m_ipa];
            u_sym = [theta; phi; thrust; roll];
            params = [m_dry; gSym; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight; ...
                      J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot; MaxMdot_d; J_d(:); TB_d];
            xdot = [qdot; v; FI/m; wdot; mdl; mdi];
            dyn_fnc = Function('dyn_fnc', {x_sym, u_sym, params}, {xdot});
        end

        %% ======================== HELPERS, PLOTTING & EXPORT ========================

        function qout = vectorToQuat(~, v_in)
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

        function R = quatToRot(~, q_in)
            w = q_in(1); x = q_in(2); y = q_in(3); z = q_in(4);
            R = [1-2*(y^2+z^2), 2*(x*y-w*z), 2*(x*z+w*y); ...
                 2*(x*y+w*z), 1-2*(x^2+z^2), 2*(y*z-w*x); ...
                 2*(x*z-w*y), 2*(y*z+w*x), 1-2*(x^2+y^2)];
        end

        function name = getVehicleName(obj)
            if obj.Vehicle == 1
                name = "TOAD";
            else
                name = "ASTRA";
            end
        end

        function fn = getFormattedFilename(obj)
            if obj.Filename ~= ""
                fn = obj.Filename;
            else
                fn = sprintf('%s_%s_v%03d', obj.getVehicleName(), obj.Maneuver, obj.Version);
            end
        end

        function figHandles = plot(obj)
            %% PLOT  Generate publication-quality 3D trajectory and control dashboard.
            % Displays 3D trajectory with takeoff glideslope cone, flaring landing funnel,
            % survey waypoints, apex marker, TVC gimbal deflections, and commanded thrust.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                error('TrajectoryOptimizerHybrid:NoSolution', 'No solution available to plot. Call solve() first.');
            end

            tv = obj.Solution.Time;
            Xs = obj.Solution.X;
            Us = obj.Solution.U;

            f1 = figure('Name', sprintf('%s 6-DoF Trajectory (%s v%03d)', ...
                        obj.Maneuver, obj.getVehicleName(), obj.Version), ...
                        'Position', [100, 100, 1300, 800], 'Visible', 'on', 'Color', 'w');

            % -------------------------------------------------------------
            % 1. 3D Flight Profile with Glideslope Cone & Landing Funnel
            % -------------------------------------------------------------
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
            p_launch = plot3(ax3d, obj.r0(1), obj.r0(2), obj.r0(3), 'o', 'MarkerSize', 9, ...
                             'MarkerFaceColor', [0.15 0.75 0.20], 'MarkerEdgeColor', [0.05 0.40 0.10], 'LineWidth', 1.8);
            p_land   = plot3(ax3d, obj.r_f(1), obj.r_f(2), obj.r_f(3), 's', 'MarkerSize', 9, ...
                             'MarkerFaceColor', [0.85 0.15 0.15], 'MarkerEdgeColor', [0.50 0.05 0.05], 'LineWidth', 1.8);

            leg_handles = [h_traj, p_launch, p_land];
            leg_labels  = {'6-DoF Trajectory', 'Launch Pad', 'Landing Target'};

            % Takeoff Glideslope Cone wireframe (GlideslopeAngle half-angle)
            z_max_traj = max(Xs(7,:));
            z_to_max = obj.getTakeoffClearanceAlt();
            z_cone = linspace(0, z_to_max, 22);
            th_c = linspace(0, 2*pi, 36);
            [TH_to, ZC_to] = meshgrid(th_c, z_cone);
            R_to = ZC_to * tand(obj.GlideslopeAngle) + 0.15;
            XC_to = obj.r0(1) + R_to .* cos(TH_to);
            YC_to = obj.r0(2) + R_to .* sin(TH_to);
            p_to = mesh(ax3d, XC_to, YC_to, ZC_to, 'FaceColor', [0.15 0.75 0.20], 'FaceAlpha', 0.05, ...
                        'EdgeColor', [0.15 0.70 0.20], 'EdgeAlpha', 0.35, 'LineStyle', ':');
            leg_handles(end+1) = p_to;
            leg_labels{end+1}  = sprintf('Takeoff Cone (%.0f°)', obj.GlideslopeAngle);

            % Flaring Landing Funnel wireframe (quadratic flare)
            z_ld_max = obj.getLandingFlareAlt();
            cf_plot = obj.FunnelCurvature;
            z_fun = linspace(0, z_ld_max, 25);
            [TH_ld, ZC_ld] = meshgrid(th_c, z_fun);
            gsa = obj.GlideslopeAngle;
            R_ld = ZC_ld * tand(gsa) + cf_plot * (ZC_ld.^2) + 0.15;
            XC_ld = obj.r_f(1) + R_ld .* cos(TH_ld);
            YC_ld = obj.r_f(2) + R_ld .* sin(TH_ld);
            p_ld = mesh(ax3d, XC_ld, YC_ld, ZC_ld, 'FaceColor', [0.85 0.15 0.15], 'FaceAlpha', 0.05, ...
                        'EdgeColor', [0.85 0.15 0.15], 'EdgeAlpha', 0.35, 'LineStyle', ':');
            leg_handles(end+1) = p_ld;
            leg_labels{end+1}  = sprintf('Landing Funnel (%.0f°)', gsa);

            % Maneuver overlays: Circle orbit or Waypoints
            if obj.Maneuver == "Circle"
                thc = linspace(0, 2*pi, 100);
                mp = obj.ManeuverParams;
                p_orb = plot3(ax3d, mp.circle_center(1) + mp.circle_radius * cos(thc), ...
                              mp.circle_center(2) + mp.circle_radius * sin(thc), ...
                              mp.circle_alt * ones(size(thc)), 'k--', 'LineWidth', 1.5);
                leg_handles(end+1) = p_orb;
                leg_labels{end+1}  = 'Nominal Orbit';
            elseif (obj.Maneuver == "Waypoint" || obj.Maneuver == "Custom") && ~isempty(obj.Waypoints)
                p_wp = plot3(ax3d, obj.Waypoints(1,:), obj.Waypoints(2,:), obj.Waypoints(3,:), 'd--', ...
                             'Color', [0.90 0.50 0.10], 'LineWidth', 1.6, 'MarkerSize', 7, ...
                             'MarkerFaceColor', [1.0 0.75 0.20], 'MarkerEdgeColor', [0.65 0.30 0.05]);
                for kw = 1:size(obj.Waypoints, 2)
                    text(ax3d, obj.Waypoints(1,kw), obj.Waypoints(2,kw), obj.Waypoints(3,kw) + 0.04*z_max_traj, ...
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

            % 1:1 physical aspect ratio with square ground footprint and non-negative altitude
            axis(ax3d, 'equal');

            % Determine square horizontal base encompassing all 3D flight elements
            x_pts = [Xs(5,:), obj.r0(1), obj.r_f(1), XC_to(:)', XC_ld(:)'];
            y_pts = [Xs(6,:), obj.r0(2), obj.r_f(2), YC_to(:)', YC_ld(:)'];
            if ~isempty(obj.Waypoints)
                x_pts = [x_pts, obj.Waypoints(1,:)];
                y_pts = [y_pts, obj.Waypoints(2,:)];
            end
            if obj.Maneuver == "Circle" && isfield(obj.ManeuverParams, 'circle_center')
                mp = obj.ManeuverParams;
                x_pts = [x_pts, mp.circle_center(1) - mp.circle_radius, mp.circle_center(1) + mp.circle_radius];
                y_pts = [y_pts, mp.circle_center(2) - mp.circle_radius, mp.circle_center(2) + mp.circle_radius];
            end

            x_min = min(x_pts); x_max = max(x_pts);
            y_min = min(y_pts); y_max = max(y_pts);
            span_x = x_max - x_min;
            span_y = y_max - y_min;
            base_span = max([span_x, span_y, 4.0]);
            half_w = 0.58 * base_span;

            x_mid = 0.5 * (x_min + x_max);
            y_mid = 0.5 * (y_min + y_max);
            xlim(ax3d, [x_mid - half_w, x_mid + half_w]);
            ylim(ax3d, [y_mid - half_w, y_mid + half_w]);

            % Strict non-negative altitude ceiling: ground starts at z = 0, no negative altitudes
            z_top = max([z_max_traj * 1.08, z_ld_max * 1.15, 5.0]);
            zlim(ax3d, [0, z_top]);
            xlabel(ax3d, 'X (North) [m]', 'FontWeight', 'bold');
            ylabel(ax3d, 'Y (East) [m]', 'FontWeight', 'bold');
            zlabel(ax3d, 'Z (Altitude) [m]', 'FontWeight', 'bold');
            title(ax3d, sprintf('%s Flight Profile (%s, Duration: %.2f s, Path: %.1f m)', ...
                obj.Maneuver, obj.getVehicleName(), obj.Solution.T_total, obj.Solution.L_path), ...
                'FontWeight', 'bold');

            % Maneuver-matched camera view
            switch obj.Maneuver
                case "Circle",    view(ax3d, -35, 26);
                case "Hop",       view(ax3d, -25, 22);
                case "Backflip",  view(ax3d, -30, 20);
                case "Waypoint",  view(ax3d, -32, 28);
                otherwise,        view(ax3d, -35, 25);
            end

            % -------------------------------------------------------------
            % 2. TVC Gimbal Deflections
            % -------------------------------------------------------------
            subplot(2, 2, 2);
            plot(tv(1:end-1), rad2deg(Us(1,:)), 'r-', 'LineWidth', 1.6); hold on;
            plot(tv(1:end-1), rad2deg(Us(2,:)), 'b-', 'LineWidth', 1.6);
            yline(rad2deg(obj.max_gimbal_angle), 'k:', 'Max Gimbal', 'LabelHorizontalAlignment', 'left');
            yline(-rad2deg(obj.max_gimbal_angle), 'k:', '-Max Gimbal', 'LabelHorizontalAlignment', 'left');
            grid on; xlabel('Time [s]', 'FontWeight', 'bold');
            ylabel('Gimbal Angle [deg]', 'FontWeight', 'bold');
            legend({'\theta (Pitch)', '\phi (Yaw)'}, 'Location', 'best');
            title('TVC Gimbal Deflections', 'FontWeight', 'bold');

            % -------------------------------------------------------------
            % 3. Commanded Thrust vs Time
            % -------------------------------------------------------------
            subplot(2, 2, 4);
            plot(tv(1:end-1), Us(3,:), 'k-', 'LineWidth', 1.8); hold on;
            hover_thrust = obj.constants.m_dry * obj.constants.g;
            yline(hover_thrust, 'r--', sprintf('Dry Hover (%.1f N)', hover_thrust), 'LabelHorizontalAlignment', 'left');
            yline(obj.constants.MaxThrust, 'b:', sprintf('Max Thrust (%.1f N)', obj.constants.MaxThrust), 'LabelHorizontalAlignment', 'left');
            grid on; xlabel('Time [s]', 'FontWeight', 'bold');
            ylabel('Thrust [N]', 'FontWeight', 'bold');
            title(sprintf('Commanded Thrust Profile (Status: %s)', obj.Solution.Status), 'FontWeight', 'bold');
            legend({'Thrust Cmd', 'Hover Ref', 'Max Limit'}, 'Location', 'best');

            figHandles = f1;
        end

        function figHandle = plotInitialGuess(obj)
            %% PLOTINITIALGUESS  Visualize analytical initial guess.
            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'X')
                obj.buildInitialGuess();
            end
            tv = obj.InitialGuess.Time; Xg = obj.InitialGuess.X; Ug = obj.InitialGuess.U;

            figHandle = figure('Name', sprintf('%s Initial Guess (%s)', obj.Maneuver, obj.getVehicleName()), ...
                               'Position', [150, 150, 950, 650], 'Visible', 'on', 'Color', 'w');
            ax_g3d = subplot(2, 1, 1);
            set(ax_g3d, 'Tag', 'InitialGuess3D');
            hold(ax_g3d, 'on'); grid(ax_g3d, 'on'); box(ax_g3d, 'on');
            plot3(ax_g3d, Xg(5,:), Xg(6,:), Xg(7,:), 'm--', 'LineWidth', 2);
            plot3(ax_g3d, obj.r0(1), obj.r0(2), obj.r0(3), 'go', 'MarkerSize', 8, 'LineWidth', 2);
            plot3(ax_g3d, obj.r_f(1), obj.r_f(2), obj.r_f(3), 'rs', 'MarkerSize', 8, 'LineWidth', 2);

            axis(ax_g3d, 'equal');
            x_pts_g = [Xg(5,:), obj.r0(1), obj.r_f(1)];
            y_pts_g = [Xg(6,:), obj.r0(2), obj.r_f(2)];
            x_min_g = min(x_pts_g); x_max_g = max(x_pts_g);
            y_min_g = min(y_pts_g); y_max_g = max(y_pts_g);
            base_g = max([x_max_g - x_min_g, y_max_g - y_min_g, 4.0]);
            half_g = 0.58 * base_g;
            x_mid_g = 0.5 * (x_min_g + x_max_g);
            y_mid_g = 0.5 * (y_min_g + y_max_g);
            xlim(ax_g3d, [x_mid_g - half_g, x_mid_g + half_g]);
            ylim(ax_g3d, [y_mid_g - half_g, y_mid_g + half_g]);
            zlim(ax_g3d, [0, max(max(Xg(7,:)) * 1.08, 5.0)]);
            xlabel('X [m]', 'FontWeight', 'bold');
            ylabel('Y [m]', 'FontWeight', 'bold');
            zlabel('Z [m]', 'FontWeight', 'bold');
            title(sprintf('%s Stage-1 QP Seed 3D Trajectory (%s)', obj.Maneuver, obj.getVehicleName()), 'FontWeight', 'bold');
            legend({'QP Seed', 'Launch', 'Landing'}, 'Location', 'best');

            switch obj.Maneuver
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

        function [tbl, filepath] = exportCSV(obj, filepath)
            %% EXPORTCSV  Export trajectory table to CSV matching TOAD format.
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
                default_name = obj.getFormattedFilename() + ".csv";
                filepath = fullfile(target_dir, default_name);
            end

            ts = obj.Solution.Time(:);
            Xr = obj.Solution.X;
            Ur = obj.Solution.U;

            Q  = Xr(1:4, :)';
            P  = Xr(5:7, :)';
            V  = Xr(8:10, :)';
            W  = Xr(11:13, :)';
            ml = Xr(14, :)';
            mf = Xr(15, :)';

            th = [Ur(1, :), Ur(1, end)]';
            ph = [Ur(2, :), Ur(2, end)]';
            T  = [Ur(3, :), Ur(3, end)]';
            ro = [Ur(4, :), Ur(4, end)]';

            tbl = table(ts, Q(:,1), Q(:,2), Q(:,3), Q(:,4), ...
                P(:,1), P(:,2), P(:,3), V(:,1), V(:,2), V(:,3), ...
                W(:,1), W(:,2), W(:,3), ml, mf, th, ph, T, ro, ...
                'VariableNames', {'Time', 'QuatW', 'QuatX', 'QuatY', 'QuatZ', ...
                                  'PosX', 'PosY', 'PosZ', 'VelX', 'VelY', 'VelZ', ...
                                  'AngRateX', 'AngRateY', 'AngRateZ', 'MassLox', 'MassFuel', ...
                                  'GimbalTheta', 'GimbalPhi', 'ThrustMag', 'RollCmd'});
            writetable(tbl, filepath);
            fprintf('  [TrajectoryOptimizerHybrid] Exported trajectory to: %s\n', filepath);
        end
    end
end
