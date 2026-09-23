classdef TrajectoryOptimizer < handle
    %% TrajectoryOptimizer  Unified 6-DoF optimal trajectory generation engine.
    % Formulates and solves vehicle trajectories using CasADi and IPOPT:
    %  1. Angle sweep formulation with bilinear cross-product progression
    %     and tight radial clamping (replaces non-convex annular corridors).
    %  2. Maneuver switch spots parameterized as decision variables.
    %  3. Minimal-time / minimal-path-length Pareto cost with slew regularization.
    %  4. Self-referencing initial guess via analytical quintic splines and
    %     differential flatness with zero-gimbal initialization.
    %
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lullo, & Antigravity)

    properties
        % System & Vehicle Properties
        constants           % Vehicle constants struct (from LoadTOADParams)
        Vehicle double = 0  % 1 for TOAD (liquid biprop), 0 for ASTRA (electric)
    end

    properties (Dependent)
        isElectric logical
    end

    properties
        % Discretization & Mesh
        N double = 120              % Number of control intervals
        T_bounds = [15, 50]         % Total duration bounds [s]
        T_initial double = 25       % Initial total duration guess [s]
        TimeParamMode string = "MultiPhase"  % "MultiPhase" (switch spots) or "Single"
        SwitchSpots double = []     % Solved transition switch times [s]

        % Formulation Settings
        CircleMode string = "AngleSweep"      % "AngleSweep" or "Corridor"
        CircleTightness double = 0.10         % Radial tolerance band [m]
        InitialGuessMode string = "FlatnessDynamic"

        % Multi-Objective Cost Function Weights (Non-Dimensionalized)
        w_time double   = 1.00      % Minimal mission duration weight
        w_length double = 0.50      % Minimal 3D spatial path length weight
        w_effort double = 0.10      % Control effort weight (hover deviation)
        w_slew double   = 0.05      % Actuator slew rate regularization (jerk)
        w_rate double   = 0.050     % Body angular rate penalty weight
        w_qz double     = 0.050     % Yaw deflection penalty weight

        % Maneuver Definition
        Maneuver string = "Circle"  % "Circle", "Backflip", "Hop", "Custom"
        ManeuverParams struct
        CustomWaypoints = []

        % Boundary Conditions
        q0 = [1; 0; 0; 0]           % Initial attitude quaternion [w; x; y; z]
        r0 = [0; 0; 0]              % Launch pad position [m]
        v0 = [0; 0; 0]              % Launch pad velocity [m/s]
        w0 = [0; 0; 0]              % Initial angular rate [rad/s]
        r_f = [0; 0; 0]             % Landing zone position [m]
        v_f_tol = 0.1               % Max allowable touchdown speed [m/s]

        % Control & Rate Limits
        thrust_margin = 0.05
        gimbal_margin = 0.15
        max_gimbal_rate = deg2rad(30)   % Gimbal slew limit [rad/s]
        max_thrust_rate = 1000          % Thrust slew limit [N/s]
        max_roll_rate = 4               % Roll torque limit [N*m]
        max_gimbal_angle = pi/15        % Max physical gimbal angle [rad] (12 deg)

        % Scaling
        Sx double
        Su double
        L_c double = 50                 % Characteristic position length [m]

        % Solver Configuration
        MaxIter double = 1500
        Tol double = 1e-3
        ConstrViolTol double = 1e-3
        PrintLevel double = 0

        % Output & Export Settings
        PlotResults logical = false
        AutoExport logical = false
        Version double = 1              % Version integer (formats as v001, v002, etc.)
        Filename string = ""            % Custom filename override
        SaveDir string = ""             % Target directory for CSV export

        % Solution Storage
        InitialGuess struct
        Solution struct
        OptiVars struct
    end

    methods
        function set.Vehicle(obj, val)
            if ischar(val) || isstring(val)
                strVal = string(val);
                if strcmpi(strVal, "ASTRA") || strcmpi(strVal, "ASTRAv2")
                    obj.Vehicle = 0;
                elseif strcmpi(strVal, "TOAD")
                    obj.Vehicle = 1;
                else
                    error('TrajectoryOptimizer:UnknownVehicle', ...
                        'Unknown vehicle: %s. Must be "TOAD" or "ASTRA".', strVal);
                end
            elseif islogical(val) || isnumeric(val)
                assert(val == 0 || val == 1, 'Vehicle flag must be 0 (ASTRA) or 1 (TOAD).');
                obj.Vehicle = double(val);
            else
                error('TrajectoryOptimizer:InvalidVehicle', 'Unsupported vehicle type: %s', class(val));
            end
        end

        function val = get.isElectric(obj), val = (obj.Vehicle == 0); end

        function obj = TrajectoryOptimizer(constants6DoF, varargin)
            %% TrajectoryOptimizer  Construct optimizer instance from vehicle constants.
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizer requires constants6DoF struct from LoadTOADParams.');
            end

            % Project and tool path resolution
            if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir')
                addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
            end
            this_file = mfilename('fullpath');
            if ~isempty(this_file)
                this_dir = fileparts(this_file);
                proj_root = this_dir;
                while ~isempty(proj_root) && ~exist(fullfile(proj_root, 'LoadTOADParams.m'), 'file')
                    parent = fileparts(proj_root);
                    if strcmp(parent, proj_root), break; end
                    proj_root = parent;
                end
                if exist(fullfile(proj_root, 'LoadTOADParams.m'), 'file')
                    addpath(proj_root);
                    addpath(this_dir, '-begin');
                    addpath(fullfile(proj_root, 'Helper'));
                    addpath(fullfile(proj_root, 'Controls', 'Time Varying LQI'));
                end
            end
            obj.constants = constants6DoF;
            obj.SaveDir = fullfile(pwd, 'Guidance', 'Trajectories');

            p = inputParser; p.KeepUnmatched = true;
            addParameter(p, 'Vehicle', [], @(x) ischar(x)||isstring(x)||isnumeric(x)||islogical(x));
            addParameter(p, 'Maneuver', obj.Maneuver, @(x) ischar(x)||isstring(x));
            addParameter(p, 'Version', obj.Version, @isnumeric);
            addParameter(p, 'N', obj.N, @isnumeric);
            addParameter(p, 'T_initial', obj.T_initial, @isnumeric);
            addParameter(p, 'T_bounds', obj.T_bounds, @isnumeric);
            addParameter(p, 'TimeParamMode', obj.TimeParamMode, @(x) ischar(x)||isstring(x));
            addParameter(p, 'CircleMode', obj.CircleMode, @(x) ischar(x)||isstring(x));
            addParameter(p, 'CircleTightness', obj.CircleTightness, @isnumeric);
            addParameter(p, 'w_time', obj.w_time, @isnumeric);
            addParameter(p, 'w_length', obj.w_length, @isnumeric);
            addParameter(p, 'w_effort', obj.w_effort, @isnumeric);
            addParameter(p, 'w_slew', obj.w_slew, @isnumeric);
            addParameter(p, 'w_rate', obj.w_rate, @isnumeric);
            addParameter(p, 'w_qz', obj.w_qz, @isnumeric);
            addParameter(p, 'PrintLevel', obj.PrintLevel, @isnumeric);
            addParameter(p, 'PlotResults', obj.PlotResults, @islogical);
            addParameter(p, 'AutoExport', obj.AutoExport, @islogical);
            addParameter(p, 'SaveDir', obj.SaveDir, @(x) ischar(x)||isstring(x));
            addParameter(p, 'Filename', obj.Filename, @(x) ischar(x)||isstring(x));
            addParameter(p, 'MaxIter', obj.MaxIter, @isnumeric);
            addParameter(p, 'Tol', obj.Tol, @isnumeric);
            addParameter(p, 'ConstrViolTol', obj.ConstrViolTol, @isnumeric);
            parse(p, varargin{:});

            if ~isempty(p.Results.Vehicle)
                obj.Vehicle = p.Results.Vehicle;
            elseif isfield(constants6DoF, 'Vehicle')
                cV = constants6DoF.Vehicle;
                if ischar(cV) || isstring(cV)
                    obj.Vehicle = double(~strcmpi(string(cV), "ASTRA"));
                else
                    obj.Vehicle = double(cV);
                end
            end

            flds = {'Maneuver','Version','N','T_initial','T_bounds','TimeParamMode', ...
                    'CircleMode','CircleTightness','w_time','w_length','w_effort', ...
                    'w_slew','w_rate','w_qz','PrintLevel','PlotResults', ...
                    'AutoExport','SaveDir','Filename','MaxIter','Tol','ConstrViolTol'};
            for i = 1:numel(flds)
                f = flds{i}; val = p.Results.(f);
                if ischar(val) || isstring(val), obj.(f) = string(val);
                else, obj.(f) = val; end
            end

            unm = p.Unmatched;
            if isfield(unm, 'max_iter'), obj.MaxIter = double(unm.max_iter); end
            if isfield(unm, 'tol'), obj.Tol = double(unm.tol); end
            if isfield(unm, 'save_dir'), obj.SaveDir = string(unm.save_dir); end

            obj.updateScaling();
            obj.setManeuver(obj.Maneuver);
        end

        function updateScaling(obj)
            %% UPDATESCALING  Build state and input characteristic scaling vectors.
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

        function setManeuver(obj, name, varargin)
            %% SETMANEUVER  Configure maneuver parameters dynamically.
            obj.Maneuver = string(name);
            p = inputParser;

            if obj.Vehicle == 0
                dR = 5.0;  dA = 7.0;  dDR = 2.5; dCB = [14.0, 35.0]; dCI = 22.0;
                dBA = 20.0; dBB = [10.0, 22.0]; dBI = 16.0;
                dHA = 15.0; dHB = [12.0, 30.0]; dHI = 18.0;
            else
                dR = 15.0; dA = 25.0; dDR = 5.0; dCB = [16.0, 45.0]; dCI = 26.0;
                dBA = 50.0; dBB = [15.0, 35.0]; dBI = 24.0;
                dHA = 50.0; dHB = [15.0, 40.0]; dHI = 25.0;
            end

            switch obj.Maneuver
                case "Circle"
                    addParameter(p, 'circle_radius', dR, @isnumeric);
                    addParameter(p, 'circle_alt', dA, @isnumeric);
                    addParameter(p, 'circle_center', [0; 0], @isnumeric);
                    addParameter(p, 'descent_glideslope', 1.0, @isnumeric);
                    addParameter(p, 'max_descent_rate', dDR, @isnumeric);
                    addParameter(p, 'wp_tol', 0.25, @isnumeric);
                    addParameter(p, 'f_orbit_start', 0.25, @isnumeric);
                    addParameter(p, 'f_orbit_end', 0.75, @isnumeric);
                    parse(p, varargin{:});

                    obj.T_bounds = dCB; obj.T_initial = dCI;
                    fs = p.Results.f_orbit_start; fe = p.Results.f_orbit_end;
                    Ns = round(fs * obj.N); Ne = round(fe * obj.N);
                    obj.ManeuverParams = struct( ...
                        'circle_radius', double(p.Results.circle_radius), ...
                        'circle_alt', double(p.Results.circle_alt), ...
                        'circle_center', double(p.Results.circle_center(:)), ...
                        'descent_glideslope', double(p.Results.descent_glideslope), ...
                        'max_descent_rate', double(p.Results.max_descent_rate), ...
                        'wp_tol', double(p.Results.wp_tol), ...
                        'f_orbit_start', fs, 'f_orbit_end', fe, ...
                        'N_orbit_start', Ns, 'N_orbit_end', Ne, ...
                        'N_ascent', Ns, 'N_ingress', Ns, 'N_orbit', Ne, 'N_egress', Ne, 'N_descent', obj.N);

                case "Backflip"
                    addParameter(p, 'apex_alt', dBA, @isnumeric);
                    addParameter(p, 'flip_start_frac', 0.30, @isnumeric);
                    addParameter(p, 'flip_end_frac', 0.70, @isnumeric);
                    addParameter(p, 'theta_tol', deg2rad(15), @isnumeric);
                    addParameter(p, 'Glideslope', tand(25), @isnumeric);
                    parse(p, varargin{:});

                    obj.T_bounds = dBB; obj.T_initial = dBI;
                    obj.ManeuverParams = struct( ...
                        'apex_alt', double(p.Results.apex_alt), ...
                        'theta_tol', double(p.Results.theta_tol), ...
                        'Glideslope', double(p.Results.Glideslope), ...
                        'N_ascent', round(p.Results.flip_start_frac * obj.N), ...
                        'N_flip', round(0.50 * obj.N), ...
                        'N_approach', round(p.Results.flip_end_frac * obj.N), ...
                        'q_inverted', [0; 0; 1; 0]);

                case "Hop"
                    addParameter(p, 'apex_alt', dHA, @isnumeric);
                    addParameter(p, 'Glideslope', tand(25), @isnumeric);
                    parse(p, varargin{:});

                    obj.T_bounds = dHB; obj.T_initial = dHI;
                    obj.ManeuverParams = struct( ...
                        'apex_alt', double(p.Results.apex_alt), ...
                        'Glideslope', double(p.Results.Glideslope), ...
                        'N_ascent', round(0.20 * obj.N), ...
                        'N_approach', round(0.80 * obj.N));

                otherwise
                    obj.ManeuverParams = struct('Name', obj.Maneuver);
            end
            obj.InitialGuess = struct();
        end

        function setBoundaries(obj, r0, r_f, varargin)
            %% SETBOUNDARIES  Configure launch pad and landing zone positions.
            obj.r0 = r0(:); obj.r_f = r_f(:);
            p = inputParser;
            addParameter(p, 'v0', [0; 0; 0], @isnumeric);
            addParameter(p, 'q0', [1; 0; 0; 0], @isnumeric);
            addParameter(p, 'v_f_tol', 0.1, @isnumeric);
            parse(p, varargin{:});
            obj.v0 = p.Results.v0(:); obj.q0 = p.Results.q0(:);
            obj.v_f_tol = double(p.Results.v_f_tol);
            obj.InitialGuess = struct();
        end

        function addWaypoint(obj, node_idx, pos, tol)
            if nargin < 4, tol = 0.5; end
            obj.CustomWaypoints = [obj.CustomWaypoints; struct('node_idx', node_idx, 'pos', pos(:), 'tol', tol)];
        end

        function clearWaypoints(obj), obj.CustomWaypoints = []; end

        %% ═══════════════════ INITIAL GUESS ═══════════════════

        function guess = generateInitialGuess(obj)
            %% GENERATEINITIALGUESS  Analytical quintic-spline self-referencing guess.
            Nn = obj.N + 1; Nc = obj.N; Ts = obj.T_initial;

            if obj.TimeParamMode == "MultiPhase" && obj.Maneuver == "Circle"
                mp = obj.ManeuverParams;
                N1 = mp.N_orbit_start; N2 = mp.N_orbit_end - N1; N3 = obj.N - mp.N_orbit_end;
                dt_vec = [repmat(0.25*Ts/N1, 1, N1), repmat(0.50*Ts/N2, 1, N2), repmat(0.25*Ts/N3, 1, N3)];
                tv = [0, cumsum(dt_vec)];
            else
                dt_vec = repmat(Ts / obj.N, 1, obj.N);
                tv = linspace(0, Ts, Nn);
            end

            m_dry = obj.constants.m_dry; g = obj.constants.g; MT = obj.constants.MaxThrust;
            thr_lo = (0.25 + obj.thrust_margin) * MT; thr_hi = (1.00 - obj.thrust_margin) * MT;

            if obj.Vehicle == 0
                m_v = m_dry * ones(1, Nn); ml = zeros(1, Nn); mi = zeros(1, Nn);
            else
                ml = linspace(obj.constants.OxMass, 0.20 * obj.constants.OxMass, Nn);
                mi = linspace(obj.constants.FuMass, 0.20 * obj.constants.FuMass, Nn);
                m_v = m_dry + ml + mi;
            end

            r_s = zeros(3, Nn); v_s = zeros(3, Nn); a_s = zeros(3, Nn);

            switch obj.Maneuver
                case "Circle"
                    mp = obj.ManeuverParams;
                    cx = mp.circle_center(1); cy = mp.circle_center(2);
                    R  = mp.circle_radius;     h  = mp.circle_alt;
                    ks = mp.N_orbit_start + 1; ke = mp.N_orbit_end + 1;
                    t1 = tv(ks); t2 = tv(ke); Torb = max(1e-2, t2 - t1); wc = 2 * pi / Torb;
                    re = [cx + R; cy; h]; ve = [0; R * wc; 0]; ae = [-R * wc^2; 0; 0];

                    for k = 1:ks
                        tau = tv(k) / max(1e-2, t1);
                        [r_s(:,k), v_s(:,k), a_s(:,k)] = obj.evalQuinticSpline( ...
                            obj.r0, [0; 0; 0.5], [0; 0; 0], re, ve, ae, tau, t1);
                    end
                    for k = (ks + 1):ke
                        th = 2 * pi * (tv(k) - t1) / Torb;
                        r_s(:,k) = [cx + R * cos(th); cy + R * sin(th); h];
                        v_s(:,k) = [-R * wc * sin(th); R * wc * cos(th); 0];
                        a_s(:,k) = [-R * wc^2 * cos(th); -R * wc^2 * sin(th); 0];
                    end
                    td = max(1e-2, tv(end) - t2);
                    for k = (ke + 1):Nn
                        tau = (tv(k) - t2) / td;
                        [r_s(:,k), v_s(:,k), a_s(:,k)] = obj.evalQuinticSpline( ...
                            re, ve, ae, obj.r_f, [0; 0; -0.2], [0; 0; 0], tau, td);
                    end

                case {"Hop", "Backflip"}
                    mp = obj.ManeuverParams;
                    apex = [0.5 * (obj.r0(1) + obj.r_f(1)); 0.5 * (obj.r0(2) + obj.r_f(2)); mp.apex_alt];
                    ta = 0.5 * Ts; vh = (obj.r_f - obj.r0) / max(1e-2, Ts); va = [vh(1); vh(2); 0];
                    for k = 1:Nn
                        t = tv(k);
                        if t <= ta
                            tau = t / ta;
                            [r_s(:,k), v_s(:,k), a_s(:,k)] = obj.evalQuinticSpline( ...
                                obj.r0, [0; 0; 0.5], [0; 0; 0], apex, va, [0; 0; 0], tau, ta);
                        else
                            tau = (t - ta) / (Ts - ta);
                            [r_s(:,k), v_s(:,k), a_s(:,k)] = obj.evalQuinticSpline( ...
                                apex, va, [0; 0; 0], obj.r_f, [0; 0; -0.2], [0; 0; 0], tau, Ts - ta);
                        end
                    end

                otherwise
                    for k = 1:Nn
                        frac = (k - 1) / obj.N;
                        r_s(:,k) = obj.r0 + (obj.r_f - obj.r0) * frac;
                        v_s(:,k) = (obj.r_f - obj.r0) / Ts;
                        a_s(:,k) = [0; 0; 0];
                    end
            end
            r_s(:,1) = obj.r0; v_s(:,1) = obj.v0; r_s(:,end) = obj.r_f; v_s(:,end) = [0; 0; 0];

            q_s = zeros(4, Nn); w_s = zeros(3, Nn); U_p = zeros(4, Nc);

            if obj.Maneuver == "Backflip"
                ks = obj.ManeuverParams.N_ascent; ke = obj.ManeuverParams.N_approach;
                Lf = max(2, ke - ks); Tf = Lf * dt_vec(1);
                for k = 1:Nn
                    if k <= ks
                        q_s(:,k) = obj.q0;
                    elseif k <= ke
                        tau = (k - ks) / Lf; th = -2*pi*(tau - sin(2*pi*tau)/(2*pi));
                        dth = -(2*pi/Tf)*(1 - cos(2*pi*tau));
                        q_s(:,k) = [cos(th/2); 0; sin(th/2); 0]; w_s(:,k) = [0; dth; 0];
                    else
                        q_s(:,k) = -obj.q0;
                    end
                end
            else
                for k = 1:Nn, q_s(:,k) = obj.vectorToQuat(m_v(k) * (a_s(:,k) + [0; 0; g])); end
                q_s(:,1) = obj.q0; q_s(:,end) = obj.q0;
                for k = 1:Nc
                    dq = (q_s(:,k+1) - q_s(:,k)) / dt_vec(k);
                    qw = q_s(1,k); qv = q_s(2:4,k);
                    w_s(:,k) = 2 * (qw * dq(2:4) - dq(1) * qv - cross(qv, dq(2:4)));
                end
            end

            for k = 1:Nc
                F_req = m_v(k) * (a_s(:,k) + [0; 0; g]);
                T_cmd = max(min(norm(F_req), thr_hi), thr_lo);
                U_p(:,k) = [0; 0; T_cmd; 0];
            end
            q_s(:,1) = obj.q0; w_s(:,1) = [0; 0; 0]; w_s(:,end) = [0; 0; 0];

            Xp = [q_s; r_s; v_s; w_s; ml; mi];
            guess = struct('Time', tv, 'X', Xp, 'U', U_p, 'T_total', Ts, ...
                           'Xhat', Xp ./ obj.Sx, 'Uhat', U_p ./ obj.Su);
            obj.InitialGuess = guess;
        end

        %% ═══════════════════ OPTIMIZATION SOLVER ═══════════════════

        function sol = solve(obj)
            %% SOLVE  Assemble CasADi Opti problem and execute IPOPT.
            import casadi.*
            N = obj.N; opti = Opti(); dyn_fnc = obj.getCasADiDynamics();
            obj.OptiVars = struct();

            % 1. Timing Parameterization
            if obj.Maneuver == "Circle" && obj.TimeParamMode == "MultiPhase"
                mp = obj.ManeuverParams;
                N1 = mp.N_orbit_start; N2 = mp.N_orbit_end - N1; N3 = N - mp.N_orbit_end;
                tmin = obj.T_bounds(1); tmax = obj.T_bounds(2);
                T_asc  = opti.variable(); opti.subject_to(0.10*tmin <= T_asc);  opti.subject_to(T_asc <= 0.45*tmax);
                T_orb  = opti.variable(); opti.subject_to(0.20*tmin <= T_orb);  opti.subject_to(T_orb <= 0.70*tmax);
                T_desc = opti.variable(); opti.subject_to(0.10*tmin <= T_desc); opti.subject_to(T_desc <= 0.45*tmax);
                opti.set_initial(T_asc,  0.25 * obj.T_initial);
                opti.set_initial(T_orb,  0.50 * obj.T_initial);
                opti.set_initial(T_desc, 0.25 * obj.T_initial);
                T_total = T_asc + T_orb + T_desc;
                opti.subject_to(tmin <= T_total); opti.subject_to(T_total <= tmax);
                dt_row = [repmat(T_asc/N1, 1, N1), repmat(T_orb/N2, 1, N2), repmat(T_desc/N3, 1, N3)];
                obj.OptiVars.T_asc = T_asc; obj.OptiVars.T_orb = T_orb; obj.OptiVars.T_desc = T_desc;
            else
                T_total = opti.variable();
                opti.subject_to(obj.T_bounds(1) <= T_total); opti.subject_to(T_total <= obj.T_bounds(2));
                opti.set_initial(T_total, obj.T_initial);
                dt_row = repmat(T_total / N, 1, N);
                obj.OptiVars.T_asc = []; obj.OptiVars.T_orb = []; obj.OptiVars.T_desc = [];
            end

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
            opti.subject_to(sum(X(1:4, :).^2, 1) == 1.00);
            opti.subject_to(-35 <= X(5, :)); opti.subject_to(X(5, :) <= 35);
            opti.subject_to(-35 <= X(6, :)); opti.subject_to(X(6, :) <= 35);
            alt_ceil = 75 * (obj.Vehicle == 0) + 150 * (obj.Vehicle == 1);
            opti.subject_to(-1 <= X(7, :)); opti.subject_to(X(7, :) <= alt_ceil);

            m_lox0 = obj.constants.OxMass; m_ipa0 = obj.constants.FuMass;
            opti.subject_to(X(:, 1) == [obj.q0; obj.r0; obj.v0; obj.w0; m_lox0; m_ipa0]);
            if obj.Maneuver == "Backflip", opti.subject_to(X(1:4, end) == -obj.q0);
            else, opti.subject_to(X(1:4, end) == obj.q0); end
            opti.subject_to(X(5:7, end) == obj.r_f);
            opti.subject_to(sum(X(8:10, end).^2) <= obj.v_f_tol^2);

            if obj.Vehicle == 1
                opti.subject_to(X(14, end) >= 0.10 * m_lox0);
                opti.subject_to(X(15, end) >= 0.10 * m_ipa0);
            else
                opti.subject_to(X(14:15, :) == 0);
            end

            % 4. Control Bounds & Slew Limits
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

            % 5. Maneuver Constraints & Multi-Objective Cost
            obj.applyManeuverConstraints(opti, X);
            obj.applyCostFunction(opti, X, Xhat, U, Uhat, T_total, dt_row);

            % 6. Initial Guess Seeding
            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'Xhat'), obj.generateInitialGuess(); end
            if ~(obj.Maneuver == "Circle" && obj.TimeParamMode == "MultiPhase")
                opti.set_initial(T_total, obj.InitialGuess.T_total);
            end
            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);
            obj.OptiVars.Xhat = Xhat; obj.OptiVars.Uhat = Uhat;
            obj.OptiVars.T_total = T_total; obj.OptiVars.dt_row = dt_row;

            % 7. Solver Setup & Execution
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, 'tol', obj.Tol, ...
                'constr_viol_tol', obj.ConstrViolTol, 'acceptable_tol', 1e-3, ...
                'acceptable_constr_viol_tol', 1e-3, 'acceptable_iter', 5, 'print_level', obj.PrintLevel);
            opti.solver('ipopt', p_opts, s_opts);

            if obj.PrintLevel > 0
                fprintf('Starting %s solve for %s (TimeMode: %s, CircleMode: %s)...\n', ...
                    obj.Maneuver, obj.getVehicleName(), obj.TimeParamMode, obj.CircleMode);
            end

            t0 = tic;
            try
                sc = opti.solve(); status = 'Success';
                X_r = obj.Sx .* sc.value(Xhat); U_r = obj.Su .* sc.value(Uhat);
                T_r = sc.value(T_total); dt_r = sc.value(dt_row); so = sc;
                stats = sc.stats(); stats.t_wall_total = toc(t0);

                if obj.TimeParamMode == "MultiPhase" && obj.Maneuver == "Circle"
                    ta = sc.value(obj.OptiVars.T_asc); to = sc.value(obj.OptiVars.T_orb);
                    obj.SwitchSpots = [ta, ta + to, T_r];
                else
                    mp = obj.ManeuverParams;
                    if isfield(mp, 'N_orbit_start') && isfield(mp, 'N_orbit_end')
                        obj.SwitchSpots = [(mp.N_orbit_start/N)*T_r, (mp.N_orbit_end/N)*T_r, T_r];
                    else
                        obj.SwitchSpots = [0.5*T_r, T_r];
                    end
                end
            catch ME
                tw = toc(t0); status = 'Infeasible';
                X_r = obj.Sx .* opti.debug.value(Xhat); U_r = obj.Su .* opti.debug.value(Uhat);
                T_r = opti.debug.value(T_total); dt_r = opti.debug.value(dt_row); so = ME;
                stats = struct('t_wall_total', tw, 'iter_count', obj.MaxIter, 'return_status', 'Infeasible');
                obj.SwitchSpots = [0.25*T_r, 0.75*T_r, T_r];
            end

            t_r = [0, cumsum(dt_r)];
            if obj.Vehicle == 0, X_r(14:15, :) = 0; end

            if obj.Vehicle == 0, mv = obj.constants.m_dry * ones(1, N);
            else, mv = obj.constants.m_dry + X_r(14, 1:N) + X_r(15, 1:N); end
            uh = mv * obj.constants.g / MT; udT = U_r(3, :) / MT - uh;
            gb = max(1e-3, (1 - obj.gimbal_margin) * mg);
            esq = udT.^2 + (U_r(1, :) / gb).^2 + (U_r(2, :) / gb).^2 + ...
                  (U_r(4, :) / max(1e-3, (1 - tm) * obj.max_roll_rate)).^2;
            CE = sum(dt_r .* esq);
            dr_r = diff(X_r(5:7, :), 1, 2); LP = sum(sqrt(sum(dr_r.^2, 1)));
            dUr = diff(U_r, 1, 2); dtsr = 0.5 * (dt_r(1:end-1) + dt_r(2:end));
            ES = sum(sum((dUr ./ dtsr).^2, 1));

            obj.Solution = struct( ...
                'Time', t_r, 't', t_r, 'X', X_r, 'x', X_r, 'U', U_r, 'u', U_r, ...
                'T_total', T_r, 'L_path', LP, 'ControlEnergy', CE, 'E_slew', ES, ...
                'Var_u', var(U_r, 0, 2)', 'Var_omega', var(X_r(11:13, :), 0, 2)', ...
                'SwitchSpots', obj.SwitchSpots, 'TimeParamMode', obj.TimeParamMode, ...
                'CircleMode', obj.CircleMode, 'Status', status, 'stats', stats, 'sol', so);
            sol = obj.Solution;

            if obj.PlotResults, obj.plot(); end
            if obj.AutoExport,  obj.exportCSV(); end
        end

        %% ═══════════════════ CONSTRAINTS & COST ═══════════════════

        function applyManeuverConstraints(obj, opti, X)
            %% APPLYMANEUVERCONSTRAINTS  Maneuver-tailored constraint injection.
            p = obj.ManeuverParams;

            switch obj.Maneuver
                case "Circle"
                    cx = p.circle_center(1); cy = p.circle_center(2);
                    R  = p.circle_radius;     h  = p.circle_alt;
                    Nos = p.N_orbit_start;   Noe = p.N_orbit_end; Norb = Noe - Nos;

                    % Stage 1: Ascent glideslope cone
                    kc = max(2, round(0.08 * obj.N));
                    opti.subject_to(sqrt(X(5, 1:kc).^2 + X(6, 1:kc).^2 + 1e-4) <= X(7, 1:kc) * tand(5) + 0.15);
                    opti.subject_to(X(5, 1:Nos) >= min(obj.r0(1), cx + R) - 0.10);
                    opti.subject_to(X(10, 1:Nos) >= -0.10);

                    % Stage 2: Circle Orbit
                    if obj.CircleMode == "AngleSweep"
                        rsq = (X(5, Nos:Noe) - cx).^2 + (X(6, Nos:Noe) - cy).^2; tr = obj.CircleTightness;
                        opti.subject_to((R - tr)^2 <= rsq); opti.subject_to(rsq <= (R + tr)^2);
                        opti.subject_to(abs(X(7, Nos:Noe) - h) <= 0.20);

                        % Bilinear cross-product angle progression
                        rx = X(5, Nos:Noe) - cx; ry = X(6, Nos:Noe) - cy;
                        cprog = rx(1:end-1) .* ry(2:end) - ry(1:end-1) .* rx(2:end);
                        opti.subject_to(cprog >= (R^2) * sin(2 * pi / Norb * 0.70));

                        % Quadrants & terminal closure
                        kq1 = Nos + round(0.25*Norb); kq2 = Nos + round(0.50*Norb); kq3 = Nos + round(0.75*Norb);
                        opti.subject_to(X(6, kq1) - cy >= 0);
                        opti.subject_to(X(5, kq2) - cx <= 0);
                        opti.subject_to(X(6, kq3) - cy <= 0);
                        opti.subject_to((X(5, Noe) - (cx + R))^2 + (X(6, Noe) - cy)^2 <= tr^2);
                        opti.subject_to((X(7, Noe) - h)^2 <= 0.10^2);
                    else
                        rsq = (X(5, Nos:Noe) - cx).^2 + (X(6, Nos:Noe) - cy).^2;
                        opti.subject_to((R - 0.50)^2 <= rsq); opti.subject_to(rsq <= (R + 0.50)^2);
                        rx = X(5, Nos:Noe) - cx; ry = X(6, Nos:Noe) - cy;
                        opti.subject_to(rx .* X(9, Nos:Noe) - ry .* X(8, Nos:Noe) >= (R^2 * 0.10));
                        kq1 = Nos + round(0.25*Norb); kq2 = Nos + round(0.50*Norb); kq3 = Nos + round(0.75*Norb);
                        opti.subject_to(X(6, kq1) - cy >= 0); opti.subject_to(X(5, kq2) - cx <= 0); opti.subject_to(X(6, kq3) - cy <= 0);
                        opti.subject_to((X(5, Noe) - (cx + R))^2 + (X(6, Noe) - cy)^2 <= 0.25^2);
                    end

                    % Stage 3: Descent glideslope
                    rf = obj.r_f;
                    rxy = sqrt((X(5, Noe:end) - rf(1)).^2 + (X(6, Noe:end) - rf(2)).^2 + 1e-4);
                    opti.subject_to(rxy <= (X(7, Noe:end) - rf(3)) * p.descent_glideslope + 0.50);
                    opti.subject_to(X(5, Noe:end) >= min(rf(1), cx + R) - 0.10);
                    opti.subject_to(-p.max_descent_rate <= X(10, Noe:end)); opti.subject_to(X(10, Noe:end) <= 0.10);
                    R33d = X(1, Noe:end).^2 - X(2, Noe:end).^2 - X(3, Noe:end).^2 + X(4, Noe:end).^2;
                    opti.subject_to(R33d >= cosd(30));

                case "Backflip"
                    Na = p.N_ascent; Nf = p.N_flip; Nap = p.N_approach; GS = p.Glideslope;
                    opti.subject_to(X(10, 1:Na) >= 0);
                    opti.subject_to(sqrt(X(5, 1:Na).^2 + X(6, 1:Na).^2 + 1e-4) <= X(7, 1:Na) * GS + 0.05);
                    opti.subject_to(X(7, :) <= p.apex_alt + 2.5);
                    opti.subject_to(X(7, Nf) >= p.apex_alt - 2.0); opti.subject_to(X(7, Na) >= 0.50 * p.apex_alt);
                    att_tol = cos(p.theta_tol / 2);
                    opti.subject_to(p.q_inverted' * X(1:4, Nf) >= att_tol);
                    opti.subject_to([-1; 0; 0; 0]' * X(1:4, Nap) >= att_tol);
                    opti.subject_to(X(12, Na:Nap) <= 0.05);
                    R33d = X(1, Nap:end).^2 - X(2, Nap:end).^2 - X(3, Nap:end).^2 + X(4, Nap:end).^2;
                    opti.subject_to(R33d >= cosd(45));

                case "Hop"
                    Na = p.N_ascent; Nap = p.N_approach; GS = p.Glideslope;
                    opti.subject_to(X(10, 1:Na) >= 0); opti.subject_to(X(10, Nap:end) <= 0.10);
                    R33 = X(1, :).^2 - X(2, :).^2 - X(3, :).^2 + X(4, :).^2; opti.subject_to(R33 >= cosd(35));
                    pd = X(5:7, Nap:end);
                    opti.subject_to(sqrt((pd(1,:) - obj.r_f(1)).^2 + (pd(2,:) - obj.r_f(2)).^2 + 1e-4) <= ...
                        (pd(3,:) - obj.r_f(3)) * GS + 0.50);
            end

            for i = 1:length(obj.CustomWaypoints)
                wp = obj.CustomWaypoints(i);
                opti.subject_to(sum((X(5:7, wp.node_idx) - wp.pos).^2) <= wp.tol^2);
            end
        end

        function applyCostFunction(obj, opti, X, Xhat, U, Uhat, T_total, dt_row)
            %% APPLYCOSTFUNCTION  Minimal-time, minimal-path-length Pareto cost.
            N = obj.N; g = obj.constants.g; MT = obj.constants.MaxThrust;

            J_time = T_total / obj.T_initial;
            dr = X(5:7, 2:end) - X(5:7, 1:end-1);
            L_path = sum(sqrt(sum(dr.^2, 1) + 1e-6));
            if obj.Maneuver == "Circle"
                L_ref = norm(obj.r_f - obj.r0) + 2 * pi * obj.ManeuverParams.circle_radius + 2 * obj.ManeuverParams.circle_alt;
            elseif obj.Maneuver == "Backflip" || obj.Maneuver == "Hop"
                L_ref = 2 * obj.ManeuverParams.apex_alt + norm(obj.r_f - obj.r0);
            else
                L_ref = max(1.0, norm(obj.r_f - obj.r0) + 20.0);
            end
            J_length = L_path / L_ref;

            if obj.Vehicle == 0, mv = obj.constants.m_dry * ones(1, N);
            else, mv = obj.constants.m_dry + X(14, 1:N) + X(15, 1:N); end
            uh = mv * g / MT; gb = max(1e-3, (1 - obj.gimbal_margin) * obj.max_gimbal_angle);
            esq = (U(3, :) / MT - uh).^2 + (U(1, :) / gb).^2 + (U(2, :) / gb).^2 + ...
                  (U(4, :) / max(1e-3, (1 - obj.thrust_margin) * obj.max_roll_rate)).^2;
            J_effort = sum(dt_row .* esq) / (0.20 * obj.T_initial);

            dUh = Uhat(:, 2:end) - Uhat(:, 1:end-1);
            J_slew = sum(sum(dUh.^2, 1)) / (N - 1);
            J_rate = sum(sum(Xhat(11:13, 1:end-1).^2, 1)) / N;
            J_qz   = sum(Xhat(4, :).^2) / (N + 1);

            J = obj.w_time * J_time + obj.w_length * J_length + ...
                    obj.w_effort * J_effort + obj.w_slew * J_slew + obj.w_rate * J_rate + obj.w_qz * J_qz;
            opti.minimize(J);
        end

        %% ═══════════════════ 6-DOF DYNAMICS ═══════════════════

        function dyn_fnc = getCasADiDynamics(obj)
            %% GETCASADIDYNAMICS  Symbolic CasADi 6-DoF dynamics function.
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

            if obj.Vehicle == 0
                mdl = 0; mdi = 0; OFH = 0; FFH = 0;
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
            dyn_fnc = Function('dynamics_fnc', {x, u, params}, {[qdot; v; FI / m; wdot; mdl; mdi]});
        end

        %% ═══════════════════ UTILITIES & EXPORT ═══════════════════

        function [rk, vk, ak] = evalQuinticSpline(~, r0, v0, a0, r1, v1, a1, tau, Ts)
            %% EVALQUINTICSPLINE  C2 quintic polynomial trajectory evaluation.
            tau = max(0, min(1, tau)); Ts = max(Ts, 1e-3);
            vs = v0 * Ts; v1s = v1 * Ts; as = a0 * Ts^2; a1s = a1 * Ts^2;
            c0 = r0; c1 = vs; c2 = 0.5 * as;
            dx = r1 - (c0 + c1 + c2); dv = v1s - (c1 + 2 * c2); da = a1s - 2 * c2;
            c3 = 10 * dx - 4 * dv + 0.5 * da; c4 = -15 * dx + 7 * dv - da; c5 = 6 * dx - 3 * dv + 0.5 * da;
            rk = c0 + c1*tau + c2*tau^2 + c3*tau^3 + c4*tau^4 + c5*tau^5;
            vk = (c1 + 2*c2*tau + 3*c3*tau^2 + 4*c4*tau^3 + 5*c5*tau^4) / Ts;
            ak = (2*c2 + 6*c3*tau + 12*c4*tau^2 + 20*c5*tau^3) / Ts^2;
        end

        function q = vectorToQuat(~, vec)
            %% VECTORTOQUAT  Unit quaternion aligning body z-axis with vector.
            nv = norm(vec);
            if nv < 1e-4, q = [1; 0; 0; 0]; return; end
            zt = vec(:) / nv; cv = cross([0; 0; 1], zt); dv = dot([0; 0; 1], zt);
            if dv < -0.9999, q = [0; 1; 0; 0];
            else, q = [1 + dv; cv]; q = q / norm(q); end
        end

        function name = getVehicleName(obj)
            if obj.Vehicle == 0, name = "ASTRA"; else, name = "TOAD"; end
        end

        function fn = getFormattedFilename(obj)
            if obj.Filename ~= ""
                fn = obj.Filename;
            else
                fn = sprintf('%s_%s_v%03d', obj.getVehicleName(), obj.Maneuver, obj.Version);
            end
        end

        function figHandles = plot(obj)
            %% PLOT  Generate clean, interactive diagnostic figures for trajectory analysis.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                error('TrajectoryOptimizer:NoSolution', 'No solution available to plot. Call solve() first.');
            end

            tv = obj.Solution.Time; Xs = obj.Solution.X; Us = obj.Solution.U;

            f1 = figure('Name', sprintf('%s 6-DoF Trajectory (%s v%03d)', ...
                        obj.Maneuver, obj.getVehicleName(), obj.Version), ...
                        'Position', [100, 100, 1200, 800], 'Visible', 'on');

            % 1. 3D Flight Path Overview with Time-Scheduled Blue Gradient
            ax3d = subplot(2, 2, [1, 3]);
            % Color-coded line: Light Sky Blue (Takeoff) -> Deep Navy (Touchdown)
            h_traj = surface([Xs(5,:); Xs(5,:)], [Xs(6,:); Xs(6,:)], [Xs(7,:); Xs(7,:)], ...
                             [tv; tv], 'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.8);
            hold on;
            plot3(Xs(5,:), Xs(6,:), zeros(size(Xs(7,:))), 'Color', [0.75 0.75 0.75], 'LineStyle', ':');
            p_launch = plot3(obj.r0(1), obj.r0(2), obj.r0(3), 'go', 'MarkerSize', 10, 'LineWidth', 2.5);
            p_land   = plot3(obj.r_f(1), obj.r_f(2), obj.r_f(3), 'rs', 'MarkerSize', 10, 'LineWidth', 2.5);

            % Blue gradient colormap
            n_c = 256;
            blue_cmap = [linspace(0.72, 0.02, n_c)', linspace(0.90, 0.12, n_c)', linspace(1.00, 0.58, n_c)'];
            colormap(ax3d, blue_cmap);
            cb = colorbar('Location', 'southoutside');
            cb.Label.String = 'Mission Time [s] (Light: Takeoff \rightarrow Dark: Landing)';
            cb.FontSize = 8.5;

            % Minimum spatial span to avoid planar collapse
            x_data = [Xs(5,:), obj.r0(1), obj.r_f(1)];
            y_data = [Xs(6,:), obj.r0(2), obj.r_f(2)];
            z_data = [Xs(7,:), obj.r0(3), obj.r_f(3)];
            if obj.Maneuver == "Circle"
                thc = linspace(0, 2*pi, 100); mp = obj.ManeuverParams;
                p_orb = plot3(mp.circle_center(1) + mp.circle_radius * cos(thc), ...
                              mp.circle_center(2) + mp.circle_radius * sin(thc), ...
                              mp.circle_alt * ones(size(thc)), 'k--', 'LineWidth', 1.5);
                x_data = [x_data, mp.circle_center(1) - mp.circle_radius, mp.circle_center(1) + mp.circle_radius];
                y_data = [y_data, mp.circle_center(2) - mp.circle_radius, mp.circle_center(2) + mp.circle_radius];
                z_data = [z_data, mp.circle_alt];
                legend([h_traj, p_launch, p_land, p_orb], ...
                       {'Flight Path (Time-Scheduled)', 'Launch Pad', 'Landing Zone', 'Nominal Orbit'}, ...
                       'Location', 'best');
            else
                legend([h_traj, p_launch, p_land], ...
                       {'Flight Path (Time-Scheduled)', 'Launch Pad', 'Landing Zone'}, ...
                       'Location', 'best');
            end

            x_span = max(x_data) - min(x_data);
            y_span = max(y_data) - min(y_data);
            z_span = max(z_data) - min(z_data);
            max_span = max([x_span, y_span, z_span, 5.0]);
            min_allowed_span = max(6.0, 0.35 * max_span);

            x_mid = 0.5 * (min(x_data) + max(x_data));
            y_mid = 0.5 * (min(y_data) + max(y_data));
            half_x = max(0.5 * x_span + 1.0, min_allowed_span / 2);
            half_y = max(0.5 * y_span + 1.0, min_allowed_span / 2);
            z_top  = max(max(z_data) + 1.5, min_allowed_span);

            grid on;
            view(-35, 25);
            xlim([x_mid - half_x, x_mid + half_x]);
            ylim([y_mid - half_y, y_mid + half_y]);
            zlim([0, z_top]);
            pbaspect([1, 1, 0.8]);
            xlabel('X (North) [m]'); ylabel('Y (East) [m]'); zlabel('Z (Altitude) [m]');
            title(sprintf('%s Flight Profile (%s, Duration: %.2f s)', ...
                obj.Maneuver, obj.getVehicleName(), obj.Solution.T_total));

            % 2. TVC Gimbal Deflections
            subplot(2, 2, 2);
            plot(tv(1:end-1), rad2deg(Us(1,:)), 'r-', 'LineWidth', 1.6); hold on;
            plot(tv(1:end-1), rad2deg(Us(2,:)), 'b-', 'LineWidth', 1.6);
            yline(rad2deg(obj.max_gimbal_angle), 'k:', 'Max Gimbal');
            yline(-rad2deg(obj.max_gimbal_angle), 'k:');
            grid on; xlabel('Time [s]'); ylabel('Gimbal Angle [deg]');
            legend({'\theta (Pitch)', '\phi (Yaw)'}, 'Location', 'best');
            title('TVC Gimbal Deflections');

            % 3. Commanded Thrust vs Time
            subplot(2, 2, 4);
            plot(tv(1:end-1), Us(3,:), 'k-', 'LineWidth', 1.8); hold on;
            hover_thrust = obj.constants.m_dry * obj.constants.g;
            yline(hover_thrust, 'r--', sprintf('Dry Hover (%.1f N)', hover_thrust));
            yline(obj.constants.MaxThrust, 'b:', sprintf('Max Thrust (%.1f N)', obj.constants.MaxThrust));
            grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
            title(sprintf('Commanded Thrust Profile (Path Length: %.2f m)', obj.Solution.L_path));
            legend({'Thrust Cmd', 'Hover Ref', 'Max Limit'}, 'Location', 'best');

            figHandles = f1;
        end

        function figHandle = plotInitialGuess(obj)
            %% PLOTINITIALGUESS  Visualize analytical initial guess.
            if isempty(obj.InitialGuess) || ~isfield(obj.InitialGuess, 'X'), obj.generateInitialGuess(); end
            tv = obj.InitialGuess.Time; Xg = obj.InitialGuess.X; Ug = obj.InitialGuess.U;

            figHandle = figure('Name', sprintf('%s Initial Guess (%s)', obj.Maneuver, obj.getVehicleName()), ...
                               'Position', [150, 150, 950, 650], 'Visible', 'on');
            subplot(2, 1, 1);
            plot3(Xg(5,:), Xg(6,:), Xg(7,:), 'm--', 'LineWidth', 2); hold on;
            plot3(obj.r0(1), obj.r0(2), obj.r0(3), 'go', 'MarkerSize', 8, 'LineWidth', 2);
            plot3(obj.r_f(1), obj.r_f(2), obj.r_f(3), 'rs', 'MarkerSize', 8, 'LineWidth', 2);

            % Minimum span to avoid planar collapse
            x_data = [Xg(5,:), obj.r0(1), obj.r_f(1)];
            y_data = [Xg(6,:), obj.r0(2), obj.r_f(2)];
            z_data = [Xg(7,:), obj.r0(3), obj.r_f(3)];
            x_span = max(x_data) - min(x_data);
            y_span = max(y_data) - min(y_data);
            z_span = max(z_data) - min(z_data);
            max_span = max([x_span, y_span, z_span, 5.0]);
            min_allowed = max(6.0, 0.40 * max_span);
            x_mid = 0.5 * (min(x_data) + max(x_data));
            y_mid = 0.5 * (min(y_data) + max(y_data));
            half_x = max(0.5 * x_span + 1.0, min_allowed / 2);
            half_y = max(0.5 * y_span + 1.0, min_allowed / 2);
            z_top  = max(max(z_data) + 1.5, min_allowed);

            grid on;
            view(-35, 25);
            xlim([x_mid - half_x, x_mid + half_x]);
            ylim([y_mid - half_y, y_mid + half_y]);
            zlim([0, z_top]);
            pbaspect([1, 1, 0.8]);
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            title(sprintf('%s Analytical Initial Guess 3D Trajectory', obj.Maneuver));
            legend({'Initial Guess', 'Launch', 'Landing'}, 'Location', 'best');

            subplot(2, 1, 2);
            plot(tv(1:end-1), Ug(3,:), 'k-', 'LineWidth', 1.5);
            grid on; xlabel('Time [s]'); ylabel('Thrust [N]');
            title('Initial Guess Commanded Thrust');
        end

        function [tbl, filepath] = exportCSV(obj, filepath)
            %% EXPORTCSV  Export trajectory table to CSV matching TOAD format.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                error('TrajectoryOptimizer:NoSolution', 'No solution available for export. Call solve() first.');
            end

            if nargin < 2 || isempty(filepath)
                if obj.SaveDir ~= ""
                    target_dir = obj.SaveDir;
                else
                    target_dir = fullfile(pwd, 'Guidance', 'Trajectories');
                end
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

            out_dir = fileparts(filepath);
            if ~exist(out_dir, 'dir') && ~isempty(out_dir)
                mkdir(out_dir);
            end
            writetable(tbl, filepath);
            fprintf('Trajectory successfully exported to:\n  %s\n', filepath);
        end
    end
end
