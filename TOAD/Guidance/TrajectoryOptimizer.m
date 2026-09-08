classdef TrajectoryOptimizer < handle
    %% TRAJECTORYOPTIMIZER  Unified 6-DoF trajectory generation and optimization engine.
    % Formulates and solves vehicle trajectories using CasADi and IPOPT, featuring
    % a decoupled 3-DoF forward guidance initial guess generator. Supports both
    % TOAD (bipropellant liquid rocket) and ASTRAv2 (electric drone prototype).
    %
    % Authors: PSP Active Controls (Pablo Plata, Andrew Lulo, & Antigravity)
    
    properties
        % System and Vehicle Properties
        constants           % Vehicle constants struct (from LoadTOADSim)
        Vehicle string      % "TOAD" or "ASTRAv2"
        isElectric logical  % True if vehicle has no propellant drain (ASTRAv2)
        
        % Discretization and Timing
        N double = 200      % Number of control intervals
        T_bounds = [15, 50] % [T_min, T_max] total trajectory duration bounds [s]
        T_initial = 35      % Initial guess for total time [s]
        
        % Maneuver Definition
        Maneuver string = "Backflip"  % "Backflip", "Circle", "Hop", "Custom"
        ManeuverParams struct         % Maneuver-specific parameters
        CustomWaypoints = []          % User-injected waypoints [N_wp x 3]
        
        % Boundary Conditions
        q0 = [1; 0; 0; 0]   % Initial attitude quaternion [w; x; y; z] (upright)
        r0 = [0; 0; 0]      % Launch pad position [m]
        v0 = [0; 0; 0]      % Launch pad velocity [m/s]
        w0 = [0; 0; 0]      % Initial angular rate [rad/s]
        r_f = [0; 0; 0]     % Landing zone position [m]
        v_f_tol = 0.1       % Max allowable landing speed [m/s]
        
        % Control and Rate Limits
        thrust_margin = 0.05
        gimbal_margin = 0.15
        max_gimbal_rate = deg2rad(30)   % Gimbal slew limit [rad/s]
        max_thrust_rate = 1000          % Thrust ramp limit [N/s]
        max_roll_rate = 4               % Roll torque limit [N*m]
        max_gimbal_angle = pi/15        % Max physical gimbal angle [rad]
        
        % Scaling Vectors
        Sx double           % State scale vector (15 x 1)
        Su double           % Input scale vector (4 x 1)
        L_c double = 50     % Characteristic position length [m]
        
        % Cost Function Weights
        w_crit = 6e-1       % Critical tilt penalty weight
        w_rate = 5e-3       % Body rate penalty weight
        w_marginGimbal = 2e-2 % Gimbal margin penalty weight
        w_qz = 6e-1         % Yaw deflection penalty weight
        w_time = 0.5        % Mission duration penalty weight
        
        % Solver Configuration
        MaxIter double = 500
        Tol double = 2e-3
        ConstrViolTol double = 1e-3
        
        % Execution & Output Flags
        PlotResults logical = false
        PlotInitialGuess logical = false
        AutoExport logical = false
        SaveFile logical = false
        Version double = 1            % Version number (formats as v001, v002, etc.)
        Filename string = ""          % Custom filename override (defaults to Vehicle_Maneuver_v###)
        SaveDir string = "sandbox/experiments/"
        
        % Storage Structures
        InitialGuess struct
        Solution struct
        OptiVars struct
    end
    
    methods
        function obj = TrajectoryOptimizer(constants6DoF, varargin)
            %% TRAJECTORYOPTIMIZER  Construct optimizer instance from vehicle constants.
            if nargin < 1 || isempty(constants6DoF)
                error('TrajectoryOptimizer requires constants6DoF struct from LoadTOADSim.');
            end
            
            % Add necessary paths
            if exist('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b', 'dir')
                addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b');
            end
            addpath(fullfile(pwd, 'Helper'));
            addpath(fullfile(pwd, 'Flight Dynamics'));
            
            obj.constants = constants6DoF;
            obj.Vehicle = constants6DoF.Vehicle;
            
            
            % Apply optional name-value pairs
            if ~isempty(varargin)
                for i = 1:2:length(varargin)
                    propName = varargin{i};
                    if isprop(obj, propName)
                        obj.(propName) = varargin{i+1};
                    end
                end
            end
            
            % Compute characteristic scales with finalized properties
            obj.updateScaling();
            
            % Set maneuver parameters with finalized N and Maneuver
            obj.setManeuver(obj.Maneuver);
        end
        
        function updateScaling(obj)
            %% UPDATESCALING  Build state and control characteristic scale vectors.
            g = obj.constants.g;
            L = obj.L_c;
            V_c = sqrt(g * L);
            W_c = sqrt(g / L);
            F_c = obj.constants.MaxThrust;
            G_c = obj.max_gimbal_angle;
            Roll_c = 10;
            
            % For electric drone, propellant mass scales default to 1 to avoid zero division
            if obj.Vehicle == 0
                Mlox_c = 1;
                Mipa_c = 1;
            else
                Mlox_c = max(obj.constants.OxMass, 1);
                Mipa_c = max(obj.constants.FuMass, 1);
            end
            
            obj.Sx = [1; 1; 1; 1; L; L; L; V_c; V_c; V_c; W_c; W_c; W_c; Mlox_c; Mipa_c];
            obj.Su = [G_c; G_c; F_c; Roll_c];
        end
        
        function setManeuver(obj, name, varargin)
            %% SETMANEUVER  Configure maneuver preset and parameter overrides.
            obj.Maneuver = string(name);
            p = struct();
            
            switch obj.Maneuver
                case "Backflip"
                    p.N_ascent = round(0.3 * obj.N);
                    p.N_flip   = round(0.5 * obj.N);
                    p.N_approach = round(0.6 * obj.N);
                    p.Glideslope = tan(deg2rad(10));
                    p.theta_tol = deg2rad(30);
                    p.q_inverted = [0; 0; -1; 0];
                    if obj.Vehicle == 0
                        p.apex_alt = 35; % Scaled for ASTRAv2 lower ceiling
                        obj.T_initial = 30;
                    else
                        p.apex_alt = 75; % TOAD full scale backflip apex
                        obj.T_initial = 35;
                    end
                    
                case "Circle"
                    p.N_ascent   = round(0.10 * obj.N);
                    p.N_c1       = round(0.25 * obj.N);
                    p.N_c2       = round(0.40 * obj.N);
                    p.N_c3       = round(0.55 * obj.N);
                    p.N_c4       = round(0.70 * obj.N);
                    p.N_approach = round(0.75 * obj.N);
                    p.circle_radius = 5.0;
                    if obj.Vehicle == 0
                        p.circle_alt = 12.0;
                    else
                        p.circle_alt = 20.0;
                    end
                    p.wp_tol = 0.5;
                    p.circle_waypoints = [ 5,  0;
                                           0,  5;
                                          -5,  0;
                                           0, -5];
                                           
                case "Hop"
                    p.N_ascent   = round(0.3 * obj.N);
                    p.N_approach = round(0.7 * obj.N);
                    p.Glideslope = tan(deg2rad(15));
                    if obj.Vehicle == 0
                        p.apex_alt = 15;
                    else
                        p.apex_alt = 30;
                    end
                    
                otherwise
                    % Custom or generic
                    p.apex_alt = 25;
            end
            
            % Override any provided parameters
            if ~isempty(varargin)
                for i = 1:2:length(varargin)
                    paramName = varargin{i};
                    p.(paramName) = varargin{i+1};
                end
            end
            obj.ManeuverParams = p;
        end
        
        function setBoundaries(obj, r0, r_f, varargin)
            %% SETBOUNDARIES  Update start and landing positions.
            obj.r0 = r0(:);
            obj.r_f = r_f(:);
            if ~isempty(varargin)
                for i = 1:2:length(varargin)
                    prop = varargin{i};
                    if isprop(obj, prop)
                        obj.(prop) = varargin{i+1};
                    end
                end
            end
        end
        
        function addWaypoint(obj, pos, tol, node_idx)
            %% ADDWAYPOINT  Inject a specific waypoint constraint at a given node.
            wp.pos = pos(:);
            wp.tol = tol;
            wp.node_idx = node_idx;
            obj.CustomWaypoints = [obj.CustomWaypoints; wp];
        end
        
        function guess = generateInitialGuess(obj)
            %% GENERATEINITIALGUESS  Decoupled 3-DoF forward guidance initial guess generator.
            % Simulates translational point-mass motion forward in time using proportional
            % waypoint guidance, then algebraically reconstructs consistent attitude
            % quaternions, angular rates, thrust, and mass consumption profiles.
            
            N_nodes = obj.N + 1;
            N_ctrl  = obj.N;
            T_sim   = obj.T_initial;
            dt_sim  = T_sim / obj.N;
            time_vec = linspace(0, T_sim, N_nodes);
            
            % Initial mass state
            m_dry = obj.constants.m_dry;
            g = obj.constants.g;
            MaxThrust = obj.constants.MaxThrust;
            min_thrust = (0.25 + obj.thrust_margin) * MaxThrust;
            max_thrust = (1.0 - obj.thrust_margin) * MaxThrust;
            
            if obj.Vehicle == 0
                m_curr = m_dry;
                m_lox_sim = zeros(1, N_nodes);
                m_ipa_sim = zeros(1, N_nodes);
            else
                m_curr = obj.constants.m_wet;
                m_lox0 = obj.constants.OxMass;
                m_ipa0 = obj.constants.FuMass;
                m_lox_sim = linspace(m_lox0, 0.20 * m_lox0, N_nodes);
                m_ipa_sim = linspace(m_ipa0, 0.20 * m_ipa0, N_nodes);
            end
            
            % 1. Kinematic Position and Velocity directly from Guidance Schedule
            r_sim = zeros(3, N_nodes);
            v_sim = zeros(3, N_nodes);
            for k = 1:N_nodes
                frac = (k - 1) / obj.N;
                [r_tgt, v_tgt] = obj.getGuidanceTarget(frac);
                r_sim(:, k) = r_tgt;
                v_sim(:, k) = v_tgt;
            end
            r_sim(:, 1)   = obj.r0;
            v_sim(:, 1)   = obj.v0;
            r_sim(:, end) = obj.r_f;
            v_sim(:, end) = [0; 0; 0];
            
            % 2. Attitude, Angular Rates, and Control Allocation
            quat_sim  = zeros(4, N_nodes);
            omega_sim = zeros(3, N_nodes);
            U_phys    = zeros(4, N_ctrl);
            
            if obj.Maneuver == "Backflip"
                k_flip = obj.ManeuverParams.N_flip;
                k_span = max(4, round(3.5 / dt_sim));
                k_flip_start = max(1, k_flip - k_span);
                k_flip_end   = min(N_ctrl, k_flip + k_span);
                L1 = max(1, k_flip - k_flip_start);
                L2 = max(1, k_flip_end - k_flip);
                T1 = L1 * dt_sim;
                T2 = L2 * dt_sim;
                
                % Prescribe continuous smooth pitch profile 0 -> -pi -> 0
                for k = 1:N_nodes
                    if k <= k_flip_start
                        quat_sim(:, k)  = obj.q0;
                        omega_sim(:, k) = [0; 0; 0];
                    elseif k <= k_flip
                        tau1 = (k - k_flip_start) / L1;
                        th_p = -pi * sin(0.5 * pi * tau1)^2;
                        dth_p = -(0.5 * pi^2 / T1) * sin(pi * tau1);
                        quat_sim(:, k)  = [cos(th_p / 2); 0; sin(th_p / 2); 0];
                        omega_sim(:, k) = [0; dth_p; 0];
                    elseif k <= k_flip_end
                        tau2 = (k - k_flip) / L2;
                        th_p = -pi * cos(0.5 * pi * tau2)^2;
                        dth_p = (0.5 * pi^2 / T2) * sin(pi * tau2);
                        quat_sim(:, k)  = [cos(th_p / 2); 0; sin(th_p / 2); 0];
                        omega_sim(:, k) = [0; dth_p; 0];
                    else
                        quat_sim(:, k)  = obj.q0;
                        omega_sim(:, k) = [0; 0; 0];
                    end
                end
                
                % Compute control inputs
                for k = 1:N_ctrl
                    if k > k_flip_start && k <= k_flip_end
                        % During flip: minimum thrust and pitch gimbal torque
                        T_cmd = min_thrust;
                        J_yy = obj.constants.J(2, 2);
                        if k <= k_flip
                            tau1 = (k - k_flip_start) / L1;
                            ddth_p = -(0.5 * pi^3 / (T1^2)) * cos(pi * tau1);
                        else
                            tau2 = (k - k_flip) / L2;
                            ddth_p = (0.5 * pi^3 / (T2^2)) * cos(pi * tau2);
                        end
                        M_pitch = J_yy * ddth_p;
                        gimbal_phi = -M_pitch / (max(obj.constants.rTB * T_cmd, 1e-2));
                        max_g = (1 - obj.gimbal_margin) * obj.max_gimbal_angle;
                        gimbal_phi = max(min(gimbal_phi, max_g), -max_g);
                        U_phys(:, k) = [0; gimbal_phi; T_cmd; 0];
                    else
                        % Outside flip: vertical thrust
                        a_acc = (v_sim(:, k+1) - v_sim(:, k)) / dt_sim;
                        F_req = m_curr * (a_acc + [0; 0; g]);
                        T_cmd = max(min(F_req(3), max_thrust), min_thrust);
                        U_phys(:, k) = [0; 0; T_cmd; 0];
                    end
                end
                
            else
                % Hop, Circle, and General maneuvers (Differential Flatness)
                for k = 1:N_ctrl
                    a_acc = (v_sim(:, k+1) - v_sim(:, k)) / dt_sim;
                    F_req = m_curr * (a_acc + [0; 0; g]);
                    F_mag = norm(F_req);
                    T_cmd = max(min(F_mag, max_thrust), min_thrust);
                    U_phys(:, k) = [0; 0; T_cmd; 0];
                    quat_sim(:, k) = obj.vectorToQuat(F_req);
                end
                quat_sim(:, end) = obj.q0;
                
                for k = 1:N_ctrl
                    dq = (quat_sim(:, k+1) - quat_sim(:, k)) / dt_sim;
                    qw = quat_sim(1, k); qv = quat_sim(2:4, k);
                    omega_sim(:, k) = 2 * (qw * dq(2:4) - dq(1) * qv - cross(qv, dq(2:4)));
                end
                omega_sim(:, end) = [0; 0; 0];
                
                % Reconstruct gimbal angles matching required torque MB = J * domega/dt
                max_g = (1 - obj.gimbal_margin) * obj.max_gimbal_angle;
                for k = 1:N_ctrl
                    domega = (omega_sim(:, k+1) - omega_sim(:, k)) / dt_sim;
                    J_xx = obj.constants.J(1, 1);
                    J_yy = obj.constants.J(2, 2);
                    lever = max(obj.constants.rTB * U_phys(3, k), 1e-2);
                    g_theta = -J_xx * domega(1) / lever;
                    g_phi   = -J_yy * domega(2) / lever;
                    U_phys(1, k) = max(min(g_theta, max_g), -max_g);
                    U_phys(2, k) = max(min(g_phi,   max_g), -max_g);
                end
            end
            
            quat_sim(:, 1)   = obj.q0;
            quat_sim(:, end) = obj.q0;
            omega_sim(:, 1)   = [0; 0; 0];
            omega_sim(:, end) = [0; 0; 0];
            
            % Assemble physical state array (15 x N_nodes)
            X_phys = [quat_sim;
                      r_sim;
                      v_sim;
                      omega_sim;
                      m_lox_sim;
                      m_ipa_sim];
                      
            % Nondimensionalize
            Xhat_guess = X_phys ./ obj.Sx;
            Uhat_guess = U_phys ./ obj.Su;
            
            % Package initial guess struct
            guess.Time = time_vec;
            guess.X = X_phys;
            guess.U = U_phys;
            guess.T_total = T_sim;
            guess.Xhat = Xhat_guess;
            guess.Uhat = Uhat_guess;
            
            obj.InitialGuess = guess;
            
            if obj.PlotInitialGuess
                obj.plotInitialGuess();
            end
        end
        
        function [r_tgt, v_tgt] = getGuidanceTarget(obj, frac)
            %% GETGUIDANCETARGET  Evaluate target waypoint for current normalized mission time.
            p = obj.ManeuverParams;
            r0_pad = obj.r0;
            rf_pad = obj.r_f;
            
            switch obj.Maneuver
                case "Backflip"
                    f_flip = p.N_flip / obj.N;
                    apex = [0.5 * (r0_pad(1) + rf_pad(1));
                            0.5 * (r0_pad(2) + rf_pad(2));
                            p.apex_alt];
                    if frac < f_flip
                        % Smooth harmonic ascent to apex: zero vertical speed at apex
                        sub_f = frac / f_flip;
                        s_pos = sin(0.5 * pi * sub_f)^2;
                        s_vel = (pi / (2 * f_flip * obj.T_initial)) * sin(pi * sub_f);
                        r_tgt = r0_pad + (apex - r0_pad) * s_pos;
                        v_tgt = (apex - r0_pad) * s_vel;
                    else
                        % Smooth harmonic descent from apex: zero speed at touchdown
                        sub_f = (frac - f_flip) / (1 - f_flip);
                        s_pos = cos(0.5 * pi * sub_f)^2;
                        s_vel = -(pi / (2 * (1 - f_flip) * obj.T_initial)) * sin(pi * sub_f);
                        r_tgt = rf_pad + (apex - rf_pad) * s_pos;
                        v_tgt = (apex - rf_pad) * s_vel;
                    end
                    
                case "Circle"
                    c_alt = p.circle_alt;
                    rad   = p.circle_radius;
                    f_asc = p.N_ascent / obj.N;
                    f_c1  = p.N_c1 / obj.N;
                    f_c4  = p.N_c4 / obj.N;
                    f_app = p.N_approach / obj.N;
                    T_init = obj.T_initial;
                    
                    % Angular velocity in circular arc (0 -> 1.5*pi)
                    omega_c = 1.5 * pi / ((f_c4 - f_c1) * T_init);
                    
                    if frac < f_asc
                        % Phase 1: Smooth vertical ascent to c_alt within 1m pad cylinder
                        s = frac / f_asc;
                        t_seg = f_asc * T_init;
                        r_start = r0_pad;
                        v_start = [0; 0; 0];
                        r_end   = [r0_pad(1); r0_pad(2); c_alt];
                        v_end   = [0; 0; 0];
                        [r_tgt, v_tgt] = obj.evalHermite(r_start, v_start, r_end, v_end, s, t_seg);
                        
                    elseif frac < f_c1
                        % Phase 2: Smooth horizontal ingress to Quadrant 1 (rad, 0, c_alt)
                        s = (frac - f_asc) / (f_c1 - f_asc);
                        t_seg = (f_c1 - f_asc) * T_init;
                        r_start = [r0_pad(1); r0_pad(2); c_alt];
                        v_start = [0; 0; 0];
                        r_end   = [rad; 0; c_alt];
                        v_end   = [0; rad * omega_c; 0];
                        [r_tgt, v_tgt] = obj.evalHermite(r_start, v_start, r_end, v_end, s, t_seg);
                        
                    elseif frac < f_c4
                        % Phase 3: Smooth circular orbit from quadrant 1 to 4 (0 to 1.5*pi)
                        tau = (frac - f_c1) / (f_c4 - f_c1);
                        theta_c = 1.5 * pi * tau;
                        r_tgt = [rad * cos(theta_c); rad * sin(theta_c); c_alt];
                        v_tgt = [-rad * omega_c * sin(theta_c); rad * omega_c * cos(theta_c); 0];
                        
                    elseif frac < f_app
                        % Phase 4: Smooth horizontal egress to descent column over landing pad
                        s = (frac - f_c4) / (f_app - f_c4);
                        t_seg = (f_app - f_c4) * T_init;
                        r_start = [0; -rad; c_alt];
                        v_start = [rad * omega_c; 0; 0];
                        r_end   = [rf_pad(1); rf_pad(2); c_alt];
                        v_end   = [0; 0; 0];
                        [r_tgt, v_tgt] = obj.evalHermite(r_start, v_start, r_end, v_end, s, t_seg);
                        
                    else
                        % Phase 5: Smooth vertical descent to landing pad with zero touchdown speed
                        s = (frac - f_app) / (1 - f_app);
                        t_seg = (1 - f_app) * T_init;
                        r_start = [rf_pad(1); rf_pad(2); c_alt];
                        v_start = [0; 0; 0];
                        r_end   = rf_pad;
                        v_end   = [0; 0; 0];
                        [r_tgt, v_tgt] = obj.evalHermite(r_start, v_start, r_end, v_end, s, t_seg);
                    end
                    
                otherwise
                    % Hop or standard
                    apex = [0.5*(r0_pad(1)+rf_pad(1)); 0.5*(r0_pad(2)+rf_pad(2)); p.apex_alt];
                    if frac < 0.5
                        r_tgt = r0_pad + (apex - r0_pad) * (frac / 0.5);
                        v_tgt = (apex - r0_pad) / (0.5 * obj.T_initial);
                    else
                        r_tgt = apex + (rf_pad - apex) * ((frac - 0.5) / 0.5);
                        v_tgt = (rf_pad - apex) / (0.5 * obj.T_initial);
                    end
            end
        end
        
        function [r, v] = evalHermite(~, r0, v0, r1, v1, s, dt_seg)
            %% EVALHERMITE  Evaluate cubic Hermite polynomial position and velocity.
            s = max(0, min(1, s));
            dt_seg = max(dt_seg, 1e-3);
            h00 = (2*s^3 - 3*s^2 + 1);
            h10 = (s^3 - 2*s^2 + s);
            h01 = (-2*s^3 + 3*s^2);
            h11 = (s^3 - s^2);
            dh00 = (6*s^2 - 6*s);
            dh10 = (3*s^2 - 4*s + 1);
            dh01 = (-6*s^2 + 6*s);
            dh11 = (3*s^2 - 2*s);
            r = r0*h00 + dt_seg*v0*h10 + r1*h01 + dt_seg*v1*h11;
            v = (r0*dh00 + dt_seg*v0*dh10 + r1*dh01 + dt_seg*v1*dh11) / dt_seg;
        end
        
        function q = vectorToQuat(~, vec)
            %% VECTORTOQUAT  Construct quaternion aligning body z-axis [0;0;1] with vector.
            n_vec = norm(vec);
            if n_vec < 1e-4
                q = [1; 0; 0; 0];
                return;
            end
            z_target = vec(:) / n_vec;
            z_body = [0; 0; 1];
            
            cross_v = cross(z_body, z_target);
            dot_v   = dot(z_body, z_target);
            
            if dot_v < -0.9999
                % Exactly opposing vectors: 180 deg rotation about X
                q = [0; 1; 0; 0];
            else
                qw = 1 + dot_v;
                qv = cross_v;
                q = [qw; qv];
                q = q / norm(q);
            end
        end
        
        function opti = buildOptimizationProblem(obj)
            %% BUILDOPTIMIZATIONPROBLEM  Assemble CasADi Opti variables and constraints.
            import casadi.*
            opti = casadi.Opti();
            N = obj.N;
            
            % Generate forward dynamics function
            dyn_fnc = obj.getCasADiDynamics();
            
            % Optimization Variables
            T_total = opti.variable();
            opti.subject_to(obj.T_bounds(1) <= T_total);
            opti.subject_to(T_total <= obj.T_bounds(2));
            dt = T_total / N;
            
            Xhat = opti.variable(15, N+1);
            Uhat = opti.variable(4, N);
            
            % Physical Unit Aliases
            X = obj.Sx .* Xhat;
            U = obj.Su .* Uhat;
            
            % System Parameter Vector for Dynamics Call
            params_val = [
                obj.constants.m_dry;
                obj.constants.g;
                obj.constants.rTB;
                obj.constants.Ox_Z;
                obj.constants.OxMass;
                obj.constants.OxHeight;
                obj.constants.Fu_Z;
                obj.constants.FuMass;
                obj.constants.FuHeight;
                obj.constants.J(:);
                obj.constants.OxRadius;
                obj.constants.FuRadius;
                obj.constants.MaxThrust;
                obj.constants.OF;
                obj.constants.MaxMdot;
                0;              % MaxMdot_d
                zeros(9, 1);    % J_d_vec
                zeros(3, 1)     % TB_d_val
            ];
            
            % Single-step RK4 Integrator Function
            xhat_sym = MX.sym('xhat', 15);
            uhat_sym = MX.sym('uhat', 4);
            dt_sym   = MX.sym('dt');
            
            x_sym = obj.Sx .* xhat_sym;
            u_sym = obj.Su .* uhat_sym;
            
            k1 = dyn_fnc(x_sym,                u_sym, params_val);
            k2 = dyn_fnc(x_sym + dt_sym/2*k1,  u_sym, params_val);
            k3 = dyn_fnc(x_sym + dt_sym/2*k2,  u_sym, params_val);
            k4 = dyn_fnc(x_sym + dt_sym*k3,    u_sym, params_val);
            x_next_sym = x_sym + dt_sym/6 * (k1 + 2*k2 + 2*k3 + k4);
            q_next_unit = x_next_sym(1:4) / sqrt(sum(x_next_sym(1:4).^2) + 1e-12);
            x_next_sym = [q_next_unit; x_next_sym(5:end)];
            
            xhat_next_sym = x_next_sym ./ obj.Sx;
            F_step = Function('F_step', {xhat_sym, uhat_sym, dt_sym}, {xhat_next_sym});
            
            % Vectorized RK4 Constraint across N intervals
            F_map = F_step.map(N);
            dt_row = repmat(dt, 1, N);
            Xhat_next_all = F_map(Xhat(:, 1:N), Uhat, dt_row);
            opti.subject_to(Xhat(:, 2:end) == Xhat_next_all);
            
            % Quaternion Unit Norm Constraint
            opti.subject_to(sum(X(1:4, :).^2, 1) == 1.00);
            
            % Flight Sandbox Constraints
            opti.subject_to(-30 <= X(5, :) <= 30); %#ok<CHAIN>
            opti.subject_to(-30 <= X(6, :) <= 30); %#ok<CHAIN>
            if obj.Vehicle == 0
                opti.subject_to(-1 <= X(7, :) <= 75); %#ok<CHAIN>  % ASTRAv2 ceiling
            else
                opti.subject_to(-1 <= X(7, :) <= 150); %#ok<CHAIN> % TOAD ceiling
            end
            
            % Initial State Constraints (On the Launch Pad)
            m_lox0 = obj.constants.OxMass;
            m_ipa0 = obj.constants.FuMass;
            opti.subject_to(X(:, 1) == [obj.q0; obj.r0; obj.v0; obj.w0; m_lox0; m_ipa0]);
            
            % Final State Constraints (Landing Zone)
            opti.subject_to(X(1:4, end) == obj.q0);
            opti.subject_to(X(5:7, end) == obj.r_f);
            opti.subject_to(sum(X(8:10, end).^2) <= obj.v_f_tol^2);
            
            % Propellant Margin (Enforced only for chemical rocket with drain dynamics)
            if ~(obj.Vehicle == 0)
                prop_margin_frac = 0.10;
                opti.subject_to(X(14, end) >= prop_margin_frac * m_lox0);
                opti.subject_to(X(15, end) >= prop_margin_frac * m_ipa0);
            else
                % Electric drone keeps mass states zero
                opti.subject_to(X(14:15, :) == 0);
            end
            
            % Control Input Bounds
            MaxThrust = obj.constants.MaxThrust;
            t_margin = obj.thrust_margin;
            g_margin = obj.gimbal_margin;
            max_gimbal = obj.max_gimbal_angle;
            
            opti.subject_to((0.25 + t_margin) * MaxThrust <= U(3, :) <= (1 - t_margin) * MaxThrust); %#ok<CHAIN>
            opti.subject_to(-(1 - g_margin) * max_gimbal <= U(1, :) <= (1 - g_margin) * max_gimbal); %#ok<CHAIN>
            opti.subject_to(-(1 - g_margin) * max_gimbal <= U(2, :) <= (1 - g_margin) * max_gimbal); %#ok<CHAIN>
            opti.subject_to(-(1 - t_margin) * obj.max_roll_rate <= U(4, :) <= (1 - t_margin) * obj.max_roll_rate); %#ok<CHAIN>
            
            % Control Rate Limits
            dU_phys = U(:, 2:end) - U(:, 1:end-1);
            opti.subject_to(-obj.max_gimbal_rate * dt <= dU_phys(1, :) <= obj.max_gimbal_rate * dt); %#ok<CHAIN>
            opti.subject_to(-obj.max_gimbal_rate * dt <= dU_phys(2, :) <= obj.max_gimbal_rate * dt); %#ok<CHAIN>
            opti.subject_to(-obj.max_thrust_rate * dt <= dU_phys(3, :) <= obj.max_thrust_rate * dt); %#ok<CHAIN>
            opti.subject_to(-obj.max_roll_rate * dt   <= dU_phys(4, :) <= obj.max_roll_rate * dt); %#ok<CHAIN>
            
            % Maneuver-Specific Constraints
            obj.applyManeuverConstraints(opti, X);
            
            % Objective Function Formulation
            obj.applyCostFunction(opti, Xhat, Uhat, T_total);
            
            % Apply Initial Guess
            if isempty(obj.InitialGuess)
                obj.generateInitialGuess();
            end
            opti.set_initial(T_total, obj.InitialGuess.T_total);
            opti.set_initial(Xhat, obj.InitialGuess.Xhat);
            opti.set_initial(Uhat, obj.InitialGuess.Uhat);
            
            % Store optimization variables on object
            obj.OptiVars = struct('Xhat', Xhat, 'Uhat', Uhat, 'T_total', T_total);
        end
        
        function applyManeuverConstraints(obj, opti, X)
            %% APPLYMANEUVERCONSTRAINTS  Inject constraints tailored to chosen maneuver.
            p = obj.ManeuverParams;
            
            switch obj.Maneuver
                case "Backflip"
                    N_ascent   = p.N_ascent;
                    N_flip     = p.N_flip;
                    N_approach = p.N_approach;
                    Glideslope = p.Glideslope;
                    
                    % Ascent glideslope cone (Second-order cone formulation)
                    pos_asc = X(5:7, 1:N_ascent);
                    opti.subject_to(X(10, 1:N_ascent) >= 0);
                    opti.subject_to(sqrt(pos_asc(1,:).^2 + pos_asc(2,:).^2 + 1e-4) <= pos_asc(3,:) * Glideslope + 0.05);
                    
                    % Flip Maneuver Target Attitude (Convex linear inner product)
                    att_tol = cos(p.theta_tol / 2);
                    q_flip = X(1:4, N_flip);
                    opti.subject_to(p.q_inverted' * q_flip >= att_tol);
                    
                    % Descent glideslope cone (Second-order cone formulation)
                    pos_desc = X(5:7, N_approach:end);
                    opti.subject_to(X(10, N_approach:end) <= 0.1);
                    opti.subject_to(sqrt((pos_desc(1,:) - obj.r_f(1)).^2 + (pos_desc(2,:) - obj.r_f(2)).^2 + 1e-4) <= ...
                                    (pos_desc(3,:) - obj.r_f(3)) * Glideslope + 0.5);
                                    
                case "Circle"
                    % Ascent
                    opti.subject_to(X(1:4, p.N_ascent) == obj.q0);
                    opti.subject_to(sum(X(5:6, 1:p.N_ascent).^2, 1) <= 1.0);
                    
                    % Quadrant waypoints (normalized)
                    circle_nodes = [p.N_c1, p.N_c2, p.N_c3, p.N_c4];
                    for i = 1:4
                        k = circle_nodes(i);
                        tgt = p.circle_waypoints(i, :)';
                        opti.subject_to(((X(5, k) - tgt(1))^2 + (X(6, k) - tgt(2))^2) / p.wp_tol^2 <= 1.0);
                    end
                    
                    % Altitude and radius boundaries during circle orbit (normalized)
                    alt_band = 3.0;
                    rad_band = 2.0;
                    for k = p.N_c1:p.N_c4
                        opti.subject_to(X(7, k) >= (p.circle_alt - alt_band));
                        opti.subject_to(X(7, k) <= (p.circle_alt + alt_band));
                        opti.subject_to(((p.circle_radius - rad_band) / p.circle_radius)^2 <= ...
                                        (X(5,k)^2 + X(6,k)^2) / p.circle_radius^2);
                        opti.subject_to((X(5,k)^2 + X(6,k)^2) / p.circle_radius^2 <= ...
                                        ((p.circle_radius + rad_band) / p.circle_radius)^2);
                    end
                    
                    % Descent
                    opti.subject_to(X(1:4, p.N_approach) == obj.q0);
                    opti.subject_to(sum((X(5:6, p.N_approach:end) - obj.r_f(1:2)).^2, 1) <= 1.0);
                    opti.subject_to(-2 <= X(10, p.N_approach:end));
                    opti.subject_to(X(10, p.N_approach:end) <= 2);
                    
                case "Hop"
                    N_ascent = p.N_ascent;
                    N_approach = p.N_approach;
                    
                    % Ascent & Descent vertical speed constraints
                    opti.subject_to(X(10, 1:N_ascent) >= 0);
                    opti.subject_to(X(10, N_approach:end) <= 0.1);
                    
                    % Keep attitude upright within 30 deg of vertical
                    R33 = X(1,:).^2 - X(2,:).^2 - X(3,:).^2 + X(4,:).^2;
                    opti.subject_to(R33 >= cosd(35));
            end
            
            % Custom Waypoint Injection
            for i = 1:length(obj.CustomWaypoints)
                wp = obj.CustomWaypoints(i);
                opti.subject_to(sum((X(5:7, wp.node_idx) - wp.pos).^2) <= wp.tol^2);
            end
        end
        
        function applyCostFunction(obj, opti, Xhat, Uhat, T_total)
            %% APPLYCOSTFUNCTION  Formulate multi-objective survivability and control cost.
            % Penalize angular rates across all maneuvers
            J_rate = sum(sum(Xhat(11:13, 1:end-1).^2, 1));
            
            % Control margins
            gimbal_bound = (1 - obj.gimbal_margin);
            J_marginGimbal = sum(sum((Uhat(1:2, :) / gimbal_bound).^2));
            J_qz = sum((Xhat(4, :)).^2);
            
            opti.minimize( ...
                1.0                * J_marginGimbal  + ...
                obj.w_rate         * J_rate          + ...
                obj.w_qz           * J_qz            + ...
                obj.w_time         * T_total);
        end
        
        function sol = solve(obj)
            %% SOLVE  Execute CasADi IPOPT optimization.
            opti = obj.buildOptimizationProblem();
            
            p_opts = struct('expand', true);
            s_opts = struct('max_iter', obj.MaxIter, ...
                            'tol', obj.Tol, ...
                            'constr_viol_tol', obj.ConstrViolTol, ...
                            'acceptable_tol', 1e-2, ...
                            'acceptable_constr_viol_tol', 1e-3, ...
                            'acceptable_iter', 10, ...
                            'print_level', 0);
            opti.solver('ipopt', p_opts, s_opts);
            
            Xhat = obj.OptiVars.Xhat;
            Uhat = obj.OptiVars.Uhat;
            T_total = obj.OptiVars.T_total;
            
            fprintf('Starting %s optimization solve for vehicle: %s\n', obj.Maneuver, obj.Vehicle);
            
            try
                sol_casadi = opti.solve();
                disp('Optimal trajectory found!');
                status = 'Success';
                
                X_res = obj.Sx .* sol_casadi.value(Xhat);
                U_res = obj.Su .* sol_casadi.value(Uhat);
                T_res = sol_casadi.value(T_total);
                sol_obj = sol_casadi;
                
            catch ME
                disp('Solver did not achieve full convergence. Retrieving debug trajectory...');
                status = 'Debug';
                X_res = obj.Sx .* opti.debug.value(Xhat);
                U_res = obj.Su .* opti.debug.value(Uhat);
                T_res = opti.debug.value(T_total);
                sol_obj = ME;
            end
            
            time_res = linspace(0, T_res, obj.N + 1);
            
            % Store in Solution struct
            obj.Solution = struct( ...
                'Time', time_res, ...
                'X', X_res, ...
                'U', U_res, ...
                'T_total', T_res, ...
                'Status', status, ...
                'sol', sol_obj);
            sol = obj.Solution;
            
            if obj.PlotResults
                obj.plot();
            end
            
            if obj.AutoExport
                out_path = fullfile(pwd, obj.SaveDir, obj.getFormattedFilename() + ".csv");
                obj.exportCSV(out_path);
            end
        end
        
        function dyn_fnc = getCasADiDynamics(obj)
            %% GETCASADIDYNAMICS  Construct or retrieve symbolic CasADi dynamics function.
            import casadi.*
            
            % Symbolic states
            q = MX.sym('q', 4);
            r = MX.sym('r', 3);
            v = MX.sym('v', 3);
            omegaB = MX.sym('omegaB', 3);
            m_lox = MX.sym('m_lox', 1);
            m_ipa = MX.sym('m_ipa', 1);
            
            % Symbolic parameters
            m_dry = MX.sym('m_dry', 1);
            g = MX.sym('g');
            rTB = MX.sym('rTB');
            Ox_Z = MX.sym('Ox_Z');
            OxMassI = MX.sym('OxMassI');
            OxHeight = MX.sym('OxHeight');
            Fu_Z = MX.sym('Fu_Z');
            FuMassI = MX.sym('FuMassI');
            FuHeight = MX.sym('FuHeight');
            J = MX.sym('J', 3, 3);
            OxRadius = MX.sym('OxRadius');
            FuRadius = MX.sym('FuRadius');
            MaxThrust = MX.sym('MaxThrust');
            OF = MX.sym('OF');
            MaxMdot = MX.sym('MaxMdot');
            MaxMdot_d = MX.sym('MaxMdot_d');
            J_d = MX.sym('J_d', 3, 3);
            TB_d = MX.sym('TB_d', 3);
            
            % Symbolic controls
            theta = MX.sym('theta', 1);
            phi = MX.sym('phi', 1);
            thrust = MX.sym('thrust', 1);
            roll = MX.sym('roll', 1);
            
            m = m_dry + m_lox + m_ipa;
            C_BI = quatRot(q);
            C_IB = C_BI.';
            
            TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            FI = C_IB * TB + [0; 0; -m*g];
            
            rdot = v;
            vdot = FI / m;
            
            % Propellant Drain Dynamics
            if obj.Vehicle == 0
                mdot_lox = 0;
                mdot_ipa = 0;
                OxFluidHeight = 0;
                FuFluidHeight = 0;
            else
                mdot_lox = -thrust / MaxThrust * OF / (1 + OF) * (MaxMdot + MaxMdot_d);
                mdot_ipa = -thrust / MaxThrust * 1 / (1 + OF) * (MaxMdot + MaxMdot_d);
                OxFluidHeight = (m_lox / OxMassI) * OxHeight * 0.9;
                FuFluidHeight = (m_ipa / FuMassI) * FuHeight * 0.9;
            end
            
            % Propellant Inertias
            J_xx = 1/12 * m_lox * (3 * OxRadius^2 + OxFluidHeight^2);
            J_zz = 1/2 * m_lox * OxRadius^2;
            J_lox = [J_xx, 0, 0; 0, J_xx, 0; 0, 0, J_zz];
            
            J_xx_fu = 1/12 * m_ipa * (3 * FuRadius^2 + FuFluidHeight^2);
            J_zz_fu = 1/2 * m_ipa * FuRadius^2;
            J_ipa = [J_xx_fu, 0, 0; 0, J_xx_fu, 0; 0, 0, J_zz_fu];
            
            % Centers of Mass
            OxFluidLocation = Ox_Z + OxFluidHeight / 2;
            FuFluidLocation = Fu_Z + FuFluidHeight / 2;
            CGz = (m_dry * rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;
            
            d_dry = rTB - CGz + TB_d(3);
            d_lox = OxFluidLocation - CGz;
            d_ipa = FuFluidLocation - CGz;
            
            J_dry = J + m_dry * diag([d_dry^2, d_dry^2, 0]);
            J_lox_tot = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
            J_ipa_tot = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);
            J_tot = J_dry + J_lox_tot + J_ipa_tot + J_d;
            
            % Body Moments
            thrustDir = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
            if obj.Vehicle == 0
                MB = zetaCross([0; 0; -CGz] + TB_d)*TB + roll * thrustDir;
            else
                MB = zetaCross([0; 0; -CGz] + TB_d)*TB + [0; 0; roll];
            end
            
            qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
            omegaBdot = (MB - zetaCross(omegaB) * J_tot * omegaB) ./ diag(J_tot);
            
            x = [q; r; v; omegaB; m_lox; m_ipa];
            u = [theta; phi; thrust; roll];
            params = [m_dry; g; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight;
                      J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot; MaxMdot_d; J_d(:); TB_d];
            xdot = [qdot; rdot; vdot; omegaBdot; mdot_lox; mdot_ipa];
            
            dyn_fnc = Function('dynamics_fnc', {x, u, params}, {xdot});
        end
        
        function figHandles = plot(obj)
            %% PLOT  Generate visualization figures matching the TOAD mission profile.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                warning('No solution available to plot. Call solve() first.');
                return;
            end
            
            t_state = obj.Solution.Time;
            t_ctrl  = linspace(0, obj.Solution.T_total - (obj.Solution.T_total/obj.N), obj.N);
            X_res   = obj.Solution.X;
            U_res   = obj.Solution.U;
            
            pos = X_res(5:7, :);
            vel = X_res(8:10, :);
            quat = X_res(1:4, :);
            
            theta_cmd = U_res(1, :);
            phi_cmd   = U_res(2, :);
            thrust    = U_res(3, :);
            
            % Plot 1: 3D Mission Trajectory
            f1 = figure('Name', sprintf('%s: 3D Mission Profile (%s)', obj.Maneuver, obj.Vehicle));
            tl = tiledlayout(f1, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            axMain = nexttile(tl, 1, [3 3]);
            hold(axMain, 'on'); grid(axMain, 'on'); axis(axMain, 'equal'); view(axMain, 3);
            xlabel(axMain, 'X [m]'); ylabel(axMain, 'Y [m]'); zlabel(axMain, 'Z (Alt) [m]');
            title(axMain, sprintf('3D Trajectory (%s - %s)', obj.Maneuver, obj.Vehicle));
            
            patch(axMain, [pos(1,:), NaN], [pos(2,:), NaN], [pos(3,:), NaN], [t_state, NaN], ...
                  'FaceColor', 'none', 'EdgeColor', 'interp', 'LineWidth', 2.5);
            cb = colorbar(axMain); cb.Label.String = 'Time [s]'; colormap(axMain, 'turbo');
            
            plot3(axMain, obj.r0(1), obj.r0(2), obj.r0(3), 'gs', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
            plot3(axMain, obj.r_f(1), obj.r_f(2), obj.r_f(3), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
            
            axTop = nexttile(tl, 4);
            hold(axTop, 'on'); grid(axTop, 'on'); axis(axTop, 'equal');
            plot(axTop, pos(1,:), pos(2,:), 'b', 'LineWidth', 1.5);
            xlabel(axTop, 'X [m]'); ylabel(axTop, 'Y [m]'); title(axTop, 'Top View (X-Y)');
            
            axSide = nexttile(tl, 8);
            hold(axSide, 'on'); grid(axSide, 'on'); axis(axSide, 'equal');
            plot(axSide, pos(1,:), pos(3,:), 'b', 'LineWidth', 1.5);
            xlabel(axSide, 'X [m]'); ylabel(axSide, 'Z [m]'); title(axSide, 'Side View (X-Z)');
            
            axFront = nexttile(tl, 12);
            hold(axFront, 'on'); grid(axFront, 'on'); axis(axFront, 'equal');
            plot(axFront, pos(2,:), pos(3,:), 'b', 'LineWidth', 1.5);
            xlabel(axFront, 'Y [m]'); ylabel(axFront, 'Z [m]'); title(axFront, 'Front View (Y-Z)');
            
            % Plot 2: Control & Attitude Performance
            f2 = figure('Name', sprintf('%s: Control & Attitude (%s)', obj.Maneuver, obj.Vehicle));
            subplot(3, 1, 1); hold on; grid on;
            yyaxis left;
            plot(t_state, quat(1, :), 'k', 'LineWidth', 1.5, 'DisplayName', 'q_w');
            plot(t_state, quat(2, :), 'r', 'LineWidth', 1.5, 'DisplayName', 'q_x');
            plot(t_state, quat(3, :), 'g', 'LineWidth', 1.5, 'DisplayName', 'q_y');
            plot(t_state, quat(4, :), 'b', 'LineWidth', 1.5, 'DisplayName', 'q_z');
            ylabel('Quaternions'); ylim([-1.1, 1.1]);
            
            yyaxis right;
            R33 = quat(1,:).^2 - quat(2,:).^2 - quat(3,:).^2 + quat(4,:).^2;
            tilt = acosd(max(min(R33, 1), -1));
            plot(t_state, tilt, '--m', 'LineWidth', 2, 'DisplayName', 'Tilt (deg)');
            ylabel('Tilt [deg]'); ylim([0, 185]);
            legend('Location', 'eastoutside'); title('Attitude Tracking');
            
            subplot(3, 1, 2); hold on; grid on;
            plot(t_ctrl, rad2deg(theta_cmd), 'Color', '#7E2F8E', 'LineWidth', 1.5);
            plot(t_ctrl, rad2deg(phi_cmd), 'Color', '#77AC30', 'LineWidth', 1.5);
            ylabel('Gimbal [deg]'); legend('\theta', '\phi', 'Location', 'best');
            title('Gimbal Commands');
            
            subplot(3, 1, 3); hold on; grid on;
            area(t_ctrl, thrust, 'FaceColor', '#EDB120', 'FaceAlpha', 0.4, 'EdgeColor', '#A27712', 'LineWidth', 1.5);
            ylabel('Thrust [N]'); xlabel('Time [s]'); title('Thrust Command');
            
            % Plot 3: Mission Kinematics
            f3 = figure('Name', sprintf('%s: Kinematics (%s)', obj.Maneuver, obj.Vehicle));
            tl_kin = tiledlayout(f3, 3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            labels_p = {'X Pos [m]', 'Y Pos [m]', 'Z Pos [m]'};
            labels_v = {'X Vel [m/s]', 'Y Vel [m/s]', 'Z Vel [m/s]'};
            cols = {'#0072BD', '#D95319', '#EDB120'};
            for i = 1:3
                nexttile(tl_kin, (i-1)*2 + 1); hold on; grid on;
                plot(t_state, pos(i, :), 'Color', cols{i}, 'LineWidth', 2);
                ylabel(labels_p{i}, 'FontWeight', 'bold');
                
                nexttile(tl_kin, (i-1)*2 + 2); hold on; grid on;
                plot(t_state, vel(i, :), 'Color', cols{i}, 'LineWidth', 2);
                ylabel(labels_v{i}, 'FontWeight', 'bold');
            end
            
            figHandles = [f1, f2, f3];
        end
        
        function f = plotInitialGuess(obj)
            %% PLOTINITIALGUESS  Visualize the 3-DoF forward guidance initial guess.
            if isempty(obj.InitialGuess)
                obj.generateInitialGuess();
            end
            
            t = obj.InitialGuess.Time;
            pos = obj.InitialGuess.X(5:7, :);
            vel = obj.InitialGuess.X(8:10, :);
            quat = obj.InitialGuess.X(1:4, :);
            thrust = obj.InitialGuess.U(3, :);
            t_ctrl = t(1:end-1);
            
            f = figure('Name', sprintf('Initial Guess: %s (%s)', obj.Maneuver, obj.Vehicle));
            tl = tiledlayout(f, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            % 3D Path
            nexttile(tl, 1); hold on; grid on; axis equal; view(3);
            plot3(pos(1,:), pos(2,:), pos(3,:), 'b-', 'LineWidth', 2);
            plot3(obj.r0(1), obj.r0(2), obj.r0(3), 'gs', 'MarkerFaceColor', 'g');
            plot3(obj.r_f(1), obj.r_f(2), obj.r_f(3), 'rs', 'MarkerFaceColor', 'r');
            title('3D Position Guess'); xlabel('X [m]'); ylabel('Y [m]'); zlabel('Z [m]');
            
            % Altitude and Velocity
            nexttile(tl, 2); hold on; grid on;
            plot(t, pos(3, :), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Z Pos');
            plot(t, vel(3, :), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Z Vel');
            title('Altitude and Vertical Speed'); xlabel('Time [s]'); legend('Location', 'best');
            
            % Attitude Quaternions
            nexttile(tl, 3); hold on; grid on;
            plot(t, quat(1,:), 'k', 'DisplayName', 'q_w');
            plot(t, quat(2,:), 'r', 'DisplayName', 'q_x');
            plot(t, quat(3,:), 'g', 'DisplayName', 'q_y');
            plot(t, quat(4,:), 'b', 'DisplayName', 'q_z');
            title('Attitude Quaternions Guess'); xlabel('Time [s]'); legend('Location', 'best');
            
            % Thrust Profile
            nexttile(tl, 4); hold on; grid on;
            plot(t_ctrl, thrust, 'Color', '#EDB120', 'LineWidth', 2);
            title('Thrust Command Guess'); xlabel('Time [s]'); ylabel('Thrust [N]');
        end
        
        function fn = getFormattedFilename(obj)
            %% GETFORMATTEDFILENAME  Standardized naming: Vehicle_ManeuverType_v###
            if strlength(obj.Filename) > 0
                fn = obj.Filename;
            else
                fn = sprintf('%s_%s_v%03d', obj.Vehicle, obj.Maneuver, obj.Version);
            end
        end
        
        function [tbl, filepath] = exportCSV(obj, filepath)
            %% EXPORTCSV  Export trajectory results table to CSV matching TOAD format.
            if isempty(obj.Solution) || ~isfield(obj.Solution, 'X')
                error('No solution available for export.');
            end
            
            if nargin < 2 || isempty(filepath)
                default_name = obj.getFormattedFilename() + ".csv";
                filepath = fullfile(pwd, 'Guidance', 'Trajectories', default_name);
            end
            
            t_state = obj.Solution.Time(:);
            X_res   = obj.Solution.X;
            U_res   = obj.Solution.U;
            
            quat     = X_res(1:4, :)';
            pos      = X_res(5:7, :)';
            vel      = X_res(8:10, :)';
            ang_rate = X_res(11:13, :)';
            m_lox    = X_res(14, :)';
            m_fuel   = X_res(15, :)';
            
            theta_out  = [U_res(1, :), U_res(1, end)]';
            phi_out    = [U_res(2, :), U_res(2, end)]';
            thrust_out = [U_res(3, :), U_res(3, end)]';
            roll_out   = [U_res(4, :), U_res(4, end)]';
            
            tbl = table(t_state, ...
                quat(:,1), quat(:,2), quat(:,3), quat(:,4), ...
                pos(:,1), pos(:,2), pos(:,3), ...
                vel(:,1), vel(:,2), vel(:,3), ...
                ang_rate(:,1), ang_rate(:,2), ang_rate(:,3), ...
                m_lox, m_fuel, ...
                theta_out, phi_out, thrust_out, roll_out, ...
                'VariableNames', {'Time', 'QuatW', 'QuatX', 'QuatY', 'QuatZ', ...
                                  'PosX', 'PosY', 'PosZ', ...
                                  'VelX', 'VelY', 'VelZ', ...
                                  'AngRateX', 'AngRateY', 'AngRateZ', ...
                                  'MassLox', 'MassFuel', ...
                                  'GimbalTheta', 'GimbalPhi', 'ThrustMag', 'RollCmd'});
            
            out_dir = fileparts(filepath);
            if ~exist(out_dir, 'dir') && ~isempty(out_dir)
                mkdir(out_dir);
            end
            writetable(tbl, filepath);
            fprintf('Trajectory successfully exported to: %s\n', filepath);
        end
    end
end
