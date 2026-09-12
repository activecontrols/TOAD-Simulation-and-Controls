function [U_cmd, SpecRad, X_err] = TOAD_TVLQI_2(GND, X_est, X_trg, U_ff, K, t, constantsTOAD, Dist_LESO)
% LESO Augmented Time Varying LQR formulation
%
%   Augmented with Dynamic Cone Clamping & Altitude Pullout Protection to
%   avoid crashes during highly dynamic manouvers, as well as integral
%   action augmentation
    persistent t_last U_last LastSpecRad I_trans
    
    % Reset persistent states when on ground or first initialization
    if isempty(t_last) || GND == 1
        t_last = t;
        U_last = [0; 0; constantsTOAD.m_wet * constantsTOAD.g; 0];
        LastSpecRad = 0.99 * ones(2,1);
        I_trans = zeros(3, 1);
        if GND == 1
            U_cmd = U_last;
            SpecRad = 0.99 * ones(2,1);
            X_err = zeros(12,1);
            debug_info = struct('crab_angle', [0;0], 'f_req', [0;0;0], 'U_dist', zeros(4,1));
            return;
        end
    end
    dT = t - t_last;
    if dT <= 0, dT = 1/500; end
    t_last = t;

    % Safeguard for missing LESO connection
    if nargin < 8 || isempty(Dist_LESO)
        Dist_LESO = zeros(6,1);
    end

    % Unpack gains and disturbance estimates
    K_trans  = K(1:3, 1:6);
    K_rot    = K(4:6, 1:6);
    a_dist   = Dist_LESO(1:3);
    ang_accel_dist = Dist_LESO(4:6);
    
    Mass = constantsTOAD.m_dry + sum(X_est(14:15));
    g_vec = [0; 0; -constantsTOAD.g];

    %% Translational Trim (Outer Loop w/ Integrator)
    pos_err = X_est(5:7) - X_trg(5:7);
    vel_err = X_est(8:10) - X_trg(8:10);
    
    % Leaked integration for the translational integral trim
    I_trans = I_trans * 0.999 + pos_err * dT;
    MaxInt = [0.35; 0.35; 0.15]; 
    I_trans = max(min(I_trans, MaxInt), -MaxInt);
    K_I = diag([0.45, 0.45, 0.15]); 
    Delta_A = -K_trans * [pos_err; vel_err] - K_I * I_trans;
    
    % Acceleration saturation limit
    MaxAccelCorr = [1.5; 1.5; 4.0];
    Delta_A = max(min(Delta_A, MaxAccelCorr), -MaxAccelCorr);

    % Nominal Trajectory Acceleration
    Q_ref = X_trg(1:4);
    Q_est = X_est(1:4);
    C_B2I_ref = quatRot(Q_ref)'; 
    C_B2I_est = quatRot(Q_est)'; 
    T_B_ff = U_ff(3) * [cos(U_ff(1))*sin(U_ff(2)); -sin(U_ff(1)); cos(U_ff(1))*cos(U_ff(2))];
    a_ff = (C_B2I_ref * T_B_ff) / Mass + g_vec;

    % Commanded inertial acceleration with LESO disturbance compensation
    a_cmd = a_ff - a_dist + Delta_A;
    f_req = a_cmd - g_vec; 

    %% Dynamic Pullout Protection & Tilt Clamping
    Z_B_trg = C_B2I_ref * (T_B_ff / max(norm(T_B_ff), 1e-4));
    z_curr  = X_est(7);
    vz_curr = X_est(10);
    
    %% Altitude & Sink-Rate Adaptive Cone Clamping
    % If current velocity is negative indicating descent, we calculate the
    % stopping distance given our current altitude target and limit the
    % tilt accordingly.
    if vz_curr < -1.0
        thrustDist = min(a_dist(3), 0);
        availAccel = max((constantsTOAD.MaxThrust / Mass + thrustDist) - constantsTOAD.g, 0.4);
        stopDist = 0.5 * vz_curr^2 / availAccel;
        
        if z_curr < stopDist * 1.35 + 2.5
            % Emergency Braking Pullout: Lock thrust strictly vertical to maximize deceleration
            MaxConeAng = deg2rad(2.0);
            f_req(3) = max(f_req(3), constantsTOAD.MaxThrust / Mass);
        elseif z_curr < stopDist * 1.80 + 6.0
            % Anticipatory Braking, squeeze cone angle smoothly
            coneRatio = (z_curr - stopDist) / max(stopDist, 1.0);
            MaxConeAng = deg2rad(15) * min(max(coneRatio, 0.15), 1.0);
            f_req(3) = max(f_req(3), constantsTOAD.g + 1.2);
        else
            MaxConeAng = deg2rad(15);
        end
    else
        MaxConeAng = deg2rad(15);
    end

    % Project the required force direction into the allowable cone relative
    % to nominal thrust axis
    norm_f = norm(f_req);
    dotVal = dot(f_req, Z_B_trg) / (norm_f * norm(Z_B_trg));
    if dotVal < cos(MaxConeAng)
        f_perp = f_req - dot(f_req, Z_B_trg) * Z_B_trg;
        if norm(f_perp) > 1e-6
            f_perp = f_perp / norm(f_perp);
            f_req = (cos(MaxConeAng) * Z_B_trg + sin(MaxConeAng) * f_perp) * norm_f;
        else
            f_req = Z_B_trg * norm_f;
        end
        norm_f = norm(f_req);
    end

    if norm_f > 1e-6
        Z_t = f_req / norm_f;
    else
        Z_t = [0; 0; 1];
    end

    %% Actuator Disturbance Trim & Intermediate Frame C_T2B
    [J_tot, lever_arm] = ComputeJtot(X_est(14), X_est(15), constantsTOAD);
    torque = ang_accel_dist' * J_tot;
    U_cmd = zeros(4,1);

    % Correct input clamping to actual limits
    U_cmd(3) = max(min(norm_f * Mass, constantsTOAD.MaxThrust), 0.2 * constantsTOAD.MaxThrust);

    % Gimbal angles required to cancel CoM offset torque:
    U_dist = zeros(4,1);
    U_dist(3) = torque(3);
    U_dist(1) = -1 * torque(1) / (lever_arm * U_cmd(3));
    U_dist(2) = -1 * torque(2) / (lever_arm * (1 - U_ff(1)^2/2) * U_cmd(3));

    % Construct Thrust-to-Body Rotation Matrix (C_T2B)
    theta_eff = U_ff(1) - U_dist(1); 
    phi_eff   = U_ff(2) - U_dist(2);
   
    R_x = [1, 0, 0; 
           0, cos(theta_eff), -sin(theta_eff); 
           0, sin(theta_eff), cos(theta_eff)];
           
    R_y = [cos(phi_eff), 0, sin(phi_eff); 
           0, 1, 0; 
           -sin(phi_eff), 0, cos(phi_eff)];
           
    C_T2B = R_y * R_x; 

    %% Triad Generation 
    % Accounts for CoM shifts.
    X_t_ref = C_B2I_ref * C_T2B(:, 1);
    X_t_ref = X_t_ref - dot(X_t_ref, Z_t) * Z_t;
    
    if norm(X_t_ref) > 1e-6
        X_t_ref = X_t_ref / norm(X_t_ref);
        Y_t = cross(Z_t, X_t_ref);
        Y_t = Y_t / norm(Y_t);
        X_t = cross(Y_t, Z_t);
    else
        Y_t = [0; 1; 0];
        X_t = cross(Y_t, Z_t);
        X_t = X_t / norm(X_t);
        Y_t = cross(Z_t, X_t);
    end

    C_T2I_cmd = [X_t, Y_t, Z_t];
    C_B2I_cmd = C_T2I_cmd * C_T2B'; 
    Q_cmd = DCM_Quat_Conversion(C_B2I_cmd);
    Q_cmd = Q_cmd / norm(Q_cmd);
    if Q_cmd(1) < 0
        Q_cmd = -Q_cmd;
    end

    %% Rotational Inner-Loop Control
    Q_cmd_Conj = [Q_cmd(1); -Q_cmd(2:4)];
    AttError = HamiltonianProd(Q_cmd_Conj) * X_est(1:4);
    if AttError(1) < 0
        AttError = -AttError;
    end
    
    omegaTRG = C_B2I_est' * C_B2I_ref * X_trg(11:13);

    n = norm(AttError(2:4));
    if n > 1e-6
        Axis = AttError(2:4) / n;
    else
        Axis = zeros(3,1);
    end
    ThetaErr = 2 * atan2(n, AttError(1)) * Axis;
    X_err_rot = [ThetaErr;
                 X_est(11:13) - omegaTRG];
    Delta_u = -K_rot * X_err_rot;

    %% Actuator Trim Integration & Disturbance Allocation
    Channels = [1, 2, 4];
    thrustMax = constantsTOAD.MaxThrust;
    gimbalMax = deg2rad(12);
    InputBounds = [-gimbalMax       gimbalMax;
                   -gimbalMax       gimbalMax;
                   0.2 * thrustMax  thrustMax;
                   -4.0             4.0];
    
    % Feedback Trim Bounds
    MaxTrim = [ones(2,1) * deg2rad(10); 3.0];
    Trim_FB = min(max(Delta_u(:), -MaxTrim), MaxTrim);
    
    % Active Gimbal Disturbance Trim Feedforward
    MaxDist   = [ones(2,1) * deg2rad(10); 3.0]; 
    U_dist_req = min(max(-U_dist(Channels), -MaxDist), MaxDist);
    
    % Combine feedback and disturbance trim
    Trim = Trim_FB + U_dist_req;
    
    % Allocate to actuator channels
    U_cmd(Channels, 1) = Trim + U_ff(Channels);
    U_cmd = min(max(U_cmd, InputBounds(:, 1)), InputBounds(:, 2));

    %% Spectral Radius Telemetry (at 10 Hz)
    if mod(round(t*100), 10) == 0
        A_est = JacobianX(X_est, U_cmd);
        B_est = JacobianU(X_est, U_cmd);
    
        T_sr = zeros(15, 12);
        T_sr(1:4, 1:3)   = 0.5 * XiMat(Q_est);
        T_sr(5:13, 4:12) = eye(9);
        A_lin_sr = pinv(T_sr) * A_est * T_sr;
        B_lin_sr = pinv(T_sr) * B_est;
    
        A_trans_c_sr = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
        B_trans_c_sr = [zeros(3,3); eye(3,3)];
        M_c_trans_sr = [A_trans_c_sr, B_trans_c_sr; zeros(3,6), zeros(3,3)];
        M_d_trans_sr = expm(M_c_trans_sr * dT);
        A_trans_d_sr = M_d_trans_sr(1:6, 1:6);
        B_trans_d_sr = M_d_trans_sr(1:6, 7:9);
    
        A_rot_c_sr = [A_lin_sr(1:3,1:3),   A_lin_sr(1:3,10:12);
                      A_lin_sr(10:12,1:3), A_lin_sr(10:12,10:12)];
        B_rot_c_sr = [B_lin_sr(1:3,  [1,2,4]);
                      B_lin_sr(10:12,[1,2,4])];
        M_c_rot_sr = [A_rot_c_sr, B_rot_c_sr; zeros(3,6), zeros(3,3)];
        M_d_rot_sr = expm(M_c_rot_sr * dT);
        A_rot_d_sr = M_d_rot_sr(1:6, 1:6);
        B_rot_d_sr = M_d_rot_sr(1:6, 7:9);
    
        eig_cl_trans = eig(A_trans_d_sr - B_trans_d_sr * K_trans);
        eig_cl_rot   = eig(A_rot_d_sr   - B_rot_d_sr   * K_rot);
        SpecRad = [max(abs(eig_cl_trans)); max(abs(eig_cl_rot))];
    else
        SpecRad = LastSpecRad;
    end
    LastSpecRad = SpecRad;

    % Diagnostics and logging
    X_err = [X_err_rot; [pos_err; vel_err]];
    U_last = U_cmd;
    
    if nargout > 3
        debug_info = struct();
        debug_info.crab_angle = [phi_eff; theta_eff];
        debug_info.f_req = f_req;
        debug_info.U_dist = U_dist;
        debug_info.Trim_Dist = U_dist_req;
    end
end
