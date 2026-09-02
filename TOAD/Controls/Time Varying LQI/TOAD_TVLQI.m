function [U_cmd, SpecRad, X_err] = TOAD_TVLQI(GND, X_est, X_trg, U_ff, K, t, constantsTOAD, Dist_LESO)
    % TOAD_TVLQI Trim Controller
    
    persistent t_last U_last LastSpecRad
    
    % Reset persistent variables if on the ground
    if isempty(t_last) || GND == 1
        t_last = t;
        U_last = [0; 0; constantsTOAD.m_wet * constantsTOAD.g; 0];
        LastSpecRad = 0.99 * ones(2,1);
        if GND == 1
            U_cmd = U_last;
            SpecRad = 0.99 * ones(2,1);
            X_err = zeros(12,1);
            return;
        end
    end
    dT = t - t_last;
    t_last = t;

    % Safeguard for missing LESO connection
    if nargin < 8 || isempty(Dist_LESO)
        Dist_LESO = zeros(6,1);
    end

    % Unpack 6x6 Gain and 6x1 Disturbance Vectors
    K_trans  = K(1:3, 1:6);
    K_rot    = K(4:6, 1:6);
    a_dist   = Dist_LESO(1:3);
    ang_accel_dist = Dist_LESO(4:6);
    
    Mass = constantsTOAD.m_dry + sum(X_est(14:15));
    g_vec = [0; 0; -constantsTOAD.g];

    %% Translational Trim (Outer Loop)
    X_err_trans = [X_est(5:7) - X_trg(5:7);
                   X_est(8:10) - X_trg(8:10)];

    Delta_A = -K_trans * X_err_trans;
    MaxAccelCorr = [1.5;1.5;4];
    Delta_A = max(min(Delta_A, MaxAccelCorr), -MaxAccelCorr);

    % Nominal NLP Acceleration
    Q_ref = X_trg(1:4);
    Q_est = X_est(1:4);
    C_B2I_ref = quatRot(Q_ref)'; 
    C_B2I_est = quatRot(Q_est)'; 
    T_B_ff = U_ff(3) * [cos(U_ff(1))*sin(U_ff(2)); -sin(U_ff(1)); cos(U_ff(1))*cos(U_ff(2))];
    a_ff = (C_B2I_ref * T_B_ff) / Mass + g_vec;

    % Torque from LESO estimate — already available, no dependency on a_cmd
    [J_tot,lever_arm] = ComputeJtot(X_est(14), X_est(15), constantsTOAD);
    torque = ang_accel_dist' * J_tot;
    a_cmd = a_ff - a_dist + Delta_A;

    %% Triad Generation & Roll Tracking
    % Section needs to be modified to align the THRUST vector with the
    % desired accel, not the Z axis. We now build an intermediate frame and
    % convert between. This is done so we can properly account for CoM
    % shifts.
    f_req = a_cmd - g_vec; 
    norm_f = norm(f_req);

    % Thrust for use in angular disturbances
    U_cmd = zeros(4,1);
    U_cmd(3) = (norm_f * Mass);

    %% Disturbance tracking
    U_dist = zeros(4,1);
    U_dist(3) = torque(3);
    U_dist(1) = -1 * torque(1)/(lever_arm*U_cmd(3));
    U_dist(2) = -1 * torque(2)/(lever_arm*(1-U_cmd(1)^2/2)*U_cmd(3));

    % Construct the Thrust-to-Body Rotation Matrix (C_T2B)
    theta_eff = U_ff(1) - U_dist(1); 
    phi_eff   = U_ff(2) - U_dist(2);
   
    R_x = [1, 0, 0; 
           0, cos(theta_eff), -sin(theta_eff); 
           0, sin(theta_eff), cos(theta_eff)];
           
    R_y = [cos(phi_eff), 0, sin(phi_eff); 
           0, 1, 0; 
           -sin(phi_eff), 0, cos(phi_eff)];
           
    C_T2B = R_y * R_x; 
    
    % Construct the Thrust-to-Inertial Frame (C_T2I_cmd)
    if norm_f > 1e-6
        Z_t = f_req / norm_f;
    else
        Z_t = [0; 0; 1];
    end

    % Reference thrust-frame X axis corresponding to target body X axis
    X_t_ref = C_B2I_ref * C_T2B(:,1);
    
    % Re-orthogonalize against requested thrust direction
    X_t_ref = X_t_ref - dot(X_t_ref, Z_t) * Z_t;
    
    if norm(X_t_ref) > 1e-6
        X_t_ref = X_t_ref / norm(X_t_ref);
        Y_t = cross(Z_t, X_t_ref);
        Y_t = Y_t / norm(Y_t);
        X_t = cross(Y_t, Z_t);
    else
        % Fallback
        Y_t = [0; 1; 0];
        X_t = cross(Y_t, Z_t);
        X_t = X_t / norm(X_t);
        Y_t = cross(Z_t, X_t);
    end
    
    C_T2I_cmd = [X_t, Y_t, Z_t];
    
    % Apply Inverse Mapping to determine Commanded Body Attitude
    C_B2I_cmd = C_T2I_cmd * C_T2B'; 
    Q_cmd = DCM_Quat_Conversion(C_B2I_cmd);
    Q_cmd = Q_cmd / norm(Q_cmd);

    %% Rotational Trim
    Q_cmd_Conj = [Q_cmd(1); -Q_cmd(2:4)];
    AttError = HamiltonianProd(Q_cmd_Conj) * X_est(1:4);
    
    if AttError(1) < 0
        AttError = -AttError;
    end
    
    % Frame correction for the angular rates. Need to be represented in the
    % actual body frame not the target frame. UNTESTED FIX for now
    omegaTRG = C_B2I_est' * C_B2I_ref * X_trg(11:13);

    %% Test implementation
    n = norm(AttError(2:4));
    if n > 1e-6
        Axis = AttError(2:4) / n;
    else
        Axis = zeros(3,1);
    end
    ThetaErr = 2 * atan2(norm(AttError(2:4)), AttError(1)) * Axis;
    X_err_rot = [ThetaErr;
                 X_est(11:13) - omegaTRG];
    Delta_u = -K_rot * X_err_rot;

    
    %% Trim Integration & clamping   
    U_cmd(1) = U_ff(1) + Delta_u(1) - U_dist(1); 
    U_cmd(2) = U_ff(2) + Delta_u(2) - U_dist(2);               
    U_cmd(4) = U_ff(4) + Delta_u(3) - U_dist(3); 

    % SpecRad = [Delta_u(1); Delta_u(2); U_cmd(3) - U_ff(3); Delta_u(3)];
    %% Spectral radius calculation for analysis (Taken from MatrixVerif)

    % Run at 10Hz for efficiency
    if mod(round(t*100), 10) == 0
        % Calculate full Jacobians
        A_est = JacobianX(X_est, U_cmd);
        B_est = JacobianU(X_est, U_cmd);
    
        % Map to reduced (attitude-error) state via quaternion kinematics
        T_sr = zeros(15, 12);
        T_sr(1:4, 1:3)   = 0.5 * XiMat(Q_est);
        T_sr(5:13, 4:12) = eye(9);
        A_lin_sr = pinv(T_sr) * A_est * T_sr;
        B_lin_sr = pinv(T_sr) * B_est;
    
        % Translational model
        A_trans_c_sr = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
        B_trans_c_sr = [zeros(3,3); eye(3,3)];
        M_c_trans_sr = [A_trans_c_sr, B_trans_c_sr; zeros(3,6), zeros(3,3)];
        M_d_trans_sr = expm(M_c_trans_sr * dT);
        A_trans_d_sr = M_d_trans_sr(1:6, 1:6);
        B_trans_d_sr = M_d_trans_sr(1:6, 7:9);
    
        % Rotational model
        A_rot_c_sr = [A_lin_sr(1:3,1:3),   A_lin_sr(1:3,10:12);
                      A_lin_sr(10:12,1:3), A_lin_sr(10:12,10:12)];
        B_rot_c_sr = [B_lin_sr(1:3,  [1,2,4]);
                      B_lin_sr(10:12,[1,2,4])];
        M_c_rot_sr = [A_rot_c_sr, B_rot_c_sr; zeros(3,6), zeros(3,3)];
        M_d_rot_sr = expm(M_c_rot_sr * dT);
        A_rot_d_sr = M_d_rot_sr(1:6, 1:6);
        B_rot_d_sr = M_d_rot_sr(1:6, 7:9);
    
        % Closed-loop spectral radius
        eig_cl_trans = eig(A_trans_d_sr - B_trans_d_sr * K_trans);
        eig_cl_rot   = eig(A_rot_d_sr   - B_rot_d_sr   * K_rot);
        SpecRad = [max(abs(eig_cl_trans)); max(abs(eig_cl_rot))];
    else
        SpecRad = LastSpecRad;
    end

    LastSpecRad = SpecRad;
    % Error computation for output & clamping 
    X_err = [X_err_rot; X_err_trans]; 

    thrustMax = constantsTOAD.MaxThrust;
    gimbalMax = pi/12;
    InputBounds = [-gimbalMax       gimbalMax;
                   -gimbalMax       gimbalMax;
                   .2 * thrustMax   thrustMax;
                   -7               7];
    U_cmd = min(max(U_cmd, InputBounds(:, 1)), InputBounds(:, 2));
    U_last = U_cmd;
end