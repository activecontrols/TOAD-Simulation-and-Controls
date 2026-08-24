function [U_cmd, U_fb, X_err] = TOAD_TVLQI(GND, X_est, X_trg, U_ff, K, t, constantsTOAD, Dist_LESO)
    % TOAD_TVLQI Trim Controller
    
    persistent t_last U_last Z_err_int
    
    % Reset persistent variables if on the ground
    if isempty(t_last) || GND == 1
        t_last = t;
        U_last = [0; 0; constantsTOAD.m_wet * constantsTOAD.g; 0];
        Z_err_int = 0;
        if GND == 1
            U_cmd = U_last;
            U_fb = zeros(4,1);
            X_err = zeros(12,1);
            return;
        end
    end
    dT = max(t - t_last, 0.001);
    t_last = t;

    % Safeguard for missing LESO connection
    if nargin < 8 || isempty(Dist_LESO)
        Dist_LESO = zeros(6,1);
    end

    % Unpack 6x6 Gain and 6x1 Disturbance Vectors
    K_trans  = K(1:3, 1:6);
    K_rot    = K(4:6, 1:6);
    a_dist   = Dist_LESO(1:3);
    U_dist = Dist_LESO(4:6);

    Mass = constantsTOAD.m_dry + sum(X_est(14:15));
    g_vec = [0; 0; -constantsTOAD.g];

    %% Translational Trim (Outer Loop)
    X_err_trans = [X_est(5:7) - X_trg(5:7);
                   X_est(8:10) - X_trg(8:10)];

    Delta_A = -K_trans * X_err_trans;
    MaxAccelCorr = 2;
    Delta_A = max(min(Delta_A, MaxAccelCorr), -MaxAccelCorr);

    % Nominal NLP Acceleration
    Q_ref = X_trg(1:4);
    C_B2I_ref = quatRot(Q_ref)'; 
    T_B_ff = U_ff(3) * [cos(U_ff(1))*sin(U_ff(2)); -sin(U_ff(1)); cos(U_ff(1))*cos(U_ff(2))];
    a_ff = (C_B2I_ref * T_B_ff) / Mass + g_vec; 

    % Thrust trim
    K_i = 40;
    Z_err_int = Z_err_int + (X_trg(7) - X_est(7)) * dT;
    ThrustTrim = K_i * Z_err_int;

    a_cmd = a_ff - a_dist + Delta_A; % + Delta_a - a_dist;

    %% Triad Generation & Roll Tracking
    f_req = a_cmd - g_vec; 
    norm_f = norm(f_req);
    
    if norm_f > 1e-6
        Z_b = f_req / norm_f;
    else
        Z_b = [0; 0; 1];
    end

    X_b_ref = C_B2I_ref(:, 1);
    Y_b_unnorm = cross(Z_b, X_b_ref);
    if norm(Y_b_unnorm) > 1e-6
        Y_b = Y_b_unnorm / norm(Y_b_unnorm);
    else
        Y_b = [0; 1; 0];
    end
    X_b = cross(Y_b, Z_b);

    C_IB_cmd = [X_b, Y_b, Z_b];
    
    Q_cmd = DCM_Quat_Conversion(C_IB_cmd);
    Q_cmd = Q_cmd / norm(Q_cmd);

    %% Rotational Trim
    Q_cmd_Conj = [Q_cmd(1); -Q_cmd(2:4)];
    AttError = HamiltonianProd(Q_cmd_Conj) * X_est(1:4);
    
    if AttError(1) < 0
        AttError = -AttError;
    end

    X_err_rot = [2 * AttError(2:4);
                 X_est(11:13) - X_trg(11:13)];

    Delta_u = -K_rot * X_err_rot;

    %% Trim Integration & clamping
    U_cmd = zeros(4,1);
    U_cmd(1) = U_ff(1) + Delta_u(1) - U_dist(1); 
    U_cmd(2) = U_ff(2) + Delta_u(2) - U_dist(2); 
    U_cmd(3) = (norm_f * Mass) + ThrustTrim;                 
    U_cmd(4) = U_ff(4) + Delta_u(3) - U_dist(3); 

    U_fb = [Delta_u(1); Delta_u(2); U_cmd(3) - U_ff(3); Delta_u(3)];
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