%% Translational Linear Extended State Observer (LESO)
% Estimates vehicle position, velocity, and unmodeled acceleration disturbances
% in the inertial frame (m/s^2).

function a_dist = LESO_Position(GND, X_est, X_trg, U_trg, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat 

    if isempty(t_last) || GND == 1
        t_last = t;
        xhat = zeros(9, 1);
        a_dist = zeros(3, 1);
        return;
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        a_dist = xhat(7:9);
        return;
    end
    
    %% Setup Dynamics & State Propagation
    idx_pos_vel = [5:7, 8:10];
    y_meas = X_est(idx_pos_vel);
    y_meas = y_meas(:);
    
    Mass = constantsTOAD.m_dry + sum(X_est(14:15));
    C_BI = quatRot(X_est(1:4));
    C_IB = C_BI';
    
    % Nominal thrust force vector in body frame
    theta = U_trg(1); 
    phi   = U_trg(2);
    thrust = U_trg(3);
    ThrustVec_B = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
    
    % Nominal acceleration in inertial frame (inc. gravity correction)
    g_vec = [0; 0; -constantsTOAD.g];
    a_nominal_I = (C_IB * ThrustVec_B) / Mass + g_vec;

    % Discrete STMs
    A_LESO_Thr = [eye(3), eye(3)*dT, eye(3)*dT^2/2;
                  zeros(3), eye(3),    eye(3)*dT;
                  zeros(3), zeros(3),  eye(3)];
                  
    B_LESO_Thr = [eye(3)*dT^2/2;
                  eye(3)*dT;
                  zeros(3,3)];

    % State Prediction 
    xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * a_nominal_I;

    % State Update 
    y_hat = xhat_pred(1:6);
    y_hat = y_hat(:);
    xhat = xhat_pred + L_Thrust * (y_meas - y_hat);

    % Disturbance Acceleration
    a_dist = xhat(7:9);
end