%% Rotational Linear Extended State Observer (LESO)
% Estimates body angular rates and disturbance torques, outputting an
% equivalent actuator disturbance trim U_dist = [theta, phi, roll].

function U_dist = LESO_Attitude(GND, X_est, X_trg, U_trg, L_Att, constantsTOAD, t)
    persistent t_last
    persistent xhat 

    if isempty(t_last) || GND == 1
        t_last = t;
        xhat = zeros(6, 1);
        U_dist = zeros(3, 1);
        return;
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        U_dist = xhat(4:6);
        return;
    end
    
    %% Setup Dynamics & State Propagation
    y_meas = X_est(11:13);
    y_meas = y_meas(:);
    C_BI = quatRot(X_est(1:4));
    
    % Nominal thrust force vector in body frame
    theta = U_trg(1); 
    phi   = U_trg(2);
    thrust = U_trg(3);
    ThrustVec_B = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

    % Compute inertia
    [J_tot,lever_arm]  = ComputeJtot(X_est(14), X_est(15), constantsTOAD);
    lever_arm_vec = [0;0;-lever_arm];
    anglAccel = (cross(lever_arm_vec, ThrustVec_B) + [0; 0; U_trg(4)])' / J_tot;

    % Substract gyroscopic coupling
    gyroCoupling = J_tot^(-1) * cross(xhat(1:3), J_tot * xhat(1:3));

    % Discrete STMs
    A_LESO_Thr = [eye(3), eye(3)*dT;
                  zeros(3), eye(3)];
                  
    B_LESO_Thr = [eye(3)*dT;
                  zeros(3)];

    % State Prediction 
    xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * (anglAccel' - gyroCoupling);

    % State Update 
    y_hat = xhat_pred(1:3);
    y_hat = y_hat(:);
    xhat = xhat_pred + L_Att * (y_meas - y_hat);

    % Disturbance Acceleration
    U_dist = xhat(4:6);
end
