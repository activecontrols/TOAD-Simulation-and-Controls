%% Rotational Linear Extended State Observer (Enhanced LESO)
% Estimates body angular rates and disturbance torques, outputting an
% equivalent angular acceleration disturbance vector u_dist (rad/s^2)
% using actual commanded plant inputs (U_cmd).

function u_dist = LESO_Attitude_2(GND, X_est, X_trg, U_cmd, L_Att, constantsTOAD, t)
    persistent t_last
    persistent xhat 

    if isempty(t_last) || GND == 1
        t_last = t;
        xhat = zeros(6, 1);
        if ~isempty(X_est) && length(X_est) >= 13
            xhat(1:3) = X_est(11:13); % Initialize with measured body rates
        end
        u_dist = zeros(3, 1);
        return;
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        u_dist = xhat(4:6);
        return;
    end
    
    %% Setup Dynamics & State Propagation
    y_meas = X_est(11:13);
    y_meas = y_meas(:);
    
    % Commanded thrust force vector in body frame
    theta  = U_cmd(1); 
    phi    = U_cmd(2);
    thrust = U_cmd(3);
    ThrustVec_B = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

    % Compute inertia & lever arm
    [J_tot, lever_arm] = ComputeJtot(X_est(14), X_est(15), constantsTOAD);
    lever_arm_vec = [0; 0; -lever_arm];
    
    % Commanded torque divided by J_tot
    tau_cmd = cross(lever_arm_vec, ThrustVec_B) + [0; 0; U_cmd(4)];
    anglAccel = J_tot \ tau_cmd;

    % Gyroscopic coupling
    gyroCoupling = J_tot \ cross(y_meas, J_tot * y_meas);

    % Discrete State Transition Matrices
    A_LESO_Thr = [eye(3), eye(3)*dT;
                  zeros(3), eye(3)];
                  
    B_LESO_Thr = [eye(3)*dT;
                  zeros(3)];

    % State Prediction 
    xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * (anglAccel - gyroCoupling);

    % State Update with Observer Gain
    y_hat = xhat_pred(1:3);
    xhat = xhat_pred + L_Att * (y_meas - y_hat);

    % Disturbance Angular Acceleration (rad/s^2)
    u_dist = xhat(4:6);
end
