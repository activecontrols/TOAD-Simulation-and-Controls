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
    y_meas = X_est(1:3);
    y_meas = y_meas(:);
    
    Mass = constantsTOAD.m_dry + sum(X_est(14:15));
    C_BI = quatRot(X_est(1:4));
    C_IB = C_BI';
    
    % Nominal thrust force vector in body frame
    theta = U_trg(1); 
    phi   = U_trg(2);
    thrust = U_trg(3);

    ThrustVec_B = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

        
    % Fluid fill heights
    OxFluidHeight = (X_est(14) / constantsTOAD.OxMass) * constantsTOAD.OxHeight * 0.9;
    FuFluidHeight = (X_est(15) / constantsTOAD.FuMass) * constantsTOAD.FuHeight * 0.9;

    % Fluid fills & location of CoM
    OxFluidLocation = constantsTOAD.Ox_Z + OxFluidHeight / 2;
    FuFluidLocation = constantsTOAD.Fu_Z + FuFluidHeight / 2;

    
    CGz = (constantsTOAD.m_dry * constantsTOAD.rTB + X_est(14) * OxFluidLocation + X_est(15) * FuFluidLocation) / Mass;
    
    % Distances to CG
    lever_arm = constantsTOAD.rTB - CGz;
    d_lox = OxFluidLocation - CGz;
    d_ipa = FuFluidLocation - CGz;

    J_xx = 1/12 * X_est(14) * (3 * constantsTOAD.OxRadius^2 + OxFluidHeight^2);
    J_zz = 1/2 * X_est(14) * constantsTOAD.OxRadius^2;
    J_lox = [J_xx, 0, 0;
             0, J_xx, 0;
             0,    0, J_zz];
    
    J_xx = 1/12 * X_est(15) * (3 * constantsTOAD.FuRadius^2 + FuFluidHeight^2);
    J_zz = 1/2 * X_est(15) * constantsTOAD.FuRadius^2;
    J_ipa = [J_xx, 0, 0;
             0, J_xx, 0;
             0,    0, J_zz];
    
    J_dry = constantsTOAD.J + constantsTOAD.m_dry * diag([lever_arm^2, lever_arm^2, 0]);
    J_lox = J_lox + X_est(14) * diag([d_lox^2, d_lox^2, 0]);
    J_ipa = J_ipa + X_est(15) * diag([d_ipa^2, d_ipa^2, 0]);

    J_tot = J_dry + J_lox + J_ipa;
    lever_arm_vec = [0;0;lever_arm];

    anglAccel = cross(lever_arm_vec, ThrustVec_B)' / J_tot;

    % Discrete STMs
    A_LESO_Thr = [eye(3), eye(3)*dT;
                  zeros(3), eye(3)];
                  
    B_LESO_Thr = [eye(3)*dT;
                  zeros(3)];

    % State Prediction 
    xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * anglAccel';

    % State Update 
    y_hat = xhat_pred(1:3);
    y_hat = y_hat(:);
    xhat = xhat_pred + L_Att * (y_meas - y_hat);

    % Disturbance Acceleration
    U_dist = xhat(4:6);
end

% Multiplicative quaternion error computation
function dtheta = QuatError(q_est, q_ref)
    Q_Ref_Conj = [q_ref(1); -q_ref(2:4)];
    q_err = HamiltonianProd(Q_Ref_Conj) * q_est;
    if q_err(1) < 0
        q_err = -q_err; 
    end
    dtheta = 2 * q_err(2:4);
end