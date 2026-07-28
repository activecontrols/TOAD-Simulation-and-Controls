%% This file contains a sequence of two parallel Linear Extended State Observers
% The main goal of LESO's is to estimate disturbances felt by the vehicle
% online and correct for them quickly. The sequence contains two parallel
% LESOs: A 3 channel coupled LESO for Attitude estimating the net
% torque-equivalent disturbance, and a single channel scalar LESO
% estimating the equivalent thrust disturbance.

function [D_Att, D_Thrust, U_corr] = LESO_Sequence(X_est, X_trg, U_trg, L_Att, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat_att   
    persistent xhat_thr  

    if isempty(t_last)
        t_last = t;
        xhat_att = zeros(9,1);
        xhat_thr = zeros(2,1);
        D_Att = zeros(3,1);
        D_Thrust = 0;
        U_corr = zeros(4,1);
        return
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        D_Att = xhat_att(7:9);
        D_Thrust = xhat_thr(2);
        U_corr = zeros(4,1);
        return
    end

    %% Jacobian Eval
    X_lin = X_trg;
    X_lin(14:15) = X_est(14:15);

    Jx = JacobianX(X_lin, U_trg);
    Ju = JacobianU(X_lin, U_trg);

    % Kinematic mapping, evaluated at the reference quaternion
    T = zeros(15, 12);
    q_ref = X_trg(1:4);
    T(1:4, 1:3) = 0.5 * XiMat(q_ref);
    T(5:13, 4:12) = eye(9);

    A_lin = pinv(T) * Jx * T;
    B_lin = pinv(T) * Ju;

    nx = size(A_lin, 1);
    A_d = eye(nx) + A_lin * dT;
    B_d = B_lin * dT;

    % Current inertia, evaluated at the *estimated* propellant state
    J_tot = ComputeJtot(X_est(14), X_est(15), constantsTOAD);

    %% Attitude LESO
        idx_att = [1:3, 10:12];
        A_d_Att = A_d(idx_att, idx_att);   
        B_d_Att = B_d(idx_att, :);      
    
        % Disturbance torque coupling
        Bd_dist_att = [zeros(3,3); J_tot \ eye(3)] * dT; 
    
        A_LESO_Att = [A_d_Att, Bd_dist_att;
                      zeros(3,6), eye(3)  ];   % 9x9
        B_LESO_Att = [B_d_Att; zeros(3,4)];    % 9x4
        C_LESO_Att = [eye(6), zeros(6,3)];     % 6x9
    
        % Measurement attitude/rate tracking error
        dtheta_meas = QuatError(X_est(1:4), X_trg(1:4));  
        domega_meas = X_est(11:13) - X_trg(11:13);       
        y_att = [dtheta_meas; domega_meas];
    
        % Predict + correct
        U_att = U_trg;  
        xhat_att_pred = A_LESO_Att * xhat_att + B_LESO_Att * U_att;
        y_hat_att = C_LESO_Att * xhat_att_pred;
        xhat_att = xhat_att_pred + L_Att * (y_att - y_hat_att);
    
        D_Att = xhat_att(7:9);

    %% Thrust LESO
        idx_th = 7:9; 
        C_BI_ref = quatRot(X_trg(1:4));
        C_IB_ref = C_BI_ref.';
        e_thrust = C_IB_ref * [0;0;1];
    
        A_d_Thr = e_thrust' * A_d(idx_th, idx_th) * e_thrust;   
        B_d_Thr = e_thrust' * B_d(idx_th, 3);                    
    
        A_LESO_Thr = [A_d_Thr, dT;
                      0,       1 ];   
        B_LESO_Thr = [B_d_Thr; 0];     
        C_LESO_Thr = [1, 0];           
    
        % Measurement velocity along the thrust axis
        y_thr = e_thrust' * X_est(8:10);
    
        U_thr = U_trg(3);   
        xhat_thr_pred = A_LESO_Thr * xhat_thr + B_LESO_Thr * U_thr;
        y_hat_thr = C_LESO_Thr * xhat_thr_pred;
        xhat_thr = xhat_thr_pred + L_Thrust * (y_thr - y_hat_thr);
    
        D_Thrust = xhat_thr(2);

    %% Corrections 
    B_att_torque = B_lin(10:12, [1, 2, 4]);   
    U_corr_att = -pinv(B_att_torque) * (J_tot \ D_Att);

    b0 = e_thrust' * B_lin(idx_th, 3);
    b0_min = 1e-6;
    if abs(b0) > b0_min
        U_corr_th = -D_Thrust / b0;
    else
        U_corr_th = 0;
    end
    
    U_corr = zeros(4,1);
    U_corr([1, 2, 4]) = U_corr_att;
    U_corr(3) = U_corr_th;

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