%% This file contains a sequence of two parallel Linear Extended State Observers
% The main goal of LESO's is to estimate disturbances felt by the vehicle
% online and correct for them quickly. The sequence contains two parallel
% LESOs: A 3 channel coupled LESO for Attitude estimating the net
% torque-equivalent disturbance, and a single channel scalar LESO
% estimating the equivalent thrust disturbance.

function [D_Att, D_Thrust, U_corr] = LESO_Thrust(GND, X_est, X_trg, U_trg, L_Att, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat_att 
    % pos, vel, disturbances
    persistent xhat_thr  

    if isempty(t_last)
        t_last = t;
        xhat_att = zeros(6,1);
        xhat_thr = zeros(2,1);
        D_Att = zeros(3,1);
        D_Thrust = 0;
        U_corr = zeros(4,1);
        return
    end

    if GND == 1
        t_last = t; 
        xhat_att(:) = 0;
        xhat_thr(:) = 0;

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

    dT_LESO = 1/1000;
    
    %% Second Order thrust LESO
    idx_thr = [4:6, 7:9];
    % body to inertial
    % C_BI_ref = quatRot(X_trg(1:4));

    % A_d_thr = A_d(idx_thr, idx_thr);   
    b_0 = B_d(idx_thr, 3);

    A_LESO_Thr = [1 dT_LESO dT_LESO^2/2;
                  0   1         dT_LESO;
                  0   0              1];
    B_LESO_Thr = [b_0 * dT_LESO^2/2;
                        b_0*dT_LESO;
                                 0];

    y_thr_abs = X_est(idx_thr);

    %xhat_thr_pred = A_LESO_Thr * xhat_thr + B_LESO_Thr * U_thr_abs;
    xhat_thr_pred = X_trg(idx_thr);
    xhat_thr = xhat_thr_pred + L_Thrust * (y_thr_abs - xhat_thr_pred)
    

    % Turn this disturbance error into a quaternion to modify attitude and
    % error
    
    U_corr = zeros(4,1);
end


function [thrust, quat]=quaterror()
%%% 
% extract the quaternion and thrust for rotation to counter the
% disturbances
end