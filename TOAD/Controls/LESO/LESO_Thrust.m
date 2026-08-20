%% This file contains a sequence of two parallel Linear Extended State Observers
% The main goal of LESO's is to estimate disturbances felt by the vehicle
% online and correct for them quickly. A 3 channel scalar LESO
% estimating the equivalent thrust disturbance.

function [D_Att, D_Thrust, U_corr] = LESO_Thrust(GND, X_est, X_trg, U_trg, L_Att, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat_att 
    % pos, vel, disturbances
    persistent xhat_thr  

    if isempty(t_last)
        t_last = t;
        xhat_att = zeros(6,1);
        % will become 9x1
        xhat_thr = zeros(9,1);
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
        D_Thrust = xhat_thr(3);
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

    pinv(T)
    A_lin = pinv(T) * Jx * T;
    B_lin = pinv(T) * Ju;

    nx = size(A_lin, 1);
    A_d = eye(nx) + A_lin * dT;
    B_d = B_lin * dT;

    dT_LESO = 1/1000;
    
    %% Second Order thrust LESO
    idx_thr = [4:6, 7:9];
    
     % body to inertial
    C_BI_ref = quatRot(X_est(1:4));
    C_IB_ref = C_BI_ref';

    % A_d_thr = A_d(idx_thr, idx_thr);   
    % b0_thr = B_d(idx_thr(1:3), 1:4);

    b_0 = C_IB_ref/constantsTOAD.m_wet;

    A_LESO_Thr = [eye(3,3) dT_LESO*eye(3,3) dT_LESO^2/2*eye(3,3);
                  zeros(3)   eye(3,3)         dT_LESO*eye(3,3);
                  zeros(3)   zeros(3)              eye(3,3)];
    B_LESO_Thr = [b_0 * dT_LESO.^2/2;
                         b_0*dT_LESO;
                         zeros(3,3)];

    y_thr_abs = X_est(idx_thr);
    theta = U_trg(1); 
    phi = U_trg(2);
    U = U_trg(3) * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

    xhat_thr_pred = A_LESO_Thr * xhat_thr + B_LESO_Thr * U;
    
    xhat_thr = xhat_thr_pred + L_Thrust * (y_thr_abs - xhat_thr_pred(1:6));


    e = [0 0 1];
   
    e_thr = e*C_IB_ref;
    D_Thrust = e_thr*xhat_thr(7:9);
    D_Att = 0;
    

    % Turn this disturbance error into a quaternion to modify attitude and
    % error
    
    U_corr = zeros(4,1);
    b0_min = 1e-6;

    
    abs(e_thr*b_0)
    
    if abs(e_thr*b_0) > b0_min
        U_corr_th = (D_Thrust*constantsTOAD.m_wet);
    else
        U_corr_th = 0;
    end
    
    U_corr_th
    U_corr(3) = U_corr_th;
    U_corr;
end


function [thrust, quat] = quaterror(thrust_corr, U_trg)
%%% 
% extract the quaternion and thrust for rotation to counter the
% disturbances
% Input:
%   Thrust: x y z components

[cos(theta) * sin(gamma), sin(theta) * cos(gamma), sin(theta) * sin(gamma)]


end