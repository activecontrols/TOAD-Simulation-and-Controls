%% This file contains a sequence of two parallel Linear Extended State Observers
% The main goal of LESO's is to estimate disturbances felt by the vehicle
% online and correct for them quickly. A 3 channel scalar LESO
% estimating the equivalent thrust disturbance.

function [D_Att, D_Thrust, U_corr] = LESO_Position(GND, X_est, X_trg, U_trg, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat 

    if isempty(t_last)
        t_last = t;
        xhat = zeros(9,1);
        D_Att = zeros(3,1);
        D_Thrust = 0;
        U_corr = zeros(4,1);
        return
    end

    if GND == 1
        t_last = t; 
        xhat(:) = 0;

        D_Att = zeros(3,1);
        D_Thrust = 0;
        U_corr = zeros(4,1);
        return
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        D_Thrust = xhat(3);
        D_Att = zeros(4,1);
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
    
    %% Second Order thrust LESO
    idx_thr = [5:7, 8:10];
    
    % Body to inertial
    C_BI_ref = quatRot(X_est(1:4));
    b_0 = C_BI_ref / (constantsTOAD.m_dry + sum(X_est(14:15)));

    % Linearized matrices
    A_LESO_Thr = [eye(3,3) dT*eye(3,3) dT^2/2*eye(3,3);
                  zeros(3)   eye(3,3)         dT*eye(3,3);
                  zeros(3)   zeros(3)              eye(3,3)];
    B_LESO_Thr = [b_0 * dT.^2/2;
                         b_0*dT;
                         zeros(3,3)];

    % Prediction
    y_abs = X_est(idx_thr);
    theta = U_trg(1); 
    phi = U_trg(2);
    ThrustVec = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
    U = U_trg(3) * ThrustVec;
    g_vec = [0; 0; -constantsTOAD.g];
    xhat_pred = A_LESO_Thr * xhat + B_LESO_Thr * U + ...
                [g_vec * dT^2/2; g_vec * dT; zeros(3,1)];

    % Update
    xhat = xhat_pred + L_Thrust * (y_abs - xhat_pred(1:6));

    e = [0; 0; 1];
    e_thr = C_BI_ref * e;
    D_Thrust = e_thr' * xhat(7:9);
    D_Att = zeros(4,1);

    %% Correction
    U_corr_th = -(C_BI_ref' * xhat(7:9))' * e * constantsTOAD.m_wet;

    U_corr = zeros(4,1);
    U_corr(3) = min(max(U_corr_th, -500), 500);
end


function [thrust, quat] = quaterror(thrust_corr, U_trg)
%%% 
% extract the quaternion and thrust for rotation to counter the
% disturbances
% Input:
%   Thrust: x y z components

[cos(theta) * sin(gamma), sin(theta) * cos(gamma), sin(theta) * sin(gamma)];


end