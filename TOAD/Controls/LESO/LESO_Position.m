%% This file contains a sequence of two parallel Linear Extended State Observers
% The main goal of LESO's is to estimate disturbances felt by the vehicle
% online and correct for them quickly. A 3 channel scalar LESO
% estimating the equivalent thrust disturbance.

function [Att_Corr, U_Corr] = LESO_Position(GND, X_est, X_trg, U_trg, L_Thrust, constantsTOAD, t)

    persistent t_last
    persistent xhat 
    persistent AltErrInt

    if isempty(t_last)
        t_last = t;
        xhat = zeros(9,1);
        AltErrInt = 0;
        Att_Corr = [1; 0; 0; 0];
        U_Corr = 0;
        return
    end

    if GND == 1
        t_last = t; 
        xhat(:) = 0;
        AltErrInt = 0;
        Att_Corr = [1; 0; 0; 0];
        U_Corr = 0;
        return
    end

    dT = t - t_last;
    t_last = t;
    if dT <= 0
        Att_Corr = zeros(4,1);
        U_Corr = 0;
        return
    end
    
    %% Second Order thrust LESO
    idx_thr = [5:7, 8:10];
    
    % Body to inertial
    C_BI_ref = quatRot(X_est(1:4));
    Mass = (constantsTOAD.m_dry + sum(X_est(14:15)));
    b_0 = C_BI_ref / Mass;

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

    % Disturbance
    e = [0; 0; 1];
    D_Thrust = (C_BI_ref' * xhat(7:9))' * e;

    %% Corrections (with added integral trim for thrust)
    Ki = 0.10;
    Clamp = 500;
    AltErr = X_trg(7) - X_est(7);
    U_corr_th = -D_Thrust * Mass + Ki * AltErrInt;
    
    % Anti-windup
    if abs(U_corr_th) < Clamp || sign(U_corr_th) ~= sign(AltErr)
        AltErrInt = AltErrInt + AltErr * dT;
    end

    U_Corr = min(max(U_corr_th, -Clamp), Clamp);
    Thrust = U_Corr + U_trg(3);
    Att_Corr = quatError(xhat(7:9), constantsTOAD, C_BI_ref, X_est, Thrust);
end


function DelQ = quatError(distVec, constantsTOAD, C_B2I, X_est, Thrust)
    %%% 
    % Extract the quaternion and thrust for rotation to counter the
    % disturbances
    
    % Disturbance transformations
    bodyZ = [0;0;1];
    DistBody = C_B2I' * distVec;
    DistLat = DistBody(1:2);
    eps_lat = 1e-4;
    if norm(DistLat) < eps_lat
        Tilt = 0;
        Axis = [0;0;0];
    else
        % Required tilt angle 
        MaxTilt = pi/24;
        Mass = constantsTOAD.m_dry + sum(X_est(14:15));
        Tilt = Mass * norm(DistLat) / Thrust;
        Tilt = max(min(Tilt, MaxTilt), -MaxTilt);
        
        % Correction
        Axis = cross(bodyZ, [-DistLat / norm(DistLat); 0]);
    end
    
    DelQ = [cos(Tilt/2); sin(Tilt/2) * Axis];

end