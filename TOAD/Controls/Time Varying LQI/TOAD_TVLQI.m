%% First generation Time Varying Linear Quadratic Integral controller for TOAD.
% This script receives an estimated state, a state target, feedforward
% input and a feedback gain matrix, and generates a corresponding control
% law using a TVLQI formualation. 

function [U_trg, X_err] = TOAD_TVLQI(X_est, X_trg, U_ff, K_f, t, constantsTOAD)

    persistent PosErrorI t_last
    if isempty(PosErrorI)
        PosErrorI = zeros(3,1);
        t_last = 0;
    end
    dT = t - t_last;
    t_last = t;

    %% Build Mutiplicative Quaternion Error
    Q_Conj = [X_est(1); -X_est(2:4, :)];
    Q_Trg  = X_trg(1:4);
    AttError = HamiltonianProd(Q_Conj) * Q_Trg;
    if AttError(1) < 0
        AttError = -AttError;
    end
    
    % Build the state error vector
    X_err = [2 * AttError(2:4); 
             X_trg(5:13) - X_est(5:13)];

    % Final control law
    U_trg = U_ff + K_f * X_err;

    % Clamping 
    thrustMax = constantsTOAD.MaxThrust;   %N
    gimbalMax = pi/12;
    InputBounds = [-gimbalMax       gimbalMax;
                   -gimbalMax       gimbalMax;
                   .2 * thrustMax   thrustMax;
                   -7               7];

    uMax = InputBounds(:, 2);
    uMin = InputBounds(:, 1);
    U_trg = min(max(U_trg, uMin), uMax);
end