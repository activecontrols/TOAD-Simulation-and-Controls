%% First generation Time Varying Linear Quadratic Integral controller for TOAD.
% This script receives an estimated state, a state target, feedforward
% input and a feedback gain matrix, and generates a corresponding control
% law using a TVLQI formualation. 

function [U_trg, U_fb, X_err] = TOAD_TVLQI(X_est, X_trg, U_ff, P_t, t, constantsTOAD, Corr)
    
    if nargin < 7 || isempty(Corr)
        Corr = [0; 1; 0; 0; 0];
    end

    persistent t_last U_last
    if isempty(t_last)
        t_last = 0;
        U_last = [0; 0; constantsTOAD.m_wet * constantsTOAD.g; 0];
    end
    dT = t - t_last;
    t_last = t;

    %% Build Mutiplicative Quaternion Error
    Q_Conj = [X_est(1); -X_est(2:4, :)];
    Q_Trg  = HamiltonianProd(X_trg(1:4)) * Corr(2:end);
    AttError = HamiltonianProd(Q_Conj) * Q_Trg;
    if AttError(1) < 0
        AttError = -AttError;
    end
    
    % Build the state error vector
    X_err = [2 * AttError(2:4); 
             X_trg(5:13) - X_est(5:13)];

    %% Gain Matrix Computation
    % Evaluate the Jacobians at the current estimated state and previous
    % input
    A_lin = JacobianX(X_est, U_last);
    B_lin = JacobianU(X_est, U_last);

    % Kinematic mapping (evaluated at target quaternion)
    T = zeros(15, 12);
    q_ref = X_trg(1:4);
    T(1:4, 1:3) = 0.5 * XiMat(q_ref);
    T(5:13, 4:12) = eye(9);

    % Reduce jacobians
    A_lin = pinv(T) * A_lin * T;
    B_lin = pinv(T) * B_lin;
    
    % Matrix discretization using ZOH
    nx = size(A_lin, 1); 

    %% Exact Matrix exponential computation, as loop >250Hz, not needed for now
    % nu = size(B_lin, 2); 
    % 
    % % Construct the continuous block matrix
    % M_c = [A_lin, B_lin; 
    %        zeros(nu, nx), zeros(nu, nu)];
    % 
    % % Discretize using the dT
    % M_d = expm(M_c * dT); 
    % A_d = M_d(1:nx, 1:nx);
    % B_d = M_d(1:nx, (nx+1):end);
    
    dT = mean(diff(constantsTOAD.Traj.Time));
    A_d = eye(nx) + A_lin * dT;
    B_d = B_lin * dT;
    R = constantsTOAD.R_Control;

    % Optimal gain matrix using the interpolated P_t passed into the function
    K_f = OptGain(A_d, B_d, R, P_t);

    % Final control law
    U_trg = U_ff + K_f * X_err;
    U_trg(3) = U_trg(3) + Corr(1);
    U_last = U_trg;
    U_fb = K_f * X_err;

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