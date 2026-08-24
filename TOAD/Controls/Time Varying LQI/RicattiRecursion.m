function [K_trans_List, K_rot_List, L_List_Att, L_List_Thr] = RicattiRecursion(Trajectory, constantsTOAD)
% RICATTIRECURSION Find decoupled optimal gain matrices for outer/inner loops
%   Separates the LQR formulation into a 3x6 translational point-mass model
%   and a 3x6 angular dynamics model.

% LESO bandwidths
omega_att = 10;
omega_thr = 2;

N = size(Trajectory.x, 2);
K_trans_List = zeros(N, 3, 6);
K_rot_List   = zeros(N, 3, 6);
L_List_Att   = zeros(N, 6, 3);
L_List_Thr   = zeros(N, 9, 6);

% Load the cascaded tuning matrices
Q_trans = constantsTOAD.Q_trans;
R_trans = constantsTOAD.R_trans;
Q_rot   = constantsTOAD.Q_rot;
R_rot   = constantsTOAD.R_rot;

% Load symbolic functions (maybe reintroduce generation)
A_fcn = str2func('JacobianX');
B_fcn = str2func('JacobianU');

% Translational Model
A_trans_c = [zeros(3,3), eye(3,3); zeros(3,3), zeros(3,3)];
B_trans_c = [zeros(3,3); eye(3,3)];

for n = N:-1:2
    x_n   = Trajectory.x(:, n);
    u_n1  = Trajectory.u(:, n-1);
    dT    = 1/1000; %Trajectory.t(n) - Trajectory.t(n-1);
    
    A_lin = A_fcn(x_n, u_n1);
    B_lin = B_fcn(x_n, u_n1);

    % Kinematic mapping 
    T = zeros(15, 12);
    q_ref = x_n(1:4);
    T(1:4, 1:3) = 0.5 * XiMat(q_ref);
    T(5:13, 4:12) = eye(9);
    A_lin = pinv(T) * A_lin * T;
    B_lin = pinv(T) * B_lin;
    
    %% Translational Discretization
    M_c_trans = [A_trans_c, B_trans_c; zeros(3,6), zeros(3,3)];
    M_d_trans = expm(M_c_trans * dT);
    A_trans_d = M_d_trans(1:6, 1:6);
    B_trans_d = M_d_trans(1:6, 7:9);
    
    %% Rotational Discretization
    A_rot_c = [A_lin(1:3, 1:3), A_lin(1:3, 10:12);
               A_lin(10:12, 1:3), A_lin(10:12, 10:12)];
    B_rot_c = [B_lin(1:3, [1, 2, 4]);
               B_lin(10:12, [1, 2, 4])];
               
    M_c_rot = [A_rot_c, B_rot_c; zeros(3,6), zeros(3,3)];
    M_d_rot = expm(M_c_rot * dT);
    A_rot_d = M_d_rot(1:6, 1:6);
    B_rot_d = M_d_rot(1:6, 7:9);

    % Initialize converged terminal costs
    if n == N
        P_trans = idare(A_trans_d, B_trans_d, Q_trans, R_trans);
        P_rot   = idare(A_rot_d, B_rot_d, Q_rot, R_rot);
    end
    
    % Optimal Gains
    K_trans = (B_trans_d'*P_trans*B_trans_d + R_trans) \ (B_trans_d'*P_trans*A_trans_d);
    K_rot   = (B_rot_d'*P_rot*B_rot_d + R_rot) \ (B_rot_d'*P_rot*A_rot_d);
    
    K_trans_List(n, :, :) = K_trans;
    K_rot_List(n, :, :) = K_rot;
    
    % Riccati Updates
    P_trans = Q_trans + A_trans_d'*P_trans*A_trans_d - A_trans_d'*P_trans*B_trans_d * K_trans;
    P_rot   = Q_rot + A_rot_d'*P_rot*A_rot_d - A_rot_d'*P_rot*B_rot_d * K_rot;
    
    % LESO Scheduler
    dT_LESO = 1/1000;
    
    % Attitude LESO 
    A_LESO_A = [eye(3), eye(3)*dT_LESO; zeros(3,3), eye(3)];
    C_LESO_A = [eye(3), zeros(3,3)];
    
    % Thrust LESO 
    A_LESO_T = [eye(3), eye(3)*dT_LESO, eye(3)*dT_LESO^2/2;
                zeros(3), eye(3), eye(3)*dT_LESO;
                zeros(3), zeros(3), eye(3)];
    C_LESO_T = [eye(6), zeros(6,3)];
    
    % Pole Placement
    s_poles_att = -omega_att * (1 + (0:5)*0.01);
    z_poles_att = exp(s_poles_att * dT_LESO);
    L_att_n = place(A_LESO_A', C_LESO_A', z_poles_att)';
    
    s_poles_th = -omega_thr * (1 + (0:8)*0.01);
    z_poles_th = exp(s_poles_th * dT_LESO);
    L_th_n = place(A_LESO_T', C_LESO_T', z_poles_th)';
    
    L_List_Att(n, :, :) = L_att_n;
    L_List_Thr(n, :, :) = L_th_n;
end

% Copy 2nd to last gain matrices to first index
K_trans_List(1, :, :) = K_trans_List(2, :, :);
K_rot_List(1, :, :)   = K_rot_List(2, :, :);
L_List_Att(1, :, :)   = L_List_Att(2, :, :);
L_List_Thr(1, :, :)   = L_List_Thr(2, :, :);
end