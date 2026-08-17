function [MatrixList, CostList, L_List_Att, L_List_Thr] = RicattiRecursion(Trajectory,Q, R, constantsTOAD)
% RICATTIRECURSION Find optimal gain matricies for sequence of waypoints
%   Backwards pass for finding optimal gain matrix at each sequence
%       Given trajectory (K, u, x, )
% evolves backwards in time from P_T = Q
% Q and R hand tuned for now, so they're inputs

% state variables, x
syms q0 q1 q2 q3
syms r1 r2 r3
syms v1 v2 v3
syms omega1 omega2 omega3
syms m_lox m_ipa
% control for u
syms theta phi thrust roll

% LESO bandwidths
omega_att = 1;
omega_thr = 1*pi/2;

MatrixList = permute(repmat([zeros(12,4)],[1,1,size(Trajectory.x,2)]), [3,1,2]);
CostList = permute(repmat([zeros(size(Q))],[1,1,size(Trajectory.x,2)]), [3,1,2]);
L_List_Att = permute(repmat(zeros(6,3), [1,1,size(Trajectory.x,2)]), [3,1,2]);
L_List_Thr = permute(repmat(zeros(9,6), [1,1,size(Trajectory.x,2)]), [3,1,2]);

q = [q0;q1;q2;q3];
r = [r1;r2;r3];
v = [v1;v2;v3];
m = [m_lox; m_ipa];
omegaB = [omega1; omega2; omega3];

xn = [q;r;v; omegaB; m]; % x_n

% U
u = [theta;phi;thrust;roll];

%% Dynamics
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
m_dry = constantsTOAD.m_dry;
m = m_dry + m_lox + m_ipa;

% Compute inertia
[J_tot, CGz] = ComputeJtot(m_lox, m_ipa, constantsTOAD);

% Angular dynamics
MB = zetaCross([0; 0; -CGz])*TB + (TB * roll)/thrust;
M = [q(1) -q(2) -q(3) -q(4);
     q(2)  q(1) -q(4)  q(3);
     q(3)  q(4)  q(1) -q(2);
     q(4) -q(3)  q(2)  q(1)];
qdot = 0.5 * M * [0; omegaB];
omegaBdot = J_tot \ (MB - zetaCross(omegaB) * J_tot * omegaB);

% Linear Dynamics
C_BI = quatRot(q);
C_IB = C_BI.';
FI = C_IB * TB + [0; 0;-m*constantsTOAD.g];

rdot = [v1;v2;v3];
vdot = FI/m;

% Mass Dynamics
mdot_lox = -thrust / constantsTOAD.MaxThrust * constantsTOAD.OF / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);
mdot_ipa = -thrust / constantsTOAD.MaxThrust * 1 / (1 + constantsTOAD.OF) * (constantsTOAD.MaxMdot);

%% State vector derivative
xdot = [qdot;rdot;vdot; omegaBdot;mdot_lox;mdot_ipa];

% Turn the jacobians into matlab functions for evaluation speed
A = jacobian(xdot, xn);
B = jacobian(xdot, u);
A_fcn = matlabFunction(A, 'Vars', {xn, u},'File', './Controls/Time Varying LQI/JacobianX');
B_fcn = matlabFunction(B, 'Vars', {xn, u},'File', './Controls/Time Varying LQI/JacobianU');

for n = size(Trajectory.x,2):-1:2
    % Extract states
    x_n   = Trajectory.x(:, n);
    u_n1  = Trajectory.u(:, n-1);
    dT    = Trajectory.t(n) - Trajectory.t(n-1);
    
    % Jacobian Eval
    A_lin = A_fcn(x_n, u_n1);
    B_lin = B_fcn(x_n, u_n1);

    % Kinematic mapping 
    T = zeros(15, 12);
    q_ref = x_n(1:4);
    T(1:4, 1:3) = 0.5 * XiMat(q_ref);
    T(5:13, 4:12) = eye(9);

    % Reduce jacobians
    A_lin = pinv(T) * A_lin * T;
    B_lin = pinv(T) * B_lin;
    
    %% Matrix discretization using ZOH
    % Dimensions
    nx = size(A_lin, 1); 
    nu = size(B_lin, 2); 
    
    % Construct the continuous block matrix
    M_c = [A_lin, B_lin; 
           zeros(nu, nx), zeros(nu, nu)];
    
    % Discretize and extract A_d and B_d
    M_d = expm(M_c * dT);
    A_d = M_d(1:nx, 1:nx);
    B_d = M_d(1:nx, (nx+1):end);

    % Initialize terminal cost as the result of the converged DARE for the
    % terminal state. 
    if n == size(Trajectory.x,2)
        P_t = idare(A_d, B_d, Q, R);
    end

    % LESO scheduler
    dT_LESO = 1/1000;

    %% Attitude LESO - Absolute Rates + Disturbance Acceleration
    % A_LESO_A assumes simple integration from acceleration to rate
    A_LESO_A = [eye(3),      eye(3) * dT_LESO;
                zeros(3,3),  eye(3)     ];         
    C_LESO_A = [eye(3), zeros(3,3)];                  
 
    %% Thrust LESO - Absolute Velocity + Disturbance Velocity Rate
    % A_LESO_T assumes simple integration from acceleration to velocity
    A_LESO_T = [eye(3) eye(3)*dT_LESO eye(3)*dT_LESO^2/2;
                  zeros(3,3)  eye(3), eye(3)*dT_LESO;
                  zeros(3,3)  zeros(3,3)  eye(3)];  
    C_LESO_T = [eye(6) zeros(6,3)];              
 
    %% Pole Placement
    s_poles_att = -omega_att * (1 + (0:5)*0.01);
    z_poles_att = exp(s_poles_att * dT_LESO);
    L_att_n = place(A_LESO_A', C_LESO_A', z_poles_att)';
 
    s_poles_th = -omega_thr * (1 + (0:8)*0.01);
    z_poles_th = exp(s_poles_th * dT_LESO);
    L_th_n = place(A_LESO_T', C_LESO_T', z_poles_th)';
 
    L_List_Att(n, :, :) = L_att_n;
    L_List_Thr(n, :, :) = L_th_n;

    % Matrix Eval and Updating Ricatti Cost
    MatrixList(n, :, :) = gain(A_d,B_d,R,P_t)';
    CostList(n, :, :)  = P_t;
    P_t = riccati(A_d,B_d, R , Q , P_t);
end

% Copy over 2nd to last gain matrix to first spot
MatrixList(1, :, :) = MatrixList(2, :, :);
CostList(1, :, :) = riccati(A_d,B_d, R , Q , P_t);
L_List_Att(1, :, :) = L_List_Att(2, :, :);
L_List_Thr(1, :, :) = L_List_Thr(2, :, :);

end

% Calculates P_(t-1)
function P = riccati(A, B, R, Q, P_t)
    P = Q + A'*P_t*A - (A'*P_t*B) / (B'*P_t*B + R) * B'*P_t*A;
end

% Optimal Gain Matrix
function K = gain(A, B, R, P_t)
    K = (B.'*P_t*B+R) \ B.'*P_t*A;
end

% P = P_(t+1)
function u_t = optimal(A, B, P, R, x_t)
    u_t = -inv(B'*P*B+R)*(B'*P*A)*x_t;
end



