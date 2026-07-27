function MatrixList = RicattiRecursion(Trajectory,Q, R, constantsTOAD)
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

% state variables for x_n-1
syms q01 q11 q21 q31
syms r11 r21 r31
syms v11 v21 v31
syms omega11 omega21 omega31
syms m_lox1 m_ipa1

MatrixList = permute(repmat([zeros(15,4)],[1,1,size(Trajectory.x,2)]), [3,1,2]);

q = [q0;q1;q2;q3];
r = [r1;r2;r3];
v = [v1;v2;v3];
m = [m_lox; m_ipa];
omegaB = [omega1; omega2; omega3];

xn = [q;r;v; omegaB; m]; % x_n

% U
u = [theta;phi;thrust;roll];

% X_n-1
qn1 = [q01;q11;q21;q31];
rn1 = [r11;r21;r31];
vn1 = [v11;v21;v31];
mn1 = [m_lox1; m_ipa1];
omegaBn1 = [omega11;omega21;omega31];
xn1 = [qn1;rn1;vn1;omegaBn1;mn1]; % x_n-1

%% Dynamics
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
m_dry = constantsTOAD.m_dry;
m = m_dry + m_lox + m_ipa;

% Propellant Fill height
OxFluidHeight = (m_lox / constantsTOAD.OxMass) * constantsTOAD.OxHeight * 0.9;
FuFluidHeight = (m_ipa / constantsTOAD.FuMass) * constantsTOAD.FuHeight * 0.9;

% Propellant inertias
J_xx = 1/12 * m_lox * (3 * constantsTOAD.OxRadius^2 + OxFluidHeight^2);
J_zz = 1/2 * m_lox * constantsTOAD.OxRadius^2;
J_lox = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

J_xx = 1/12 * m_ipa * (3 * constantsTOAD.FuRadius^2 + FuFluidHeight^2);
J_zz = 1/2 * m_ipa * constantsTOAD.FuRadius^2;
J_ipa = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

% Fluid fills & location of CoM
OxFluidLocation = constantsTOAD.Ox_Z + OxFluidHeight / 2;
FuFluidLocation = constantsTOAD.Fu_Z + FuFluidHeight / 2;
CGz = (constantsTOAD.m_dry * constantsTOAD.rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

% Distances to CG
d_dry = constantsTOAD.rTB - CGz;
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Shifted Inertias
J_dry = constantsTOAD.J + m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);
% Bring everything to instantaneous CG (w.r.t Engine attachment frame)
J_tot = J_dry + J_lox + J_ipa;

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
A_fcn = matlabFunction(A, 'Vars', {xn, xn1, u});
B_fcn = matlabFunction(B, 'Vars', {xn, xn1, u});

for n = size(Trajectory.x,2):-1:2
    % Extract states
    x_n   = Trajectory.x(:, n);
    x_n1  = Trajectory.x(:, n-1);
    u_n1  = Trajectory.u(:, n-1);
    dT = Trajectory.t(n) - Trajectory.t(n-1);
    
    % Jacobian Eval
    A_lin = A_fcn(x_n, x_n1, u_n1);
    B_lin = B_fcn(x_n, x_n1, u_n1);

    % Kinematic mapping 
    T = zeros(15, 12);
    q_ref = x_n(1:4);
    T(1:4, 1:3) = 0.5 * XiMat(q_ref);
    T(5:13, 4:12) = eye(9);

    % Reduce jacobians
    A_lin = pinv(T) * A_lin * T;
    B_lin = pinv(T) * B_lin;

    % Matrix augmentation for LQI
    C_mat = zeros(3, 12);       C_mat(1:3, 4:6) = eye(3);
    A_aug = [A_lin, zeros(12, 3); 
             C_mat, zeros(3, 3)];         
    B_aug = [B_lin; 
             zeros(3, 4)];
    
    %% Matrix discretization using ZOH
    % Dimensions
    nx = size(A_aug, 1); 
    nu = size(B_aug, 2); 
    
    % Construct the continuous block matrix
    M_c = [A_aug, B_aug; 
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

    % Matrix Eval and Updating Ricatti Cost
    MatrixList(n, :, :) = gain(A_d,B_d,R,P_t)';
    P_t = riccati(A_d,B_d, R , Q , P_t);
end

% Copy over 2nd to last gain matrix to first spot
MatrixList(1, :, :) = MatrixList(2, :, :);

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



