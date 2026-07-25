function MatrixList = RicattiRecursion(TrajectoryList,Q, R, constantsTOAD)
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

P_t = Q;
size(P_t)
size(Q)
size(TrajectoryList.x, 2)
MatrixList = permute(repmat([zeros(12,4)],[1,1,size(TrajectoryList.x,2)]), [3,1,2]);
size(MatrixList)

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


xdot = [qdot;rdot;vdot; omegaBdot;mdot_lox;mdot_ipa];

% Ignore the constrained state
T = [zeros(1,12); eye(12);zeros(2,12)];
T(1:4,1:3) = [zeros(1,3); eye(3)];

A = jacobian(xdot, xn);
B = jacobian(xdot, u);

A = pinv(T) * A * T;
B = pinv(T) * B;

for n = size(TrajectoryList.x,2):-1:2
    A_lin = double(subs(A, [xn; xn1; u], [TrajectoryList.x(:,n); TrajectoryList.x(:,n-1); TrajectoryList.u(:,n-1)]));
    B_lin = double(subs(B, [xn; xn1; u],[TrajectoryList.x(:,n); TrajectoryList.x(:,n-1); TrajectoryList.u(:,n-1)]));
    MatrixList(n, :, :) = gain(A_lin,B_lin,R,P_t)';
    P_t = riccati(A_lin,B_lin,R,Q,P_t);
end
SolveLQR(A_lin, B_lin, Q, R)

end

% Calculates P_(t-1)
function P = riccati(A, B, R, Q, P_t)

P=Q+A.'*P_t*A - (A.'*P_t*B)/(B.'*P_t*B+R)*B.'*P_t*A;
end

function K = gain(A, B, R, P_t)

 K = (B.'*P_t*B+R)\B.'*P_t*A;
end

% P = P_(t+1)
function u_t = optimal(A, B, P, R, x_t)
    u_t = -inv(B'*P*B+R)*(B'*P*A)*x_t;
end

