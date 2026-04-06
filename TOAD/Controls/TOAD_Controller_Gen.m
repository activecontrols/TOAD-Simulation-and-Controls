function [K, lin] = TOAD_Controller_Gen(constants, OxMass, FuMass)
%% Generation script for Toad Controls (LQRi Step, Third Loop)
% Symbolic Variables
syms q0 q1 q2 q3            % earth-body quaternion
syms omega1 omega2 omega3   % earth-body angular velocity
syms theta phi;             % thrust angles
syms thrust                 % thrust magnitude
syms roll                   % roll input (rad/sec^2)
syms m g rTB              % system constants
syms J [3 3]                % Intertia Matrix

%% Equations of Motion
%Full quaternion vector
q = [q0; q1; q2; q3];

%Angular velocity (taken wrt E-frame, represented in B-frame) 
omegaB = [omega1; omega2; omega3];

%Body frame thrust vector
T_vec = [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];
TB = thrust * T_vec;

% Total mass at linearization point
m_tot = constants.m_dry + OxMass + FuMass;

% Fluid heights at linearization fill level
OxFluidHeight = (OxMass / constants.OxMass) * constants.OxHeight * 0.9;
FuFluidHeight = (FuMass / constants.FuMass) * constants.FuHeight * 0.9;

% Propellant CoM locations
OxFluidLocation = constants.Ox_Z + OxFluidHeight / 2;
FuFluidLocation = constants.Fu_Z + FuFluidHeight / 2;

% CG at linearization point
CGz = (constants.m_dry * constants.rTB + ...
       OxMass * OxFluidLocation + ...
       FuMass * FuFluidLocation) / m_tot;

% Propellant inertias
Jxx_lox = 1/12 * OxMass * (3*constants.OxRadius^2 + OxFluidHeight^2);
Jzz_lox = 1/2  * OxMass * constants.OxRadius^2;
J_lox = diag([Jxx_lox, Jxx_lox, Jzz_lox]);

Jxx_ipa = 1/12 * FuMass * (3*constants.FuRadius^2 + FuFluidHeight^2);
Jzz_ipa = 1/2  * FuMass * constants.FuRadius^2;
J_ipa = diag([Jxx_ipa, Jxx_ipa, Jzz_ipa]);

% Parallel-axis shift to linearization-point CG
d_dry = constants.rTB   - CGz;
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Compute inertias
J_dry = constants.J + constants.m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox       + OxMass          * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa       + FuMass          * diag([d_ipa^2, d_ipa^2, 0]);
J_tot = J_dry + J_lox + J_ipa;

% Dynamics
% Torque vector
MB = zetaCross([0; 0; -CGz])*TB + (T_vec * roll);
M = [q(1) -q(2) -q(3) -q(4);
     q(2)  q(1) -q(4)  q(3);
     q(3)  q(4)  q(1) -q(2);
     q(4) -q(3)  q(2)  q(1)];
qdot = 0.5 * M * [0; omegaB];
omegaBdot = J_tot \ (MB - zetaCross(omegaB) * J_tot * omegaB);

% State vectors
x = [q; omegaB];
xdot = [qdot; omegaBdot];
u = [theta; phi; roll];

% Substitude constants in
constVal = [m_tot; constants.g; CGz; J_tot(:); m_tot * constants.g];
constVec = [m; g; rTB; J(:); thrust];

% Subsitute numeric values of constants into xdot
xdot = subs(xdot, constVec, constVal);

%% Linearization and Controller Generaiton
lin.A = jacobian(xdot, x);
lin.B = jacobian(xdot, u);

% Map to 6 states.
T = [zeros(1,6); eye(6)];
T(1:4,1:3) = [zeros(1,3); eye(3)];
lin.A = pinv(T) * lin.A * T;
lin.B = pinv(T) * lin.B;

% Augment system with Error-Integral
lin.A = [lin.A zeros(6,3);
         -eye(3) zeros(3,6)];
lin.B = [lin.B; zeros(3)];

% Linearization point
delx = [1; 0; 0; 0; 0; 0; 0];
delu = [0; 0; 0];
lin.A = double(subs(lin.A, [x; u], [delx; delu]));
lin.B = double(subs(lin.B, [x; u], [delx; delu]));

% Hand tuning for Q for now
a_weights = ones(6,1);
a_weights = a_weights / norm(a_weights);
max_x = [0.26, 0.26, 0.25, 18, 18, 1.0];
Q = eye(6) .* a_weights ./ max_x.^2;
R = diag([18, 18, 2]);

% Augment Q with integral states
Qi = diag([1.8, 1.8, 5]);
Q = [Q zeros(6,3);
     zeros(3,6) Qi];

% Initial Controller Solution
K = lqr(lin.A, lin.B, Q, R);
end