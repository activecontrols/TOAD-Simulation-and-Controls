addpath('C:\MATLAB Tools\casadi-3.7.2-windows64-matlab2018b')
import casadi.*

% States
q = MX.sym('q', 4);
r = MX.sym('r', 3);
v = MX.sym('v', 3);
omegaB = MX.sym('omegaB', 3);
m_lox = MX.sym('m_lox', 1);
m_ipa = MX.sym('m_ipa', 1);

% Constants
J = MX.sym('J', 3, 3);
J_d = MX.sym('J_d', 3, 3);
TB_d = MX.sym('TB_d', 3);
m_dry = MX.sym('m_dry', 1);
g = MX.sym('g');
rTB = MX.sym('rTB');
MaxThrust = MX.sym('MaxThrust');
MaxMdot = MX.sym('MaxMdot');
OF = MX.sym('OF');
OxMassI = MX.sym('OxMassI');
FuMassI = MX.sym('FuMassI');
OxHeight = MX.sym('OxHeight');
FuHeight = MX.sym('FuHeight');
OxRadius = MX.sym('OxRadius');
FuRadius = MX.sym('FuRadius');
Ox_Z = MX.sym('Ox_Z');
Fu_Z = MX.sym('Fu_Z');
 % = MX.sym('OF');
 % = MX.sym('OF');

% Monte Carlo Params
MaxMdot_d = MX.sym('MaxMdot_d');

% Controls
theta = MX.sym('theta', 1);
phi = MX.sym('phi', 1);
thrust = MX.sym('thrust', 1);
roll = MX.sym('roll', 1);

% Total Mass
m = m_dry + m_lox + m_ipa;

%DCM's
C_BI = quatRot(q);
C_IB = C_BI.';

%Body frame thrust vector
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

%Earth frame net force
FI = C_IB * TB + [0; 0;-m*g];

%Rdot and Vdot
rdot = v;
vdot = FI / m;

% Tank drain dynamics
mdot_lox = -thrust / MaxThrust * OF / (1 + OF) * (MaxMdot + MaxMdot_d);
mdot_ipa = -thrust / MaxThrust * 1 / (1 + OF) * (MaxMdot + MaxMdot_d);

% Propellant Fill height
OxFluidHeight = (m_lox / OxMassI) * OxHeight * 0.9;
FuFluidHeight = (m_ipa / FuMassI) * FuHeight * 0.9;

% Propellant inertias
J_xx = 1/12 * m_lox * (3 * OxRadius^2 + OxFluidHeight^2);
J_zz = 1/2 * m_lox * OxRadius^2;
J_lox = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

J_xx = 1/12 * m_ipa * (3 * FuRadius^2 + FuFluidHeight^2);
J_zz = 1/2 * m_ipa * FuRadius^2;
J_ipa = [J_xx, 0, 0;
         0, J_xx, 0;
         0,    0, J_zz];

% Fluid fills & location of CoM
OxFluidLocation = Ox_Z + OxFluidHeight / 2;
FuFluidLocation = Fu_Z + FuFluidHeight / 2;
CGz = (m_dry * rTB + m_lox * OxFluidLocation + m_ipa * FuFluidLocation) / m;

% Distances to CG
d_dry = rTB - CGz + TB_d(3);
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Shifted Inertias
J_dry = J + m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);

% Bring everything to instantaneous CG (w.r.t Engine attachment frame)
J_tot = J_dry + J_lox + J_ipa + J_d;

% Body frame moment (Roll torque along the thrust vector axis due to contra
% EDF) (TODO: UPDATE for use with RCS)
% Off center moments for Simulation
MB = zetaCross([0; 0; -CGz] + TB_d)*TB + (TB * roll) / thrust;

% Dynamics
qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
netTau = (MB - zetaCross(omegaB) * J_tot * omegaB);

% Function inputs
x = [q; r; v; omegaB; m_lox; m_ipa];
u = [theta; phi; thrust; roll];
params = [m_dry; g; rTB; Ox_Z; OxMassI; OxHeight; Fu_Z; FuMassI; FuHeight; 
    J(:); OxRadius; FuRadius; MaxThrust; OF; MaxMdot; MaxMdot_d; J_d(:); TB_d];

%State derivatives
xdot = [qdot; rdot; vdot; netTau; mdot_lox; mdot_ipa];
dynamics_fnc = Function('dynamics_fnc', {x, u, params}, {xdot}); 