function [x, u, xdot] = FlightDynamicsGen(constants, mode)
% Symbolic defiitions
syms r1 r2 r3               % earth frame position (R3 up, NWU frame)
syms v1 v2 v3               % earth frame velocity
syms q0 q1 q2 q3            % earth-body quaternion
syms omega1 omega2 omega3   % earth-body angular velocity
syms m_lox m_ipa            % Propellant masses
syms theta phi;             % thrust angles
syms thrust                 % thrust magnitude
syms roll                   % roll input (Nm)
syms m_dry l g rTB          % system constants
syms m_lox_i m_ipa_i        % constants
syms OxHeight FuHeight      % constants
syms OxRadius FuRadius      % constants
syms Ox_Z Fu_Z              % constants
syms MaxThrust MaxMdot OF   % constants

syms J [3 3]                % Dry Intertia Matrix

% Total Mass
m = m_dry + m_lox + m_ipa;

%Earth-frame position & velocity
r = [r1; r2; r3];
v = [v1; v2; v3];

%Full quaternion vector
q = [q0; q1; q2; q3];

%DCM's
C_BI = quatRot(q);
C_IB = C_BI.';

%Angular velocities
omegaB = [omega1; omega2; omega3];

%Body frame thrust vector
TB = thrust * [cos(theta)*sin(phi); -sin(theta); cos(theta)*cos(phi)];

%Earth frame net force
FI = C_IB * TB + [0; 0;-m*g];

%Rdot and Vdot
rdot = v;
vdot = FI / m;

% Tank drain dynamics
mdot_lox = -thrust / MaxThrust * OF / (1 + OF) * MaxMdot;
mdot_ipa = -thrust / MaxThrust * 1 / (1 + OF) * MaxMdot;

% Propellant Fill height
OxFluidHeight = (m_lox / m_lox_i) * OxHeight * 0.9;
FuFluidHeight = (m_ipa / m_ipa_i) * FuHeight * 0.9;

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
d_dry = rTB - CGz;
d_lox = OxFluidLocation - CGz;
d_ipa = FuFluidLocation - CGz;

% Shifted Inertias
J_dry = J + m_dry * diag([d_dry^2, d_dry^2, 0]);
J_lox = J_lox + m_lox * diag([d_lox^2, d_lox^2, 0]);
J_ipa = J_ipa + m_ipa * diag([d_ipa^2, d_ipa^2, 0]);

% Bring everything to instantaneous CG (w.r.t Engine attachment frame)
J_tot = J_dry + J_lox + J_ipa;

% Body frame moment (Roll torque along the thrust vector axis due to contra
% EDF) (TODO: UPDATE for use with RCS)
% Off center moments for Simulation
TBx = -0.008 * (mode == 1);
TBy = 0.0120 * (mode == 1);
MB = zetaCross([TBx; TBy; -CGz])*TB + (TB * roll) / thrust;

% Dynamics
qdot = 0.5 * HamiltonianProd(q) * [0; omegaB];
omegaBdot = J_tot \ (MB - zetaCross(omegaB) * J_tot * omegaB);

%State vector
x = [q; r; v; omegaB; m_lox; m_ipa];

%Input vector (Assuming concentric EDF design)
u = [theta; phi; thrust; roll];
    
%State derivative
xdot = [qdot; rdot; vdot; omegaBdot; mdot_lox; mdot_ipa];

% Create symbolic and numeric constant arrays
% Perturb true m, rTB, and J values (TODO: Update with disturbance vector &
% Monte Carlo workflow)
if mode == 1
    constants.rTB = constants.rTB;
    constants.m = constants.m * 0.97;
    constants.J(:) = constants.J(:) * 1.3;
end

constVal = [constants.m; constants.l; constants.g; constants.rTB; constants.J(:)];
constVec = [m_dry; l; g; rTB; J(:)];

% Subsitute numeric values of constants into xdot
xdot = subs(xdot, constVec, constVal);
end
