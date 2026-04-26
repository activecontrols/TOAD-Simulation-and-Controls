function [PropMass, FlightTime] = SKIPPER_3DoF(MaxThrustINPUT, WetMass)

%Variable set-up
syms r_1 r_2 v_1 v_2 theta theta_dot m m_dot thrust alpha_ang
%Constants
syms g init_m radius lengthTOAD max_thrust
%Time variable 
syms t 

%Declare constants 
consts = transpose([g init_m radius lengthTOAD max_thrust]);
r = [r_1; r_2];
v = [v_1; v_2];
dictionary = [r,     v;
              theta, theta_dot;
              m,    m_dot];

%Set up assumptions for symbolic variables
assume(imag(dictionary) == 0);
assume(imag(consts) == 0);
assume(consts > 0);
assume(thrust > 0);
assume(alpha_ang > 0);
assume(m > 0);

%State Vector
x = [r; v; theta; theta_dot; m];

%Input vector
u = [thrust; alpha_ang];

%BF Thrust
BF = thrust*[sin(alpha_ang); cos(alpha_ang)];

%Rotation matrix from Inertial frame to Body frame
I2B = [cos(theta), -sin(theta);
       sin(theta),  cos(theta)];

%Rotation matrix from Body frame to Inertial frame
B2I = [cos(-theta), -sin(-theta);
       sin(-theta),  cos(-theta)];

%Inertial frame net force
IF = B2I*BF - [0; m*g];

%Derivative of position
r_dot = v;

%Derivative of velocity
v_dot = IF/m;

%Derivative of mass
Isp_MAX = 184.77; % (s)
Isp_MIN = 159.29; % (s)
Isp = Isp_MIN + 2 * (thrust / max_thrust - 0.5) * (Isp_MAX - Isp_MIN);
m_dot = -thrust / (Isp * g);

%Torque
tau = lengthTOAD*thrust*sin(alpha_ang);

%Moment of inertia
I = (1/3)*m*radius^2+(1/4)*m*2*lengthTOAD^2;

%Derivative of angular velocity
theta_ddot = tau/I;

%Derivative of state vector
x_dot = [r_dot; v_dot; theta_dot; theta_ddot; m_dot];

%Create numerical functions
matlabFunction(x_dot, 'File', './odefcn2.m', 'Vars', {t, x, u, consts});
%%
%Define constants and initial conditions

g = 9.8066;                         % m/s^2
rad = 0.0254;                       % m
len = 4;                            % m
constants = [g; WetMass; rad; len; MaxThrustINPUT];
x_0 = [0; 0; 0; 0; 0; 0; WetMass];

%Calculate jacobians for X and U
J_x = jacobian(x_dot, x);
J_u = jacobian(x_dot, u);

%System critical points
x_crit = [0;0;0;0;0;0;WetMass];
u_crit = [constants(2)*constants(1);0];

%Generate A and B matrices using jacobians, constants, and critical points
A = double(subs(J_x, [x; u; consts], [x_crit; u_crit; constants]));
B = double(subs(J_u, [x; u; consts], [x_crit; u_crit; constants]));

%C and D matrices
C = eye(size(x, 1));
D = 0;

%Set up LQR controller
a_weights = [5; 8; 1; 5; 2; 2; 0];
b_weights = [1; 1];
maxA = 0.8*g;

a_weights = a_weights / norm(a_weights);
b_weights = b_weights / norm(b_weights);

max_x = [2000, 2000, 0.5, 2, pi/6, pi/6, 1000];
max_u = [(maxA / g + 1) * WetMass* 0.5 * g , pi/24];

Q = eye(size(x,1)) .* a_weights ./ max_x.^2;
R = eye(size(u,1)) .* b_weights ./ max_u.^2;

% Gain Matrix
[K, ~, ~] = lqr(A, B, Q, R);

%%
% Simulation
clear ref_generator;
clear inputfcn;
tspan = [0 60];

% Solve using ODE113 and nonlinear dynamics
[tsim, xsim] = ode45(@(tsim, xsim) odefcn2(tsim, xsim, inputfcn(K, xsim, tsim, MaxThrustINPUT), constants), tspan, x_0);

% Solve for inputs
u = zeros(size(xsim,1 ), 2);
clear ref_generator;
clear inputfcn;
for i = 1:1:size(xsim, 1)
    x_i = xsim(i,1:7).';
    t_i = tsim(i);
    u(i,1:2) = inputfcn(K, x_i, t_i, MaxThrustINPUT);
end
%% OUTPUTS
EndIndex = find(xsim(100:end,2) < 0.1 & abs(xsim(100:end,4)) < 10);
if isempty(EndIndex)
    PropMass = 0;
    FlightTime = 0;
    ThrustDev = 0;
else
    EndIndex = EndIndex(1) + 100;
    PropMass = WetMass - xsim(EndIndex, 7);
    FlightTime = tsim(EndIndex);
    AccelVert = (xsim(2:end,2) - xsim(1:end-1,2)) ./ (tsim(2:end) - tsim(1:end-1));
    AccelLat = (xsim(2:end,1) - xsim(1:end-1,1)) ./ (tsim(2:end) - tsim(1:end-1));
    AccelTot = (AccelVert.^2 + AccelLat.^2).^(1/2);

    ThrustDev = u(1:EndIndex,1) - 0.7 * MaxThrustINPUT;
    ThrustDev = sum(ThrustDev.^2);
end



