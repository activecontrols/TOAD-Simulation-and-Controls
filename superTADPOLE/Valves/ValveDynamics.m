%% Valve Actuator Model for superTADPOLE.

function xdot = ValveDynamics(x, AngleTarget, Mode)

% Define the motor parameters
J_Motor = 3e-5;         % Motor Inertia
D_Motor = 1e-4;         % Motor Damping
Kp_Motor = 50;          % Motor Stiffness
Kd_Motor = 0.04;        % Motor Damping

% Define the gearbox parameters
N = 10;                 % Gear Ratio
B = 3 * pi / 180;       % Half Angle of backlash
Kp_Gear = 1e4;          % Gearbox Stiffness
Kd_Gear = 10;           % Gearbox damping

% Valve Parameters
J_Valve = 9.4e-6;       % Valve ball Inertia
D_Valve = 0.001;        % Valve ball damping
F_Valve = 8;            % Friction limit (8Nm for Cryo Valve)

% Outer Loop Parameters
Kp_Outer = 30;          % Cascaded positional gain
B_Comp = 5 * pi / 180;  % Estimated backlash for use in Controls.

%% Dynamics
% Extract States
MotorPos = x(1);
MotorVel = x(2);
ValvePos = x(3);
ValveVel = x(4);

% Control Modes
if Mode == "Internal"
    % 1. Standard Control: Motor blindly tracks its own scaled target
    MotorTarget = AngleTarget * N;
    T_Motor = Kp_Motor * (MotorTarget - MotorPos) - Kd_Motor * MotorVel;

elseif Mode == "Closed"
    % Cascaded Control: Outer loop reads valve, Inner loop executes
    ValveError = AngleTarget - ValvePos;
    
    % Fast Pacing (Backlash Compensation)
    Backlash_Jump = (B_Comp * N) * tanh(100 * ValveError);
    
    % Final Motor Target = Base + Fast Pacing + Outer Correction
    MotorTarget = (AngleTarget * N) + Backlash_Jump + (Kp_Outer * ValveError);
    T_Motor = Kp_Motor * (MotorTarget - MotorPos) - Kd_Motor * MotorVel;

else
    T_Motor = 0;
end

% Physics
RelDisp = MotorPos / N - ValvePos;
RelVel = MotorVel / N - ValveVel;

% Transmitted torque to valve 
if abs(RelDisp) < B
    T_Valve = 0;
elseif RelDisp > B
    T_Valve = Kp_Gear * (RelDisp - B) + Kd_Gear * RelVel;
else 
    T_Valve = Kp_Gear * (RelDisp + B) + Kd_Gear * RelVel;
end

% Valve Static Friction (Smoothed)
T_Fric = F_Valve * tanh(100 * ValveVel);

%% State Derivatives
xdot = zeros(4,1);
xdot(1) = MotorVel;
xdot(2) = (T_Motor - D_Motor * MotorVel - T_Valve / N) / J_Motor;
xdot(3) = ValveVel;
xdot(4) = (T_Valve - D_Valve * ValveVel - T_Fric) / J_Valve;

end