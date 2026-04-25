%% Valve Actuator Model for superTADPOLE.

function xdot = ValveDynamics(x, AngleTarget, Mode)

% Define the motor parameters
J_Motor = 3e-5;         % Motor Inertia
D_Motor = 1e-4;         % Motor Damping
Kp_Motor = 50;          % Motor Stiffness
Kd_Motor = 0.04;        % Motor Damping
T_Hold = 2;             % Holding Torque for NEMA-23
CornerVel = 40;         % rad/sec (Torque dropoff)

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
Kp_Outer = 12;          % Cascaded positional gain
B_Comp = 2.5 * pi / 180;

%% Dynamics
% Extract States
MotorPos = x(1);
MotorVel = x(2);
ValvePos = x(3);
ValveVel = x(4);
Filter = x(5);

% Control Modes
if Mode == "Internal"
    % 1. Standard Control: Motor blindly tracks its own scaled target
    MotorTarget = AngleTarget * N;
    T_Motor = Kp_Motor * (MotorTarget - MotorPos) - Kd_Motor * MotorVel;
    FDot = 0;

elseif Mode == "Closed"
    % Cascaded Control: Outer loop reads valve, Inner loop executes
    ValveError = AngleTarget - ValvePos;

    % Gearbox disconnect with a deadband
    RawDisconnect = MotorPos / N - ValvePos;
    Deadband = 10 * pi / 1800;
    if abs(RawDisconnect) < Deadband
        Disconnect = 0;
    elseif RawDisconnect > 0
        Disconnect = RawDisconnect - Deadband; 
    else
        Disconnect = RawDisconnect + Deadband; 
    end
    
    % Backlash compensation + filtering
    TargetComp = max(-B_Comp, min(B_Comp, Disconnect));
    FDot = (TargetComp - Filter) / 0.05;

    % Final Motor Target = Base + Fast Pacing + Outer Correction
    MotorTarget = (AngleTarget * N) + (Kp_Outer * ValveError) + 1 * (Filter * N);
    T_Motor = Kp_Motor * (MotorTarget - MotorPos) - Kd_Motor * MotorVel;

else
    T_Motor = 0;
    FDot = 0;
end

% Calculate maximum available torque at current velocity
T_Max_Available = T_Hold / (1 + (abs(MotorVel) / CornerVel));

% Clamp the requested control torque to the motor's physical capability
T_Motor = max(-T_Max_Available, min(T_Max_Available, T_Motor));

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
xdot = zeros(5,1);
xdot(1) = MotorVel;
xdot(2) = (T_Motor - D_Motor * MotorVel - T_Valve / N) / J_Motor;
xdot(3) = ValveVel;
xdot(4) = (T_Valve - D_Valve * ValveVel - T_Fric) / J_Valve;
xdot(5) = FDot;
end