% Version 2 Controller formulation for TOAD. Structure consists of 3
% cascaded loops.
%
% First Loop is a P Loop in charge of Pos -> Vel commands.
%
% Second loop is a PI Loop in charge of Vel -> Acceleration commands, fitted
% with anti-windup and soft-gatingcbased on attitude error.
%
% Second loop is an LQRi Loop in charge of Target Attitude -> Gimbal +
% Torque commands. Integral action is fitted with anti-windup clamps and is
% self soft-gated.
% Tuning of the second loop is done by using a Genetic algorithm to
% maximize a Crossover Frequency vs. Disk Margin tradeoff on all actuator
% channels.
%
% By: Pablo Plata   -   11/27/25 (Happy Thanksgiving!)
function [U, State_ERR] = TOAD_Control_FCN(PosTarget, X, constantsTOAD, t, MaxVel, VelFF, HDGRef, GND)

% Time Counter
persistent lastT VelErrorI AttErrorI lastAttError lastAccelZ LatTrim
if isempty(lastT)
    lastT = 0;
    VelErrorI = zeros(3,1);
    AttErrorI = zeros(3,1);
    lastAttError = zeros(3,1);
    LatTrim = zeros(3,1);
end
dT = t - lastT;
lastT = t;
K_Att = constantsTOAD.K_Att_Dry + (constantsTOAD.K_Att_Wet - constantsTOAD.K_Att_Dry) * (X(14) / constantsTOAD.OxMass);

% Controller Limits
thrustMax = constantsTOAD.MaxThrust;   %N
gimbalMax = pi/18;
InputBounds = [-gimbalMax       gimbalMax;
               -gimbalMax       gimbalMax;
               .4 * thrustMax   thrustMax;
               -7               7];
U = zeros(4,1);

%% First Loop (P Loop)
    % Position Error Vector
    PosError = PosTarget - X(5:7, :);
    
    % Velocity Command
    K_P = [0.43; 0.43; 0.8];
    VelTarget = K_P .* PosError + VelFF;

    % Velocity Saturation Step
    VelTarget = max(min(VelTarget, MaxVel), -MaxVel);

%% Second Loop (PI Loop)
    % Velocity Error Vector
    VelError = VelTarget - X(8:10, :);

    % Steady-State Trim Observer (Handles MEKF/CG bias)
    SteadyState = norm(VelTarget(1:2)) < 0.5 && norm(lastAttError(2:3)) < 0.1;
    K_Trim = [0.05; 0.05; 0];
    
    LatTrim = LatTrim + (K_Trim .* VelError .* dT .* (1 - GND) .* SteadyState);
    LatTrim(1:2) = max(min(LatTrim(1:2), 0.7), -0.7);

    % Dynamic Integrator Gating (Lorentzian)
    Leak = 0.1;
    GateWidth = 0.5;
    ErrorMag = norm(VelError) + 5.0 * norm(lastAttError(2:3));
    Gate = Leak + (1 - Leak) * (1 / (1 + (ErrorMag / GateWidth)^2));

    % Dynamic Integrator
    K_I = [0.25; 0.25; 0.3];
    Clamp = [5; 5; 5];
    
    K_I = K_I .* Gate;
    VelErrorI = VelErrorI + K_I .* VelError .* dT .* (1 - GND);
    VelErrorI = max(min(VelErrorI, Clamp), -Clamp);
    
    % Acceleration Target
    K_P = [2.4; 2.4; 2.0];
    AccelTarget = K_P .* VelError + VelErrorI + [0; 0; constantsTOAD.g];

    % Acceleration Saturation Step
    MaxAccelUp = [2.3 2.3 15]';
    MaxAccelDown = [-2.3 -2.3 4]';
    AccelTarget = max(min(AccelTarget, MaxAccelUp), MaxAccelDown);

    % Acceleration Rate Limit
    MaxJerkZ = 10;
    if isempty(lastAccelZ)
        lastAccelZ = AccelTarget(3);
    end
    MaxDeltaZ = MaxJerkZ * dT;
    DeltaZ = AccelTarget(3) - lastAccelZ;
    AccelTarget(3) = lastAccelZ + max(min(DeltaZ, MaxDeltaZ), -MaxDeltaZ);
    lastAccelZ = AccelTarget(3);

    % Inject Trim Acceleration from the Low Rate Integrator, scaled to
    % maintain physical tilt
    TrimAccel = LatTrim .* (AccelTarget(3) / constantsTOAD.g);
    AccelTarget(1:2) = AccelTarget(1:2) + TrimAccel(1:2);
    
%% Kinematics Step
    % Compute thrust target (Update to use estimated mass, for now I gain takes care)
    TargetForce_I = (constantsTOAD.m_dry + X(14) + X(15)) * AccelTarget;
    U(3) = norm(TargetForce_I);
    
    % Compute target attitude via GSP.
    Z_b = AccelTarget / norm(AccelTarget);

    % Heading reference (+X axis rolled to north)
    Y_b = cross(Z_b, HDGRef);
    Y_b = Y_b / norm(Y_b);

    % Complete the triad
    X_b = cross(Y_b, Z_b);

    % Create DCM and convert to quaternion
    DCM = [X_b Y_b Z_b];
    TargetAtt = DCM_Quat_Conversion(DCM);

%% Third Loop (LQRi)
    % Attitude Error computation
    Q_Conj = [X(1); -X(2:4, :)];
    AttError = HamiltonianProd(Q_Conj) * TargetAtt;
    if AttError(1) < 0
        AttError = -AttError;
    end
    lastAttError = AttError(2:4);

    % Error accumulation and clamping
    Clamp = [0.0; 0.0; 0.4];
    AttErrorI = AttErrorI + AttError(2:4) .* dT .* (1 - GND);
    AttErrorI = max(min(AttErrorI, Clamp), -Clamp);

    % State vector and error
    State_ERR = [AttError(2:4); PosError; VelError].^2;
    MaxRollErr = 0.25;
    if abs(AttError(4)) > MaxRollErr
        AttError(4) = sign(AttError(4)) * MaxRollErr;
    end
    X_Err = [-AttError(2:4, :); X(11:13, :); AttErrorI];
    
    % LQR Controller
    U([1 2 4]) = -K_Att * X_Err;

%% Controls Saturation
uMax = InputBounds(:, 2);
uMin = InputBounds(:, 1);
U = min(max(U, uMin), uMax);
    

