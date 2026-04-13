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
function [U, State_ERR] = TOAD_Control_FCN(PosTarget, X, constantsTOAD, t, MaxVel, VelFF)

% Time Counter
persistent lastT VelErrorI AttErrorI lastAttError
if isempty(lastT)
    lastT = 0;
    VelErrorI = zeros(3,1);
    AttErrorI = zeros(3,1);
    lastAttError = zeros(3,1);
end
dT = t - lastT;
lastT = t;
K_Att = constantsTOAD.K_Att;

% Controller Limits
thrustMax = constantsTOAD.MaxThrust;   %N
gimbalMax = pi/18;
InputBounds = [-gimbalMax       gimbalMax;
               -gimbalMax       gimbalMax;
               .4 * thrustMax   thrustMax;
               -10              10];
U = zeros(4,1);

%% First Loop (P Loop)
    % Position Error Vector
    PosError = PosTarget - X(5:7, :);
    
    % Velocity Command
    K_P = [0.57; 0.57; 0.75];
    VelTarget = K_P .* PosError + VelFF;

    % Velocity Saturation Step
    VelTarget = max(min(VelTarget, MaxVel), -MaxVel);

%% Second Loop (PI Loop)
    % Velocity Error Vector
    VelError = VelTarget - X(8:10, :);

    % Integral Accumulator
    K_I = [0.8; 0.8; 2];
    Leak = 0.05;
    Clamp = [5; 5; 5];

    % Normalize errors (0 to 1 scale)
    MaxAttError = [0.15; 0.15; 0.1];
    MaxVelError = [0.2; 0.2; 0.2];
    NormAttErr = abs(lastAttError) ./ MaxAttError;
    NormVelErr = abs(VelError) ./ MaxVelError;
    
    % Combine errors (Vector magnitude)
    TotalErrorMetric = NormAttErr + NormVelErr;
    
    % Calculate Gate using Gaussian function
    Gate = (1 - Leak) * exp(-2 * TotalErrorMetric.^2) + Leak;

    % Integrator
    K_I = K_I .* Gate;
    VelErrorI = VelErrorI + K_I .* VelError .* dT;
    VelErrorI = max(min(VelErrorI, Clamp), -Clamp);
    K_P = [1.7; 1.7; 7];

    % Acceleration Target
    AccelTarget = K_P .* VelError + VelErrorI  + [0; 0; constantsTOAD.g];

    % Acceleration Saturation Step
    MaxAccelUp = [2 2 15]';
    MaxAccelDown = [-2 -2 4]';
    AccelTarget = max(min(AccelTarget, MaxAccelUp), MaxAccelDown);

%% Kinematics Step
    % Compute thrust target (Update to use estimated mass, for now I gain takes care)
    TargetForce_I = constantsTOAD.m_dry * AccelTarget;
    U(3) = norm(TargetForce_I);
    
    % Compute target attitude via GSP.
    AccelTarget(3) = max(AccelTarget(3), constantsTOAD.g);
    Z_b = AccelTarget / norm(AccelTarget);

    % Heading reference (+X axis rolled to north)
    HDGRef = [1; 0; 0];
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
    Clamp = [0.4; 0.4; 0.4];
    AttErrorI = AttErrorI + AttError(2:4) .* dT;
    AttErrorI = max(min(AttErrorI, Clamp), -Clamp);

    % State vector and error
    X_Err = [-AttError(2:4, :); X(11:13, :); AttErrorI];
    State_ERR = [AttError(2:4); PosError; VelError].^2;

    % LQR Controller
    U([1 2 4]) = -K_Att * X_Err;

%% Controls Saturation
uMax = InputBounds(:, 2);
uMin = InputBounds(:, 1);
U = min(max(U, uMin), uMax);
    

