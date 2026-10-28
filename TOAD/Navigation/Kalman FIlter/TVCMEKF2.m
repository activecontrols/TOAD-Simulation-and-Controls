function [x_est, lastP, logRow] = TVCMEKF2(x_est, constantsASTRA, z, dT, U, mass, P0)
%% TVC-aided M-EKF
% Estimates the first 15 error states; only P(1:15,1:15) is passed in/out.
%
% Inputs
%   x_est  state vector; estimated: q(1:4), p(5:7), v(8:10), b_g(11:13), b_a(14:16)
%          (b_m in 17:19 is read for bias removal but not estimated)
%   z      15x1: [accel(1:3); gyro(4:6); mag(7:9); gps_pos(10:12); gps_vel(13:15)]
%   P0     15x15 initial covariance
%   U      [theta; phi; T; roll] applied gimbal angles, thrust
%   mass   current vehicle mass [kg]
%
% Error state (15): [dtheta(1:3) dp(4:6) dv(7:9) dbg(10:12) dba(13:15)]
 
nx = 15;
zRaw = z;

persistent P lastZ
if isempty(P)
    P = P0;
    P(13:15,1:12) = 0;
    P(1:12,13:15) = 0;
    P(13:15,13:15) = 0.3^2 * eye(3);
    lastZ = zeros(15,1);
end

%%  Thrust model (body frame) 
theta = U(1);
phi = U(2);
T = U(3);

a_thrust = zeros(3,1);

a_thrust(1) = T * cos(theta) * sin(phi) / (constantsASTRA.m_dry+ mass(1) + mass(2));
a_thrust(2) = -T * sin(theta) / (constantsASTRA.m_dry+ mass(1) + mass(2));
a_thrust(3) = T * cos(theta) * cos(phi) / (constantsASTRA.m_dry+ mass(1) + mass(2));

%% Remove bias from IMU
%z(1:3) = z(1:3) - x_est(14:16);
z(4:6) = z(4:6) - x_est(11:13);
z(7:9) = z(7:9) - x_est(17:19);
S = norm(z(7:9));
z(7:9) = z(7:9) / S;

%% Quaternion propagation (gyro)
dx = zeros(nx,1);
q = x_est(1:4);
if sqrt(sum(z(4:6).^2)) < 1e-5 || dT == 0
    x_est(1:4) = q;
else
    x_est(1:4) = (HamiltonianProd(q) * [cos(sqrt(sum(z(4:6).^2)).*dT./2); (z(4:6)/sqrt(sum(z(4:6).^2))).*sin(sqrt(sum(z(4:6).^2)).*dT./2)]);
end
q = x_est(1:4) / norm(x_est(1:4));
 
% A-priori quaternion estimate and rotation matrix
q = q / norm(q);
R_b2i = quatRot(q)';

% GPS velocity lever-arm correction
rGPS = [0 0 0.31]';
z(13:15) = z(13:15) - R_b2i * cross(z(4:6), rGPS);

a_model = a_thrust + x_est(14:16);%without this change ba effectively feeds nothing and velo propagates only on the thrust model


%% State transition matrix
% Thrust-model specific force replaces the accelerometer. The 12-state block
% comes from StateTransitionMat; accel bias is a random walk and does not
% drive velocity (propagation uses the thrust model, not the accelerometer)
F = zeros(nx);
F(1:12, 1:12) = StateTransitionMat(a_model, z(4:6), R_b2i, 0);
F(7:9,13:15) = R_b2i; %fixing yaw required repurposing accel biases as "thrust model error" so this rotation term just puts the error into the inertial frame

% Propagate rest of state using THRUST acceleration
x_est(8:10) = x_est(8:10) + (R_b2i * a_model - [0; 0; constantsASTRA.g]) * dT;
x_est(5:7)  = x_est(5:7) + x_est(8:10) * dT;

% Discrete state transition matrix
Phi = expm(F * dT);
 
% Magnetometer normalization factor
MagMatrix = (eye(3) - z(7:9) * z(7:9)') / S;
 
% Noise matrices
Q = constantsASTRA.Q2(1:nx, 1:nx);

% Thrust scaled velocity noise
aT   = norm(a_thrust);
qThr = aT^2 * (constantsASTRA.tvc.sigDir^2 + constantsASTRA.tvc.sigMag^2) * dT;
Q(7:9,7:9) = Q(7:9,7:9) + qThr*eye(3);

Rc = constantsASTRA.R2;

R = zeros(6);
R(1:3, 1:3) = Rc(1:3, 1:3);% accel vs thrust model
R(4:6, 4:6) = 1e-2 * MagMatrix + 1e-8 * (z(7:9) * z(7:9)');% mag

% Process noise covariance and a-priori propagation step
P = Phi * P * Phi' + Q;
P = (P + P') / 2;

dxI = zeros(nx,1);  dxG = zeros(nx,1);
resI = zeros(6,1);  resG = zeros(6,1);

%%  IMU updates 
if any(lastZ(1:9) ~= zRaw(1:9))
    % Measurement matrix
    % Accel row: prediction is a_thrust
    % so only accel bias is observed. Mag row: attitude only 
    H = zeros(6, nx);
    H(1:3, 13:15) = eye(3);
    H(4:6, 1:3)   = zetaCross(R_b2i' * constantsASTRA.mag);
 
    % Predicted measurements
    z_hat = [a_model;
             R_b2i' * constantsASTRA.mag];
 
    % Kalman gain
    L = P * H' / (H * P * H' + R);

    %L(1:12,1:3) = 0; %makes sure that the accel update only updates the accel bias, not the other states
 
    ILH = (eye(nx) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z([1:3 7:9]) - z_hat);
    inn = L*residual;

    dxI = inn;  resI = residual;

    dx = dx + inn;
end

%%  GPS update 
if any(lastZ(10:15) ~= zRaw(10:15))
    % Measurement matrix
    H = zeros(6, nx);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);
 
    % Measurement covariance matrix
    GyroCovar = eye(3) * 1e-2;
    R = zeros(6);
    R(1:3, 1:3) = 0.3 * eye(3);
    R(4:6, 4:6) = .1 * eye(3) + R_b2i * zetaCross(rGPS) * GyroCovar * (R_b2i * zetaCross(rGPS))';
 
    % Kalman gain
    L = P * H' / (H * P * H' + R);
 
    % Predicted measurements
    z_hat = [x_est(5:7);
             x_est(8:10)];
 
    ILH = (eye(nx) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z(10:15) - z_hat) - H * dx;
    inn = L * residual;

    dxG = inn;  resG = residual;

    dx = dx + inn;
end

logRow = [dxI(1:3)', dxI(10:12)', ...   % 1-6   IMU-step att, bg correction
          dxG(1:3)', dxG(10:12)', ...   % 7-12  GPS-step att, bg correction
          resI(1:3)', resI(4:6)', ...   % 13-18 accel residual, mag residual
          resG(1:3)', resG(4:6)', ...   % 19-24 GPS pos, vel residual
          P(10,13), P(11,14), P(12,15)];% 25-27 cov(dbg_i, dba_i)

% Output
lastP = P;
 
% Update full-state estimates
dq = [1; dx(1:3) / 2];
dq = dq / norm(dq);
 
q_nom = quatmultiply(q', dq');
q_nom = q_nom / norm(q_nom);
x_est(1:4) = q_nom';
x_est(5:16) = x_est(5:16) + dx(4:15);   % pos, vel, b_g, b_a
lastZ = zRaw;
end
