function [x_est, lastP] = TVCMEKF2(x_est, constantsASTRA, z, dT, U, mass)
%% M-EKF with thrust-model propagation, force-error state d_f, and
%% TVC-aided M-EKF (15-state replacement for FlightEstimator2)
% Estimates the first 15 error states; only P(1:15,1:15) is passed in/out.
%   [x_est, lastP(1:15,1:15)] = TVCMEKF2(x_est, TOAD, z, dt, lastP(1:15,1:15), U, mass);
%
% Inputs
%   x_est  state vector; estimated: q(1:4), p(5:7), v(8:10), b_g(11:13), b_a(14:16)
%          (b_m in 17:19 is read for bias removal but not estimated)
%   z      15x1: [accel(1:3); gyro(4:6); mag(7:9); gps_pos(10:12); gps_vel(13:15)]
%   P0     15x15 initial covariance
%   U      [theta; phi; T; ...] applied gimbal angles [rad], thrust [N]
%   mass   current vehicle mass [kg]
%
% Error state (15): [dtheta(1:3) dp(4:6) dv(7:9) dbg(10:12) dba(13:15)]
 
nx = 15;

persistent P lastZ
if isempty(P)
    P = P0;
    lastZ = zeros(15,1);
end

%%  Thrust model (body frame) 
theta = U(1);  phi = U(2);  T = U(3);
T_b = T * [cos(theta)*sin(phi);
          -sin(theta);
           cos(theta)*cos(phi)];
a_thrust = T_b / mass;

%% Remove bias from IMU
z(1:3) = z(1:3) - x_est(14:16);
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


%% State transition matrix
% Thrust-model specific force replaces the accelerometer. The 12-state block
% comes from StateTransitionMat; accel bias is a random walk and does not
% drive velocity (propagation uses the thrust model, not the accelerometer)
F = zeros(nx);
F(1:12, 1:12) = StateTransitionMat(a_thrust, z(4:6), R_b2i, 0);


% Propagate rest of state using THRUST acceleration
x_est(8:10) = x_est(8:10) + (R_b2i * a_thrust - [0; 0; constantsASTRA.g]) * dT;
x_est(5:7)  = x_est(5:7) + x_est(8:10) * dT;

% Discrete state transition matrix
Phi = expm(F * dT);
 
% Magnetometer normalization factor
MagMatrix = (eye(3) - z(7:9) * z(7:9)') / S;
 
% Noise matrices
Q = constantsASTRA.Q(1:nx, 1:nx);
 
Rc = constantsASTRA.R;
R = zeros(6);
R(1:3, 1:3) = Rc(1:3, 1:3);                                   % accel vs thrust model
R(4:6, 4:6) = 1e-2 * MagMatrix + 1e-8 * (z(7:9) * z(7:9)');   % mag
 
% Process noise covariance and a-priori propagation step
P = Phi * P * Phi' + Q;
P = (P + P') / 2;

%%  IMU updates 
if any(lastZ(1:9) ~= zRaw(1:9))
    % Measurement matrix
    % Accel row: prediction is a_thrust
    % so only accel bias is observed. Mag row: attitude only (no mag bias state)
    H = zeros(6, nx);
    H(1:3, 13:15) = eye(3);
    H(4:6, 1:3)   = zetaCross(R_b2i' * constantsASTRA.mag);
 
    % Predicted measurements
    z_hat = [a_thrust;
             R_b2i' * constantsASTRA.mag];
 
    % Kalman gain
    L = P * H' / (H * P * H' + R);
 
    ILH = (eye(nx) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z([1:3 7:9]) - z_hat);
    inn = L*residual;
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
    R(1:3, 1:3) = 0.1 * eye(3);
    R(4:6, 4:6) = 0.3 * eye(3) + R_b2i * zetaCross(rGPS) * GyroCovar * (R_b2i * zetaCross(rGPS))';
 
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
    dx = dx + inn;
end


% Output
lastP = P;
 
% Update full-state estimates
dq = [1; dx(1:3) / 2];
dq = dq / norm(dq);
 
q_nom = quatmultiply(q', dq');
q_nom = q_nom / norm(q_nom);
x_est(1:4) = q_nom';
x_est(5:16) = x_est(5:16) + dx(4:15);   % pos, vel, b_g, b_a
lastZ = z;
end
