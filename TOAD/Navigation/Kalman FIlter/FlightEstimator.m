function [x_est, lastP] = FlightEstimator(x_est,constantsASTRA,z,dT,P0)
%% M-EKF Implementation
% Remove bias from IMU
z(1:3) = z(1:3) - x_est(14:16);
z(4:6) = z(4:6) - x_est(11:13);
z(7:9) = z(7:9) - x_est(17:19);

% Extract quaternion
dx = zeros(9,1);
q = x_est(1:4);
qdot = 0.5 * HamiltonianProd(q) * [0; z(4:6)]; 
x_est(1:4) = q + qdot * dT;
q = x_est(1:4);

% A-priori quaternion estimate and rotation matrix
q = q / norm(q);
R_b2i = quatRot(q)';

% GPS velocity correction
rGPS = [0 0 0.31]';
z(13:15) = z(13:15) - R_b2i * cross(z(4:6), rGPS);

% Process Covariance Matrix
persistent P lastZ 
if isempty(P)
    P = P0;  
    lastZ = zeros(15,1);
end

% State Transition Matrix
F = StateTransitionMat(z(1:3), z(4:6), R_b2i, 0);

% Propagate rest of state using IMU
x_est(8:10) = x_est(8:10) + (R_b2i * z(1:3) - [0; 0; constantsASTRA.g]) * dT;
x_est(5:7) = x_est(5:7) + x_est(8:10) * dT;

% Discrete STM
Phi = expm(F * dT);

% Extract Matrices
Q = constantsASTRA.Q(1:9,1:9);

% Process Noise Covariance and a-priori propagation step
Q = 0.5 * Q;
P = Phi * P * Phi' + Q;
RTK = 1;

%% GPS Update
if sum(lastZ(10:15) - z(10:15)) ~=0

    % Measurement matrix
    H = zeros(6,9);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);

    % Measurement Covariance Matrix
    gps_pos_covar = 0.2 * RTK + 10 * (1 - RTK);
    gps_vel_covar = 0.05; %gps_pos_covar * 5;
    R = diag([gps_pos_covar^2 * ones(3,1); gps_vel_covar^2 * ones(3,1)]);

    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);

    % Predicted measurements 
    z_hat = [x_est(5:7);
             x_est(8:10)];

    % Update Error State
    ILH = (eye(9) - L * H);
    P = ILH * P * ILH' + L * R * L';
    residual = (z(10:15) - z_hat);
    inn = L * residual;
    dx = dx + inn;
end

% Update full-state estimates
dq = [1; dx(1:3) / 2];
dq = dq / norm(dq);

q_nom = quatmultiply(q', dq');
q_nom = q_nom / norm(q_nom); 
x_est(1:4) = q_nom';
x_est(5:10) = x_est(5:10) + dx(4:9);
lastZ = z;
lastP = P;
end