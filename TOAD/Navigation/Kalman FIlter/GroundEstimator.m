function [x_est, lastP] = GroundEstimator(x_est,constantsASTRA,z,dT)
%% M-EKF Implementation
% Remove bias from IMU
z(1:3) = z(1:3) - x_est(14:16);
z(4:6) = z(4:6) - x_est(11:13);

% Mag bias
z(7:9) = z(7:9) - x_est(17:19);
S = norm(z(7:9));
z(7:9) = z(7:9) / S;

% Extract quaternion
dx = zeros(18,1);
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
    P = 1 * eye(18); 
    P(1:3,1:3) = 0.5;
    P(10:12, 10:12) = 0.03;
    P(13:15, 13:15) = 0.25;
    P(16:18, 16:18) = 10;
    lastZ = zeros(15,1);
end

% State Transition Matrix
F = StateTransitionMat(z(1:3), z(4:6), R_b2i, 1);

% Propagate rest of state using IMU
x_est(8:10) = x_est(8:10) + (R_b2i * z(1:3) - [0; 0; constantsASTRA.g]) * dT;
x_est(5:7) = x_est(5:7) + x_est(8:10) * dT;

% Discrete STM
Phi = expm(F * dT);

% Magnetometer normalization factor
MagMatrix = (eye(3) - z(7:9) * z(7:9)') / S;

% Extract Matrices
Q = constantsASTRA.Q;
R = constantsASTRA.R;

% Matrix updates for the mags
R(4:6, 4:6) = 2e-3 * MagMatrix + 1e-5 * (z(7:9) * z(7:9)');

% Process Noise Covariance and a-priori propagation step
P = Phi * P * Phi' + Q;
P = (P + P') / 2;
RTK = 1;

%% IMU Update
if any(lastZ(1:9) ~= z(1:9))

    % Measurement matrix
    H = zeros(6,18);
    H(1:3, 1:3) = zetaCross(R_b2i' * [0; 0; constantsASTRA.g]);
    H(1:3, 13:15) = eye(3);
    H(4:6, 1:3) = zetaCross(R_b2i' * constantsASTRA.mag);
    H(4:6, 16:18) = MagMatrix;
        
    % Predicted measurements 
    z_hat = [R_b2i' * [0; 0; constantsASTRA.g];
             R_b2i' * constantsASTRA.mag];

    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);
    
    % Kalman Gain Weighting based on predicted acceleration
    ILH = (eye(18) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z([1:3 7:9]) - z_hat);
    dx = dx + L * residual;
end

%% GPS Update
if any(lastZ(10:15) ~= z(10:15))

    % Measurement matrix
    H = zeros(6,18);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);

    % Measurement Covariance Matrix
    gps_pos_covar = 1 * RTK + 10 * (1 - RTK);
    gps_vel_covar = gps_pos_covar * 1;
    R = diag([gps_pos_covar^2 * ones(3,1); gps_vel_covar^2 * ones(3,1)]);

    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);
    L(10:18, :) = 0;

    % Predicted measurements 
    z_hat = [x_est(5:7);
             x_est(8:10)];

    % Kalman Gain Weighting based on predicted acceleration
    ILH = (eye(18) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z(10:15) - z_hat) - H * dx;
    inn = L * residual;
    dx = dx + inn;
else
    %% Section purely for output on testing
    % % Measurement matrix
    % H = zeros(6,18);
    % H(1:3, 1:3) = zetaCross(R_b2i' * [0; 0; constantsASTRA.g]);
    % H(1:3, 13:15) = eye(3);
    % H(4:6, 1:3) = zetaCross(R_b2i' * constantsASTRA.mag);
    % H(4:6, 16:18) = MagMatrix;
    
    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);
    L(10:18, :) = 0;
end

% Output
lastP = P;

% Update full-state estimates
dq = [1; dx(1:3) / 2];
dq = dq / norm(dq);

q_nom = quatmultiply(q', dq');
q_nom = q_nom / norm(q_nom); 
x_est(1:4) = q_nom';
x_est(5:19) = x_est(5:19) + dx(4:18);
lastZ = z;
% lastP = P;
end