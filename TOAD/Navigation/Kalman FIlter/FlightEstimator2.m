function [x_est, lastP] = FlightEstimator2(x_est,constantsASTRA,z,dT,P0)
%% M-EKF Implementation
% Remove bias from IMU
z(1:3) = z(1:3) - x_est(14:16);
z(4:6) = z(4:6) - x_est(11:13);
z(7:9) = z(7:9) - x_est(17:19);
S = norm(z(7:9));
z(7:9) = z(7:9) / S;

% Extract quaternion
dx = zeros(12,1);
q = x_est(1:4);
if sqrt(sum(z(4:6).^2))<1e-5 || dT == 0
    x_est(1:4) = q;
else
    x_est(1:4) = (HamiltonianProd(q) * [cos(sqrt(sum(z(4:6).^2)).*dT./2); (z(4:6)/sqrt(sum(z(4:6).^2))).*sin(sqrt(sum(z(4:6).^2)).*dT./2)]);
end
q = x_est(1:4)/norm(x_est(1:4));

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

% Magnetometer normalization factor
MagMatrix = (eye(3) - z(7:9) * z(7:9)') / S;

% Extract Matrices
Q = constantsASTRA.Q(1:12, 1:12);
R = constantsASTRA.R(4:6, 4:6);

% Matrix updates for the mags
R(1:3, 1:3) = 5e-1 * MagMatrix + 6e-9 * (z(7:9) * z(7:9)');

% Process Noise Covariance and a-priori propagation step
P = Phi * P * Phi' + Q;
P = (P + P') / 2;

%% IMU Update
if any(lastZ(1:9) ~= z(1:9))

    % Measurement matrix
    H = zeros(3,12);
    H(1:3, 1:3) = zetaCross(R_b2i' * constantsASTRA.mag);
        
    % Predicted measurements 
    z_hat = R_b2i' * constantsASTRA.mag;

    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);
    
    % Kalman Gain Weighting based on predicted acceleration
    ILH = (eye(12) - L * H);
    P = ILH * P * ILH' + L * R * L';
    P = (P + P') / 2;
    residual = (z(7:9) - z_hat);
    dx = dx + L * residual;
end

%% GPS Update
if any(lastZ(10:15) ~= z(10:15))

    % Measurement matrix
    H = zeros(6,12);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);

    % Measurement Covariance Matrix
    GyroCovar = eye(3) * 5e-2;
    R = zeros(6);
    R(1:3, 1:3) = 0.1 * eye(3);
    R(4:6, 4:6) = 0.3 * eye(3) + R_b2i * zetaCross(rGPS) * GyroCovar * (R_b2i * zetaCross(rGPS))';

    % A priori covariance and Kalman gain
    L = P * H' / (H * P * H' + R);

    % Predicted measurements 
    z_hat = [x_est(5:7);
             x_est(8:10)];

    % Kalman Gain Weighting based on predicted acceleration
    ILH = (eye(12) - L * H);
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
x_est(5:13) = x_est(5:13) + dx(4:12);
lastZ = z;
% lastP = P;
end