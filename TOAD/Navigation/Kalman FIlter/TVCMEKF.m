function [x_est, lastP] = TVCMEKF(x_est, constantsASTRA, z, dT, U, mass, a_ref)
%% M-EKF with thrust-model propagation, force-error state d_f, and
%  reference-acceleration update.
%
% Inputs
%   x_est    22x1: [q(1:4); p(5:7); v(8:10); b_g(11:13); b_a(14:16); b_m(17:19); d_f(20:22)]
%            d_f = body-frame specific-force error [m/s^2] (new state)
%   z        15x1: [accel(1:3); gyro(4:6); mag(7:9); gps_pos(10:12); gps_vel(13:15)]
%   U        [theta; phi; T; roll] applied gimbal angles, thrust
%   mass     current vehicle mass 
%   a_ref    3x1 planned INERTIAL acceleration (e.g. a_ff from the trajectory)
%
% Error state (21): [dtheta(1:3) dp(4:6) dv(7:9) dbg(10:12) dba(13:15) dbm(16:18) ddf(19:21)]
% dtheta is a BODY-frame error (q_true = q_est (x) dq), so
%   R_true = R_est*(I + [dtheta]x)   and   R_true'*w = R'*w + [R'*w]x*dtheta.


%%  Tuning 
sigDf    = 0.05;    % d_f random-walk sigma (tune)
sigTrack = 0.5;     % [m/s^2] trajectory-tracking error treated as accel noise (tune)
RTK      = 1;

g    = constantsASTRA.g;
gvec = [0; 0; -g];
nx   = 21;

zRaw = z;           % raw copy for new-sample detection

%%  Thrust model (body frame) 
theta = U(1);  phi = U(2);  T = U(3);
T_b = T * [cos(theta)*sin(phi);
          -sin(theta);
           cos(theta)*cos(phi)];
a_thrust = T_b / mass;
d_f = x_est(20:22);

%%  Bias removal 
z(1:3) = z(1:3) - x_est(14:16);
z(4:6) = z(4:6) - x_est(11:13);
z(7:9) = z(7:9) - x_est(17:19);
S = norm(z(7:9));
z(7:9) = z(7:9) / S;

%%  Quaternion propagation (gyro) 
dx = zeros(nx,1);
q  = x_est(1:4);
w  = z(4:6);
wn = norm(w);
if wn < 1e-5 || dT == 0
    x_est(1:4) = q;
else
    x_est(1:4) = HamiltonianProd(q) * [cos(wn*dT/2); (w/wn)*sin(wn*dT/2)];
end
q = x_est(1:4) / norm(x_est(1:4));
R_b2i = quatRot(q)';

% GPS velocity lever-arm correction
rGPS = [0 0 0.31]';
z(13:15) = z(13:15) - R_b2i * cross(z(4:6), rGPS);

%%  Covariance init / reset 
persistent P lastZ
if isempty(P)
    P = initP();
    lastZ = zeros(15,1);
end
if dT <= 0
    P = initP();
    lastZ = zeros(15,1);
    dT = 0;
end

%%  Velocity / position propagation 
% Thrust model + force error
f_prop = a_thrust + d_f;
x_est(8:10) = x_est(8:10) + (R_b2i * f_prop + gvec) * dT;
x_est(5:7)  = x_est(5:7)  + x_est(8:10) * dT;

%%  Error-state dynamics 
F = zeros(nx);
F(1:3, 1:3)   = -zetaCross(w);
F(1:3, 10:12) = -eye(3);
F(4:6, 7:9)   = eye(3);
F(7:9, 1:3)   = -R_b2i * zetaCross(f_prop);
F(7:9, 19:21) = R_b2i;
Phi = expm(F * dT);

%%  Noise matrices 
Q = zeros(nx);
Q(1:18, 1:18) = constantsASTRA.Q;
Q(19:21, 19:21) = (sigDf * max(norm(a_thrust), 1))^2 * dT * eye(3);

Rc = constantsASTRA.R;
MagMatrix = (eye(3) - z(7:9)*z(7:9)') / S;
Rmag = 2e-1 * MagMatrix + 1e-6 * (z(7:9)*z(7:9)');

P = Phi * P * Phi' + Q;
P = (P + P') / 2;

%%  IMU updates 
if any(lastZ(1:9) ~= zRaw(1:9))

    % Reference-acceleration update (accel vs planned inertial accel, estimated attitude)
    f_pred = R_b2i' * (a_ref(:) + [0; 0; g]);

    H = zeros(3, nx);
    H(:, 1:3)   = zetaCross(f_pred);       % pitch/yaw observbility fix
    H(:, 13:15) = eye(3);
    Racc = Rc(1:3,1:3) + sigTrack^2 * eye(3);
    [P, dx] = seqUpdate(P, dx, H, Racc, z(1:3) - f_pred);

    % Thrust-model update: observes d_f + b_a
    H = zeros(3, nx);
    H(:, 13:15) = eye(3);
    H(:, 19:21) = eye(3);
    [P, dx] = seqUpdate(P, dx, H, Rc(1:3,1:3), z(1:3) - (a_thrust + d_f));

    % Magnetometer
    H = zeros(3, nx);
    H(:, 1:3)   = zetaCross(R_b2i' * constantsASTRA.mag);
    H(:, 16:18) = MagMatrix;
    [P, dx] = seqUpdate(P, dx, H, Rmag, z(7:9) - R_b2i' * constantsASTRA.mag);
end

%%  GPS update 
if any(lastZ(10:15) ~= zRaw(10:15))
    H = zeros(6, nx);
    H(1:3, 4:6) = eye(3);
    H(4:6, 7:9) = eye(3);

    gps_pos_covar = 1 * RTK + 10 * (1 - RTK);
    gps_vel_covar = gps_pos_covar * 1;
    Rg = diag([gps_pos_covar^2 * ones(3,1); gps_vel_covar^2 * ones(3,1)]);

    z_hat = [x_est(5:7); x_est(8:10)];
    [P, dx] = seqUpdate(P, dx, H, Rg, z(10:15) - z_hat);
end

lastP = P;

%%  Inject error state 
dq = [1; dx(1:3)/2];
dq = dq / norm(dq);
q_nom = quatmultiply(q', dq');
q_nom = q_nom / norm(q_nom);
x_est(1:4)  = q_nom';
x_est(5:22) = x_est(5:22) + dx(4:21);   % pos, vel, b_g, b_a, b_m, d_f
lastZ = zRaw;
end

%% Local functions 
function [P, dx] = seqUpdate(P, dx, H, Rm, r0)
    L = P * H' / (H * P * H' + Rm);
    L(16:18, :) = 0;                     % mag-bias gain frozen 
    ILH = eye(size(P)) - L * H;
    P = ILH * P * ILH' + L * Rm * L';
    P = (P + P') / 2;
    dx = dx + L * (r0 - H * dx);
end

function P = initP()
    P = eye(21);
    P(1:3,   1:3)   = 0.05 * eye(3);   % Attitude
    P(10:12, 10:12) = 0.01 * eye(3);   % Gyro bias
    P(13:15, 13:15) = 0.01 * eye(3);   % Accel bias
    P(16:18, 16:18) = 0.01 * eye(3);   % Mag bias
    P(19:21, 19:21) = 0.01 * eye(3);   % Force error d_f
end