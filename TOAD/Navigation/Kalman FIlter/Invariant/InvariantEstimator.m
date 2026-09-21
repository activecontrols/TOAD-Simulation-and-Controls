function [x_est, lastP] = GroundIEKF(x_est, constantsTOAD, z, U, Mass, dT)
% GroundIEKF
%
% TVC-aided iterated invariant EKF adapted to the existing TOAD estimator.
%
% Nominal state:
%   x_est(1:4)    quaternion, scalar first
%   x_est(5:7)    position, NWU
%   x_est(8:10)   velocity, NWU
%   x_est(11:13)  gyro bias
%   x_est(14:16)  accelerometer bias
%   x_est(17:19)  magnetometer bias
%
% Error state:
%   1:3     attitude
%   4:6     position
%   7:9     velocity
%   10:12   gyro bias
%   13:15   accelerometer bias
%   16:18   magnetometer bias
%
% Measurements:
%   z(1:3)    accelerometer
%   z(4:6)    gyroscope
%   z(7:9)    magnetometer
%   z(10:12)  GNSS position
%   z(13:15)  GNSS velocity
%
% Controller:
%   U(1)  theta gimbal
%   U(2)  phi gimbal
%   U(3)  thrust
%   U(4)  roll torque
%
% The same controller command sent to the vehicle is used by the IEKF as
% a known deterministic input.
%
% d_f and d_tau are treated as stochastic disturbances, not filter states.

persistent P lastZ

if isempty(P)

    P = eye(18);

    P(1:3,1:3)     = 0.05 * eye(3);
    P(4:6,4:6)     = 1.00 * eye(3);
    P(7:9,7:9)     = 1.00 * eye(3);
    P(10:12,10:12) = 0.01 * eye(3);
    P(13:15,13:15) = 0.01 * eye(3);
    P(16:18,16:18) = 0.01 * eye(3);

    lastZ = zeros(15,1);

end

if dT <= 0
    lastP = P;
    lastZ = z;
    return
end

%% Extract nominal state

q     = x_est(1:4);
pHat  = x_est(5:7);
vHat  = x_est(8:10);
bgHat = x_est(11:13);
baHat = x_est(14:16);
bmHat = x_est(17:19);

q = q / norm(q);

% Existing quatRot returns inertial-to-body.
R_b2i = quatRot(q)';

%% Extract measurements

accelMeas = z(1:3);
gyroMeas  = z(4:6);
magMeas   = z(7:9);

accelHat = accelMeas - baHat;
omegaHat = gyroMeas  - bgHat;

%% Determine vehicle mass

if numel(Mass) >= 2
    mass = constantsTOAD.m_dry + Mass(1) + Mass(2);
else
    mass = Mass;
end

mass = max(mass, constantsTOAD.m_dry);

%% Nominal attitude propagation

omegaNorm = norm(omegaHat);

if omegaNorm > 1e-10

    angle = omegaNorm * dT;

    dq = [cos(angle/2);
          (omegaHat/omegaNorm) * sin(angle/2)];

    q = HamiltonianProd(q) * dq;
    q = q / norm(q);

end

R_b2i = quatRot(q)';

%% Nominal translational propagation

% NWU gravity.
g_I = [0;
       0;
      -constantsTOAD.g];

% Accelerometer is specific force.
%
% TVC thrust is NOT added here because the accelerometer already measures
% the vehicle's thrust-induced specific acceleration.
a_I = R_b2i * accelHat + g_I;

pHat = pHat + vHat*dT + 0.5*a_I*dT^2;
vHat = vHat + a_I*dT;

%% Construct continuous invariant error dynamics

A = StateTransitionMat_IEKF( ...
    omegaHat, ...
    R_b2i, ...
    pHat, ...
    vHat, ...
    U, ...
    mass, ...
    constantsTOAD);

Phi = expm(A*dT);

%% Process covariance

Qd = BuildIEKFProcessNoise( ...
    constantsTOAD, ...
    mass, ...
    dT);

Ppred = Phi * P * Phi' + Qd;
Ppred = 0.5 * (Ppred + Ppred');

%% Save predicted state
%
% Every measurement iteration is referenced back to this prediction.
% The iterations are NOT separate Kalman measurement updates.

ChiPred = [R_b2i,      vHat, pHat;
           zeros(1,3), 1,    0;
           zeros(1,3), 0,    1];

bgPred = bgHat;
baPred = baHat;
bmPred = bmHat;

%% Determine which measurements are new

newMag = any(lastZ(7:9) ~= z(7:9));
newGPS = any(lastZ(10:15) ~= z(10:15));

if newMag || newGPS

    %% Initialize iterated state at prediction

    ChiIter = ChiPred;

    bgIter = bgPred;
    baIter = baPred;
    bmIter = bmPred;

    maxIter = 5;

    correctionTolerance = 1e-6;

    Kfinal = [];
    Hfinal = [];
    Rfinal = [];

    for iter = 1:maxIter

        %% Extract current iteration state

        Rj = ChiIter(1:3,1:3);
        vj = ChiIter(1:3,4);
        pj = ChiIter(1:3,5);

        %% Determine displacement of current iterate from prediction
        %
        % eta_j = Log(Chi_j * Chi_pred^-1)
        %
        % This is essential. Each iteration solves for the total
        % correction relative to the original predicted state rather than
        % treating every iteration as another independent measurement.

        ChiRelative = ChiIter / ChiPred;

        etaGroup = LogSE23(ChiRelative);

        eta18 = zeros(18,1);

        % LogSE23 ordering:
        %   attitude
        %   velocity
        %   position
        %
        % TOAD error ordering:
        %   attitude
        %   position
        %   velocity

        eta18(1:3) = etaGroup(1:3);
        eta18(4:6) = etaGroup(7:9);
        eta18(7:9) = etaGroup(4:6);

        eta18(10:12) = bgIter - bgPred;
        eta18(13:15) = baIter - baPred;
        eta18(16:18) = bmIter - bmPred;

        %% Build measurement model at current iterate

        Hblocks = {};
        residualBlocks = {};
        Rblocks = {};

        %% Magnetometer

        if newMag

            magRaw = magMeas - bmIter;
            magMagnitude = norm(magRaw);

            if magMagnitude > 1e-10

                magUnit = magRaw / magMagnitude;

                MagMatrix = ...
                    (eye(3) - magUnit*magUnit') / magMagnitude;

                magPred = Rj' * constantsTOAD.mag;
                magPred = magPred / norm(magPred);

                Hmag = zeros(3,18);

                Hmag(:,1:3) = zetaCross(magPred);

                % Magnetometer bias is a TOAD-specific extension to the
                % state used by the paper.
                Hmag(:,16:18) = MagMatrix;

                RmagRaw = constantsTOAD.R(4:6,4:6);

                Rmag = ...
                    MagMatrix * RmagRaw * MagMatrix' ...
                    + 1e-9*eye(3);

                rMag = magUnit - magPred;

                Hblocks{end+1} = Hmag;
                residualBlocks{end+1} = rMag;
                Rblocks{end+1} = Rmag;

            end

        end

        %% GNSS

        if newGPS

            Hgps = zeros(6,18);

            % Right-invariant GNSS Jacobian.
            Hgps(1:3,1:3) = ...
                -zetaCross(Rj' * pj);

            Hgps(1:3,4:6) = eye(3);

            Hgps(4:6,1:3) = ...
                -zetaCross(Rj' * vj);

            Hgps(4:6,7:9) = eye(3);

            zGPS = [z(10:12);
                    z(13:15)];

            zHatGPS = [pj;
                       vj];

            rGPS = zGPS - zHatGPS;

            % Preserve the original estimator's current RTK assumption.
            RTK = 1;

            gps_pos_covar = ...
                1*RTK + 10*(1-RTK);

            gps_vel_covar = gps_pos_covar;

            Rgps = diag([ ...
                gps_pos_covar^2 * ones(3,1);
                gps_vel_covar^2 * ones(3,1)]);

            Hblocks{end+1} = Hgps;
            residualBlocks{end+1} = rGPS;
            Rblocks{end+1} = Rgps;

        end

        %% Stack available measurements

        if isempty(Hblocks)
            break
        end

        H = vertcat(Hblocks{:});
        residual = vertcat(residualBlocks{:});

        Rmeas = blkdiag(Rblocks{:});

        %% Iterated Kalman gain
        %
        % Ppred remains fixed during all iterations.
        %
        % Updating P inside this loop would incorrectly count the same
        % physical measurement multiple times.

        S = H * Ppred * H' + Rmeas;

        K = (Ppred * H') / S;

        %% Compute new total correction
        %
        % The H*eta term accounts for the fact that the nonlinear
        % measurement has been relinearized away from the original
        % predicted state.
        %
        % dx is therefore the new estimate of the TOTAL correction from
        % ChiPred, not an incremental correction from ChiIter.

        dx = eta18 + ...
             K * (residual - H*eta18);

        %% Check convergence

        if norm(dx - eta18) < correctionTolerance

            Kfinal = K;
            Hfinal = H;
            Rfinal = Rmeas;

            break

        end

        %% Construct next iterate relative to prediction

        dtheta = dx(1:3);
        dp     = dx(4:6);
        dv     = dx(7:9);

        dChi = ExpSE23( ...
            dtheta, ...
            dv, ...
            dp);

        % IMPORTANT:
        %
        % Use ChiPred here, NOT ChiIter.
        %
        % dx represents the total correction relative to the original
        % prediction.
        ChiIter = dChi * ChiPred;

        % Euclidean states are also reconstructed from the prediction,
        % rather than incrementally accumulated.
        bgIter = bgPred + dx(10:12);
        baIter = baPred + dx(13:15);
        bmIter = bmPred + dx(16:18);

        Kfinal = K;
        Hfinal = H;
        Rfinal = Rmeas;

    end

    %% Use final iterated state

    R_b2i = ChiIter(1:3,1:3);
    vHat  = ChiIter(1:3,4);
    pHat  = ChiIter(1:3,5);

    bgHat = bgIter;
    baHat = baIter;
    bmHat = bmIter;

    %% Project rotation back onto SO(3)

    [Ur,~,Vr] = svd(R_b2i);

    R_b2i = Ur*Vr';

    if det(R_b2i) < 0
        Ur(:,3) = -Ur(:,3);
        R_b2i = Ur*Vr';
    end

    %% Covariance update
    %
    % This happens ONCE after all iterations have converged/completed.

    if ~isempty(Kfinal)

        I18 = eye(18);

        IKH = I18 - Kfinal*Hfinal;

        P = ...
            IKH*Ppred*IKH' ...
            + Kfinal*Rfinal*Kfinal';

        P = 0.5*(P + P');

    else

        P = Ppred;

    end

else

    % No new external observation.
    P = Ppred;

end

%% Convert corrected attitude back to quaternion

qNew = DCM_Quat_Conversion(R_b2i');

qNew = qNew(:);
qNew = qNew / norm(qNew);

% Prevent equivalent q/-q representation from causing apparent jumps.
if dot(qNew, q) < 0
    qNew = -qNew;
end

%% Repack nominal state

x_est(1:4)   = qNew;
x_est(5:7)   = pHat;
x_est(8:10)  = vHat;
x_est(11:13) = bgHat;
x_est(14:16) = baHat;
x_est(17:19) = bmHat;

lastP = P;
lastZ = z;

end