function A = StateTransitionMat_IEKF( ...
    omegaHat, R_b2i, pHat, vHat, U, mass, constantsTOAD)
% StateTransitionMat_IEKF
%
% Continuous-time error dynamics for the 18-state TVC-aided IEKF.
%
% Error ordering:
%   1:3     attitude
%   4:6     position
%   7:9     velocity
%   10:12   gyro bias
%   13:15   accelerometer bias
%   16:18   magnetometer bias

A = zeros(18);

theta  = U(1);
phi    = U(2);
thrust = U(3);

% Exact TOAD thrust-vector model.
%
% This intentionally replaces the paper's small-angle TVC approximation
% with the force equation actually used by TOAD's plant.
T_B = thrust * [ ...
    cos(theta)*sin(phi);
   -sin(theta);
    cos(theta)*cos(phi)];

g_I = [0;
       0;
      -constantsTOAD.g];

% Attitude error dynamics.
A(1:3,1:3) = -zetaCross(omegaHat);

% Gyro-bias error corrupts attitude propagation.
A(1:3,10:12) = -eye(3);

% Position derivative is velocity.
A(4:6,7:9) = eye(3);

% Right-invariant gyro-bias / position coupling.
A(4:6,10:12) = -zetaCross(pHat);

% TVC coupling.
%
% This is the main addition from the TVC-aided estimator:
%
% delta(v_dot)_TVC =
%     -(1/m) R [T_B]x delta(theta)
%
% An attitude error therefore predicts a translational error whenever the
% vehicle is producing thrust.
A_TVC = ...
    -(1/mass) ...
    * R_b2i ...
    * zetaCross(T_B);

% Gravity / attitude coupling.
A_gravity = ...
    zetaCross(R_b2i' * g_I);

A(7:9,1:3) = ...
    A_TVC + A_gravity;

% Right-invariant gyro-bias / velocity coupling.
A(7:9,10:12) = ...
    -zetaCross(vHat);

% Accelerometer bias corrupts navigation-frame acceleration.
A(7:9,13:15) = ...
    -R_b2i;

% Bias deterministic dynamics are zero because the biases are modeled as
% random walks.

end