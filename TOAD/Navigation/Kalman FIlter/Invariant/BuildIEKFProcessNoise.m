function Qd = BuildIEKFProcessNoise(constantsTOAD, mass, dT)
% BuildIEKFProcessNoise
%
% Builds the discrete process-noise covariance for the 18-state filter.
%
% The existing TOAD Q matrix already includes discretized gyro noise,
% accelerometer noise, gyro-bias random walk, accelerometer-bias random
% walk, magnetometer-bias random walk, and the appropriate position /
% velocity cross-covariance terms.
%
% Rather than discard that tested structure, this function preserves it
% and adds stochastic force and torque disturbances.
%
% d_f and d_tau are NOT filter states. Their uncertainty is represented
% only through process covariance.

% constantsTOAD.Q was constructed using dt = 0.005.
%
% If the estimator continues running at exactly that rate, using it
% directly is appropriate as a first implementation.
%
% If dT varies significantly, this matrix should eventually be rebuilt
% from continuous noise spectral densities or discretized using Van Loan.
Qd = constantsTOAD.Q;

% Stochastic unmodeled force disturbance.
%
% sigma_df has units that should correspond to force-noise intensity.
% This is a tuning parameter because the paper does not provide a TOAD-
% specific value.
if isfield(constantsTOAD,'sigma_df')
    sigma_df = constantsTOAD.sigma_df;
else
    sigma_df = 5.0;
end

Q_df = sigma_df^2 * eye(3);

% Force disturbance produces acceleration:
%
%     delta(v_dot) = d_f / m
%
% For white force noise, first-order covariance injection into velocity
% is sufficient for the initial implementation.
Gdf = zeros(18,3);
Gdf(7:9,:) = eye(3)/mass;

Qd_df = Gdf * Q_df * Gdf' * dT;

% Stochastic unmodeled torque disturbance.
%
% There is an unavoidable approximation here because angular velocity is
% not part of the 18-state estimate
%
% Physically:
%
%     d_tau -> omega_dot -> omega -> attitude
%
% A fully rigorous continuous-time noise injection would therefore
% require angular-rate error states. Since the existing TOAD filter uses
% gyro measurements directly instead, torque uncertainty is mapped into
% attitude over one sample using:
%
%     delta(theta) ~= 0.5 * J^-1 * d_tau * dt^2
%
% This preserves the existing 18-state architecture at the cost of an
% approximate torque-disturbance model.
if isfield(constantsTOAD,'sigma_dtau')
    sigma_dtau = constantsTOAD.sigma_dtau;
else
    sigma_dtau = 1.0;
end

Q_tau = sigma_dtau^2 * eye(3);

Btau = zeros(18,3);
Btau(1:3,:) = 0.5 * (constantsTOAD.J \ eye(3)) * dT^2;

Qd_tau = Btau * Q_tau * Btau';

Qd = Qd + Qd_df + Qd_tau;

% Numerical symmetrization prevents roundoff from slowly making P
% asymmetric.
Qd = 0.5*(Qd + Qd');

end