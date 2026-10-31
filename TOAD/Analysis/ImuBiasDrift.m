% Brownian motion with drift, per axis:  dX = mu dt + sigma dZ
%   X_t = X_0 + mu*t + sigma*Z_t,   Z_t ~ N(0, t), with X_0 = 0
% (the constant turn-on bias is already accelBias/gyroBias)
% mu is the expected drift RATE: E[X_t] = X_0 + mu*t

function [AccelDrift, GyroDrift] = ImuBiasDrift(dt, IMUDriftParams)

persistent Za Zg t
if isempty(Za)
    rng(IMUDriftParams(9), 'twister');   % per-run seed, set by MonteCarlo_v5
    Za = zeros(3,1);  Zg = zeros(3,1);  t = 0;
end

mu_a = IMUDriftParams(1:3)';  sig_a = IMUDriftParams(4);
mu_g = IMUDriftParams(5:7)';  sig_g = IMUDriftParams(8);

t  = t + dt;
Za = Za + sqrt(dt) * randn(3,1);         % Brownian increments ~ N(0, dt)
Zg = Zg + sqrt(dt) * randn(3,1);

AccelDrift = mu_a * t + sig_a * Za;      % m/s^2
GyroDrift  = mu_g * t + sig_g * Zg;      % rad/s
end