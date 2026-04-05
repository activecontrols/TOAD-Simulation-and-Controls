function GPS = GPS_Sim(PosVel, Rate, R_b2i, t)

% Integrate rate limiting directly into GPS measurements instead of through
% SensorRates function. Keep Sensor Rates for IMU only.

GPS_Rate = 20;      %Hz
persistent error_pos lastGPS lastTime
if isempty(error_pos)
    error_pos = 0;
    lastGPS = zeros(6,1);
    lastTime = 0;  % Update last GPS timestamp
end

dT = t - lastTime;

if t - lastTime > 1 / GPS_Rate
    % Correlated noise measurements (GMP)
    RTK = 1;
    gps_pos_covar = 0.01 * RTK + 0.5 * (1 - RTK);
    gps_vel_covar = 0.005 * RTK + 0.05 * (1 - RTK);
    decay = 50;
    beta = exp(-dT / decay);
    error_pos = beta * error_pos + sqrt(1 - beta^2) * gps_pos_covar * randn;
    
    % Fake GPS Measurements (derivate velocity from GPS_Pos)
    GPS = zeros(6,1);
    GPS(1:3) = PosVel(1:3) + error_pos;
    GPS(4:6) = (PosVel(4:6) + gps_vel_covar * randn) + R_b2i * cross(Rate, [0 0 0.31]');
    lastGPS = GPS;
    lastTime = t;
else
    GPS = lastGPS;
end
