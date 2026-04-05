function z = SensorRates(Accel, Gyros, Mag, GPS, constantsASTRA, THRUST, t)

% Simulates magnetometer, accelerometer, and gyros given the true state and
% a covariance vector. No biases for now.

% Update rates
IMU_Rate = 1000;     %Hz
GPS_Rate = 15;      %Hz
z = zeros(15,1);

persistent lastGPS lastIMU lastZ phase

if isempty(lastGPS)
    lastGPS = 0;  % Update last GPS timestamp
    lastIMU = 0;  % Update last IMU timestamp
    lastZ = [0; 0; 9.81; zeros(3,1); constantsASTRA.mag; zeros(6,1)];
    phase = 0;
end

%% Fake IMU measurements (add simulated vibrations based on thrust data)
% Thrust track
Track = [1.3525    42.6278];
if t - lastIMU > 1 / IMU_Rate
    % Thrust vibration generator
    w = 2 * pi * Track(1) * THRUST + Track(2) + 5 * randn();
    phase = phase + w * (t - lastIMU);
    MaxAmplAccel = 10;
    MaxAmplGyro = 0.1;
    z(1:3) = Accel + MaxAmplAccel * sin(phase);
    z(4:6) = Gyros + MaxAmplGyro * sin(phase);
    z(7:9) = Mag;
    lastIMU = t;
else
    z(1:9) = lastZ(1:9);
end

% Fake GPS measurements
z(10:15) = GPS;
lastZ = z;