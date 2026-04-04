function F = StateTransitionMat(accel, gyro, R_b2i, GND)
% Remove angular rates from error-state. Make safety copies of all relevant
% files into Archive
if GND
    F = zeros(18);
    F(1:3,1:3) = -zetaCross(gyro);
    F(1:3,10:12) = -eye(3);
    F(4:6,7:9) = eye(3);
    F(7:9,1:3) = -R_b2i * zetaCross(accel);
    F(7:9,13:15) = -R_b2i;
else
    F = zeros(9);
    F(1:3,1:3) = -zetaCross(gyro);
    F(4:6,7:9) = eye(3);
    F(7:9,1:3) = -R_b2i * zetaCross(accel);
end
