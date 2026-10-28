
%Taken from example implementation
p2.est_cov = 1.0;
gyro_bias_cov = 1e-9;
accel_proc_cov = 1e-2;
accel_bias_cov = 1e-12;
mag_proc_cov = 1e-1;
mag_bias_cov = 1e-12; % * 3000;

gyro_cov = 5e-5;
accel_obs_cov = 3e-2;
mag_obs_cov = 2e-1;

%Init constants
dt = 0.005;

%body vectors, defined from flight computer notation
p2.f_global_body = [0 0 -9.81]; %Expected gravity vector from 0 (standing upright)
p2.m_global_body = [1 0 0]; %Expected normalized field vector from 0 (standing upright)

p2.obsv_cov_mat = eye(6);
p2.obsv_cov_mat(1:3,1:3) = accel_obs_cov * eye(3);
p2.obsv_cov_mat(4:6,4:6) = mag_obs_cov * eye(3);

p2.gyro_cov_mat = gyro_cov*eye(3);
p2.gyro_bias_cov_mat = gyro_bias_cov*eye(3);
p2.accel_cov_mat = accel_proc_cov*eye(3);
p2.accel_bias_cov_mat = accel_bias_cov*eye(3);
p2.mag_cov_mat = mag_proc_cov*eye(3);
p2.mag_bias_cov_mat = mag_bias_cov*eye(3);

p2.G = zeros(15);
p2.G(1:3,10:12) = -eye(3);
p2.G(7:9,4:6) = eye(3);

p2.Q = zeros(18);
p2.Q(1:3, 1:3) = p2.gyro_cov_mat*dt + p2.gyro_bias_cov_mat*(dt^3)/3.0;
p2.Q(1:3, 10:12) = -p2.gyro_bias_cov_mat*(dt^2)/2.0;
p2.Q(4:6, 4:6) = p2.accel_cov_mat*(dt^3)/3.0 + p2.accel_bias_cov_mat*(dt^5)/20.0;
p2.Q(4:6, 7:9) = p2.accel_cov_mat*(dt^2)/2.0 + p2.accel_bias_cov_mat*(dt^4)/8.0;
p2.Q(4:6, 13:15) = -p2.accel_bias_cov_mat*(dt^3)/6.0;
p2.Q(7:9, 4:6) = p2.accel_cov_mat*(dt^2)/2.0 + p2.accel_bias_cov_mat*(dt^4)/8.0;
p2.Q(7:9, 7:9) = p2.accel_cov_mat*dt + p2.accel_bias_cov_mat*(dt^3)/3.0;
p2.Q(7:9, 13:15) = -p2.accel_bias_cov_mat*(dt^2)/2.0;
p2.Q(10:12, 1:3) = -p2.gyro_bias_cov_mat*(dt^2)/2.0;
p2.Q(10:12, 10:12) = p2.gyro_bias_cov_mat*dt;
p2.Q(13:15, 4:6) = -p2.accel_bias_cov_mat*(dt^2)/2.0;
p2.Q(13:15, 7:9) = -p2.accel_bias_cov_mat*(dt^3)/6.0;
p2.Q(13:15, 13:15) = p2.accel_bias_cov_mat*dt;
p2.Q(16:18, 16:18) = p2.mag_bias_cov_mat*dt;

%% TVCMEKF2 Constants
dt2 = 0.002;

% Tunable scalars
t2.qAtt = 500;
t2.gyroBiasCov = 1e-9;
t2.velProcCov = 1e-2;
t2.accelBiasCov= 1e-2;

% Thrust-model uncertainty scaled by thrust
p2.tvc.sigDir = deg2rad(1);
p2.tvc.sigMag = 0.1;

% Error state: [theta position velocity gyro bias accel bias]
p2.Q2 = zeros(18);

p2.Q2(1:3,1:3) = t2.qAtt*(gyro_cov*dt2 + t2.gyroBiasCov*dt2^3/3)*eye(3);
p2.Q2(1:3,10:12) = -t2.gyroBiasCov*dt2^2/2*eye(3);
p2.Q2(10:12,1:3) = -t2.gyroBiasCov*dt2^2/2*eye(3);
p2.Q2(10:12,10:12) = t2.gyroBiasCov*dt2*eye(3);


p2.Q2(4:6,4:6) = t2.velProcCov*dt2^3/3*eye(3);
p2.Q2(4:6,7:9) = t2.velProcCov*dt2^2/2*eye(3);
p2.Q2(7:9,4:6) = t2.velProcCov*dt2^2/2*eye(3);
p2.Q2(7:9,7:9) = t2.velProcCov*dt2*eye(3);

% Accel bias no coupling to pos/vel only random walk
p2.Q2(13:15,13:15) = t2.accelBiasCov*dt2*eye(3);

% R for tvc filter -> accel row uses a_thrust
p2.R2 = eye(6);
p2.R2(1:3,1:3) = 3e-2*eye(3); %tune
p2.R2(4:6,4:6) = 1e-2*eye(3);