% test script cooked up by Claude not me
clear astra_mex;  % wipes the persistent C++ state so every run starts fresh
clc;

%% Settings
dT        = 0.001;   % 1 kHz loop (Q_gen and DigitalNF assume this)
t_ground  = 5.0;     % seconds on the pad (estimator settles)
t_flight  = 15.0;     % seconds of "flight"
g         = 9.8015;
mass      = 1.169;
noise_on  = true;    % set false for a perfectly clean run

N_g = round(t_ground / dT);
N_f = round(t_flight / dT);
N   = N_g + N_f;

%% Initialize the C++ code
astra_mex('init');

%% Constant target: upright, at the origin, not moving
X_trg = zeros(15,1);
X_trg(1) = 1;                   % quaternion [w x y z] = identity
U_ff  = [0; 0; mass*g; 0];      % [gimbal pitch; gimbal yaw; thrust N; roll accel], hover thrust

%% Sensor input template (vehicle sitting upright and still)
ci.GND_val = true;
ci.new_mag_packet = false;
ci.new_gps_packet = false;
ci.imu = struct('accel_x',0,'accel_y',0,'accel_z',g, ...
    'gyro_yaw',0,'gyro_pitch',0,'gyro_roll',0);
ci.mag = struct('mag_x',0.38535,'mag_y',0.03198,'mag_z',-0.9222);
ci.gps.pos = struct('north',0,'west',0,'up',0);
ci.gps.vel = struct('north',0,'west',0,'up',0);
for f = ["posCovNN","posCovNE","posCovND","posCovEE","posCovED","posCovDD", ...
        "velCovNN","velCovNE","velCovND","velCovEE","velCovED","velCovDD"]
    ci.gps.(f) = 1;
end

%% Logs
t      = (0:N-1)' * dT;
log_co = zeros(N,4);    % thrust, roll accel, gimbal pitch, gimbal yaw
log_q  = zeros(N,4);
log_p  = zeros(N,3);
log_v  = zeros(N,3);
log_gnd = false(N,1);

%% Run
tic;
for k = 1:N
    ci.GND_val = (k <= N_g);
    ci.new_mag_packet = mod(k,10) == 0;    % mag at 100 Hz
    ci.new_gps_packet = mod(k,100) == 0;   % GPS at 10 Hz

    s = struct('ax',0,'ay',0,'az',0,'gx',0,'gy',0,'gz',0);
    if noise_on
        s.ax = 0.05*randn; s.ay = 0.05*randn; s.az = 0.05*randn;
        s.gx = 0.002*randn; s.gy = 0.002*randn; s.gz = 0.002*randn;
    end
    ci.imu.accel_x = 0 + s.ax;
    ci.imu.accel_y = 0 + s.ay;
    ci.imu.accel_z = g + s.az;
    ci.imu.gyro_yaw   = s.gx;
    ci.imu.gyro_pitch = s.gy;
    ci.imu.gyro_roll  = s.gz;

    [co, in] = astra_mex('step', ci, X_trg, U_ff, dT, dT);

    log_co(k,:) = [co.thrust_N, co.roll_rad_sec_squared, co.gimbal_pitch_deg, co.gimbal_yaw_deg];
    log_q(k,:)  = [in.x_est.q_vec_w, in.x_est.q_vec_x, in.x_est.q_vec_y, in.x_est.q_vec_z];
    log_p(k,:)  = [in.x_est.est_pos_north, in.x_est.est_pos_west, in.x_est.est_pos_up];
    log_v(k,:)  = [in.x_est.est_vel_north, in.x_est.est_vel_west, in.x_est.est_vel_up];
    log_gnd(k)  = ci.GND_val;
end
fprintf('Ran %d steps in %.2f s\n', N, toc);

%% Sanity checks
fprintf('\n--- Sanity checks ---\n');
fprintf('Any NaN/Inf in outputs:      %d (want 0)\n', any(~isfinite([log_co(:); log_q(:); log_p(:); log_v(:)])));
fprintf('Quaternion norm range:       %.4f to %.4f (want ~1)\n', min(vecnorm(log_q,2,2)), max(vecnorm(log_q,2,2)));
fprintf('Ground thrust (want 0):      %.3f N\n', max(abs(log_co(log_gnd,1))));
fprintf('Final quaternion:            [%s]\n', num2str(log_q(end,:), '%.3f '));
fprintf('Final position estimate:     [%s] m\n', num2str(log_p(end,:), '%.3f '));
fprintf('Final thrust command:        %.2f N (hover ~ %.2f N)\n', log_co(end,1), mass*g);

%% Plots
figure('Name','ASTRA MEX test');
subplot(2,2,1); plot(t, log_q); grid on; xline(t_ground,'--k');
title('Estimated quaternion'); legend('w','x','y','z'); xlabel('s');
subplot(2,2,2); plot(t, log_p); grid on; xline(t_ground,'--k');
title('Estimated position (N, W, Up)'); legend('N','W','Up'); xlabel('s'); ylabel('m');
subplot(2,2,3); plot(t, log_co(:,1)); grid on; xline(t_ground,'--k');
title('Thrust command'); xlabel('s'); ylabel('N');
subplot(2,2,4); plot(t, log_co(:,3:4)); grid on; xline(t_ground,'--k');
title('Gimbal commands'); legend('pitch','yaw'); xlabel('s'); ylabel('deg');