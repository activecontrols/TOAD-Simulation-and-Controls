%% Parallel Monte Carlo Trajectory Setup & Extraction (V5)
% V5 additions (vehicle-agnostic; adapted for TOAD_Simulation_CPP):
%   - Physical mean wind in m/s (speed + direction) with Dryden turbulence intensity
%   - Symmetric inertia dispersion (J_d)
%   - Center of mass / lever arm dispersion (TB_d)
%   - IMU constant biases (accelBias, gyroBias), dispersed per run
%   - Magnetometer soft iron (magDistMatrix)
%   - Pass/fail per run (position and attitude error after pad hold)
%   - SSE_Wind extracted as timeseries RMSE
%   - Direct constants6DoF parameter ingestion
%   - Retains only verified, physically connected disturbances

%% Configuration
model_name = 'TOAD_Simulation_CPP';
num_sims = 250;

% Set to true for a new random batch every run, or false to lock a specific seed
use_random_master_seed = true; 

if use_random_master_seed
    rng('shuffle', 'twister');
    rng_state = rng;
    seed = rng_state.Seed;
else
    seed = uint32(1788146097);
    rng(seed, 'twister');
end

% Fast Restart skips recompiling between runs.
use_fast_restart = true;
clear simIn out

%% Nominal parameters
% Ensure constants6DoF is loaded from LoadTOADSim
J_nom = constants6DoF.J;
GrommetIDX = 2;
G = GrommetSelect(GrommetIDX);
m_FC = 0.1;
K_nom = G.K;
B_nom = G.C / (2 * sqrt(K_nom * m_FC));

%% V5 disturbance settings (all dispersion sizes in one place)
% Inertia: symmetric +/- fraction of nominal per principal axis
J_frac = 0.10;

% Mean wind, at the 6 m reference height used by the Dryden block
wind_sweep         = false;   % true: uniform 0..wind_max (envelope search)
wind_max           = 12;      % m/s, upper end of the envelope sweep
wind_weibull_scale = 4.0;     % m/s, operational campaign
wind_weibull_shape = 2.0;     % Rayleigh-like; matches measured surface wind
dryden_W6          = 15;      % m/s, "Wind speed at 6 m" set in the Dryden block mask
turb_mult_range    = [0.5, 2.0]; % multiplier on spec turbulence intensity

% IMU constant biases, perturbed around nominal values
accelBias_nom   = accelBias;              % m/s^2
gyroBias_nom    = gyroBias;               % rad/s
sigma_accelBias = 0.02;                   % m/s^2 (~2 mg) 1-sigma per axis
sigma_gyroBias  = deg2rad(0.1);           % rad/s 1-sigma per axis

% Magnetometer soft iron, perturbed around nominal values
magDist_nom     = magDistMatrix;
sigma_softIron  = 0.01;                   % per element of the (symmetric) matrix

% Pass/fail thresholds (after the pad hold)
pad_hold_time   = 15;      % s, PadHoldTime in TrajectoryManager
pos_err_fail    = 1.0;     % m, matches max_x_trans position limit
att_err_fail    = 15;      % deg

%% Preallocate Monte Carlo disturbance arrays
J_d_vals  = cell(1, num_sims);
TB_d_vals = cell(1, num_sims);

Wind_Gain_vals  = cell(1, num_sims);
Wind_Covar_vals = cell(1, num_sims);

Wind_Mean_vals     = cell(1, num_sims);   % m/s
Wind_Dir_vals      = cell(1, num_sims);   % rad
windMerid_vals     = cell(1, num_sims);
windZonal_vals     = cell(1, num_sims);
accelBias_vals     = cell(1, num_sims);
gyroBias_vals      = cell(1, num_sims);
magDistMatrix_vals = cell(1, num_sims);

% Seeds table for reproducibility (9 physically active channels)
num_seed_cols = 9;
rng(double(seed), 'twister');
noiseSeeds = randi(double(intmax('uint32')) - 1, num_sims, num_seed_cols, 'uint32');

%% Preallocate arrays for simulation outputs
RMSE_Filter_all = zeros(3, num_sims);
RMSE_Wind_all   = zeros(1, num_sims);

%% Preallocate interpolated state histories
t_sim = 60;
t_common = (0:0.1:t_sim)';

pos_all   = nan(num_sims, length(t_common), 3);
vel_all   = nan(num_sims, length(t_common), 3);
quat_all  = nan(num_sims, length(t_common), 4);
rates_all = nan(num_sims, length(t_common), 3);

target_state_all = [];

fprintf('Master Seed: %u\n', seed);
fprintf('Generating focused disturbances for %d runs with master seed %u...\n', num_sims, seed);

%% Generate Monte Carlo disturbances
for i = 1:num_sims

    %% 1. Moment of Inertia Disturbances (symmetric)
    rng(double(noiseSeeds(i, 1)), 'twister');
    dI_xx = (J_frac * J_nom(1,1)) * (2 * rand() - 1);
    dI_yy = (J_frac * J_nom(2,2)) * (2 * rand() - 1);
    dI_zz = (J_frac * J_nom(3,3)) * (2 * rand() - 1);

    dI_xy = 0;
    dI_xz = 0;
    dI_yz = 0;

    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];

    %% 2. Center-of-Mass / Lever-Arm Disturbance
    rng(double(noiseSeeds(i, 2)), 'twister');
    sigma_lever = [0.01; 0.01; 0.01];
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;

    %% 3. Wind Disturbances
    % 3a. Mean wind speed (m/s) and direction
    rng(double(noiseSeeds(i, 5)), 'twister');
    if wind_sweep
        Wind_Mean_vals{i} = wind_max * rand();
    else
        Wind_Mean_vals{i} = wind_weibull_scale * (-log(rand())).^(1/wind_weibull_shape);
    end
    rng(double(noiseSeeds(i, 6)), 'twister');
    Wind_Dir_vals{i} = 2 * pi * rand();
    windMerid_vals{i} = Wind_Mean_vals{i} * cos(Wind_Dir_vals{i});
    windZonal_vals{i} = Wind_Mean_vals{i} * sin(Wind_Dir_vals{i});

    % 3b. Turbulence scaling
    rng(double(noiseSeeds(i, 3)), 'twister');
    turb_mult = turb_mult_range(1) + diff(turb_mult_range) * rand();
    Wind_Gain_vals{i} = turb_mult * Wind_Mean_vals{i} / dryden_W6;

    rng(double(noiseSeeds(i, 4)), 'twister');
    Wind_Covar_vals{i} = 8 * rand();

    %% 4. IMU constant biases
    rng(double(noiseSeeds(i, 7)), 'twister');
    accelBias_vals{i} = accelBias_nom(:) + sigma_accelBias * randn(3, 1);
    rng(double(noiseSeeds(i, 8)), 'twister');
    gyroBias_vals{i}  = gyroBias_nom(:)  + sigma_gyroBias  * randn(3, 1);

    %% 5. Magnetometer soft iron (symmetric)
    rng(double(noiseSeeds(i, 9)), 'twister');
    E = sigma_softIron * randn(3);
    magDistMatrix_vals{i} = magDist_nom + (E + E') / 2;
end

%% Setup Simulation Inputs for Parallel Execution
disp('Configuring Parallel Simulation Inputs...');
simIn(1:num_sims) = Simulink.SimulationInput(model_name);

for i = 1:num_sims
    simIn(i) = simIn(i).setVariable('J_d', J_d_vals{i});
    simIn(i) = simIn(i).setVariable('TB_d', TB_d_vals{i});
    simIn(i) = simIn(i).setVariable('Wind_Gain', Wind_Gain_vals{i});
    simIn(i) = simIn(i).setVariable('Wind_Covar', Wind_Covar_vals{i});

    %% V5 physically connected disturbances
    simIn(i) = simIn(i).setVariable('windMerid', windMerid_vals{i});
    simIn(i) = simIn(i).setVariable('windZonal', windZonal_vals{i});
    simIn(i) = simIn(i).setVariable('magDistMatrix', magDistMatrix_vals{i});

    %% Logging configuration for TOAD_Simulation_CPP
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/state_log'], 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/state_log'], 'SampleTime', '0.1');
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/target_pos_log'], 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/target_pos_log'], 'SampleTime', '0.1');

    % Disable unused high-rate logs to minimize parsim worker transfer overhead
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/meas_log'], 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/inputCMD'], 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter( ...
        [model_name '/input_log'], 'Commented', 'on');
end

%% Base simulation variables
for i = 1:num_sims
    simIn(i) = simIn(i).setVariable('IMU_Rate', IMU_Rate);
    simIn(i) = simIn(i).setVariable('accelBias', accelBias_vals{i});
    simIn(i) = simIn(i).setVariable('gyroBias', gyroBias_vals{i});
    simIn(i) = simIn(i).setVariable('gyroNoisePower', gyroNoisePower);
    simIn(i) = simIn(i).setVariable('dt_SIM', dt_SIM);
    simIn(i) = simIn(i).setVariable('dt', dt);
end

%% Execute Parallel Simulations
disp('Starting Parallel Monte Carlo Trajectory Simulations (parsim)...');
if use_fast_restart, fast_restart = 'on'; else, fast_restart = 'off'; end
out = parsim(simIn, ...
    'ShowProgress', 'on', ...
    'UseFastRestart', fast_restart, ...
    'TransferBaseWorkspaceVariables', 'on');

ok = arrayfun(@(s) isempty(s.ErrorMessage), out);

if ~any(ok)
    error('MonteCarlo:AllRunsFailed', ...
        'No simulations succeeded. First error:\n%s', ...
        out(1).ErrorMessage);
end

%% Extract Metrics and Interpolate Trajectories
disp('Simulations complete. Extracting and interpolating trajectories...');

for i = 1:num_sims
    if isempty(out(i).ErrorMessage)
        % Filter Attitude RMSE from SSE_Filter
        RMSE_Filter_all(:, i) = sqrt(out(i).SSE_Filter(:) ./ t_sim);

        % Wind metric: SSE_Wind as a timeseries of wind velocity vector
        w_data = squeeze(out(i).SSE_Wind(:));
        RMSE_Wind_all(i) = sqrt(mean(sum(w_data.^2, 2)));

        % Actual state extraction
        ts_state = out(i).state_log;
        t_raw = ts_state.Time;
        data_raw = squeeze(ts_state.Data);

        if size(data_raw, 1) ~= length(t_raw)
            data_raw = data_raw';
        end

        if size(data_raw,2) < 13
            error('Simulation %d: state_log contains only %d columns; at least 13 are required for [q, pos, vel, omega].', ...
                i, size(data_raw,2));
        end

        quat_raw = data_raw(:, 1:4);
        pos_raw  = data_raw(:, 5:7);
        vel_raw  = data_raw(:, 8:10);
        rates_raw = data_raw(:, 11:13);

        for dim = 1:3
            pos_all(i,:,dim) = interp1(t_raw, pos_raw(:,dim), ...
                t_common, 'linear', 'extrap');
            vel_all(i,:,dim) = interp1(t_raw, vel_raw(:,dim), ...
                t_common, 'linear', 'extrap');
            rates_all(i,:,dim) = interp1(t_raw, rates_raw(:,dim), ...
                t_common, 'linear', 'extrap');
        end

        for dim = 1:4
            quat_all(i,:,dim) = interp1(t_raw, quat_raw(:,dim), ...
                t_common, 'linear', 'extrap');
        end

        % Full-state target trajectory
        ts_target = out(i).target_pos_log;
        t_target = ts_target.Time;
        target_raw = squeeze(ts_target.Data);

        if size(target_raw, 1) ~= length(t_target)
            target_raw = target_raw';
        end

        % Allocate target state history
        if isempty(target_state_all)
            n_target_states = size(target_raw, 2);
            target_state_all = nan(num_sims, length(t_common), n_target_states);
        end

        if size(target_raw, 2) ~= size(target_state_all, 3)
            error('Simulation %d: target trajectory state dimension changed from %d to %d.', ...
                i, size(target_state_all,3), size(target_raw,2));
        end

        for dim = 1:size(target_raw,2)
            target_state_all(i,:,dim) = interp1(t_target, target_raw(:,dim), t_common, 'linear', 'extrap'); %#ok<SAGROW>
        end
    else
        warning('Simulation %d failed: %s', i, out(i).ErrorMessage);
        RMSE_Filter_all(:, i) = NaN;
        RMSE_Wind_all(:, i)   = NaN;
    end
end

%% Derived disturbance metrics for sensitivity analysis
disp('Calculating disturbance magnitudes...');

Lever_Radial = nan(num_sims,1);
Lever_Axial  = nan(num_sims,1);

J_Trans_Scale = nan(num_sims,1);
J_Axial_Scale = nan(num_sims,1);
J_Wobble_Coup = nan(num_sims,1);
J_Trans_Coup  = nan(num_sims,1);

Wind_Gain_all  = cell2mat(Wind_Gain_vals);
Wind_Covar_all = cell2mat(Wind_Covar_vals);

Wind_Mean_all      = cell2mat(Wind_Mean_vals);
Wind_Dir_all       = cell2mat(Wind_Dir_vals);
AccelBias_Mag      = cellfun(@norm, accelBias_vals);
GyroBias_Mag       = cellfun(@norm, gyroBias_vals);
SoftIron_Dev       = cellfun(@(M) norm(M - magDist_nom, 'fro'), magDistMatrix_vals);

for i = 1:num_sims
    Lever_Radial(i) = norm(TB_d_vals{i}(1:2));
    Lever_Axial(i)  = abs(TB_d_vals{i}(3));

    dI_xx = J_d_vals{i}(1,1);
    dI_yy = J_d_vals{i}(2,2);
    dI_zz = J_d_vals{i}(3,3);
    dI_xy = J_d_vals{i}(1,2);
    dI_xz = J_d_vals{i}(1,3);
    dI_yz = J_d_vals{i}(2,3);

    J_Trans_Scale(i) = norm([dI_xx, dI_yy]);
    J_Axial_Scale(i) = abs(dI_zz);
    J_Wobble_Coup(i) = norm([dI_xz, dI_yz]);
    J_Trans_Coup(i)  = abs(dI_xy);
end

%% Target-relative trajectory metrics
if isempty(target_state_all) || size(target_state_all,3) < 7
    error('The full-state target trajectory was not logged with the expected position entries.');
end

target_pos_all = target_state_all(:,:,5:7);
target_vel_all = target_state_all(:,:,8:10);

pos_error_all = pos_all - target_pos_all;
vel_error_all = vel_all - target_vel_all;

pos_error_rmse = squeeze(sqrt(mean(pos_error_all.^2, 2, 'omitnan')));
vel_error_rmse = squeeze(sqrt(mean(vel_error_all.^2, 2, 'omitnan')));

Pos_RMSE_total = sqrt(sum(pos_error_rmse.^2, 2));
Vel_RMSE_total = sqrt(sum(vel_error_rmse.^2, 2));

%% V5 pass/fail per run
% Errors measured against target trajectory after pad hold
eval_idx = t_common >= pad_hold_time;

pos_err_norm = sqrt(sum(pos_error_all.^2, 3));                 % num_sims x T
Pos_Err_Max  = max(pos_err_norm(:, eval_idx), [], 2);          % m

q_sim = quat_all ./ sqrt(sum(quat_all.^2, 3));
q_trg = target_state_all(:,:,1:4);
q_trg = q_trg ./ sqrt(sum(q_trg.^2, 3));
q_dot = abs(sum(q_sim .* q_trg, 3));
q_dot(q_dot > 1) = 1;                                          % keeps NaN as NaN
att_err_deg = 2 * acosd(q_dot);                                % num_sims x T
Att_Err_Max = max(att_err_deg(:, eval_idx), [], 2);            % deg

Fail_Error    = ~ok(:);
Fail_NoData   = isnan(Pos_Err_Max) | isnan(Att_Err_Max);
Fail_Position = Pos_Err_Max > pos_err_fail;
Fail_Attitude = Att_Err_Max > att_err_fail;
MC_Fail = Fail_Error | Fail_NoData | Fail_Position | Fail_Attitude;

fprintf('\nV5 pass/fail: %d of %d runs failed (%.1f%%)\n', ...
    nnz(MC_Fail), num_sims, 100 * nnz(MC_Fail) / num_sims);
fprintf('   sim error      %d\n', nnz(Fail_Error));
fprintf('   no state data  %d\n', nnz(Fail_NoData & ~Fail_Error));
fprintf('   position > %.1f m   %d\n', pos_err_fail, nnz(Fail_Position));
fprintf('   attitude > %.0f deg  %d\n', att_err_fail, nnz(Fail_Attitude));

disp('Data successfully calculated.');

%% Auto-Save
save_data = 1;
if save_data
    disp('Saving workspace data...');
    save_dir = fullfile(pwd, 'Analysis', 'Monte Carlo Runs');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    time_str = char(datetime('now', 'Format', 'yyyyMMdd_HHmm'));
    mat_filename = fullfile(save_dir, sprintf('MC5_%druns_%s.mat', num_sims, time_str));

    save(mat_filename, '-v7.3');
    disp(['Data saved successfully to: ', mat_filename]);
end

%% Render Plots
PlotMC_v5();
