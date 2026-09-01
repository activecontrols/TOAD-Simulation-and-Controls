%% Parallel Monte Carlo Trajectory Setup & Extraction (V4)
% V4 focus:
%   - Wind disturbances
%   - Center-of-mass / lever-arm disturbances
%   - Moment-of-inertia disturbances
%
% Legacy vibration / estimator / actuator uncertainties are retained below
% as commented-out code for easy restoration, but are not generated or
% passed into the simulations in this version.
%
% The new controller uses a full-state trajectory target directly; there is
% no Waypoint/TrajectoryBuilder dependency in this script.

%% Configuration
model_name = 'TOAD_Simulation';
num_sims = 50;
seed = 1788146097; % int64(seconds(datetime('now','Timezone','UTC')-datetime('1970-01-01','Timezone','UTC')));

clear simIn out

%% Nominal parameters
GrommetIDX = 2;
J_nom = constantsTOAD.J;
G = GrommetSelect(GrommetIDX);
m_FC = 0.1;
K_nom = G.K;
B_nom = G.C / (2 * sqrt(K_nom * m_FC));

%% Preallocate Monte Carlo disturbance arrays
J_d_vals  = cell(1, num_sims);
TB_d_vals = cell(1, num_sims);

Wind_Gain_vals  = cell(1, num_sims);
Wind_Covar_vals = cell(1, num_sims);

%% Legacy / optional uncertainty arrays
% Leave these commented for data-transfer efficiency. Uncomment the
% corresponding generation and setVariable calls if needed.
%
% GyroNoisePower_vals = cell(1, num_sims);
% G_RMAX_vals         = cell(1, num_sims);
% kGrom_vals          = cell(1, num_sims);
% bGrom_vals          = cell(1, num_sims);
% Kg2_vals            = cell(1, num_sims);
% PSD_Low_vals        = cell(1, num_sims);
% PSD_High_vals       = cell(1, num_sims);
%
% % Optional explicit CoM disturbance, if the model exposes a separate
% % CoM_d workspace variable. Keep disabled unless that variable exists.
% % CoM_d_vals = cell(1, num_sims);

%% Preallocate arrays for simulation outputs
RMSE_Controls_all = zeros(12, num_sims);
RMSE_Filter_all   = zeros(3, num_sims);
RMSE_Wind_all     = zeros(1, num_sims);
MaxLESODist_all   = zeros(6, num_sims);

%% Preallocate interpolated state histories
t_sim = 60;
t_common = (0:0.1:t_sim)';

pos_all  = nan(num_sims, length(t_common), 3);
vel_all  = nan(num_sims, length(t_common), 3);
quat_all = nan(num_sims, length(t_common), 4);

% The target is now a full-state trajectory. Its width is determined after
% parsim from the actual target logger rather than assuming a legacy
% waypoint/position-only target format.
target_state_all = [];
disp(seed)
disp(['Generating focused disturbances for ', num2str(num_sims), ' runs with seed ', seed]);

%% Generate Monte Carlo disturbances
rng(seed)
for i = 1:num_sims

    %% 1. Moment of Inertia Disturbances
    dI_xx = (0.1 * J_nom(1,1)) * rand();
    dI_yy = (0.1 * J_nom(2,2)) * rand();
    dI_zz = (0.1 * J_nom(3,3)) * rand();

    % Keep products diagonal to match the V3 disturbance model.
    dI_xy = 0;
    dI_xz = 0;
    dI_yz = 0;

    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];

    %% 2. Center-of-Mass / Lever-Arm Disturbance
    % V3 used TB_d as the 3-axis lever-arm / off-center-CG disturbance.
    % Preserve that model-facing variable and keep its radial/axial
    % components available for sensitivity analysis.
    sigma_lever = [0.01; 0.01; 0.01];
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;

    %% 3. Wind Disturbances
    % V3 wind gain distribution and covariance range are retained.
    a = 0.35;
    b = 2;
    Wind_Gain_vals{i} = a * (-log(rand(1,1))).^(1/b);
    Wind_Covar_vals{i} = 8 * rand();

    %% Legacy uncertainty generation (disabled)
    % GyroNoisePower_vals{i} = LogNormal(10^-7, 1);
    % G_RMAX_vals{i} = (8 - 3) * rand() + 3;
    % kGrom_vals{i} = K_nom * (0.5 * (1 + 2 * rand()));
    % bGrom_vals{i} = B_nom * (0.5 * (1 + 2 * rand()));
    % Kg2_vals{i} = (0.1 - 0.002) * rand() + 0.002;
    % PSD_Low_vals{i} = 20 + (100 - 20) * rand();
    % PSD_High_vals{i} = 500 + (1000 - 500) * rand();
    %
    % % Explicit CoM disturbance, if supported by the model:
    % % CoM_d_vals{i} = randn(3,1) .* sigma_com;
end

%% Setup Simulation Inputs for Parallel Execution
disp('Configuring Parallel Simulation Inputs...');
simIn(1:num_sims) = Simulink.SimulationInput(model_name);

for i = 1:num_sims
    simIn(i) = simIn(i).setVariable('J_d', J_d_vals{i});
    simIn(i) = simIn(i).setVariable('TB_d', TB_d_vals{i});
    simIn(i) = simIn(i).setVariable('Wind_Gain', Wind_Gain_vals{i});
    simIn(i) = simIn(i).setVariable('Wind_Covar', Wind_Covar_vals{i});

    %% Optional legacy variables (disabled)
    % simIn(i) = simIn(i).setVariable('gyroNoisePower', GyroNoisePower_vals{i});
    % simIn(i) = simIn(i).setVariable('G_RMAX', G_RMAX_vals{i});
    % simIn(i) = simIn(i).setVariable('kGrom', kGrom_vals{i});
    % simIn(i) = simIn(i).setVariable('bGrom', bGrom_vals{i});
    % simIn(i) = simIn(i).setVariable('Kg2', Kg2_vals{i});
    % simIn(i) = simIn(i).setVariable('lowEnd', PSD_Low_vals{i});
    % simIn(i) = simIn(i).setVariable('highEnd', PSD_High_vals{i});
    %
    % % Explicit CoM disturbance, if supported:
    % simIn(i) = simIn(i).setVariable('CoM_d', CoM_d_vals{i});

    %% Logging: retain only the state and full target trajectory
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/state_log', 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/state_log', 'SampleTime', '0.5');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/Controller & LESOs/MaxLESODist', 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/target_pos_log', 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/target_pos_log', 'SampleTime', '0.5');

    % Disable unused logs to reduce data movement / RAM use.
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/meas_log', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/State Estimator & Filters/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_state', ...
        'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/State Estimator & Filters/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_P', ...
        'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter( ...
        'TOAD_Simulation/inputCMD', 'Commented', 'on');
end

%% Execute Parallel Simulations
disp('Starting Parallel Monte Carlo Trajectory Simulations (parsim)...');
out = sim(simIn, 'ShowProgress', 'on', 'UseFastRestart', 'on');

%% Extract Metrics and Interpolate Trajectories
disp('Simulations complete. Extracting and interpolating trajectories...');

for i = 1:num_sims
    if isempty(out(i).ErrorMessage)
        % Existing scalar/vector SSE metrics converted to RMSE in the same
        % manner as V3. Wind is now an RMSE metric rather than max-gust.
        RMSE_Controls_all(:, i) = sqrt(out(i).SSE_Controls(:) ./ t_sim);
        RMSE_Filter_all(:, i)   = sqrt(out(i).SSE_Filter(:) ./ t_sim);

        % SSE_Wind is treated as an integrated squared-error metric, just
        % like the other SSE outputs. If it contains multiple wind
        % channels, combine them into one RMS wind-error magnitude.
        wind_sse = out(i).SSE_Wind(:);
        RMSE_Wind_all(i) = sqrt(sum(wind_sse) / t_sim);
       MaxLESODist_all(:, i) = out(i).MaxLESODist.Data(:);

        % Actual state
        ts_state = out(i).state_log;
        t_raw = ts_state.Time;
        data_raw = squeeze(ts_state.Data);

        if size(data_raw, 1) ~= length(t_raw)
            data_raw = data_raw';
        end

        if size(data_raw,2) < 10
            error('Simulation %d: state_log contains only %d columns; at least 10 are required.', ...
                i, size(data_raw,2));
        end

        quat_raw = data_raw(:, 1:4);
        pos_raw  = data_raw(:, 5:7);
        vel_raw  = data_raw(:, 8:10);

        for dim = 1:3
            pos_all(i,:,dim) = interp1(t_raw, pos_raw(:,dim), ...
                t_common, 'linear', 'extrap');
            vel_all(i,:,dim) = interp1(t_raw, vel_raw(:,dim), ...
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

        % Allocate after observing the actual target-state width.
        if isempty(target_state_all)
            n_target_states = size(target_raw, 2);
            target_state_all = nan(num_sims, length(t_common), n_target_states);
        end

        if size(target_raw, 2) ~= size(target_state_all, 3)
            error(['Simulation %d: target trajectory state dimension changed ' ...
                   'from %d to %d.'], ...
                i, size(target_state_all,3), size(target_raw,2));
        end

        for dim = 1:size(target_raw,2)
            target_state_all(i,:,dim) = interp1( ...
                t_target, target_raw(:,dim), t_common, 'linear', 'extrap');
        end
    else
        warning('Simulation %d failed: %s', i, out(i).ErrorMessage);
        RMSE_Controls_all(:, i) = NaN;
        RMSE_Filter_all(:, i)   = NaN;
        RMSE_Wind_all(:, i)     = NaN;
    end
end

% Derived disturbance metrics for sensitivity analysis
disp('Calculating disturbance magnitudes...');

Lever_Radial = nan(num_sims,1);
Lever_Axial  = nan(num_sims,1);

J_Trans_Scale = nan(num_sims,1);
J_Axial_Scale = nan(num_sims,1);
J_Wobble_Coup = nan(num_sims,1);
J_Trans_Coup  = nan(num_sims,1);

Wind_Gain_all  = cell2mat(Wind_Gain_vals);
Wind_Covar_all = cell2mat(Wind_Covar_vals);

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

%%Target-relative trajectory metrics
% No waypoint assumptions are used. Compare each simulated trajectory to
% the full target trajectory directly. The state ordering used here follows
% the logged vehicle state convention from V3:
%   [q(1:4), position(1:3), velocity(1:3), ...]
if isempty(target_state_all) || size(target_state_all,3) < 10
    error(['The full-state target trajectory was not logged with the expected ' ...
           'position/velocity entries (states 5:10).']);
end

target_pos_all = target_state_all(:,:,5:7);
target_vel_all = target_state_all(:,:,8:10);

pos_error_all = pos_all - target_pos_all;
vel_error_all = vel_all - target_vel_all;

pos_error_rmse = squeeze(sqrt(mean(pos_error_all.^2, 2, 'omitnan')));
vel_error_rmse = squeeze(sqrt(mean(vel_error_all.^2, 2, 'omitnan')));

% Aggregate position / velocity errors for compact plotting.
Pos_RMSE_total = sqrt(sum(pos_error_rmse.^2, 2));
Vel_RMSE_total = sqrt(sum(vel_error_rmse.^2, 2));

disp('Data successfully calculated.');

% Auto-Save
save_data = 1;
if save_data
    disp('Saving workspace data...');
    save_dir = fullfile(pwd, 'Analysis', 'Monte Carlo Runs');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end

    mat_filename = fullfile(save_dir, ...
        sprintf('MC4_%druns_%s.mat', num_sims, datestr(now, 'yyyymmdd_HHMM')));

    save(mat_filename, '-v7.3');
    disp(['Data saved successfully to: ', mat_filename]);
end

% Render Plots
PlotMC_v4();

%% Local helper
function samples = LogNormal(target_mode, sigma)
    mu_normal = log(target_mode) + sigma^2;
    samples = exp(mu_normal + sigma * randn());
end
