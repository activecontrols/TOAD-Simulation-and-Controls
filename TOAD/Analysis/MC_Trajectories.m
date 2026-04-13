%% Parallel Monte Carlo Trajectory Setup & Extraction
% --- Configuration ---
model_name = 'TOAD_Simulation';
num_sims = 100; % Reduced for trajectory logging
clear simIn out

% Nominal parameters (Ensure constantsTOAD is loaded in base workspace first)
GrommetIDX = 1;
J_nom = constantsTOAD.J;
G = GrommetSelect(GrommetIDX);
m_FC = 0.1;
K_nom = G.K;
B_nom = G.C / (2 * sqrt(K_nom * m_FC));

% Preallocate arrays for parameters
J_d_vals  = cell(1, num_sims);
TB_d_vals = cell(1, num_sims);
GyroNoisePower_vals = cell(1, num_sims);
G_RMAX_vals = cell(1, num_sims);
kGrom_vals = cell(1, num_sims);
bGrom_vals = cell(1, num_sims);
Kg2_vals = cell(1, num_sims);

% Preallocate arrays for the final SSE vectors
RMSE_Controls_all = zeros(9, num_sims); 
RMSE_Filter_all   = zeros(3, num_sims); 

disp(['Generating disturbances for ', num2str(num_sims), ' runs...']);

for i = 1:num_sims
    % 1. Moment of Inertia Disturbances (Delta J)
    dI_xx = (0.1 * J_nom(1,1)) * rand();
    dI_yy = (0.1 * J_nom(2,2)) * rand();
    dI_zz = (0.1 * J_nom(3,3)) * rand();
    dI_xy = 0; dI_xz = 0; dI_yz = 0;
    
    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];
               
    % 2. Lever Arm Disturbances (Delta Lever Arm)
    sigma_lever = [0.005; 0.005; 0.005]; 
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;

    % 3. Gyro Biases & Params
    GyroNoisePower_vals{i} = LogNormal(10^-7, 1);
    G_RMAX_vals{i} = (8 - 3) * rand() + 3;
    kGrom_vals{i} = K_nom * (1 + 0.150 * randn());
    bGrom_vals{i} = B_nom * (1 + 0.150 * randn());
    Kg2_vals{i} = (0.05-0.0005) * rand() + 0.0005;
end

%% Setup Simulation Inputs for Parallel Execution
disp('Configuring Parallel Simulation Inputs...');

% Properly preallocate the SimulationInput object array to suppress warnings
simIn(1:num_sims) = Simulink.SimulationInput(model_name);

for i = 1:num_sims
    % Securely attach this iteration's variables
    simIn(i) = simIn(i).setVariable('J_d', J_d_vals{i});
    simIn(i) = simIn(i).setVariable('TB_d', TB_d_vals{i});
    simIn(i) = simIn(i).setVariable('gyroNoisePower', GyroNoisePower_vals{i});
    simIn(i) = simIn(i).setVariable('kGrom', kGrom_vals{i});
    simIn(i) = simIn(i).setVariable('bGrom', bGrom_vals{i});
    simIn(i) = simIn(i).setVariable('Kg2', Kg2_vals{i});
    simIn(i) = simIn(i).setVariable('G_RMAX', G_RMAX_vals{i});

    % Enable state logging at 10Hz (SampleTime = 0.1)
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/state_log', 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/state_log', 'SampleTime', '0.1');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/target_pos_log', 'Commented', 'off');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/target_pos_log', 'SampleTime', '0.1');
    
    % Keep other large data blocks commented out to save memory
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/meas_log', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/Subsystem Reference2/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_state', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/Subsystem Reference2/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_P', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/inputCMD', 'Commented', 'on');
end

%% Execute Parallel Simulations
disp('Starting Parallel Monte Carlo Trajectory Simulations (parsim)...');
myCluster = parcluster('Processes');
delete(myCluster.Jobs);

out = parsim(simIn, 'ShowProgress', 'on', 'UseFastRestart', 'on');

%% Extract Data 
disp('Simulations Complete. Extracting and Interpolating Trajectories...');
t_sim = 75;
t_common = (0:0.1:t_sim)'; % Common time vector at 10Hz
pos_all = nan(num_sims, length(t_common), 3);
vel_all = nan(num_sims, length(t_common), 3);

for i = 1:num_sims
    if isempty(out(i).ErrorMessage)
        % Metric extraction
        RMSE_Controls_all(:, i) = sqrt(out(i).SSE_Controls(:) ./ t_sim); 
        RMSE_Filter_all(:, i)   = sqrt(out(i).SSE_Filter(:) ./ t_sim);
        
        % State extraction (Position is indices 5:7)
        ts_state = out(i).state_log;
        t_raw = ts_state.Time;
        data_raw = squeeze(ts_state.Data);
        if size(data_raw, 1) ~= length(t_raw)
            data_raw = data_raw';
        end
        pos_raw = data_raw(:, 5:7);
        
        % Interpolate onto common time grid for average calculation
        pos_all(i, :, 1) = interp1(t_raw, pos_raw(:,1), t_common, 'linear', 'extrap');
        pos_all(i, :, 2) = interp1(t_raw, pos_raw(:,2), t_common, 'linear', 'extrap');
        pos_all(i, :, 3) = interp1(t_raw, pos_raw(:,3), t_common, 'linear', 'extrap');

        vel_raw = data_raw(:, 8:10); 
        vel_all(i, :, 1) = interp1(t_raw, vel_raw(:,1), t_common, 'linear', 'extrap');
        vel_all(i, :, 2) = interp1(t_raw, vel_raw(:,2), t_common, 'linear', 'extrap');
        vel_all(i, :, 3) = interp1(t_raw, vel_raw(:,3), t_common, 'linear', 'extrap');
    else
        warning('Simulation %d failed: %s', i, out(i).ErrorMessage);
        RMSE_Controls_all(:, i) = NaN;
        RMSE_Filter_all(:, i)   = NaN;
    end
end

%% Analysis & Averages
disp('Analyzing Data...');
pos_avg = squeeze(mean(pos_all, 1, 'omitnan'));

% Calculate Symmetry Metrics
Lever_Radial = zeros(num_sims, 1); Lever_Axial = zeros(num_sims, 1);
J_Trans_Scale = zeros(num_sims, 1); J_Axial_Scale = zeros(num_sims, 1);

for i = 1:num_sims
    Lever_Radial(i) = norm(TB_d_vals{i}(1:2)); 
    Lever_Axial(i)  = abs(TB_d_vals{i}(3));    
    J_Trans_Scale(i) = norm([J_d_vals{i}(1,1), J_d_vals{i}(2,2)]); 
    J_Axial_Scale(i) = abs(J_d_vals{i}(3,3));           
end

disp('Data successfully calculated.');

%% Render Plots 
PlotMCTrajectories();

function samples = LogNormal(target_mode, sigma)
    mu_normal = log(target_mode) + (sigma^2);
    samples = exp(mu_normal + sigma * randn());
end