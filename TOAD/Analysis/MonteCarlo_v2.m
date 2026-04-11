%% Parallel Monte Carlo Setup & Disturbance Generation
% --- Configuration ---
model_name = 'TOAD_Simulation';
num_sims = 200;
clear simIn out

% Nominal parameters (Ensure constantsTOAD is loaded in base workspace first)
J_nom = constantsTOAD.J;

% Preallocate arrays for parameters
J_d_vals  = cell(1, num_sims);
TB_d_vals = cell(1, num_sims);
GyroNoisePower_vals = cell(1, num_sims);
G_RMAX_vals = call(1, num_sims);
kGrom_vals = cell(1, num_sims);
bGrom_vals = cell(1, num_sims);
Kg2_vals = cell(1, num_sims);

% Preallocate arrays for the final SSE vectors to save memory
RMSE_Controls_all = zeros(9, num_sims); 
RMSE_Filter_all   = zeros(3, num_sims); % 3x1 for Roll, Pitch, Yaw Error

disp(['Generating disturbances for ', num2str(num_sims), ' runs...']);

for i = 1:num_sims
    % 1. Moment of Inertia Disturbances (Delta J)
    dI_xx = (0.12 * J_nom(1,1)) * randn();
    dI_yy = (0.12 * J_nom(2,2)) * randn();
    dI_zz = (0.12 * J_nom(3,3)) * randn();
    dI_xy = 1 * randn();
    dI_xz = 1 * randn();
    dI_yz = 1 * randn();
    
    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];
               
    % 2. Lever Arm Disturbances (Delta Lever Arm)
    sigma_lever = [0.01; 0.01; 0.02]; 
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;

    % 3. Gyro Biases
    GyroNoisePower_vals{i} = LogNormal(10^-6, 1.2);

    G_RMAX_vals(i) = (20 - 4) * rand() + 4;
    
    kGrom_nom = Grom.k;
    kGrom = kGrom_nom * (1 + 0.05 * randn());
    bGrom_nom = Grom.c / 1 / sqrt(kGrom_nom * 0.1);
    bGrom = bGrom_nom * (1 + 0.05 * randn());

    Kg2_vals(i) = (0.08-0.005) * rand() + 0.005;
end

%% Setup Simulation Inputs for Parallel Execution
disp('Configuring Parallel Simulation Inputs...');

% Create an array of SimulationInput objects
simIn(1:num_sims) = Simulink.SimulationInput(model_name);

for i = 1:num_sims
    % Securely attach this iteration's variables
    simIn(i) = simIn(i).setVariable('J_d', J_d_vals{i});
    simIn(i) = simIn(i).setVariable('TB_d', TB_d_vals{i});
    simIn(i) = simIn(i).setVariable('gyroNoisePower', GyroNoisePower_vals{i});

    % Save RAM & Comment out other timeseries logs
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/state_log', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/meas_log', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/Subsystem Reference2/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_state', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/Subsystem Reference2/Multiplicative Extended Kalman Filter [M-EKF]/MEKF_P', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/inputCMD', 'Commented', 'on');
    simIn(i) = simIn(i).setBlockParameter('TOAD_Simulation/target_pos_log', 'Commented', 'on');
end

%% Execute Parallel Simulations
disp('Starting Parallel Monte Carlo Simulations (parsim)...');

myCluster = parcluster('Processes');
delete(myCluster.Jobs);

out = parsim(simIn, 'ShowProgress', 'on', 'UseFastRestart', 'on');

%% Extract Data 
disp('Simulations Complete. Extracting Data...');
t_sim = 60;
for i = 1:num_sims
    if isempty(out(i).ErrorMessage)
        
        % Calculate running RMSE from the SSE Simulink outputs
        RMSE_Controls_all(:, i) = sqrt(out(i).SSE_Controls(:) ./ t_sim); 
        RMSE_Filter_all(:, i)   = sqrt(out(i).SSE_Filter(:) ./ t_sim);
        
    else
        warning('Simulation %d failed: %s', i, out(i).ErrorMessage);
        RMSE_Controls_all(:, i) = NaN;
        RMSE_Filter_all(:, i)   = NaN;
    end
end

%% Analysis
disp('Analyzing Data...');

% --- New Isolated Symmetry Metrics ---
Lever_Radial  = zeros(num_sims, 1);
Lever_Axial   = zeros(num_sims, 1);
J_Trans_Scale = zeros(num_sims, 1);
J_Axial_Scale = zeros(num_sims, 1);
J_Wobble_Coup = zeros(num_sims, 1);
J_Trans_Coup  = zeros(num_sims, 1);

metric_Ctrl_Att = zeros(num_sims, 1);
metric_Ctrl_Pos = zeros(num_sims, 1);
metric_Ctrl_Vel = zeros(num_sims, 1);
metric_Filter   = zeros(num_sims, 1);

for i = 1:num_sims
    % 1. Lever Arm: Axial vs. Radial
    Lever_Radial(i) = norm(TB_d_vals{i}(1:2)); % sqrt(dX^2 + dY^2)
    Lever_Axial(i)  = abs(TB_d_vals{i}(3));    % |dZ|
    
    % 2. Inertia Tensor: Dynamic Symmetry Grouping
    dI_xx = J_d_vals{i}(1,1);
    dI_yy = J_d_vals{i}(2,2);
    dI_zz = J_d_vals{i}(3,3);
    dI_xy = J_d_vals{i}(1,2);
    dI_xz = J_d_vals{i}(1,3);
    dI_yz = J_d_vals{i}(2,3);
    
    J_Trans_Scale(i) = norm([dI_xx, dI_yy]); % Pitch/Yaw inertia error
    J_Axial_Scale(i) = abs(dI_zz);           % Roll inertia error
    J_Wobble_Coup(i) = norm([dI_xz, dI_yz]); % Roll -> Pitch/Yaw coupling
    J_Trans_Coup(i)  = abs(dI_xy);           % Pitch <-> Yaw coupling
    
    % Split the 9x1 Controller SSE into its physical domains
    metric_Ctrl_Att(i) = norm(RMSE_Controls_all(1:3, i)); 
    metric_Ctrl_Pos(i) = norm(RMSE_Controls_all(4:6, i)); 
    metric_Ctrl_Vel(i) = norm(RMSE_Controls_all(7:9, i)); 
    
    % Filter Metric
    metric_Filter(i) = norm(RMSE_Filter_all(:, i));
end

disp('-----------------------------------');
disp('Data successfully calculated.');
disp('-----------------------------------');

%% Auto-Save for Large Datasets (Data Only)
save_data = 1;
if save_data
    disp('Saving workspace data...');
    
    script_dir = pwd; 
    save_dir = fullfile(script_dir, 'Analysis', 'Monte Carlo Runs');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    timestamp = datestr(now, 'yyyymmdd_HHMM');
    mat_filename = fullfile(save_dir, sprintf('MC_%druns_%s.mat', num_sims, timestamp));
    
    % Save Workspace (using -v7.3 to support files larger than 2GB)
    save(mat_filename, '-v7.3');
    
    disp(['Data saved successfully to: ', mat_filename]);
else
    disp('Save bypassed (save_data flag set to false).');
end

%% Render Plots using the active workspace
PlotMonteCarlo();

function samples = LogNormal(target_mode, sigma)
    % Calculate the underlying normal mean (mu) based on the target mode
    mu_normal = log(target_mode) + (sigma^2);
    
    % Generate using standard normal (randn), scale, shift, and exponentiate
    samples = exp(mu_normal + sigma * randn());
end