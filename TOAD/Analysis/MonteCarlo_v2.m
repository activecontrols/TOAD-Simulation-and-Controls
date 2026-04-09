%% Parallel Monte Carlo Setup & Disturbance Generation
% --- Configuration ---
model_name = 'TOAD_Simulation';
num_sims = 200;

% Nominal parameters (Ensure constantsTOAD is loaded in base workspace first)
J_nom = constantsTOAD.J;

% Preallocate arrays for parameters
J_d_vals  = cell(1, num_sims);
TB_d_vals = cell(1, num_sims);

% Preallocate arrays for the final SSE vectors to save memory
SSE_Controls_all = zeros(9, num_sims); 
SSE_Filter_all   = zeros(3, num_sims); % 3x1 for Roll, Pitch, Yaw Error

disp(['Generating disturbances for ', num2str(num_sims), ' runs...']);

for i = 1:num_sims
    % 1. Moment of Inertia Disturbances (Delta J)
    dI_xx = (0.10 * J_nom(1,1)) * randn();
    dI_yy = (0.10 * J_nom(2,2)) * randn();
    dI_zz = (0.10 * J_nom(3,3)) * randn();
    dI_xy = 5.0 * randn();
    dI_xz = 5.0 * randn();
    dI_yz = 5.0 * randn();
    
    J_d_vals{i} = [dI_xx, dI_xy, dI_xz;
                   dI_xy, dI_yy, dI_yz;
                   dI_xz, dI_yz, dI_zz];
               
    % 2. Lever Arm Disturbances (Delta Lever Arm)
    sigma_lever = [0.02; 0.02; 0.1]; 
    TB_d_vals{i} = randn(3, 1) .* sigma_lever;
end

%% Setup Simulation Inputs for Parallel Execution
disp('Configuring Parallel Simulation Inputs...');

% Create an array of SimulationInput objects
simIn(1:num_sims) = Simulink.SimulationInput(model_name);

for i = 1:num_sims
    % Securely attach this iteration's variables
    simIn(i) = simIn(i).setVariable('J_d', J_d_vals{i});
    simIn(i) = simIn(i).setVariable('TB_d', TB_d_vals{i});
    
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

out = parsim(simIn, 'ShowProgress', 'on');

%% Extract Data 
disp('Simulations Complete. Extracting Data...');

for i = 1:num_sims
    if isempty(out(i).ErrorMessage)
        
        % The timeseries .Data property holds the exact vector we need.
        % The (:) forces it into a strict column vector, absolutely 
        % guaranteeing MATLAB won't corrupt the dimensions.
        SSE_Controls_all(:, i) = out(i).SSE_Controls.Data(:); 
        SSE_Filter_all(:, i)   = out(i).SSE_Filter.Data(:);
        
    else
        warning('Simulation %d failed: %s', i, out(i).ErrorMessage);
        SSE_Controls_all(:, i) = NaN;
        SSE_Filter_all(:, i)   = NaN;
    end
end

%% Analysis
disp('Analyzing Data...');

J_mag           = zeros(num_sims, 1);
Lever_mag       = zeros(num_sims, 1);

metric_Ctrl_Att = zeros(num_sims, 1);
metric_Ctrl_Pos = zeros(num_sims, 1);
metric_Ctrl_Vel = zeros(num_sims, 1);
metric_Filter   = zeros(num_sims, 1);

for i = 1:num_sims
    % Input magnitudes
    J_mag(i)     = norm(J_d_vals{i}); 
    Lever_mag(i) = norm(TB_d_vals{i});
    
    % Split the 9x1 Controller SSE into its physical domains
    metric_Ctrl_Att(i) = norm(SSE_Controls_all(1:3, i)); 
    metric_Ctrl_Pos(i) = norm(SSE_Controls_all(4:6, i)); 
    metric_Ctrl_Vel(i) = norm(SSE_Controls_all(7:9, i)); 
    
    % Filter Metric
    metric_Filter(i) = norm(SSE_Filter_all(:, i));
end

disp('-----------------------------------');
disp('Data successfully calculated.');
disp('-----------------------------------');

%% Auto-Save for Large Datasets (Data Only)
save_data = 1;
if save_data
    disp('Saving workspace data...');
    
    % Use Present Working Directory (pwd) to bypass the AppData Temp bug
    % This relies on your Current Folder pane being set to the Analysis folder
    script_dir = pwd; 
    
    % Simple folder structure inside your project directory
    save_dir = fullfile(script_dir, 'Analysis', 'Monte Carlo Runs');
    if ~exist(save_dir, 'dir')
        mkdir(save_dir);
    end
    
    % Clean, simple naming convention: MC_500runs_YYYYMMDD_HHMM.mat
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