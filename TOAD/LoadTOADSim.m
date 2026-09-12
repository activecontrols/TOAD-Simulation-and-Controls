%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file loads in all the constants and parameters for the Simulink into
% workspace. Please always run this file before running a full-scale
% simulation if you've made any changes to trajectory, controls, filtering,
% or others.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clean reset of Simulink models, cached buses, and workspace
% if bdIsLoaded('TOAD_Simulation')
%     bdclose('TOAD_Simulation');
% end
% bdclose('all');
% clear classes;
% clear functions;
% clear slBus* Vehicle_Bus;

%% Select vehicle (1 for TOAD, 0 for ASTRA)
Vehicle = 0; % 1 for TOAD, 0 for ASTRA (default ASTRA)
constants6DoF = LoadTOADParams(Vehicle);
FlightDynamicsGen(constants6DoF);
x0 = [1; zeros(12,1); constants6DoF.OxMass; constants6DoF.FuMass];
u0 = [0; 0; constants6DoF.g * constants6DoF.m_wet; 0];

%% Guidance, Navigation, Control parameters

%% Kalman Filter params
    MEKF_Constants;
    % Kalman Filter & Control Parameters
    constants6DoF.Q = p2.Q;
    constants6DoF.R = p2.obsv_cov_mat;
    constants6DoF.BSigma = 5e-2;
    constants6DoF.BBias = 1e-8;
    
    % Magnetometer base measurement
    constants6DoF.mag = [0.385202; 0.030609; -0.922324];

%% Controller params
    % Outer Loop
    if Vehicle == 1
        % TOAD Tuning
        max_x_trans = 1.4 * ones(1,6);
    else
        % ASTRAv2 Tuning
        max_x_trans = 1.5 * ones(1,6);
    end
    constants6DoF.Q_trans = diag(1 ./ max_x_trans.^2);
    max_a_trans = 1.2; 
    constants6DoF.R_trans = eye(3) .* (1 / max_a_trans^2);
    constants6DoF.OmegaThr = 2.2;
    
    % Inner Loop
    if Vehicle == 1
        % TOAD Tuning
        max_x_rot = [0.12, 0.12, 0.12, 0.22, 0.22, 0.21];
        constants6DoF.R_rot = diag([35, 35, 1/4^2]);
    else
        % ASTRAv2 Tuning
        max_x_rot = [0.15, 0.15, 0.12, 0.5, 0.5, 0.7];
        constants6DoF.R_rot = diag([60, 60, 1/0.1^2]);
    end

    constants6DoF.Q_rot = diag(1 ./ max_x_rot.^2);
    constants6DoF.OmegaAtt = 3.2;

%% Trajectory Params
    % Pick a trajectory filename (e.g. "TOAD_Backflip_v001",
    % "ASTRA_Circle_v001")
    filename = "ASTRA_Backflip_v001";
    
    % Resolve trajectory CSV file path
    traj_dir = fullfile(pwd, 'Guidance', 'Trajectories');
    [~, name_stem, ext] = fileparts(filename);
    if isempty(ext)
        traj_file = fullfile(traj_dir, name_stem + ".csv");
        if ~exist(traj_file, 'file')
            traj_file = fullfile(traj_dir, filename);
        end
    else
        traj_file = fullfile(traj_dir, filename) + '.csv';
    end
    
    % Check existance
    if ~exist(traj_file, 'file')
        warning('LoadTOADSim:TrajectoryNotFound', ...
            ['Trajectory file not found: %s\n' ...
             'Please generate the trajectory using TrajectoryGenerator.m.'], traj_file);
    end
    
    % Load trajectory matrices into constants6DoF.Traj
    if exist(traj_file, 'file')
        Data = readmatrix(traj_file);
        constants6DoF.Traj.Time   = Data(:, 1);
        constants6DoF.Traj.States = Data(:, 2:16);
        constants6DoF.Traj.Inputs = Data(:, 17:20);

        % Always regenerate gains, even if a file exists.
        SaveGains(name_stem, constants6DoF);
        [constants6DoF.Traj.KTGain, constants6DoF.Traj.KRGain, ...
            constants6DoF.Traj.LAGain, constants6DoF.Traj.LTGain] = ReadGains(name_stem);
        fprintf('Gain files successfully generated for %s.\n', name_stem);
    end
    
clear slBus* 
busInfo = Simulink.Bus.createObject(constants6DoF);
topLevelBusName = busInfo(end).busName;
Vehicle_Bus = evalin('base', topLevelBusName);

Waypoints = TrajectoryBuilder;
J_d = constants6DoF.J * 0.1;
MaxMdot_d = 0;
TB_d = [0.01, 0.01, 0]';

% Constant vars (varied usage)
[windMerid, windZonal] = atmoshwm(40.4258686, -86.9080655, 186 + 50);
accelBias = 0.02 * ones(3,1);
gyroBias = 0.0 * ones(3, 1);
distMode = 0;

% Sim parameters
dM_xx = 0.015;      % 3.5% Scaling from SS Rods
dM_zz = 0.010;      % 6.0% Scaling from crown
dM_xz = 0.010;      % 1.00% Coupling
dM_xy = 0.005;      % 0.50% Coupling
magDistMatrix = [dM_xx, dM_xy, dM_xz;
    dM_xy, dM_xx, dM_xz;
    dM_xz, dM_xz, dM_zz] * 0.1 + eye(3);

covar_vec = [accel_proc_cov; gyro_cov; mag_proc_cov];
IMU_Rate = 500;
Checkpoints =  [0, 5, 5,  5;
    0, 5, 10, 10;
    0, 50, 0, 0];
HoldTimeReqs = [20, 10, 2, 3];
dt_SIM = 1/500;

% MC Variables
gyroNoisePower = 10^-6;
GrommetIDX = 1;
G = GrommetSelect(GrommetIDX);
m_FC = 0.1;
kGrom = G.K;
bGrom = G.C / (2 * sqrt(kGrom * m_FC));
Kg2 = 0.03;
G_RMAX = 4;
Wind_Gain = 0.3;
Wind_Covar = 7;
lowEnd = 50;
highEnd = 800;
