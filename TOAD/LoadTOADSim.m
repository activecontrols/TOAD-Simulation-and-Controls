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
        max_x_trans = 1.3 * ones(1,6);
    else
        % ASTRAv2 Tuning
        max_x_trans = 1.25 * ones(1,6);
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
        max_x_rot = [0.15, 0.15, 0.12, 0.75, 0.75, 0.7];
        constants6DoF.R_rot = diag([70, 70, 1/0.10^2]);
    end

    constants6DoF.Q_rot = diag(1 ./ max_x_rot.^2);
    constants6DoF.OmegaAtt = 2.2;

    constants6DoF.K_Trans = [9.583831e-01, 0, 0, 1.683824e+00, 0, 0;
        0, 9.583831e-01, 0, 0, 1.683824e+00, 0;
        0, 0, 9.583831e-01, 0, 0, 1.683824e+00];
    constants6DoF.K_Rot = [-7.541335e-01, -8.115214e-05, -1.362233e-04, -2.208310e-01, 2.617803e-06, -2.971590e-05;
        8.246603e-05, -7.541662e-01, -1.921107e-03, -2.577455e-06, -2.208326e-01, 6.048713e-04;
        -5.516150e-05, -5.104331e-03, 8.236951e-01, 2.019352e-05, -4.129795e-04, 2.299734e-01];
    constants6DoF.L_Thrust = [5.616296e-03, 0, 0, 2.005830e-03, 0, 0; 
        0, 4.727018e-03, 4.929804e-06 0, 2.010921e-03, 2.808143e-08; 
        0, 5.023681e-06, 5.088132e-03, 0, -1.601761e-07, 2.006986e-03; 
        1.889238e-04, -9.886904e-10, 0, 1.123218e-02, 0, 0;
        0, 9.157307e-04, -2.406860e-05, 0, 9.845408e-03, -4.535823e-07;
        0, 3.624535e-06, 2.974015e-04, 0, -3.739723e-07, 9.309009e-03;
        5.350414e-04, -2.769941e-09, 0, 1.575495e-02, 1.030446e-10, 0;
        0, 2.273589e-03, -5.977943e-05, 0, 1.205917e-02, -1.064738e-06;
        0, 8.297828e-06, 6.894103e-04, 0, -8.723946e-07, 1.079742e-02];
    constants6DoF.L_Att = [6.808373e-03, 0, 0;
        0, 8.123445e-03, 0;
        0, 0, 8.123349e-03;
        5.788234e-03, 0, 0;
        0, 8.242795e-03, 0;
        0, 0, 8.194591e-03];

%% Trajectory Params
    % Pick a trajectory filename (e.g. "TOAD_Backflip_v001",
    % "ASTRA_Circle_v001")
    filename = "ASTRA_Circle_v003";
    
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

J_d = constants6DoF.J * 0.1;
MaxMdot_d = 0;
TB_d = [0.01, 0.01, 0]';

% Constant vars (varied usage)
[windMerid, windZonal] = atmoshwm(40.4258686, -86.9080655, 186 + 50);
accelBias = 0.01 * ones(3,1);
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

%V5 MC variables

gimbalBias = zeros(4,1);
magHardIron = 0.05 * ones(3,1);
sigma_accel_mu = 1e-4;              % TODO: Allan variance / warm-up log
accel_sigma_range = [0.5e-4, 2e-4];
sigma_gyro_mu = deg2rad(1e-3);
gyro_sigma_range = [0, 0];   
IMUDriftParams = [zeros(1,8),1];
seedDryden = [23341 23342 23343 23344];
seedAccel = [23093 23094 23095];
seedGyro = [23093 23094 23095];
seedMagNoise = [23341,23342,23343];
