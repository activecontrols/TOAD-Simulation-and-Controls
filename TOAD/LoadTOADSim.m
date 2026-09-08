%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file loads in all the constants and parameters for the Simulink into
% workspace. Please always run this file before running a full-scale
% simulation if you've made any changes to trajectory, controls, filtering,
% or others.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize parameters and clear functions
% Initial conditions for state
clear functions;
MEKF_Constants;

%% Select vehicle
constants6DoF.Vehicle = 0; % 1 for TOAD, 0 for ASTRA

%% Create constants struct for vehicle (Approximate values, all metric)
if constants6DoF.Vehicle == 1
    % TOAD Parameters
    constants6DoF.m_dry = 141.521;
    constants6DoF.g = 9.80145; 
    constants6DoF.rTB = 0.75;
    constants6DoF.J = diag([110 110 20]);
    constants6DoF.MaxThrust = 2446.52;
    constants6DoF.MaxMdot = 1.3204;
    constants6DoF.OF = 1;
    constants6DoF.OxMass = 20.78;   constants6DoF.FuMass = 19.79;
    constants6DoF.OxHeight = 0.377; constants6DoF.FuHeight = 0.495;
    constants6DoF.OxRadius = 0.146; constants6DoF.FuRadius = 0.146;
    constants6DoF.Ox_Z = 0.85;      constants6DoF.Fu_Z = 1.35;
    constants6DoF.m_wet = constants6DoF.m_dry + constants6DoF.OxMass + constants6DoF.FuMass;
elseif constants6DoF.Vehicle == 0
    constants6DoF.m_dry = 1.275;
    constants6DoF.g = 9.80145; 
    constants6DoF.rTB = 0.26;
    constants6DoF.J = diag([0.067 0.067 0.02]);
    constants6DoF.MaxThrust = constants6DoF.m_dry * constants6DoF.g * 1.30;  % Check properly
    constants6DoF.MaxMdot = 0;
    constants6DoF.OF = 1;
    constants6DoF.OxMass = 0;       constants6DoF.FuMass = 0;
    constants6DoF.OxHeight = 1;     constants6DoF.FuHeight = 1;
    constants6DoF.OxRadius = 1;     constants6DoF.FuRadius = 1; 
    constants6DoF.Ox_Z = 1;         constants6DoF.Fu_Z = 1;
    constants6DoF.m_wet = constants6DoF.m_dry; 
end

% Dynamic Files Generation & Control
FlightDynamicsGen(constants6DoF);
x0 = [1; zeros(12,1); constants6DoF.OxMass; constants6DoF.FuMass];
u0 = [0; 0; constants6DoF.g * constants6DoF.m_wet; 0];

% Kalman Filter & Control Parameters
constants6DoF.Q = p2.Q;
constants6DoF.R = p2.obsv_cov_mat;
constants6DoF.BSigma = 5e-2;
constants6DoF.BBias = 1e-8;

% Magnetometer
constants6DoF.mag = [0.385202; 0.030609; -0.922324];
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

%% Trajectory Load
% Controller gains
% Outer Loop
max_x_trans = 1.4 * ones(1,6);
constants6DoF.Q_trans = diag(1 ./ max_x_trans.^2);
max_a_trans = 1.2; 
constants6DoF.R_trans = eye(3) .* (1 / max_a_trans^2);
constants6DoF.OmegaThr = 2.2;

% Inner Loop
max_x_rot = [0.12, 0.12, 0.12, 0.22, 0.22, 0.21];
constants6DoF.Q_rot = diag(1 ./ max_x_rot.^2);
constants6DoF.R_rot = diag([32, 32, 1/4^2]);
constants6DoF.OmegaAtt = 3.0;

% Pick a trajectory filename (e.g. "TOAD_Backflip_v001", "ASTRAv2_Circle_v001", or "Backflip_v3")
filename = "ASTRAv2_Backflip_v001";

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
    error('LoadTOADSim:TrajectoryNotFound', ...
        ['Trajectory file not found: %s\n' ...
         'Please generate the trajectory using TrajectoryGenerator.m before loading.'], traj_file);
end

% Load trajectory matrices into constants6DoF.Traj
Data = readmatrix(traj_file);
constants6DoF.Traj.Time   = Data(:, 1);
constants6DoF.Traj.States = Data(:, 2:16);
constants6DoF.Traj.Inputs = Data(:, 17:20);

% Check for gain files
gains_dir = fullfile(traj_dir, 'Gains');
gain_prefix = fullfile(gains_dir, "K_trans_" + name_stem);
gain_exists = exist(gain_prefix + ".txt", 'file') || exist(gain_prefix + ".csv", 'file') || exist(gain_prefix, 'file');

if ~gain_exists
    fprintf('Gain files not found for %s. Automatically calling SaveGains to generate TV-LQI tracking gains...\n', name_stem);
    SaveGains(name_stem, constants6DoF);
    fprintf('Gain files successfully generated for %s.\n', name_stem);
end

% Load tracking gains into constants6DoF.Traj
[constants6DoF.Traj.KTGain, constants6DoF.Traj.KRGain, ...
 constants6DoF.Traj.LAGain, constants6DoF.Traj.LTGain] = ReadGains(name_stem);

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

% MC Variables
gyroNoisePower = 10^-6;
GrommetIDX = 1;
G = GrommetSelect(GrommetIDX);
m_FC = 0.1;
kGrom = G.K;
bGrom = G.C / (2 * sqrt(kGrom * m_FC));
Kg2 = 0.03;
G_RMAX = 4;
Wind_Gain = 0.6;
Wind_Covar = 7;
lowEnd = 50;
highEnd = 800;
