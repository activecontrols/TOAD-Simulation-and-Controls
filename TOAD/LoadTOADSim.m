%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file loads in all the constants and parameters for the Simulink into
% workspace. Please always run this file before running a full-scale
% simulation if you've made any changes to trajectory, controls, filtering,
% or others.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize parameters and clear functions
% Initial conditions for state
clear functions; %#ok<CLFUNC>

%% Ensure required project paths are added, should be redundant but oh well
sim_dir = fileparts(mfilename('fullpath'));
if ~isempty(sim_dir)
    p_proj = sim_dir;
    while ~isempty(p_proj) && ~exist(fullfile(p_proj, 'Navigation'), 'dir')
        p_parent = fileparts(p_proj);
        if strcmp(p_parent, p_proj), break; end
        p_proj = p_parent;
    end
    if exist(fullfile(p_proj, 'Navigation'), 'dir')
        addpath(fullfile(p_proj, 'Navigation', 'Kalman FIlter'));
        addpath(fullfile(p_proj, 'Helper'));
        addpath(fullfile(p_proj, 'Flight Dynamics'));
        addpath(fullfile(p_proj, 'Analysis'));
        addpath(fullfile(p_proj, 'Controls'));
        addpath(fullfile(p_proj, 'Guidance'));
    end
end

MEKF_Constants;
LoadTOADParams;


%% Select vehicle (1 for TOAD, 0 for ASTRA)
if isempty(constants6DoF.Vehicle)|| ~isfield(constants6DoF, 'Vehicle')
    constants6DoF.Vehicle = 0; % 1 for TOAD, 0 for ASTRA (default ASTRA)
else
    % Normalize strings if passed in
    if ischar(constants6DoF.Vehicle) || isstring(constants6DoF.Vehicle)
        if any(strcmpi(string(constants6DoF.Vehicle), ["ASTRA", "ASTRAv2"]))
            constants6DoF.Vehicle = 0;
        elseif strcmpi(string(constants6DoF.Vehicle), "TOAD")
            constants6DoF.Vehicle = 1;
        end
    end
end

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

% Pick a trajectory filename (e.g. "TOAD_Backflip_v001", "ASTRA_Circle_v001", or "Backflip_v3")
filename = "ASTRA_Hop_v001";

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
end
% Check for gain files
    gains_dir = fullfile(traj_dir, 'Gains');
    gain_prefix = fullfile(gains_dir, "K_trans_" + name_stem);
    gain_exists = exist(gain_prefix + ".txt", 'file') || exist(gain_prefix + ".csv", 'file') || exist(gain_prefix, 'file');


if ~gain_exists && exist(traj_file, 'file')
    fprintf('Gain files not found for %s. Automatically calling SaveGains to generate TV-LQI tracking gains...\n', name_stem);
    SaveGains(name_stem, constants6DoF);
    fprintf('Gain files successfully generated for %s.\n', name_stem);
end

% Load tracking gains into constants6DoF.Traj
if gain_exists || exist(traj_file, 'file')
[constants6DoF.Traj.KTGain, constants6DoF.Traj.KRGain, ...
 constants6DoF.Traj.LAGain, constants6DoF.Traj.LTGain] = ReadGains(name_stem);
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
