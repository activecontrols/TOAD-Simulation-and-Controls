%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file loads in all the constants and parameters for the Simulink into
% workspace. Please always run this file before running a full-scale
% simulation if you've made any changes to trajectory, controls, filtering,
% or others.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize parameters and clear functions
% Initial conditions for state
clear ref_generator3;
clear inputfcn3;
clear EstimateStateFCN;
clear SensorSimulation;
clear GPS_Sim;
clear DigitalNF;
clear slBus1
MEKF_Constants;

%% Create constants struct for TOAD (Approximate values, all metric)
% Vehicle Parameters
constantsTOAD.m_dry = 132.5;
constantsTOAD.g = 9.80145; 
constantsTOAD.rTB = 1.1;
constantsTOAD.J = diag([50 50 15]);
constantsTOAD.MaxThrust = 2446.52;
constantsTOAD.MaxMdot = 1.3204;
constantsTOAD.OF = 1;
constantsTOAD.OxMass = 20.78;   constantsTOAD.FuMass = 19.79;
constantsTOAD.OxHeight = 0.377; constantsTOAD.FuHeight = 0.495;
constantsTOAD.OxRadius = 0.146; constantsTOAD.FuRadius = 0.146;
constantsTOAD.Ox_Z = 0.85;      constantsTOAD.Fu_Z = 1.35;
constantsTOAD.m_wet = constantsTOAD.m_dry + constantsTOAD.OxMass + constantsTOAD.FuMass;

% Dynamic Files Generation & Control
% FlightDynamicsGen(constantsTOAD);
[K_Att, ~] = TOAD_Controller_Gen(constantsTOAD, constantsTOAD.OxMass, constantsTOAD.FuMass);
x0 = [1; zeros(12,1); constantsTOAD.OxMass; constantsTOAD.FuMass];
u0 = [0; 0; constantsTOAD.g * constantsTOAD.m_wet; 0];

% Kalman Filter & Control Parameters
constantsTOAD.Q = p2.Q;
constantsTOAD.R = p2.obsv_cov_mat;
constantsTOAD.BSigma = 5e-2;
constantsTOAD.BBias = 1e-8;
constantsTOAD.mag = [cos(pi/6); 0; -sin(pi/6)];
magDistMatrix = eye(3) + 0.02 * randn(3);
magDistMatrix = (magDistMatrix + magDistMatrix') / 2;
magDistMatrix = eye(3);
constantsTOAD.K_Att = K_Att;
covar_vec = [accel_proc_cov; gyro_cov; mag_proc_cov];
IMU_Rate = 1000;
Checkpoints =  [0, 5, 5,  5;
                0, 5, 10, 10;
                0, 50, 0, 0];
HoldTimeReqs = [20, 10, 2, 3];
dt_SIM = 1/250;

% Create slBus
TOAD = Simulink.Bus.createObject(constantsTOAD);
Waypoints = TrajectoryBuilder;
J_d = zeros(3);
MaxMdot_d = 0;
TB_d = zeros(3,1);

% Constant vars (varied usage)
[windMerid, windZonal] = atmoshwm(40.4258686, -86.9080655, 186 + 50);
accelBias = 0.09 * ones(3,1);
gyroBias = 0.05 * ones(3, 1);
distMode = 0;

%% Load the data dictionary
% dictObj = Simulink.data.dictionary.open('Model_Vars_2024b.sldd');
% importFromBaseWorkspace(dictObj, 'varList', {'accelBias', 'constantsTOAD', 'distMode', 'dt', ...
%     'dt_SIM', 'gyroBias', 'IMU_Rate', 'J_d', 'magDistMatrix', 'MaxMdot_d', 'slBus1', 'TB_d', ...
%     'u0', 'x0', 'Waypoints', 'windMerid', 'windZonal'}, 'existingVarsAction', 'overwrite');
% saveChanges(dictObj);