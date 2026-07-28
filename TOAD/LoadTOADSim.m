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
clear GroundEstimator;
clear FlightEstimator2;
clear SensorSimulation;
clear GPS_Sim;
clear DigitalNF;
clear MEKF;
clear slBus1
MEKF_Constants;

%% Create constants struct for TOAD (Approximate values, all metric)
% Vehicle Parameters
constantsTOAD.m_dry = 141.521;
constantsTOAD.g = 9.80145; 
constantsTOAD.rTB = 0.75;
constantsTOAD.J = diag([110 110 20]);
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
[K_Att_Wet, ~] = TOAD_Controller_Gen(constantsTOAD, constantsTOAD.OxMass, constantsTOAD.FuMass);
[K_Att_Dry, ~] = TOAD_Controller_Gen(constantsTOAD, 0, 0);
x0 = [1; zeros(12,1); constantsTOAD.OxMass; constantsTOAD.FuMass];
u0 = [0; 0; constantsTOAD.g * constantsTOAD.m_wet; 0];

% Kalman Filter & Control Parameters
constantsTOAD.Q = p2.Q;
constantsTOAD.R = p2.obsv_cov_mat;
constantsTOAD.BSigma = 5e-2;
constantsTOAD.BBias = 1e-8;

% Magnetometer
constantsTOAD.mag = [0.385202; 0.030609; -0.922324];
dM_xx = 0.035;      % 3.5% Scaling from SS Rods
dM_zz = 0.060;      % 6.0% Scaling from crown
dM_xz = 0.010;      % 1.00% Coupling
dM_xy = 0.005;      % 0.50% Coupling
magDistMatrix = [dM_xx, dM_xy, dM_xz;
                 dM_xy, dM_xx, dM_xz;
                 dM_xz, dM_xz, dM_zz] * 0.1 + eye(3);

constantsTOAD.K_Att_Wet = K_Att_Wet;
constantsTOAD.K_Att_Dry = K_Att_Dry;
covar_vec = [accel_proc_cov; gyro_cov; mag_proc_cov];
IMU_Rate = 1000;
Checkpoints =  [0, 5, 5,  5;
                0, 5, 10, 10;
                0, 50, 0, 0];
HoldTimeReqs = [20, 10, 2, 3];
dt_SIM = 1/1000;

%% Trajectory Load
% Pick a trajectory filename 
filename = "Backflip_v1.csv";
Data = readmatrix("Guidance\Trajectories\"+filename);
constantsTOAD.Traj.Time = Data(:, 1);
constantsTOAD.Traj.States = Data(:, 2:16);
constantsTOAD.Traj.Inputs = Data(:, 17:20);
[constantsTOAD.Traj.FBGain, constantsTOAD.Traj.FBCost] = ReadGains(filename);

% Controller gains
max_x = [0.4, 0.4, 0.1, 2, 2, 2, 3, 3, 3, 0.8, 0.8, 0.5];
constantsTOAD.Q_Control = eye(12) .* 1 ./ max_x.^2;
constantsTOAD.R_Control = diag([12, 12, 1/100^2, 5]);

clear slBus* 
busInfo = Simulink.Bus.createObject(constantsTOAD);
topLevelBusName = busInfo(end).busName;
TOAD_Bus = evalin('base', topLevelBusName);

Waypoints = TrajectoryBuilder;
J_d = zeros(3);
MaxMdot_d = 0;
TB_d = zeros(3,1);

% Constant vars (varied usage)
[windMerid, windZonal] = atmoshwm(40.4258686, -86.9080655, 186 + 50);
accelBias = 0.09 * ones(3,1);
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
Wind_Gain = 1;
Wind_Covar = 1;
lowEnd = 50;
highEnd = 800;
