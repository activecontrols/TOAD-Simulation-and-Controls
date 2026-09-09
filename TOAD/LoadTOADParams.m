constants6DoF.Vehicle = 0;

MEKF_Constants;
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
    % ASTRA Parameters
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
% Ensure working directory is project root for relative file writes, once
% again should be redundant but needed. Failing when working on diffrent
% machines
if ~exist(fullfile(pwd, 'Flight Dynamics'), 'dir')
    p_check = fileparts(mfilename('fullpath'));
    while ~isempty(p_check) && ~exist(fullfile(p_check, 'Flight Dynamics'), 'dir')
        parent_check = fileparts(p_check);
        if strcmp(parent_check, p_check), break; end
        p_check = parent_check;
    end
    if exist(fullfile(p_check, 'Flight Dynamics'), 'dir')
        cd(p_check);
    end
end
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