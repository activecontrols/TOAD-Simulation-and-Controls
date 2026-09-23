function constants6DoF = LoadTOADParams(Vehicle)

constants6DoF.Vehicle = Vehicle;
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, 'Navigation')));
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

%% Controller params (Gains and Tuning)
% Outer Loop
if constants6DoF.Vehicle == 1
    % TOAD Tuning
    max_x_trans = 1.4 * ones(1,6);
else
    % ASTRA Tuning
    max_x_trans = 1.5 * ones(1,6);
end
constants6DoF.Q_trans = diag(1 ./ max_x_trans.^2);
max_a_trans = 1.2; 
constants6DoF.R_trans = eye(3) .* (1 / max_a_trans^2);
constants6DoF.OmegaThr = 2.2;

% Inner Loop
if constants6DoF.Vehicle == 1
    % TOAD Tuning
    max_x_rot = [0.12, 0.12, 0.12, 0.22, 0.22, 0.21];
    constants6DoF.R_rot = diag([35, 35, 1/4^2]);
else
    % ASTRA Tuning
    max_x_rot = [0.15, 0.15, 0.12, 0.5, 0.5, 0.7];
    constants6DoF.R_rot = diag([60, 60, 1/0.1^2]);
end

constants6DoF.Q_rot = diag(1 ./ max_x_rot.^2);
constants6DoF.OmegaAtt = 3.2;


end
