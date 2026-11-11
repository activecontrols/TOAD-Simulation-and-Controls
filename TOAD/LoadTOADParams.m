function constants6DoF = LoadTOADParams(Vehicle)

constants6DoF.Vehicle = Vehicle;
this_dir = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(this_dir, 'Navigation')));
MEKF_Constants;
%% Create constants struct for vehicle (Approximate values, all metric)
if constants6DoF.Vehicle == 1
    % TOAD Parameters
    constants6DoF.m_dry = 181.437;
    constants6DoF.g = 9.80145; 
    constants6DoF.rTB = 0.75;
    constants6DoF.J = diag([110 110 20]);
    constants6DoF.MaxThrust = 2935.83;
    constants6DoF.MaxMdot = 1.5513;
    constants6DoF.OF = 1.2;
    constants6DoF.OxMass = 33.42;   constants6DoF.FuMass = 26.181;
    constants6DoF.OxHeight = 0.377; constants6DoF.FuHeight = 0.495;
    constants6DoF.OxRadius = 0.146; constants6DoF.FuRadius = 0.146;
    constants6DoF.Ox_Z = 0.85;      constants6DoF.Fu_Z = 1.35;
    constants6DoF.m_wet = constants6DoF.m_dry + constants6DoF.OxMass + constants6DoF.FuMass;
elseif constants6DoF.Vehicle == 0
    % ASTRA Parameters
    constants6DoF.m_dry = 1.1690;
    constants6DoF.g = 9.80145; 
    constants6DoF.rTB = 0.26;
    constants6DoF.J = diag([0.067 0.067 0.02]);
    constants6DoF.MaxThrust = 16;  
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

constants6DoF.K_trans = [9.583831e-01, 0, 0, 1.683824e+00, 0, 0, 0, 9.583831e-01, 0, 0, 1.683824e+00, 0, 0, 0, 9.583831e-01, 0, 0, 1.683824e+00];
constants6DoF.K_rot = [-7.541335e-01, -8.115214e-05, -1.362233e-04, -2.208310e-01, 2.617803e-06, -2.971590e-05, 8.246603e-05, -7.541662e-01, -1.921107e-03, -2.577455e-06, -2.208326e-01, 6.048713e-04, -5.516150e-05, -5.104331e-03, 8.236951e-01, 2.019352e-05, -4.129795e-04, 2.299734e-01];
constants6DoF.L_Thrust = [5.616296e-03, 0, 0, 2.005830e-03, 0, 0; 
    0, 4.727018e-03, 4.929804e-06 0, 2.010921e-03, 2.808143e-08; 
    0, 5.023681e-06, 5.088132e-03, 0, -1.601761e-07, 2.006986e-03; 
    1.889238e-04, -9.886904e-10, 0, 1.123218e-02, 0, 0;
    0, 9.157307e-04, -2.406860e-05, 0, 9.845408e-03, -4.535823e-07;
    0, 3.624535e-06, 2.974015e-04, 0, -3.739723e-07, 9.309009e-03;
    5.350414e-04, -2.769941e-09, 0, 1.575495e-02, 1.030446e-10, 0;
    0, 2.273589e-03, -5.977943e-05, 0, 1.205917e-02, -1.064738e-06;
    0, 8.297828e-06, 6.894103e-04, 0, -8.723946e-07, 1.079742e-02];
constants6DoF.L_Att = [6.808373e-03, 0, 0, 0, 8.123445e-03, 0, 0, 0, 8.123349e-03, 5.788234e-03, 0, 0, 0, 8.242795e-03, 0, 0, 0, 8.194591e-03];


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
