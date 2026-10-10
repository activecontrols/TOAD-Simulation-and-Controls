clear slBus1
addpath('Fluids\');
addpath('Thrust\');
addpath('Valves\');


% *************************************************************************
% Fluid dynamics initial conditions

constantsSTADPOLE.init_cond = [0.8795; 0.7329; 1723690; 2240797; 2068428; 2e6; 2e6]; 
% 1: mass flow oxygen across pipe [kg/s]
% 2: mass flow fuel across pipe [kg/s]
% 3: chamber pressure [Pa]
% 4: oxygen pressure at injector [Pa]
% 5: fuel pressure at injector [Pa]
% 6: oxygen pressure at valve [Pa] ------------ NEED VALUE
% 7: fuel pressure at valve [Pa] -------------- NEED VALUE


% *************************************************************************
% Oxygen related constants

constantsSTADPOLE.dens_o = 1130;                    % [kg/m^3], liquid oxygen density, at 90 K and 1 atm
constantsSTADPOLE.a_i_o = 3.8791e-5;                % [m^2], oxygen injector area
constantsSTADPOLE.d_coeff_ox = 0.66;                % [unitless], oxygen discharge coefficient (injector)

% Oxygen tank pressure = 3.7921E+06 [Pa]
% Oxygen temperature = 90.17 [K]


% *************************************************************************
% Fuel related constants

constantsSTADPOLE.dens_f = 691.4185;                % [kg/m^3], fuel density (at injector), at 293 K and 1 atm
constantsSTADPOLE.a_i_f = 4.4991e-5;                % [m^2], fuel injector area
constantsSTADPOLE.d_coeff_fu = 0.749;               % [unitless], fuel discharge coefficient (injector)

% Fuel tank pressure = 3.7921E+06 [Pa]
% Fuel temperature = 297 [K]


% *************************************************************************
% Other constants

constantsSTADPOLE.c_f = 1.334;                      % [unitless], thrust coefficient
constantsSTADPOLE.a_t = 1.2809e-3;                  % [m^2], throat area
constantsSTADPOLE.dens_w = 997;                     % [kg/m^3], water density
constantsSTADPOLE.g = 9.8066;                       % [m/s^2], acceleration due to gravity
constantsSTADPOLE.r = 442.234043;                   % [J/kg-K], specific gas constant of exhaust (molecular weight 18.8 g/mol)
constantsSTADPOLE.c_star = 1419;                    % [m/s] characteristic velocity
constantsSTADPOLE.OF_target = 1.2;                  % [unitless], oxygen:fuel target ratio


% constantsSTADPOLE.OF_engine = 1.2; calculate it in the engine with mass
% of ox over mass of fuel


% *************************************************************************
STADPOLE = Simulink.Bus.createObject(constantsSTADPOLE);