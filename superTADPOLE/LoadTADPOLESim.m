% fluid mechanics conditions
constantsSTADPOLE.init_cond = [0.4; 0.4; 1e+05; 3.447e+6; 3.447e+6]; % initial conditions for fluid mechanics
constantsSTADPOLE.c_f = 1.1; % thrust coefficient
constantsSTADPOLE.a_t = 1.0948365 * 1e-3; % [m^2], throat area

constantsSTADPOLE.dens_o = 1141; % [kg/m^3], liquid oxygen density, at 90 K and 1 atm
constantsSTADPOLE.dens_f = 786; % [kg/m^3], fuel density, at 293 K and 1 atm
constantsSTADPOLE.dens_w = 1000; % [kg/m^3], water density
constantsSTADPOLE.g = 9.81; % [m/s^2], acceleration due to gravity
constantsSTADPOLE.a_i_o = 3.14577787810908e-5; % [m^2] area injector ox
constantsSTADPOLE.a_i_f = 2.60064441031368e-5; % [m^2] area injector fuel
constantsSTADPOLE.r = 442.234043; % [J/kg-K], specific gas constant of exhaust (molecular weight 18.8 g/mol)
constantsSTADPOLE.c_star = 1750; % [m/s] characteristic velocity
constantsSTADPOLE.OF_target = 1;
% constantsSTADPOLE.OF_engine = 1.2; calculate it in the engine with mass
% of ox over mass of fuel
STADPOLE = Simulink.Bus.createObject(constantsSTADPOLE);