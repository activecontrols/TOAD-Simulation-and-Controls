function [tankV, tankP0, tankT0, copvV, copvP0, copvT0, timeStep, regSetP] = initializeModel()
% sets up the model parameters

% tank params
tankV = 6 * 0.05; %m^3
tankP0 = 4.1369e+7; %Pa (6000 psia)
tankT0 = 297.039; %Kelvin (75 deg. F)

% copv params
copvV = 0.036; %m^3
copvP0 = 101325; %Pa (ambient)
copvT0 = 297.039; %Kelvin (75 deg. F); 

% timestep
timeStep = 1; %seconds

% reg set pressure

regSetP = 3.1026e+7; %Pa (4500 psia)



end