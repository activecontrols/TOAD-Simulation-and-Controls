function [tankp, tankm, tankh, tankU, copvp, copvm, copvh, copvU, alumT, carbFibT, fibGlassT] = ...
    initializeState(tankV0, tankP0, tankT0, copvV0, copvP0, copvT0)
% initializes state of each control volume

%% tank

% density [kg/m^3]
tankp = py.CoolProp.CoolProp.PropsSI('D', 'T', tankT0, 'P', tankP0, 'Nitrogen'); 

% mass [kg]
tankm = tankp * tankV0;

% specific enthalpy [J/kg]
tankh = py.CoolProp.CoolProp.PropsSI('Hmass', 'T', tankT0, 'P', tankP0, 'Nitrogen');

% internal energy [J]
tankU = (py.CoolProp.CoolProp.PropsSI('U', 'T', tankT0, 'P', tankP0, 'Nitrogen')) * tankm;


%% copv

% density [kg/m^3]
copvp = py.CoolProp.CoolProp.PropsSI('D', 'T', copvT0, 'P', copvP0, 'Nitrogen');

% mass [kg]
copvm = copvp * copvV0;

% specific enthalpy [J/kg]
copvh = py.CoolProp.CoolProp.PropsSI('Hmass', 'T', copvT0, 'P', copvP0, 'Nitrogen');

% internal energy [J]
copvU = (py.CoolProp.CoolProp.PropsSI('U', 'T', copvT0, 'P', copvP0, 'Nitrogen')) * copvm;

% layer temperatures [K]

alumT = 297.039; 

carbFibT = 297.039;

fibGlassT = 297.039;

end

