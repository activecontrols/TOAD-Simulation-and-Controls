function [tankm1, tankU1, tankP1, tankT1, tankh1, tankp1, copvm1, copvU1, copvP1, copvT1, copvh1, copvp1] ...
    = updateState(tankm0, tankU0, tankV, tankh0, copvm0, copvU0, copvV, mdot, timeStep)
% update state of tank and copv (mass, internal energy, density, spec.
% internal energy, pressure, and temperature. (in order))

% tank
tankm1 = tankm0 - mdot * timeStep; 
tankU1 = tankU0 - (mdot * tankh0) * timeStep;
tankp1 = tankm1 / tankV;
tanku1 = tankU1 / tankm1;

tankP1 = py.CoolProp.CoolProp.PropsSI('P','Umass', tanku1,'Dmass', tankp1, 'Nitrogen');
tankT1 = py.CoolProp.CoolProp.PropsSI('T','Umass', tanku1,'Dmass', tankp1, 'Nitrogen');
tankh1 = py.CoolProp.CoolProp.PropsSI('Hmass','Umass', tanku1,'Dmass', tankp1, 'Nitrogen');

% copv
copvm1 = copvm0 + mdot * timeStep;
copvU1 = copvU0 + (mdot * tankh0) * timeStep;
copvp1 = copvm1 / copvV;
copvu1 = copvU1 / copvm1;

copvP1 = py.CoolProp.CoolProp.PropsSI('P','Umass', copvu1,'Dmass', copvp1, 'Nitrogen');
copvT1 = py.CoolProp.CoolProp.PropsSI('T','Umass', copvu1,'Dmass', copvp1, 'Nitrogen');
copvh1 = py.CoolProp.CoolProp.PropsSI('Hmass','Umass', copvu1,'Dmass', copvp1, 'Nitrogen');



end
