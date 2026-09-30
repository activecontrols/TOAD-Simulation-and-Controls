function [copvU1, copvP1, copvT1, copvh1] = updateState2(outQ, copvU0, copvp, copvm)
% updates pressure, temperature, internal energy, and enthalpy of copv due
% to heat transfer

copvU1 = copvU0 - outQ; 
copvu1 = copvU1 / copvm;

copvP1 = py.CoolProp.CoolProp.PropsSI('P','Umass', copvu1,'Dmass', copvp, 'Nitrogen');
copvT1 = py.CoolProp.CoolProp.PropsSI('T','Umass', copvu1,'Dmass', copvp, 'Nitrogen');
copvh1 = py.CoolProp.CoolProp.PropsSI('Hmass','Umass', copvu1,'Dmass', copvp, 'Nitrogen');

end