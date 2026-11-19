function [mdot, isChoked, regOutP, pressureEq] = modelRegulator(tankP, tankT, copvP, regSetP, tankP0)
%determines effects due to regulator (by far the lowest Cv out of any
%instrument -> we know that this will set the mass flow.

pressureEq = 0;

% reg CdA (according to specs)

CdA = 2.4828e-6; %0.07 inch diameter -> area in m^2

% reg outlet pressure increase due to change in supply pressure drop (according to specs)

regOutP = regSetP + 0.07 * (tankP0 - tankP);

% mass flow and choking through reg

if tankP > copvP

    [mdot, isChoked] = calculateMassFlow(tankP, tankT, copvP, CdA);

else

    pressureEq = 1;
    mdot = 0;
    isChoked = 0;

end

end