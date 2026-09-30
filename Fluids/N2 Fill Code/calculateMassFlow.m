function [mdot, isChoked] = calculateMassFlow(Pup, Tup, Pdown, CdA)
% calculates mass flow

cp = py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', Pup, 'T', Tup, 'Nitrogen');

cv = py.CoolProp.CoolProp.PropsSI('Cvmass', 'P', Pup, 'T', Tup, 'Nitrogen');

gamma = cp / cv; %specific heat ratio

R = py.CoolProp.CoolProp.PropsSI('GAS_CONSTANT', 'Nitrogen') / ...
    py.CoolProp.CoolProp.PropsSI('MOLAR_MASS', 'Nitrogen'); %specific gas constant

isChoked = (Pdown / Pup) <= ((2 / (gamma + 1)) ^ (gamma / (gamma - 1)));

if isChoked
    mdot = CdA * Pup * sqrt(gamma / (R * Tup)) * ...
        (2 / (gamma + 1))^((gamma + 1) / (2 * (gamma - 1)));
else
    mdot = CdA * Pup * sqrt((2 * gamma / (R * Tup * (gamma - 1))) * ...
        ((Pdown/Pup)^(2 / gamma) - (Pdown/Pup)^((gamma + 1) / gamma)));
end

end