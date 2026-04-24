function x = ThrustPlant(t)
%THRUSTPLANT Summary of this function goes here
%   random plant function for the thrust in the engine

x = cosint(t*pi/180);

end