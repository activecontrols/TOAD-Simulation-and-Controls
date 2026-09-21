function Chi = ExpSE23(phi, dv, dp)
% ExpSE23
%
% Exponential map for the SE_2(3) correction.
%
% Input ordering here is:
%   phi  attitude correction
%   dv   velocity correction
%   dp   position correction

theta = norm(phi);
Phi = zetaCross(phi);

if theta < 1e-8

    R = eye(3) ...
        + Phi ...
        + 0.5*Phi^2;

    J = eye(3) ...
        + 0.5*Phi ...
        + (1/6)*Phi^2;

else

    R = eye(3) ...
        + (sin(theta)/theta)*Phi ...
        + ((1-cos(theta))/theta^2)*Phi^2;

    J = eye(3) ...
        + ((1-cos(theta))/theta^2)*Phi ...
        + ((theta-sin(theta))/theta^3)*Phi^2;

end

Chi = [R,          J*dv, J*dp;
       zeros(1,3), 1,    0;
       zeros(1,3), 0,    1];

end