function eta = LogSE23(Chi)
% LogSE23
%
% Logarithm map for SE_2(3).
%
% Output:
%   eta(1:3)  attitude
%   eta(4:6)  velocity
%   eta(7:9)  position

R = Chi(1:3,1:3);
v = Chi(1:3,4);
p = Chi(1:3,5);

% Compute SO(3) logarithm.
cosTheta = (trace(R)-1)/2;

% Protect acos against numerical values slightly outside [-1,1].
cosTheta = max(-1,min(1,cosTheta));

theta = acos(cosTheta);

if theta < 1e-8

    Phi = 0.5*(R-R');

    phi = [Phi(3,2);
           Phi(1,3);
           Phi(2,1)];

else

    Phi = ...
        theta/(2*sin(theta)) ...
        * (R-R');

    phi = [Phi(3,2);
           Phi(1,3);
           Phi(2,1)];

end

% Inverse SO(3) left Jacobian.
theta = norm(phi);
Phi = zetaCross(phi);

if theta < 1e-6

    Jinv = ...
        eye(3) ...
        - 0.5*Phi ...
        + (1/12)*Phi^2;

else

    coefficient = ...
        1/theta^2 ...
        - (1+cos(theta)) ...
          /(2*theta*sin(theta));

    Jinv = ...
        eye(3) ...
        - 0.5*Phi ...
        + coefficient*Phi^2;

end

dv = Jinv*v;
dp = Jinv*p;

eta = [phi;
       dv;
       dp];

end