function xdot = TOAD_Dynamics(x, u)
[xdot_6DoF, xdot_Mass, J_tot, netTau] = RawDynamics(x, u);

% Numerical calculation of angular acceleration
omegaBdot = J_tot \ netTau;

% Full state derivative
xdot = [xdot_6DoF; omegaBdot; xdot_Mass];
end