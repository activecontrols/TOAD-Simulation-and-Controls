function xdot = ActuatorDynamics(x, u)
    xdot = zeros(6, 1);
    
    % 1st-order gimbal dynamics (tau = 25 ms)
    tau_gimbal = 0.025;
    xdot(1) = (1 / tau_gimbal) * (u(1) - x(1));  % dTheta/dt
    xdot(2) = 0;                                 % Unused 
    
    xdot(3) = (1 / tau_gimbal) * (u(2) - x(3));  % dPhi/dt
    xdot(4) = 0;                                 % Unused 
    
    % Thrust and Roll dynamics
    tau = 0.03;
    xdot(5) = (1 / tau) * (u(3) - x(5));
    xdot(6) = (1 / (tau / 2)) * (u(4) - x(6));
end