function xdot = ActuatorDynamics(x, u, GND)
    xdot = zeros(6, 1);
    % Simple actuator model
    % x(1) Theta 
    % x(2) Rate 1
    % x(3) Phi 2
    % x(4) Rate 2
    % x(5) Thrust
    % x(6) Roll

    % Assume identical servos for gimbal
    a = 400;
    b = 35;

    k1 = 1;
    k2 = 1;

    theta_error = k1 * (1 - cos(x(1)));
    phi_error = k2 * (1 - cos(x(3)));

    x(3) = x(3) - theta_error;
    x(1) = x(1) - phi_error;

    xdot(1:4) = [x(2);
                a*(u(1) - x(1)) - b*x(2); 
                x(4);
                a*(u(2) - x(3)) - b*x(4)];

    % Assume first order response for thrust
    tau = 0.15;
    xdot(5) = (1/tau) * (u(3) - x(5));

    % Since torque depends on thrust changes on two fans, response time is
    % half that of thrust. We artificially add an error to the applied
    % input, proportional to the squared magnitude of the thrust pct
    roll_error = 0.05 * (u(3) / (1.5 * 9.8))^2 * (1 - GND);
    xdot(6) = (1/(tau/2)) * (u(4) + roll_error - x(6));

end