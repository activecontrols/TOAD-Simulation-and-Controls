function x_next = RK4(f, x_k, System, dT)
    % RK4 integrator for x_dot = f(x, u) with fixed step size dt.
    % Inputs:
    %   f: function handle for dynamics, f(x, u)
    %   x_k: current state (column vector)
    %   u: input (assumed constant over the step)
    %   dt: time step (scalar)
    % Output:
    %   x_next: next state (column vector)
    
    k1 = f(x_k, System);
    k2 = f(x_k + (1/2) * k1 * dT, System);
    k3 = f(x_k + (1/2) * k2 * dT, System);
    k4 = f(x_k + dT * k3, System);
    
    x_next = x_k + (dT/6) * (k1 + 2*k2 + 2*k3 + k4);
end