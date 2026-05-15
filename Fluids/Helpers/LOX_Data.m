function [h, Cp, rho, A, B] = LOX_Data(T, P)
    % Linear Cp fit for Liquid Oxygen
    A = 1600; 
    B = 1.0;  
    
    Cp = A + B .* T;
    h  = A .* T + 0.5 .* B .* T.^2; % Absolute enthalpy ref at 0K
    
    % Linear density fit for cryogenic LOX (bounded to prevent negative values at high T)
    rho = max(1141 - 4.5 .* (T - 90), 10); 
end