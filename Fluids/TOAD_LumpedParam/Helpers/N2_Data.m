function [h, Cp, rho, A, B] = N2_Data(T, P)
    % N2 Gas (kept B small to ensure stability up to 3000K in the combustor)
    A = 1000;
    B = 0.1;
    
    Cp = A + B .* T;
    h  = A .* T + 0.5 .* B .* T.^2;
    
    R_N2 = 296.8;
    rho = P ./ (R_N2 .* T); % Ideal gas density
end