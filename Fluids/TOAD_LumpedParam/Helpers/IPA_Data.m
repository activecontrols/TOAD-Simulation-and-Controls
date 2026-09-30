function [h, Cp, rho, A, B] = IPA_Data(T, P)
    % Linear Cp fit for Isopropyl Alcohol
    A = 1500;
    B = 4.0;
    
    Cp = A + B .* T;
    h  = A .* T + 0.5 .* B .* T.^2;
    
    rho = max(786 - 0.8 .* (T - 293), 10); 
end