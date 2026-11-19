function [mu, k, rho, cp, saturation_temp, state] = LOX_properties(T, P)
    % Thermophysical properties of Oxygen
    % Inputs: T [Kelvin], P [Pascals]
    
    %% Constants
    % Molar mass of O2: 31.999 kg/kmol
    % Critical Point: Tc = 154.58 K, Pc = 5.043 MPa
    R = 259.8; % Gas constant for O2 [J/kgK]
    
    %% Saturation Curve (Antoine-style approximation for O2)
    % P in Pascals, T in Kelvin
    P_bar = P / 1e5;
    if P_bar < 50.43
        saturation_temp = 100 / (1 - 0.1 * log(P_bar / 1.01325)) * 0.901; % Approx curve
    else
        saturation_temp = 154.58; % Supercritical
    end
    
    %% State Determination
    if T >= saturation_temp || P < 1000
        state = "gas";
    else
        state = "liquid";
    end
    
    %% Properties
    if strcmp(state, "liquid")
        % Liquid properties near 90K (1 atm boiling point)
        rho = 1141; % [kg/m^3] - relatively incompressible
        cp = 1.70e3; % [J/kgK]
        mu = 1.9e-4; % Viscosity [Pa*s]
        k = 0.15; % Thermal conductivity [W/mK]
    else
        % Gas properties (Ideal Gas assumption for real-time speed)
        rho = P / (R * T); % [kg/m^3]
        cp = 0.918e3; % [J/kgK]
        mu = 2.0e-5 * (T / 293)^0.7; % Sutherland's approx for viscosity
        k = 0.026 * (T / 293)^0.8; % Approx thermal conductivity
    end
end