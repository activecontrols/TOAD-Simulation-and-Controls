function [Q_gas, COPV_Temps_new] = COPVHeatTransfer(fillMode, dT, P_gas, T_gas, COPV_Temps, WindVel)
% Resistance network for the COPV heat transfer. Most of this code is taken
% straight from the v1 script and is essentially just cleaned up and sanity
% checked from what I did during regen 
% COPV_Temps = [Al, Carbon Fiber, Fiber Glass
    
    % Extract the temps for ease of use
    T_al = COPV_Temps(1);
    T_cf = COPV_Temps(2);
    T_fg = COPV_Temps(3);
    
    T_amb = 297.04;     % Ambient temp
    L_copv = 0.57;      % COPV Length, need to check this
    D_copv = 0.172;
    sigma = 5.670374e-8; % Stefan Boltzmann
    epsilon_fg = 0.90;   % Emissivity of outer fiberglass
    
    %% Gas & Convection Properties
    A_in = 1.1289; 
    % Equivalent diameter using the length value from the main script. Get
    % an actual value eventually
    D_in = 0.162;

    % Eval at film temp
    T_film = (T_al + T_gas) / 2;
    deltaT_in = max(abs(T_gas - T_al), 0.01);
    
    % CoolProp transport properties at current COPV state
    rho = double(py.CoolProp.CoolProp.PropsSI('D', 'P', P_gas, 'T', T_film, 'Nitrogen'));
    cp  = double(py.CoolProp.CoolProp.PropsSI('Cpmass', 'P', P_gas, 'T', T_film, 'Nitrogen'));
    k   = double(py.CoolProp.CoolProp.PropsSI('L', 'P', P_gas, 'T', T_film, 'Nitrogen'));
    mu  = double(py.CoolProp.CoolProp.PropsSI('V', 'P', P_gas, 'T', T_film, 'Nitrogen'));
    beta = double(py.CoolProp.CoolProp.PropsSI('isobaric_expansion_coefficient', ...
        'P', P_gas, 'T', T_film, 'Nitrogen'));
    
    % Rayleigh number
    g = 9.81;
    Ra_in = (g * beta * deltaT_in * (D_in^3) * (rho^2) * cp) / (mu * k);
    
    % Rayleigh Nusselt correlation. Far as i understand it, it's of the
    % form Nu = C * Ra^b, where C and B are chosen based on geometry, flow
    % conditions and such. Using turbulent flow one rn, need to read up.
    Nu_in = 0.10 * (Ra_in^(1/3));
    h_nat_in = (k / D_in) * Nu_in;

    if fillMode == 1
        h_in = max(150, h_nat_in); % Forced mixing during active fill
    else
        h_in = h_nat_in;           % Natural convection
    end
    
    %% Wall Layer Specifications
    % Aluminum Liner, value from v1
    A_al = 1.1468; t_al = 0.00200; k_al = 167; rho_al = 2700; cp_al = 900;  
    m_al = rho_al * A_al * t_al;  
    
    % Carbon Fiber Composite, values from v1
    A_cf = 1.2040; t_cf = 0.00423; k_cf = 0.5; rho_cf = 1600; cp_cf = 900;  
    m_cf = rho_cf * A_cf * t_cf;  
    
    % Fiberglass Outer Wrap, values from v1
    A_fg = 1.2480; t_fg = 0.00050; k_fg = 0.3; rho_fg = 1900; cp_fg = 800;  
    m_fg = rho_fg * A_fg * t_fg;  
    
    %% Outside Natural Convection
    % Most of this is taken blindly from v1
    A_out = 1.2529;   
    k_air = 0.024;  
    nu_air = 1.56e-5;
    deltaT = max(abs(T_fg - T_amb), 0.01);  
    Ra_L = 10^8 * deltaT * (L_copv^3);  
    Pr = 0.71;  

    % Convection to outside air, using the Churchill-Bernstein equation,
    % using the est
    Re_D = WindVel * D_copv / nu_air;
    Nu_D = 0.3 + (0.62 * Re_D^(1/2) * Pr^(1/3)) / ...
                 (1 + (0.4/Pr)^(2/3))^(1/4) * (1 + (Re_D / 282e3)^(5/8))^(4/5);
    h_conv_forced = (k_air / D_copv) * Nu_D;

    % Compute stagnant air convection as a floor
    h_conv_nat = (k_air / L_copv) * (0.68 + (0.67 * Ra_L^(1/4)) / (1 + (0.492 / Pr)^(9/16))^(4/9));
    
    % Pick the highest
    h_conv = max(h_conv_forced, h_conv_nat);

    % Radiation cooling, apparently not negligible lmao 
    h_rad = epsilon_fg * sigma * (T_fg + T_amb) * (T_fg^2 + T_amb^2);
    h_out = h_rad + h_conv;

    
    %% Thermal Resistances (Conduction & Convection)
    R_in  = 1 / (h_in * A_in);  
    R_al  = t_al / (k_al * A_al);  
    R_cf  = t_cf / (k_cf * A_cf);  
    R_fg  = t_fg / (k_fg * A_fg);  
    R_out = 1 / (h_out * A_out);  
    
    % Resistance between cell midpoints
    R_gas_al = R_in + R_al / 2;  
    R_al_cf  = R_al / 2 + R_cf / 2;  
    R_cf_fg  = R_cf / 2 + R_fg / 2;  
    R_fg_amb = R_fg / 2 + R_out;  
    
    %% Heat Fluxes
    Qdot_gas_al = (T_gas - T_al) / R_gas_al;  
    Qdot_al_cf  = (T_al - T_cf)  / R_al_cf;  
    Qdot_cf_fg  = (T_cf - T_fg)  / R_cf_fg;  
    Qdot_fg_amb = (T_fg - T_amb) / R_fg_amb;  
    
    %% Tempe updates for the nodes
    dU_al = (Qdot_gas_al - Qdot_al_cf) * dT;  
    dU_cf = (Qdot_al_cf  - Qdot_cf_fg) * dT;  
    dU_fg = (Qdot_cf_fg  - Qdot_fg_amb) * dT;  
    
    T_al_new = T_al + dU_al / (m_al * cp_al);  
    T_cf_new = T_cf + dU_cf / (m_cf * cp_cf);  
    T_fg_new = T_fg + dU_fg / (m_fg * cp_fg);  
    
    COPV_Temps_new = [T_al_new, T_cf_new, T_fg_new];
    
    % Total heat energy removed during the step.
    Q_gas = Qdot_gas_al * dT;  
end