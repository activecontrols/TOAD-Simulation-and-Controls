function [Q_gas, ThermalState_new] = COPVThermalModel(P_gas, T_gas, ThermalState, dt, isFlowing, WindVel)
% COPVTHERMALMODEL 3-layer resistance network for COPV thermal transient.
% Adapted from N2 Fill Code v2 (Aluminum liner, Carbon Fiber, Fiberglass).
%
% Inputs:
%   P_gas            - Gas pressure [Pa]
%   T_gas            - Gas temperature [K]
%   ThermalState     - Struct with field .T_layers = [T_al, T_cf, T_fg]
%   dt               - Time step [s]
%   isFlowing        - Boolean (true if active blowdown/fill)
%   WindVel          - Wind velocity [m/s] (default 1.5)
%
% Outputs:
%   Q_gas            - Net heat added to fluid gas [J] (negative if heat leaves gas)
%   ThermalState_new - Updated thermal state struct

    if nargin < 6 || isempty(WindVel), WindVel = 1.5; end

    T_al = ThermalState.T_layers(1);
    T_cf = ThermalState.T_layers(2);
    T_fg = ThermalState.T_layers(3);

    T_amb = 297.04;
    L_copv = 0.57;
    D_copv = 0.172;
    sigma = 5.670374e-8;
    epsilon_fg = 0.90;

    % Internal Convection
    A_in = 1.1289;
    D_in = 0.162;
    T_film = (T_al + T_gas) / 2;
    deltaT_in = max(abs(T_gas - T_al), 0.01);

    % Gas properties at film temp (single call, zero CoolProp overhead)
    props_film = FluidProperties('Nitrogen', 'From_P_T', P_gas, T_film);
    rho = props_film.rho;
    cp  = props_film.cp;
    k   = props_film.k;
    mu  = props_film.mu;

    beta = 1.0 / max(T_film, 50);
    g = 9.81;
    Ra_in = (g * beta * deltaT_in * (D_in^3) * (rho^2) * cp) / max(mu * k, 1e-12);
    Nu_in = 0.10 * (max(Ra_in, 1e2)^(1/3));
    h_nat_in = (k / D_in) * Nu_in;

    if isFlowing
        h_in = max(150, h_nat_in); % Enhanced forced convection during active flow
    else
        h_in = h_nat_in;
    end

    % Wall Layer Specs
    A_al = 1.1468; t_al = 0.00200; k_al = 167; rho_al = 2700; cp_al = 900;
    m_al = rho_al * A_al * t_al;

    A_cf = 1.2040; t_cf = 0.00423; k_cf = 0.5; rho_cf = 1600; cp_cf = 900;
    m_cf = rho_cf * A_cf * t_cf;

    A_fg = 1.2480; t_fg = 0.00050; k_fg = 0.3; rho_fg = 1900; cp_fg = 800;
    m_fg = rho_fg * A_fg * t_fg;

    % Outside Air Convection & Radiation
    A_out = 1.2529; k_air = 0.024; nu_air = 1.56e-5; Pr = 0.71;
    deltaT_out = max(abs(T_fg - T_amb), 0.01);
    Ra_L = 10^8 * deltaT_out * (L_copv^3);

    Re_D = WindVel * D_copv / nu_air;
    Nu_D = 0.3 + (0.62 * sqrt(Re_D) * (Pr^(1/3))) / ...
                 (1 + (0.4/Pr)^(2/3))^(1/4) * (1 + (Re_D / 282e3)^(5/8))^(4/5);
    h_conv_forced = (k_air / D_copv) * Nu_D;
    h_conv_nat = (k_air / L_copv) * (0.68 + (0.67 * (Ra_L^(1/4))) / (1 + (0.492 / Pr)^(9/16))^(4/9));
    h_conv = max(h_conv_forced, h_conv_nat);

    h_rad = epsilon_fg * sigma * (T_fg + T_amb) * (T_fg^2 + T_amb^2);
    h_out = h_rad + h_conv;

    % Resistances
    R_in  = 1 / (h_in * A_in);
    R_al  = t_al / (k_al * A_al);
    R_cf  = t_cf / (k_cf * A_cf);
    R_fg  = t_fg / (k_fg * A_fg);
    R_out = 1 / (h_out * A_out);

    R_gas_al = R_in + R_al / 2;
    R_al_cf  = R_al / 2 + R_cf / 2;
    R_cf_fg  = R_cf / 2 + R_fg / 2;
    R_fg_amb = R_fg / 2 + R_out;

    % Heat Fluxes [W]
    Qdot_gas_al = (T_gas - T_al) / R_gas_al;
    Qdot_al_cf  = (T_al - T_cf)  / R_al_cf;
    Qdot_cf_fg  = (T_cf - T_fg)  / R_cf_fg;
    Qdot_fg_amb = (T_fg - T_amb) / R_fg_amb;

    % Dynamic Wall Temperature Updates
    dU_al = (Qdot_gas_al - Qdot_al_cf) * dt;
    dU_cf = (Qdot_al_cf  - Qdot_cf_fg) * dt;
    dU_fg = (Qdot_cf_fg  - Qdot_fg_amb) * dt;

    T_al_new = T_al + dU_al / (m_al * cp_al);
    T_cf_new = T_cf + dU_cf / (m_cf * cp_cf);
    T_fg_new = T_fg + dU_fg / (m_fg * cp_fg);

    ThermalState_new.T_layers = [T_al_new, T_cf_new, T_fg_new];

    % Net heat added to gas [J] (negative if gas is hotter than wall and cools)
    Q_gas = -Qdot_gas_al * dt;
end
