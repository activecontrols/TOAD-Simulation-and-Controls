function [Q_ull, Q_liq, h_ull, h_int] = TankThermalModel_v2(TankNode, mdot_n2_in, dt)
%% Tank Thermals
% Uses forced convection w/ diffuser params from CDR. Mixes in nayural
% convection, heat transfer into constant temp 90K wall in the LOX tank
% (hella conservative). returns out the thermal energy transfered and the
% convection coefs.

    if nargin < 2 || isempty(mdot_n2_in), mdot_n2_in = 0.0; end
    if nargin < 3 || isempty(dt), dt = 0.010; end

    % Extract States
    T_ull = TankNode.Ullage.T;
    rho_ull = TankNode.Ullage.rho;

    T_liq = TankNode.Liquid.T;
    V_liq = TankNode.Liquid.V;
    V_total = TankNode.V;

    % Geometry Specifications
    D_tank = TankNode.D;
    A_int = (pi / 4.0) * (D_tank^2);           % Liquid surface interface area [m^2]

    A_total = 2.0 * A_int + (4.0 * V_total / D_tank);

    liqFrac = max(0.0, min(1.0, V_liq / V_total));
    A_dry   = A_total * (1.0 - liqFrac);       % Surface area in contact with ullage gas

    % Fixed Wall Temperature
    isLOX = contains(lower(TankNode.Fluid), 'oxygen') || contains(lower(TankNode.Name), 'o2');
    if isLOX
        T_wall = 90.17;                        % Constant steady-state LOX temperature [K]
    else
        T_wall = 293.15;                       % Constant ambient temperature for fuel [K]
    end

    % Nitrogen Transport Properties at Mean Film Temperature
    T_film = max(60.0, (T_ull + T_wall) * 0.5);

    % Dense Nitrogen transport calibrations
    mu_gas = 1.78e-5 * ((T_film / 293.15)^0.75);
    k_gas  = 0.0250 * ((T_film / 293.15)^0.85);
    cp_gas = 1040.0;
    Pr_gas = (mu_gas * cp_gas) / max(k_gas, 1e-4);

    % Forced Convection from Diffuser Jet (Length scale = D_diff)
    v_jet     = max(0.0, mdot_n2_in) / max(rho_ull * A_diff, 1e-6);
    Re_d      = (rho_ull * v_jet * D_diff) / max(mu_gas, 1e-7);
    Nu_forced = 0.023 * (Re_d^0.8) * (Pr_gas^0.4);
    h_forced  = (k_gas / D_diff) * Nu_forced;

    % Natural Convection Floor (Length scale = D_tank)
    dT_wall  = max(abs(T_ull - T_wall), 0.1);
    beta_gas = 1.0 / T_film;
    g        = 9.81;
    Ra       = (g * beta_gas * dT_wall * (D_tank^3) * (rho_ull^2) * cp_gas) / max(mu_gas * k_gas, 1e-9);
    Nu_nat   = 0.10 * (max(Ra, 0.0)^(1.0 / 3.0));
    h_nat    = (k_gas / D_tank) * Nu_nat;

    % Churchill-Usagi blend in dimensional space
    h_ull = (h_forced^3.0 + h_nat^3.0)^(1.0 / 3.0);

    % Heat Transfer into Tank Wall
    Qdot_wall = h_ull * A_dry * (T_wall - T_ull); 

    % Liquid-Gas Interface Heat Transfer
    % Gas recirculation sweeps over the interface; scaled by diffuser/tank diameter ratio
    h_int = h_ull * sqrt(D_diff / D_tank);
    dT_int = T_liq - T_ull;
    Qdot_int = h_int * A_int * dT_int;            

    % Integrated Energy Terms [J]
    Q_ull = (Qdot_wall + Qdot_int) * dt;
    Q_liq = (-Qdot_int) * dt;
end
