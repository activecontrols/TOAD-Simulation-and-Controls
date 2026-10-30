function [Q_ull, Q_liq, ThermalState_new] = TankThermalModel(TankState, ThermalState, dt, T_amb)
% TANKTHERMALMODEL Models convective & radiative heat inleak into a propellant tank.
%
% Inputs:
%   TankState        - Struct containing .Ullage and .Liquid sub-nodes
%   ThermalState     - Struct with field .T_wall
%   dt               - Time step [s]
%   T_amb            - Ambient temperature [K] (default 293.15)
%
% Outputs:
%   Q_ull            - Heat added to ullage gas [J]
%   Q_liq            - Heat added to liquid propellant [J]
%   ThermalState_new - Updated thermal state struct

    if nargin < 4 || isempty(T_amb), T_amb = 293.15; end

    T_wall = ThermalState.T_wall;
    T_ull  = TankState.Ullage.T;
    T_liq  = TankState.Liquid.T;

    % Tank Geometry & Wall Specs (Aluminum tank shell)
    A_total = 0.65;      % Total surface area [m^2]
    V_total = TankState.V;
    V_liq   = TankState.Liquid.V;
    liqFrac = max(0.0, min(1.0, V_liq / V_total));

    A_liq = A_total * liqFrac;
    A_ull = A_total * (1 - liqFrac);

    m_wall  = 4.5;       % Wall mass [kg]
    cp_wall = 900;       % J/(kg-K) Aluminum

    % Outside Convection + Radiation to Wall
    h_ext = 15.0;        % W/(m^2-K) natural convection to ambient
    epsilon_wall = 0.2;  % Bare/polished aluminum emissivity
    sigma = 5.670374e-8;
    h_rad = epsilon_wall * sigma * (T_wall + T_amb) * (T_wall^2 + T_amb^2);
    Qdot_ext_in = (h_ext + h_rad) * A_total * (T_amb - T_wall);

    % Internal Convection from Wall to Fluids
    h_int_liq = 300.0;   % High convection / pool boiling for liquid
    h_int_ull = 25.0;    % Gas natural convection

    Qdot_wall_to_liq = h_int_liq * A_liq * (T_wall - T_liq);
    Qdot_wall_to_ull = h_int_ull * A_ull * (T_wall - T_ull);

    % Wall Thermal Energy Balance
    dU_wall = (Qdot_ext_in - Qdot_wall_to_liq - Qdot_wall_to_ull) * dt;
    T_wall_new = T_wall + dU_wall / (m_wall * cp_wall);

    ThermalState_new.T_wall = T_wall_new;
    Q_ull = Qdot_wall_to_ull * dt;
    Q_liq = Qdot_wall_to_liq * dt;
end
