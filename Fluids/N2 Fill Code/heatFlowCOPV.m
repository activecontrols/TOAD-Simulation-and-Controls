function [Q_gas, T_al_new, T_cf_new, T_fg_new, h_out, Ra_L, deltaT] = ...
    heatFlowCOPV(fillMode, timeStep, T_gas, T_al, T_cf, T_fg)

k_N2 = double(py.CoolProp.CoolProp.PropsSI('L', ...
    'P', P_gas, 'T', T_al, 'Nitrogen'));

rho_N2 = double(py.CoolProp.CoolProp.PropsSI('D', ...
    'P', P_gas, 'T', T_al, 'Nitrogen'));

cp_N2 = double(py.CoolProp.CoolProp.PropsSI('Cpmass', ...
    'P', Pcopv, 'T', Tfilm, 'Nitrogen'));




% Dynamic COPV wall heat-transfer model
% Positive Q_gas = energy leaving nitrogen and entering wall
deltaT = abs(T_fg - T_amb);
L_copv = .57; %m
alpha_N2 = k_N2 / (rho_N2 * cp_N2);
D_copv

%% Ambient
T_amb = 297.04; % K

%% Inside convection
A_in = 1.1289;
h_in_fill = 150; % W / (m^2 * K)

%calculate Raleigh number
Ra = (g * beta_N2 * cp_N2 * rho_N2^2 * Lc^3 * ...
    abs(Tcopv - Tliner)) / (mu_N2 * k_N2);

Nu_nat = 0.104 * Ra^0.352;

h_in_pause = Nu_nat * k_N2 / Lc;

if fillMode == 0
    h_in = h_in_pause;
else
    h_in = h_in_fill;
end

%% Aluminum liner
A_al = 1.1468;
t_al = 0.002;
k_al = 167;

rho_al = 2700;      % kg/m^3
cp_al  = 900;       % J/(kg*K)

m_al = rho_al * A_al * t_al;

%% Carbon fiber
A_cf = 1.2040;
t_cf = 0.00423;
k_cf = 0.5;

rho_cf = 1600;      % approximate
cp_cf  = 900;       % approximate

m_cf = rho_cf * A_cf * t_cf;

%% Fiberglass
A_fg = 1.2480;
t_fg = 0.0005;
k_fg = 0.3;

rho_fg = 1900;      % approximate
cp_fg  = 800;       % approximate

m_fg = rho_fg * A_fg * t_fg;

%% Outside convection
A_out = 1.2529; %m^2
k_air = 0.024; %W / (m * K)
Ra_L = 10^8 * (deltaT) * L_copv^3;
Pr = 0.71;

%h_out using wiki: vertical plane external flow.
h_out = (k_air/L_copv) * (0.68 + (0.67 * Ra_L^(1/4)) / ...
    (1 + (0.492/Pr)^(9/16))^(4/9));

%% Full-layer resistances
R_in  = 1 / (h_in * A_in);

R_al = t_al / (k_al * A_al);
R_cf = t_cf / (k_cf * A_cf);
R_fg = t_fg / (k_fg * A_fg);

R_out = 1 / (h_out * A_out);

%% Resistance between thermal nodes
% Each wall temperature represents approximately the CENTER of its layer

R_gas_al = R_in + R_al/2;

R_al_cf = R_al/2 + R_cf/2;

R_cf_fg = R_cf/2 + R_fg/2;

R_fg_amb = R_fg/2 + R_out;

%% Heat-transfer rates
Qdot_gas_al = (T_gas - T_al) / R_gas_al;

Qdot_al_cf = (T_al - T_cf) / R_al_cf;

Qdot_cf_fg = (T_cf - T_fg) / R_cf_fg;

Qdot_fg_amb = (T_fg - T_amb) / R_fg_amb;

%% Energy balance on each wall layer

dU_al = (Qdot_gas_al - Qdot_al_cf) * timeStep;

dU_cf = (Qdot_al_cf - Qdot_cf_fg) * timeStep;

dU_fg = (Qdot_cf_fg - Qdot_fg_amb) * timeStep;

%% Update layer temperatures

T_al_new = T_al + dU_al / (m_al * cp_al);

T_cf_new = T_cf + dU_cf / (m_cf * cp_cf);

T_fg_new = T_fg + dU_fg / (m_fg * cp_fg);

%% Energy removed from nitrogen
Q_gas = Qdot_gas_al * timeStep;


end