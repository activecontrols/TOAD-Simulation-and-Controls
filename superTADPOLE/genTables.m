function genTables(I_sp, cstar, A_t,  g0, CdA_fu, CdA_ox, K_fu, K_ox)

warning ('off','all');

f_thrust = linspace(200, 590, 50);

mdot = f_thrust / I_sp;

P_c = cstar/(A_t * g0) .* mdot;

% table_Pc = [f_thrust;P_c];

% writematrix(table_Pc, "Tables/thrustToPc.csv");

mdot_fu = mdot * 6/11;
mdot_ox = mdot * 5/11;

KGM3_TO_LBIN3 = 3.61273e-5;
PSI_TO_PA = 6894.7;
ipa_temp = 293.17; %%linspace(290,400, 50);
lox_temp =  90.17; %%linspace(90,150,50);
P_inj_ox = P_c * 2.28571428571 - 78.57; %% basic linear interpolation first guess
P_inj_fu = P_c * 1.64204545455 - 47.20; %% basic linear interpolation

delta_rho_ox = 10^10;
delta_rho_fu = 10^10;

rho_ox = FastFluidProps("Oxygen", 'From_P_T', P_inj_ox * PSI_TO_PA, lox_temp).rho * KGM3_TO_LBIN3;
rho_fu = FastFluidProps("IPA", 'From_P_T', P_inj_fu * PSI_TO_PA, ipa_temp).rho * KGM3_TO_LBIN3;

while abs(delta_rho_fu(end)) > 10^-5
    rho_fu_old = rho_fu;
    rho_fu;
    inj_fu = P_c + (mdot_fu/g0).^2./(2*rho_fu*(CdA_fu)^2);

    rho_fu = FastFluidProps("IPA", 'From_P_T', inj_fu * PSI_TO_PA, ipa_temp).rho * KGM3_TO_LBIN3;

    delta_rho_fu = rho_fu_old - rho_fu;
end

while abs(delta_rho_ox(end)) > 10^-5
    rho_ox_old = rho_ox;

    inj_ox = P_c + (mdot_ox/g0).^2./(2*rho_ox*(CdA_ox)^2);

    rho_ox = FastFluidProps("Oxygen", 'From_P_T', inj_ox * PSI_TO_PA, lox_temp).rho * KGM3_TO_LBIN3;
    delta_rho_ox = rho_ox_old - rho_ox;
end


p_inj_ox = P_c + (mdot_ox/g0).^2./(2*rho_ox*(CdA_ox)^2);

max(p_inj_ox)

p_inj_fu = P_c + (mdot_fu/g0).^2./(2*rho_fu*(CdA_fu)^2);

p_valve_fu = p_inj_fu + K_fu * mdot_fu.^2;
p_valve_ox = p_inj_ox + K_ox * mdot_ox.^2;


rho_water = 13/360;
p_tank = 550;

max(p_valve_fu)
min(p_tank - p_valve_fu)

cv_fu = mdot_fu./(rho_fu.* rho_water .* (p_tank - p_valve_fu)).^0.5 * 60 / 231;
cv_ox = mdot_ox./(rho_ox.* rho_water .* (p_tank - p_valve_ox)).^0.5 * 60 / 231;

table_cv_fu = [f_thrust; cv_fu];
table_cv_ox = [f_thrust; cv_ox];

writematrix(table_cv_fu, "Tables/thrustCvFu.csv");
writematrix(table_cv_ox, "Tables/thrustCvOx.csv");
end
