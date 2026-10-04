function extractTables(constantsTADPOLE)
    thrust_cmd = out.thrust_cmd(:);
    ox_angle_cmd = out.angle_ox_cmd(:);
    pc_cmd = out.pc_cmd(:);
    fu_angle_cmd = out.angle_fu_cmd(:);
    pfu_cmd = out.pfu_cmd(:);
    pox_cmd = out.pox_cmd(:);
    mdot_fu = out.mdot_fu(:);
    mdot_ox = out.mdot_ox(:);

    rho_valve_ox = rhoToImperial(constantsTADPOLE.dens_o);
    rho__valve_fu = rhoToImperial(constantsTADPOLE.dens_fu);

    rho_inj_ox = rhoToImperial(constantsTADPOLE.dens_o); %% Will end up being table based, changes based on time
    rho_inj_fu = rhoToImperial(constantsTADPOLE.dens_fu); %% Will be table based, changes based on throttle

    
end

function cv = extract_cv(P_tank, P_valve, mdot, rho)
    
    rho_water = 13/360; % lb/in^3

    cv = mdot/sqrt(rho * rho_water * (P_tank - P_valve));
end

function cstar = extract_cstar(P_c, A_t, mdot)
    cstar = P_c * A_t / mdot;
end

function I_sp = extract_I_sp(F_thrust, mdot, g0)
    %% specific impulse

    I_sp = F_thrust / (mdot * g0);
end

function CdA = extract_CdA(P_inj, P_c, mdot, rho)
    CdA =  mdot/sqrt((P_inj - P_c)*2*rho);

end

function K = extract_K(P_inj, P_valve, mdot)
    K = (P_valve - P_inj) / mdot^2; 
end