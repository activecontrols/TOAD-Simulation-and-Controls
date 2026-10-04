function extractTables

end

function cv = extract_cv(P_up, P_down, mdot, rho)
    
    rho_water = 13/360; % lb/in^3

    cv = mdot/sqrt(rho * rho_water * (P_up - P_down));
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

function K = extract_K(P_inj, P_down, mdot)
    K = (P_down - P_inj) / mdot^2; 
end