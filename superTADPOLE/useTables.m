function rho_ox = injOxDens(time)

end

function rho_fu = injFuDens(throttle)
    throttle = [];
    rho_list = [];
end

function cv = angleToCv(angle)
    cv_list = [0, 0.01, 0.06, 0.21, 0.46, 0.83, 1.35, 2.04, 2.93];
    angle_list = [90*0.2, 90*0.3, 90*0.4, 90*0.5, 90*0.6, 90*0.7, 90*0.8, 90*0.9, 90*1];
    if angle >= 90
        cv = 2.93/(1.316e6);
    elseif angle <= 18 || angle == 1/0
        cv = 0;
    else
        cv = interp1(angle_list, cv_list, angle);
    end
end

function angle = cvToAngle(cv)
    cv_list = [0, 0.01, 0.06, 0.21, 0.46, 0.83, 1.35, 2.04, 2.93];
    angle_list = [90*0.2, 90*0.3, 90*0.4, 90*0.5, 90*0.6, 90*0.7, 90*0.8, 90*0.9, 90*1];
    if cv >= 2.93
        angle = 90;
    elseif cv < 0
        angle = 0;
    else
        angle = interp1(cv_list, angle_list, cv);
    end
end

function Cv = commandedCv(p_tank, p_valve, mdot, rho)
    rho_water = 13/360; % lb/in^3
    Cv = mdot/sqrt(rho * rho_water * (p_tank - p_valve));
end

function p_valve = psiValve(p_inj, mdot, K)
    p_valve = p_inj + K * mdot^2;
end

function p_inj = psiInj(p_chamber, mdot, rho, CdA)
    p_inj = p_chamber + mdot^2 / (2 * rho * CdA^2);
end

function p_chamber = psiChamber(Isp, thrust, cstar , constantsTADPOLE)
    mdot = thrust / (Isp * constantsTADPOLE.g);
    p_chamber = cstar * mdot / constantsTADPOLE.a_t;
end

