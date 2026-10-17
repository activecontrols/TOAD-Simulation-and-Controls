function extractTables(constantsSTADPOLE)
    att2 = readmatrix("Hotfire Data (from Throttle Box DAQ) - att5.csv");

    help = att2(1:end, 4) * 360;
    startup = find(att2(1:end,2) > 1, 1); % bc tadpole is scuffed find((att2(1:end, 4)*360)> 20, 1)+1
    time = att2(startup:end,1);

    
    % sampled at 200 hz
    [b,a] = butter(2, 0.5);
    
    thrust_cmd = filtfilt(b,a, att2(startup:end,3));
    
    ox_commanded =  filtfilt(b,a, att2(startup:end, 4)*360);
    fu_cmd = filtfilt(b,a, att2(startup:end, 5)*360);
    fu_pos = filtfilt(b,a, att2(startup:end, 6)*360);
    fu_valve_down = filtfilt(b,a, att2(startup:end, 21));
    fu_valve_up = filtfilt(b,a, att2(startup:end, 20));
    ox_valve_down =  filtfilt(b,a, att2(startup:end, 16));
    ox_valve_up =  filtfilt(b,a, att2(startup:end, 15));

    p_chamber = filtfilt(b,a, att2(startup:end, 14));
    mdot_ox =  filtfilt(b,a, att2(startup:end, 29));
    mdot_fu = filtfilt(b,a, att2(startup:end, 30));

    p_tank = 550;
    a_t = cdaToImperial(constantsSTADPOLE.a_t);

    g_0 = cstarToImperial(constantsSTADPOLE.g);

    o_dens_valve = rhoToImperial(constantsSTADPOLE.dens_o); % constant
    fu_dens_valve = rhoToImperial(constantsSTADPOLE.dens_f); % constant

    o_dens_inj = rhoToImperial(constantsSTADPOLE.dens_o); % varying w.r.t. time
    fu_dens_inj = rhoToImperial(constantsSTADPOLE.dens_f); % varying w.r.t. throttle 
    % size(ox_valve_up)
    % size(ox_valve_down)
    % size(mdot_ox)
    cvs_ox = extract_cv(ox_valve_up, ox_valve_down, mdot_ox, o_dens_valve);
    cvs_fu = extract_cv(fu_valve_up, fu_valve_down, mdot_fu, fu_dens_valve);
    find(abs(cvs_fu) > 10, 50);
    mdot_fu(find(abs(cvs_fu) > 5, 50));
    fu_valve_up(find(abs(cvs_fu) > 5, 50))-fu_valve_down(find(abs(cvs_fu) > 5, 50));
   
    
    % disp(cvs_ox - arrayfun(@(x) AngleToCoefficient(x),ox_commanded))
    figure;
    hold on;
    scatter(ox_commanded, real(cvs_ox));
    scatter(ox_commanded, arrayfun(@(x) AngleToCoefficient(x),ox_commanded));
    % time, arrayfun(@(x) AngleToCoefficient(x),ox_commanded)
    legend("calculated ox cv", "datasheet");

    hold off;
    figure;
    hold on;
    scatter(fu_pos, real(cvs_fu));
    scatter(fu_pos, arrayfun(@(x) AngleToCoefficient(x),fu_pos));
    legend("calculated fu cv", "datasheet");

    figure;
    cstar = extract_cstar(p_chamber, a_t, mdot_ox + mdot_fu, g_0);
    scatter(p_chamber, cstar);
    title("cstar wrt. chamber pressure");

    figure;
    cstar = extract_cstar(p_chamber, a_t, mdot_ox + mdot_fu, g_0);
    scatter(time, cstar);
    %vplot(time, mdot_fu)
    title("cstar wrt. time");

    figure;
    CdA = extract_CdA(ox_valve_down, p_chamber, mdot_ox, o_dens_inj);
    scatter(p_chamber, CdA);
    title("CdA wrt. chamber pressure");
    find((CdA) < 0, 50);
    mdot_ox(find((CdA) < 0, 50));

    figure;
    I_sp = extract_I_sp(thrust_cmd, mdot_ox + mdot_fu);
    scatter(thrust_cmd, I_sp);
    title("I_{sp} wrt thrust");
    find((I_sp) > 10, 50);
    mdot_ox(find((I_sp) > 10, 50));
end

function cv = extract_cv(P_tank, P_valve, mdot, rho)
    
    rho_water = 13/360; % lb/in^3
    cv = mdot./(rho * rho_water * (P_tank - P_valve)).^0.5 * 60 / 231;
end

function cstar = extract_cstar(P_c, A_t, mdot, g0)
    cstar = P_c * A_t ./ mdot .* g0;
end

function I_sp = extract_I_sp(F_thrust, mdot)
    %% specific impulse

    I_sp = F_thrust ./mdot;
end

function CdA = extract_CdA(P_inj, P_c, mdot, rho)
    CdA =  mdot./((P_inj - P_c)*2*rho).^0.5;
end

function K = extract_K(P_inj, P_valve, mdot)
    K = (P_valve - P_inj) ./ mdot.^2; 
end

function psi = PaToImperial(Pa)
psi = Pa * 0.0001450377;
end

function Pa = psiToMetric(psi)
Pa = psi / 0.0001450377;
end

% rho is kg/m^3 to lb/in^3
function imprho = rhoToImperial(rho)
imprho = rho * 13/360000;
end

% cstar is m/s to ft/s
function impcstar = cstarToImperial(cstar)
impcstar = cstar * 3.28084;
end

% cda is m^2 to in^2
function impcda = cdaToImperial(cda)
impcda = cda * 39.37007874^2;
end

function impcv = cvToImperial(cv)
impcv = cv / 15850.3 * sqrt(1.450377e-4);
end

function impmdot = mdotToImperial(mdot)
impmdot = mdot * 2.204623;
end

function mdot = mdotToMetric(impmdot)
mdot = impmdot * 2.204623;
end
