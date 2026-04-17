function valve_coeff_ox_cmd = oxygenCircuit(thrust, chamber_pressure_m, tank_pressure_ox_m, t, constants)
    persistent lastT trimOx
    if isempty(lastT)
        lastT = 0;
        trimOx = 0;
    
    end
    % if you have a better way to do this please do
    dT = t - lastT;
    lastT = t;
    
    % thrust [N]
    % chamber_pressure_m [Pa]
    throat_area = constantsSTADPOLE.a_t; % [m^2]
    thrust_coeff = constantsSTADPOLE.c_f; % [unitless]
    char_vel = constants.c_star; % [m/s]
    of_ratio = constants.OF_target; % [unitless]
    rho_ox = ; % [kg/m^3]
    discharge_coeff_ox = []; % [unitless]
    orifice_area_ox = []; % [m^2]
    friction_pressure_drop_ox = ; % [Pa]
    rho_water = ; % [kg/m^3]
    k_ox = ; % [unitless]

    chamber_pressure_cmd = thrust / (throat_area * thrust_coeff); % [Pa]
    massflow_tot_cmd = chamber_pressure_cmd * throat_area / char_vel; % [kg/s]
    massflow_ox_cmd = massflow_tot_cmd * of_ratio / (of_ratio + 1); % [kg/s]
    term_1 = (massflow_ox_cmd ^ 2) / (2 * rho_ox * (discharge_coeff_ox * orifice_area_ox) ^ 2); % [Pa]
    injector_pressure_ox_cmd = chamber_pressure_cmd + term_1; % [Pa]

    term_2 = sqrt(1 / (rho_ox * rho_water * (tank_pressure_ox_m - friction_pressure_drop_ox - injector_pressure_ox_cmd))); % [unknown]
    valve_coeff_ox_cmd = massflow_ox_cmd * term_2; % [unitless?]
    
    % Add integral trim
    % unsure how to do integral (based on time)
    
    error_ox = chamber_pressure_cmd - chamber_pressure_m;
    trim_ox = trimOx + k_ox * (error_ox) * dT; % less wrong probably??
    trimOx = trim_ox;
    valve_coeff_ox_cmd = valve_coeff_ox_cmd + trim_ox; % [unitless?]

    % Technically, the function should output valve angle, but idk what the
    % function is since its based on valve coefficient. also the fluid
    % mechanics function only take sin valve coeff rn. 