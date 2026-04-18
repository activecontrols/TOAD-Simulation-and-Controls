function valve_coeff_ox_cmd = oxygenCircuit(thrust, chamber_pressure_m, tank_pressure_ox_m)
    
    % Notes:
    % - Oxygen valve is the lead valve
    % - The function should output valve angle, but the fluid mechanics 
    %   function only accepts valve coefficients, so this functin will 
    %   be updated once the fluid mechanics function is updated.

    % Inputs:
    % - thrust: commanded thrust [N]
    % - chamber_pressure_m: measured chamber pressure [Pa]
    % - tank_pressure_ox_m: measured oxygen tank pressure [Pa]
    
    % Output
    % - valve_coeff_ox_cmd: commanded valve coefficient [unitless]

    % Variables
    throat_area = constantsSTADPOLE.a_t; % [m^2]
    thrust_coeff = constantsSTADPOLE.c_f; % [unitless]
    char_vel = constants.c_star; % [m/s]
    of_ratio = constants.OF_target; % [unitless]
    rho_ox = ; % [kg/m^3]
    discharge_coeff_ox = ; % [unitless]
    orifice_area_ox = ; % [m^2]
    friction_pressure_drop_ox = ; % [Pa]
    rho_water = ; % [kg/m^3]
    
    % Feedback trim variables
    k_ox = ; % [unitless], integral gain for oxygen trim
    time_step = ; % [s], based on loop time

    persistent integral_error_ox
   
    if isempty(integral_error_ox)
        integral_error_ox = 0;
    end

    % Valve coefficient equations
    chamber_pressure_cmd = thrust / (throat_area * thrust_coeff); % [Pa]
    massflow_tot_cmd = chamber_pressure_cmd * throat_area / char_vel; % [kg/s]
    massflow_ox_cmd = massflow_tot_cmd * of_ratio / (of_ratio + 1); % [kg/s]
    term_1 = (massflow_ox_cmd ^ 2) / (2 * rho_ox * (discharge_coeff_ox * orifice_area_ox) ^ 2); % [Pa]
    injector_pressure_ox_cmd = chamber_pressure_cmd + term_1; % [Pa]
    term_2 = sqrt(1 / (rho_ox * rho_water * (tank_pressure_ox_m - friction_pressure_drop_ox - injector_pressure_ox_cmd))); % [m^3.5s/kg^1.5]
    valve_coeff_ox_cmd = massflow_ox_cmd * term_2; % [m^3.5/kg^0.5]
    
    % Feedback trim equations
    error_ox = chamber_pressure_cmd - chamber_pressure_m; % [Pa]
    integral_error_ox = integral_error_ox + error_ox * time_step;
    trim_ox = k_ox * integral_error_ox;
    valve_coeff_ox_cmd = valve_coeff_ox_cmd + trim_ox; % [m^3.5/kg^0.5]

