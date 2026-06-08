function [valve_coeff_ox_cmd,error_ox] = oxygenCircuit(thrust_cmd, chamber_pressure_m, tank_pressure_ox_m, constantsSTADPOLE)
    
    % Notes:
    % - Oxygen valve is the lead valve
    % - The function should output valve angle, but the fluid mechanics 
    %   function only accepts valve coefficients, so this functin will 
    %   be updated once the fluid mechanics function is updated.

    % Inputs:
    % - thrust: commanded thrust [N]
    % - chamber_pressure_m: measured chamber pressure [Pa]
    % - tank_pressure_ox_m: measured oxygen tank pressure [Pa]
    % - constantsSTADPOLE: constant values
    
    % Output
    % - valve_coeff_ox_cmd: commanded valve coefficient [unitless]

    % Variables
    throat_area = constantsSTADPOLE.a_t; % [m^2]
    thrust_coeff = constantsSTADPOLE.c_f; % [unitless]
    char_vel = constantsSTADPOLE.c_star; % [m/s]
    of_ratio = constantsSTADPOLE.OF_target; % [unitless]
    rho_ox = constantsSTADPOLE.dens_o; % [kg/m^3]
    discharge_coeff_ox = constantsSTADPOLE.d_coeff_ox; % [unitless],  -------------Random value just to get the code to run
    orifice_area_ox = constantsSTADPOLE.a_i_o; % [m^2],  -------------Random value just to get the code to run
    rho_water = constantsSTADPOLE.dens_w; % [kg/m^3]
    
    % Feedback trim variables
    k_ox = 1E-11; % [unitless], integral gain for oxygen trim,  -------------Random value just to get the code to run
    time_step = 0.001; % [s], based on loop time,  -------------Random value just to get the code to run

    persistent integral_error_ox
   
    if isempty(integral_error_ox)
        integral_error_ox = 0;
    end

    % Valve coefficient equations
    chamber_pressure_cmd = thrust_cmd / (throat_area * thrust_coeff); % [Pa]
    massflow_tot_cmd = chamber_pressure_cmd * throat_area / char_vel; % [kg/s]
    massflow_ox_cmd = massflow_tot_cmd * of_ratio / (of_ratio + 1); % [kg/s]
    friction_pressure_drop_ox = 96402.2 * massflow_ox_cmd^2;
    injector_pressure_ox_cmd = chamber_pressure_cmd + (massflow_ox_cmd ^ 2) / (2 * rho_ox * (discharge_coeff_ox * orifice_area_ox) ^ 2); % [Pa]
    valve_coeff_ox_cmd = massflow_ox_cmd * sqrt(max(1 / (rho_ox * rho_water * (tank_pressure_ox_m - friction_pressure_drop_ox - injector_pressure_ox_cmd)), 0)); % [m^3.5/kg^0.5]
    
    % Feedback trim equations
    error_ox = chamber_pressure_cmd - chamber_pressure_m; % [Pa]
    integral_error_ox = integral_error_ox + error_ox * time_step;
    % this is very scuffed, needs to be fixed
    trim_ox = (k_ox * integral_error_ox)* massflow_ox_cmd/max(massflow_tot_cmd-massflow_ox_cmd, 1);
    valve_coeff_ox_cmd = valve_coeff_ox_cmd + trim_ox; % [m^3.5/kg^0.5]

