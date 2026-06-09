function [injPress, valve_coeff_fu_cmd] = fuelCircuit(chamber_pressure_m, injector_pressure_ox_m, injector_pressure_fu_m, tank_pressure_fu_m, constantsSTADPOLE)
    
    % Notes:
    % - Fuel valve is the follow valve
    % - The function should output valve angle, but the fluid mechanics 
    %   function only accepts valve coefficients, so this functin will 
    %   be updated once the fluid mechanics function is updated.

    % Inputs:
    % - chamber_pressure_m: measured chamber pressure [Pa]
    % - injector_pressure_ox_m: measured oxygen injector pressure [Pa]
    % - tank_pressure_fu_m: measured fuel tank pressure [Pa]
    % - constantsSTADPOLE: constant values
    
    % Output
    % - valve_coeff_ox_cmd: commanded valve coefficient [unitless]

    % Variables
    of_ratio = constantsSTADPOLE.OF_target; % [unitless]
    rho_ox = constantsSTADPOLE.dens_o; % [kg/m^3]
    rho_fu = constantsSTADPOLE.dens_f; % [kg/m^3]
    discharge_coeff_ox = constantsSTADPOLE.d_coeff_ox; % [unitless],  -------------Random value just to get the code to run
    discharge_coeff_fu = constantsSTADPOLE.d_coeff_fu; % [unitless],  -------------Random value just to get the code to run
    orifice_area_ox = constantsSTADPOLE.a_i_o; % [m^2],  -------------Random value just to get the code to run
    orifice_area_fu = constantsSTADPOLE.a_i_f; % [m^2],  -------------Random value just to get the code to run
    
    % 150 psi at max thrust, quarter at 50%
    rho_water = constantsSTADPOLE.dens_w; % [kg/m^3]
    
    % Feedback trim variables
    k_fu = 1E-11; % [unitless], integral gain for fuel trim,  -------------Random value just to get the code to run
    time_step = 0.001; % [s], based on loop time,  -------------Random value just to get the code to run

    persistent integral_error_fu
   
    if isempty(integral_error_fu)
        integral_error_fu = 0;
    end

    % Valve coefficient equations
    delta_p_injector_ox_m = max(injector_pressure_ox_m - chamber_pressure_m, 0); % [Pa]
    term_1 = sqrt(2 * delta_p_injector_ox_m * rho_ox); % [kg/m^2-s]
    massflow_ox_m = discharge_coeff_ox * orifice_area_ox * term_1; % [kg/s]
    massflow_fu_cmd = massflow_ox_m / of_ratio; % [kg/s]
    injector_pressure_fu_cmd = chamber_pressure_m + (massflow_fu_cmd ^ 2) / (2 * rho_fu * (discharge_coeff_fu * orifice_area_fu) ^ 2); % [Pa]
    friction_pressure_drop_fu = massflow_fu_cmd^2 * 2082581.1;
    valve_coeff_fu_cmd = massflow_fu_cmd *  sqrt(max(1 / (rho_fu * rho_water * (tank_pressure_fu_m - friction_pressure_drop_fu - injector_pressure_fu_cmd)), 0)); % [m^3.5/kg^0.5]
    injPress = injector_pressure_fu_cmd;
    
    % Feedback trim equations
    error_fu = injector_pressure_fu_cmd - injector_pressure_fu_m; % [Pa]
    integral_error_fu = integral_error_fu + error_fu * time_step;
    trim_fu = k_fu * integral_error_fu;
    % valve_coeff_fu_cmd = valve_coeff_fu_cmd + trim_fu; % [m^3.5/kg^0.5]

