function linear_system = fluidMechanics(t,y,u, constantsSTADPOLE)
    
    % State vector y
    
    % y(1) = mass flow oxygen [kg/s]
    % y(2) = mass flow fuel [kg/s]
    % y(3) = chamber pressure [Pa]
    % y(4) = oxygen pressure [Pa]
    % y(5) = fuel pressure [Pa]
    
    % Derivative of state vector y'

    % y(1)' = time rate of mass flow oxygen [kg/s^2]
    % y(2)' = time rate of mass flow fuel [kg/s^2]
    % y(3)' = time rate of chamber pressure [Pa/s]
    % y(4)' = time rate of oxygen pressure [Pa/s]
    % y(5)' = time rate of fuel pressure [Pa/s]


    a_t = constantsSTADPOLE.a_t ; % [m^2], throat area
    
    dens_o = constantsSTADPOLE.dens_o; % [kg/m^3], liquid oxygen density, at 90 K and 1 atm
    dens_f = constantsSTADPOLE.dens_f ; % [kg/m^3], fuel density, at 293 K and 1 atm
    dens_w = constantsSTADPOLE.dens_w ; % [kg/m^3], water density
    g = constantsSTADPOLE.g ; % [m/s^2], acceleration due to gravity
    a_i_o = constantsSTADPOLE.a_i_o; % [m^2] area injector ox
    a_i_f = constantsSTADPOLE.a_i_f ; % [m^2] area injector fuel
    r = constantsSTADPOLE.r ; % [J/kg-K], specific gas constant of exhaust (molecular weight 18.8 g/mol)
    c_star = constantsSTADPOLE.c_star ; % [m/s] characteristic velocity
    of = constantsSTADPOLE.OF_target;

    c_fric = 0.03; % [unitless], line friction coefficient
    c_d_o = 0.445; % [unitless], oxygen injector coefficient
    c_d_f = 0.7; % [unitless], fuel injector coefficient
    l_eq_o = 7 * 0.0254; % [m], oxygen equivalent line length 
    l_eq_f = 7 * 0.0254; % [m], fuel equivalent line length 
    l_o = 1200 * 0.0254; % [m], oxygen line length
    l_f = 600 * 0.0254; % [m], fuel line length
    d = (0.5 - 2 * 0.065) * 0.0254; % [m], line diameter
    
   
 
    v_c = 1.103170594 * 1e-3; % [m^3], chamber volume
    a_l = pi * (d / 2) ^ 2; % [m^2], line area
  
    v_l_o = a_l * l_eq_o;
    v_l_f = a_l * l_eq_f;
    t_o = 90; 
    t_f = 293;

    % Assume constant?
    c_v_o = u(1,1); %1.231e-06; % [m^3.5/kg^0.5], oxygen valve coefficient
    c_v_f = u(2,1); %1.231e-06; % [], fuel valve coefficient
    p_tank_o = 4.137e+6; % [Pa]
    p_tank_f = 4.137e+6; % [Pa]
    t_c = 3570; % [K], chamber temperature
    

    del_p_i_o = 1 / (2 * dens_o * (c_d_o * a_i_o) ^ 2); 
    del_p_l_o = c_fric * l_eq_o / (2 * d * dens_o * a_l ^ 2);
    del_p_i_f = 1 / (2 * dens_f * (c_d_f * a_i_f) ^ 2); 
    del_p_l_f = c_fric * l_eq_f / (2 * d * dens_f * a_l ^ 2);
    
    eqn_1 = (y(4) - y(3) - (del_p_i_o + del_p_l_o) * (y(1) ^ 2)) * a_l / l_o;
    eqn_2 = (y(5) - y(3) - (del_p_i_f + del_p_l_f) * (y(2) ^ 2)) * a_l / l_f;

    eqn_3 = (r * t_c / v_c) * (y(1) + y(2) - a_t * y(3) / c_star);
    eqn_4 = (r * t_o / v_l_o) * (c_v_o * sqrt(dens_o * dens_w * (p_tank_o - y(4))) - y(1));
    eqn_5 = (r * t_f / v_l_f) * (c_v_f * sqrt(dens_f * dens_w * (p_tank_f - y(5))) - y(2));
    
    linear_system = [eqn_1; eqn_2; eqn_3; eqn_4; eqn_5];
end