function [Node_new, mdot_nozzle, isLit] = SKIPPER_Model(Node, Inflow, Nozzle, dt, IGN)

%% SKIPPER_Model uses Quasi-Steady LRE equations to model the state response of the engine 
% Has a very approximate 1st order delay on the chamber pressure response.
% Can handle both acting as a standard node for coldflows or ignited if IGN
% == true.
% Inputs:
%   Node    - Struct holding current chamber state (.P, .T, .V, .m, .U, .u, etc.)
%   Inflow  - Struct with incoming flow rates (.mdot_ox, .mdot_fu, .mdot_n2)
%   Nozzle  - Struct with nozzle parameters (.A_throat, .Cd, .P_back)
%   dt      - Discrete simulation time step [s]
%   IGN     - Boolean ignition command / spark state (true/false)
%
% Outputs:
%   Node_new    - Updated chamber node state struct
%   mdot_nozzle - Mass flow rate leaving nozzle [kg/s]
%   isLit       - Boolean indicating active combustion status

    %% Combustion Properties (feel free to make them state dependent by calling specific UDFs)
    cstar_hot = 1419.15; % [m/s]
    Tc_hot    = 2900.0; % [K]
    gamma_hot = 1.21;   % Specific heat ratio
    R_hot     = 370.0;  % [J/(kg*K)]
    tau       = 0.05;   % [s]

   
    %% Sum incoming massflows
    mdot_ox = max(0.0, Inflow.mdot_ox);
    mdot_fu = max(0.0, Inflow.mdot_fu);
    mdot_prop = mdot_ox + mdot_fu;

    
    A_t = 1.2810e-3;
    Cd = 0.95;
    P_back = 101325.0;

    % =========================================================================
    % 4. FIRST-ORDER LAG DYNAMICS & STATE UPDATE
    % =========================================================================
    V_ch = 0.001302;
    if isfield(Node, 'V') && Node.V > 0, V_ch = Node.V; end
    P_curr = P_back;
    if isfield(Node, 'P') && Node.P > 0, P_curr = Node.P; end

    if IGN
        isLit  = true;
        cstar  = cstar_hot;
        Tc     = Tc_hot;
        gamma  = gamma_hot;
        R_gas  = R_hot;
        cv_gas = R_gas / (gamma - 1.0);
        T_new  = Tc;

        % Analytical rocket chamber pressure: P_c = (mdot_prop * cstar) / (A_t * Cd)
        Pc_ss = max(P_back, (mdot_prop * cstar) / (A_t * Cd));
        decay = exp(-dt / tau);
        P_new = P_curr * decay + Pc_ss * (1.0 - decay);
    else
        isLit  = false;
        cstar  = 1100.0;
        gamma  = 1.40;
        R_gas  = 296.8; % [J/(kg*K)]
        cv_gas = R_gas / (gamma - 1.0);
        T_new  = 293.15;

        decay = exp(-dt / max(tau, 0.02));
        P_new = P_curr * decay + P_back * (1.0 - decay);
    end
    P_new = max(P_new, P_back);

    % =========================================================================
    % 6. CHOKED & UNCHOKED NOZZLE MASS FLOW
    % =========================================================================
    PR = P_back / P_new;
    PR_crit = (2.0 / (gamma + 1.0))^(gamma / (gamma - 1.0));

    if P_new <= P_back
        mdot_nozzle = 0.0;
    elseif PR <= PR_crit
        % Choked flow: sonic velocity at throat
        mdot_nozzle = (P_new * A_t * Cd) / cstar;
    else
        % Subsonic / unchoked compressible orifice flow
        mdot_nozzle = Cd * A_t * sqrt(2.0 * P_new * (P_new / (R_gas * T_new)) * ...
                      (gamma / (gamma - 1.0)) * max(0.0, PR^(2.0 / gamma) - PR^((gamma + 1.0) / gamma)));
    end

    % =========================================================================
    % 7. THERMODYNAMIC STATE CLOSING
    % =========================================================================
    rho_new = P_new / (R_gas * T_new);
    m_new   = rho_new * V_ch;
    u_new   = cv_gas * T_new;
    U_new   = m_new * u_new;
    h_new   = u_new + P_new / rho_new;

    Node_new       = Node;
    Node_new.V     = V_ch;
    Node_new.m     = m_new;
    Node_new.U     = U_new;
    Node_new.u     = u_new;
    Node_new.rho   = rho_new;
    Node_new.T     = T_new;
    Node_new.P     = P_new;
    Node_new.h     = h_new;
    Node_new.gamma = gamma;
    Node_new.cstar = cstar;
    Node_new.isLit = isLit;
    if isfield(Node, 'm_liq'), Node_new.m_liq = 0.0; end
    if isfield(Node, 'm_gas'), Node_new.m_gas = m_new; end
end
