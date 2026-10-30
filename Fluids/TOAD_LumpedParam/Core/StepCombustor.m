function [Node_new, mdot_nozzle, Q_rxn, isFired] = StepCombustor(Node, Inflow, Nozzle, dt, SparkActive, TorchActive, Q_wall)
% STEPCOMBUSTOR Generic, modular combustion chamber state stepper.
% Evaluates thermochemical state, reaction heat release, and sonic nozzle outflow.
% Supports custom combustor tables or default CEA tables for LOX/IPA.
%
% Inputs:
%   Node        - Struct holding current chamber state (.m, .U, .P, .T, .V, .isLit)
%   Inflow      - Struct with .mdot_ox, .mdot_fu, .mdot_n2, .h_ox, .h_fu, .h_n2
%   Nozzle      - Struct with .A_throat, .Cd, .P_back
%   dt          - Time step [s]
%   SparkActive - Boolean (true if spark plug energized)
%   TorchActive - Boolean (true if hot torch flame entering from igniter)
%   Q_wall      - Thermal heat loss to chamber wall [J]
%
% Outputs:
%   Node_new    - Updated chamber state struct
%   mdot_nozzle - Mass flow through nozzle [kg/s]
%   Q_rxn       - Chemical heat release energy in this step [J]
%   isFired     - Boolean indicating active combustion

    if nargin < 7 || isempty(Q_wall), Q_wall = 0.0; end
    if nargin < 6 || isempty(TorchActive), TorchActive = false; end
    if nargin < 5 || isempty(SparkActive), SparkActive = false; end

    % Standard / Default CEA Baseline Data for LOX/IPA
    if isfield(Node, 'CEA_OF') && ~isempty(Node.CEA_OF)
        CEA_OF    = Node.CEA_OF;
        CEA_Temp  = Node.CEA_Temp;
        CEA_Cstar = Node.CEA_Cstar;
    else
        CEA_OF = [0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, ...
                  1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, ...
                  2.8, 2.9, 3.0, 3.1, 3.2, 3.3];
        CEA_Temp = [1562, 1879, 2172, 2436, 2668, 2862, 3017, 3135, 3219, 3277, ...
                    3315, 3338, 3352, 3359, 3361, 3360, 3356, 3349, 3341, 3332, ...
                    3322, 3311, 3300, 3288, 3275, 3262];
        CEA_Cstar = [1350, 1460, 1540, 1600, 1640, 1670, 1690, 1700, 1705, 1700, ...
                     1690, 1680, 1665, 1650, 1630, 1610, 1590, 1570, 1550, 1530, ...
                     1510, 1490, 1470, 1450, 1430, 1410];
    end

    mdot_ox = max(Inflow.mdot_ox, 0.0);
    mdot_fu = max(Inflow.mdot_fu, 0.0);
    mdot_n2 = 0.0;
    if isfield(Inflow, 'mdot_n2'), mdot_n2 = max(Inflow.mdot_n2, 0.0); end
    if isfield(Inflow, 'mdot_inert'), mdot_n2 = mdot_n2 + max(Inflow.mdot_inert, 0.0); end

    mdot_prop = mdot_ox + mdot_fu;
    mdot_tot  = mdot_prop + mdot_n2;

    % Flammability and Ignition
    isLit = false;
    if isfield(Node, 'isLit'), isLit = Node.isLit; end
    
    OF_in = mdot_ox / max(mdot_fu, 1e-6);
    isFlammable = (mdot_ox > 0.005) && (mdot_fu > 0.003) && ...
                  (OF_in >= 0.5) && (OF_in <= 4.0);

    if isFlammable && (SparkActive || TorchActive || isLit)
        isLit = true;
    elseif mdot_prop < 0.005 || ~isFlammable
        isLit = false;
    end
    isFired = isLit;

    % Gas Properties
    if isLit
        OF_clamped = max(CEA_OF(1), min(CEA_OF(end), OF_in));
        OF_step = CEA_OF(2) - CEA_OF(1);
        idx_f = (OF_clamped - CEA_OF(1)) / OF_step + 1.0;
        i_low = max(1, min(length(CEA_OF)-1, floor(idx_f)));
        f_hi  = idx_f - i_low;
        T_flame = CEA_Temp(i_low) * (1.0 - f_hi) + CEA_Temp(i_low + 1) * f_hi;
        cstar   = CEA_Cstar(i_low) * (1.0 - f_hi) + CEA_Cstar(i_low + 1) * f_hi;
        gamma   = 1.21;
        if isfield(Node, 'gamma') && Node.gamma > 1.05, gamma = Node.gamma; end
        R_gas   = 370.0; % J/(kg-K)
        cp_gas  = gamma * R_gas / (gamma - 1);
        cv_gas  = R_gas / (gamma - 1);

        % Wall heat loss correction
        deltaT_wall = Q_wall / (max(mdot_tot, 1e-4) * cp_gas * dt);
        T_ch = max(500.0, T_flame - deltaT_wall);
        Q_rxn = mdot_prop * cp_gas * (T_flame - 298.15) * dt;
    else
        T_ch   = 293.15;
        cstar  = 1100.0;
        gamma  = 1.4;
        R_gas  = 296.8;
        cv_gas = 743.0;
        Q_rxn  = 0.0;
    end

    % Nozzle Throat Parameters
    A_t    = Nozzle.A_throat;
    Cd     = 0.95;
    if isfield(Nozzle, 'Cd') && Nozzle.Cd > 0, Cd = Nozzle.Cd; end
    P_back = max(Nozzle.P_back, 101325);

    % Unconditionally Stable Chamber Filling Dynamics
    tau_res = max(Node.V / (A_t * Cd * cstar), 1e-4);
    tau_eff = max(tau_res, 0.025);

    % Steady-state target pressure for current inflow
    P_ss = max(P_back, (mdot_tot * cstar) / (A_t * Cd));

    % Exponential transition: P(t+dt) = P_curr * exp(-dt/tau) + P_ss * (1 - exp(-dt/tau))
    decay = exp(-dt / tau_eff);
    P_new = Node.P * decay + P_ss * (1.0 - decay);
    P_new = max(P_new, P_back);

    % Choked nozzle mass flow at current chamber pressure
    PR = P_back / P_new;
    PR_crit = (2.0 / (gamma + 1.0))^(gamma / (gamma - 1.0));

    if P_new <= P_back
        mdot_nozzle = 0.0;
    elseif PR <= PR_crit
        mdot_nozzle = (P_new * A_t * Cd) / cstar;
    else
        mdot_nozzle = Cd * A_t * sqrt(2.0 * P_new * (P_new / (R_gas * T_ch)) * ...
                      (gamma / (gamma - 1.0)) * (PR^(2.0 / gamma) - PR^((gamma + 1.0) / gamma)));
    end

    % Thermodynamic State Closing
    rho_new = P_new / (R_gas * T_ch);
    m_new   = rho_new * Node.V;
    u_new   = cv_gas * T_ch;
    U_new   = m_new * u_new;

    Node_new       = Node;
    Node_new.m     = m_new;
    Node_new.U     = U_new;
    Node_new.u     = u_new;
    Node_new.rho   = rho_new;
    Node_new.T     = T_ch;
    Node_new.P     = P_new;
    Node_new.h     = u_new + P_new / rho_new;
    Node_new.gamma = gamma;
    Node_new.cstar = cstar;
    Node_new.isLit = isLit;
end
