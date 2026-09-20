function [Q_nodes, State] = StepThermalNetwork(State, dt, Env)
% STEPTHERMALNETWORK Master thermal network solver.
% Computes heat exchange for all registered thermal sub-models and updates
% internal wall thermal states.
%
% Inputs:
%   State   - Full simulation State struct
%   dt      - Time step [s]
%   Env     - Ambient conditions struct (.T_amb, .WindVel, .P_atm)
%
% Outputs:
%   Q_nodes - Struct of net thermal energy [J] transferred into each node
%   State   - Updated State struct with advanced thermal wall states

    if nargin < 3 || isempty(Env)
        Env.T_amb   = 293.15;
        Env.WindVel = 1.5;
        Env.P_atm   = 101325;
    end

    Q_nodes = struct();

    % 1. COPV Heat Transfer
    if isfield(State.Nodes, 'TK_N2') && isfield(State.Thermal, 'COPV')
        P_copv = State.Nodes.TK_N2.P;
        T_copv = State.Nodes.TK_N2.T;
        isFlowing = isfield(State, 'COPV_Flowing') && State.COPV_Flowing;

        [Q_copv, State.Thermal.COPV] = COPVThermalModel(P_copv, T_copv, ...
            State.Thermal.COPV, dt, isFlowing, Env.WindVel);
        Q_nodes.TK_N2 = Q_copv;
    else
        Q_nodes.TK_N2 = 0.0;
    end

    % 2. LOX Tank Heat Inleak
    if isfield(State.Nodes, 'TK_O2_01') && isfield(State.Thermal, 'TK_O2_01')
        [Q_ull_ox, Q_liq_ox, State.Thermal.TK_O2_01] = TankThermalModel(...
            State.Nodes.TK_O2_01, State.Thermal.TK_O2_01, dt, Env.T_amb);
        Q_nodes.TK_O2_01_ull = Q_ull_ox;
        Q_nodes.TK_O2_01_liq = Q_liq_ox;
    else
        Q_nodes.TK_O2_01_ull = 0.0;
        Q_nodes.TK_O2_01_liq = 0.0;
    end

    % 3. Fuel Tank Heat Inleak
    if isfield(State.Nodes, 'TK_FU_01') && isfield(State.Thermal, 'TK_FU_01')
        [Q_ull_fu, Q_liq_fu, State.Thermal.TK_FU_01] = TankThermalModel(...
            State.Nodes.TK_FU_01, State.Thermal.TK_FU_01, dt, Env.T_amb);
        Q_nodes.TK_FU_01_ull = Q_ull_fu;
        Q_nodes.TK_FU_01_liq = Q_liq_fu;
    else
        Q_nodes.TK_FU_01_ull = 0.0;
        Q_nodes.TK_FU_01_liq = 0.0;
    end

    % 4. Regenerative Chamber Cooling
    if isfield(State.Nodes, 'SKIPPER') && isfield(State.Thermal, 'Regen')
        T_cham = State.Nodes.SKIPPER.T;
        P_cham = State.Nodes.SKIPPER.P;
        mdot_fu = 0.0;
        T_fu_in = 293.15;
        if isfield(State, 'Regen_mdot'), mdot_fu = State.Regen_mdot; end
        if isfield(State, 'Regen_T_in'), T_fu_in = State.Regen_T_in; end
        OF = 1.2;
        if isfield(State, 'OF'), OF = State.OF; end

        [Q_cham_loss, Q_fuel_gain, State.Thermal.Regen] = RegenThermalModel(...
            T_cham, P_cham, mdot_fu, T_fu_in, State.Thermal.Regen, dt, OF);
        Q_nodes.SKIPPER_loss = Q_cham_loss;
        Q_nodes.Regen_gain   = Q_fuel_gain;
    else
        Q_nodes.SKIPPER_loss = 0.0;
        Q_nodes.Regen_gain   = 0.0;
    end
end
