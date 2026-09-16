function [State, FlowRates, System] = StepSimulation(State, dt, System, LinkStates)
% STEPSIMULATION Advances the TOAD lumped parameter fluid system by one time step (dt)
% using discrete Mass & Internal Energy (m, U) tracking, quasi-steady branch hydraulics,
% and precomputed thermodynamic property tables.
%
% Inputs:
%   State      - Simulation State struct holding all named nodes, thermal states, Time
%   dt         - Discrete time step [s]
%   System     - System struct holding .Links, .Link.State, .Env
%   LinkStates - (Optional) Commanded valve states struct (fractions 0.0 to 1.0 of MaxCv),
%                e.g. from TOADSchedule(t) or direct struct overrides.
%
% Outputs:
%   State      - Updated State struct
%   FlowRates  - Struct of link mass flow rates [kg/s]
%   System     - (Optional) Updated System struct with synchronized .Links and .Link.State

    if nargin < 4, LinkStates = []; end

    %% 1. Actuator Dynamics, Valve States, and In-Loop Overrides
    % Ensure System.Link.State exists for direct in-loop access/overrides
    if ~isfield(System, 'Link') || ~isfield(System.Link, 'State')
        System.Link.State = struct();
        linkNames = fieldnames(System.Links);
        for k = 1:length(linkNames)
            System.Link.State.(linkNames{k}) = System.Links.(linkNames{k}).State;
        end
    end

    % 1a. Detect direct in-loop overrides (e.g. System.Link.State.(valve) = val)
    linkNames = fieldnames(System.Links);
    for k = 1:length(linkNames)
        fn = linkNames{k};
        if isfield(System.Link.State, fn)
            if System.Link.State.(fn) ~= System.Links.(fn).State
                System.Links.(fn).State = max(0.0, min(1.0, System.Link.State.(fn)));
            end
        else
            System.Link.State.(fn) = System.Links.(fn).State;
        end
    end

    % 1b. Apply commanded LinkStates (fraction 0.0 to 1.0 of MaxCv) with actuator lag
    if ~isempty(LinkStates)
        if isstruct(LinkStates) && isscalar(LinkStates)
            cmdFields = fieldnames(LinkStates);
            for k = 1:length(cmdFields)
                rawName = cmdFields{k};
                fn = strrep(rawName, '-', '_');
                if isfield(System.Links, fn)
                    target = max(0.0, min(1.0, LinkStates.(rawName)));
                    if isfield(System.Links.(fn), 'Tau') && System.Links.(fn).Tau > dt
                        tau = System.Links.(fn).Tau;
                        alpha = 1.0 - exp(-dt / tau);
                        cur = System.Links.(fn).State;
                        new_state = cur + alpha * (target - cur);
                        if abs(new_state - target) < 1e-4
                            new_state = target;
                        end
                        System.Links.(fn).State = new_state;
                    else
                        System.Links.(fn).State = target;
                    end
                end
            end
        elseif iscell(LinkStates) || (isstruct(LinkStates) && ~isscalar(LinkStates))
            for k = 1:length(LinkStates)
                if iscell(LinkStates)
                    entry = LinkStates{k};
                else
                    entry = LinkStates(k);
                end
                if isstruct(entry)
                    if isfield(entry, 'Name')
                        vName = entry.Name;
                    elseif isfield(entry, 'Valve')
                        vName = entry.Valve;
                    else
                        continue;
                    end
                    if isfield(entry, 'Value')
                        vVal = entry.Value;
                    elseif isfield(entry, 'State')
                        vVal = entry.State;
                    elseif isfield(entry, 'TargetCv')
                        vVal = entry.TargetCv;
                    else
                        continue;
                    end
                elseif iscell(entry) && numel(entry) >= 2
                    vName = entry{1};
                    vVal = entry{2};
                else
                    continue;
                end
                fn = strrep(vName, '-', '_');
                if isfield(System.Links, fn)
                    if isfield(System.Links.(fn), 'MaxCv') && System.Links.(fn).MaxCv > 0 && vVal > 1.0
                        target = max(0.0, min(1.0, vVal / System.Links.(fn).MaxCv));
                    else
                        target = max(0.0, min(1.0, vVal));
                    end
                    if isfield(System.Links.(fn), 'Tau') && System.Links.(fn).Tau > dt
                        tau = System.Links.(fn).Tau;
                        alpha = 1.0 - exp(-dt / tau);
                        cur = System.Links.(fn).State;
                        new_state = cur + alpha * (target - cur);
                        if abs(new_state - target) < 1e-4
                            new_state = target;
                        end
                        System.Links.(fn).State = new_state;
                    else
                        System.Links.(fn).State = target;
                    end
                end
            end
        end
    end

    % 1c. Update active Cv = MaxCv * State for all valves and keep alias synchronized
    for k = 1:length(linkNames)
        fn = linkNames{k};
        if isfield(System.Links.(fn), 'MaxCv')
            System.Links.(fn).Cv = System.Links.(fn).MaxCv * System.Links.(fn).State;
        end
        System.Link.State.(fn) = System.Links.(fn).State;
    end
    State.LinkStates = System.Link.State;


    FlowRates = struct();

    %% 2. Ground Bulk N2 Fill Link (if present)
    if isfield(State.Nodes, 'TK_N2_BULK') && isfield(System.Links, 'BV_N2_FILL')
        if System.Links.BV_N2_FILL.Cv > 0.001
            [P_bulk, Props_bulk] = getNodeLinkState(State.Nodes, 'TK_N2_BULK', 'Up');
            [P_copv, Props_copv] = getNodeLinkState(State.Nodes, 'TK_N2', 'Down');
            [mdot_fill, ~, ~] = CalculateLinkFlow(System.Links.BV_N2_FILL, P_bulk, P_copv, Props_bulk, Props_copv, dt);
            FlowRates.BV_N2_FILL = mdot_fill;
        else
            FlowRates.BV_N2_FILL = 0.0;
        end
    else
        FlowRates.BV_N2_FILL = 0.0;
    end

    %% 3. Regulated N2 Pressurization Network (Quasi-Steady Equilibrium)
    % Physics: Dome-loaded regulator maintains downstream manifold pressure matching demand
    P_copv = State.Nodes.TK_N2.P;
    P_set  = System.Links.REG_873D.P_set;
    Droop  = System.Links.REG_873D.Droop;
    SPE    = System.Links.REG_873D.SPE;
    P_reg_target = P_set - SPE * (P_copv - P_set);

    % Solenoid & check valve conductances
    Cv_n2_02    = System.Links.BV_N2_02.Cv;
    K_press_ox  = System.Links.Press_OX.Cv * 2.402e-5;
    K_press_fu  = System.Links.Press_FU.Cv * 2.402e-5;
    P_ox_ull    = State.Nodes.TK_O2_01.Ullage.P;
    P_fu_ull    = State.Nodes.TK_FU_01.Ullage.P;
    rho_n2_pman = max(State.Nodes.Press_Manifold.rho, 1.0);

    % Purge line conductances
    K_pur_ox = System.Links.SV_N2_05.Cv * 2.402e-5;
    K_pur_fu = System.Links.SV_N2_06.Cv * 2.402e-5;
    K_pur_dt = System.Links.SV_N2_07.Cv * 2.402e-5;
    P_skip   = State.Nodes.SKIPPER.P;

    % Demand calculations when valves open
    mdot_demand_ox = 0.0;
    mdot_demand_fu = 0.0;
    if Cv_n2_02 > 0.01
        if P_reg_target > P_ox_ull
            mdot_demand_ox = K_press_ox * sqrt(rho_n2_pman * (P_reg_target - P_ox_ull));
        end
        if P_reg_target > P_fu_ull
            mdot_demand_fu = K_press_fu * sqrt(rho_n2_pman * (P_reg_target - P_fu_ull));
        end
    end
    mdot_pur_ox = 0.0;
    mdot_pur_fu = 0.0;
    mdot_pur_dt = 0.0;
    if K_pur_ox > 1e-7 && P_reg_target > P_skip
        mdot_pur_ox = K_pur_ox * sqrt(rho_n2_pman * (P_reg_target - P_skip));
    end
    if K_pur_fu > 1e-7 && P_reg_target > P_skip
        mdot_pur_fu = K_pur_fu * sqrt(rho_n2_pman * (P_reg_target - P_skip));
    end
    if K_pur_dt > 1e-7 && P_reg_target > P_skip
        mdot_pur_dt = K_pur_dt * sqrt(rho_n2_pman * (P_reg_target - P_skip));
    end

    total_demand = mdot_demand_ox + mdot_demand_fu + mdot_pur_ox + mdot_pur_fu + mdot_pur_dt;
    mdot_reg_max = max(0.85 * (P_copv / 3e7), 0.01);
    openFrac     = min(1.0, total_demand / mdot_reg_max);
    P_man_equil  = max(101325, P_reg_target - Droop * openFrac);

    % Actual flows under equilibrium pressure (Courant-stable ullage flow matching)
    if Cv_n2_02 > 0.01
        V_ull_ox = State.Nodes.TK_O2_01.Ullage.V;
        V_ull_fu = State.Nodes.TK_FU_01.Ullage.V;
        R_n2 = 296.8;
        T_n2_in = max(State.Nodes.TK_N2.T, 100);

        % Maximum inflow per step dt to prevent overshooting P_man_equil (eliminates 3.5 psi sawtooth)
        cap_flow_ox = max(0.0, (V_ull_ox / (R_n2 * T_n2_in * dt)) * max(0, P_man_equil - P_ox_ull));
        cap_flow_fu = max(0.0, (V_ull_fu / (R_n2 * T_n2_in * dt)) * max(0, P_man_equil - P_fu_ull));

        raw_flow_ox = K_press_ox * sqrt(rho_n2_pman * max(0, P_man_equil - P_ox_ull));
        raw_flow_fu = K_press_fu * sqrt(rho_n2_pman * max(0, P_man_equil - P_fu_ull));

        FlowRates.Press_OX = min(raw_flow_ox, cap_flow_ox);
        FlowRates.Press_FU = min(raw_flow_fu, cap_flow_fu);
    else
        FlowRates.Press_OX = 0.0;
        FlowRates.Press_FU = 0.0;
    end
    FlowRates.SV_N2_05 = K_pur_ox * sqrt(rho_n2_pman * max(0, P_man_equil - P_skip));
    FlowRates.SV_N2_06 = K_pur_fu * sqrt(rho_n2_pman * max(0, P_man_equil - P_skip));
    FlowRates.SV_N2_07 = K_pur_dt * sqrt(rho_n2_pman * max(0, P_man_equil - P_skip));

    FlowRates.REG_873D = FlowRates.Press_OX + FlowRates.Press_FU + ...
                         FlowRates.SV_N2_05 + FlowRates.SV_N2_06 + FlowRates.SV_N2_07;
    FlowRates.BV_N2_02 = FlowRates.Press_OX + FlowRates.Press_FU;

    % Update Manifold Nodal States (Smooth, non-oscillatory)
    State.Nodes.Purge_Manifold.P   = P_man_equil;
    State.Nodes.Purge_Manifold.T   = State.Nodes.TK_N2.T;
    State.Nodes.Purge_Manifold.rho = P_man_equil / (296.8 * State.Nodes.TK_N2.T);
    State.Nodes.Purge_Manifold.m   = State.Nodes.Purge_Manifold.rho * State.Nodes.Purge_Manifold.V;
    State.Nodes.Purge_Manifold.h   = State.Nodes.TK_N2.h;
    State.Nodes.Purge_Manifold.U   = State.Nodes.Purge_Manifold.m * State.Nodes.TK_N2.u;

    State.Nodes.Press_Manifold.P   = P_man_equil;
    State.Nodes.Press_Manifold.T   = State.Nodes.TK_N2.T;
    State.Nodes.Press_Manifold.rho = State.Nodes.Purge_Manifold.rho;
    State.Nodes.Press_Manifold.m   = State.Nodes.Press_Manifold.rho * State.Nodes.Press_Manifold.V;
    State.Nodes.Press_Manifold.h   = State.Nodes.TK_N2.h;
    State.Nodes.Press_Manifold.U   = State.Nodes.Press_Manifold.m * State.Nodes.TK_N2.u;

    %% 4. Exact Series Liquid Feed Branch Solver (Zero Chatter)
    % A. OXIDIZER FEED BRANCH (LOX)
    % TK-O2-01.Liquid -> OX_Line_1 -> BV-02-03 -> Inter_OX -> BV-02-04 -> Post_Throttle_OX -> Inj_OX -> SKIPPER
    P_tank_ox  = State.Nodes.TK_O2_01.Liquid.P;
    rho_ox     = max(State.Nodes.TK_O2_01.Liquid.rho, 900.0);
    T_lox      = State.Nodes.TK_O2_01.Liquid.T;
    h_lox      = State.Nodes.TK_O2_01.Liquid.h;
    u_lox      = State.Nodes.TK_O2_01.Liquid.u;

    Cv_main_ox = System.Links.BV_02_03.Cv;
    Cv_thrt_ox = System.Links.BV_02_04.Cv;
    K_pipe_ox  = 1.581e-4; % Pipe conductance
    K_inj_ox   = System.Links.Inj_OX.Cd * System.Links.Inj_OX.A;

    mdot_ox = 0.0;
    if Cv_main_ox > 0.01 && Cv_thrt_ox > 0.01
        K_main_ox = Cv_main_ox * 2.402e-5 * 0.7071;
        K_thrt_ox = Cv_thrt_ox * 2.402e-5 * 0.7071;
        invK2_ox  = (1.0 / K_pipe_ox^2) + (1.0 / K_main_ox^2) + (1.0 / K_thrt_ox^2) + (1.0 / K_inj_ox^2);
        K_eq_ox   = 1.0 / sqrt(invK2_ox);
        DP_ox     = max(0.0, P_tank_ox - P_skip);
        mdot_ox   = K_eq_ox * sqrt(2.0 * rho_ox * DP_ox);

        P_pre_ox  = P_tank_ox - (mdot_ox^2) / (2.0 * rho_ox * K_pipe_ox^2);
        P_intr_ox = P_pre_ox  - (mdot_ox^2) / (2.0 * rho_ox * K_main_ox^2);
        P_post_ox = P_intr_ox - (mdot_ox^2) / (2.0 * rho_ox * K_thrt_ox^2);
    else
        if Cv_main_ox > 0.01
            P_pre_ox  = P_tank_ox;
            P_intr_ox = P_tank_ox;
            P_post_ox = P_skip;
        else
            P_pre_ox  = P_tank_ox;
            P_intr_ox = 101325;
            P_post_ox = P_skip;
        end
    end

    FlowRates.OX_Line_1   = mdot_ox;
    FlowRates.BV_02_03    = mdot_ox;
    FlowRates.BV_02_04    = mdot_ox;
    FlowRates.OX_Inj_Line = mdot_ox + FlowRates.SV_N2_05;
    FlowRates.Inj_OX      = mdot_ox + FlowRates.SV_N2_05;

    % Direct O(1) state assignment for OX line nodes
    State.Nodes.Pre_Main_OX.P   = P_pre_ox;
    State.Nodes.Pre_Main_OX.T   = T_lox;
    State.Nodes.Pre_Main_OX.rho = rho_ox;
    State.Nodes.Pre_Main_OX.h   = h_lox;
    State.Nodes.Pre_Main_OX.m   = rho_ox * State.Nodes.Pre_Main_OX.V;
    State.Nodes.Pre_Main_OX.U   = State.Nodes.Pre_Main_OX.m * u_lox;

    State.Nodes.Inter_OX.P   = P_intr_ox;
    State.Nodes.Inter_OX.T   = T_lox;
    State.Nodes.Inter_OX.rho = rho_ox;
    State.Nodes.Inter_OX.h   = h_lox;
    State.Nodes.Inter_OX.m   = rho_ox * State.Nodes.Inter_OX.V;
    State.Nodes.Inter_OX.U   = State.Nodes.Inter_OX.m * u_lox;

    State.Nodes.Post_Throttle_OX.P   = P_post_ox;
    State.Nodes.Post_Throttle_OX.T   = T_lox;
    State.Nodes.Post_Throttle_OX.rho = rho_ox;
    State.Nodes.Post_Throttle_OX.h   = h_lox;
    State.Nodes.Post_Throttle_OX.m   = rho_ox * State.Nodes.Post_Throttle_OX.V;
    State.Nodes.Post_Throttle_OX.U   = State.Nodes.Post_Throttle_OX.m * u_lox;

    State.Nodes.OX_Manifold.P   = P_post_ox;
    State.Nodes.OX_Manifold.T   = T_lox;
    State.Nodes.OX_Manifold.rho = rho_ox;
    State.Nodes.OX_Manifold.h   = h_lox;
    State.Nodes.OX_Manifold.m   = rho_ox * State.Nodes.OX_Manifold.V;
    State.Nodes.OX_Manifold.U   = State.Nodes.OX_Manifold.m * u_lox;

    % B. FUEL FEED BRANCH (IPA with Regenerative Jacket)
    % TK-FU-01.Liquid -> FU_Line_1 -> BV-FU-03 -> Inter_FU -> BV-FU-04 -> Post_Throttle_FU -> Regen -> FU_Manifold -> Inj_FU -> SKIPPER
    P_tank_fu  = State.Nodes.TK_FU_01.Liquid.P;
    rho_fu     = max(State.Nodes.TK_FU_01.Liquid.rho, 700.0);
    T_ipa      = State.Nodes.TK_FU_01.Liquid.T;
    h_ipa      = State.Nodes.TK_FU_01.Liquid.h;
    u_ipa      = State.Nodes.TK_FU_01.Liquid.u;

    Cv_main_fu = System.Links.BV_FU_03.Cv;
    Cv_thrt_fu = System.Links.BV_FU_04.Cv;
    K_pipe_fu  = 1.581e-4;
    K_inj_fu   = System.Links.Inj_FU.Cd * System.Links.Inj_FU.A;

    DP_regen = 0.0;
    if isfield(State.Thermal, 'Regen') && isfield(State.Thermal.Regen, 'DeltaP')
        DP_regen = State.Thermal.Regen.DeltaP;
    end

    mdot_fu = 0.0;
    if Cv_main_fu > 0.01 && Cv_thrt_fu > 0.01
        K_main_fu = Cv_main_fu * 2.402e-5 * 0.7071;
        K_thrt_fu = Cv_thrt_fu * 2.402e-5 * 0.7071;
        invK2_fu  = (1.0 / K_pipe_fu^2) + (1.0 / K_main_fu^2) + (1.0 / K_thrt_fu^2) + (1.0 / K_inj_fu^2);
        K_eq_fu   = 1.0 / sqrt(invK2_fu);
        DP_fu_net = max(0.0, P_tank_fu - P_skip - DP_regen);
        mdot_fu   = K_eq_fu * sqrt(2.0 * rho_fu * DP_fu_net);

        % Smooth fixed-point relaxation with regen jacket impedance (eliminates 100 Hz limit cycle):
        if mdot_fu > 0.01 && DP_regen > 500 && isfield(State, 'Regen_mdot') && State.Regen_mdot > 0.01
            DP_regen_est = DP_regen * (mdot_fu / State.Regen_mdot)^1.75;
            DP_regen_damped = 0.5 * DP_regen + 0.5 * DP_regen_est;
            DP_fu_net = max(0.0, P_tank_fu - P_skip - DP_regen_damped);
            mdot_fu   = K_eq_fu * sqrt(2.0 * rho_fu * DP_fu_net);
        end

        P_pre_fu  = P_tank_fu - (mdot_fu^2) / (2.0 * rho_fu * K_pipe_fu^2);
        P_intr_fu = P_pre_fu  - (mdot_fu^2) / (2.0 * rho_fu * K_main_fu^2);
        P_post_fu = P_intr_fu - (mdot_fu^2) / (2.0 * rho_fu * K_thrt_fu^2);
        P_man_fu  = max(P_skip, P_post_fu - DP_regen);
    else
        if Cv_main_fu > 0.01
            P_pre_fu  = P_tank_fu;
            P_intr_fu = P_tank_fu;
            P_post_fu = P_skip;
            P_man_fu  = P_skip;
        else
            P_pre_fu  = P_tank_fu;
            P_intr_fu = 101325;
            P_post_fu = P_skip;
            P_man_fu  = P_skip;
        end
    end

    FlowRates.FU_Line_1   = mdot_fu;
    FlowRates.BV_FU_03    = mdot_fu;
    FlowRates.BV_FU_04    = mdot_fu;
    FlowRates.FU_Inj_Line = mdot_fu + FlowRates.SV_N2_06;
    FlowRates.Inj_FU      = mdot_fu + FlowRates.SV_N2_06;

    % Direct O(1) state assignment for Fuel line nodes
    State.Nodes.Pre_Main_FU.P   = P_pre_fu;
    State.Nodes.Pre_Main_FU.T   = T_ipa;
    State.Nodes.Pre_Main_FU.rho = rho_fu;
    State.Nodes.Pre_Main_FU.h   = h_ipa;
    State.Nodes.Pre_Main_FU.m   = rho_fu * State.Nodes.Pre_Main_FU.V;
    State.Nodes.Pre_Main_FU.U   = State.Nodes.Pre_Main_FU.m * u_ipa;

    State.Nodes.Inter_FU.P   = P_intr_fu;
    State.Nodes.Inter_FU.T   = T_ipa;
    State.Nodes.Inter_FU.rho = rho_fu;
    State.Nodes.Inter_FU.h   = h_ipa;
    State.Nodes.Inter_FU.m   = rho_fu * State.Nodes.Inter_FU.V;
    State.Nodes.Inter_FU.U   = State.Nodes.Inter_FU.m * u_ipa;

    State.Nodes.Post_Throttle_FU.P   = P_post_fu;
    State.Nodes.Post_Throttle_FU.T   = T_ipa;
    State.Nodes.Post_Throttle_FU.rho = rho_fu;
    State.Nodes.Post_Throttle_FU.h   = h_ipa;
    State.Nodes.Post_Throttle_FU.m   = rho_fu * State.Nodes.Post_Throttle_FU.V;
    State.Nodes.Post_Throttle_FU.U   = State.Nodes.Post_Throttle_FU.m * u_ipa;

    % C. DART TORCH IGNITER BRANCH
    Cv_dart_ox = System.Links.SV_DART_OX.Cv;
    Cv_dart_fu = System.Links.SV_DART_FU.Cv;
    K_dart_ox  = Cv_dart_ox * 2.402e-5 * 0.7071;
    K_dart_fu  = Cv_dart_fu * 2.402e-5 * 0.7071;
    P_dart     = State.Nodes.DART_Chamber.P;

    FlowRates.SV_DART_OX = 0.0;
    FlowRates.SV_DART_FU = 0.0;
    if Cv_dart_ox > 0.001 && P_intr_ox > P_dart
        FlowRates.SV_DART_OX = K_dart_ox * sqrt(2.0 * rho_ox * (P_intr_ox - P_dart));
    end
    if Cv_dart_fu > 0.001 && P_intr_fu > P_dart
        FlowRates.SV_DART_FU = K_dart_fu * sqrt(2.0 * rho_fu * (P_intr_fu - P_dart));
    end

    %% 5. Thermal Network & Regen Heat Transfer
    State.COPV_Flowing = (FlowRates.REG_873D > 1e-4) || (FlowRates.BV_N2_FILL > 1e-4);
    State.Regen_mdot   = abs(FlowRates.Inj_FU);
    if FlowRates.Inj_FU > 1e-4
        State.OF = FlowRates.Inj_OX / FlowRates.Inj_FU;
    else
        State.OF = 1.2;
    end
    [Q_nodes, State] = StepThermalNetwork(State, dt, System.Env);

    % FU Manifold (warmed by regenerative cooling)
    h_man_fu = h_ipa;
    if mdot_fu > 1e-4 && isfield(Q_nodes, 'Regen_gain')
        h_man_fu = h_ipa + Q_nodes.Regen_gain / max(mdot_fu * dt, 1e-6);
    end
    State.Nodes.FU_Manifold.P   = P_man_fu;
    State.Nodes.FU_Manifold.T   = T_ipa + (h_man_fu - h_ipa) / 2600.0;
    State.Nodes.FU_Manifold.rho = rho_fu;
    State.Nodes.FU_Manifold.h   = h_man_fu;
    State.Nodes.FU_Manifold.m   = rho_fu * State.Nodes.FU_Manifold.V;
    State.Nodes.FU_Manifold.U   = State.Nodes.FU_Manifold.m * (h_man_fu - P_man_fu / rho_fu);

    %% 6. Advance Combustors (DART & SKIPPER)
    DART_Inflow.mdot_ox = FlowRates.SV_DART_OX;
    DART_Inflow.mdot_fu = FlowRates.SV_DART_FU;
    DART_Inflow.mdot_n2 = FlowRates.SV_N2_07;
    DART_Inflow.h_ox    = State.Nodes.Inter_OX.h;
    DART_Inflow.h_fu    = State.Nodes.Inter_FU.h;
    DART_Inflow.h_n2    = State.Nodes.Purge_Manifold.h;

    DART_Nozzle.A_throat = System.Links.DART_Nozzle.A;
    DART_Nozzle.Cd       = System.Links.DART_Nozzle.Cd;
    DART_Nozzle.P_back   = State.Nodes.SKIPPER.P;

    SparkActive = (System.Links.Spark.Cv > 0.5);
    [State.Nodes.DART_Chamber, FlowRates.DART_Nozzle, ~, ~] = StepCombustor(...
        State.Nodes.DART_Chamber, DART_Inflow, DART_Nozzle, dt, SparkActive, false, 0.0);

    SKIPPER_Inflow.mdot_ox = FlowRates.Inj_OX;
    SKIPPER_Inflow.mdot_fu = FlowRates.Inj_FU;
    SKIPPER_Inflow.mdot_n2 = FlowRates.DART_Nozzle;
    SKIPPER_Inflow.h_ox    = State.Nodes.OX_Manifold.h;
    SKIPPER_Inflow.h_fu    = State.Nodes.FU_Manifold.h;
    SKIPPER_Inflow.h_n2    = State.Nodes.DART_Chamber.h;

    SKIPPER_Nozzle.A_throat = System.Links.Nozzle.A;
    SKIPPER_Nozzle.Cd       = System.Links.Nozzle.Cd;
    SKIPPER_Nozzle.P_back   = State.Nodes.Atmosphere.P;

    TorchActive = State.Nodes.DART_Chamber.isLit && (FlowRates.DART_Nozzle > 0.001);
    [State.Nodes.SKIPPER, FlowRates.Nozzle, ~, ~] = StepCombustor(...
        State.Nodes.SKIPPER, SKIPPER_Inflow, SKIPPER_Nozzle, dt, false, TorchActive, Q_nodes.SKIPPER_loss);

    %% 7. Advance Conserved States (Mass & Energy) for Capacitances
    % A. Bulk N2 Supply
    if isfield(State.Nodes, 'TK_N2_BULK')
        mdot_bulk_out = FlowRates.BV_N2_FILL;
        if abs(mdot_bulk_out) > 1e-6
            h_bulk = State.Nodes.TK_N2_BULK.h;
            dm_bulk = -mdot_bulk_out * dt;
            dU_bulk = -(mdot_bulk_out * h_bulk) * dt;

            m_bulk_new = max(State.Nodes.TK_N2_BULK.m + dm_bulk, 1e-4);
            U_bulk_new = State.Nodes.TK_N2_BULK.U + dU_bulk;
            rho_bulk   = m_bulk_new / State.Nodes.TK_N2_BULK.V;
            u_bulk     = U_bulk_new / m_bulk_new;

            props_bulk = FluidProperties('Nitrogen', 'From_u_rho', u_bulk, rho_bulk);
            State.Nodes.TK_N2_BULK.m   = m_bulk_new;
            State.Nodes.TK_N2_BULK.U   = U_bulk_new;
            State.Nodes.TK_N2_BULK.u   = u_bulk;
            State.Nodes.TK_N2_BULK.rho = rho_bulk;
            State.Nodes.TK_N2_BULK.P   = props_bulk.P;
            State.Nodes.TK_N2_BULK.T   = props_bulk.T;
            State.Nodes.TK_N2_BULK.h   = props_bulk.h;
        end
    end

    % B. COPV (TK-N2)
    mdot_copv_in  = FlowRates.BV_N2_FILL;
    h_in_copv     = State.Nodes.TK_N2.h;
    if isfield(State.Nodes, 'TK_N2_BULK'), h_in_copv = State.Nodes.TK_N2_BULK.h; end
    mdot_copv_out = FlowRates.REG_873D;
    h_copv        = State.Nodes.TK_N2.h;

    dm_copv = (mdot_copv_in - mdot_copv_out) * dt;
    dU_copv = (mdot_copv_in * h_in_copv - mdot_copv_out * h_copv) * dt + Q_nodes.TK_N2;

    m_copv_new = max(State.Nodes.TK_N2.m + dm_copv, 1e-4);
    U_copv_new = State.Nodes.TK_N2.U + dU_copv;
    rho_copv   = m_copv_new / State.Nodes.TK_N2.V;
    u_copv     = U_copv_new / m_copv_new;

    props_copv = FluidProperties('Nitrogen', 'From_u_rho', u_copv, rho_copv);
    State.Nodes.TK_N2.m   = m_copv_new;
    State.Nodes.TK_N2.U   = U_copv_new;
    State.Nodes.TK_N2.u   = u_copv;
    State.Nodes.TK_N2.rho = rho_copv;
    State.Nodes.TK_N2.P   = props_copv.P;
    State.Nodes.TK_N2.T   = props_copv.T;
    State.Nodes.TK_N2.h   = props_copv.h;

    % C. LOX Tank (TK-O2-01, Two-Zone Stratified)
    mdot_ox_press_in = FlowRates.Press_OX;
    mdot_ox_feed_out = mdot_ox;

    % Ullage Zone (N2 gas)
    dm_ull_ox = mdot_ox_press_in * dt;
    dU_ull_ox = (mdot_ox_press_in * State.Nodes.Press_Manifold.h) * dt + Q_nodes.TK_O2_01_ull;
    m_ull_ox  = max(State.Nodes.TK_O2_01.Ullage.m + dm_ull_ox, 1e-4);
    U_ull_ox  = State.Nodes.TK_O2_01.Ullage.U + dU_ull_ox;

    % Liquid Zone (LOX liquid)
    dm_liq_ox = -mdot_ox_feed_out * dt;
    dU_liq_ox = -(mdot_ox_feed_out * h_lox) * dt + Q_nodes.TK_O2_01_liq;
    m_liq_ox  = max(State.Nodes.TK_O2_01.Liquid.m + dm_liq_ox, 0.0);
    U_liq_ox  = State.Nodes.TK_O2_01.Liquid.U + dU_liq_ox;

    u_liq_ox   = U_liq_ox / max(m_liq_ox, 1e-4);
    props_lox  = FluidProperties('Oxygen', 'From_u_rho', u_liq_ox, 1141);
    V_liq_ox   = m_liq_ox / max(props_lox.rho, 100);
    V_ull_ox   = max(State.Nodes.TK_O2_01.V - V_liq_ox, 1e-4);

    rho_ull_ox = m_ull_ox / V_ull_ox;
    u_ull_ox   = U_ull_ox / m_ull_ox;
    props_ull_ox = FluidProperties('Nitrogen', 'From_u_rho', u_ull_ox, rho_ull_ox);

    P_tank_ox = props_ull_ox.P;
    State.Nodes.TK_O2_01.P = P_tank_ox;
    State.Nodes.TK_O2_01.Ullage.m   = m_ull_ox;
    State.Nodes.TK_O2_01.Ullage.U   = U_ull_ox;
    State.Nodes.TK_O2_01.Ullage.V   = V_ull_ox;
    State.Nodes.TK_O2_01.Ullage.P   = P_tank_ox;
    State.Nodes.TK_O2_01.Ullage.T   = props_ull_ox.T;
    State.Nodes.TK_O2_01.Ullage.rho = rho_ull_ox;
    State.Nodes.TK_O2_01.Ullage.h   = props_ull_ox.h;

    State.Nodes.TK_O2_01.Liquid.m   = m_liq_ox;
    State.Nodes.TK_O2_01.Liquid.U   = U_liq_ox;
    State.Nodes.TK_O2_01.Liquid.V   = V_liq_ox;
    State.Nodes.TK_O2_01.Liquid.P   = P_tank_ox;
    State.Nodes.TK_O2_01.Liquid.T   = props_lox.T;
    State.Nodes.TK_O2_01.Liquid.rho = props_lox.rho;
    State.Nodes.TK_O2_01.Liquid.h   = props_lox.h;

    % D. Fuel Tank (TK-FU-01, Two-Zone Stratified)
    mdot_fu_press_in = FlowRates.Press_FU;
    mdot_fu_feed_out = mdot_fu;

    % Ullage Zone (N2 gas)
    dm_ull_fu = mdot_fu_press_in * dt;
    dU_ull_fu = (mdot_fu_press_in * State.Nodes.Press_Manifold.h) * dt + Q_nodes.TK_FU_01_ull;
    m_ull_fu  = max(State.Nodes.TK_FU_01.Ullage.m + dm_ull_fu, 1e-4);
    U_ull_fu  = State.Nodes.TK_FU_01.Ullage.U + dU_ull_fu;

    % Liquid Zone (IPA liquid)
    dm_liq_fu = -mdot_fu_feed_out * dt;
    dU_liq_fu = -(mdot_fu_feed_out * h_ipa) * dt + Q_nodes.TK_FU_01_liq;
    m_liq_fu  = max(State.Nodes.TK_FU_01.Liquid.m + dm_liq_fu, 0.0);
    U_liq_fu  = State.Nodes.TK_FU_01.Liquid.U + dU_liq_fu;

    u_liq_fu   = U_liq_fu / max(m_liq_fu, 1e-4);
    props_ipa  = FluidProperties('IPA', 'From_u_rho', u_liq_fu, 786);
    V_liq_fu   = m_liq_fu / max(props_ipa.rho, 100);
    V_ull_fu   = max(State.Nodes.TK_FU_01.V - V_liq_fu, 1e-4);

    rho_ull_fu = m_ull_fu / V_ull_fu;
    u_ull_fu   = U_ull_fu / m_ull_fu;
    props_ull_fu = FluidProperties('Nitrogen', 'From_u_rho', u_ull_fu, rho_ull_fu);

    P_tank_fu = props_ull_fu.P;
    State.Nodes.TK_FU_01.P = P_tank_fu;
    State.Nodes.TK_FU_01.Ullage.m   = m_ull_fu;
    State.Nodes.TK_FU_01.Ullage.U   = U_ull_fu;
    State.Nodes.TK_FU_01.Ullage.V   = V_ull_fu;
    State.Nodes.TK_FU_01.Ullage.P   = P_tank_fu;
    State.Nodes.TK_FU_01.Ullage.T   = props_ull_fu.T;
    State.Nodes.TK_FU_01.Ullage.rho = rho_ull_fu;
    State.Nodes.TK_FU_01.Ullage.h   = props_ull_fu.h;

    State.Nodes.TK_FU_01.Liquid.m   = m_liq_fu;
    State.Nodes.TK_FU_01.Liquid.U   = U_liq_fu;
    State.Nodes.TK_FU_01.Liquid.V   = V_liq_fu;
    State.Nodes.TK_FU_01.Liquid.P   = P_tank_fu;
    State.Nodes.TK_FU_01.Liquid.T   = props_ipa.T;
    State.Nodes.TK_FU_01.Liquid.rho = props_ipa.rho;
    State.Nodes.TK_FU_01.Liquid.h   = props_ipa.h;

    %% 8. Advance Time
    State.Time = State.Time + dt;
end

%% --- Helper: Extract Node Pressure and Properties for Link Flow ---
function [P, Props] = getNodeLinkState(Nodes, nodeName, endType)
    if ~isfield(Nodes, nodeName)
        error('getNodeLinkState: Unknown node name "%s"', nodeName);
    end
    N = Nodes.(nodeName);

    if strcmp(N.Type, 'TwoZoneTank')
        if strcmp(endType, 'Up')
            P = N.Liquid.P;
            Props = N.Liquid;
        else
            P = N.Ullage.P;
            Props = N.Ullage;
        end
    else
        P = N.P;
        Props = N;
    end
end
