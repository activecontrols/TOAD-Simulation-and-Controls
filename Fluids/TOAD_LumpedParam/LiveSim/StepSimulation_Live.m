function [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates)
% STEPSIMULATION_LIVE Advances the merged TOAD fluid system by one timestep dt
% using a first-principles Control Volume (CV) nodal massflow dynamics solver
% with internal sub-stepping for unconditional stability and physical fidelity.

    if nargin < 4, LinkStates = []; end

    %% 1. Actuator Dynamics, Valve States, and In-Loop Overrides
    if ~isfield(System, 'Link') || ~isfield(System.Link, 'State')
        System.Link.State = struct();
        linkNames = fieldnames(System.Links);
        for k = 1:length(linkNames)
            System.Link.State.(linkNames{k}) = System.Links.(linkNames{k}).State;
        end
    end

    % Direct in-loop overrides
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

    % Commanded LinkStates with actuator lag
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
                        if abs(new_state - target) < 1e-2
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

    % Update effective Cv and synchronize Link.State for all links
    for k = 1:length(linkNames)
        fn = linkNames{k};
        if isfield(System.Links.(fn), 'MaxCv')
            System.Links.(fn).Cv = System.Links.(fn).MaxCv * System.Links.(fn).State;
        end
        System.Link.State.(fn) = System.Links.(fn).State;
    end
    State.LinkStates = System.Link.State;

    %% 2. Environmental & Boundary Setup
    P_atm = System.Env.P_atm;
    gamma_n2 = 1.4;
    R_n2 = 296.8;

    %% 3. Thermal Network Step
    Q_nodes = struct();
    Q_nodes.TK_N2 = 0.0;
    Q_nodes.TK_O2_01_ull = 0.0;
    Q_nodes.TK_O2_01_liq = 0.0;
    Q_nodes.TK_FU_01_ull = 0.0;
    Q_nodes.TK_FU_01_liq = 0.0;
    Q_nodes.Regen = 0.0;
    Q_nodes.SKIPPER_loss = 0.0;

    if isfield(State, 'Thermal')
        [Q_nodes, State] = StepThermalNetwork(State, dt, System.Env);
    end

    %% 4. Ensure Dynamic Node Initialization
    if ~isfield(State.Nodes, 'Inter_N2_Reg')
        V_inr = 5.0e-4;
        props_inr = FluidProperties('Nitrogen', 'From_P_T', P_atm, 293.15);
        m_inr0 = props_inr.rho * V_inr;
        State.Nodes.Inter_N2_Reg = struct(...
            'Name', 'Inter N2 Reg', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
            'V', V_inr, 'm', m_inr0, 'U', props_inr.u * m_inr0, ...
            'u', props_inr.u, 'rho', props_inr.rho, ...
            'P', P_atm, 'T', 293.15, 'h', props_inr.h, ...
            'gamma', 1.4, 'Fixed', false);
    end

    if ~isfield(State.Nodes, 'Press_Manifold')
        State.Nodes.Press_Manifold = State.Nodes.Purge_Manifold;
    end

    %% 5. Sub-Stepped Control Volume (CV) Dynamic Solver
    % Sub-stepping: N_sub = 5 provides 2.0 ms dt_sub, ensuring 100% stability
    % with CV volumes >= 0.5L while running in < 4 ms per step (ample headroom for 100 Hz real-time)
    N_sub = 5;
    dt_sub = dt / N_sub;

    % Initialize flow rate accumulators
    FlowAccum = struct();
    linkList = fieldnames(System.Links);
    for k = 1:length(linkList)
        FlowAccum.(linkList{k}) = 0.0;
    end

    for sub = 1:N_sub
        %% --- A. Pneumatics & Pressurization Network ---
        % 1. GSE Bulk N2 Fill (BV-N2-FILL)
        mdot_fill = 0.0;
        if isfield(System.Links, 'BV_N2_FILL') && isfield(State.Nodes, 'TK_N2_BULK')
            if System.Links.BV_N2_FILL.Cv > 1e-4
                P_bulk = State.Nodes.TK_N2_BULK.P;
                P_copv = State.Nodes.TK_N2.P;
                if P_bulk > P_copv + 1.0
                    [mdot_raw, ~] = CalculateLinkFlow(System.Links.BV_N2_FILL, ...
                        P_bulk, P_copv, State.Nodes.TK_N2_BULK, State.Nodes.TK_N2, dt_sub);
                    T_b = max(State.Nodes.TK_N2_BULK.T, 100);
                    T_c = max(State.Nodes.TK_N2.T, 100);
                    V_b = State.Nodes.TK_N2_BULK.V;
                    V_c = State.Nodes.TK_N2.V;
                    cap_fill = max(0.0, (P_bulk - P_copv) / (gamma_n2 * R_n2 * ((T_b / V_b) + (T_c / V_c)) * dt_sub));
                    mdot_fill = min(max(0.0, mdot_raw), cap_fill);
                end
            end
        end

        % 2. COPV Dump Valve (BV-N2-01)
        mdot_copv_dump = 0.0;
        if isfield(System.Links, 'BV_N2_01') && System.Links.BV_N2_01.Cv > 1e-6
            P_copv = State.Nodes.TK_N2.P;
            if P_copv > P_atm + 1.0
                mdot_copv_dump = CalculateLinkFlow(System.Links.BV_N2_01, P_copv, P_atm, ...
                    State.Nodes.TK_N2, State.Nodes.Atmosphere, dt_sub);
            end
        end

        % 3. Dome Regulator (REG-873D: TK_N2 -> Inter_N2_Reg)
        P_copv = State.Nodes.TK_N2.P;
        P_inr  = State.Nodes.Inter_N2_Reg.P;
        P_set  = System.Links.REG_873D.P_set;
        SPE    = System.Links.REG_873D.SPE;
        Droop  = System.Links.REG_873D.Droop;

        P_reg_target = P_set - SPE * (P_copv - P_set);
        mdot_reg = 0.0;
        if P_inr < P_reg_target && P_copv > P_inr + 1.0
            x_reg = min(1.0, max(0.0, (P_reg_target - P_inr) / Droop));
            reg_link = System.Links.REG_873D;
            reg_link.Cv = reg_link.MaxCv * x_reg;
            if reg_link.Cv > 1e-4
                mdot_reg_raw = CalculateLinkFlow(reg_link, P_copv, P_inr, ...
                    State.Nodes.TK_N2, State.Nodes.Inter_N2_Reg, dt_sub);
                % Courant cap to prevent overshooting P_reg_target
                T_inr = max(State.Nodes.Inter_N2_Reg.T, 100);
                V_inr = State.Nodes.Inter_N2_Reg.V;
                cap_reg = max(0.0, (V_inr / (gamma_n2 * R_n2 * T_inr * dt_sub)) * (P_reg_target - P_inr));
                mdot_reg = min(max(0.0, mdot_reg_raw), cap_reg);
            end
        end

        % 4. Main Vehicle N2 Isolation Valve (BV-N2-02: Inter_N2_Reg -> Purge_Manifold)
        P_pm = State.Nodes.Purge_Manifold.P;
        Cv_n2_02 = System.Links.BV_N2_02.Cv;
        mdot_iso = 0.0;
        if Cv_n2_02 > 1e-4 && P_inr > P_pm + 1.0
            mdot_iso_raw = CalculateLinkFlow(System.Links.BV_N2_02, P_inr, P_pm, ...
                State.Nodes.Inter_N2_Reg, State.Nodes.Purge_Manifold, dt_sub);
            T_pm = max(State.Nodes.Purge_Manifold.T, 100);
            V_pm = State.Nodes.Purge_Manifold.V;
            cap_iso = max(0.0, (V_pm / (gamma_n2 * R_n2 * T_pm * dt_sub)) * (P_inr - P_pm));
            mdot_iso = min(max(0.0, mdot_iso_raw), cap_iso);
        end

        % 5. Propellant Pressurization Check Valves (Press_OX, Press_FU)
        P_ox_ull = State.Nodes.TK_O2_01.Ullage.P;
        P_fu_ull = State.Nodes.TK_FU_01.Ullage.P;

        mdot_pr_ox = 0.0;
        if P_pm > P_ox_ull + 1.0
            mdot_pr_ox_raw = CalculateLinkFlow(System.Links.Press_OX, P_pm, P_ox_ull, ...
                State.Nodes.Purge_Manifold, State.Nodes.TK_O2_01.Ullage, dt_sub);
            V_ull_ox = State.Nodes.TK_O2_01.Ullage.V;
            T_ull_ox = max(State.Nodes.TK_O2_01.Ullage.T, 100);
            cap_ox = max(0.0, (V_ull_ox / (gamma_n2 * R_n2 * T_ull_ox * dt_sub)) * (P_pm - P_ox_ull));
            mdot_pr_ox = min(max(0.0, mdot_pr_ox_raw), cap_ox);
        end

        mdot_pr_fu = 0.0;
        if P_pm > P_fu_ull + 1.0
            mdot_pr_fu_raw = CalculateLinkFlow(System.Links.Press_FU, P_pm, P_fu_ull, ...
                State.Nodes.Purge_Manifold, State.Nodes.TK_FU_01.Ullage, dt_sub);
            V_ull_fu = State.Nodes.TK_FU_01.Ullage.V;
            T_ull_fu = max(State.Nodes.TK_FU_01.Ullage.T, 100);
            cap_fu = max(0.0, (V_ull_fu / (gamma_n2 * R_n2 * T_ull_fu * dt_sub)) * (P_pm - P_fu_ull));
            mdot_pr_fu = min(max(0.0, mdot_pr_fu_raw), cap_fu);
        end

        % 6. Manifold & Chamber Purge Solenoids (SV-N2-05, SV-N2-06, SV-N2-07)
        P_post_ox = State.Nodes.Post_Throttle_OX.P;
        P_post_fu = State.Nodes.Post_Throttle_FU.P;
        P_dart    = State.Nodes.DART_Chamber.P;

        mdot_pur_ox = 0.0;
        if System.Links.SV_N2_05.Cv > 1e-4 && P_pm > P_post_ox + 1.0
            mdot_pur_ox_raw = CalculateLinkFlow(System.Links.SV_N2_05, P_pm, P_post_ox, ...
                State.Nodes.Purge_Manifold, State.Nodes.Post_Throttle_OX, dt_sub);
            V_pt_ox = State.Nodes.Post_Throttle_OX.V;
            T_pt_ox = max(State.Nodes.Post_Throttle_OX.T, 100);
            cap_p_ox = max(0.0, (V_pt_ox / (gamma_n2 * R_n2 * T_pt_ox * dt_sub)) * (P_pm - P_post_ox));
            mdot_pur_ox = min(max(0.0, mdot_pur_ox_raw), cap_p_ox);
        end

        mdot_pur_fu = 0.0;
        if System.Links.SV_N2_06.Cv > 1e-4 && P_pm > P_post_fu + 1.0
            mdot_pur_fu_raw = CalculateLinkFlow(System.Links.SV_N2_06, P_pm, P_post_fu, ...
                State.Nodes.Purge_Manifold, State.Nodes.Post_Throttle_FU, dt_sub);
            V_pt_fu = State.Nodes.Post_Throttle_FU.V;
            T_pt_fu = max(State.Nodes.Post_Throttle_FU.T, 100);
            cap_p_fu = max(0.0, (V_pt_fu / (gamma_n2 * R_n2 * T_pt_fu * dt_sub)) * (P_pm - P_post_fu));
            mdot_pur_fu = min(max(0.0, mdot_pur_fu_raw), cap_p_fu);
        end

        mdot_pur_dt = 0.0;
        if System.Links.SV_N2_07.Cv > 1e-4 && P_pm > P_dart + 1.0
            mdot_pur_dt_raw = CalculateLinkFlow(System.Links.SV_N2_07, P_pm, P_dart, ...
                State.Nodes.Purge_Manifold, State.Nodes.DART_Chamber, dt_sub);
            V_dt = State.Nodes.DART_Chamber.V;
            T_dt = max(State.Nodes.DART_Chamber.T, 100);
            cap_p_dt = max(0.0, (V_dt / (gamma_n2 * R_n2 * T_dt * dt_sub)) * (P_pm - P_dart));
            mdot_pur_dt = min(max(0.0, mdot_pur_dt_raw), cap_p_dt);
        end

        % 7. RCS Cold-Gas Thrusters (SV-N2-01 .. 04: Purge_Manifold -> Atmosphere)
        mdot_rcs_1 = 0.0; mdot_rcs_2 = 0.0; mdot_rcs_3 = 0.0; mdot_rcs_4 = 0.0;
        if isfield(System.Links, 'SV_N2_01') && System.Links.SV_N2_01.Cv > 1e-4 && P_pm > P_atm + 1.0
            mdot_rcs_1 = CalculateLinkFlow(System.Links.SV_N2_01, P_pm, P_atm, ...
                State.Nodes.Purge_Manifold, State.Nodes.Atmosphere, dt_sub);
        end
        if isfield(System.Links, 'SV_N2_02') && System.Links.SV_N2_02.Cv > 1e-4 && P_pm > P_atm + 1.0
            mdot_rcs_2 = CalculateLinkFlow(System.Links.SV_N2_02, P_pm, P_atm, ...
                State.Nodes.Purge_Manifold, State.Nodes.Atmosphere, dt_sub);
        end
        if isfield(System.Links, 'SV_N2_03') && System.Links.SV_N2_03.Cv > 1e-4 && P_pm > P_atm + 1.0
            mdot_rcs_3 = CalculateLinkFlow(System.Links.SV_N2_03, P_pm, P_atm, ...
                State.Nodes.Purge_Manifold, State.Nodes.Atmosphere, dt_sub);
        end
        if isfield(System.Links, 'SV_N2_04') && System.Links.SV_N2_04.Cv > 1e-4 && P_pm > P_atm + 1.0
            mdot_rcs_4 = CalculateLinkFlow(System.Links.SV_N2_04, P_pm, P_atm, ...
                State.Nodes.Purge_Manifold, State.Nodes.Atmosphere, dt_sub);
        end

        % 8. Propellant Vents & Drains
        mdot_vent_ox = 0.0; mdot_drain_ox = 0.0;
        mdot_vent_fu = 0.0; mdot_drain_fu = 0.0;

        if isfield(System.Links, 'BV_O2_01') && System.Links.BV_O2_01.Cv > 1e-6
            if P_ox_ull > P_atm + 1.0
                mdot_vent_ox = CalculateLinkFlow(System.Links.BV_O2_01, P_ox_ull, P_atm, ...
                    State.Nodes.TK_O2_01.Ullage, State.Nodes.Atmosphere, dt_sub);
            end
        end
        if isfield(System.Links, 'BV_O2_02') && System.Links.BV_O2_02.Cv > 1e-6
            P_ox_liq = State.Nodes.TK_O2_01.Liquid.P;
            if State.Nodes.TK_O2_01.Liquid.m > 1e-3 && P_ox_liq > P_atm + 1.0
                mdot_drain_ox = CalculateLinkFlow(System.Links.BV_O2_02, P_ox_liq, P_atm, ...
                    State.Nodes.TK_O2_01.Liquid, State.Nodes.Atmosphere, dt_sub);
            end
        end

        if isfield(System.Links, 'BV_FU_01') && System.Links.BV_FU_01.Cv > 1e-6
            if P_fu_ull > P_atm + 1.0
                mdot_vent_fu = CalculateLinkFlow(System.Links.BV_FU_01, P_fu_ull, P_atm, ...
                    State.Nodes.TK_FU_01.Ullage, State.Nodes.Atmosphere, dt_sub);
            end
        end
        if isfield(System.Links, 'BV_FU_02') && System.Links.BV_FU_02.Cv > 1e-6
            P_fu_liq = State.Nodes.TK_FU_01.Liquid.P;
            if State.Nodes.TK_FU_01.Liquid.m > 1e-3 && P_fu_liq > P_atm + 1.0
                mdot_drain_fu = CalculateLinkFlow(System.Links.BV_FU_02, P_fu_liq, P_atm, ...
                    State.Nodes.TK_FU_01.Liquid, State.Nodes.Atmosphere, dt_sub);
            end
        end

        %% --- B. Propellant Feed Branches & Combustors ---
        P_skip = State.Nodes.SKIPPER.P;

        % LOX Feed Branch
        P_tank_ox  = State.Nodes.TK_O2_01.Liquid.P;
        rho_ox     = max(State.Nodes.TK_O2_01.Liquid.rho, 900.0);
        T_lox      = State.Nodes.TK_O2_01.Liquid.T;
        h_lox      = State.Nodes.TK_O2_01.Liquid.h;
        u_lox      = State.Nodes.TK_O2_01.Liquid.u;
        m_liq_ox   = State.Nodes.TK_O2_01.Liquid.m;

        Cv_main_ox = System.Links.BV_02_03.Cv;
        Cv_thrt_ox = System.Links.BV_02_04.Cv;
        K_pipe_ox  = 1.581e-4;
        K_inj_ox   = System.Links.Inj_OX.Cd * System.Links.Inj_OX.A;

        mdot_ox = 0.0;
        if Cv_main_ox > 0.01 && Cv_thrt_ox > 0.01 && m_liq_ox > 1e-3
            K_main_ox = Cv_main_ox * 2.402e-5 * 0.7071;
            K_thrt_ox = Cv_thrt_ox * 2.402e-5 * 0.7071;
            invK2_ox  = (1.0 / K_pipe_ox^2) + (1.0 / K_main_ox^2) + (1.0 / K_thrt_ox^2) + (1.0 / K_inj_ox^2);
            K_eq_ox   = 1.0 / sqrt(invK2_ox);
            DP_ox     = max(0.0, P_tank_ox - P_skip);
            mdot_ox   = K_eq_ox * sqrt(2.0 * rho_ox * DP_ox);

            P_pre_ox  = P_tank_ox - (mdot_ox^2) / (2.0 * rho_ox * K_pipe_ox^2);
            P_intr_ox = P_pre_ox  - (mdot_ox^2) / (2.0 * rho_ox * K_main_ox^2);
            P_post_ox = P_intr_ox - (mdot_ox^2) / (2.0 * rho_ox * K_thrt_ox^2);
            P_man_ox  = P_post_ox;
        else
            P_pre_ox  = P_tank_ox;
            if Cv_main_ox > 0.01
                P_intr_ox = P_tank_ox;
            else
                P_intr_ox = P_atm;
            end
            if mdot_pur_ox > 1e-5
                A_inj = max(K_inj_ox, 1e-6);
                P_man_purge = P_skip + (mdot_pur_ox / (A_inj * sqrt(gamma_n2 / (R_n2 * 293.15)) * 0.6847));
                P_post_ox = P_man_purge;
                P_man_ox  = P_man_purge;
            else
                P_post_ox = P_skip;
                P_man_ox  = P_skip;
            end
        end

        % FUEL Feed Branch (IPA with Regen Jacket)
        P_tank_fu  = State.Nodes.TK_FU_01.Liquid.P;
        rho_fu     = max(State.Nodes.TK_FU_01.Liquid.rho, 750.0);
        T_fu       = State.Nodes.TK_FU_01.Liquid.T;
        h_ipa      = State.Nodes.TK_FU_01.Liquid.h;
        u_ipa      = State.Nodes.TK_FU_01.Liquid.u;
        m_liq_fu   = State.Nodes.TK_FU_01.Liquid.m;

        Cv_main_fu = System.Links.BV_FU_03.Cv;
        Cv_thrt_fu = System.Links.BV_FU_04.Cv;
        K_pipe_fu  = 1.581e-4;
        K_inj_fu   = System.Links.Inj_FU.Cd * System.Links.Inj_FU.A;
        K_regen    = 4.67e-5;

        mdot_fu = 0.0;
        if Cv_main_fu > 0.01 && Cv_thrt_fu > 0.01 && m_liq_fu > 1e-3
            K_main_fu = Cv_main_fu * 2.402e-5 * 0.7071;
            K_thrt_fu = Cv_thrt_fu * 2.402e-5 * 0.7071;
            invK2_fu  = (1.0 / K_pipe_fu^2) + (1.0 / K_main_fu^2) + (1.0 / K_thrt_fu^2) + (1.0 / K_regen^2) + (1.0 / K_inj_fu^2);
            K_eq_fu   = 1.0 / sqrt(invK2_fu);
            DP_fu     = max(0.0, P_tank_fu - P_skip);
            mdot_fu   = K_eq_fu * sqrt(2.0 * rho_fu * DP_fu);

            P_pre_fu  = P_tank_fu - (mdot_fu^2) / (2.0 * rho_fu * K_pipe_fu^2);
            P_intr_fu = P_pre_fu  - (mdot_fu^2) / (2.0 * rho_fu * K_main_fu^2);
            P_post_fu = P_intr_fu - (mdot_fu^2) / (2.0 * rho_fu * K_thrt_fu^2);
            P_man_fu  = P_post_fu - (mdot_fu^2) / (2.0 * rho_fu * K_regen^2);
        else
            P_pre_fu  = P_tank_fu;
            if Cv_main_fu > 0.01
                P_intr_fu = P_tank_fu;
            else
                P_intr_fu = P_atm;
            end
            if mdot_pur_fu > 1e-5
                A_inj = max(K_inj_fu, 1e-6);
                P_man_purge = P_skip + (mdot_pur_fu / (A_inj * sqrt(gamma_n2 / (R_n2 * 293.15)) * 0.6847));
                P_post_fu = P_man_purge;
                P_man_fu  = P_man_purge;
            else
                P_post_fu = P_skip;
                P_man_fu  = P_skip;
            end
        end

        % DART Torch Igniter
        mdot_dart_ox = 0.0;
        mdot_dart_fu = 0.0;
        if System.Links.SV_DART_OX.Cv > 1e-6
            mdot_dart_ox = CalculateLinkFlow(System.Links.SV_DART_OX, P_intr_ox, P_dart, ...
                State.Nodes.Inter_OX, State.Nodes.DART_Chamber, dt_sub);
        end
        if System.Links.SV_DART_FU.Cv > 1e-6
            mdot_dart_fu = CalculateLinkFlow(System.Links.SV_DART_FU, P_intr_fu, P_dart, ...
                State.Nodes.Inter_FU, State.Nodes.DART_Chamber, dt_sub);
        end

        SparkActive = false;
        Inflow_DART = struct('mdot_ox', mdot_dart_ox, 'mdot_fu', mdot_dart_fu, ...
                             'mdot_n2', mdot_pur_dt, 'h_ox', h_lox, 'h_fu', h_ipa, 'h_n2', State.Nodes.Purge_Manifold.h);
        Nozzle_DART = struct('A_throat', System.Links.DART_Nozzle.A, 'Cd', System.Links.DART_Nozzle.Cd, 'P_back', P_skip);

        [State.Nodes.DART_Chamber, mdot_dart_nozzle, ~, isDartLit] = StepCombustor(...
            State.Nodes.DART_Chamber, Inflow_DART, Nozzle_DART, dt_sub, SparkActive, false, 0.0);

        % SKIPPER Main Combustor
        Inflow_SKIPPER = struct('mdot_ox', mdot_ox, 'mdot_fu', mdot_fu, ...
                                'mdot_n2', mdot_pur_ox + mdot_pur_fu, ...
                                'h_ox', h_lox, 'h_fu', h_ipa, 'h_n2', State.Nodes.Purge_Manifold.h);
        Nozzle_SKIPPER = struct('A_throat', System.Links.Nozzle.A, 'Cd', System.Links.Nozzle.Cd, 'P_back', P_atm);

        [State.Nodes.SKIPPER, mdot_nozzle, ~, ~] = StepCombustor(...
            State.Nodes.SKIPPER, Inflow_SKIPPER, Nozzle_SKIPPER, dt_sub, false, isDartLit, Q_nodes.SKIPPER_loss / N_sub);

        %% --- C. Control Volume Continuity & State Updates ---
        % 1. Ground Bulk Supply Tank (TK-N2-BULK)
        if isfield(State.Nodes, 'TK_N2_BULK')
            if mdot_fill > 0
                dm_bulk = -mdot_fill * dt_sub;
                dU_bulk = -mdot_fill * State.Nodes.TK_N2_BULK.h * dt_sub;
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

        % 2. COPV (TK-N2)
        h_in_copv = State.Nodes.TK_N2.h;
        if isfield(State.Nodes, 'TK_N2_BULK'), h_in_copv = State.Nodes.TK_N2_BULK.h; end
        dm_copv = (mdot_fill - mdot_reg - mdot_copv_dump) * dt_sub;
        dU_copv = (mdot_fill * h_in_copv - (mdot_reg + mdot_copv_dump) * State.Nodes.TK_N2.h) * dt_sub + (Q_nodes.TK_N2 / N_sub);

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

        % 3. Intermediate Regulator Node (Inter_N2_Reg)
        dm_inr = (mdot_reg - mdot_iso) * dt_sub;
        m_inr_new = max(1e-5, State.Nodes.Inter_N2_Reg.m + dm_inr);
        State.Nodes.Inter_N2_Reg.m = m_inr_new;
        State.Nodes.Inter_N2_Reg.rho = m_inr_new / State.Nodes.Inter_N2_Reg.V;
        State.Nodes.Inter_N2_Reg.P = State.Nodes.Inter_N2_Reg.rho * R_n2 * State.Nodes.Inter_N2_Reg.T;

        % 4. Main Vehicle N2 Distribution Manifold (Purge_Manifold)
        mdot_pm_out = mdot_pr_ox + mdot_pr_fu + mdot_pur_ox + mdot_pur_fu + mdot_pur_dt + ...
                      mdot_rcs_1 + mdot_rcs_2 + mdot_rcs_3 + mdot_rcs_4;
        dm_pm = (mdot_iso - mdot_pm_out) * dt_sub;
        m_pm_new = max(1e-5, State.Nodes.Purge_Manifold.m + dm_pm);
        State.Nodes.Purge_Manifold.m = m_pm_new;
        State.Nodes.Purge_Manifold.rho = m_pm_new / State.Nodes.Purge_Manifold.V;
        State.Nodes.Purge_Manifold.P = State.Nodes.Purge_Manifold.rho * R_n2 * State.Nodes.Purge_Manifold.T;

        % Mirror to Press_Manifold for backward compatibility
        State.Nodes.Press_Manifold = State.Nodes.Purge_Manifold;

        % 5. Two-Zone Stratified LOX Tank (TK-O2-01)
        dm_ull_ox = (mdot_pr_ox - mdot_vent_ox) * dt_sub;
        dU_ull_ox = (mdot_pr_ox * State.Nodes.Purge_Manifold.h - mdot_vent_ox * State.Nodes.TK_O2_01.Ullage.h) * dt_sub + (Q_nodes.TK_O2_01_ull / N_sub);
        m_ull_ox  = max(State.Nodes.TK_O2_01.Ullage.m + dm_ull_ox, 1e-4);
        U_ull_ox  = State.Nodes.TK_O2_01.Ullage.U + dU_ull_ox;

        mdot_ox_liq_out = mdot_ox + mdot_dart_ox + mdot_drain_ox;
        dm_liq_ox = -mdot_ox_liq_out * dt_sub;
        dU_liq_ox = -(mdot_ox_liq_out * h_lox) * dt_sub + (Q_nodes.TK_O2_01_liq / N_sub);
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
        State.Nodes.TK_O2_01.P          = P_tank_ox;
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

        % 6. Two-Zone Stratified Fuel Tank (TK-FU-01)
        dm_ull_fu = (mdot_pr_fu - mdot_vent_fu) * dt_sub;
        dU_ull_fu = (mdot_pr_fu * State.Nodes.Purge_Manifold.h - mdot_vent_fu * State.Nodes.TK_FU_01.Ullage.h) * dt_sub + (Q_nodes.TK_FU_01_ull / N_sub);
        m_ull_fu  = max(State.Nodes.TK_FU_01.Ullage.m + dm_ull_fu, 1e-4);
        U_ull_fu  = State.Nodes.TK_FU_01.Ullage.U + dU_ull_fu;

        mdot_fu_liq_out = mdot_fu + mdot_dart_fu + mdot_drain_fu;
        dm_liq_fu = -mdot_fu_liq_out * dt_sub;
        dU_liq_fu = -(mdot_fu_liq_out * h_ipa) * dt_sub + (Q_nodes.TK_FU_01_liq / N_sub);
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
        State.Nodes.TK_FU_01.P          = P_tank_fu;
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

        % 7. Update Feed Line Node States
        State.Nodes.Pre_Main_OX.P = P_pre_ox;
        State.Nodes.Inter_OX.P    = P_intr_ox;
        State.Nodes.Post_Throttle_OX.P = P_post_ox;
        State.Nodes.OX_Manifold.P = P_man_ox;

        State.Nodes.Pre_Main_FU.P = P_pre_fu;
        State.Nodes.Inter_FU.P    = P_intr_fu;
        State.Nodes.Post_Throttle_FU.P = P_post_fu;
        State.Nodes.FU_Manifold.P = P_man_fu;

        % --- Accumulate Flows for 10 ms Step Output ---
        FlowAccum.BV_N2_FILL   = FlowAccum.BV_N2_FILL + mdot_fill;
        FlowAccum.BV_N2_01     = FlowAccum.BV_N2_01 + mdot_copv_dump;
        FlowAccum.REG_873D     = FlowAccum.REG_873D + mdot_reg;
        FlowAccum.BV_N2_02     = FlowAccum.BV_N2_02 + mdot_iso;
        FlowAccum.Press_OX     = FlowAccum.Press_OX + mdot_pr_ox;
        FlowAccum.Press_FU     = FlowAccum.Press_FU + mdot_pr_fu;
        FlowAccum.SV_N2_05     = FlowAccum.SV_N2_05 + mdot_pur_ox;
        FlowAccum.SV_N2_06     = FlowAccum.SV_N2_06 + mdot_pur_fu;
        FlowAccum.SV_N2_07     = FlowAccum.SV_N2_07 + mdot_pur_dt;

        if isfield(FlowAccum, 'SV_N2_01'), FlowAccum.SV_N2_01 = FlowAccum.SV_N2_01 + mdot_rcs_1; end
        if isfield(FlowAccum, 'SV_N2_02'), FlowAccum.SV_N2_02 = FlowAccum.SV_N2_02 + mdot_rcs_2; end
        if isfield(FlowAccum, 'SV_N2_03'), FlowAccum.SV_N2_03 = FlowAccum.SV_N2_03 + mdot_rcs_3; end
        if isfield(FlowAccum, 'SV_N2_04'), FlowAccum.SV_N2_04 = FlowAccum.SV_N2_04 + mdot_rcs_4; end

        FlowAccum.BV_O2_01     = FlowAccum.BV_O2_01 + mdot_vent_ox;
        FlowAccum.BV_O2_02     = FlowAccum.BV_O2_02 + mdot_drain_ox;
        FlowAccum.BV_FU_01     = FlowAccum.BV_FU_01 + mdot_vent_fu;
        if isfield(FlowAccum, 'BV_FU_02'), FlowAccum.BV_FU_02 = FlowAccum.BV_FU_02 + mdot_drain_fu; end

        FlowAccum.OX_Line_1    = FlowAccum.OX_Line_1 + mdot_ox;
        FlowAccum.BV_02_03     = FlowAccum.BV_02_03 + mdot_ox;
        FlowAccum.BV_02_04     = FlowAccum.BV_02_04 + mdot_ox;
        FlowAccum.OX_Inj_Line  = FlowAccum.OX_Inj_Line + mdot_ox + mdot_pur_ox;
        FlowAccum.Inj_OX       = FlowAccum.Inj_OX + mdot_ox + mdot_pur_ox;

        FlowAccum.FU_Line_1    = FlowAccum.FU_Line_1 + mdot_fu;
        FlowAccum.BV_FU_03     = FlowAccum.BV_FU_03 + mdot_fu;
        FlowAccum.BV_FU_04     = FlowAccum.BV_FU_04 + mdot_fu;
        FlowAccum.FU_Inj_Line  = FlowAccum.FU_Inj_Line + mdot_fu + mdot_pur_fu;
        FlowAccum.Inj_FU       = FlowAccum.Inj_FU + mdot_fu + mdot_pur_fu;

        FlowAccum.SV_DART_OX   = FlowAccum.SV_DART_OX + mdot_dart_ox;
        FlowAccum.SV_DART_FU   = FlowAccum.SV_DART_FU + mdot_dart_fu;
        FlowAccum.DART_Nozzle  = FlowAccum.DART_Nozzle + mdot_dart_nozzle;
        FlowAccum.Nozzle       = FlowAccum.Nozzle + mdot_nozzle;
    end

    %% 6. Final FlowRates Normalization & Time Step Sync
    FlowRates = struct();
    accFields = fieldnames(FlowAccum);
    for k = 1:length(accFields)
        FlowRates.(accFields{k}) = FlowAccum.(accFields{k}) / N_sub;
    end

    State.Time = State.Time + dt;
    linkNames = fieldnames(System.Links);
    for k = 1:length(linkNames)
        fn = linkNames{k};
        System.Link.State.(fn) = System.Links.(fn).State;
    end
    State.LinkStates = System.Link.State;
end
