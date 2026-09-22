function [State, FlowRates, System] = StepSimulation(State, dt, System, LinkStates)
% STEPSIMULATION Universal, system-blind Control Volume network solver.
% Operates strictly on Nodes, Links, and LinkMap defined by the system builder.
% Features:
% - Vectorized link hydraulics via CalculateLinkFlow (PhysicsEngine2 equations)
% - Courant-bounded outflow mass limiter eliminating pressure chatter
% - Vectorized mass and energy accumulation via single unified LinkMap
% - Discrete (m, U) Control Volume conservation across all nodes
% - Stratified TwoZoneTank and modular pluggable Combustor stepping
% - Zero hardcoded component names; 100% C++ compileable data layout
% - Backward-compatible named struct outputs for UI, telemetry, and scripts
%
% Inputs:
%   State      - State struct (.Nodes, .Time, etc.)
%   dt         - Discrete time step [s]
%   System     - System struct (.Nodes, .Links, .LinkMap, .PortMap, etc.)
%   LinkStates - (Optional) Commanded valve states (struct, cell, or vector U)
%
% Outputs:
%   State      - Updated State struct
%   FlowRates  - Struct of link mass flow rates [kg/s]
%   System     - Updated System struct with synchronized actuator positions

    if nargin < 4, LinkStates = []; end

    %% 1. Synchronize System Links & Nodes to Flat Arrays
    [Nodes, Links, LinkMap, PortMap, NodeNames, LinkNames, System] = unpackSystemTopology(System, State);
    N_N = numel(Nodes);
    N_L = numel(Links);

    %% 2. Vectorized Actuator Dynamics (First-Order Lag Filter)
    if ~isempty(LinkStates)
        Links = applyActuatorCommands(Links, LinkNames, LinkStates, dt);
    end

    % Update active Cv = MaxCv * State for all links
    for j = 1:N_L
        if isfield(Links(j), 'MaxCv') && Links(j).MaxCv > 0
            st = 1.0;
            if isfield(Links(j), 'State'), st = Links(j).State; end
            Links(j).Cv = Links(j).MaxCv * st;
        end
    end

    %% 3. Extract Thermodynamic States for All Links via LinkMap & PortMap
    UpProps   = repmat(struct('P', 0, 'T', 293.15, 'rho', 1.18, 'gamma', 1.4, 'h', 0, 'u', 0, 'm', 0), N_L, 1);
    DownProps = UpProps;
    P_up      = zeros(N_L, 1);
    P_down    = zeros(N_L, 1);

    for j = 1:N_L
        u_idx = LinkMap(j, 1);
        d_idx = LinkMap(j, 2);
        u_port = PortMap(j, 1);
        d_port = PortMap(j, 2);

        % Fetch Upstream node/port state
        [P_up(j), UpProps(j)] = extractPortState(Nodes(u_idx), u_port);

        % Fetch Downstream node/port state
        [P_down(j), DownProps(j)] = extractPortState(Nodes(d_idx), d_port);
    end

    %% 4. Vectorized Link Flow & Heat Transfer Evaluation
    [mdot_links, isChoked_links, flowDir_links, Qdot_links] = CalculateLinkFlow(...
        Links, P_up, P_down, UpProps, DownProps, dt);

    %% 5. Courant Stability / Dynamic Outflow Mass Limiter
    % Prevents any node or port from draining more than 50% of available mass in dt
    mdot_links = applyCourantOutflowLimiter(mdot_links, Links, LinkMap, PortMap, Nodes, dt);

    %% 6. Vectorized Flux Accumulation over Single Unified LinkMap
    % Flow positive: Up -> Down (leaves Up, enters Down)
    % Flow negative: Down -> Up (leaves Down, enters Up)
    mdot_pos = max(mdot_links, 0.0);
    mdot_neg = max(-mdot_links, 0.0);

    h_up_vec   = [UpProps.h]';
    h_down_vec = [DownProps.h]';

    % Advected enthalpy and heat transfer
    Q_pos = max(Qdot_links, 0.0);
    Q_neg = max(-Qdot_links, 0.0);

    E_pos = mdot_pos .* h_up_vec   + Q_pos;
    E_neg = mdot_neg .* h_down_vec + Q_neg;

    % Accumulate fluxes per node
    u_nodes = LinkMap(:, 1);
    d_nodes = LinkMap(:, 2);

    dm_in_node  = accumarray(d_nodes, mdot_pos, [N_N, 1]) + accumarray(u_nodes, mdot_neg, [N_N, 1]);
    dm_out_node = accumarray(u_nodes, mdot_pos, [N_N, 1]) + accumarray(d_nodes, mdot_neg, [N_N, 1]);
    dE_in_node  = accumarray(d_nodes, E_pos,    [N_N, 1]) + accumarray(u_nodes, E_neg,    [N_N, 1]);
    dE_out_node = accumarray(u_nodes, E_pos,    [N_N, 1]) + accumarray(d_nodes, E_neg,    [N_N, 1]);

    % Separate port accumulation for TwoZoneTanks (Port 1 = Ullage, Port 2 = Liquid)
    [dm_in_ull, dm_out_ull, dE_in_ull, dE_out_ull, ...
     dm_in_liq, dm_out_liq, dE_in_liq, dE_out_liq] = accumulateTwoZonePortFluxes(...
        LinkMap, PortMap, mdot_pos, mdot_neg, E_pos, E_neg, N_N);

    %% 7. Discrete (m, U) Conservation Stepping for All Nodes (System-Blind)
    for i = 1:N_N
        if isfield(Nodes(i), 'Fixed') && Nodes(i).Fixed
            continue; % Boundary node (fixed P, T)
        end

        nodeType = lower(Nodes(i).Type);

        switch nodeType
            case {'gas', 'controlvolume'}
                % Single-phase compressible gas control volume
                dm = (dm_in_node(i) - dm_out_node(i)) * dt;
                dE = (dE_in_node(i) - dE_out_node(i)) * dt;

                m_new = max(1e-5, Nodes(i).m + dm);
                U_new = Nodes(i).U + dE;
                rho_new = m_new / Nodes(i).V;
                u_new   = U_new / m_new;

                props = FluidProperties(Nodes(i).Fluid, 'From_u_rho', u_new, rho_new);
                Nodes(i).m   = m_new;
                Nodes(i).U   = U_new;
                Nodes(i).u   = u_new;
                Nodes(i).rho = rho_new;
                Nodes(i).P   = props.P;
                Nodes(i).T   = props.T;
                Nodes(i).h   = props.h;
                Nodes(i).gamma = props.gamma;

            case 'twozonetank'
                % Stratified tank with Ullage gas and Liquid zones
                % Ullage zone (gas)
                dm_u = (dm_in_ull(i) - dm_out_ull(i)) * dt;
                dE_u = (dE_in_ull(i) - dE_out_ull(i)) * dt;
                m_ull_new = max(1e-4, Nodes(i).Ullage.m + dm_u);
                U_ull_new = Nodes(i).Ullage.U + dE_u;

                % Liquid zone (propellant)
                dm_l = (dm_in_liq(i) - dm_out_liq(i)) * dt;
                dE_l = (dE_in_liq(i) - dE_out_liq(i)) * dt;
                m_liq_new = max(0.0, Nodes(i).Liquid.m + dm_l);
                U_liq_new = Nodes(i).Liquid.U + dE_l;

                % Liquid volume and remaining gas volume
                u_liq_new = U_liq_new / max(m_liq_new, 1e-4);
                ref_rho = 1000.0;
                if isfield(Nodes(i).Liquid, 'rho') && Nodes(i).Liquid.rho > 100
                    ref_rho = Nodes(i).Liquid.rho;
                end
                props_liq = FluidProperties(Nodes(i).Fluid, 'From_u_rho', u_liq_new, ref_rho);
                V_liq_new = m_liq_new / max(props_liq.rho, 100.0);
                V_ull_new = max(Nodes(i).V - V_liq_new, 1e-5);

                % Ullage gas state
                rho_ull_new = m_ull_new / V_ull_new;
                u_ull_new   = U_ull_new / m_ull_new;
                ullFluid = 'Nitrogen';
                if isfield(Nodes(i).Ullage, 'Fluid'), ullFluid = Nodes(i).Ullage.Fluid; end
                props_ull = FluidProperties(ullFluid, 'From_u_rho', u_ull_new, rho_ull_new);

                P_tank = props_ull.P;
                Nodes(i).P = P_tank;

                % Update Ullage sub-struct
                Nodes(i).Ullage.m   = m_ull_new;
                Nodes(i).Ullage.U   = U_ull_new;
                Nodes(i).Ullage.V   = V_ull_new;
                Nodes(i).Ullage.P   = P_tank;
                Nodes(i).Ullage.T   = props_ull.T;
                Nodes(i).Ullage.rho = rho_ull_new;
                Nodes(i).Ullage.u   = u_ull_new;
                Nodes(i).Ullage.h   = props_ull.h;

                % Update Liquid sub-struct
                Nodes(i).Liquid.m   = m_liq_new;
                Nodes(i).Liquid.U   = U_liq_new;
                Nodes(i).Liquid.V   = V_liq_new;
                Nodes(i).Liquid.P   = P_tank;
                Nodes(i).Liquid.T   = props_liq.T;
                Nodes(i).Liquid.rho = props_liq.rho;
                Nodes(i).Liquid.u   = u_liq_new;
                Nodes(i).Liquid.h   = props_liq.h;

            case {'line', 'liquidline'}
                % Universal Control Volume Line
                % 1. Incompressible liquid hydraulic line (with dead-end lockup and flow-through)
                % 2. Pneumatic gas purge displacement (N2 entering liquid line, expelling liquid)
                % 3. Pure compressible gas flow
                
                dm_out = dm_out_node(i);
                dE_in  = dE_in_node(i);
                dE_out = dE_out_node(i);
                
                % Reference liquid density
                rho_l_ref = 1140.0;
                if contains(lower(Nodes(i).Fluid), 'ipa') || contains(lower(Nodes(i).Fluid), 'fu')
                    rho_l_ref = 786.0;
                end
                
                % Sub-mass tracking (Liquid vs Gas)
                if ~isfield(Nodes(i), 'm_liq') || isempty(Nodes(i).m_liq)
                    if strcmpi(Nodes(i).Fluid, 'nitrogen') || Nodes(i).rho < 600.0
                        Nodes(i).m_liq = 0.0;
                        Nodes(i).m_gas = Nodes(i).m;
                    else
                        Nodes(i).m_liq = Nodes(i).m;
                        Nodes(i).m_gas = 1e-5;
                    end
                end
                
                % Incoming flows: distinguish gas purge from liquid flow
                inLinks = find(LinkMap(:, 2) == i);
                dm_gas_in = 0.0;
                dm_liq_in = 0.0;
                for k = 1:numel(inLinks)
                    l_idx = inLinks(k);
                    md = mdot_pos(l_idx);
                    if md > 1e-8
                        if UpProps(l_idx).rho < 600.0
                            dm_gas_in = dm_gas_in + md * dt;
                        else
                            dm_liq_in = dm_liq_in + md * dt;
                        end
                    end
                end
                
                % Outgoing flow: liquid leaves first, then gas
                dm_tot_out = dm_out * dt;
                if Nodes(i).m_liq > 1e-4
                    dm_liq_out = min(Nodes(i).m_liq, dm_tot_out);
                    dm_gas_out = dm_tot_out - dm_liq_out;
                else
                    dm_liq_out = 0.0;
                    dm_gas_out = dm_tot_out;
                end
                
                m_liq_new = max(0.0, Nodes(i).m_liq + dm_liq_in - dm_liq_out);
                m_gas_new = max(1e-5, Nodes(i).m_gas + dm_gas_in - dm_gas_out);
                m_tot_new = m_liq_new + m_gas_new;
                
                U_new = max(100.0, Nodes(i).U + (dE_in - dE_out) * dt);
                u_new = U_new / m_tot_new;
                
                % Liquid volume and gas volume
                V_liq = m_liq_new / rho_l_ref;
                V_gas = max(1e-7, Nodes(i).V - V_liq);
                
                % State evaluation based on phase state:
                if m_liq_new < 1e-3
                    % Pure gas line (purged or N2 line)
                    props = FluidProperties('Nitrogen', 'From_u_rho', u_new, m_gas_new / Nodes(i).V);
                    P_line = props.P;
                    T_line = props.T;
                    rho_bulk = m_gas_new / Nodes(i).V;
                    gamma_line = props.gamma;
                elseif dm_gas_in > 1e-8 || m_gas_new > 1e-3
                    % Active gas purge displacing liquid:
                    R_gas = 296.8;
                    T_gas = max(90.0, min(350.0, Nodes(i).T));
                    P_line = (m_gas_new * R_gas * T_gas) / V_gas;
                    P_line = max(101325.0, min(5e7, P_line));
                    rho_bulk = m_tot_new / Nodes(i).V;
                    T_line = T_gas;
                    gamma_line = 1.4;
                else
                    % Liquid-dominated line: unconditionally stable hydraulic nodal relaxation
                    % Balances inflow and outflow conductances (Kirchhoff's nodal law)
                    sum_C_P = 0.0;
                    sum_C   = 0.0;
                    
                    % Incoming links to node i
                    for k = 1:numel(inLinks)
                        l_idx = inLinks(k);
                        L = Links(l_idx);
                        if strcmpi(L.Type, 'thermal') || strcmpi(L.Type, 'signal'), continue; end
                        P_adj = P_up(l_idx);
                        
                        % Link hydraulic capacity K
                        K_lnk = 0.0;
                        if isfield(L, 'MaxCv') && L.MaxCv > 0
                            st = 1.0;
                            if isfield(L, 'State'), st = L.State; end
                            K_lnk = L.MaxCv * st * 2.402e-5 * sqrt(rho_l_ref);
                        elseif isfield(L, 'A') && L.A > 0
                            if isfield(L, 'Zeta') && L.Zeta > 0
                                K_lnk = L.A * sqrt(2.0 * rho_l_ref / max(L.Zeta, 1.0));
                            else
                                cd = 0.70;
                                if isfield(L, 'Cd') && L.Cd > 0, cd = L.Cd; end
                                K_lnk = L.A * cd * sqrt(2.0 * rho_l_ref);
                            end
                        end
                        
                        % Check valve masking: reverse flow blocked
                        if isfield(L, 'Type') && strcmpi(L.Type, 'check') && (P_adj <= Nodes(i).P)
                            K_lnk = 0.0;
                        end
                        
                        if K_lnk > 0
                            C_lnk = K_lnk / sqrt(max(abs(P_adj - Nodes(i).P), 5000.0));
                            sum_C_P = sum_C_P + C_lnk * P_adj;
                            sum_C   = sum_C   + C_lnk;
                        end
                    end
                    
                    % Outgoing links from node i
                    outLinks = find(LinkMap(:, 1) == i);
                    for k = 1:numel(outLinks)
                        l_idx = outLinks(k);
                        L = Links(l_idx);
                        if strcmpi(L.Type, 'thermal') || strcmpi(L.Type, 'signal'), continue; end
                        P_adj = P_down(l_idx);
                        
                        % Link hydraulic capacity K
                        K_lnk = 0.0;
                        if isfield(L, 'MaxCv') && L.MaxCv > 0
                            st = 1.0;
                            if isfield(L, 'State'), st = L.State; end
                            K_lnk = L.MaxCv * st * 2.402e-5 * sqrt(rho_l_ref);
                        elseif isfield(L, 'A') && L.A > 0
                            if isfield(L, 'Zeta') && L.Zeta > 0
                                K_lnk = L.A * sqrt(2.0 * rho_l_ref / max(L.Zeta, 1.0));
                            else
                                cd = 0.70;
                                if isfield(L, 'Cd') && L.Cd > 0, cd = L.Cd; end
                                K_lnk = L.A * cd * sqrt(2.0 * rho_l_ref);
                            end
                        end
                        
                        % Check valve masking: reverse flow blocked
                        if isfield(L, 'Type') && strcmpi(L.Type, 'check') && (Nodes(i).P <= P_adj)
                            K_lnk = 0.0;
                        end
                        
                        if K_lnk > 0
                            C_lnk = K_lnk / sqrt(max(abs(Nodes(i).P - P_adj), 5000.0));
                            sum_C_P = sum_C_P + C_lnk * P_adj;
                            sum_C   = sum_C   + C_lnk;
                        end
                    end
                    
                    if sum_C > 1e-12
                        P_target = sum_C_P / sum_C;
                        % Damped relaxation toward hydraulic equilibrium (eliminates acoustic stiffness)
                        P_line = Nodes(i).P + 0.8 * (P_target - Nodes(i).P);
                    else
                        P_line = Nodes(i).P;
                    end
                    P_line = max(101325.0, min(5e7, P_line));
                    m_liq_new = rho_l_ref * Nodes(i).V;
                    m_tot_new = m_liq_new + m_gas_new;
                    rho_bulk = rho_l_ref;
                    T_line = Nodes(i).T;
                    gamma_line = 1.15;
                end
                
                Nodes(i).m_liq = m_liq_new;
                Nodes(i).m_gas = m_gas_new;
                Nodes(i).m     = m_tot_new;
                Nodes(i).U     = U_new;
                Nodes(i).u     = u_new;
                Nodes(i).rho   = rho_bulk;
                Nodes(i).P     = P_line;
                Nodes(i).T     = T_line;
                Nodes(i).h     = u_new + P_line / rho_bulk;
                Nodes(i).gamma = gamma_line;

            case 'combustor'
                % Pluggable combustion chamber node
                % Identify all incoming links to this combustor
                inLinks = find(LinkMap(:, 2) == i);
                Inflow = struct('mdot_ox', 0.0, 'mdot_fu', 0.0, 'mdot_n2', 0.0, ...
                                'h_ox', 0.0, 'h_fu', 0.0, 'h_n2', 0.0);

                for k = 1:numel(inLinks)
                    l_idx = inLinks(k);
                    md = mdot_links(l_idx);
                    if md > 1e-6
                        src_node = Nodes(LinkMap(l_idx, 1));
                        src_fluid = lower(src_node.Fluid);
                        if contains(src_fluid, 'ox') || contains(src_fluid, 'oxygen')
                            Inflow.mdot_ox = Inflow.mdot_ox + md;
                            Inflow.h_ox    = UpProps(l_idx).h;
                        elseif contains(src_fluid, 'fu') || contains(src_fluid, 'ipa') || contains(src_fluid, 'methane')
                            Inflow.mdot_fu = Inflow.mdot_fu + md;
                            Inflow.h_fu    = UpProps(l_idx).h;
                        else
                            Inflow.mdot_n2 = Inflow.mdot_n2 + md;
                            Inflow.h_n2    = UpProps(l_idx).h;
                        end
                    end
                end

                % Identify outgoing nozzle link
                outLinks = find(LinkMap(:, 1) == i);
                Nozzle = struct('A_throat', 1e-4, 'Cd', 0.95, 'P_back', 101325);
                nozzle_link_idx = 0;
                for k = 1:numel(outLinks)
                    l_idx = outLinks(k);
                    if isfield(Links(l_idx), 'A') && Links(l_idx).A > 0
                        Nozzle.A_throat = Links(l_idx).A;
                        if isfield(Links(l_idx), 'Cd') && Links(l_idx).Cd > 0
                            Nozzle.Cd = Links(l_idx).Cd;
                        end
                        Nozzle.P_back = P_down(l_idx);
                        nozzle_link_idx = l_idx;
                        break;
                    end
                end

                % Ignition triggers
                SparkActive = false;
                if isfield(Nodes(i), 'SparkActive') && Nodes(i).SparkActive
                    SparkActive = true;
                end
                for j_sig = 1:N_L
                    if strcmpi(Links(j_sig).Type, 'signal') && LinkMap(j_sig, 2) == i
                        if Links(j_sig).Cv > 0.5 || Links(j_sig).State > 0.5
                            SparkActive = true;
                        end
                    end
                end

                TorchActive = false;
                if Inflow.mdot_n2 > 0.001 && isfield(Nodes(i), 'TorchLit') && Nodes(i).TorchLit
                    TorchActive = true;
                end

                % Advance Combustor
                [Nodes(i), mdot_nozzle, ~, ~] = StepCombustor(...
                    Nodes(i), Inflow, Nozzle, dt, SparkActive, TorchActive, 0.0);

                % Update nozzle link flow if present
                if nozzle_link_idx > 0
                    mdot_links(nozzle_link_idx) = mdot_nozzle;
                end
        end
    end

    %% 8. Pack Synchronized Output Structs (Backward Compatibility)
    FlowRates = struct();
    for j = 1:N_L
        fn = LinkNames{j};
        is_therm = isfield(Links(j), 'Type') && strcmpi(Links(j).Type, 'thermal');
        if is_therm
            FlowRates.(fn) = Qdot_links(j);
        else
            FlowRates.(fn) = mdot_links(j);
        end
        Links(j).Flow = mdot_links(j);
        Links(j).Q    = Qdot_links(j);
        Links(j).IsChoked = isChoked_links(j);
        Links(j).FlowDir = flowDir_links(j);
        System.Links.(fn) = Links(j);
        System.Link.State.(fn) = Links(j).State;
    end

    for i = 1:N_N
        fn = NodeNames{i};
        State.Nodes.(fn) = Nodes(i);
    end
    State.NodeArray = Nodes;
    System.NodeArray = Nodes;
    System.LinkArray = Links;
    State.LinkStates = System.Link.State;
    State.Time = State.Time + dt;
end

%% --- Helper: Unpack System Topology into Fast Numerical Arrays ---
function [Nodes, Links, LinkMap, PortMap, NodeNames, LinkNames, System] = unpackSystemTopology(System, State)
    % Extract or build NodeArray
    if isfield(State, 'NodeArray') && ~isempty(State.NodeArray)
        Nodes = State.NodeArray;
        NodeNames = {Nodes.Name};
    elseif isfield(System, 'NodeArray') && ~isempty(System.NodeArray)
        Nodes = System.NodeArray;
        NodeNames = {Nodes.Name};
    elseif isfield(State, 'Nodes') && isstruct(State.Nodes)
        NodeNames = fieldnames(State.Nodes);
        N_N = numel(NodeNames);
        Nodes = repmat(State.Nodes.(NodeNames{1}), N_N, 1);
        for i = 1:N_N
            Nodes(i) = State.Nodes.(NodeNames{i});
            if ~isfield(Nodes(i), 'ID') || isempty(Nodes(i).ID)
                Nodes(i).ID = i;
            end
        end
        System.NodeArray = Nodes;
    else
        NodeNames = fieldnames(System.Nodes);
        N_N = numel(NodeNames);
        Nodes = repmat(System.Nodes.(NodeNames{1}), N_N, 1);
        for i = 1:N_N
            Nodes(i) = System.Nodes.(NodeNames{i});
            if ~isfield(Nodes(i), 'ID') || isempty(Nodes(i), 'ID')
                Nodes(i).ID = i;
            end
        end
        System.NodeArray = Nodes;
    end

    % Unconditionally synchronize any overrides/changes from State.Nodes
    if isfield(State, 'Nodes') && isstruct(State.Nodes)
        for i = 1:numel(Nodes)
            fn = NodeNames{i};
            if isfield(State.Nodes, fn)
                src = State.Nodes.(fn);
                src_f = fieldnames(src);
                for k = 1:numel(src_f)
                    fk = src_f{k};
                    if isfield(Nodes(i), fk)
                        Nodes(i).(fk) = src.(fk);
                    end
                end
            end
        end
    end

    % Extract or build LinkArray
    if isfield(System, 'LinkArray') && ~isempty(System.LinkArray)
        Links = System.LinkArray;
        LinkNames = {Links.Name};
    else
        LinkNames = fieldnames(System.Links);
        N_L = numel(LinkNames);
        Links = repmat(System.Links.(LinkNames{1}), N_L, 1);
        for j = 1:N_L
            Links(j) = System.Links.(LinkNames{j});
            if ~isfield(Links(j), 'ID') || isempty(Links(j).ID)
                Links(j).ID = j;
            end
        end
        System.LinkArray = Links;
    end

    % Extract or build LinkMap and PortMap
    if isfield(System, 'LinkMap') && ~isempty(System.LinkMap)
        LinkMap = System.LinkMap;
    else
        N_L = numel(Links);
        LinkMap = zeros(N_L, 2);
        for j = 1:N_L
            u_id = resolveNodeID(Links(j).Up, Nodes, NodeNames);
            d_id = resolveNodeID(Links(j).Down, Nodes, NodeNames);
            LinkMap(j, :) = [u_id, d_id];
        end
        System.LinkMap = LinkMap;
    end

    if isfield(System, 'PortMap') && ~isempty(System.PortMap)
        PortMap = System.PortMap;
    else
        N_L = numel(Links);
        PortMap = zeros(N_L, 2);
        for j = 1:N_L
            u_p = 0;
            if isfield(Links(j), 'UpPort'), u_p = parsePortCode(Links(j).UpPort); end
            d_p = 0;
            if isfield(Links(j), 'DownPort'), d_p = parsePortCode(Links(j).DownPort); end
            PortMap(j, :) = [u_p, d_p];
        end
        System.PortMap = PortMap;
    end
end

%% --- Helper: Resolve Node Identifier to Integer Index ---
function id = resolveNodeID(ref, Nodes, NodeNames)
    if isnumeric(ref)
        id = ref;
    elseif ischar(ref) || isstring(ref)
        id = find(strcmp(NodeNames, ref), 1);
        if isempty(id)
            % Try searching by .Name field inside Nodes
            for k = 1:numel(Nodes)
                if strcmp(Nodes(k).Name, ref)
                    id = k;
                    return;
                end
            end
            error('resolveNodeID: Unknown node reference "%s"', ref);
        end
    else
        id = 1;
    end
end

%% --- Helper: Parse Port Selector (0 = Bulk, 1 = Ullage, 2 = Liquid) ---
function code = parsePortCode(portStr)
    if isnumeric(portStr)
        code = portStr;
    elseif ischar(portStr) || isstring(portStr)
        pLower = lower(strtrim(portStr));
        if strcmp(pLower, 'ullage') || strcmp(pLower, 'gas') || strcmp(pLower, 'top')
            code = 1;
        elseif strcmp(pLower, 'liquid') || strcmp(pLower, 'bottom')
            code = 2;
        else
            code = 0;
        end
    else
        code = 0;
    end
end

%% --- Helper: Extract Port State from a Node ---
function [P, Props] = extractPortState(Node, portCode)
    if strcmpi(Node.Type, 'twozonetank')
        if portCode == 1
            src = Node.Ullage;
        else
            src = Node.Liquid;
        end
    else
        src = Node;
    end
    P = src.P;
    Props.P = src.P;
    Props.T = src.T;
    Props.rho = src.rho;
    Props.gamma = 1.4;
    if isfield(src, 'gamma') && src.gamma > 1.0, Props.gamma = src.gamma; end
    Props.h = src.h;
    Props.u = src.u;
    Props.m = src.m;
end

%% --- Helper: Apply Actuator Commands with Lag Filter ---
function Links = applyActuatorCommands(Links, LinkNames, LinkStates, dt)
    N_L = numel(Links);
    if isstruct(LinkStates) && isscalar(LinkStates)
        fn_cmd = fieldnames(LinkStates);
        for k = 1:numel(fn_cmd)
            raw = fn_cmd{k};
            fn = strrep(raw, '-', '_');
            j = find(strcmp(LinkNames, fn) | strcmp(LinkNames, raw), 1);
            if ~isempty(j)
                val = LinkStates.(raw);
                Links(j) = filterActuatorState(Links(j), val, dt);
            end
        end
    elseif isnumeric(LinkStates)
        nCmd = min(numel(LinkStates), N_L);
        for j = 1:nCmd
            Links(j) = filterActuatorState(Links(j), LinkStates(j), dt);
        end
    end
end

%% --- Helper: Filter Single Actuator Position ---
function L = filterActuatorState(L, targetVal, dt)
    if isfield(L, 'MaxCv') && L.MaxCv > 0 && targetVal > 1.0
        target = max(0.0, min(1.0, targetVal / L.MaxCv));
    else
        target = max(0.0, min(1.0, targetVal));
    end

    if isfield(L, 'Tau') && L.Tau > dt
        alpha = 1.0 - exp(-dt / L.Tau);
        cur = 0.0;
        if isfield(L, 'State'), cur = L.State; end
        new_st = cur + alpha * (target - cur);
        if abs(new_st - target) < 1e-4, new_st = target; end
        L.State = new_st;
    else
        L.State = target;
    end
end

%% --- Helper: Courant Outflow Limiter ---
function mdot = applyCourantOutflowLimiter(mdot, Links, LinkMap, PortMap, Nodes, dt)
    N_L = numel(Links);
    N_N = numel(Nodes);

    % Separate demand for standard nodes, ullage ports, and liquid ports
    demand_node = zeros(N_N, 1);
    demand_ull  = zeros(N_N, 1);
    demand_liq  = zeros(N_N, 1);

    for j = 1:N_L
        md = mdot(j);
        if md > 0
            u = LinkMap(j, 1);
            u_p = PortMap(j, 1);
            if u_p == 1, demand_ull(u) = demand_ull(u) + md;
            elseif u_p == 2, demand_liq(u) = demand_liq(u) + md;
            else, demand_node(u) = demand_node(u) + md;
            end
        elseif md < 0
            d = LinkMap(j, 2);
            d_p = PortMap(j, 2);
            if d_p == 1, demand_ull(d) = demand_ull(d) + abs(md);
            elseif d_p == 2, demand_liq(d) = demand_liq(d) + abs(md);
            else, demand_node(d) = demand_node(d) + abs(md);
            end
        end
    end

    % Compute scale factors
    scale_node = ones(N_N, 1);
    scale_ull  = ones(N_N, 1);
    scale_liq  = ones(N_N, 1);

    for i = 1:N_N
        if isfield(Nodes(i), 'Fixed') && Nodes(i).Fixed, continue; end
        if strcmpi(Nodes(i).Type, 'twozonetank')
            % Ullage
            avail_u = 0.5 * Nodes(i).Ullage.m / max(dt, 1e-4);
            if demand_ull(i) > avail_u && demand_ull(i) > 0
                scale_ull(i) = avail_u / demand_ull(i);
            end
            % Liquid
            avail_l = 0.5 * Nodes(i).Liquid.m / max(dt, 1e-4);
            if demand_liq(i) > avail_l && demand_liq(i) > 0
                scale_liq(i) = avail_l / demand_liq(i);
            end
        else
            avail_m = 0.5 * Nodes(i).m / max(dt, 1e-4);
            if demand_node(i) > avail_m && demand_node(i) > 0
                scale_node(i) = avail_m / demand_node(i);
            end
        end
    end

    % Apply scale factor to links
    for j = 1:N_L
        md = mdot(j);
        if md > 0
            u = LinkMap(j, 1);
            u_p = PortMap(j, 1);
            if u_p == 1, s = scale_ull(u);
            elseif u_p == 2, s = scale_liq(u);
            else, s = scale_node(u);
            end
            mdot(j) = md * s;
        elseif md < 0
            d = LinkMap(j, 2);
            d_p = PortMap(j, 2);
            if d_p == 1, s = scale_ull(d);
            elseif d_p == 2, s = scale_liq(d);
            else, s = scale_node(d);
            end
            mdot(j) = md * s;
        end
    end
end

%% --- Helper: Accumulate Port Fluxes for TwoZoneTanks ---
function [dm_in_ull, dm_out_ull, dE_in_ull, dE_out_ull, ...
          dm_in_liq, dm_out_liq, dE_in_liq, dE_out_liq] = accumulateTwoZonePortFluxes(...
          LinkMap, PortMap, mdot_pos, mdot_neg, E_pos, E_neg, N_N)

    dm_in_ull  = zeros(N_N, 1);
    dm_out_ull = zeros(N_N, 1);
    dE_in_ull  = zeros(N_N, 1);
    dE_out_ull = zeros(N_N, 1);

    dm_in_liq  = zeros(N_N, 1);
    dm_out_liq = zeros(N_N, 1);
    dE_in_liq  = zeros(N_N, 1);
    dE_out_liq = zeros(N_N, 1);

    N_L = size(LinkMap, 1);
    for j = 1:N_L
        u = LinkMap(j, 1);
        d = LinkMap(j, 2);
        u_p = PortMap(j, 1);
        d_p = PortMap(j, 2);

        % Outflow from Upstream (forward flow leaves Up)
        if mdot_pos(j) > 0 || E_pos(j) > 0
            if u_p == 1
                dm_out_ull(u) = dm_out_ull(u) + mdot_pos(j);
                dE_out_ull(u) = dE_out_ull(u) + E_pos(j);
            elseif u_p == 2
                dm_out_liq(u) = dm_out_liq(u) + mdot_pos(j);
                dE_out_liq(u) = dE_out_liq(u) + E_pos(j);
            end
        end

        % Inflow into Downstream (forward flow enters Down)
        if mdot_pos(j) > 0 || E_pos(j) > 0
            if d_p == 1
                dm_in_ull(d) = dm_in_ull(d) + mdot_pos(j);
                dE_in_ull(d) = dE_in_ull(d) + E_pos(j);
            elseif d_p == 2
                dm_in_liq(d) = dm_in_liq(d) + mdot_pos(j);
                dE_in_liq(d) = dE_in_liq(d) + E_pos(j);
            end
        end

        % Outflow from Downstream (reverse flow leaves Down)
        if mdot_neg(j) > 0 || E_neg(j) > 0
            if d_p == 1
                dm_out_ull(d) = dm_out_ull(d) + mdot_neg(j);
                dE_out_ull(d) = dE_out_ull(d) + E_neg(j);
            elseif d_p == 2
                dm_out_liq(d) = dm_out_liq(d) + mdot_neg(j);
                dE_out_liq(d) = dE_out_liq(d) + E_neg(j);
            end
        end

        % Inflow into Upstream (reverse flow enters Up)
        if mdot_neg(j) > 0 || E_neg(j) > 0
            if u_p == 1
                dm_in_ull(u) = dm_in_ull(u) + mdot_neg(j);
                dE_in_ull(u) = dE_in_ull(u) + E_neg(j);
            elseif u_p == 2
                dm_in_liq(u) = dm_in_liq(u) + mdot_neg(j);
                dE_in_liq(u) = dE_in_liq(u) + E_neg(j);
            end
        end
    end
end
