function [System, InitialState] = SKIPPER_StandA()
%% SKIPPER_StandA.m: Test Stand A Fluid Topology for SKIPPER Hot-Fire
% Fixed 550 psi supply LOX and IPA run tanks feeding SKIPPER combustor node
% via variable throttle valves, line drop orifices, and injector orifices.
% Outflow exits through choked sonic nozzle into atmospheric boundary sink.
%
% Topology:
%   Nodes:
%     1. TK_OX          - Fixed Boundary LOX at 550 psi (90.17 K)
%     2. TK_FU          - Fixed Boundary IPA at 550 psi (293.15 K)
%     3. OX_Valve_Down  - Line/manifold volume downstream of OX throttle (0.5 L)
%     4. OX_Manifold    - Line/manifold volume before OX injector (0.5 L)
%     5. FU_Valve_Down  - Line/manifold volume downstream of FU throttle (0.5 L)
%     6. FU_Regen_Inlet - Line/manifold volume at regen jacket inlet (0.3 L)
%     7. FU_Manifold    - Line/manifold volume before FU injector (0.5 L)
%     8. SKIPPER        - Main combustor chamber node (1.8 L, cstar = 1550 m/s)
%     9. Atmosphere     - Ambient boundary sink (101,325 Pa, 293.15 K)
%
%   Links:
%     1. Valve_OX       - Main throttle valve on LOX feed (MaxCv = 2.5)
%     2. OX_Line        - Line pressure drop between valve and manifold (Cv = 8.0)
%     3. Inj_OX         - Injector orifice into SKIPPER (A = 3.25e-5 m^2, Cd = 0.70)
%     4. Valve_FU       - Main throttle valve on Fuel feed (MaxCv = 2.5)
%     5. FU_Line        - Line pressure drop between valve and regen (Cv = 10.0)
%     6. Regen_Line     - Regen cooling jacket pressure drop (Cv = 4.5)
%     7. Inj_FU         - Injector orifice into SKIPPER (A = 3.65e-5 m^2, Cd = 0.75)
%     8. Nozzle         - Sonic throat nozzle into atmosphere (A = 7.14e-4 m^2, Cd = 0.98)
%
% Standardized for StepSimulation universal lumped-parameter solver.

    psi2Pa = 6894.757;
    P_atm  = 101325.0;
    P_sup  = 550.0 * psi2Pa;

    %% 1. Node Geometric & Baseline State Definitions
    State = struct();
    State.Nodes = struct();

    % Node 1: TK_OX (Fixed Boundary LOX at 550 psi, 90.17 K)
    V_tk_ox = 1.0; % Arbitrary boundary reservoir volume
    props_tk_ox = FluidProperties('Oxygen', 'From_P_T', P_sup, 90.17);
    m_tk_ox0 = props_tk_ox.rho * V_tk_ox;
    State.Nodes.TK_OX = struct(...
        'ID', 1, 'Name', 'TK_OX', 'Type', 'Boundary', 'Fluid', 'Oxygen', ...
        'V', V_tk_ox, 'm', m_tk_ox0, 'U', props_tk_ox.u * m_tk_ox0, ...
        'u', props_tk_ox.u, 'rho', props_tk_ox.rho, ...
        'P', P_sup, 'T', 90.17, 'h', props_tk_ox.h, ...
        'gamma', 1.1, 'Fixed', true, 'm_liq', m_tk_ox0, 'm_gas', 1e-5);

    % Node 2: TK_FU (Fixed Boundary IPA at 550 psi, 293.15 K)
    V_tk_fu = 1.0; % Arbitrary boundary reservoir volume
    props_tk_fu = FluidProperties('IPA', 'From_P_T', P_sup, 293.15);
    m_tk_fu0 = props_tk_fu.rho * V_tk_fu;
    State.Nodes.TK_FU = struct(...
        'ID', 2, 'Name', 'TK_FU', 'Type', 'Boundary', 'Fluid', 'IPA', ...
        'V', V_tk_fu, 'm', m_tk_fu0, 'U', props_tk_fu.u * m_tk_fu0, ...
        'u', props_tk_fu.u, 'rho', props_tk_fu.rho, ...
        'P', P_sup, 'T', 293.15, 'h', props_tk_fu.h, ...
        'gamma', 1.15, 'Fixed', true, 'm_liq', m_tk_fu0, 'm_gas', 1e-5);

    % Node 3: OX_Manifold (Line volume before OX injector, V = 0.1 L, Fluid = 'Oxygen')
    V_ox_man = 0.1e-3; % 100 mL
    props_ox_init = FluidProperties('Oxygen', 'From_P_T', P_atm, 90.17);
    m_ox_man0 = props_ox_init.rho * V_ox_man;
    State.Nodes.OX_Manifold = struct(...
        'ID', 3, 'Name', 'OX_Manifold', 'Type', 'Line', 'Fluid', 'Oxygen', ...
        'V', V_ox_man, 'm', m_ox_man0, 'U', props_ox_init.u * m_ox_man0, ...
        'u', props_ox_init.u, 'rho', props_ox_init.rho, ...
        'P', P_atm, 'T', 90.17, 'h', props_ox_init.h, ...
        'gamma', 1.1, 'Fixed', false, 'm_liq', m_ox_man0, 'm_gas', 1e-5);

    % Node 4: FU_Regen_Inlet (Line volume at regen jacket inlet, V = 0.1 L, Fluid = 'IPA')
    V_fu_reg = 0.1e-3; % 100 mL
    props_fu_init = FluidProperties('IPA', 'From_P_T', P_atm, 293.15);
    m_fu_reg0 = props_fu_init.rho * V_fu_reg;
    State.Nodes.FU_Regen_Inlet = struct(...
        'ID', 4, 'Name', 'FU_Regen_Inlet', 'Type', 'Line', 'Fluid', 'IPA', ...
        'V', V_fu_reg, 'm', m_fu_reg0, 'U', props_fu_init.u * m_fu_reg0, ...
        'u', props_fu_init.u, 'rho', props_fu_init.rho, ...
        'P', P_atm, 'T', 293.15, 'h', props_fu_init.h, ...
        'gamma', 1.15, 'Fixed', false, 'm_liq', m_fu_reg0, 'm_gas', 1e-5);

    % Node 5: FU_Manifold (Line volume before FU injector, V = 0.1 L, Fluid = 'IPA')
    V_fu_man = 2.72e-5; 
    m_fu_man0 = props_fu_init.rho * V_fu_man;
    State.Nodes.FU_Manifold = struct(...
        'ID', 5, 'Name', 'FU_Manifold', 'Type', 'Line', 'Fluid', 'IPA', ...
        'V', V_fu_man, 'm', m_fu_man0, 'U', props_fu_init.u * m_fu_man0, ...
        'u', props_fu_init.u, 'rho', props_fu_init.rho, ...
        'P', P_atm, 'T', 293.15, 'h', props_fu_init.h, ...
        'gamma', 1.15, 'Fixed', false, 'm_liq', m_fu_man0, 'm_gas', 1e-5);

    % Node 6: SKIPPER (Combustor node, V = 0.0018 m^3, Fluid = 'CombustionGas', cstar = 1550)
    V_skip = 0.001302;
    props_gas_init = FluidProperties('Nitrogen', 'From_P_T', P_atm, 293.15);
    m_skip0 = props_gas_init.rho * V_skip;
    State.Nodes.SKIPPER = struct(...
        'ID', 6, 'Name', 'SKIPPER', 'Type', 'Combustor', 'Fluid', 'CombustionGas', ...
        'V', V_skip, 'm', m_skip0, 'U', props_gas_init.u * m_skip0, ...
        'u', props_gas_init.u, 'rho', props_gas_init.rho, ...
        'P', P_atm, 'T', 293.15, 'h', props_gas_init.h, ...
        'gamma', 1.21, 'cstar', 1550.0, 'isLit', false, ...
        'SparkActive', false, 'TorchLit', false, 'Fixed', false, ...
        'm_liq', 0.0, 'm_gas', m_skip0);

    % Node 7: Atmosphere (Fixed Boundary at 101325 Pa, 293.15 K)
    State.Nodes.Atmosphere = struct(...
        'ID', 7, 'Name', 'Atmosphere', 'Type', 'Boundary', 'Fluid', 'Nitrogen', ...
        'V', 1.0, 'm', 1.18, 'U', 2e5, 'u', 1.7e5, 'rho', 1.18, ...
        'P', P_atm, 'T', 293.15, 'h', 2.8e5, 'gamma', 1.4, ...
        'Fixed', true, 'm_liq', 0.0, 'm_gas', 1.18);

    %% 2. Link Definitions
    Links = struct();

    % Link 1: Valve_OX (Throttle valve: TK_OX -> OX_Manifold)
    Links.Valve_OX = struct('ID', 1, 'Name', 'Valve_OX', 'Type', 'Throttle', ...
        'Up', 'TK_OX', 'Down', 'OX_Manifold', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'MaxCv', 2.5, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05, 'IsOneWay', true);

    % Link 2: Inj_OX (Injector orifice: OX_Manifold -> SKIPPER)
    Links.Inj_OX = struct('ID', 2, 'Name', 'Inj_OX', 'Type', 'Orifice', ...
        'Up', 'OX_Manifold', 'Down', 'SKIPPER', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'A', 2.56012e-5, 'Cd', 1, 'MaxCv', 0.0, 'Cv', 0.0, 'State', 1.0, 'Tau', 0.0, 'IsOneWay', true);

    % Link 3: Valve_FU (Throttle valve: TK_FU -> FU_Regen_Inlet)
    Links.Valve_FU = struct('ID', 3, 'Name', 'Valve_FU', 'Type', 'Throttle', ...
        'Up', 'TK_FU', 'Down', 'FU_Regen_Inlet', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'MaxCv', 2.5, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05, 'IsOneWay', true);

    % Link 4: Regen_Line (Regen pressure drop: FU_Regen_Inlet -> FU_Manifold, tuned for 155 psi dP at target flow)
    Links.Regen_Line = struct('ID', 4, 'Name', 'Regen_Line', 'Type', 'Orifice', ...
        'Up', 'FU_Regen_Inlet', 'Down', 'FU_Manifold', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'MaxCv', 1.053, 'Cv', 1.053, 'State', 1.0, 'Tau', 0.0, 'IsOneWay', false);

    % Link 5: Inj_FU (Injector orifice: FU_Manifold -> SKIPPER)
    Links.Inj_FU = struct('ID', 5, 'Name', 'Inj_FU', 'Type', 'Orifice', ...
        'Up', 'FU_Manifold', 'Down', 'SKIPPER', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'A', 3.32903e-5, 'Cd', 1, 'MaxCv', 0.0, 'Cv', 0.0, 'State', 1.0, 'Tau', 0.0, 'IsOneWay', true);

    % Link 6: Nozzle (Sonic throat orifice: SKIPPER -> Atmosphere)
    Links.Nozzle = struct('ID', 6, 'Name', 'Nozzle', 'Type', 'Nozzle', ...
        'Up', 'SKIPPER', 'Down', 'Atmosphere', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'A', 7.14e-4, 'Cd', 0.98, 'MaxCv', 0.0, 'Cv', 0.0, 'State', 1.0, 'Tau', 0.0, 'IsOneWay', true);

    %% 3. Assemble Homogeneous Arrays & Connection Maps
    nodeNames = fieldnames(State.Nodes);
    N_N = numel(nodeNames);
    firstNode = standardizeNode(State.Nodes.(nodeNames{1}));
    NodeArray = repmat(firstNode, N_N, 1);
    for i = 1:N_N
        NodeArray(i) = standardizeNode(State.Nodes.(nodeNames{i}));
        NodeArray(i).ID = i;
    end

    linkNames = fieldnames(Links);
    N_L = numel(linkNames);
    firstLink = standardizeLink(Links.(linkNames{1}));
    LinkArray = repmat(firstLink, N_L, 1);
    LinkMap   = zeros(N_L, 2);
    PortMap   = zeros(N_L, 2);

    for j = 1:N_L
        fn = linkNames{j};
        LinkArray(j) = standardizeLink(Links.(fn));
        LinkArray(j).ID = j;

        u_id = find(strcmp(nodeNames, LinkArray(j).Up), 1);
        d_id = find(strcmp(nodeNames, LinkArray(j).Down), 1);
        if isempty(u_id), error('SKIPPER_StandA: Unknown Up node "%s"', LinkArray(j).Up); end
        if isempty(d_id), error('SKIPPER_StandA: Unknown Down node "%s"', LinkArray(j).Down); end

        LinkMap(j, :) = [u_id, d_id];
        PortMap(j, :) = [parsePortCode(LinkArray(j).UpPort), parsePortCode(LinkArray(j).DownPort)];
    end

    %% 4. Pack System & Initial State
    System.Nodes     = State.Nodes;
    System.Links     = Links;
    System.NodeArray = NodeArray;
    System.LinkArray = LinkArray;
    System.LinkMap   = LinkMap;
    System.PortMap   = PortMap;
    System.Env       = struct('T_amb', 293.15, 'WindVel', 1.5, 'P_atm', P_atm);

    System.Link.State = struct();
    for k = 1:N_L
        fn = linkNames{k};
        System.Link.State.(fn) = Links.(fn).State;
    end

    State.Time = 0.0;
    State.LinkStates = System.Link.State;
    State.NodeArray  = NodeArray;
    State.LinkArray  = LinkArray;
    InitialState = State;
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

%% --- Helper: Standardize Node to Homogeneous Schema ---
function n_std = standardizeNode(n)
    defZone = struct('m', 0.0, 'U', 0.0, 'V', 0.0, 'P', 0.0, 'T', 293.15, ...
                     'rho', 0.0, 'u', 0.0, 'h', 0.0, 'gamma', 1.4);
    defNode = struct(...
        'ID', 0, 'Name', '', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', 1e-3, 'D', 0.0, 'm', 0.0, 'U', 0.0, 'u', 0.0, 'rho', 0.0, ...
        'P', 101325.0, 'T', 293.15, 'h', 0.0, 'gamma', 1.4, ...
        'cstar', 0.0, 'isLit', false, 'Fixed', false, ...
        'SparkActive', false, 'TorchLit', false, ...
        'm_liq', 0.0, 'm_gas', 0.0, ...
        'Ullage', defZone, 'Liquid', defZone);
    n_std = defNode;
    f = fieldnames(n);
    for k = 1:numel(f)
        n_std.(f{k}) = n.(f{k});
    end
end

%% --- Helper: Standardize Link to Homogeneous Schema ---
function l_std = standardizeLink(l)
    defLink = struct(...
        'ID', 0, 'Name', '', 'Type', 'Valve', ...
        'Up', '', 'Down', '', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'MaxCv', 0.0, 'Cv', 0.0, 'State', 1.0, 'Tau', 0.05, ...
        'A', 0.0, 'Cd', 1.0, 'Zeta', 0.0, ...
        'P_set', 0.0, 'Droop', 0.0, 'SPE', 0.0, 'P_crack', 0.0, ...
        'IsOneWay', false, 'G', 0.0, 'Flow', 0.0, 'IsChoked', false, 'FlowDir', 1.0);
    l_std = defLink;
    f = fieldnames(l);
    for k = 1:numel(f)
        l_std.(f{k}) = l.(f{k});
    end
end
