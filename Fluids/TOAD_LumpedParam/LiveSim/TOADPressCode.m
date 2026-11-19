function [System, InitialState] = TOADPressCode()
%% Press Code System Def
% Topology:
%   Nodes:
%     1. TK_N2          - 36 L vehicle COPV (4,500 psi initial)
%     2. Press_Manifold - Regulated N2 distribution manifold (~1 L)
%     3. TK_O2_01       - Two-zone LOX tank (14.7 psi, 90 K liquid, 10% ullage)
%     4. TK_FU_01       - Two-zone Fuel tank (14.7 psi, 293.15 K liquid, 10% ullage)
%     5. Atmosphere     - Ambient boundary sink (14.7 psi)
%   Links:
%     1. REG_873D       - Dome-loaded pressure regulator (TK_N2 -> Press_Manifold)
%     2. Press_OX       - Pressurization check valve (Press_Manifold -> TK_O2_01 Ullage)
%     3. Press_FU       - Pressurization check valve (Press_Manifold -> TK_FU_01 Ullage)
%     4. OX_Drain       - Propellant liquid outflow link (TK_O2_01 Liquid -> Atmosphere)
%     5. FU_Drain       - Propellant liquid outflow link (TK_FU_01 Liquid -> Atmosphere)

    psi2Pa = 6894.757;
    P_atm  = 101325.0;

    %% Node Geometric & Baseline Definitions
    TankVolOX = 0.03075;   % m^3 (~30.75 L)
    TankVolFU = 0.03746;   % m^3 (~37.46 L)
    TankDia   = 11.5 * 0.0254;
    TargetUllage = 0.10;

    % Node 1: TK_N2 (Vehicle COPV)
    V_copv = 0.036;        % m^3 (36 L)
    P_copv0 = 4500.0 * psi2Pa;
    T_copv0 = 297.15;
    props_copv = FluidProperties('Nitrogen', 'From_P_T', P_copv0, T_copv0);
    m_copv0 = props_copv.rho * V_copv;
    State.Nodes.TK_N2 = struct(...
        'ID', 1, 'Name', 'TK_N2', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_copv, 'm', m_copv0, 'U', props_copv.u * m_copv0, ...
        'u', props_copv.u, 'rho', props_copv.rho, ...
        'P', P_copv0, 'T', T_copv0, 'h', props_copv.h, ...
        'gamma', props_copv.gamma, 'Fixed', false);

    % Pad unpressurized atmospheric pressure
    P_pad = P_atm; % 101325.0 Pa (14.696 psia)

    % Node 2: Press_Manifold (Regulated N2 Distribution Manifold)
    % (not physically representative but needed for sim stability, so end
    % consumption value is overestimated slightly as it includes this
    % volume. 2lbm maybe??)
    V_pman = 3.0e-3;       % 3 L manifold volume
    props_pman = FluidProperties('Nitrogen', 'From_P_T', P_pad, 293.15);
    m_pman0 = props_pman.rho * V_pman;
    State.Nodes.Press_Manifold = struct(...
        'ID', 2, 'Name', 'Press_Manifold', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_pman, 'm', m_pman0, 'U', props_pman.u * m_pman0, ...
        'u', props_pman.u, 'rho', props_pman.rho, ...
        'P', P_pad, 'T', 293.15, 'h', props_pman.h, ...
        'gamma', props_pman.gamma, 'Fixed', false);

    % Node 3: TK_O2_01 (Two-Zone Stratified LOX Tank)
    V_ull_ox = TankVolOX * TargetUllage;
    V_liq_ox = TankVolOX * (1.0 - TargetUllage);
    T_lox0   = 90.17; % K (LOX saturation at pad)
    props_ox_ull = FluidProperties('Nitrogen', 'From_P_T', P_pad, T_lox0);
    m_ull_ox0    = props_ox_ull.rho * V_ull_ox;
    U_ull_ox0    = props_ox_ull.u * m_ull_ox0;

    props_ox_liq = FluidProperties('Oxygen', 'From_P_T', P_pad, T_lox0);
    m_liq_ox0    = props_ox_liq.rho * V_liq_ox;
    U_liq_ox0    = props_ox_liq.u * m_liq_ox0;

    State.Nodes.TK_O2_01 = struct(...
        'ID', 3, 'Name', 'TK_O2_01', 'Type', 'TwoZoneTank', 'Fluid', 'Oxygen', ...
        'V', TankVolOX, 'D', TankDia, 'P', P_pad, 'Fixed', false, ...
        'Ullage', struct('m', m_ull_ox0, 'U', U_ull_ox0, 'V', V_ull_ox, ...
                         'P', P_pad, 'T', T_lox0, 'rho', props_ox_ull.rho, ...
                         'u', props_ox_ull.u, 'h', props_ox_ull.h, 'gamma', 1.4), ...
        'Liquid', struct('m', m_liq_ox0, 'U', U_liq_ox0, 'V', V_liq_ox, ...
                         'P', P_pad, 'T', T_lox0, 'rho', props_ox_liq.rho, ...
                         'u', props_ox_liq.u, 'h', props_ox_liq.h, 'gamma', 1.1));

    % Node 4: TK_FU_01 (Two-Zone Stratified Fuel Tank)
    V_ull_fu = TankVolFU * TargetUllage;
    V_liq_fu = TankVolFU * (1.0 - TargetUllage);
    T_fu0    = 293.15; % K
    props_fu_ull = FluidProperties('Nitrogen', 'From_P_T', P_pad, T_fu0);
    m_ull_fu0    = props_fu_ull.rho * V_ull_fu;
    U_ull_fu0    = props_fu_ull.u * m_ull_fu0;

    props_fu_liq = FluidProperties('IPA', 'From_P_T', P_pad, T_fu0);
    m_liq_fu0    = props_fu_liq.rho * V_liq_fu;
    U_liq_fu0    = props_fu_liq.u * m_liq_fu0;

    State.Nodes.TK_FU_01 = struct(...
        'ID', 4, 'Name', 'TK_FU_01', 'Type', 'TwoZoneTank', 'Fluid', 'IPA', ...
        'V', TankVolFU, 'D', TankDia, 'P', P_pad, 'Fixed', false, ...
        'Ullage', struct('m', m_ull_fu0, 'U', U_ull_fu0, 'V', V_ull_fu, ...
                         'P', P_pad, 'T', T_fu0, 'rho', props_fu_ull.rho, ...
                         'u', props_fu_ull.u, 'h', props_fu_ull.h, 'gamma', 1.4), ...
        'Liquid', struct('m', m_liq_fu0, 'U', U_liq_fu0, 'V', V_liq_fu, ...
                         'P', P_pad, 'T', T_fu0, 'rho', props_fu_liq.rho, ...
                         'u', props_fu_liq.u, 'h', props_fu_liq.h, 'gamma', 1.15));

    % Node 5: Atmosphere
    State.Nodes.Atmosphere = struct(...
        'ID', 5, 'Name', 'Atmosphere', 'Type', 'Boundary', 'Fluid', 'Nitrogen', ...
        'V', 1.0, 'm', 1.18, 'U', 2e5, 'u', 1.7e5, 'rho', 1.18, ...
        'P', P_atm, 'T', 293.15, 'h', 2.8e5, 'gamma', 1.4, 'Fixed', true);


    %% Link Definitions
    Links = struct();

    % 1. REG_873D: COPV -> Press_Manifold
    Links.REG_873D = struct('ID', 1, 'Name', 'REG_873D', 'Type', 'Regulator', ...
        'Up', 'TK_N2', 'Down', 'Press_Manifold', 'UpPort', 'Bulk', 'DownPort', 'Bulk', ...
        'MaxCv', 0.70, 'Cv', 0.70, 'State', 1.0, ...
        'P_set', 563.5 * psi2Pa, 'Droop', 30 * psi2Pa, 'SPE', 0.003, 'IsOneWay', true);

    % Press_OX: Press_Manifold -> TK_O2_01 Ullage Check Valve (Cv = 0.70)
    Links.Press_OX = struct('ID', 2, 'Name', 'Press_OX', 'Type', 'Check', ...
        'Up', 'Press_Manifold', 'Down', 'TK_O2_01', 'UpPort', 'Bulk', 'DownPort', 'Ullage', ...
        'MaxCv', 0.70, 'Cv', 0.70, 'State', 1.0, 'P_crack', 0.0, 'IsOneWay', true);

    % Press_FU: Press_Manifold -> TK_FU_01 Ullage Check Valve (Cv = 0.70)
    Links.Press_FU = struct('ID', 3, 'Name', 'Press_FU', 'Type', 'Check', ...
        'Up', 'Press_Manifold', 'Down', 'TK_FU_01', 'UpPort', 'Bulk', 'DownPort', 'Ullage', ...
        'MaxCv', 0.70, 'Cv', 0.70, 'State', 1.0, 'P_crack', 0.0, 'IsOneWay', true);

    % OX_Drain: TK_O2_01 Liquid Outflow Link -> Atmosphere
    Links.OX_Drain = struct('ID', 4, 'Name', 'OX_Drain', 'Type', 'Throttle', ...
        'Up', 'TK_O2_01', 'Down', 'Atmosphere', 'UpPort', 'Liquid', 'DownPort', 'Bulk', ...
        'MaxCv', 10.0, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.0, 'IsOneWay', true);

    % FU_Drain: TK_FU_01 Liquid Outflow Link -> Atmosphere
    Links.FU_Drain = struct('ID', 5, 'Name', 'FU_Drain', 'Type', 'Throttle', ...
        'Up', 'TK_FU_01', 'Down', 'Atmosphere', 'UpPort', 'Liquid', 'DownPort', 'Bulk', ...
        'MaxCv', 10.0, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.0, 'IsOneWay', true);

    %% Thermal State Initialization
    State.Thermal.COPV     = struct('T_layers', [293.15, 293.15, 293.15]);
    State.Thermal.TK_O2_01 = struct('T_wall', 90.17);
    State.Thermal.TK_FU_01 = struct('T_wall', 293.15);

    %% Assemble Numerical Arrays and Maps for StepSimulation
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
        if isempty(u_id), error('TOADPressCode: Unknown Up node "%s"', LinkArray(j).Up); end
        if isempty(d_id), error('TOADPressCode: Unknown Down node "%s"', LinkArray(j).Down); end

        LinkMap(j, :) = [u_id, d_id];
        PortMap(j, :) = [parsePortCode(LinkArray(j).UpPort), parsePortCode(LinkArray(j).DownPort)];
    end

    %% Pack System Definition & Initial Link States
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

%% --- Helper: Parse Port Code ---
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
