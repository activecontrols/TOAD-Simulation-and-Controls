function [System, InitialState] = BuildTOADSystem_Merged()
% BUILDTOADSYSTEM_MERGED Assembles the unified TOAD fluid system topology,
% fully merged between the Simulation Lumped Parameter model and Ground Control UI.
%
% Initial Conditions: Always Flight-Ready
%   - TK-N2 (COPV): Pre-charged to 4,500 psi
%   - TK-N2-BULK: Ground N2 supply at 6,000 psi
%   - TK-O2-01 & TK-FU-01: Propellant tanks filled (10% ullage) at 1 atm
%   - RCS Solenoids (SV-N2-01..04) and Vent/Drain Valves included

    psi2Pa = 6894.757;
    P_atm  = 14.7 * psi2Pa;

    %% 1. Node Geometric & Baseline Definitions
    TankVolOX = 0.03075;
    TankVolFU = 0.03746;
    StartP_Pa = P_atm;
    TargetUllage = 0.10; % 10% ullage fraction

    % --- Node 0: TK-N2-BULK (Ground N2 Supply, 6000 psi, 0.3 m^3) ---
    V_bulk = 0.30;
    P_bulk0 = 6000 * psi2Pa;
    T_bulk0 = 297.15;
    props_bulk = FluidProperties('Nitrogen', 'From_P_T', P_bulk0, T_bulk0);
    m_bulk0 = props_bulk.rho * V_bulk;
    U_bulk0 = props_bulk.u * m_bulk0;

    State.Nodes.TK_N2_BULK = struct(...
        'Name', 'TK-N2-BULK', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_bulk, 'm', m_bulk0, 'U', U_bulk0, ...
        'u', props_bulk.u, 'rho', props_bulk.rho, ...
        'P', P_bulk0, 'T', T_bulk0, 'h', props_bulk.h, ...
        'gamma', props_bulk.gamma, 'Fixed', false);

    % --- Node 1: TK-N2 (Vehicle COPV) - Flight Ready at 4,500 psi ---
    V_copv = 0.036;
    P_copv0 = 4500 * psi2Pa;
    T_copv0 = 297.15;
    props_copv = FluidProperties('Nitrogen', 'From_P_T', P_copv0, T_copv0);
    m_copv0 = props_copv.rho * V_copv;
    U_copv0 = props_copv.u * m_copv0;

    State.Nodes.TK_N2 = struct(...
        'Name', 'TK-N2', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_copv, 'm', m_copv0, 'U', U_copv0, ...
        'u', props_copv.u, 'rho', props_copv.rho, ...
        'P', P_copv0, 'T', T_copv0, 'h', props_copv.h, ...
        'gamma', props_copv.gamma, 'Fixed', false);

    % --- Node 16: Purge Manifold (Main Vehicle N2 Distribution Manifold) ---
    V_pm = 3.0e-3; % 3.0 L
    P_pm0 = P_atm;
    T_pm0 = 293.15;
    props_pm = FluidProperties('Nitrogen', 'From_P_T', P_pm0, T_pm0);
    m_pm0 = props_pm.rho * V_pm;
    U_pm0 = props_pm.u * m_pm0;

    State.Nodes.Purge_Manifold = struct(...
        'Name', 'Purge Manifold', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_pm, 'm', m_pm0, 'U', U_pm0, ...
        'u', props_pm.u, 'rho', props_pm.rho, ...
        'P', P_pm0, 'T', T_pm0, 'h', props_pm.h, ...
        'gamma', props_pm.gamma, 'Fixed', false);

    % --- Node 2: Press Manifold (Maintained for backward compatibility) ---
    V_pman = 1.0e-3;
    P_pman0 = P_atm;
    T_pman0 = 293.15;
    props_pman = FluidProperties('Nitrogen', 'From_P_T', P_pman0, T_pman0);
    m_pman0 = props_pman.rho * V_pman;
    U_pman0 = props_pman.u * m_pman0;

    State.Nodes.Press_Manifold = struct(...
        'Name', 'Press Manifold', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_pman, 'm', m_pman0, 'U', U_pman0, ...
        'u', props_pman.u, 'rho', props_pman.rho, ...
        'P', P_pman0, 'T', T_pman0, 'h', props_pman.h, ...
        'gamma', props_pman.gamma, 'Fixed', false);

    % --- Node 18: Inter N2 Reg (Between REG-873D and BV-N2-02) ---
    V_inr = 5.0e-4; % 0.5 L
    P_inr0 = 550 * psi2Pa; % Pre-charged to regulated pressure in flight-ready state
    T_inr0 = 293.15;
    props_inr = FluidProperties('Nitrogen', 'From_P_T', P_inr0, T_inr0);
    m_inr0 = props_inr.rho * V_inr;
    U_inr0 = props_inr.u * m_inr0;

    State.Nodes.Inter_N2_Reg = struct(...
        'Name', 'Inter N2 Reg', 'Type', 'Gas', 'Fluid', 'Nitrogen', ...
        'V', V_inr, 'm', m_inr0, 'U', U_inr0, ...
        'u', props_inr.u, 'rho', props_inr.rho, ...
        'P', P_inr0, 'T', T_inr0, 'h', props_inr.h, ...
        'gamma', props_inr.gamma, 'Fixed', false);

    % --- Node 3: TK-O2-01 (LOX Tank, Two-Zone Stratified) ---
    V_ull_ox = TankVolOX * TargetUllage;
    V_liq_ox = TankVolOX * (1 - TargetUllage);
    T_lox0   = 90.0; % K

    props_ox_ull = FluidProperties('Nitrogen', 'From_P_T', StartP_Pa, 293.15);
    m_ull_ox0    = props_ox_ull.rho * V_ull_ox;
    U_ull_ox0    = props_ox_ull.u * m_ull_ox0;

    props_ox_liq = FluidProperties('Oxygen', 'From_P_T', StartP_Pa, T_lox0);
    m_liq_ox0    = props_ox_liq.rho * V_liq_ox;
    U_liq_ox0    = props_ox_liq.u * m_liq_ox0;

    State.Nodes.TK_O2_01 = struct(...
        'Name', 'TK-O2-01', 'Type', 'TwoZoneTank', 'Fluid', 'Oxygen', ...
        'V', TankVolOX, 'P', StartP_Pa, 'Fixed', false, ...
        'Ullage', struct('m', m_ull_ox0, 'U', U_ull_ox0, 'V', V_ull_ox, ...
                         'P', StartP_Pa, 'T', 293.15, 'rho', props_ox_ull.rho, ...
                         'u', props_ox_ull.u, 'h', props_ox_ull.h, 'gamma', 1.4), ...
        'Liquid', struct('m', m_liq_ox0, 'U', U_liq_ox0, 'V', V_liq_ox, ...
                         'P', StartP_Pa, 'T', T_lox0, 'rho', props_ox_liq.rho, ...
                         'u', props_ox_liq.u, 'h', props_ox_liq.h, 'gamma', 1.1));

    % --- Node 4: TK-FU-01 (Fuel Tank, Two-Zone Stratified) ---
    V_ull_fu = TankVolFU * TargetUllage;
    V_liq_fu = TankVolFU * (1 - TargetUllage);
    T_fu0    = 293.15; % K

    props_fu_ull = FluidProperties('Nitrogen', 'From_P_T', StartP_Pa, T_fu0);
    m_ull_fu0    = props_fu_ull.rho * V_ull_fu;
    U_ull_fu0    = props_fu_ull.u * m_ull_fu0;

    props_fu_liq = FluidProperties('IPA', 'From_P_T', StartP_Pa, T_fu0);
    m_liq_fu0    = props_fu_liq.rho * V_liq_fu;
    U_liq_fu0    = props_fu_liq.u * m_liq_fu0;

    State.Nodes.TK_FU_01 = struct(...
        'Name', 'TK-FU-01', 'Type', 'TwoZoneTank', 'Fluid', 'IPA', ...
        'V', TankVolFU, 'P', StartP_Pa, 'Fixed', false, ...
        'Ullage', struct('m', m_ull_fu0, 'U', U_ull_fu0, 'V', V_ull_fu, ...
                         'P', StartP_Pa, 'T', T_fu0, 'rho', props_fu_ull.rho, ...
                         'u', props_fu_ull.u, 'h', props_fu_ull.h, 'gamma', 1.4), ...
        'Liquid', struct('m', m_liq_fu0, 'U', U_liq_fu0, 'V', V_liq_fu, ...
                         'P', StartP_Pa, 'T', T_fu0, 'rho', props_fu_liq.rho, ...
                         'u', props_fu_liq.u, 'h', props_fu_liq.h, 'gamma', 1.15));

    % --- Node 5: Pre Main OX (Primed Liquid LOX) ---
    V_pmo = 1e-5;
    m_pmo0 = props_ox_liq.rho * V_pmo;
    State.Nodes.Pre_Main_OX = struct(...
        'Name', 'Pre Main OX', 'Type', 'LiquidLine', 'Fluid', 'Oxygen', ...
        'V', V_pmo, 'm', m_pmo0, 'U', props_ox_liq.u * m_pmo0, ...
        'u', props_ox_liq.u, 'rho', props_ox_liq.rho, ...
        'P', StartP_Pa, 'T', T_lox0, 'h', props_ox_liq.h, ...
        'gamma', 1.1, 'Fixed', false);

    % --- Node 6: Pre Main FU (Primed Liquid IPA) ---
    V_pmf = 1e-5;
    m_pmf0 = props_fu_liq.rho * V_pmf;
    State.Nodes.Pre_Main_FU = struct(...
        'Name', 'Pre Main FU', 'Type', 'LiquidLine', 'Fluid', 'IPA', ...
        'V', V_pmf, 'm', m_pmf0, 'U', props_fu_liq.u * m_pmf0, ...
        'u', props_fu_liq.u, 'rho', props_fu_liq.rho, ...
        'P', StartP_Pa, 'T', T_fu0, 'h', props_fu_liq.h, ...
        'gamma', 1.15, 'Fixed', false);

    % --- Line Nodes: Inter OX, Inter FU, Post Throttle OX/FU, Manifolds ---
    lineDefs = {
        'Inter_OX',          'Inter OX',          1.0e-3, 'Oxygen';  % 1.0 L
        'Inter_FU',          'Inter FU',          1.0e-3, 'IPA';     % 1.0 L
        'Post_Throttle_OX',  'Post Throttle OX',  1.0e-3, 'Oxygen';  % 1.0 L
        'Post_Throttle_FU',  'Post Throttle FU',  1.0e-3, 'IPA';     % 1.0 L
        'OX_Manifold',       'OX Manifold',       2.0e-3, 'Oxygen';  % 2.0 L
        'FU_Manifold',       'FU Manifold',       2.0e-3, 'IPA';     % 2.0 L
    };

    props_gas_init = FluidProperties('Nitrogen', 'From_P_T', P_atm, 293.15);
    for k = 1:size(lineDefs, 1)
        fName = lineDefs{k, 1};
        uName = lineDefs{k, 2};
        vol   = lineDefs{k, 3};
        propF = lineDefs{k, 4};

        m_init = props_gas_init.rho * vol;
        State.Nodes.(fName) = struct(...
            'Name', uName, 'Type', 'Line', 'Fluid', propF, ...
            'V', vol, 'm_prop', 0.0, 'm_n2', m_init, 'm', m_init, ...
            'U', props_gas_init.u * m_init, 'u', props_gas_init.u, ...
            'rho', props_gas_init.rho, 'P', P_atm, 'T', 293.15, ...
            'h', props_gas_init.h, 'gamma', 1.4, 'Fixed', false);
    end

    % --- Node 13: SKIPPER (Main Combustor) ---
    V_skip = 0.001302;
    m_skip0 = props_gas_init.rho * V_skip;
    State.Nodes.SKIPPER = struct(...
        'Name', 'SKIPPER', 'Type', 'Combustor', 'Fluid', 'CombustionGas', ...
        'V', V_skip, 'm', m_skip0, 'U', props_gas_init.u * m_skip0, ...
        'u', props_gas_init.u, 'rho', props_gas_init.rho, ...
        'P', P_atm, 'T', 293.15, 'h', props_gas_init.h, ...
        'gamma', 1.21, 'cstar', 1600.0, 'isLit', false, 'Fixed', false);

    % --- Node 17: DART Chamber (Torch Igniter) ---
    V_dart = 5.0e-4; % 0.5 L
    m_dart0 = props_gas_init.rho * V_dart;
    State.Nodes.DART_Chamber = struct(...
        'Name', 'DART Chamber', 'Type', 'Combustor', 'Fluid', 'CombustionGas', ...
        'V', V_dart, 'm', m_dart0, 'U', props_gas_init.u * m_dart0, ...
        'u', props_gas_init.u, 'rho', props_gas_init.rho, ...
        'P', P_atm, 'T', 293.15, 'h', props_gas_init.h, ...
        'gamma', 1.21, 'cstar', 1600.0, 'isLit', false, 'Fixed', false);

    % --- Node 14: Atmosphere (Fixed Boundary) ---
    State.Nodes.Atmosphere = struct(...
        'Name', 'Atmosphere', 'Type', 'Boundary', 'Fluid', 'Nitrogen', ...
        'V', 1.0, 'm', 1.18, 'U', 2e5, 'u', 1.7e5, 'rho', 1.18, ...
        'P', P_atm, 'T', 293.15, 'h', 2.8e5, 'gamma', 1.4, 'Fixed', true);

    % --- Node 15: Virtual (Spark Signal Source) ---
    State.Nodes.Virtual = struct(...
        'Name', 'Virtual', 'Type', 'Boundary', 'Fluid', 'Nitrogen', ...
        'V', 1.0, 'm', 1.0, 'U', 0.0, 'u', 0.0, 'rho', 1.0, ...
        'P', 0.0, 'T', 0.0, 'h', 0.0, 'gamma', 1.4, 'Fixed', true);

    %% 2. Link Definitions (Merged System)
    Links = struct();

    % 1. REG-873D (Regulator: TK_N2 -> Inter_N2_Reg)
    Links.REG_873D = struct('ID', 1, 'Name', 'REG-873D', 'Type', 'Regulator', ...
        'Up', 'TK_N2', 'Down', 'Inter_N2_Reg', 'MaxCv', 0.70, 'Cv', 0.70, 'State', 1.0, ...
        'P_set', 550 * psi2Pa, 'Droop', 30 * psi2Pa, 'SPE', 0.003);

    % 20. BV-N2-02 (Main Vehicle N2 Isolation Valve: Inter_N2_Reg -> Purge_Manifold)
    Links.BV_N2_02 = struct('ID', 20, 'Name', 'BV-N2-02', 'Type', 'Throttle', ...
        'Up', 'Inter_N2_Reg', 'Down', 'Purge_Manifold', 'MaxCv', 2.0, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 2. Press OX Check (Purge_Manifold -> TK_O2_01 Ullage)
    Links.Press_OX = struct('ID', 2, 'Name', 'Press OX', 'Type', 'Check', ...
        'Up', 'Purge_Manifold', 'Down', 'TK_O2_01', 'MaxCv', 25.0, 'Cv', 25.0, 'State', 1.0, 'P_crack', 0.0);

    % 3. Press FU Check (Purge_Manifold -> TK_FU_01 Ullage)
    Links.Press_FU = struct('ID', 3, 'Name', 'Press FU', 'Type', 'Check', ...
        'Up', 'Purge_Manifold', 'Down', 'TK_FU_01', 'MaxCv', 25.0, 'Cv', 25.0, 'State', 1.0, 'P_crack', 0.0);

    % 4. OX Line 1
    Links.OX_Line_1 = struct('ID', 4, 'Name', 'OX Line 1', 'Type', 'Pipe', ...
        'Up', 'TK_O2_01', 'Down', 'Pre_Main_OX', 'A', 5e-4, 'Zeta', 20.0, 'Cv', 0, 'State', 1.0);

    % 5. FU Line 1
    Links.FU_Line_1 = struct('ID', 5, 'Name', 'FU Line 1', 'Type', 'Pipe', ...
        'Up', 'TK_FU_01', 'Down', 'Pre_Main_FU', 'A', 5e-4, 'Zeta', 20.0, 'Cv', 0, 'State', 1.0);

    % 6. BV-02-03 (Main OX)
    Links.BV_02_03 = struct('ID', 6, 'Name', 'BV-02-03', 'Type', 'Solenoid', ...
        'Up', 'Pre_Main_OX', 'Down', 'Inter_OX', 'MaxCv', 2.9, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 7. BV-FU-03 (Main FU)
    Links.BV_FU_03 = struct('ID', 7, 'Name', 'BV-FU-03', 'Type', 'Solenoid', ...
        'Up', 'Pre_Main_FU', 'Down', 'Inter_FU', 'MaxCv', 2.9, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 8. BV-02-04 (Throttle OX)
    Links.BV_02_04 = struct('ID', 8, 'Name', 'BV-02-04', 'Type', 'Throttle', ...
        'Up', 'Inter_OX', 'Down', 'Post_Throttle_OX', 'MaxCv', 1.18, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.07);

    % 9. BV-FU-04 (Throttle FU)
    Links.BV_FU_04 = struct('ID', 9, 'Name', 'BV-FU-04', 'Type', 'Throttle', ...
        'Up', 'Inter_FU', 'Down', 'Post_Throttle_FU', 'MaxCv', 1.18, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.07);

    % 10. OX Inj Line
    Links.OX_Inj_Line = struct('ID', 10, 'Name', 'OX Inj Line', 'Type', 'Pipe', ...
        'Up', 'Post_Throttle_OX', 'Down', 'OX_Manifold', 'A', 5e-4, 'Zeta', 20.0, 'Cv', 0, 'State', 1.0);

    % 11. FU Inj Line
    Links.FU_Inj_Line = struct('ID', 11, 'Name', 'FU Inj Line', 'Type', 'Pipe', ...
        'Up', 'Post_Throttle_FU', 'Down', 'FU_Manifold', 'A', 5e-4, 'Zeta', 20.0, 'Cv', 0, 'State', 1.0);

    % 12. Inj OX
    Links.Inj_OX = struct('ID', 12, 'Name', 'Inj OX', 'Type', 'Orifice', ...
        'Up', 'OX_Manifold', 'Down', 'SKIPPER', 'A', 3.25e-5, 'Cd', 0.45, 'Cv', 0, 'State', 1.0);

    % 13. Inj FU
    Links.Inj_FU = struct('ID', 13, 'Name', 'Inj FU', 'Type', 'Orifice', ...
        'Up', 'FU_Manifold', 'Down', 'SKIPPER', 'A', 3.65e-5, 'Cd', 0.77, 'Cv', 0, 'State', 1.0);

    % 14. Nozzle
    Links.Nozzle = struct('ID', 14, 'Name', 'Nozzle', 'Type', 'Orifice', ...
        'Up', 'SKIPPER', 'Down', 'Atmosphere', 'A', 0.00129717, 'Cd', 0.95, 'Cv', 0, 'State', 1.0);

    % 15. SV-N2-05 (OX Purge)
    Links.SV_N2_05 = struct('ID', 15, 'Name', 'SV-N2-05', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Post_Throttle_OX', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 16. SV-N2-06 (FU Purge)
    Links.SV_N2_06 = struct('ID', 16, 'Name', 'SV-N2-06', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Post_Throttle_FU', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 22. SV-N2-07 (DART Purge)
    Links.SV_N2_07 = struct('ID', 22, 'Name', 'SV-N2-07', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'DART_Chamber', 'MaxCv', 0.10, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.10);

    % 17. SV-DART-OX (UI: SV-O2-01)
    Links.SV_DART_OX = struct('ID', 17, 'Name', 'SV-DART-OX', 'Type', 'Solenoid', ...
        'Up', 'Inter_OX', 'Down', 'DART_Chamber', 'MaxCv', 0.014, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.04);

    % 18. SV-DART-FU (UI: SV-FU-01)
    Links.SV_DART_FU = struct('ID', 18, 'Name', 'SV-DART-FU', 'Type', 'Solenoid', ...
        'Up', 'Inter_FU', 'Down', 'DART_Chamber', 'MaxCv', 0.0195, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.04);

    % 19. Spark (Signal Link - Held Off in manual mode)
    Links.Spark = struct('ID', 19, 'Name', 'Spark', 'Type', 'Signal', ...
        'Up', 'Virtual', 'Down', 'Virtual', 'MaxCv', 1.0, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.001);

    % 21. DART Nozzle
    Links.DART_Nozzle = struct('ID', 21, 'Name', 'DART Nozzle', 'Type', 'Orifice', ...
        'Up', 'DART_Chamber', 'Down', 'SKIPPER', 'A', 2.62e-5, 'Cd', 1.0, 'Cv', 0, 'State', 1.0);

    % 23. BV-N2-FILL (Ground Bulk N2 Fill Valve)
    Links.BV_N2_FILL = struct('ID', 23, 'Name', 'BV-N2-FILL', 'Type', 'Ball', ...
        'Up', 'TK_N2_BULK', 'Down', 'TK_N2', 'MaxCv', 1.5, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    %% --- MERGED P&ID COMPONENTS ADDED FROM UI ---
    % 24..27. RCS Cold-Gas Thrusters (SV-N2-01 .. SV-N2-04)
    Links.SV_N2_01 = struct('ID', 24, 'Name', 'SV-N2-01', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Atmosphere', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.02);
    Links.SV_N2_02 = struct('ID', 25, 'Name', 'SV-N2-02', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Atmosphere', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.02);
    Links.SV_N2_03 = struct('ID', 26, 'Name', 'SV-N2-03', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Atmosphere', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.02);
    Links.SV_N2_04 = struct('ID', 27, 'Name', 'SV-N2-04', 'Type', 'Solenoid', ...
        'Up', 'Purge_Manifold', 'Down', 'Atmosphere', 'MaxCv', 0.20, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.02);

    % 28. BV-N2-01: COPV High-Pressure Manual Dump/Vent
    Links.BV_N2_01 = struct('ID', 28, 'Name', 'BV-N2-01', 'Type', 'Solenoid', ...
        'Up', 'TK_N2', 'Down', 'Atmosphere', 'MaxCv', 0.50, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 29. BV-O2-01: LOX Tank Ullage Vent
    Links.BV_O2_01 = struct('ID', 29, 'Name', 'BV-O2-01', 'Type', 'Solenoid', ...
        'Up', 'TK_O2_01', 'Down', 'Atmosphere', 'MaxCv', 1.00, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 30. BV-O2-02: LOX Tank Liquid Drain
    Links.BV_O2_02 = struct('ID', 30, 'Name', 'BV-O2-02', 'Type', 'Solenoid', ...
        'Up', 'TK_O2_01', 'Down', 'Atmosphere', 'MaxCv', 1.00, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 31. BV-FU-01: Fuel Tank Ullage Vent
    Links.BV_FU_01 = struct('ID', 31, 'Name', 'BV-FU-01', 'Type', 'Solenoid', ...
        'Up', 'TK_FU_01', 'Down', 'Atmosphere', 'MaxCv', 1.00, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    % 32. BV-FU-02: Fuel Tank Liquid Drain
    Links.BV_FU_02 = struct('ID', 32, 'Name', 'BV-FU-02', 'Type', 'Solenoid', ...
        'Up', 'TK_FU_01', 'Down', 'Atmosphere', 'MaxCv', 1.00, 'Cv', 0.0, 'State', 0.0, 'Tau', 0.05);

    %% 3. Thermal State Initialization
    State.Thermal.COPV     = struct('T_layers', [293.15, 293.15, 293.15]);
    State.Thermal.TK_O2_01 = struct('T_wall', 120.0);
    State.Thermal.TK_FU_01 = struct('T_wall', 293.15);
    State.Thermal.Regen    = struct('T_wall', 293.15, 'DeltaP', 0.0, 'Qdot', 0.0);

    %% 4. Pack System Definition & Initial Link States
    System.Links = Links;
    System.Env = struct('T_amb', 293.15, 'WindVel', 1.5, 'P_atm', P_atm);

    System.Link.State = struct();
    linkNames = fieldnames(Links);
    for k = 1:length(linkNames)
        fn = linkNames{k};
        System.Link.State.(fn) = Links.(fn).State;
    end

    State.Time = 0.0;
    State.LinkStates = System.Link.State;
    InitialState = State;
end
