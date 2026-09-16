function Diagnose_Pressure_Equalization()
    fprintf('=====================================================\n');
    fprintf('  DIAGNOSING PRESSURE EQUALIZATION IN STEPSIMULATION_LIVE\n');
    fprintf('=====================================================\n');

    baseDir = fileparts(mfilename('fullpath'));
    addpath(fullfile(baseDir, '..', 'Core'));
    addpath(fullfile(baseDir, '..', 'Config'));
    addpath(fullfile(baseDir, '..', 'Thermal'));
    addpath(fullfile(baseDir, '..', 'Fluid Properties'));
    addpath(baseDir);

    psi2Pa = 6894.757;
    dt = 0.010;

    %% ---------------------------------------------------------------
    %% TEST 1: BULK N2 FILL (BV_N2_FILL)
    %% ---------------------------------------------------------------
    fprintf('\n--- TEST 1: BV_N2_FILL with COPV at 2500 psi and Bulk at 6000 psi ---\n');
    [System, State] = BuildTOADSystem_Merged();
    
    % Set COPV to 2500 psi partially depleted
    P_copv_init = 2500 * psi2Pa;
    props_copv = FluidProperties('Nitrogen', 'From_P_T', P_copv_init, 297.15);
    State.Nodes.TK_N2.P = P_copv_init;
    State.Nodes.TK_N2.rho = props_copv.rho;
    State.Nodes.TK_N2.m = props_copv.rho * State.Nodes.TK_N2.V;
    State.Nodes.TK_N2.U = props_copv.u * State.Nodes.TK_N2.m;
    State.Nodes.TK_N2.u = props_copv.u;
    State.Nodes.TK_N2.h = props_copv.h;

    fprintf('Initial Bulk P: %.1f psi | Initial COPV P: %.1f psi\n', ...
        State.Nodes.TK_N2_BULK.P / psi2Pa, State.Nodes.TK_N2.P / psi2Pa);

    LinkStates = struct();
    LinkStates.BV_N2_FILL = 1.0;

    for s = 1:500 % 5.0 seconds
        [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates);
        if mod(s, 50) == 0 || s == 1
            fprintf('  t=%.2fs | mdot_fill=%.4f kg/s | Bulk P=%.1f psi | COPV P=%.1f psi\n', ...
                State.Time, FlowRates.BV_N2_FILL, State.Nodes.TK_N2_BULK.P / psi2Pa, State.Nodes.TK_N2.P / psi2Pa);
        end
    end

    %% ---------------------------------------------------------------
    %% TEST 2: PROPELLANT TANK PRESSURIZATION (BV_N2_02)
    %% ---------------------------------------------------------------
    fprintf('\n--- TEST 2: BV_N2_02 Pressurization (COPV -> Reg -> Tanks) ---\n');
    [System, State] = BuildTOADSystem_Merged();
    fprintf('Initial: COPV: %.1f psi | PurgeMan: %.1f psi | PressMan: %.1f psi | OX Ull: %.1f psi | FU Ull: %.1f psi\n', ...
        State.Nodes.TK_N2.P / psi2Pa, State.Nodes.Purge_Manifold.P / psi2Pa, ...
        State.Nodes.Press_Manifold.P / psi2Pa, State.Nodes.TK_O2_01.P / psi2Pa, State.Nodes.TK_FU_01.P / psi2Pa);

    LinkStates = struct();
    LinkStates.BV_N2_02 = 1.0;

    for s = 1:500 % 5.0 seconds
        [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates);
        if mod(s, 50) == 0 || s <= 10
            fprintf('  t=%.2fs | mdot_ox=%.4f | mdot_fu=%.4f | Reg=%.1f | PressMan=%.1f | OX=%.1f | FU=%.1f psi\n', ...
                State.Time, FlowRates.Press_OX, FlowRates.Press_FU, ...
                System.Links.REG_873D.P_set / psi2Pa, State.Nodes.Press_Manifold.P / psi2Pa, ...
                State.Nodes.TK_O2_01.P / psi2Pa, State.Nodes.TK_FU_01.P / psi2Pa);
        end
    end

    %% ---------------------------------------------------------------
    %% TEST 3: PROPELLANT TANK VENTING (BV_O2_01 & BV_FU_01)
    %% ---------------------------------------------------------------
    fprintf('\n--- TEST 3: Tank Venting with BV_N2_02 Closed ---\n');
    LinkStates = struct();
    LinkStates.BV_N2_02 = 0.0;
    LinkStates.BV_O2_01 = 1.0;
    LinkStates.BV_FU_01 = 1.0;

    for s = 1:500 % 5.0 seconds
        [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates);
        if mod(s, 50) == 0 || s <= 10
            fprintf('  t=%.2fs | vent_ox=%.4f | vent_fu=%.4f | OX=%.1f | FU=%.1f psi\n', ...
                State.Time, FlowRates.BV_O2_01, FlowRates.BV_FU_01, ...
                State.Nodes.TK_O2_01.P / psi2Pa, State.Nodes.TK_FU_01.P / psi2Pa);
        end
    end

    %% ---------------------------------------------------------------
    %% TEST 4: RCS FIRING (SV_N2_01)
    %% ---------------------------------------------------------------
    fprintf('\n--- TEST 4: RCS Thruster Firing (SV_N2_01) ---\n');
    [System, State] = BuildTOADSystem_Merged();
    LinkStates = struct();
    LinkStates.SV_N2_01 = 1.0;

    for s = 1:200 % 2.0 seconds
        [State, FlowRates, System] = StepSimulation_Live(State, dt, System, LinkStates);
        if mod(s, 20) == 0 || s <= 5
            fprintf('  t=%.2fs | mdot_rcs=%.4f | PurgeMan=%.1f | COPV=%.1f psi\n', ...
                State.Time, FlowRates.SV_N2_01, ...
                State.Nodes.Purge_Manifold.P / psi2Pa, State.Nodes.TK_N2.P / psi2Pa);
        end
    end
end
