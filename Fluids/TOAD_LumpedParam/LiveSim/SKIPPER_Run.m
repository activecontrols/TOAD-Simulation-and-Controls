%% SKIPPER_Run.m: Master Simulation Script for SKIPPER Test Stand A
% Executes full 15-second hot-fire test sequence:
%   - 0.0 to 5.0 s:  Cold Dwell (Valves closed, IGN = 0)
%   - 5.0 to 10.0 s: Hot-Fire Burn (Valves open, IGN = 1)
%   - 10.0 to 15.0 s: Cutoff Dwell (Valves closed, IGN = 0)
%
% Architecture & schemas modeled after PressCode_Run.m and TOADPressCode.m.
% Includes clean modular valve actuation & angle-to-Cv controller interface.
% Saves run outputs to sandbox/experiments/SKIPPER_Run_Out.mat and
% multi-panel diagnostics figure to sandbox/experiments/SKIPPER_Run_Results.png.

clear; clc; close all;
fprintf('SKIPPER Hotfire Sim');

%% Path Management
scriptDir = fileparts(mfilename('fullpath'));
candidates = {
    scriptDir, ...
    fullfile(scriptDir, '..', 'Core'), ...
    fullfile(scriptDir, '..', 'Thermal'), ...
    fullfile(scriptDir, '..', 'Helpers'), ...
    fullfile(scriptDir, '..', 'Config'), ...
    fullfile(scriptDir, '..', 'Fluid Properties'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'LiveSim'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Core'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Thermal'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Helpers'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Fluid Properties'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'LiveSim'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Core'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Fluid Properties'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Helpers')
};
for c = 1:numel(candidates)
    if exist(candidates{c}, 'dir'), addpath(candidates{c}); end
end

%% System Initialization
[System, State] = SKIPPER_StandA();
psi2Pa = 6894.757;
P_atm  = System.Env.P_atm;

% Simulation Timing Parameters
totalTime = 15.0;       % Total duration: 15.0 s
dt        = 0.002;      % Timestep
N_steps   = ceil(totalTime / dt) + 1;

%% Pre-allocate Logging History Arrays
hist_t             = zeros(1, N_steps);
hist_Pc            = zeros(1, N_steps); % Chamber pressure (psi)
hist_Tc            = zeros(1, N_steps); % Chamber temperature (K)
hist_cstar         = zeros(1, N_steps); % Characteristic velocity (m/s)
hist_isLit         = zeros(1, N_steps); % Combustion flag
hist_mdot_ox       = zeros(1, N_steps); % Ox injector mass flow (kg/s)
hist_mdot_fu       = zeros(1, N_steps); % Fuel injector mass flow (kg/s)
hist_mdot_noz      = zeros(1, N_steps); % Sonic throat nozzle mass flow (kg/s)
hist_mdot_v_ox     = zeros(1, N_steps); % Throttle valve OX flow (kg/s)
hist_mdot_v_fu     = zeros(1, N_steps); % Throttle valve FU flow (kg/s)
hist_OF            = zeros(1, N_steps); % O/F mixture ratio
hist_P_ox_man      = zeros(1, N_steps); % OX manifold pressure (psi)
hist_P_fu_reg      = zeros(1, N_steps); % FU regen inlet pressure (psi)
hist_P_fu_man      = zeros(1, N_steps); % FU manifold pressure (psi)
hist_cmd_ox        = zeros(1, N_steps); % Commanded Valve_OX State (0..1)
hist_cmd_fu        = zeros(1, N_steps); % Commanded Valve_FU State (0..1)
hist_theta_ox      = zeros(1, N_steps); % Commanded Valve_OX Angle (deg)
hist_theta_fu      = zeros(1, N_steps); % Commanded Valve_FU Angle (deg)
hist_ign           = zeros(1, N_steps); % Commanded IGN trigger (0/1)

%% Master Simulation Loop
t_sim = 0.0;
step  = 0;
CmdStates = System.Link.State;

while t_sim <= totalTime + 1e-9
    step = step + 1;

    %% Controls Hookup
    [theta_ox_deg, theta_fu_deg, ign_cmd] = ActuatorController(t_sim, State);
    cmd_ox_state = valveAngleToState(theta_ox_deg);
    cmd_fu_state = valveAngleToState(theta_fu_deg);

    % Set internal states to commands.
    CmdStates.Valve_OX = cmd_ox_state;
    CmdStates.Valve_FU = cmd_fu_state;
    CmdStates.IGN      = (ign_cmd > 0);

    % Ignition trigger interface
    State.Nodes.SKIPPER.SparkActive = (ign_cmd > 0);
    if ign_cmd == 0 && t_sim >= 10.0
        State.Nodes.SKIPPER.isLit = false; 
    end

    %% Physics Simulation Step
    [State, FlowRates, System] = StepSimulation(State, dt, System, CmdStates);

    %% Extract & Calculate Telemetry
    Pc_psi       = State.Nodes.SKIPPER.P / psi2Pa;
    Tc_K         = State.Nodes.SKIPPER.T;
    cstar_val    = State.Nodes.SKIPPER.cstar;
    isLit_val    = double(State.Nodes.SKIPPER.isLit);

    mdot_ox_val  = max(0.0, FlowRates.Inj_OX);
    mdot_fu_val  = max(0.0, FlowRates.Inj_FU);
    mdot_noz_val = max(0.0, FlowRates.Nozzle);
    mdot_vox_val = max(0.0, FlowRates.Valve_OX);
    mdot_vfu_val = max(0.0, FlowRates.Valve_FU);

    if mdot_fu_val > 1e-4
        OF_val = mdot_ox_val / mdot_fu_val;
    else
        OF_val = NaN;
    end

    P_ox_man_psi = State.Nodes.OX_Manifold.P / psi2Pa;
    P_fu_reg_psi = State.Nodes.FU_Regen_Inlet.P / psi2Pa;
    P_fu_man_psi = State.Nodes.FU_Manifold.P / psi2Pa;

    %% Log Step Data
    hist_t(step)         = t_sim;
    hist_Pc(step)        = Pc_psi;
    hist_Tc(step)        = Tc_K;
    hist_cstar(step)     = cstar_val;
    hist_isLit(step)     = isLit_val;
    hist_mdot_ox(step)   = mdot_ox_val;
    hist_mdot_fu(step)   = mdot_fu_val;
    hist_mdot_noz(step)  = mdot_noz_val;
    hist_mdot_v_ox(step) = mdot_vox_val;
    hist_mdot_v_fu(step) = mdot_vfu_val;
    hist_OF(step)        = OF_val;
    hist_P_ox_man(step)  = P_ox_man_psi;
    hist_P_fu_reg(step)  = P_fu_reg_psi;
    hist_P_fu_man(step)  = P_fu_man_psi;
    hist_cmd_ox(step)    = cmd_ox_state;
    hist_cmd_fu(step)    = cmd_fu_state;
    hist_theta_ox(step)  = theta_ox_deg;
    hist_theta_fu(step)  = theta_fu_deg;
    hist_ign(step)       = ign_cmd;

    %% Live Console Update
    if mod(step - 1, 100) == 0 || step == 1 || step == N_steps
        of_str = sprintf('%5.2f', OF_val);
        if isnan(OF_val), of_str = '  N/A'; end
        fprintf('  [t = %5.2f s] Pc = %5.1f psi | Tc = %6.1f K | OF = %s | mdot_ox = %5.3f kg/s | mdot_fu = %5.3f kg/s | mdot_noz = %5.3f kg/s | Lit = %d\n', ...
            t_sim, Pc_psi, Tc_K, of_str, mdot_ox_val, mdot_fu_val, mdot_noz_val, isLit_val);
    end

    t_sim = t_sim + dt;
end

% Trim arrays to executed steps
idx_valid = 1:step;
hist_t         = hist_t(idx_valid);
hist_Pc        = hist_Pc(idx_valid);
hist_Tc        = hist_Tc(idx_valid);
hist_cstar     = hist_cstar(idx_valid);
hist_isLit     = hist_isLit(idx_valid);
hist_mdot_ox   = hist_mdot_ox(idx_valid);
hist_mdot_fu   = hist_mdot_fu(idx_valid);
hist_mdot_noz  = hist_mdot_noz(idx_valid);
hist_mdot_v_ox = hist_mdot_v_ox(idx_valid);
hist_mdot_v_fu = hist_mdot_v_fu(idx_valid);
hist_OF        = hist_OF(idx_valid);
hist_P_ox_man  = hist_P_ox_man(idx_valid);
hist_P_fu_reg  = hist_P_fu_reg(idx_valid);
hist_P_fu_man  = hist_P_fu_man(idx_valid);
hist_cmd_ox    = hist_cmd_ox(idx_valid);
hist_cmd_fu    = hist_cmd_fu(idx_valid);
hist_theta_ox  = hist_theta_ox(idx_valid);
hist_theta_fu  = hist_theta_fu(idx_valid);
hist_ign       = hist_ign(idx_valid);
fprintf('\nSimulation completed (%d steps).\n', step);

%% Multi-Panel Diagnostic Plotting
fig = figure('Color', 'w', 'Position', [60, 60, 1400, 920], 'Visible', 'on');

% --- Subplot 1: Chamber Pressure Pc vs. Time ---
subplot(3, 2, 1);
plot(hist_t, hist_Pc, 'b-', 'LineWidth', 2.0); hold on;
yline(550.0, 'k--', 'Supply Tank (550 psi)', 'LineWidth', 1.2, 'FontSize', 8, 'LabelHorizontalAlignment', 'left');
yline(14.7, 'k:', 'Atmosphere (14.7 psi)', 'LineWidth', 1.0, 'FontSize', 8, 'LabelHorizontalAlignment', 'left');
xline(5.0, 'g--', 'Ignition (5.0s)', 'LineWidth', 1.2, 'FontSize', 8);
xline(10.0, 'r--', 'Cutoff (10.0s)', 'LineWidth', 1.2, 'FontSize', 8);
grid on; box on;
xlim([0, totalTime]);
ylim([0, 600]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Chamber Pressure Pc (psia)', 'FontSize', 9, 'FontWeight', 'bold');
title('SKIPPER Chamber Pressure Response', 'FontSize', 11, 'FontWeight', 'bold');

% --- Subplot 2: Mass Flow Rates vs. Time ---
subplot(3, 2, 2);
plot(hist_t, hist_mdot_ox, 'b-', 'LineWidth', 1.8); hold on;
plot(hist_t, hist_mdot_fu, 'r-', 'LineWidth', 1.8);
plot(hist_t, hist_mdot_noz, 'k-', 'LineWidth', 2.0);
xline(5.0, 'g--', 'LineWidth', 1.0);
xline(10.0, 'r--', 'LineWidth', 1.0);
grid on; box on;
xlim([0, totalTime]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Mass Flow Rate (kg/s)', 'FontSize', 9, 'FontWeight', 'bold');
title('Propellant Injector & Throat Mass Flow Rates', 'FontSize', 11, 'FontWeight', 'bold');
legend('LOX Injector (\itmdot_{ox}\rm)', 'Fuel Injector (\itmdot_{fu}\rm)', 'Nozzle Throat (\itmdot_{noz}\rm)', 'Location', 'northeast');

% --- Subplot 3: Feed Line Pressures (Delta-P Cascades) ---
subplot(3, 2, 3);
plot(hist_t, hist_P_ox_man, 'b-', 'LineWidth', 1.8); hold on;
plot(hist_t, hist_P_fu_reg, 'Color', [0.85, 0.45, 0.1], 'LineWidth', 1.8);
plot(hist_t, hist_P_fu_man, 'r-', 'LineWidth', 1.8);
plot(hist_t, hist_Pc, 'k--', 'LineWidth', 1.5);
yline(550.0, 'k:', '550 psi Supply', 'LineWidth', 1.0, 'FontSize', 8);
xline(5.0, 'g--', 'LineWidth', 1.0);
xline(10.0, 'r--', 'LineWidth', 1.0);
grid on; box on;
xlim([0, totalTime]);
ylim([0, 600]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Line Pressures (psia)', 'FontSize', 9, 'FontWeight', 'bold');
title('Feed Line & Manifold Pressures (Hydraulic Drops)', 'FontSize', 11, 'FontWeight', 'bold');
legend('P_{OX,Manifold}', 'P_{FU,RegenIn}', 'P_{FU,Manifold}', 'P_c (Chamber)', 'Location', 'southeast', 'FontSize', 8);

% --- Subplot 4: Hydraulic Delta-P Drops Across Components ---
subplot(3, 2, 4);
dp_ox_valve = max(0, 550.0 - hist_P_ox_man);
dp_ox_inj   = max(0, hist_P_ox_man - hist_Pc);
dp_fu_valve = max(0, 550.0 - hist_P_fu_reg);
dp_fu_regen = max(0, hist_P_fu_reg - hist_P_fu_man);
dp_fu_inj   = max(0, hist_P_fu_man - hist_Pc);

plot(hist_t, dp_ox_valve, 'b--', 'LineWidth', 1.5); hold on;
plot(hist_t, dp_ox_inj,   'b-',  'LineWidth', 2.0);
plot(hist_t, dp_fu_valve, 'r--', 'LineWidth', 1.5);
plot(hist_t, dp_fu_regen, 'Color', [0.85, 0.45, 0.1], 'LineWidth', 1.8);
plot(hist_t, dp_fu_inj,   'r-',  'LineWidth', 2.0);
xline(5.0, 'g--', 'LineWidth', 1.0);
xline(10.0, 'r--', 'LineWidth', 1.0);
grid on; box on;
xlim([0, totalTime]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('\DeltaP (psi)', 'FontSize', 9, 'FontWeight', 'bold');
title('Component Pressure Drop Breakdown (\DeltaP)', 'FontSize', 11, 'FontWeight', 'bold');
legend('\DeltaP OX Valve', '\DeltaP OX Inj', ...
       '\DeltaP FU Valve', '\DeltaP FU Regen', '\DeltaP FU Inj', ...
       'Location', 'northeast', 'FontSize', 8);

% --- Subplot 5: Chamber Gas Temperature Tc vs. Time ---
subplot(3, 2, 5);
plot(hist_t, hist_Tc, 'r-', 'LineWidth', 2.0); hold on;
yline(293.15, 'k:', 'Ambient (293.15 K)', 'LineWidth', 1.0, 'FontSize', 8, 'LabelHorizontalAlignment', 'left');
xline(5.0, 'g--', 'Ignition', 'LineWidth', 1.0);
xline(10.0, 'r--', 'Cutoff', 'LineWidth', 1.0);
grid on; box on;
xlim([0, totalTime]);
ylim([200, 3600]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
ylabel('Temperature Tc (K)', 'FontSize', 9, 'FontWeight', 'bold');
title('Chamber Gas Temperature Evolution', 'FontSize', 11, 'FontWeight', 'bold');

% --- Subplot 6: Mixture Ratio (O/F) & Actuator State Commands ---
subplot(3, 2, 6);
yyaxis left;
plot(hist_t, hist_OF, 'Color', [0.1, 0.6, 0.2], 'LineWidth', 2.0); hold on;
ylabel('Mixture Ratio (O/F)', 'FontSize', 9, 'FontWeight', 'bold');
ylim([0, 3.5]);

yyaxis right;
plot(hist_t, hist_cmd_ox, 'b-', 'LineWidth', 1.5); hold on;
plot(hist_t, hist_cmd_fu, 'r:', 'LineWidth', 1.5);
ylabel('Valve Command [0..1]', 'FontSize', 9, 'FontWeight', 'bold');
ylim([-0.05, 1.15]);

xline(5.0, 'g--', 'LineWidth', 1.0);
xline(10.0, 'r--', 'LineWidth', 1.0);
grid on; box on;
xlim([0, totalTime]);
xlabel('Time (s)', 'FontSize', 9, 'FontWeight', 'bold');
title('Mixture Ratio (O/F) & Valve Actuation Interface', 'FontSize', 11, 'FontWeight', 'bold');
legend('O/F Ratio', 'Target O/F', 'Valve\_OX Cmd', 'Valve\_FU Cmd', 'Location', 'northeast', 'FontSize', 8);

sgtitle('SKIPPER Test Stand A - Universal 15.0 s Hot-Fire Dynamic Response', 'FontSize', 14, 'FontWeight', 'bold');

%% ======================================================================
%% Modular Valve Actuation & Controller Interface
%% ======================================================================

function [theta_ox_deg, theta_fu_deg, ign_cmd] = ActuatorController(t_sim, State)

    if t_sim >= 5.0 && t_sim < 10.0
        % Phase 2: Hot-Fire Burn (5.0s to 10.0s)
        % Tuned to match: Pc = 250 psi, OF = 1.20, and Regen dP = 155 psi
        theta_ox_deg = 13.62; % Commanded angle for OX valve [deg]
        theta_fu_deg = 20.03; % Commanded angle for FU valve [deg]
        ign_cmd      = 1;     % Energize igniter / enable combustion
    else
        % Phase 1 (0.0s to 5.0s): Cold Dwell
        % Phase 3 (10.0s to 15.0s): Cutoff Dwell
        theta_ox_deg = 0.0;  % Fully closed (0 deg)
        theta_fu_deg = 0.0;  % Fully closed (0 deg)
        ign_cmd      = 0;    % Cutoff / de-energize
    end
end

function valve_state = valveAngleToState(angle_deg)
    valve_state = max(0.0, min(1.0, angle_deg / 90.0));
end
