%% PressCode_Run.m: TOAD Vehicle Flight Pressurization
% Master script for flight press analysis.
% Integrates FlightProfile_1DoF trajectory with a 4 Node Fluid Network &
% Thermal Network
%
% Pablo Plata

clear; clc; close all;
fprintf('=================================================================\n');
fprintf('   TOAD Vehicle Flight Pressurization & N2 Consumption Analysis   \n');
fprintf('=================================================================\n\n');

%% Path Management
scriptDir = fileparts(mfilename('fullpath'));
candidates = {
    fullfile(scriptDir, '..', 'Core'), ...
    fullfile(scriptDir, '..', 'Thermal'), ...
    fullfile(scriptDir, '..', 'Helpers'), ...
    fullfile(scriptDir, '..', 'Config'), ...
    fullfile(scriptDir, '..', 'Fluid Properties'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Core'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Thermal'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Helpers'), ...
    fullfile(pwd, 'sandbox', 'experiments', 'Fluid Properties'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Core'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Fluid Properties'), ...
    fullfile(pwd, 'TOAD_LumpedParam', 'Helpers')
};
for c = 1:numel(candidates)
    if exist(candidates{c}, 'dir'), addpath(candidates{c}); end
end

%% Flight Trajectory & Outflow Profile Generation
fprintf('Generating 1-DoF Flight Profile...\n');
max_isp          = 185.631;  % Max throttle ISP (s) 
min_isp          = 166.629;  % Min throttle ISP (s)
total_OF         = 1.2;      % Total OF Ratio
min_throttle     = 0.5;      % Minimum throttle fraction
max_thrust       = 660.0;    % Max thrust (lbf)
TWR              = 1.25;     % Thrust to weight ratio
Propellant_mass  = 127.0;    % Propellant mass (lbm)

TOAD_mass = max_thrust / TWR;
selection = 1;

[fuel_mass_lbm, ox_mass_lbm, totalTime, fuel_mdot_lbm, ox_mdot_lbm] = ...
    FlightProfile_1DoF(TOAD_mass, total_OF, max_isp, min_isp, min_throttle, ...
                       max_thrust, Propellant_mass, selection);
close all; % Suppress initial 1DoF standalone figures

lbm_to_kg = 0.45359237;
kg_to_lbm = 1.0 / lbm_to_kg;

ox_mdot_kg   = ox_mdot_lbm * lbm_to_kg;
fuel_mdot_kg = fuel_mdot_lbm * lbm_to_kg;
N_pts = length(ox_mdot_kg);
profile_time = linspace(0.0, totalTime, N_pts);

fprintf('Trajectory generated: Flight Duration = %.2f s, Target LOX = %.2f kg (%.1f lbm), Target Fuel = %.2f kg (%.1f lbm)\n\n', ...
    totalTime, sum(ox_mdot_kg) * (totalTime / N_pts), sum(ox_mdot_lbm) * (totalTime / N_pts), ...
    sum(fuel_mdot_kg) * (totalTime / N_pts), sum(fuel_mdot_lbm) * (totalTime / N_pts));

%% System Initialization 
[System, State] = TOADPressCode();

psi2Pa = 6894.757;
P_atm  = System.Env.P_atm;
dt     = 0.010; % 100 Hz simulation rate

% Standard conditions for SCFM (N2 at 14.696 psia, 60 deg F)
R_spec_n2   = 296.8; % J/(kg*K)
T_std       = 288.71; % K (60 F)
P_std       = 101325.0; % Pa
rho_std     = P_std / (R_spec_n2 * T_std); % ~1.182 kg/m^3
m3s_to_scfm = 2118.88; % 1 m^3/s = 2118.88 ft^3/min

%% Phase 1: Pad Pre-Pressurization & Thermal Settling
fprintf('-----------------------------------------------------------------\n');
fprintf('Phase 1: Pre-Pressurizing Propellant Tanks to 550 psia + Thermal Settling Delay...\n');

m_copv_initial = State.Nodes.TK_N2.m;
t_sim = 0.0;

CmdStates = System.Link.State;
CmdStates.REG_873D = 1.0;
CmdStates.Press_OX = 1.0;
CmdStates.Press_FU = 1.0;
CmdStates.OX_Drain = 0.0;
CmdStates.FU_Drain = 0.0;

State.mdot_pr_ox   = 0.0;
State.mdot_pr_fu   = 0.0;
State.COPV_Flowing = false;

% Step until both tanks reach operating pressure (>= 540 psi), then hold for settling
isPrepressed = false;
t_prepress_done = 0.0;
t_settle_hold   = 150.0; 

max_pad_steps = ceil(20.0 / dt);
log_pad_t        = zeros(1, max_pad_steps);
log_pad_P_copv   = zeros(1, max_pad_steps);
log_pad_P_ox     = zeros(1, max_pad_steps);
log_pad_P_fu     = zeros(1, max_pad_steps);
log_pad_P_man    = zeros(1, max_pad_steps);
log_pad_T_ox_ull = zeros(1, max_pad_steps);
log_pad_T_fu_ull = zeros(1, max_pad_steps);
log_pad_T_copv   = zeros(1, max_pad_steps);
log_pad_flow_ox  = zeros(1, max_pad_steps);
log_pad_flow_fu  = zeros(1, max_pad_steps);
log_pad_m_copv   = zeros(1, max_pad_steps);
pad_step = 0;

while true
    pad_step = pad_step + 1;
    t_sim = t_sim + dt;

    % Thermal step
    [~, State] = StepThermalNetwork(State, dt, System.Env);

    % Fluid flow & CV state step
    [State, FlowRates, System] = StepSimulation(State, dt, System, CmdStates);

    % Update flows for thermal network
    pr_ox_flow = 0.0; if isfield(FlowRates, 'Press_OX'), pr_ox_flow = max(0.0, FlowRates.Press_OX); end
    pr_fu_flow = 0.0; if isfield(FlowRates, 'Press_FU'), pr_fu_flow = max(0.0, FlowRates.Press_FU); end
    reg_flow   = 0.0; if isfield(FlowRates, 'REG_873D'), reg_flow   = max(0.0, FlowRates.REG_873D); end

    State.mdot_pr_ox   = pr_ox_flow;
    State.mdot_pr_fu   = pr_fu_flow;
    State.COPV_Flowing = (reg_flow > 1e-4);

    % Log pad transient
    log_pad_t(pad_step)        = t_sim;
    log_pad_P_copv(pad_step)   = State.Nodes.TK_N2.P / psi2Pa;
    log_pad_P_ox(pad_step)     = State.Nodes.TK_O2_01.P / psi2Pa;
    log_pad_P_fu(pad_step)     = State.Nodes.TK_FU_01.P / psi2Pa;
    log_pad_P_man(pad_step)    = State.Nodes.Press_Manifold.P / psi2Pa;
    log_pad_T_ox_ull(pad_step) = State.Nodes.TK_O2_01.Ullage.T;
    log_pad_T_fu_ull(pad_step) = State.Nodes.TK_FU_01.Ullage.T;
    log_pad_T_copv(pad_step)   = State.Nodes.TK_N2.T;
    log_pad_flow_ox(pad_step)  = pr_ox_flow * kg_to_lbm;
    log_pad_flow_fu(pad_step)  = pr_fu_flow * kg_to_lbm;
    log_pad_m_copv(pad_step)   = State.Nodes.TK_N2.m;

    % Check if 540 psi reached
    if ~isPrepressed && State.Nodes.TK_O2_01.P >= 540.0 * psi2Pa && State.Nodes.TK_FU_01.P >= 540.0 * psi2Pa
        isPrepressed = true;
        t_prepress_done = t_sim;
        fprintf('  Target pressure (540 psi) reached at t = %.2f s. Starting %.1f s thermal settling hold...\n', ...
            t_prepress_done, t_settle_hold);
    end

    if mod(pad_step, 100) == 0
        fprintf('  [Pad t = %5.2f s] P_COPV = %6.1f psi | P_OX = %5.1f psi | P_FU = %5.1f psi | T_OX_Ull = %5.1f K\n', ...
            t_sim, log_pad_P_copv(pad_step), log_pad_P_ox(pad_step), log_pad_P_fu(pad_step), log_pad_T_ox_ull(pad_step));
    end

    % Exit once settling delay is complete
    if isPrepressed && (t_sim >= t_prepress_done + t_settle_hold)
        t_pad_total = t_sim;
        break;
    end
end

m_copv_after_pad = State.Nodes.TK_N2.m;
N2_pad_used_kg   = m_copv_initial - m_copv_after_pad;
N2_pad_used_lbm  = N2_pad_used_kg * kg_to_lbm;

fprintf('Pad pre-press & settling complete at t = %.2f s.\n', t_pad_total);
fprintf('  Pre-press N2 used: %.3f kg (%.3f lbm)\n', N2_pad_used_kg, N2_pad_used_lbm);
fprintf('  Settled state: P_OX = %.1f psi, T_OX_Ull = %.1f K | P_FU = %.1f psi, T_FU_Ull = %.1f K | P_COPV = %.1f psi\n\n', ...
    State.Nodes.TK_O2_01.P / psi2Pa, State.Nodes.TK_O2_01.Ullage.T, ...
    State.Nodes.TK_FU_01.P / psi2Pa, State.Nodes.TK_FU_01.Ullage.T, ...
    State.Nodes.TK_N2.P / psi2Pa);

%% Phase 2: Flight Burn Profile Execution
fprintf('-----------------------------------------------------------------\n');
fprintf('Phase 2: Executing Flight Burn Profile (Duration: %.2f s, dt = %.3f s)...\n', totalTime, dt);

t_flight = 0.0;
step = 0;
N_flight_steps = ceil(totalTime / dt) + 1;

hist_t_fl          = zeros(1, N_flight_steps);
hist_P_copv        = zeros(1, N_flight_steps);
hist_P_ox          = zeros(1, N_flight_steps);
hist_P_fu          = zeros(1, N_flight_steps);
hist_P_pman        = zeros(1, N_flight_steps);
hist_T_ox_ull      = zeros(1, N_flight_steps);
hist_T_fu_ull      = zeros(1, N_flight_steps);
hist_T_copv        = zeros(1, N_flight_steps);
hist_m_copv        = zeros(1, N_flight_steps);
hist_m_ox_liq      = zeros(1, N_flight_steps);
hist_m_fu_liq      = zeros(1, N_flight_steps);
hist_mdot_ox_prop  = zeros(1, N_flight_steps);
hist_mdot_fu_prop  = zeros(1, N_flight_steps);
hist_mdot_n2_ox    = zeros(1, N_flight_steps);
hist_mdot_n2_fu    = zeros(1, N_flight_steps);
hist_mdot_n2_tot   = zeros(1, N_flight_steps);
hist_scfm_ox       = zeros(1, N_flight_steps);
hist_scfm_fu       = zeros(1, N_flight_steps);
hist_scfm_tot      = zeros(1, N_flight_steps);
hist_collapse_ox   = zeros(1, N_flight_steps);
hist_collapse_fu   = zeros(1, N_flight_steps);

while t_flight <= totalTime
    step = step + 1;

    % Target propellant outflow from 1-DoF trajectory
    mdot_tgt_ox = interp1(profile_time, ox_mdot_kg, min(t_flight, totalTime), 'linear', 'extrap');
    mdot_tgt_fu = interp1(profile_time, fuel_mdot_kg, min(t_flight, totalTime), 'linear', 'extrap');

    % Direct algebraic drain Cv modulation: mdot = Cv * 2.402e-5 * sqrt(rho * DP)
    % For LOX drain:
    P_tank_ox = State.Nodes.TK_O2_01.Liquid.P;
    DP_ox     = max(P_tank_ox - P_atm, 1.0 * psi2Pa);
    rho_ox    = max(State.Nodes.TK_O2_01.Liquid.rho, 900.0);
    Cv_ox_req = mdot_tgt_ox / (2.402e-5 * sqrt(rho_ox * DP_ox));
    CmdStates.OX_Drain = max(0.0, min(1.0, Cv_ox_req / System.Links.OX_Drain.MaxCv));

    % For Fuel drain:
    P_tank_fu = State.Nodes.TK_FU_01.Liquid.P;
    DP_fu     = max(P_tank_fu - P_atm, 1.0 * psi2Pa);
    rho_fu    = max(State.Nodes.TK_FU_01.Liquid.rho, 750.0);
    Cv_fu_req = mdot_tgt_fu / (2.402e-5 * sqrt(rho_fu * DP_fu));
    CmdStates.FU_Drain = max(0.0, min(1.0, Cv_fu_req / System.Links.FU_Drain.MaxCv));

    % Thermal network step
    [~, State] = StepThermalNetwork(State, dt, System.Env);

    % Fluid flow & CV state step
    [State, FlowRates, System] = StepSimulation(State, dt, System, CmdStates);

    % Update pressurant flow rates for thermal model
    pr_ox_flow = max(0.0, FlowRates.Press_OX);
    pr_fu_flow = max(0.0, FlowRates.Press_FU);
    reg_flow   = max(0.0, FlowRates.REG_873D);

    State.mdot_pr_ox   = pr_ox_flow;
    State.mdot_pr_fu   = pr_fu_flow;
    State.COPV_Flowing = (reg_flow > 1e-4);

    actual_ox_drain = FlowRates.OX_Drain;
    actual_fu_drain = FlowRates.FU_Drain;

    % Pressurant flows in lbm/s
    pr_ox_lbm  = pr_ox_flow * kg_to_lbm;
    pr_fu_lbm  = pr_fu_flow * kg_to_lbm;
    pr_tot_lbm = pr_ox_lbm + pr_fu_lbm;

    % Pressurant SCFM
    scfm_ox_val  = (pr_ox_flow / rho_std) * m3s_to_scfm;
    scfm_fu_val  = (pr_fu_flow / rho_std) * m3s_to_scfm;
    scfm_tot_val = scfm_ox_val + scfm_fu_val;

    % Collapse Factors: C = mdot_pr / mdot_pr_ideal
    % where mdot_pr_ideal = rho_N2(P_tank, T_inlet) * Vdot_displaced
    T_inlet = State.Nodes.TK_N2.T; % Nominal ambient N2 supply temperature [K]
    props_id_ox = FluidProperties('Nitrogen', 'From_P_T', P_tank_ox, T_inlet);
    props_id_fu = FluidProperties('Nitrogen', 'From_P_T', P_tank_fu, T_inlet);
    rho_n2_ox_ideal = props_id_ox.rho;
    rho_n2_fu_ideal = props_id_fu.rho;

    vdot_ox_displaced = actual_ox_drain / rho_ox;
    vdot_fu_displaced = actual_fu_drain / rho_fu;

    mdot_ideal_ox = rho_n2_ox_ideal * vdot_ox_displaced;
    mdot_ideal_fu = rho_n2_fu_ideal * vdot_fu_displaced;

    % Guard against division-by-zero before valve opens
    if mdot_ideal_ox > 1e-4
        c_ox = pr_ox_flow / mdot_ideal_ox;
    else
        c_ox = 1.0;
    end
    
    if mdot_ideal_fu > 1e-4
        c_fu = pr_fu_flow / mdot_ideal_fu;
    else
        c_fu = 1.0;
    end

    % Log flight burn
    hist_t_fl(step)         = t_flight;
    hist_P_copv(step)       = State.Nodes.TK_N2.P / psi2Pa;
    hist_P_ox(step)         = State.Nodes.TK_O2_01.P / psi2Pa;
    hist_P_fu(step)         = State.Nodes.TK_FU_01.P / psi2Pa;
    hist_P_pman(step)       = State.Nodes.Press_Manifold.P / psi2Pa;
    hist_T_ox_ull(step)     = State.Nodes.TK_O2_01.Ullage.T;
    hist_T_fu_ull(step)     = State.Nodes.TK_FU_01.Ullage.T;
    hist_T_copv(step)       = State.Nodes.TK_N2.T;
    hist_m_copv(step)       = State.Nodes.TK_N2.m;
    hist_m_ox_liq(step)     = State.Nodes.TK_O2_01.Liquid.m;
    hist_m_fu_liq(step)     = State.Nodes.TK_FU_01.Liquid.m;
    hist_mdot_ox_prop(step) = actual_ox_drain * kg_to_lbm;
    hist_mdot_fu_prop(step) = actual_fu_drain * kg_to_lbm;
    hist_mdot_n2_ox(step)   = pr_ox_lbm;
    hist_mdot_n2_fu(step)   = pr_fu_lbm;
    hist_mdot_n2_tot(step)  = pr_tot_lbm;
    hist_scfm_ox(step)      = scfm_ox_val;
    hist_scfm_fu(step)      = scfm_fu_val;
    hist_scfm_tot(step)     = scfm_tot_val;
    hist_collapse_ox(step)  = c_ox;
    hist_collapse_fu(step)  = c_fu;

    if mod(step, 200) == 0 || step == 1
        fprintf('  [Flight t = %5.2f s / %5.2f s] P_COPV = %6.1f psi | P_OX = %5.1f psi | T_OX_Ull = %5.1f K | N2_Mdot = %5.3f lb/s\n', ...
            t_flight, totalTime, hist_P_copv(step), hist_P_ox(step), hist_T_ox_ull(step), pr_tot_lbm);
    end

    t_flight = t_flight + dt;
    t_sim    = t_sim + dt;
end

% Trim flight history arrays
idx_valid = 1:step;
hist_t_fl         = hist_t_fl(idx_valid);
hist_P_copv       = hist_P_copv(idx_valid);
hist_P_ox         = hist_P_ox(idx_valid);
hist_P_fu         = hist_P_fu(idx_valid);
hist_P_pman       = hist_P_pman(idx_valid);
hist_T_ox_ull     = hist_T_ox_ull(idx_valid);
hist_T_fu_ull     = hist_T_fu_ull(idx_valid);
hist_T_copv       = hist_T_copv(idx_valid);
hist_m_copv       = hist_m_copv(idx_valid);
hist_m_ox_liq     = hist_m_ox_liq(idx_valid);
hist_m_fu_liq     = hist_m_fu_liq(idx_valid);
hist_mdot_ox_prop = hist_mdot_ox_prop(idx_valid);
hist_mdot_fu_prop = hist_mdot_fu_prop(idx_valid);
hist_mdot_n2_ox   = hist_mdot_n2_ox(idx_valid);
hist_mdot_n2_fu   = hist_mdot_n2_fu(idx_valid);
hist_mdot_n2_tot  = hist_mdot_n2_tot(idx_valid);
hist_scfm_ox      = hist_scfm_ox(idx_valid);
hist_scfm_fu      = hist_scfm_fu(idx_valid);
hist_scfm_tot     = hist_scfm_tot(idx_valid);
hist_collapse_ox  = hist_collapse_ox(idx_valid);
hist_collapse_fu  = hist_collapse_fu(idx_valid);

%% Final Results & N2 Consumption Summary
m_copv_final       = State.Nodes.TK_N2.m;
N2_flight_used_kg  = m_copv_after_pad - m_copv_final;
N2_flight_used_lbm = N2_flight_used_kg * kg_to_lbm;
N2_total_used_kg   = m_copv_initial - m_copv_final;
N2_total_used_lbm  = N2_total_used_kg * kg_to_lbm;

fprintf('\n=================================================================\n');
fprintf('                     FLIGHT RESULTS SUMMARY                      \n');
fprintf('=================================================================\n');
fprintf('Initial COPV Pressure:             %.1f psia\n', log_pad_P_copv(1));
fprintf('COPV at Liftoff:                   %.1f psia\n', hist_P_copv(1));
fprintf('Final COPV Pressure:               %.1f psia (Delta_burn = %.1f psia)\n', hist_P_copv(end), hist_P_copv(1) - hist_P_copv(end));
fprintf('Initial COPV Temperature:          %.2f K\n', State.Thermal.COPV.T_layers(1));
fprintf('Final COPV Temperature:            %.2f K\n', hist_T_copv(end));
fprintf('Pad Pre-Press & Settling Time:     %.2f s\n', t_pad_total);
fprintf('Total Flight Burn Time:            %.2f s\n', totalTime);
fprintf('-----------------------------------------------------------------\n');
fprintf('N2 Consumed during Pad Pre-Press:  %.3f kg  (%.3f lbm)\n', N2_pad_used_kg, N2_pad_used_lbm);
fprintf('N2 Consumed during Flight Burn:    %.3f kg  (%.3f lbm)\n', N2_flight_used_kg, N2_flight_used_lbm);
fprintf('TOTAL VEHICLE N2 CONSUMPTION:      %.3f kg  (%.3f lbm)\n', N2_total_used_kg, N2_total_used_lbm);
fprintf('-----------------------------------------------------------------\n');
fprintf('LOX Ullage Temperature Range:      %.1f K (liftoff) to %.1f K (final)\n', hist_T_ox_ull(1), hist_T_ox_ull(end));
fprintf('Fuel Ullage Temperature Range:     %.1f K (liftoff) to %.1f K (final)\n', hist_T_fu_ull(1), hist_T_fu_ull(end));
fprintf('LOX Tank Regulated Pressure Range: %.1f psi to %.1f psi\n', min(hist_P_ox), max(hist_P_ox));
fprintf('Fuel Tank Regulated Pressure Range:%.1f psi to %.1f psi\n', min(hist_P_fu), max(hist_P_fu));
fprintf('Mean LOX Press Collapse Factor:    %.2f\n', mean(hist_collapse_ox(hist_t_fl > 1.0)));
fprintf('Mean Fuel Press Collapse Factor:   %.2f\n', mean(hist_collapse_fu(hist_t_fl > 1.0)));
fprintf('=================================================================\n');

%% 4-Panel Output Figure (Matching Baseline Format)
fig = figure('Color', 'w', 'Position', [60, 60, 1350, 880]);

% --- Subplot 1: COPV Pressure vs. Time ---
subplot(2, 2, 1);
plot(hist_t_fl, hist_P_copv, 'LineWidth', 1.8, 'Color', [0.15, 0.45, 0.8]);
ylabel('COPV Pressure (psia)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (sec.)', 'FontSize', 10);
title('COPV Pressure vs. Time', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0, ceil(totalTime)]);
ylim([floor(min(hist_P_copv)/100)*100 - 100, 4600]);
grid on;

% --- Subplot 2: Ullage Volume Temperature vs. Time ---
% subplot(2, 2, 2);
% plot(hist_t_fl, hist_T_fu_ull, 'r-', 'LineWidth', 1.8); hold on;
% plot(hist_t_fl, hist_T_ox_ull, 'b-', 'LineWidth', 1.8);
% ylabel('Ullage Volume Temperature (K)', 'FontSize', 10, 'FontWeight', 'bold');
% xlabel('Time (sec.)', 'FontSize', 10);
% title('Ullage Volume Temperature vs. Time', 'FontSize', 11, 'FontWeight', 'bold');
% legend('Fuel Ullage', 'Ox Ullage', 'Location', 'best', 'FontSize', 9);
% xlim([0, ceil(totalTime)]);
% ylim([70, max([hist_T_fu_ull, hist_T_ox_ull]) + 20]);
% grid on;

% --- Subplot 3: Pressurant SCFM vs. Time ---
subplot(2, 2, 2);
plot(hist_t_fl, hist_scfm_fu, 'r-', 'LineWidth', 1.6); hold on;
plot(hist_t_fl, hist_scfm_ox, 'b-', 'LineWidth', 1.6);
plot(hist_t_fl, hist_scfm_tot, 'g-', 'LineWidth', 1.8);
ylabel('Press Flow Rates (ft^3/min)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (sec.)', 'FontSize', 10);
title('Pressurant SCFM vs. Time', 'FontSize', 11, 'FontWeight', 'bold');
legend('Fuel Press', 'Ox Press', 'Combined', 'Location', 'best', 'FontSize', 9);
xlim([0, ceil(totalTime)]);
grid on;

% --- Subplot 4: COPV Temperature vs. Time ---
subplot(2, 2, 3);
plot(hist_t_fl, hist_T_copv, 'LineWidth', 1.8, 'Color', [0.15, 0.45, 0.8]);
ylabel('COPV Temperature (K)', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (sec.)', 'FontSize', 10);
title('COPV Temperature vs. Time', 'FontSize', 11, 'FontWeight', 'bold');
xlim([0, ceil(totalTime)]);
grid on;

% --- Subplot 6: Collapse Factor vs. Time ---
subplot(2, 2, 4);
c_fu_plot = min(hist_collapse_fu, 5.0);
c_ox_plot = min(hist_collapse_ox, 10.0);
plot(hist_t_fl, c_fu_plot, 'r-', 'LineWidth', 1.8); hold on;
plot(hist_t_fl, c_ox_plot, 'b-', 'LineWidth', 1.8);
ylabel('Tank Press Collapse Factors', 'FontSize', 10, 'FontWeight', 'bold');
xlabel('Time (sec.)', 'FontSize', 10);
title('Collapse Factor vs. Time', 'FontSize', 11, 'FontWeight', 'bold');
legend('Fuel Press', 'Ox Press', 'Location', 'best', 'FontSize', 9);
xlim([0, ceil(totalTime)]);
ylim([0.5, 8.5]);
grid on;

%% 8. Standalone 6-Panel Pre-Press Diagnostics Figure
fig_prepress = figure('Color', 'w', 'Position', [80, 80, 1100, 800]);

% Subplot 1: Pressure Evolution
subplot(2, 2, 1);
plot(log_pad_t(1:pad_step), log_pad_P_ox(1:pad_step), 'b-', 'LineWidth', 2.0); hold on;
plot(log_pad_t(1:pad_step), log_pad_P_fu(1:pad_step), 'r-', 'LineWidth', 2.0);
plot(log_pad_t(1:pad_step), log_pad_P_man(1:pad_step), 'm--', 'LineWidth', 1.0);
yline(550.0, 'k:', '550 psi Setpoint', 'LineWidth', 1.2, 'FontSize', 9);
if t_prepress_done > 0
    xline(t_prepress_done, 'g--', sprintf('540 psi @ %.2fs', t_prepress_done), 'LineWidth', 1.2, 'FontSize', 9);
end
grid on; box on;
xlabel('Time from Valve Open (s)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Pressure (psia)', 'FontSize', 10, 'FontWeight', 'bold');
title('Tank & Manifold Pressure Evolution (14.7 -> 550 psia)', 'FontSize', 11, 'FontWeight', 'bold');
legend('LOX Tank Ullage', 'Fuel Tank Ullage', 'Press Manifold', 'Location', 'best');
ylim([0, 600]);
xlim([0, 5 * t_prepress_done]);

% Subplot 2: Temperature Evolution
subplot(2, 2, 2);
plot(log_pad_t(1:pad_step), log_pad_T_ox_ull(1:pad_step), 'b-', 'LineWidth', 2.0); hold on;
plot(log_pad_t(1:pad_step), log_pad_T_fu_ull(1:pad_step), 'r-', 'LineWidth', 2.0);
plot(log_pad_t(1:pad_step), log_pad_T_copv(1:pad_step), 'k-.', 'LineWidth', 1.5);
yline(90.17, 'c--', 'LOX Wall Saturation (90.17 K)', 'LineWidth', 1.2, 'FontSize', 9);
yline(293.15, 'm:', 'Ambient (293.15 K)', 'LineWidth', 1.2, 'FontSize', 9);
if t_prepress_done > 0
    xline(t_prepress_done, 'g--', 'LineWidth', 1.2);
end
grid on; box on;
xlabel('Time from Valve Open (s)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Temperature (K)', 'FontSize', 10, 'FontWeight', 'bold');
title('Ullage & Supply Temperature Evolution', 'FontSize', 11, 'FontWeight', 'bold');
legend('LOX Ullage T', 'Fuel Ullage T', 'COPV Gas T', 'Location', 'best');
ylim([70, 500]);

% Subplot 3: Inflow Rates
subplot(2, 2, 3);
plot(log_pad_t(1:pad_step), log_pad_flow_ox(1:pad_step), 'b-', 'LineWidth', 1.8); hold on;
plot(log_pad_t(1:pad_step), log_pad_flow_fu(1:pad_step), 'r-', 'LineWidth', 1.8);
if t_prepress_done > 0
    xline(t_prepress_done, 'g--', 'LineWidth', 1.2);
end
grid on; box on;
xlabel('Time from Valve Open (s)', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('Inflow Rate (lbm/s)', 'FontSize', 10, 'FontWeight', 'bold');
title('Pressurant Flow Rates into Tanks', 'FontSize', 11, 'FontWeight', 'bold');
legend('LOX Check Valve Flow', 'Fuel Check Valve Flow', 'Location', 'best');
xlim([0, 5 * t_prepress_done]);

% Subplot 4: Cumulative N2 Consumed & COPV Pressure
subplot(2, 2, 4);
m_pad_n2_used = (m_copv_initial - log_pad_m_copv(1:pad_step)) * kg_to_lbm;
yyaxis left;
plot(log_pad_t(1:pad_step), m_pad_n2_used, 'k-', 'LineWidth', 2.0);
ylabel('N2 Used (lbm)', 'FontSize', 10, 'FontWeight', 'bold');
ylim([0, max(m_pad_n2_used)*1.2 + 0.1]);

yyaxis right;
plot(log_pad_t(1:pad_step), log_pad_P_copv(1:pad_step), 'b--', 'LineWidth', 1.8);
ylabel('COPV Pressure (psia)', 'FontSize', 10, 'FontWeight', 'bold');
ylim([3800, 4600]);

if t_prepress_done > 0
    xline(t_prepress_done, 'g--', 'LineWidth', 1.2);
end
grid on; box on;
xlabel('Time from Valve Open (s)', 'FontSize', 10, 'FontWeight', 'bold');
title('Pre-Press N2 Mass Used & COPV Blowdown', 'FontSize', 11, 'FontWeight', 'bold');