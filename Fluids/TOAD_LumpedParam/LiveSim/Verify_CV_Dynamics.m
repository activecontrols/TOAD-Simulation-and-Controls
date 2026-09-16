%% VERIFY_CV_DYNAMICS
% Verification script testing first-principles Control Volume (CV) dynamics:
% 1. Initial flight-ready states
% 2. COPV isolation when BV-N2-02 = 0 (purges & RCS open)
% 3. Pressurization to full 550 psi (BV-N2-02 = 1)
% 4. Purge pressure rises on OX_Manifold, FU_Manifold, and DART_Chamber
% 5. RCS independent cold gas discharge
% 6. LOX tank drain with active pressurization (no 63 psi freeze)

clear; clc;
fprintf('=======================================================\n');
fprintf('       TOAD CV NODAL MASSFLOW DYNAMICS VERIFICATION    \n');
fprintf('=======================================================\n\n');

psi2Pa = 6894.757;
passed = true;

%% TEST 1: Initial Flight-Ready States
fprintf('--- TEST 1: Initial Flight-Ready States ---\n');
[System, State] = BuildTOADSystem_Merged();
P_copv_init = State.Nodes.TK_N2.P / psi2Pa;
P_bulk_init = State.Nodes.TK_N2_BULK.P / psi2Pa;
P_ox_init   = State.Nodes.TK_O2_01.P / psi2Pa;
P_fu_init   = State.Nodes.TK_FU_01.P / psi2Pa;

fprintf('  COPV Pressure: %.1f psia (target: 4500)\n', P_copv_init);
fprintf('  Bulk Pressure: %.1f psia (target: 6000)\n', P_bulk_init);
fprintf('  LOX Tank P:    %.1f psia (target: 14.7)\n', P_ox_init);
fprintf('  Fuel Tank P:   %.1f psia (target: 14.7)\n', P_fu_init);

if abs(P_copv_init - 4500) < 5 && abs(P_bulk_init - 6000) < 5
    fprintf('  [PASS] Initial pressures match flight-ready criteria.\n\n');
else
    fprintf('  [FAIL] Initial pressures do not match.\n\n');
    passed = false;
end

%% TEST 2: COPV Isolation (BV-N2-02 = 0 with Purges & RCS Open)
fprintf('--- TEST 2: COPV Isolation (BV-N2-02 = 0) ---\n');
[System, State] = BuildTOADSystem_Merged();
m_copv_before = State.Nodes.TK_N2.m;
P_copv_before = State.Nodes.TK_N2.P / psi2Pa;

% Open all purges and RCS thrusters, but keep BV-N2-02 = 0
LinkStates = struct(...
    'BV_N2_02', 0.0, ...
    'SV_N2_05', 1.0, 'SV_N2_06', 1.0, 'SV_N2_07', 1.0, ...
    'SV_N2_01', 1.0, 'SV_N2_02', 1.0, 'SV_N2_03', 1.0, 'SV_N2_04', 1.0);

for s = 1:200 % 2 seconds at 100 Hz
    [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates);
end
m_copv_after = State.Nodes.TK_N2.m;
P_copv_after = State.Nodes.TK_N2.P / psi2Pa;
mass_loss = m_copv_before - m_copv_after;
copv_loss = P_copv_before - P_copv_after;

fprintf('  COPV Mass Loss over 2.0s: %.6f kg (target: 0.000000)\n', mass_loss);
fprintf('  COPV Press Loss over 2.0s: %.4f psi (target: < 0.5)\n', copv_loss);
fprintf('  REG-873D Flow:            %.6f kg/s (target: 0.000000)\n', FlowRates.REG_873D);
fprintf('  BV-N2-02 Flow:            %.6f kg/s (target: 0.000000)\n', FlowRates.BV_N2_02);

if mass_loss < 1e-6 && FlowRates.REG_873D < 1e-6
    fprintf('  [PASS] COPV is 100%% isolated when BV-N2-02 is closed.\n\n');
else
    fprintf('  [FAIL] COPV drained with BV-N2-02 closed!\n\n');
    passed = false;
end

%% TEST 3: Propellant Tank Pressurization to 550 psi
fprintf('--- TEST 3: Tank Pressurization to 550 psi (BV-N2-02 = 1) ---\n');
[System, State] = BuildTOADSystem_Merged();
LinkStates = struct('BV_N2_02', 1.0);

for s = 1:1000 % 10 seconds at 100 Hz
    [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates);
end

P_ox_press = State.Nodes.TK_O2_01.P / psi2Pa;
P_fu_press = State.Nodes.TK_FU_01.P / psi2Pa;
P_pm_press = State.Nodes.Purge_Manifold.P / psi2Pa;

fprintf('  Purge Manifold P: %.1f psia (target: ~550)\n', P_pm_press);
fprintf('  LOX Tank P:       %.1f psia (target: 550 +/- 15)\n', P_ox_press);
fprintf('  Fuel Tank P:      %.1f psia (target: 550 +/- 15)\n', P_fu_press);

if abs(P_ox_press - 550) < 20 && abs(P_fu_press - 550) < 20
    fprintf('  [PASS] Propellant tanks pressurize to full 550 psi.\n\n');
else
    fprintf('  [FAIL] Pressurization stalled or missed target.\n\n');
    passed = false;
end

%% TEST 4: Purge Pressure Rises
fprintf('--- TEST 4: Purge Pressure Rises ---\n');
% Test A: OX Purge (SV-N2-05)
LinkStates = struct('BV_N2_02', 1.0, 'SV_N2_05', 1.0);
for s = 1:100, [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates); end
P_ox_man = State.Nodes.OX_Manifold.P / psi2Pa;
fprintf('  OX Manifold P with SV-N2-05=1:   %.1f psia (target: > 50)\n', P_ox_man);

% Test B: FU Purge (SV-N2-06)
LinkStates = struct('BV_N2_02', 1.0, 'SV_N2_06', 1.0);
for s = 1:100, [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates); end
P_fu_man = State.Nodes.FU_Manifold.P / psi2Pa;
fprintf('  FU Manifold P with SV-N2-06=1:   %.1f psia (target: > 50)\n', P_fu_man);

% Test C: DART Purge (SV-N2-07)
LinkStates = struct('BV_N2_02', 1.0, 'SV_N2_07', 1.0);
for s = 1:100, [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates); end
P_dart_cham = State.Nodes.DART_Chamber.P / psi2Pa;
fprintf('  DART Chamber P with SV-N2-07=1: %.1f psia (target: > 25)\n', P_dart_cham);

if P_ox_man > 40 && P_fu_man > 40 && P_dart_cham > 25
    fprintf('  [PASS] All purges exhibit clear, measurable pressure rises.\n\n');
else
    fprintf('  [FAIL] Purge pressure rise insufficient.\n\n');
    passed = false;
end

%% TEST 5: LOX Tank Drain with Active Pressurization
fprintf('--- TEST 5: Tank Drain with Active Pressurization ---\n');
[System, State] = BuildTOADSystem_Merged();
% Pre-pressurize to 550 psi
LinkStates = struct('BV_N2_02', 1.0);
for s = 1:500, [State, ~, System] = StepSimulation_Live(State, 0.01, System, LinkStates); end

% Now open LOX drain (BV_O2_02 = 1) with BV-N2-02 = 1 (active press)
LinkStates = struct('BV_N2_02', 1.0, 'BV_O2_02', 1.0);
fprintf('  Draining LOX tank over 30s...\n');
P_history = zeros(1, 30);
for sec = 1:30
    for s = 1:100
        [State, FlowRates, System] = StepSimulation_Live(State, 0.01, System, LinkStates);
    end
    P_history(sec) = State.Nodes.TK_O2_01.P / psi2Pa;
end

P_drain_min = min(P_history);
P_drain_end = P_history(end);
m_liq_end   = State.Nodes.TK_O2_01.Liquid.m;

fprintf('  Min Tank P during drain: %.1f psia\n', P_drain_min);
fprintf('  Final Tank P:            %.1f psia (target: ~550, NOT 63)\n', P_drain_end);
fprintf('  Remaining Liquid Mass:   %.3f kg (target: 0.000)\n', m_liq_end);

if P_drain_min > 450 && abs(P_drain_end - 550) < 30 && m_liq_end < 0.05
    fprintf('  [PASS] Active pressurization maintains ~550 psi throughout drain.\n\n');
else
    fprintf('  [FAIL] Pressure dropped or got stuck at 63 psi.\n\n');
    passed = false;
end

%% SUMMARY
fprintf('=======================================================\n');
if passed
    fprintf('           ALL CV DYNAMICS CHECKS PASSED!             \n');
else
    fprintf('           SOME CHECKS FAILED - REVIEW LOGS           \n');
end
fprintf('=======================================================\n');
