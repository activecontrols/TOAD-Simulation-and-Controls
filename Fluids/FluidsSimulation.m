%% This file is test Simulation run of the P&ID Model Simulation.
% Load in PID Structure
addpath('.\Fluids\Helpers');
addpath('.\Fluids\Links and Nodes');
PID_Structure;

%% Package initial state vector
PLength = size(System.Nodes, 2);
MDLength = size(System.Links.Dynamic, 2);
MALength = size(System.Links.Algebraic, 2);

YIdx_Start = PLength + MDLength + 1;
NumSpecies = 3;

% NEW: Enthalpy Index Tracking
HIdx_Start = YIdx_Start + (PLength * NumSpecies);

% Expand X0 to include Enthalpies (PLength nodes)
X0 = zeros(HIdx_Start + PLength - 1, 1);

% Package Pressures, Mass Fractions, and Initial Enthalpies
for i = 1:PLength
    T0 = System.Nodes(i).Temp;
    P0 = System.Nodes(i).P0 * 6895;
    Y0 = System.Nodes(i).Y0';
    
    X0(i) = P0;
    
    CurrentIDX = YIdx_Start + (i-1)*NumSpecies : YIdx_Start + (i-1)*NumSpecies + (NumSpecies-1);
    X0(CurrentIDX) = Y0;
    
    % --- SMART ENTHALPY INITIALIZATION ---
    [h_ox, ~, ~, ~, ~] = LOX_Data(T0, P0);
    [h_fu, ~, ~, ~, ~] = IPA_Data(T0, P0);
    [h_n2, ~, ~, ~, ~] = N2_Data(T0, P0);
    
    % Mass-weighted initial node enthalpy
    H0 = Y0(1)*h_ox + Y0(2)*h_fu + Y0(3)*h_n2;
    X0(HIdx_Start + i - 1) = H0;
end

% Package Dynamic Massflows 
for i = 1:MDLength
    X0(PLength + i) = 0;
end

%% Simulation Loop 2 (MATLAB Native ode15s Solver - Instantaneous Combustion)
simTime = 50;
dT = 0.05;
MaxSteps = round(simTime / dT);

% Logs
X_LOG = zeros(size(X0, 1), MaxSteps);
U_LOG = zeros(MALength, MaxSteps);
T_LOG = zeros(1, MaxSteps);
X_cur = X0;
t = 0;
T_LOG(1) = t;
X_LOG(:, 1) = X_cur;
stepCount = 1;

% Setup ode15s Options
AbsTol_Vec = ones(size(X0)) * 5e-4; % Base 1e-4 tolerance for mass flows and fractions
AbsTol_Vec(1:PLength) = 1000;       % Allow 1000 Pa (~0.15 psi) absolute tolerance for pressures
options = odeset('RelTol', 1e-2, 'AbsTol', AbsTol_Vec, 'MaxOrder', 2); 

% Warm-start
fprintf('Warming up system states...\n');
for k = 1:10
    XDOT = PhysicsEngine2(X_cur, System, zeros(MALength,1));
    X_cur = X_cur + XDOT * 5e-5;
end

fprintf('Warm up complete. Starting ode15s Loop.\n');
tic;

while t < simTime
    % Proposed next time
    t_next = t + dT;

    %% Valve Control Logic
    U = Scheduler.Step(t_next, dT);

    %% Integration Step via ode15s
    try
        [~, X_step] = ode23t(@(t_dummy, X_dummy) PhysicsEngine2(X_dummy, System, U), [t, t_next], X_cur, options);
        X_new = X_step(end, :)'; % Extract the final state of the step
    catch ME
        fprintf('Solver failed at t = %.3f. Reason: %s\n', t, ME.message);
        break;
    end

    % Advance State
    t = t_next;      
    X_cur = X_new;   
    stepCount = stepCount + 1;

    if stepCount > MaxSteps, break; end 
    
    % Log Data
    X_LOG(:, stepCount) = X_cur;
    U_LOG(:, stepCount) = U;
    T_LOG(stepCount) = t; 
    
    if mod(stepCount, 100) == 0
        fprintf('Time: %.2f s | Chamber P: %.2f psi\n', t, X_cur(10)/6895);
    end
end
toc;

%% Trim logs to actual size after loop
X_LOG = X_LOG(:, 1:stepCount-1);
U_LOG = U_LOG(:, 1:stepCount-1);
T_LOG = T_LOG(1:stepCount-1);

% VISUALIZATION (DARK MODE)
close all;

% Convert to common engineering units
Time = T_LOG;
PSI = 6895;
Pressures_PSI = X_LOG(1:PLength, :) / PSI;

% Post-Process Temperatures from Enthalpy
T_Nodes_LOG = zeros(PLength, length(Time));

% Fetch base coefficients (Pressure doesn't affect A and B in our data functions)
[~, ~, ~, A_ox, B_ox] = LOX_Data(293, 101325);
[~, ~, ~, A_fu, B_fu] = IPA_Data(293, 101325);
[~, ~, ~, A_n2, B_n2] = N2_Data(293, 101325);

% StartIdx is already PLength + MDLength (Index right before Y starts)
StartIdx = PLength + MDLength;

for i = 1:PLength
    % Extract Enthalpy history for Node i
    H_history = X_LOG(HIdx_Start + i - 1, :);

    % Extract Mass Fractions for Node i
    Idx_OX = StartIdx + (i-1)*NumSpecies + 1;
    Idx_FU = StartIdx + (i-1)*NumSpecies + 2;
    Idx_N2 = StartIdx + (i-1)*NumSpecies + 3;

    Y_OX_hist = X_LOG(Idx_OX, :);
    Y_FU_hist = X_LOG(Idx_FU, :);
    Y_N2_hist = X_LOG(Idx_N2, :);

    % Calculate Mixture Coefficients over time
    A_mix = Y_OX_hist.*A_ox + Y_FU_hist.*A_fu + Y_N2_hist.*A_n2;
    B_mix = Y_OX_hist.*B_ox + Y_FU_hist.*B_fu + Y_N2_hist.*B_n2;

    % Solve Quadratic for Temperature
    T_Nodes_LOG(i, :) = (-A_mix + sqrt(A_mix.^2 + 2 .* B_mix .* H_history)) ./ (B_mix + 1e-9);
end

% Dynamic Massflows (Pipe Links)
Pipe_Flows = X_LOG(PLength+1 : PLength+MDLength, :);

% Reconstruct Algebraic Massflows (Injectors/Valves) for plotting
Inj_Flows = zeros(2, length(Time)); % 1: OX Inj, 2: FU Inj

% Dynamically fetch parameters from the System Struct
% Algebraic Links: [1:PressIso, 2:MainOX, 3:MainFU, 4:ThrotOX, 5:ThrotFU, 6:InjOX, 7:InjFU...]
L_OX = System.Links.Algebraic(6); % Alg Link 6 (Inj OX)
L_FU = System.Links.Algebraic(7); % Alg Link 7 (Inj FU)

for k = 1:length(Time)
    % Extract State (UPDATED NODE MAP: 11=OX Man, 12=FU Man, 13=SKIPPER Chamber)
    P_Man_OX = X_LOG(11, k); 
    P_Man_FU = X_LOG(12, k);
    P_Cham   = X_LOG(13, k); 
    
    % Fetch N2 Mass Fractions to determine if Manifold is Liquid or Gas
    StartIdx = PLength + MDLength;
    Y_N2_OX = X_LOG(StartIdx + (11-1)*NumSpecies + 3, k); 
    Y_N2_FU = X_LOG(StartIdx + (12-1)*NumSpecies + 3, k); 
    
    % Dynamically calculate density (Transitions from Liquid to Gas)
    Rho_Gas_OX = P_Man_OX / (296.8 * 293); % Ideal Gas Law for N2
    Rho_Gas_FU = P_Man_FU / (296.8 * 293);
    
    Rho_Eff_OX = (1 - Y_N2_OX) * 1141 + (Y_N2_OX) * Rho_Gas_OX;
    Rho_Eff_FU = (1 - Y_N2_FU) * 786  + (Y_N2_FU) * Rho_Gas_FU;
    
    % OX Injector
    DP_OX = P_Man_OX - P_Cham;
    Inj_Flows(1, k) = L_OX.Cv * L_OX.A * sqrt(2 * Rho_Eff_OX * abs(DP_OX)) * sign(DP_OX);
    
    % FU Injector
    DP_FU = P_Man_FU - P_Cham;
    Inj_Flows(2, k) = L_FU.Cv * L_FU.A * sqrt(2 * Rho_Eff_FU * abs(DP_FU)) * sign(DP_FU);
end

% Dark Mode Settings
DarkBg = [0.15 0.15 0.15];
AxColor = [0.9 0.9 0.9];
GridAlpha = 0.2;

% FIG 1: Tank Blowdown (Supply Side)
fig1 = figure('Name', 'Supply Pressures', 'NumberTitle', 'off', 'Color', DarkBg);
    plot(Time, Pressures_PSI(2,:), 'w--', 'LineWidth', 1.5); hold on; % N2 Source (Node 1)
    plot(Time, Pressures_PSI(3,:), 'c', 'LineWidth', 1.5);   % OX Tank (Node 3)
    plot(Time, Pressures_PSI(4,:), 'm', 'LineWidth', 1.5);   % FU Tank (Node 4)
    grid on;
    legend('N2 Source', 'LOX Tank', 'IPA Tank', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    xlabel('Time [s]'); ylabel('Pressure [psi]');
    title('Tank Blowdown Curves');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 2: Feed System Drops (Manifold Dynamics)
% Trace: Tank -> Post-Throttle Node -> Manifold
fig2 = figure('Name', 'Feed System Gradients', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(3,:), 'c--'); hold on;          % Tank OX (Node 3)
    plot(Time, Pressures_PSI(9,:), 'c');                     % Post-Throttle OX (Node 9)
    plot(Time, Pressures_PSI(11,:), 'g', 'LineWidth', 1.5);  % Manifold OX (Node 11)
    grid on;
    legend('Tank', 'Post-Throttle', 'Manifold', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Oxidizer Feed Line Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);
    
subplot(2,1,2);
    plot(Time, Pressures_PSI(4,:), 'm--'); hold on;          % Tank FU (Node 4)
    plot(Time, Pressures_PSI(10,:), 'm');                    % Post-Throttle FU (Node 10)
    plot(Time, Pressures_PSI(12,:), 'y', 'LineWidth', 1.5);  % Manifold FU (Node 12)
    grid on;
    legend('Tank', 'Post-Throttle', 'Manifold', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Fuel Feed Line Pressures');
    ylabel('Pressure [psi]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 3: Engine Performance (Chamber & Injection)
fig3 = figure('Name', 'Engine Operation', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(11,:), 'g'); hold on;           % OX Man (Node 11)
    plot(Time, Pressures_PSI(12,:), 'y');                    % FU Man (Node 12)
    plot(Time, Pressures_PSI(13,:), 'w', 'LineWidth', 2);    % SKIPPER (Node 13)
    yline(14.7, 'w--');
    grid on;
    legend('OX Manifold', 'FU Manifold', 'Chamber', 'Atm', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Injection Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    % Compare Pipe Massflow (Link 5/6 - Manifold In) vs Injector Massflow (Demand)
    % Pipe Flow Indices: 1/2=Press Lines, 3/4=Pre-Mains, 5=OX Inj Line, 6=FU Inj Line
    plot(Time, Pipe_Flows(5,:), 'c'); hold on;               % OX Manifold In
    plot(Time, Inj_Flows(1,:), 'g--');                       % OX Inj
    plot(Time, Pipe_Flows(6,:), 'm');                        % FU Manifold In
    plot(Time, Inj_Flows(2,:), 'y--');                       % FU Inj
    grid on;
    legend('OX Feed (Manifold In)', 'OX Inj', 'FU Feed (Manifold In)', 'FU Inj', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Massflow Dynamics');
    ylabel('Flow [kg/s]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 4: Node Compositions
fig4 = figure('Name', 'Combustion Components', 'NumberTitle', 'off', 'Color', DarkBg);
StartIdx = PLength + MDLength;
NumSpecies = 3; 

subplot(2,2,1);
    % Node 11 (OX Manifold)
    Idx_LOX_Liq = StartIdx + (11-1)*NumSpecies + 1; 
    Idx_LOX_N2  = StartIdx + (11-1)*NumSpecies + 3;
    plot(Time, X_LOG(Idx_LOX_Liq, :), 'c', 'LineWidth', 1.5); hold on;
    plot(Time, X_LOG(Idx_LOX_N2, :), 'w--');
    grid on;
    ylabel('Mass Fraction');
    legend('LOX Liquid', 'Gases', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Oxidizer Manifold Composition');
    ylim([-0.05 1.05]); 
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,2,3);
    % Node 12 (FU Manifold)
    Idx_FU_Liq = StartIdx + (12-1)*NumSpecies + 2; 
    Idx_FU_N2  = StartIdx + (12-1)*NumSpecies + 3;
    plot(Time, X_LOG(Idx_FU_Liq, :), 'm', 'LineWidth', 1.5); hold on;
    plot(Time, X_LOG(Idx_FU_N2, :), 'w--');
    grid on;
    ylabel('Mass Fraction');
    xlabel('Time [s]');
    legend('IPA Liquid', 'Gases', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Fuel Manifold Composition');
    ylim([-0.05 1.05]);
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,2,2);
    % Node 13 (Chamber)
    Idx_C_OX = StartIdx + (13-1)*NumSpecies + 1; 
    Idx_C_FU = StartIdx + (13-1)*NumSpecies + 2; 
    Idx_C_N2 = StartIdx + (13-1)*NumSpecies + 3;
    plot(Time, X_LOG(Idx_C_OX, :), 'c', 'LineWidth', 1.5); hold on;
    plot(Time, X_LOG(Idx_C_FU, :), 'm', 'LineWidth', 1.5);
    plot(Time, X_LOG(Idx_C_N2, :), 'w--');
    grid on;
    ylabel('Mass Fraction');
    xlabel('Time [s]');
    legend('Oxidizer', 'Fuel', 'Gases', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('SKIPPER Chamber Composition');
    ylim([-0.05 1.05]);
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,2,4);
    % OF Ratio Trace
    OF_C = X_LOG(Idx_C_OX, :) ./ (X_LOG(Idx_C_FU, :) + 1e-9);
    plot(Time, OF_C, 'Color', [0.7 0 1], 'LineWidth', 1.5);
    grid on;
    ylabel('O/F Ratio');
    xlabel('Time [s]');
    legend('Chamber OF Ratio', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('SKIPPER Chamber O/F Ratio');
    ylim([0.5 1.9]);
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

   % FIG 6: COPV Dynamics (Pressure & Temperature)
fig6 = figure('Name', 'COPV Dynamics', 'NumberTitle', 'off', 'Color', DarkBg);

subplot(2,1,1);
    plot(Time, Pressures_PSI(1,:), 'w', 'LineWidth', 1.5); 
    grid on;
    legend('COPV Pressure', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('COPV Blowdown Pressure');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    plot(Time, T_Nodes_LOG(1,:), 'r', 'LineWidth', 1.5); 
    grid on;
    legend('COPV Temperature', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('COPV Internal Gas Temperature');
    ylabel('Temperature [K]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 5: AutoSequence
    Scheduler.PlotSequence(round(T_LOG(end)));