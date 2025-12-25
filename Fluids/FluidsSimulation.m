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

X0 = zeros(PLength + MDLength + PLength * NumSpecies, 1);
% Package Pressures and Transports
for i = 1:PLength
    X0(i) = System.Nodes(i).P0 * 6895;
    CurrentIDX = YIdx_Start + (i-1)*NumSpecies : YIdx_Start + (i-1)*NumSpecies + (NumSpecies-1);
    X0(CurrentIDX) = System.Nodes(i).Y0';
end
% Package Dynamic Massflows 
for i = 1:MDLength
    X0(PLength + i) = 0;
end

%% Simulation Loop 2 (Implicit Euler w/ Variable timestep)
% Initialize Solver 2 parameters
simTime = 14;
dT = 1e-3;
MaxdT = 5e-2;
MindT = 5e-5;
MaxSteps = 1e5;

% Logs
X_LOG = zeros(size(X0, 1), MaxSteps);
U_LOG = zeros(MALength, MaxSteps);
T_LOG = zeros(1, MaxSteps);
X_cur = X0;
t = 0;
T_LOG(1) = t;
X_LOG(:, 1) = X_cur;
stepCount = 1;
lastError = 0;

% Warm-start the implicit solver
fprintf('Warming up system states...\n');
for k = 1:10
    XDOT = PhysicsEngine(X_cur, System, zeros(MALength,1), true);
    X_cur = X_cur + XDOT * 5e-5;
end

fprintf('Warm up complete. Starting Implicit Loop.\n');
tic;
warning('off', 'MATLAB:nearlySingularMatrix');
while t < simTime
    X_old = X_cur;
    isAccepted = false;
    while ~isAccepted
        % Use proposed step
        t_next = t + dT;

%% Valve Control Logic
        ValveTiming = [5 9];
        RunValves = [2.7, 1.05]; % [Cv_OX, Cv_FU]
        ValveRamp = 0.03;        
        OxDelay   = 0.06;        
        
        U = zeros(MALength, 1);
        U(1:2) = 2;              % Keep Pressurization Lines Open
        
        % FU Valve
        if t_next < ValveTiming(1)
            U(4) = 0;
        elseif t_next < (ValveTiming(1) + ValveRamp)
            Progress = (t_next - ValveTiming(1)) / ValveRamp;
            U(4) = RunValves(2) * Progress;
        elseif t_next < ValveTiming(2)
            U(4) = RunValves(2);
        elseif t_next < (ValveTiming(2) + ValveRamp)
            Progress = (t_next - ValveTiming(2)) / ValveRamp;
            U(4) = RunValves(2) * (1 - Progress);
        else
            U(4) = 0;
        end

        % OX Valve
        OxStart = ValveTiming(1) + OxDelay;
        
        if t_next < OxStart
            U(3) = 0;
        elseif t_next < (OxStart + ValveRamp)
            Progress = (t_next - OxStart) / ValveRamp;
            U(3) = RunValves(1) * Progress;
        elseif t_next < ValveTiming(2)
            U(3) = RunValves(1);
        elseif t_next < (ValveTiming(2) + ValveRamp)
            Progress = (t_next - ValveTiming(2)) / ValveRamp;
            U(3) = RunValves(1) * (1 - Progress);
        else
            U(3) = 0;
        end
    
        %% Implicit Solver Step via Newton-Ralphson
        X_new = X_old;
        isConverged = false;

        for iter = 1:35
            % Evaluate Dynamics and Jacobian
            F = PhysicsEngine(X_new, System, U, true);
            J = NumericalJacobian(@(x) PhysicsEngine(x, System, U, true), X_new);
            M = eye(length(X0)) - dT * J;
            [LD, UD] = lu(M);
    
            % Residual
            R = X_new - X_old - dT * F;
    
            % Convergence check
            if max(abs(R)) < 1e-5
                isConverged = true;
                break;
            end
    
            % Newton Update: Delta = - (I - dt*J)^-1 * R
            Delta = - (UD \ (LD \ R));
            X_new = X_new + Delta;
        end

        %% Step Acceptance for Variable Timestep
        if isConverged
            Y_Start_Idx = PLength + MDLength + 1;
            X_new(Y_Start_Idx:end) = max(min(X_new(Y_Start_Idx:end), 1.0), 0.0);
            % Metric: Max Relative Change
            RelChange = abs(X_new - X_old) ./ (abs(X_new) + 1e-2);
            MaxChange = max(RelChange);
            Tolerance = 0.15;  
            Target    = 0.05;  
            
            if MaxChange > Tolerance && dT > MindT
                isAccepted = false;
                dT = max(dT * 0.5, MindT); 
            else
                isAccepted = true;

                % Lowpass for timestep
                Error = (Target / (MaxChange + 1e-9));
                Factor = (Error)^0.5;
                Beta = 0.03;
                dT = (1 - Beta) * dT + Beta * Factor * dT;
                
                % Solver Difficulty Check (Override)
                if iter >= 10
                    dT = min(dT, dT * 0.4);
                end
            end
        else
            % Diverged (Newton failed)
            if dT <= MindT
                isAccepted = true; 
            else
                dT = max(dT * 0.4, MindT);
                WasRejected = true;
            end
        end
        dT = min(max(dT, MindT), MaxdT);
    end
    % Advance State
    t = t_next;      % Update time
    X_cur = X_new;   % Update state
    stepCount = stepCount + 1;
    
    % Buffer Safety
    if stepCount > MaxSteps, break; end 
    
    % Log Data
    X_LOG(:, stepCount) = X_cur;
    U_LOG(:, stepCount) = U;
    T_LOG(stepCount) = t; 
    
    % Progress Print
    if mod(stepCount, 100) == 0
        fprintf('Time: %.2f s | Chamber P: %.2f psi | dT: %.1e\n', ...
            t, X_cur(10)/6895, dT);
    end
end
warning('on', 'MATLAB:nearlySingularMatrix');

% Trim logs to actual size after loop
X_LOG = X_LOG(:, 1:stepCount);
U_LOG = U_LOG(:, 1:stepCount);
T_LOG = T_LOG(1:stepCount);
toc;

%% VISUALIZATION (DARK MODE)
close all
% Convert to common engineering units
Time = T_LOG;
PSI = 6895;
Pressures_PSI = X_LOG(1:PLength, :) / PSI;

% Dynamic Massflows (Pipe Links)
Pipe_Flows = X_LOG(PLength+1 : PLength+MDLength, :);

% Reconstruct Algebraic Massflows (Injectors/Valves) for plotting
Inj_Flows = zeros(2, length(Time)); % 1: OX Inj, 2: FU Inj

% Dynamically fetch parameters from the System Struct
% Algebraic Links: [1:PressOX, 2:PressFU, 3:MainOX, 4:MainFU, 5:InjOX, 6:InjFU, 7:Noz]
L_OX = System.Links.Algebraic(5); % Link 9 (Inj OX)
L_FU = System.Links.Algebraic(6); % Link 10 (Inj FU)

for k = 1:length(Time)
    % Extract State (UPDATED NODE MAP: 8=OX Man, 9=FU Man, 10=Chamber)
    P_Man_OX = X_LOG(8, k); 
    P_Cham   = X_LOG(10, k); 
    P_Man_FU = X_LOG(9, k);
    
    % OX Injector
    DP_OX = P_Man_OX - P_Cham;
    Inj_Flows(1, k) = L_OX.Cv * L_OX.A * sqrt(2 * 1141 * abs(DP_OX)) * sign(DP_OX);
    
    % FU Injector
    DP_FU = P_Man_FU - P_Cham;
    Inj_Flows(2, k) = L_FU.Cv * L_FU.A * sqrt(2 * 786 * abs(DP_FU)) * sign(DP_FU);
end

% Dark Mode Settings
DarkBg = [0.15 0.15 0.15];
AxColor = [0.9 0.9 0.9];
GridAlpha = 0.2;

% FIG 1: Tank Blowdown (Supply Side)
fig1 = figure('Name', 'Supply Pressures', 'NumberTitle', 'off', 'Color', DarkBg);
    plot(Time, Pressures_PSI(1,:), 'w--', 'LineWidth', 1.5); hold on; % Regulator (Node 1)
    plot(Time, Pressures_PSI(2,:), 'c', 'LineWidth', 1.5);   % OX Tank (Node 2)
    plot(Time, Pressures_PSI(3,:), 'm', 'LineWidth', 1.5);   % FU Tank (Node 3)
    grid on;
    legend('Regulator', 'LOX Tank', 'IPA Tank', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    xlabel('Time [s]'); ylabel('Pressure [psi]');
    title('Tank Blowdown Curves');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 2: Feed System Drops (Manifold Dynamics)
% Trace: Tank -> Post-Valve Node -> Manifold
fig2 = figure('Name', 'Feed System Gradients', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(2,:), 'c--'); hold on;          % Tank OX (Node 2)
    plot(Time, Pressures_PSI(6,:), 'c');                     % Post-Valve OX (Node 6)
    plot(Time, Pressures_PSI(8,:), 'g', 'LineWidth', 1.5);   % Manifold OX (Node 8)
    grid on;
    legend('Tank', 'Line (Post-Valve)', 'Manifold', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Oxidizer Feed Line Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);
    
subplot(2,1,2);
    plot(Time, Pressures_PSI(3,:), 'm--'); hold on;          % Tank FU (Node 3)
    plot(Time, Pressures_PSI(7,:), 'm');                     % Post-Valve FU (Node 7)
    plot(Time, Pressures_PSI(9,:), 'y', 'LineWidth', 1.5);   % Manifold FU (Node 9)
    grid on;
    legend('Tank', 'Line (Post-Valve)', 'Manifold', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Fuel Feed Line Pressures');
    ylabel('Pressure [psi]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 3: Engine Performance (Chamber & Injection)
fig3 = figure('Name', 'Engine Operation', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(8,:), 'g'); hold on;            % OX Man (Node 8)
    plot(Time, Pressures_PSI(9,:), 'y');                     % FU Man (Node 9)
    plot(Time, Pressures_PSI(10,:), 'w', 'LineWidth', 2);    % Chamber (Node 10)
    yline(14.7, 'w--');
    grid on;
    legend('OX Manifold', 'FU Manifold', 'Chamber', 'Atm', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Injection Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    % Compare Pipe Massflow (Link 7/8 - Manifold In) vs Injector Massflow (Demand)
    % Pipe Flow Indices: 1=Link3, 2=Link4, 3=Link7, 4=Link8
    plot(Time, Pipe_Flows(3,:), 'c'); hold on;               % OX Pipe 2 (Link 7)
    plot(Time, Inj_Flows(1,:), 'g--');                       % OX Inj
    plot(Time, Pipe_Flows(4,:), 'm');                        % FU Pipe 2 (Link 8)
    plot(Time, Inj_Flows(2,:), 'y--');                       % FU Inj
    grid on;
    legend('OX Feed (Manifold In)', 'OX Inj', 'FU Feed (Manifold In)', 'FU Inj', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Massflow Dynamics');
    ylabel('Flow [kg/s]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 4: Tank Composition (Node 2 & 3 - Indices Unchanged)
fig4 = figure('Name', 'Tank Composition', 'NumberTitle', 'off', 'Color', DarkBg);
StartIdx = PLength + MDLength;
NumSpecies = 3; 

subplot(2,1,1);
    % Node 2 (OX Tank)
    Idx_LOX_Liq = StartIdx + (2-1)*NumSpecies + 1; 
    Idx_LOX_N2  = StartIdx + (2-1)*NumSpecies + 3;
    plot(Time, X_LOG(Idx_LOX_Liq, :), 'c', 'LineWidth', 1.5); hold on;
    plot(Time, X_LOG(Idx_LOX_N2, :), 'w--');
    grid on;
    ylabel('Mass Fraction');
    legend('LOX Liquid', 'N2 Gas', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Oxidizer Tank Composition');
    ylim([-0.05 1.05]); 
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    % Node 3 (FU Tank)
    Idx_FU_Liq = StartIdx + (3-1)*NumSpecies + 2; 
    Idx_FU_N2  = StartIdx + (3-1)*NumSpecies + 3;
    plot(Time, X_LOG(Idx_FU_Liq, :), 'm', 'LineWidth', 1.5); hold on;
    plot(Time, X_LOG(Idx_FU_N2, :), 'w--');
    grid on;
    ylabel('Mass Fraction');
    xlabel('Time [s]');
    legend('IPA Liquid', 'N2 Gas', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Fuel Tank Composition');
    ylim([-0.05 1.05]);
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);