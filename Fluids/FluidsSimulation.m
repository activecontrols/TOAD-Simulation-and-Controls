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
simTime = 20;
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
        ValveTiming = [5 70];
        ValveRamp = 0.12;
        U = zeros(MALength, 1);

        % Throttle Valve Open
        if t_next < ValveTiming(1)
            U = zeros(MALength, 1);
        elseif t_next < (ValveTiming(1) + ValveRamp)
            Progress = (t_next - ValveTiming(1)) / ValveRamp;
            U(3:4) = 2.9 * Progress;
        else
            U(3:4) = 2.9;
        end

        % COPV Valve Close
        if t_next < ValveTiming(2)
            U(1:2) = 2.9;
        elseif t_next < (ValveTiming(2) + ValveRamp)
            Progress = (t_next - ValveTiming(1)) / ValveRamp;
            U(1:2) = 0;
        else
            U(1:2) = 0;
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
        fprintf('Time: %.2f s | Tank P: %.2f psi | dT: %.1e\n', ...
            t, X_cur(2)/6895, dT);
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
% Note: Adjust indices (7 and 8) if your Link IDs change
L_OX = System.Links.Algebraic(5); 
L_FU = System.Links.Algebraic(6);

for k = 1:length(Time)
    % Extract State
    P_Man_OX = X_LOG(6, k); P_Cham = X_LOG(8, k); P_Man_FU = X_LOG(7, k);
    
    % OX Injector (Link 7)
    DP_OX = P_Man_OX - P_Cham;
    % Uses L_OX.Cv (which is Cd) and L_OX.Area automatically
    Inj_Flows(1, k) = L_OX.Cv * L_OX.A * sqrt(2 * 1141 * abs(DP_OX)) * sign(DP_OX);
    
    % FU Injector (Link 8)
    DP_FU = P_Man_FU - P_Cham;
    Inj_Flows(2, k) = L_FU.Cv * L_FU.A * sqrt(2 * 786 * abs(DP_FU)) * sign(DP_FU);
end

% Dark Mode Settings
DarkBg = [0.15 0.15 0.15]; % Dark Grey Background
AxColor = [0.9 0.9 0.9];   % Off-White Axes/Text
GridAlpha = 0.2;

% FIG 1: Tank Blowdown (Supply Side)
fig1 = figure('Name', 'Supply Pressures', 'NumberTitle', 'off', 'Color', DarkBg);
    plot(Time, Pressures_PSI(1,:), 'w--', 'LineWidth', 1.5); hold on; % Regulator
    plot(Time, Pressures_PSI(2,:), 'c', 'LineWidth', 1.5);   % OX (Cyan)
    plot(Time, Pressures_PSI(3,:), 'm', 'LineWidth', 1.5);   % FU (Magenta)
    grid on;
    legend('Regulator', 'LOX Tank', 'IPA Tank', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    xlabel('Time [s]'); ylabel('Pressure [psi]');
    title('Tank Blowdown Curves');
    
    % Apply Dark Theme
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 2: Feed System Drops (Manifold Dynamics)
fig2 = figure('Name', 'Feed System Gradients', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(2,:), 'c--'); hold on;          % Tank (Cyan Dash)
    plot(Time, Pressures_PSI(4,:), 'c');                     % Line (Cyan)
    plot(Time, Pressures_PSI(6,:), 'g', 'LineWidth', 1.5);   % Manifold (Green)
    grid on;
    legend('Tank', 'Line (Post-Valve)', 'Manifold (Pre-Inj)', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Oxidizer Feed Line Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);
    
subplot(2,1,2);
    plot(Time, Pressures_PSI(3,:), 'm--'); hold on;          % Tank (Magenta Dash)
    plot(Time, Pressures_PSI(5,:), 'm');                     % Line (Magenta)
    plot(Time, Pressures_PSI(7,:), 'y', 'LineWidth', 1.5);   % Manifold (Yellow)
    grid on;
    legend('Tank', 'Line (Post-Valve)', 'Manifold (Pre-Inj)', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Fuel Feed Line Pressures');
    ylabel('Pressure [psi]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 3: Engine Performance (Chamber & Injection)
fig3 = figure('Name', 'Engine Operation', 'NumberTitle', 'off', 'Color', DarkBg);
subplot(2,1,1);
    plot(Time, Pressures_PSI(6,:), 'g'); hold on;            % OX Man (Green)
    plot(Time, Pressures_PSI(7,:), 'y');                     % FU Man (Yellow)
    plot(Time, Pressures_PSI(8,:), 'w', 'LineWidth', 2);     % Chamber (White)
    yline(14.7, 'w--');                                      % Atm (White Dash)
    grid on;
    legend('OX Manifold', 'FU Manifold', 'Chamber', 'Atm', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Injection Pressures');
    ylabel('Pressure [psi]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    % Compare Pipe Massflow (Supply) vs Injector Massflow (Demand)
    plot(Time, Pipe_Flows(1,:), 'c'); hold on;               % OX Pipe (Cyan)
    plot(Time, Inj_Flows(1,:), 'g--');                       % OX Inj (Green Dash)
    plot(Time, Pipe_Flows(2,:), 'm');                        % FU Pipe (Magenta)
    plot(Time, Inj_Flows(2,:), 'y--');                       % FU Inj (Yellow Dash)
    grid on;
    legend('OX Pipe', 'OX Inj', 'FU Pipe', 'FU Inj', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Massflow Dynamics (Lag Check)');
    ylabel('Flow [kg/s]');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

% FIG 4: Fluid Composition (Purge/Mixing Check)
fig4 = figure('Name', 'Fluid Composition', 'NumberTitle', 'off', 'Color', DarkBg);
StartIdx = PLength + MDLength;

subplot(2,1,1);
    % Plot OX Manifold Gas Fraction
    Idx_ManOX = StartIdx + (6-1)*3 + 3; % N2 index
    plot(Time, X_LOG(Idx_ManOX, :), 'c'); hold on;
    Idx_ManFU = StartIdx + (7-1)*3 + 3; 
    plot(Time, X_LOG(Idx_ManFU, :), 'm');
    grid on;
    ylabel('N2 Mass Fraction');
    legend('OX Manifold', 'FU Manifold', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Quality Check (1.0 = Pure Gas, 0.0 = Pure Liquid)');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);

subplot(2,1,2);
    % Chamber Mixture Ratio
    Idx_ChamOX = StartIdx + (8-1)*3 + 1;
    Idx_ChamFU = StartIdx + (8-1)*3 + 2;
    plot(Time, X_LOG(Idx_ChamOX, :), 'c'); hold on;
    plot(Time, X_LOG(Idx_ChamFU, :), 'm');
    grid on;
    legend('Oxidizer', 'Fuel', 'TextColor', 'w', 'Color', 'none', 'EdgeColor', 'none');
    title('Chamber Composition');
    ylabel('Mass Fraction');
    xlabel('Time [s]');
    set(gca, 'Color', DarkBg, 'XColor', AxColor, 'YColor', AxColor, 'GridColor', 'w', 'GridAlpha', GridAlpha);