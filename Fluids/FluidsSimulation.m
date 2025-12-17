%% This file is test Simulation run of the P&ID Model Simulation.
% Load in PID Structure
addpath('.\Helpers');
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
simTime = 50;
dT = 1e-5;
MaxdT = 1e-2;
MindT = 1e-5;
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

% Generate Symbolic Files
FluidsSolver2([], System, 'Analytic');

% Warm-start the implicit solver
fprintf('Warming up system states...\n');
for k = 1:1000
    XDOT = FluidsSolver2(X_cur, System, 'Numerical');
    X_cur = X_cur + XDOT * 1e-5;
end

fprintf('Warm up complete. Starting Implicit Loop.\n');
tic;
while t < simTime
    X_old = X_cur;
    isAccepted = false;
    while ~isAccepted
        % Use proposed step
        t_next = t + dT;

        %% Valve Control Logic
        ValveTiming = 10;
        ValveRamp = 0.12;
        U = zeros(MALength, 1);
        if t_next < ValveTiming
            U = zeros(MALength, 1);
        elseif t_next < (ValveTiming + ValveRamp)
            Progress = (t_next - ValveTiming) / ValveRamp;
            U(1) = 0.5 * Progress;
        else
            U(1) = 0.5;
        end
    
        %% Implicit Solver Step via Newton-Ralphson
        X_new = X_old;
        isConverged = false;

        for iter = 1:20
            % Evaluate Dynamics and Jacobian
            F = FluidDynamics(X_new, U);
            J = FluidJacobian(X_new, U);
    
            % Residual
            R = X_new - X_old - dT * F;
    
            % Convergence check
            if max(abs(R)) < 1e-5
                isConverged = true;
                break;
            end
    
            % Newton Update: Delta = - (I - dt*J)^-1 * R
            M = eye(length(X0)) - dT * J;
            Delta = -M \ R;
            
            X_new = X_new + Delta;
        end

        % Step Acceptance
        if isConverged
            % Metric: Max Relative Change
            RelChange = abs(X_new - X_old) ./ (abs(X_new) + 1e-2);
            MaxChange = max(RelChange);
            Tolerance = 0.15;  
            Target    = 0.04;  
            
            if MaxChange > Tolerance && dT > MindT
                isAccepted = false;
                dT = max(dT * 0.5, MindT); 
            else
                isAccepted = true;

                % ~PI Control for timestep
                Error = (Target / (MaxChange + 1e-9));
                Factor = (Error)^0.4;
                Beta = 0.02;
                dT = (1 - Beta) * dT + Beta * Factor * dT;
                
                % Solver Difficulty Check (Override)
                if iter >= 10
                    dT = min(dT, dT * 0.9);
                end
            end
        else
            % Diverged (Newton failed)
            if dT <= MindT
                isAccepted = true; 
            else
                dT = max(dT * 0.5, MindT);
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
    if mod(stepCount, 500) == 0
        fprintf('Time: %.2f s | Tank P: %.2f psi | dT: %.1e\n', ...
            t, X_cur(1)/6895, dT);
    end
end

% Trim logs to actual size after loop
X_LOG = X_LOG(:, 1:stepCount);
U_LOG = U_LOG(:, 1:stepCount);
T_LOG = T_LOG(1:stepCount);
toc;




