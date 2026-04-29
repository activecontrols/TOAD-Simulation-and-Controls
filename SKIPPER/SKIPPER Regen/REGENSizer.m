%% Code to automatically size Regen channel dimensions to optimize for engine life
% Pablo Plata
% Code operates via Bayesian Optimization on channel dimensions and with
% objective of engine life. Current design shall operate at max throttle
% only, eventually goal will be to increase median engine lifespan
% including min throttle aswell.
% 
% INPUTS:   WallThickness -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelHeight -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelWidth  -> 1x2 array, [Chamber, Nozzle]
%
% Design vector for solver:
% X = [WallThickness, ChannelHeight, ChannelWidth];

clear;
clc;
close all;

%% Parameter sampling (inches)
addpath('BayesianOpt\');
Thickness = [0.05; 0.1];    Height = [0.01; 0.125];  Width = [0.01; 0.125];
LowerBounds = [ones(1, 3) * Thickness(1), ones(1, 3) * Height(1), ones(1, 2) * Width(1), 40];
UpperBounds = [ones(1, 3) * Thickness(2), ones(1, 3) * Height(2), ones(1, 2) * Width(2), 85];
MaxDP = 150; % psi

% Latin Hypercube Sampling for parameter space for GP training
NumDims = length(LowerBounds);
NumSamples = 150;
LHS = HyperSampl(NumSamples, NumDims);
Geometries = LowerBounds + LHS .* (UpperBounds - LowerBounds);

%% Feed sample space through SKIPPERRegen.m to evaluate physics
Lifespan = zeros(NumSamples, 1);
PressDrop = zeros(NumSamples, 1);
Invalid = 0;
for i = 1:NumSamples
    NC = Geometries(i, 9);
    WT = Geometries(i, 1:3);
    CH = Geometries(i, 4:6);
    CW = Geometries(i, 7:8);

    % Evaluate regen
    try
        [Lifespan(i), PressDrop(i)] = SKIPPERRegen(NC, WT, CH, CW, 0);
    catch
        % Invalid geometry
        Lifespan(i) = NaN;
        PressDrop(i) = NaN;
        Invalid = Invalid + 1;
        warning('Sample Number %i Produced invalid geometry', i);
    end
    if ~isnan(Lifespan(i))
        fprintf('Evaluation #%i Complete!, %.2f Cycles & %.2f psi drop\n', i, Lifespan(i), PressDrop(i));
    end
end
fprintf('Initial Evaluation complete!, %.2f%% of sampled geometries were invalid\n', Invalid / NumSamples * 100);

%% Gaussian Process Training 
% Output scaling
LifespanDEV = std(Lifespan, 'omitnan'); LifespanMEAN = mean(Lifespan, 'omitnan');
PressDEV = std(PressDrop, 'omitnan');   PressMEAN = mean(PressDrop, 'omitnan');
LifespanSCALED = (Lifespan - LifespanMEAN) / LifespanDEV;
PressDropSCALED = (PressDrop - PressMEAN) / PressDEV;

% NaN handling
WorstLifespan = min(LifespanSCALED);
WorstPressDrop = max(PressDropSCALED);
PenaltyLife = WorstLifespan - 0.5;
PenaltyPress = WorstPressDrop + 0.5;
LifespanSCALED(isnan(LifespanSCALED)) = PenaltyLife;
PressDropSCALED(isnan(PressDropSCALED)) = PenaltyPress;

% Parameter initialization for Kernel
% Define the kernel parameters for Gaussian Process
lengthScales = zeros(NumDims, 1);
signalVar = 0;
noiseVar = log(0.01);
logThetaInit = [lengthScales; signalVar; noiseVar];

% Set up optimization for parameters
options = optimset('Display', 'off', 'MaxFunEvals', 1000);
logThetaOpt_LIFE = fminsearch(@(t) NLML(t, Geometries, LifespanSCALED), logThetaInit, options);
logThetaOpt_PRESS = fminsearch(@(t) NLML(t, Geometries, PressDropSCALED), logThetaInit, options);
ThetaOpt_LIFE = exp(logThetaOpt_LIFE);
ThetaOpt_PRESS = exp(logThetaOpt_PRESS);

%% Optimizer
NumIterations = 80;
for iter = 1:NumIterations
    % Loop parameters
    MaxDP_Scaled = (MaxDP - PressMEAN) / PressDEV;
    y_train = [LifespanSCALED, PressDropSCALED];
    
    % Find indices that didn't crash AND are under 150 psi
    FeasibleIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP);
    if any(FeasibleIdx)
        BestValidScaled = max(LifespanSCALED(FeasibleIdx));
    else
        BestValidScaled = min(LifespanSCALED); 
    end

    % Global rough search using grid search
    % Precompute L matrices 
    NumTrain = size(Geometries, 1);
    
    % Life GP
    len_LIFE = ThetaOpt_LIFE(1:end-2)'; sig_LIFE = ThetaOpt_LIFE(end-1); noise_LIFE = ThetaOpt_LIFE(end);
    K_LIFE = MaternKernel(Geometries, Geometries, len_LIFE, sig_LIFE) + noise_LIFE * eye(NumTrain);
    L_LIFE = chol(K_LIFE, 'lower');
    alpha_LIFE = L_LIFE' \ (L_LIFE \ LifespanSCALED);
    
    % Pressure GP
    len_PRESS = ThetaOpt_PRESS(1:end-2)'; sig_PRESS = ThetaOpt_PRESS(end-1); noise_PRESS = ThetaOpt_PRESS(end);
    K_PRESS = MaternKernel(Geometries, Geometries, len_PRESS, sig_PRESS) + noise_PRESS * eye(NumTrain);
    L_PRESS = chol(K_PRESS, 'lower');
    alpha_PRESS = L_PRESS' \ (L_PRESS \ PressDropSCALED);

    % Run rough search
    NumGrid = 10000;
    GridGeometries = LowerBounds + rand(NumGrid, NumDims) .* (UpperBounds - LowerBounds);
    GridScores = zeros(NumGrid, 1);
    for g = 1:NumGrid
        GridScores(g) = Acquisition(GridGeometries(g,:), Geometries, LowerBounds, UpperBounds, L_LIFE, alpha_LIFE, ThetaOpt_LIFE, L_PRESS, alpha_PRESS, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled);
    end
    % CurrentMaxEI = -min(GridScores);
    % if iter > 10 && CurrentMaxEI < 1e-5
    %     fprintf('\nConvergence detected.\n');
    %     break;
    % end

    % Refine the search using fminsearch
    [~, bestIdx] = min(GridScores);
    x0 = GridGeometries(bestIdx, :);
    AcqObj = @(x) Acquisition(x, Geometries, LowerBounds, UpperBounds, L_LIFE, alpha_LIFE, ThetaOpt_LIFE, L_PRESS, alpha_PRESS, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled);
    x_next = fminsearch(AcqObj, x0, options);

    % Evaluate the physics
    WT = x_next(1:3); CH = x_next(4:6); CW = x_next(7:8); NC = x_next(9);
    fprintf('Testing new geometry (#%i out of %i): \n NC = %i, WT=[%.3f, %.3f, %.3f], CH=[%.3f, %.3f, %.3f], CW=[%.3f, %.3f]\n', iter, NumIterations, round(NC), x_next(1:8));
    try
        [Life_new, Drop_new] = SKIPPERRegen(NC, WT, CH, CW, 0);
    catch
        Life_new = NaN; Drop_new = NaN;
        warning('Geometry failed in physics solver.');
    end
    
    if ~isnan(Life_new)
        fprintf('Result: %.2f Cycles, %.2f psi drop\n', Life_new, Drop_new);
    end

    % Append data to arrays
    Geometries = [Geometries; x_next];
    Lifespan = [Lifespan; Life_new];
    PressDrop = [PressDrop; Drop_new];

    % Repenalize and rescale
    LifespanDEV = std(Lifespan, 'omitnan'); LifespanMEAN = mean(Lifespan, 'omitnan');
    PressDEV = std(PressDrop, 'omitnan');   PressMEAN = mean(PressDrop, 'omitnan');
    LifespanSCALED = (Lifespan - LifespanMEAN) / LifespanDEV;
    PressDropSCALED = (PressDrop - PressMEAN) / PressDEV;
    
    % NaN handling
    WorstLifespan = min(LifespanSCALED);
    WorstPressDrop = max(PressDropSCALED);
    PenaltyLife = WorstLifespan - 0.5;
    PenaltyPress = WorstPressDrop + 0.5;
    LifespanSCALED(isnan(LifespanSCALED)) = PenaltyLife;
    PressDropSCALED(isnan(PressDropSCALED)) = PenaltyPress;

    % Retrain GPs using the optimal Thetas as the guess
    logThetaOpt_LIFE = fminsearch(@(t) NLML(t, Geometries, LifespanSCALED), logThetaOpt_LIFE, options);
    logThetaOpt_PRESS = fminsearch(@(t) NLML(t, Geometries, PressDropSCALED), logThetaOpt_PRESS, options);
    ThetaOpt_LIFE = exp(logThetaOpt_LIFE);
    ThetaOpt_PRESS = exp(logThetaOpt_PRESS);
end

%% Plotting
% Filter for Valid & Feasible Results
% Must not be NaN (physics failed) AND must be under the MaxDP limit
ValidIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP);

if ~any(ValidIdx)
    warning('No valid geometries were found that satisfied the pressure constraint!');
else
    % Extract only valid data
    ValidLifespans = Lifespan(ValidIdx);
    ValidPressDrops = PressDrop(ValidIdx);
    ValidGeometries = Geometries(ValidIdx, :);
    
    % Find the absolute best lifespan among the feasible points
    [BestLife, relative_idx] = max(ValidLifespans);
    BestPress = ValidPressDrops(relative_idx);
    ChampionGeom = ValidGeometries(relative_idx, :);
    
    % Console Printout of the Champion
    fprintf('CHAMPION PERFORMANCE:\n');
    fprintf('  Lifespan:      %.2f Cycles\n', BestLife);
    fprintf('  Pressure Drop: %.2f psi\n\n', BestPress);
    
    fprintf('OPTIMAL GEOMETRY [inches]:\n');
    fprintf('  Wall Thickness (Chamber, Throat, Nozzle): [%.4f, %.4f, %.4f]\n', ChampionGeom(1:3));
    fprintf('  Channel Height (Chamber, Throat, Nozzle): [%.4f, %.4f, %.4f]\n', ChampionGeom(4:6));
    fprintf('  Channel Width  (Chamber, Nozzle):         [%.4f, %.4f]\n', ChampionGeom(7:8));
    fprintf('  Channel Count:                             %i \n', ChampionGeom(9));
    fprintf('======================================================\n');

    % Calculate Convergence Tracking
    TotalEvals = length(Lifespan);
    BestSoFar = zeros(TotalEvals, 1);
    CurrentMax = 0;
    
    for i = 1:TotalEvals
        % Update current max only if it's valid, feasible, and better
        if ~isnan(Lifespan(i)) && (PressDrop(i) <= MaxDP) && (Lifespan(i) > CurrentMax)
            CurrentMax = Lifespan(i);
        end
        % If no valid points found yet, keep it at 0
        BestSoFar(i) = CurrentMax;
    end

    % Convergence & History Plots
    figure('Name', 'Bayesian Optimization Convergence', 'Position', [100, 100, 800, 600]);
    
    % Lifespan Convergence
    subplot(2, 1, 1);
    hold on; set(gca, 'FontName', 'Times New Roman');
    scatter(1:TotalEvals, Lifespan, 20, 'k', 'filled', 'MarkerFaceAlpha', 0.3);
    plot(1:TotalEvals, BestSoFar, 'b-', 'LineWidth', 2.5);
    xline(NumSamples, 'k--', 'End of LHS Init', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
    
    title('Objective Convergence: Engine Lifespan');
    ylabel('Lifespan [Cycles]');
    legend('Sampled Geometries', 'Optimal Valid Lifespan', 'Location', 'northwest');
    grid on;
    ylim([0, max(Lifespan) * 1.1]);
    
    % Pressure Drop Tracking
    subplot(2, 1, 2);
    hold on; set(gca, 'FontName', 'Times New Roman');
    scatter(1:TotalEvals, PressDrop, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
    yline(MaxDP, 'r-', 'Max Constraint (150 psi)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
    xline(NumSamples, 'k--', 'LineWidth', 1.5);
    
    title('Constraint Tracking: Coolant Pressure Drop');
    xlabel('Evaluation Number (Initialization + Active Search)');
    ylabel('Pressure Drop [psi]');
    grid on;
    % Dynamically set Y-limits so massive 500+ psi failures don't ruin the scale
    ylim([0, max(MaxDP + 50, prctile(PressDrop, 85))]); 

    % Call SKIPPERRegen for Final Native Plots
    
    % SKIPPERRegen(ChampionGeom(9), ChampionGeom(1:3), ChampionGeom(4:6), ChampionGeom(7:8), 1);
end



