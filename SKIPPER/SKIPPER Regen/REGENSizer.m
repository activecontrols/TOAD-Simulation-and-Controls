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
addpath("cea\");
addpath("IPA Data\");
addpath("Material Data\");
addpath("Contours\");
Data = LoadData();

%% Parameter sampling (inches)
addpath('BayesianOpt\');
Thickness = [0.05; 0.1];    AR = [0.2; 3.0];  Width = [0.01; 0.125];
LowerBounds = [ones(1, 3) * Thickness(1), ones(1, 3) * AR(1), ones(1, 2) * Width(1), 40];
UpperBounds = [ones(1, 3) * Thickness(2), ones(1, 3) * AR(2), ones(1, 2) * Width(2), 100];
InputRange  = UpperBounds - LowerBounds;    % FIX 1: used to normalize all GP inputs to [0,1]^D
MaxDP = 150 * 0.5^2; % psi

% Latin Hypercube Sampling for parameter space for GP training
NumDims = length(LowerBounds);
NumSamples = 130;
LHS = HyperSampl(NumSamples, NumDims);
Geometries = LowerBounds + LHS .* (UpperBounds - LowerBounds);

%% Feed sample space through SKIPPERRegen.m to evaluate physics
Lifespan = zeros(NumSamples, 1);
PressDrop = zeros(NumSamples, 1);
Invalid = 0;
for i = 1:NumSamples
    NC = Geometries(i, 9);
    WT = Geometries(i, 1:3);
    AR = Geometries(i, 4:6);
    CW = Geometries(i, 7:8);
    try
        [Lifespan(i), PressDrop(i)] = SKIPPERRegen(Data, NC, WT, AR, CW, 0);
        if ~isreal(Lifespan(i)) || ~isreal(PressDrop(i))
            error('Complex physics output detected.'); 
        end
    catch
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
LifespanDEV  = std(Lifespan,  'omitnan'); 
LifespanMEAN = mean(Lifespan, 'omitnan');
LifespanSCALED = (Lifespan - LifespanMEAN) / LifespanDEV;

% Log-transform the pressure drop
LogPressDrop = log(PressDrop);
PressDEV  = std(LogPressDrop,  'omitnan');   
PressMEAN = mean(LogPressDrop, 'omitnan');
PressDropSCALED = (LogPressDrop - PressMEAN) / PressDEV;
MaxDP_Scaled = (log(MaxDP) - PressMEAN) / PressDEV;

% FIX 2: NaN penalty only for lifespan GP; pressure GP is trained on valid
%         points only, so crashed samples are simply excluded from it.
ValidMask = ~isnan(Lifespan);
WorstLifespan = min(LifespanSCALED(ValidMask));
PenaltyLife   = WorstLifespan - 0.1;
LifespanSCALED(isnan(LifespanSCALED)) = PenaltyLife;
% (PressDropSCALED NaNs are left in place; they are masked out below)

% FIX 1: Normalize training inputs to [0,1]^D before building any GP.
%         Raw ranges span 0.05 to 60, so with log-length-scales initialised
%         at zero (=> length scale = 1) the kernel sees nearly all raw pairs
%         as uncorrelated, giving a flat GP and therefore flat EI.
GeometriesNorm = (Geometries - LowerBounds) ./ InputRange;

% Parameter initialization for Kernel
% FIX 1: Start length scales at log(0.5) rather than 0 so the initial
%         guess is "half the [0,1] input range", which is sensible.
lengthScales = log(0.5) * ones(NumDims, 1);
signalVar = 0;
noiseVar  = log(0.01);
logThetaInit = [lengthScales; signalVar; noiseVar];

% Set up optimization options
options_GP  = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton', 'MaxFunctionEvaluations', 1000);
options_ACQ = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp',          'MaxFunctionEvaluations', 1000);

% FIX 1+2: Pass normalized inputs; pressure GP trained on valid rows only.
logThetaOpt_LIFE  = fminunc(@(t) NLML(t, GeometriesNorm,              LifespanSCALED),             logThetaInit, options_GP);
logThetaOpt_PRESS = fminunc(@(t) NLML(t, GeometriesNorm(ValidMask,:), PressDropSCALED(ValidMask)), logThetaInit, options_GP);
ThetaOpt_LIFE  = exp(logThetaOpt_LIFE);
ThetaOpt_PRESS = exp(logThetaOpt_PRESS);

%% Optimizer
NumIterations = 100;

% FIX 1: Acquisition function and fmincon both operate in normalized [0,1]^D
%         space; results are de-normalized before calling the physics solver.
LB_norm = zeros(1, NumDims);
UB_norm = ones(1, NumDims);

for iter = 1:NumIterations

    % Recompute valid mask and normalized inputs for this iteration
    ValidMask      = ~isnan(Lifespan);
    GeometriesNorm = (Geometries - LowerBounds) ./ InputRange;
    GeometriesNorm_Valid = GeometriesNorm(ValidMask, :);
    NumTrain = size(GeometriesNorm, 1);
    NumValid = sum(ValidMask);

    % Find best feasible lifespan seen so far
    FeasibleIdx = ValidMask & (PressDrop <= MaxDP);
    if any(FeasibleIdx)
        BestValidScaled = max(LifespanSCALED(FeasibleIdx));
    else
        BestValidScaled = min(LifespanSCALED(ValidMask)); 
    end

    % Precompute Cholesky factors and alpha vectors
    % Life GP  -- all points (invalid ones carry penalty value)
    len_LIFE   = ThetaOpt_LIFE(1:end-2)';  sig_LIFE   = ThetaOpt_LIFE(end-1);  noise_LIFE   = ThetaOpt_LIFE(end);
    K_LIFE     = MaternKernel(GeometriesNorm, GeometriesNorm, len_LIFE, sig_LIFE) + noise_LIFE * eye(NumTrain);
    L_LIFE     = chol(K_LIFE, 'lower');
    alpha_LIFE = L_LIFE' \ (L_LIFE \ LifespanSCALED);

    % Pressure GP -- valid points only (FIX 2)
    len_PRESS   = ThetaOpt_PRESS(1:end-2)'; sig_PRESS   = ThetaOpt_PRESS(end-1); noise_PRESS   = ThetaOpt_PRESS(end);
    K_PRESS     = MaternKernel(GeometriesNorm_Valid, GeometriesNorm_Valid, len_PRESS, sig_PRESS) + noise_PRESS * eye(NumValid);
    L_PRESS     = chol(K_PRESS, 'lower');
    alpha_PRESS = L_PRESS' \ (L_PRESS \ PressDropSCALED(ValidMask));

    MaxDP_Scaled = (log(MaxDP) - PressMEAN) / PressDEV;

    % Global grid search -- sample directly in [0,1]^D (FIX 1)
    NumGrid = 50000;
    GridGeometriesNorm = rand(NumGrid, NumDims);

    % Predict on the grid
    [muLIFE_g,  stdLIFE_g]  = PredictGP(GridGeometriesNorm, GeometriesNorm,       L_LIFE,  alpha_LIFE,  ThetaOpt_LIFE);
    [muPRESS_g, stdPRESS_g] = PredictGP(GridGeometriesNorm, GeometriesNorm_Valid,  L_PRESS, alpha_PRESS, ThetaOpt_PRESS);
    stdLIFE_g  = max(real(stdLIFE_g),  1e-6);
    stdPRESS_g = max(real(stdPRESS_g), 1e-6);

    % Vectorised Expected Improvement
    % FIX 3: xi must be >= 0. Negative xi makes EI non-zero almost
    %         everywhere, which washes out acquisition discrimination.
    xi = 0.01;
    Z  = (muLIFE_g - BestValidScaled - xi) ./ stdLIFE_g;
    EI = stdLIFE_g .* (Z .* NormCDF(Z) + NormPDF(Z));

    % Vectorised Probability of Feasibility
    PoF_g = NormCDF((MaxDP_Scaled - muPRESS_g) ./ stdPRESS_g);

    GridScores = -EI .* PoF_g;

    % Multi-start fmincon -- starts and optimises in normalised space (FIX 1)
    [~, SortedIdx] = sort(GridScores);

    NumStarts      = 3;
    BestCandidates = zeros(NumStarts, NumDims);
    BestAcqVals    = zeros(NumStarts, 1);

    AcqObj = @(x_norm) Acquisition(x_norm, GeometriesNorm, GeometriesNorm_Valid, ...
                                    LB_norm, UB_norm, ...
                                    L_LIFE,  alpha_LIFE,  ThetaOpt_LIFE, ...
                                    L_PRESS, alpha_PRESS, ThetaOpt_PRESS, ...
                                    MaxDP_Scaled, BestValidScaled);

    for k = 1:NumStarts
        x0_norm = GridGeometriesNorm(SortedIdx(k), :);
        [x_opt_norm, fval] = fmincon(AcqObj, x0_norm, [], [], [], [], LB_norm, UB_norm, [], options_ACQ);
        BestCandidates(k, :) = x_opt_norm;
        BestAcqVals(k)       = fval;
    end

    % Select best candidate and de-normalize back to physical space
    [~, trueBestIdx] = min(BestAcqVals);
    x_next_norm = BestCandidates(trueBestIdx, :);
    x_next      = LowerBounds + x_next_norm .* InputRange;

    % Evaluate the physics
    WT = x_next(1:3); AR = x_next(4:6); CW = x_next(7:8); NC = x_next(9);
    fprintf('Testing new geometry (#%i out of %i): \n NC = %i, WT=[%.3f, %.3f, %.3f], AR=[%.3f, %.3f, %.3f], CW=[%.3f, %.3f]\n', iter, NumIterations, round(NC), x_next(1:8));
    try
        [Life_new, Drop_new] = SKIPPERRegen(Data, NC, WT, AR, CW, 0);
        if ~isreal(Life_new) || ~isreal(Drop_new)
            error('Complex physics output detected.'); 
        end
    catch
        Life_new = NaN; Drop_new = NaN;
        warning('Geometry failed in physics solver.');
    end

    if ~isnan(Life_new)
        fprintf('Result: %.2f Cycles, %.2f psi drop\n', Life_new, Drop_new);
    end

    % Append new observation
    Geometries = [Geometries; x_next];
    Lifespan   = [Lifespan;   Life_new];
    PressDrop  = [PressDrop;  Drop_new];

    % Re-scale outputs over the full dataset
    LifespanDEV  = std(Lifespan,  'omitnan'); 
    LifespanMEAN = mean(Lifespan, 'omitnan');
    LifespanSCALED = (Lifespan - LifespanMEAN) / LifespanDEV;

    LogPressDrop = log(PressDrop);
    PressDEV  = std(LogPressDrop,  'omitnan');   
    PressMEAN = mean(LogPressDrop, 'omitnan');
    PressDropSCALED = (LogPressDrop - PressMEAN) / PressDEV;
    MaxDP_Scaled    = (log(MaxDP)  - PressMEAN) / PressDEV;

    % NaN penalty for lifespan GP only (FIX 2)
    ValidMask     = ~isnan(Lifespan);
    WorstLifespan = min(LifespanSCALED(ValidMask));
    PenaltyLife   = WorstLifespan - 0.1;
    LifespanSCALED(isnan(LifespanSCALED)) = PenaltyLife;

    % Retrain both GPs (warm-started from previous optimum)
    GeometriesNorm = (Geometries - LowerBounds) ./ InputRange;
    logThetaOpt_LIFE  = fminunc(@(t) NLML(t, GeometriesNorm,              LifespanSCALED),             logThetaOpt_LIFE,  options_GP);
    logThetaOpt_PRESS = fminunc(@(t) NLML(t, GeometriesNorm(ValidMask,:), PressDropSCALED(ValidMask)), logThetaOpt_PRESS, options_GP);
    ThetaOpt_LIFE  = exp(logThetaOpt_LIFE);
    ThetaOpt_PRESS = exp(logThetaOpt_PRESS);
end

%% Plotting
ValidIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP);

if ~any(ValidIdx)
    warning('No valid geometries were found that satisfied the pressure constraint!');
else
    ValidLifespans  = Lifespan(ValidIdx);
    ValidPressDrops = PressDrop(ValidIdx);
    ValidGeometries = Geometries(ValidIdx, :);

    [BestLife, relative_idx] = max(ValidLifespans);
    BestPress    = ValidPressDrops(relative_idx);
    ChampionGeom = ValidGeometries(relative_idx, :);

    fprintf('CHAMPION PERFORMANCE:\n');
    fprintf('  Lifespan:      %.2f Cycles\n', BestLife);
    fprintf('  Pressure Drop: %.2f psi\n\n', BestPress);

    fprintf('OPTIMAL GEOMETRY [inches]:\n');
    fprintf('  Wall Thickness (Chamber, Throat, Nozzle): [%.4f, %.4f, %.4f]\n', ChampionGeom(1:3));
    fprintf('  Aspect Ratio (Chamber, Throat, Nozzle): [%.4f, %.4f, %.4f]\n',   ChampionGeom(4:6));
    fprintf('  Channel Width  (Chamber, Nozzle):         [%.4f, %.4f]\n',        ChampionGeom(7:8));
    fprintf('  Channel Count:                             %i \n',                 ChampionGeom(9));
    fprintf('======================================================\n');

    TotalEvals = length(Lifespan);
    BestSoFar  = zeros(TotalEvals, 1);
    CurrentMax = 0;

    for i = 1:TotalEvals
        if ~isnan(Lifespan(i)) && (PressDrop(i) <= MaxDP) && (Lifespan(i) > CurrentMax)
            CurrentMax = Lifespan(i);
        end
        BestSoFar(i) = CurrentMax;
    end

    figure('Name', 'Bayesian Optimization Convergence', 'Position', [100, 100, 800, 600]);

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

    subplot(2, 1, 2);
    hold on; set(gca, 'FontName', 'Times New Roman');
    scatter(1:TotalEvals, PressDrop, 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
    yline(MaxDP, 'r-', 'Max Constraint (37.5 psi)', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');
    xline(NumSamples, 'k--', 'LineWidth', 1.5);
    title('Constraint Tracking: Coolant Pressure Drop');
    xlabel('Evaluation Number (Initialization + Active Search)');
    ylabel('Pressure Drop [psi]');
    grid on;
    ylim([0, max(MaxDP + 50, prctile(PressDrop, 85))]);

    % Call SKIPPERRegen for Final Native Plots
    % SKIPPERRegen(Data, ChampionGeom(9), ChampionGeom(1:3), ChampionGeom(4:6), ChampionGeom(7:8), 1);
end