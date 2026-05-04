%% Code to automatically size Regen channel dimensions to optimize for engine life
% Pablo Plata
% Optimized using Relaxation & Fine-Tuning Mixed-Integer BO
% 
% INPUTS:   WallThickness -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelHeight -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelWidth  -> 1x2 array, [Chamber, Nozzle]
%
% Design vector for solver:
% X = [WallThickness, ChannelHeight, ChannelWidth, ChannelCount];

clear;
clc;
close all;
addpath("cea\");
addpath("IPA Data\");
addpath("Material Data\");
addpath("Contours\");
addpath('BayesianOpt\');
Data = LoadData();
rng(42);

%% Parameter sampling (inches)
Thickness = [0.05; 0.1];    AR = [0.3; 3.0];  Width = [0.01; 0.125];
LowerBounds = [ones(1, 3) * Thickness(1), ones(1, 3) * AR(1), ones(1, 2) * Width(1), 40];
UpperBounds = [ones(1, 3) * Thickness(2), ones(1, 3) * AR(2), ones(1, 2) * Width(2), 80];
InputRange  = UpperBounds - LowerBounds;    
MaxDP = 150; % * 0.5^2; % psi

% Latin Hypercube Sampling for parameter space for GP training
NumDims = length(LowerBounds);
NumSamples = 150;
LHS = HyperSampl(NumSamples, NumDims);
Geometries = LowerBounds + LHS .* InputRange;

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

%% Phase 1: Continuous Global Optimization
NumPhase1 = 100;
fprintf('\nStarting Phase 1: Continuous Search\n');
[Geometries, Lifespan, PressDrop, logT_L, logT_P] = BOSearch( ...
    NumPhase1, Geometries, Lifespan, PressDrop, ...
    LowerBounds, UpperBounds, LowerBounds, UpperBounds, ...
    Data, MaxDP, [], [], []);

% Identify Continuous Champion
ValidIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP);
if ~any(ValidIdx)
    error('Phase 1 found no valid designs satisfying constraints.');
end
ValidLifespans = Lifespan(ValidIdx);
ValidGeometries = Geometries(ValidIdx, :);
[~, bestIdx] = max(ValidLifespans);
ContChampion = ValidGeometries(bestIdx, :);
fprintf('\nGlobal Optimization Complete. Continuous Champion\n');
fprintf('  WT=[%.4f, %.4f, %.4f], AR=[%.4f, %.4f, %.4f], CW=[%.4f, %.4f], NC=%i\n', ...
        ContChampion(1:8), round(ContChampion(9)));

%% Phase 2: Bracketing & Fine-Tuning
fprintf('\nStarting Phase 2: Discrete Fine-Tuning\n');
% Expanded standard slitting saw thicknesses (inches)
CW_Avail = [
    0.0100; % 0.01" (10 thou)
    0.0120; % 0.012"
    0.0130; % 0.013"
    0.0140; % 0.014"
    0.0156; % 1/64"
    0.0160; % 0.016"
    0.0180; % 0.018"
    0.0200; % 0.02"
    0.0230; % 0.023"
    0.0250; % 0.025"
    0.0280; % 0.028"
    0.0312; % 1/32"
    0.0320; % 0.032"
    0.0350; % 0.035"
    0.0360; % 0.036"
    0.0400; % 0.04"
    0.0450; % 0.045"
    0.0469; % 3/64"
    0.0510; % 0.051"
    0.0570; % 0.057"
    0.0625; % 1/16"
    0.0640; % 0.064"
    0.0720; % 0.072"
    0.0781; % 5/64"
    0.0810; % 0.081"
    0.0910; % 0.091"
    0.0938; % 3/32"
    0.1020; % 0.102"
    0.1140; % 0.114"
    0.1250  % 1/8"
];

% Bracket Chamber Width
diffC = CW_Avail - ContChampion(7);
idxC_L = find(diffC <= 0, 1, 'last');  
if isempty(idxC_L)
    idxC_L = 1; 
end
idxC_H = find(diffC >= 0, 1, 'first'); 
if isempty(idxC_H)
    idxC_H = length(CW_Avail);
end
C_Opts = unique([CW_Avail(idxC_L), CW_Avail(idxC_H)]);

% Bracket Nozzle / Throat Width
diffN = CW_Avail - ContChampion(8);
idxN_L = find(diffN <= 0, 1, 'last');
if isempty(idxN_L)
    idxN_L = 1;
end
idxN_H = find(diffN >= 0, 1, 'first'); 
if isempty(idxN_H)
    idxN_H = length(CW_Avail);
end
N_Opts = unique([CW_Avail(idxN_L), CW_Avail(idxN_H)]);

% Generate configurations
[C_Grid, N_Grid] = ndgrid(C_Opts, N_Opts);
FineTuneCWs = [C_Grid(:), N_Grid(:)];

% Restrict search space for fine-tuning
SearchRadius = 0.50 * InputRange;
FT_LB = max(LowerBounds, ContChampion - SearchRadius);
FT_UB = min(UpperBounds, ContChampion + SearchRadius);

NumPhase2_Init = 70;
NumPhase2_Search = 70;
FT_EndIndices = []; 
CurrentTotalEvals = NumSamples + NumPhase1;

for k = 1:size(FineTuneCWs, 1)
    FixedCW = FineTuneCWs(k, :);
    fprintf('\nFine-Tuning config %i/%i: CW = [%.4f, %.4f]\n', k, size(FineTuneCWs, 1), FixedCW(1), FixedCW(2));
    
    % Generate Local LHS
    LocalLHS = HyperSampl(NumPhase2_Init, NumDims);
    LocalGeom = FT_LB + LocalLHS .* (FT_UB - FT_LB);
    LocalGeom(:, 7:8) = repmat(FixedCW, NumPhase2_Init, 1); 
    
    LocalLife = zeros(NumPhase2_Init, 1);
    LocalPress = zeros(NumPhase2_Init, 1);
    
    fprintf('  Evaluating %i local LHS samples... ', NumPhase2_Init);
    for i = 1:NumPhase2_Init
        try
            [LocalLife(i), LocalPress(i)] = SKIPPERRegen(Data, LocalGeom(i, 9), LocalGeom(i, 1:3), LocalGeom(i, 4:6), FixedCW, 0);
        catch
            LocalLife(i) = NaN; LocalPress(i) = NaN;
        end
    end
    fprintf('Done.\n');
    
    % Cold-Start Local BO Search (FT bounds act as absolute normalization bounds)
    [OptGeom, OptLife, OptPress, ~, ~] = BOSearch( ...
        NumPhase2_Search, LocalGeom, LocalLife, LocalPress, ...
        FT_LB, FT_UB, FT_LB, FT_UB, ... 
        Data, MaxDP, FixedCW, [], []);
        
    % Append to master array for final plotting
    Geometries = [Geometries; OptGeom];
    Lifespan = [Lifespan; OptLife];
    PressDrop = [PressDrop; OptPress];
    CurrentTotalEvals = CurrentTotalEvals + NumPhase2_Init + NumPhase2_Search;
    FT_EndIndices = [FT_EndIndices, CurrentTotalEvals];
end

%% Plotting & Final Selection
TotalEvals = length(Lifespan);

% 1. Identify MANUFACTURABLE Champion (Phase 2 Only)
Phase2Start = NumSamples + NumPhase1 + 1;
Phase2Mask = false(TotalEvals, 1);
Phase2Mask(Phase2Start:end) = true;

% Define all valid indices for text output and plotting
ValidMfgIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP) & Phase2Mask;
ValidGlobalIdx = ~isnan(Lifespan) & (PressDrop <= MaxDP);

if ~any(ValidMfgIdx)
    warning('No valid manufacturable geometries found in Phase 2.');
else
    [BestLife, relative_idx] = max(Lifespan(ValidMfgIdx));
    
    % Map back to absolute index to grab correct geometry
    ValidMfgIndices = find(ValidMfgIdx);
    AbsBestIdx = ValidMfgIndices(relative_idx);
    
    BestPress    = PressDrop(AbsBestIdx);
    ChampionGeom = Geometries(AbsBestIdx, :);

    fprintf('\n======================================================\n');
    fprintf(' CHAMPION MANUFACTURABLE PERFORMANCE:\n');
    fprintf('  Lifespan:      %.2f Cycles\n', BestLife);
    fprintf('  Pressure Drop: %.2f psi\n\n', BestPress);
    fprintf(' OPTIMAL GEOMETRY [inches]:\n');
    fprintf('  Wall Thickness (C, T, N): [%.4f, %.4f, %.4f]\n', ChampionGeom(1:3));
    fprintf('  Aspect Ratio (C, T, N):   [%.4f, %.4f, %.4f]\n', ChampionGeom(4:6));
    fprintf('  Channel Width (C, N):     [%.4f, %.4f]\n',       ChampionGeom(7:8));
    fprintf('  Channel Count:            %i \n',                ChampionGeom(9));
    fprintf('======================================================\n');
end

% 2. Calculate "Best So Far" globally for the trend line
BestSoFar  = zeros(TotalEvals, 1);
CurrentMax = 0;
for i = 1:TotalEvals
    if ValidGlobalIdx(i) && (Lifespan(i) > CurrentMax)
        CurrentMax = Lifespan(i);
    end
    BestSoFar(i) = CurrentMax;
end

% 3. Build Figures
figure('Name', 'Bayesian Optimization Convergence', 'Position', [100, 100, 800, 600]);

% Segregate points: Crashes (NaN), Constraint Fails (>MaxDP), and Valid
CrashIdx = isnan(Lifespan) | isnan(PressDrop);
ConstraintFailIdx = ~CrashIdx & (PressDrop > MaxDP);
FailedIdx = CrashIdx | ConstraintFailIdx;

% Force crashes to y=0 for visibility, leave constraint fails at actual values
PlotLife = Lifespan; PlotLife(CrashIdx) = 0;
PlotPress = PressDrop; PlotPress(CrashIdx) = 0;

% Plot 1: Lifespan
subplot(2, 1, 1);
hold on; set(gca, 'FontName', 'Times New Roman');

% Plot failed geometries (Gray)
scatter(find(FailedIdx), PlotLife(FailedIdx), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
% Plot valid geometries (Black)
scatter(find(ValidGlobalIdx), PlotLife(ValidGlobalIdx), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
plot(1:TotalEvals, BestSoFar, 'b-', 'LineWidth', 2.5);

xline(NumSamples, 'k-', 'End LHS Init', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
xline(NumSamples + NumPhase1, 'g-', 'End Phase 1', 'LineWidth', 1.5, 'LabelVerticalAlignment', 'bottom');
for p = 1:length(FT_EndIndices)
    xline(FT_EndIndices(p), 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
end

title('Objective Convergence: Engine Lifespan');
ylabel('Lifespan [Cycles]');
legend('Failed (Crashed or > Max DP)', 'Valid Geometries', 'Optimal Valid Lifespan', 'Location', 'northwest');
grid on; ylim([0, max(PlotLife) * 1.1]);

% Plot 2: Pressure Drop
subplot(2, 1, 2);
hold on; set(gca, 'FontName', 'Times New Roman');

% Plot failed geometries (Gray)
scatter(find(FailedIdx), PlotPress(FailedIdx), 20, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
% Plot valid geometries (Red)
scatter(find(ValidGlobalIdx), PlotPress(ValidGlobalIdx), 20, 'r', 'filled', 'MarkerFaceAlpha', 0.5);

yline(MaxDP, 'r-', 'Max Constraint', 'LineWidth', 2, 'LabelHorizontalAlignment', 'left');

xline(NumSamples, 'k-', 'LineWidth', 1.5);
xline(NumSamples + NumPhase1, 'g-', 'LineWidth', 1.5);
for p = 1:length(FT_EndIndices)
    xline(FT_EndIndices(p), 'k--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.0);
end

title('Constraint Tracking: Coolant Pressure Drop');
xlabel('Evaluation Number');
ylabel('Pressure Drop [psi]');
grid on; ylim([0, max(MaxDP + 50, prctile(PlotPress(~CrashIdx), 85))]);
%SKIPPERRegen(Data, ChampionGeom(9), ChampionGeom(1:3), ChampionGeom(4:6), ChampionGeom(7:8), 1);
%% Local Functions
function [Geometries, Lifespan, PressDrop, logT_L, logT_P] = BOSearch(NumIter, Geometries, Lifespan, PressDrop, SearchLB, SearchUB, GlobalLB, GlobalUB, Data, MaxDP, FixedCW, logT_L, logT_P)
    
    NumDims = length(GlobalLB);
    GlobalRange = GlobalUB - GlobalLB;

    % Fix dimensions 7 and 8 if provided
    if ~isempty(FixedCW)
        SearchLB(7:8) = FixedCW;
        SearchUB(7:8) = FixedCW;
    end

    LB_norm = (SearchLB - GlobalLB) ./ GlobalRange;
    UB_norm = (SearchUB - GlobalLB) ./ GlobalRange;

    % Initialize hyperparams if not warm-starting
    if isempty(logT_L)
        logT_L = [log(0.5) * ones(NumDims, 1); 0; log(0.01)];
        logT_P = logT_L;
    end

    options_GP  = optimoptions('fminunc', 'Display', 'off', 'Algorithm', 'quasi-newton', 'MaxFunctionEvaluations', 1000);
    options_ACQ = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp',          'MaxFunctionEvaluations', 1000);

    for iter = 1:NumIter
        % Output scaling
        LifespanDEV  = std(Lifespan,  'omitnan'); 
        LifespanMEAN = mean(Lifespan, 'omitnan');
        LifespanSCALED = (Lifespan - LifespanMEAN) / LifespanDEV;

        LogPressDrop = log(PressDrop);
        PressDEV  = std(LogPressDrop,  'omitnan');   
        PressMEAN = mean(LogPressDrop, 'omitnan');
        PressDropSCALED = (LogPressDrop - PressMEAN) / PressDEV;
        MaxDP_Scaled = (log(MaxDP) - PressMEAN) / PressDEV;

        ValidMask = ~isnan(Lifespan);
        WorstLifespan = min(LifespanSCALED(ValidMask));
        LifespanSCALED(isnan(LifespanSCALED)) = WorstLifespan - 0.1;

        GeometriesNorm = (Geometries - GlobalLB) ./ GlobalRange;
        GeometriesNorm_Valid = GeometriesNorm(ValidMask, :);
        NumValid = sum(ValidMask);

        FeasibleIdx = ValidMask & (PressDrop <= MaxDP);
        if any(FeasibleIdx)
            BestValidScaled = max(LifespanSCALED(FeasibleIdx));
        else
            BestValidScaled = min(LifespanSCALED(ValidMask)); 
        end

        % GP Training
        logT_L = fminunc(@(t) NLML(t, GeometriesNorm,              LifespanSCALED),             logT_L, options_GP);
        logT_P = fminunc(@(t) NLML(t, GeometriesNorm_Valid, PressDropSCALED(ValidMask)), logT_P, options_GP);
        ThetaOpt_LIFE  = exp(logT_L);
        ThetaOpt_PRESS = exp(logT_P);

        len_L = ThetaOpt_LIFE(1:end-2)';  sig_L = ThetaOpt_LIFE(end-1);  noise_L = ThetaOpt_LIFE(end);
        K_L   = MaternKernel(GeometriesNorm, GeometriesNorm, len_L, sig_L) + noise_L * eye(size(GeometriesNorm, 1));
        
        % Progressive jitter for Life GP
        [L_L, pL] = chol(K_L, 'lower');
        jitterL = 1e-6;
        while pL > 0
            [L_L, pL] = chol(K_L + jitterL * max(diag(K_L)) * eye(size(K_L,1)), 'lower');
            jitterL = jitterL * 10;
        end
        a_L   = L_L' \ (L_L \ LifespanSCALED);

        len_P = ThetaOpt_PRESS(1:end-2)'; sig_P = ThetaOpt_PRESS(end-1); noise_P = ThetaOpt_PRESS(end);
        K_P   = MaternKernel(GeometriesNorm_Valid, GeometriesNorm_Valid, len_P, sig_P) + noise_P * eye(NumValid);
        
        % Progressive jitter for Pressure GP
        [L_P, pP] = chol(K_P, 'lower');
        jitterP = 1e-6;
        while pP > 0
            [L_P, pP] = chol(K_P + jitterP * max(diag(K_P)) * eye(size(K_P,1)), 'lower');
            jitterP = jitterP * 10;
        end
        a_P   = L_P' \ (L_P \ PressDropSCALED(ValidMask));

        % Grid Search
        NumGrid = 50000;
        GridGeometriesNorm = LB_norm + rand(NumGrid, NumDims) .* (UB_norm - LB_norm);

        [muLIFE_g,  stdLIFE_g]  = PredictGP(GridGeometriesNorm, GeometriesNorm,       L_L, a_L, ThetaOpt_LIFE);
        [muPRESS_g, stdPRESS_g] = PredictGP(GridGeometriesNorm, GeometriesNorm_Valid, L_P, a_P, ThetaOpt_PRESS);
        stdLIFE_g  = max(real(stdLIFE_g),  1e-6);
        stdPRESS_g = max(real(stdPRESS_g), 1e-6);

        xi = 0.01;
        Z  = (muLIFE_g - BestValidScaled - xi) ./ stdLIFE_g;
        EI = stdLIFE_g .* (Z .* NormCDF(Z) + NormPDF(Z));
        PoF_g = NormCDF((MaxDP_Scaled - muPRESS_g) ./ stdPRESS_g);
        GridScores = -EI .* PoF_g;

        [~, SortedIdx] = sort(GridScores);

        % Optimizer
        NumStarts = 2;
        BestCandidates = zeros(NumStarts, NumDims);
        BestAcqVals    = zeros(NumStarts, 1);

        AcqObj = @(x_norm) Acquisition(x_norm, GeometriesNorm, GeometriesNorm_Valid, ...
                                        LB_norm, UB_norm, L_L, a_L, ThetaOpt_LIFE, L_P, a_P, ThetaOpt_PRESS, MaxDP_Scaled, BestValidScaled);

        for j = 1:NumStarts
            x0_norm = GridGeometriesNorm(SortedIdx(j), :);
            [x_opt_norm, fval] = fmincon(AcqObj, x0_norm, [], [], [], [], LB_norm, UB_norm, [], options_ACQ);
            BestCandidates(j, :) = x_opt_norm;
            BestAcqVals(j)       = fval;
        end

        [~, trueBestIdx] = min(BestAcqVals);
        x_next = GlobalLB + BestCandidates(trueBestIdx, :) .* GlobalRange;

        % Physical evaluation
        WT = x_next(1:3); AR = x_next(4:6); CW = x_next(7:8); NC = x_next(9);
        try
            [Life_new, Drop_new] = SKIPPERRegen(Data, NC, WT, AR, CW, 0);
            if ~isreal(Life_new) || ~isreal(Drop_new), error('Complex output.'); end
        catch
            Life_new = NaN; Drop_new = NaN;
        end
        if ~isnan(Life_new)
            fprintf('Eval Complete! %.2f Cycles & %.2f psi drop\n', Life_new, Drop_new);
        end

        Geometries = [Geometries; x_next];
        Lifespan   = [Lifespan;   Life_new];
        PressDrop  = [PressDrop;  Drop_new];
    end
end