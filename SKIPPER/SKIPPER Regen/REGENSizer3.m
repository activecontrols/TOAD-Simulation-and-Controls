%% Code to automatically size Regen channel dimensions to optimize for engine life
% Pablo Plata
% Optimized using Differential Evolution (global) + fmincon polish (local)
%
% Replaces the GP-based Bayesian Optimization approach now that physics
% calls run at ~40ms instead of ~2s. At this call speed the GP surrogate's
% own overhead (O(n^3) Cholesky, repeated fminunc hyperparameter fits,
% 25k-point acquisition grids) exceeds the cost of just running more
% physics evaluations directly, so a population-based global search is
% used instead. DE handles the non-smooth constraint landscape (solver
% crashes -> NaN, hard OD limit) more robustly than a GP acquisition
% function, requires no surrogate maintenance, and its wall-clock cost is
% dominated almost entirely by the physics calls themselves.
%
% INPUTS:   WallThickness -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelHeight -> 1x3 array, [Chamber, Throat, Nozzle]
%           ChannelWidth  -> 1x2 array, [Chamber, Nozzle]  (3rd slot unused, kept for LinerODCalc compat)
%           ChannelCount  -> scalar, now part of the search space
%
% Design vector for solver:
% X = [WallThickness(3), ChannelHeight/AR(3), ChannelWidth(3), ChannelCount(1)]   -> 10D

clear;
clc;
close all;
set(0, 'DefaultFigureWindowStyle', 'docked');  % dock all figures instead of separate windows
addpath("cea\");
addpath("IPA Data\");
addpath("Material Data\");
addpath("Contours\");
addpath('BayesianOpt\');   % HyperSampl lives here
Data = LoadData();
rng(41);

%% Parameter bounds (inches, and channel count)
Thickness = [0.05; 0.1];    AR = [0.1; 5];  Width = [0.02; 0.125];
NC_Bounds = [50; 65]; 

LowerBounds = [ones(1, 3) * Thickness(1), ones(1, 3) * AR(1), ones(1, 3) * Width(1), NC_Bounds(1)];
UpperBounds = [ones(1, 3) * Thickness(2), ones(1, 3) * AR(2), ones(1, 3) * Width(2), NC_Bounds(2)];
InputRange  = UpperBounds - LowerBounds;
NumDims     = length(LowerBounds);     % 10
IntDim      = NumDims;                 % index of the integer (channel count) dimension

MaxDP = 50;     % psi, hard pressure-drop constraint (soft-penalized in cost)
ODLimit = 3.95; % in, max allowable liner OD for manufacturability (hard constraint)
R_chamber = Data.Contour(1,2);

%% Algorithm / budget settings
SecondaryTempWeight = 0.35;  % relative weight of the non-throat cooling objective vs life
ThroatExcludeFrac   = 0.06;  % fraction of profile length excluded on EACH SIDE of the peak (throat) station
                              % before averaging, so the "secondary" temp reflects chamber+nozzle
                              % cooling only, not the throat spike that Life already accounts for
NumInit          = 1200;   % LHS samples used both to seed DE and to set objective scaling
PopSize          = 75;     % DE population size
MaxGenerations   = 600;    % soft ceiling; wall-clock budget below is the real stopping rule
F_weight         = 0.8;    % DE differential weight
CR               = 0.9;    % DE crossover probability
TimeBudget_DE    = 300;    % seconds allotted to the DE phase
NumPolishStarts  = 4;      % baseline multi-start fmincon polish after DE (auto-increased if DE finishes with time left)
SecPerPolishStart = 4;     % rough time cost (s) of one extra polish multi-start, used to size the bonus below
MaxExtraPolishStarts = 1;  % cap on how many bonus polish starts freed-up DE time can buy

%% Quality-Diversity archive settings (merged from BoundaryExplorer.m)
% Tracks, for each of the 10 design dimensions, the best (lowest) cost
% seen so far within each of NumBins slices of that dimension's range.
% Any evaluation that breaks a bin record is a "boundary pusher" and gets
% injected into the active DE population in place of its worst member.
NumBins = 20;

fprintf('Generating %i LHS samples for DE seeding + objective scaling...\n', NumInit);
LHS = HyperSampl(NumInit, NumDims);
InitGeom = LowerBounds + LHS .* InputRange;

InitLife = nan(NumInit, 1);
InitDrop = nan(NumInit, 1);
InitMeanTemp = nan(NumInit, 1);
InitPeakTemp = nan(NumInit, 1);
InitSecondaryTemp = nan(NumInit, 1);

for i = 1:NumInit
    [InitLife(i), InitDrop(i), InitMeanTemp(i), InitPeakTemp(i), InitSecondaryTemp(i)] = ...
        EvaluatePhysicsRaw(InitGeom(i, :), Data, R_chamber, ODLimit, ThroatExcludeFrac);
end

ValidInit = ~isnan(InitLife);
if ~any(ValidInit)
    error('No valid geometries found in the initial LHS batch — check bounds/physics function.');
end

% Scaling references (median of valid samples; robust to outliers, no
% need for full standardization like the GP version required)
LifeScale = median(InitLife(ValidInit));
SecondaryScale = median(InitSecondaryTemp(ValidInit));
if LifeScale <= 0 || isnan(LifeScale); LifeScale = 1; end
if SecondaryScale <= 0 || isnan(SecondaryScale); SecondaryScale = 1; end

fprintf('Init complete. Valid: %i/%i. LifeScale=%.2f, SecondaryScale=%.2f (throat-excluded mean)\n', ...
    sum(ValidInit), NumInit, LifeScale, SecondaryScale);

%% Set up QD boundary archives now that bounds/dims are known
Archives = cell(NumDims, 1);
for d = 1:NumDims
    Archives{d} = struct('MinCost', inf(NumBins, 1), ...
                          'Geom', nan(NumBins, NumDims), ...
                          'Life', nan(NumBins, 1), ...
                          'Drop', nan(NumBins, 1), ...
                          'Temp', nan(NumBins, 1), ...
                          'Peak', nan(NumBins, 1), ...
                          'Secondary', nan(NumBins, 1));
    Archives{d}.Edges = linspace(LowerBounds(d), UpperBounds(d), NumBins + 1);
    Archives{d}.Edges(1) = -inf; Archives{d}.Edges(end) = inf; % catch float rounding at the extremes
end

%% Pre-allocate evaluation log (DE fills this; polish appends after)
MaxLogRows = NumInit + PopSize * MaxGenerations + (NumPolishStarts + MaxExtraPolishStarts) * 200;
LogGeom     = nan(MaxLogRows, NumDims);
LogLife     = nan(MaxLogRows, 1);
LogDrop     = nan(MaxLogRows, 1);
LogMeanTemp = nan(MaxLogRows, 1);
LogPeakTemp = nan(MaxLogRows, 1);
LogSecondaryTemp = nan(MaxLogRows, 1);
LogCost     = nan(MaxLogRows, 1);
LogCount    = 0;

% Log the initial LHS batch, seeding the archives as we go
for i = 1:NumInit
    LogCount = LogCount + 1;
    LogGeom(LogCount, :) = InitGeom(i, :);
    LogLife(LogCount)    = InitLife(i);
    LogDrop(LogCount)    = InitDrop(i);
    LogMeanTemp(LogCount)= InitMeanTemp(i);
    LogPeakTemp(LogCount)= InitPeakTemp(i);
    LogSecondaryTemp(LogCount) = InitSecondaryTemp(i);
    LogCost(LogCount)    = CostFromRaw(InitLife(i), InitDrop(i), InitSecondaryTemp(i), MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);

    if ~isnan(InitLife(i)) && InitDrop(i) <= MaxDP
        Archives = UpdateArchives(InitGeom(i, :), LogCost(LogCount), InitLife(i), InitDrop(i), InitMeanTemp(i), InitPeakTemp(i), InitSecondaryTemp(i), NumDims, Archives);
    end
end

DE_EndIndex = LogCount; %#ok<NASGU> (kept for reference, set properly after DE below)

%% Phase 1: Differential Evolution (global search) with QD diversity injection
fprintf('\nStarting Differential Evolution (PopSize=%i, up to %i generations, %.0fs budget)\n', ...
    PopSize, MaxGenerations, TimeBudget_DE);

% Seed initial population from the best NumInit samples (fill any shortfall randomly)
[~, sortIdx] = sort(LogCost(1:NumInit));
SeedIdx = sortIdx(1:min(PopSize, NumInit));
Pop = InitGeom(SeedIdx, :);
PopCost = LogCost(SeedIdx);
PopLife = InitLife(SeedIdx);
PopDrop = InitDrop(SeedIdx);
PopTemp = InitMeanTemp(SeedIdx);
PopPeak = InitPeakTemp(SeedIdx);
PopSecondary = InitSecondaryTemp(SeedIdx);

if size(Pop, 1) < PopSize
    NumFill = PopSize - size(Pop, 1);
    FillLHS = HyperSampl(NumFill, NumDims);
    FillGeom = LowerBounds + FillLHS .* InputRange;
    FillLife = nan(NumFill, 1); FillDrop = nan(NumFill, 1); FillTemp = nan(NumFill, 1); FillPeak = nan(NumFill, 1); FillSecondary = nan(NumFill, 1); FillCost = nan(NumFill, 1);
    for i = 1:NumFill
        [FillLife(i), FillDrop(i), FillTemp(i), FillPeak(i), FillSecondary(i)] = ...
            EvaluatePhysicsRaw(FillGeom(i, :), Data, R_chamber, ODLimit, ThroatExcludeFrac);
        FillCost(i) = CostFromRaw(FillLife(i), FillDrop(i), FillSecondary(i), MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);
        LogCount = LogCount + 1;
        LogGeom(LogCount, :) = FillGeom(i, :);
        LogLife(LogCount) = FillLife(i); LogDrop(LogCount) = FillDrop(i);
        LogMeanTemp(LogCount) = FillTemp(i); LogPeakTemp(LogCount) = FillPeak(i);
        LogSecondaryTemp(LogCount) = FillSecondary(i); LogCost(LogCount) = FillCost(i);

        if ~isnan(FillLife(i)) && FillDrop(i) <= MaxDP
            Archives = UpdateArchives(FillGeom(i, :), FillCost(i), FillLife(i), FillDrop(i), FillTemp(i), FillPeak(i), FillSecondary(i), NumDims, Archives);
        end
    end
    Pop = [Pop; FillGeom]; PopCost = [PopCost; FillCost];
    PopLife = [PopLife; FillLife]; PopDrop = [PopDrop; FillDrop]; PopTemp = [PopTemp; FillTemp];
    PopPeak = [PopPeak; FillPeak]; PopSecondary = [PopSecondary; FillSecondary];
end

StartTimer = tic;
BestCostHistory = nan(MaxGenerations, 1);
GensRun = 0;

for gen = 1:MaxGenerations
    if toc(StartTimer) > TimeBudget_DE
        fprintf('Time budget reached after %i generations.\n', gen - 1);
        break;
    end

    ArchiveInjections = 0;

    for i = 1:PopSize
        % Mutation: pick 3 distinct individuals != i
        idxPool = setdiff(1:PopSize, i);
        r = idxPool(randperm(length(idxPool), 3));
        Donor = Pop(r(1), :) + F_weight * (Pop(r(2), :) - Pop(r(3), :));
        Donor = min(max(Donor, LowerBounds), UpperBounds);

        % Binomial crossover
        Trial = Pop(i, :);
        CrossMask = rand(1, NumDims) < CR;
        CrossMask(randi(NumDims)) = true; % ensure at least one dim from donor
        Trial(CrossMask) = Donor(CrossMask);

        [TrialLife, TrialDrop, TrialTemp, TrialPeak, TrialSecondary] = EvaluatePhysicsRaw(Trial, Data, R_chamber, ODLimit, ThroatExcludeFrac);
        TrialCost = CostFromRaw(TrialLife, TrialDrop, TrialSecondary, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);

        LogCount = LogCount + 1;
        LogGeom(LogCount, :) = Trial;
        LogLife(LogCount) = TrialLife; LogDrop(LogCount) = TrialDrop;
        LogMeanTemp(LogCount) = TrialTemp; LogPeakTemp(LogCount) = TrialPeak;
        LogSecondaryTemp(LogCount) = TrialSecondary; LogCost(LogCount) = TrialCost;

        InjectedByArchive = false;
        if ~isnan(TrialLife) && TrialDrop <= MaxDP
            [Archives, IsBoundaryPusher] = UpdateArchives(Trial, TrialCost, TrialLife, TrialDrop, TrialTemp, TrialPeak, TrialSecondary, NumDims, Archives);
            if IsBoundaryPusher
                % Diversity injection: this design set a new best-cost
                % record in some 1D slice of some dimension. Force it
                % into the population by replacing the current worst
                % member, regardless of whether it beats its own target
                % index i. This is what keeps the population spread out
                % and productive for the full time budget, replacing the
                % old stall-based early exit.
                [~, worstIdx] = max(PopCost);
                Pop(worstIdx, :) = Trial;
                PopCost(worstIdx) = TrialCost;
                PopLife(worstIdx) = TrialLife; PopDrop(worstIdx) = TrialDrop;
                PopTemp(worstIdx) = TrialTemp; PopPeak(worstIdx) = TrialPeak; PopSecondary(worstIdx) = TrialSecondary;
                ArchiveInjections = ArchiveInjections + 1;
                InjectedByArchive = true;
            end
        end

        % Standard greedy selection (skip if this trial was already
        % placed via archive injection above, to avoid double-counting)
        if ~InjectedByArchive && TrialCost < PopCost(i)
            Pop(i, :) = Trial;
            PopCost(i) = TrialCost;
            PopLife(i) = TrialLife; PopDrop(i) = TrialDrop; PopTemp(i) = TrialTemp;
            PopPeak(i) = TrialPeak; PopSecondary(i) = TrialSecondary;
        end
    end

    GensRun = gen;
    BestCostHistory(gen) = min(PopCost);
    if mod(gen, 5) == 0 || gen == 1
        [bc, bi] = min(PopCost);
        fprintf('Gen %3i | Best cost %.4f | Life %.1f cyc | Drop %.1f psi | ThroatT %.1f K | SecondaryT %.1f K | Archive injections %3i | t=%.0fs\n', ...
            gen, bc, PopLife(bi), PopDrop(bi), PopPeak(bi), PopSecondary(bi), ArchiveInjections, toc(StartTimer));
    end
end

[~, bestPopIdx] = min(PopCost);
DE_Best = Pop(bestPopIdx, :);
fprintf('\nDE complete after %i generations (%.0fs). Best cost = %.4f\n', GensRun, toc(StartTimer), PopCost(bestPopIdx));

% Roll any DE time budget left unused (e.g. MaxGenerations was hit before
% the wall-clock limit) into a more thorough local polish.
TimeLeftFromDE = TimeBudget_DE - toc(StartTimer);
if TimeLeftFromDE > SecPerPolishStart
    BonusStarts = min(MaxExtraPolishStarts, floor(TimeLeftFromDE / SecPerPolishStart));
    NumPolishStarts = NumPolishStarts + BonusStarts;
    fprintf('DE finished with %.0fs of its budget unused -> adding %i bonus polish starts (now %i total).\n', ...
        TimeLeftFromDE, BonusStarts, NumPolishStarts);
end

%% Phase 2: Local polish (fmincon, continuous dims only, NC fixed)
fprintf('\nStarting local polish around DE best...\n');
NC_Fixed = round(DE_Best(IntDim));
NC_Fixed = min(max(NC_Fixed, NC_Bounds(1)), NC_Bounds(2));

x9_LB = LowerBounds(1:9); x9_UB = UpperBounds(1:9);
ODNonlcon = @(x9) deal(LinerODCalc(x9, R_chamber) - ODLimit, []);
CostWrapper = @(x9) CostWrapperFcn(x9, NC_Fixed, Data, R_chamber, ODLimit, ThroatExcludeFrac, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);

options_Polish = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp', ...
    'StepTolerance', 1e-6, 'OptimalityTolerance', 1e-6, 'MaxFunctionEvaluations', 200);

PolishStarts = [DE_Best(1:9); ...
    min(max(DE_Best(1:9) + 0.05*InputRange(1:9).*(2*rand(NumPolishStarts-1, 9)-1), x9_LB), x9_UB)];

PolishResults = nan(NumPolishStarts, 9);
PolishCosts   = nan(NumPolishStarts, 1);

for k = 1:NumPolishStarts
    try
        [x_opt, fval] = fmincon(CostWrapper, PolishStarts(k, :), [], [], [], [], x9_LB, x9_UB, ODNonlcon, options_Polish);
        PolishResults(k, :) = x_opt;
        PolishCosts(k) = fval;
    catch
        PolishResults(k, :) = PolishStarts(k, :);
        PolishCosts(k) = CostWrapper(PolishStarts(k, :));
    end

    % Log the final physics evaluation of this polish result
    FullX = [PolishResults(k, :), NC_Fixed];
    [PLife, PDrop, PTemp, PPeak, PSecondary] = EvaluatePhysicsRaw(FullX, Data, R_chamber, ODLimit, ThroatExcludeFrac);
    LogCount = LogCount + 1;
    LogGeom(LogCount, :) = FullX;
    LogLife(LogCount) = PLife; LogDrop(LogCount) = PDrop; LogMeanTemp(LogCount) = PTemp;
    LogPeakTemp(LogCount) = PPeak; LogSecondaryTemp(LogCount) = PSecondary;
    LogCost(LogCount) = CostFromRaw(PLife, PDrop, PSecondary, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);
end

% Trim log to actual size
LogGeom = LogGeom(1:LogCount, :);
LogLife = LogLife(1:LogCount);
LogDrop = LogDrop(1:LogCount);
LogMeanTemp = LogMeanTemp(1:LogCount);
LogPeakTemp = LogPeakTemp(1:LogCount);
LogSecondaryTemp = LogSecondaryTemp(1:LogCount);
LogCost = LogCost(1:LogCount);

%% Final Selection
LinerODAll = LinerODCalc(LogGeom(:, 1:9), R_chamber);
ValidIdx = ~isnan(LogLife) & (LogDrop <= MaxDP) & (LinerODAll <= ODLimit);

if ~any(ValidIdx)
    error('No valid designs satisfying constraints were found.');
end

% Champion by combined objective (life + mean-temp + peak-temp tradeoff)
FeasibleCost = LogCost; FeasibleCost(~ValidIdx) = Inf;
[~, ChampIdx] = min(FeasibleCost);
Champion = LogGeom(ChampIdx, :);

% Diagnostic: pure max-life champion among valid designs (for comparison)
ValidLifeOnly = LogLife; ValidLifeOnly(~ValidIdx) = -Inf;
[~, LifeChampIdx] = max(ValidLifeOnly);
LifeChampion = LogGeom(LifeChampIdx, :);

fprintf('\n======================================================\n');
fprintf(' CHAMPION (combined life / secondary-cooling objective):\n');
fprintf('  Lifespan:           %.2f Cycles\n', LogLife(ChampIdx));
fprintf('  Pressure Drop:      %.2f psi\n', LogDrop(ChampIdx));
fprintf('  Throat Temp (peak, drives Life): %.2f K\n', LogPeakTemp(ChampIdx));
fprintf('  Secondary Temp (non-throat, optimized): %.2f K\n', LogSecondaryTemp(ChampIdx));
fprintf('  Whole-Profile Mean Temp (reference only): %.2f K\n', LogMeanTemp(ChampIdx));
fprintf(' GEOMETRY [inches]:\n');
fprintf('  Wall Thickness (C, T, N): [%.4f, %.4f, %.4f]\n', Champion(1:3));
fprintf('  Aspect Ratio (C, T, N):   [%.4f, %.4f, %.4f]\n', Champion(4:6));
fprintf('  Channel Width (C, N):     [%.4f, %.4f, %.4f]\n', Champion(7:9));
fprintf('  Channel Count:            %i\n', round(Champion(10)));
fprintf('  Liner OD:                 %.3f in\n', LinerODCalc(Champion(1:9), R_chamber));
fprintf('======================================================\n');

fprintf('\n[Diagnostic] Best pure-life design (ignoring the secondary cooling objective):\n');
fprintf('  Lifespan: %.2f Cycles | Drop: %.2f psi | ThroatT: %.2f K | SecondaryT: %.2f K | NC: %i\n', ...
    LogLife(LifeChampIdx), LogDrop(LifeChampIdx), LogPeakTemp(LifeChampIdx), LogSecondaryTemp(LifeChampIdx), round(LifeChampion(10)));

fprintf('\nTotal physics evaluations: %i\n', LogCount);

% Fetch the champion's full 100-station temperature profile (mean/peak
% only were kept in the log; this one extra call is just for the
% diagnostic plot below)
[~, ~, ChampionTempProfile] = SKRegen2_ElectricBoogalo(Data, round(Champion(10)), ...
    Champion(1:3), Champion(4:6), Champion(7:9), 1);

% Run SKR2 in its own verbose/plotting mode (last arg = 1) on the champion
% geometry to inspect whatever internal diagnostics it exposes. Commented
% out by default since it opens SKR2's own figures/console output.
% SKRegen2_ElectricBoogalo(Data, round(Champion(10)), Champion(1:3), Champion(4:6), Champion(7:9), 1);

%% Plotting
TotalEvals = LogCount;
CrashIdx = isnan(LogLife) | isnan(LogDrop);
ConstraintFailIdx = ~CrashIdx & (LogDrop > MaxDP | LinerODAll > ODLimit);
FailedIdx = CrashIdx | ConstraintFailIdx;

BestSoFarLife = nan(TotalEvals, 1);
CurrentMax = -Inf;
for i = 1:TotalEvals
    if ValidIdx(i) && LogLife(i) > CurrentMax
        CurrentMax = LogLife(i);
    end
    BestSoFarLife(i) = CurrentMax;
end

figure('Name', 'DE Optimization Convergence', 'WindowStyle', 'docked');

subplot(3, 1, 1); hold on; set(gca, 'FontName', 'Times New Roman');
scatter(find(FailedIdx), zeros(sum(FailedIdx),1), 15, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
scatter(find(ValidIdx), LogLife(ValidIdx), 15, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
plot(1:TotalEvals, BestSoFarLife, 'b-', 'LineWidth', 2);
xline(NumInit, 'k--', 'End LHS Init');
title('Lifespan across evaluations'); ylabel('Cycles'); grid on;

subplot(3, 1, 2); hold on; set(gca, 'FontName', 'Times New Roman');
scatter(find(FailedIdx), LogMeanTemp(FailedIdx), 15, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
scatter(find(ValidIdx), LogPeakTemp(ValidIdx), 15, [0.6 0.0 0.6], 'filled', 'MarkerFaceAlpha', 0.25);
scatter(find(ValidIdx), LogSecondaryTemp(ValidIdx), 15, [0.85 0.33 0.1], 'filled', 'MarkerFaceAlpha', 0.5);
legend('Failed', 'Throat (Peak) Temp - drives Life', 'Secondary (non-throat) Temp - optimized', 'Location', 'best');
title('Throat vs. non-throat wall temperature across evaluations'); ylabel('Temp [K]'); grid on;

subplot(3, 1, 3); hold on; set(gca, 'FontName', 'Times New Roman');
scatter(find(FailedIdx), LogDrop(FailedIdx), 15, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
scatter(find(ValidIdx), LogDrop(ValidIdx), 15, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
yline(MaxDP, 'r-', 'Max Constraint', 'LineWidth', 2);
title('Pressure drop across evaluations'); xlabel('Evaluation Number'); ylabel('Pressure Drop [psi]'); grid on;

%% Diagnostic Figure 2: DE convergence, life/temp tradeoff, channel count, champion profile
figure('Name', 'DE Diagnostics', 'WindowStyle', 'docked');

% (a) DE best-cost convergence history
subplot(2, 2, 1); hold on; set(gca, 'FontName', 'Times New Roman');
plot(1:GensRun, BestCostHistory(1:GensRun), 'b-', 'LineWidth', 2);
title('DE Convergence'); xlabel('Generation'); ylabel('Best Cost (lower = better)'); grid on;

% (b) Life vs. Secondary (non-throat) Temp tradeoff (Pareto-style view), colored by pressure drop
% Note: plotting against throat/peak temp here would be near-degenerate,
% since peak temp is what Life is already computed from. Secondary temp
% is the actual second objective being traded off against Life.
subplot(2, 2, 2); hold on; set(gca, 'FontName', 'Times New Roman');
scatter(LogSecondaryTemp(ValidIdx), LogLife(ValidIdx), 20, LogDrop(ValidIdx), 'filled', 'MarkerFaceAlpha', 0.6);
cb = colorbar; cb.Label.String = 'Pressure Drop [psi]';
plot(LogSecondaryTemp(ChampIdx), LogLife(ChampIdx), 'p', 'MarkerSize', 16, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k', 'DisplayName', 'Champion');
plot(LogSecondaryTemp(LifeChampIdx), LogLife(LifeChampIdx), 'd', 'MarkerSize', 12, ...
    'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', 'DisplayName', 'Max-Life Only');
legend('Location', 'best');
title('Life vs. Secondary (non-throat) Temp Tradeoff'); xlabel('Secondary Wall Temp [K]'); ylabel('Lifespan [Cycles]'); grid on;

% (c) Channel count sensitivity
subplot(2, 2, 3); hold on; set(gca, 'FontName', 'Times New Roman');
scatter(round(LogGeom(FailedIdx, 10)), zeros(sum(FailedIdx), 1), 15, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.2);
scatter(round(LogGeom(ValidIdx, 10)), LogLife(ValidIdx), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.5);
plot(round(Champion(10)), LogLife(ChampIdx), 'p', 'MarkerSize', 16, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'k');
title('Channel Count vs. Lifespan'); xlabel('Channel Count'); ylabel('Lifespan [Cycles]'); grid on;

% (d) Champion's full axial temperature profile
subplot(2, 2, 4); hold on; set(gca, 'FontName', 'Times New Roman');
plot(1:length(ChampionTempProfile), ChampionTempProfile, 'r-', 'LineWidth', 1.5);
yline(mean(ChampionTempProfile), 'b--', 'Mean', 'LineWidth', 1.5);
yline(max(ChampionTempProfile), 'k--', 'Max', 'LineWidth', 1);
title('Champion Axial Temperature Profile'); xlabel('Axial Station (1-100)'); ylabel('Wall Temp [K]'); grid on;

%% Diagnostic Figure 3: Parameter sensitivity grid (each design var vs. lifespan)
DimLabels = {'WT Chamber', 'WT Throat', 'WT Nozzle', ...
             'AR Chamber', 'AR Throat', 'AR Nozzle', ...
             'CW Chamber', 'CW Mid', 'CW Nozzle', 'Channel Count'};

figure('Name', 'Parameter Sensitivity vs. Lifespan', 'WindowStyle', 'docked');
for d = 1:NumDims
    subplot(2, 5, d); hold on; set(gca, 'FontName', 'Times New Roman');
    xVals = LogGeom(ValidIdx, d);
    if d == NumDims
        xVals = round(xVals);
    end
    scatter(xVals, LogLife(ValidIdx), 12, LogSecondaryTemp(ValidIdx), 'filled', 'MarkerFaceAlpha', 0.6);
    title(DimLabels{d}, 'FontSize', 9);
    if d == 1 || d == 6
        ylabel('Lifespan [Cycles]');
    end
    grid on;
end
sgtitle('Design Variable vs. Lifespan (color = secondary/non-throat wall temp)');

%% Local Functions

function [Life, Drop, MeanTemp, PeakTemp, SecondaryTemp] = EvaluatePhysicsRaw(x, Data, R_chamber, ODLimit, ThroatExcludeFrac)
% Wraps SKRegen2_ElectricBoogalo, enforcing the OD constraint cheaply
% before spending a physics call. Collapses the 100-station temp profile
% into three metrics:
%   MeanTemp      - whole-profile average, kept for reference/plots only,
%                   NOT used in the cost function (see CostFromRaw).
%   PeakTemp      - the single hottest station. This always lands at the
%                   throat and already drives Life via the peak-strain
%                   calc inside SKRegen2_ElectricBoogalo, so it is
%                   likewise NOT re-penalized in the cost function -
%                   doing so would just double-count what Life already
%                   optimizes.
%   SecondaryTemp - the mean temperature with a window around the throat
%                   (peak) station excluded on both sides. This is what
%                   actually goes into the cost: it targets "the rest of
%                   the temperatures" (chamber + nozzle cooling quality
%                   at minimum coolant flow) without being dominated by
%                   the throat spike Life is already responsible for.
    WT = x(1:3); AR = x(4:6); CW = x(7:9); NC = round(x(10));

    if LinerODCalc(x(1:9), R_chamber) > ODLimit
        Life = NaN; Drop = NaN; MeanTemp = NaN; PeakTemp = NaN; SecondaryTemp = NaN;
        return;
    end

    try
        [Life, Drop, TempProfile] = SKRegen2_ElectricBoogalo(Data, NC, WT, AR, CW, 0);
        if ~isreal(Life) || ~isreal(Drop) || any(~isreal(TempProfile)) || isnan(Life)
            Life = NaN; Drop = NaN; MeanTemp = NaN; PeakTemp = NaN; SecondaryTemp = NaN;
            return;
        end
        TempProfile = TempProfile(:);
        NStations = numel(TempProfile);
        [PeakTemp, PeakIdx] = max(TempProfile);
        MeanTemp = mean(TempProfile);

        HalfWidth = max(2, round(ThroatExcludeFrac * NStations));
        ExcludeMask = false(NStations, 1);
        ExcludeMask(max(1, PeakIdx - HalfWidth):min(NStations, PeakIdx + HalfWidth)) = true;

        RemainingStations = TempProfile(~ExcludeMask);
        if isempty(RemainingStations)
            % Exclusion window swallowed the whole profile (only possible
            % for a pathologically short profile) - fall back to the
            % overall mean rather than returning NaN.
            SecondaryTemp = MeanTemp;
        else
            SecondaryTemp = mean(RemainingStations);
        end
    catch
        Life = NaN; Drop = NaN; MeanTemp = NaN; PeakTemp = NaN; SecondaryTemp = NaN;
    end
end

function cost = CostFromRaw(Life, Drop, SecondaryTemp, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight)
% Combines life (maximize), secondary/non-throat temp (minimize), and a
% soft pressure-drop penalty into a single scalar cost to minimize.
%
% Peak (throat) temp is deliberately NOT in this cost: it always occurs
% at the throat and already dictates Life via the peak-strain station, so
% penalizing it again here would just double-count the same physics and
% pull the search away from the best-life designs. SecondaryTemp (the
% throat-excluded mean from EvaluatePhysicsRaw) is what represents the
% true secondary goal - good cooling everywhere else at minimum coolant
% flow - and is what actually gets minimized here.
    BigPenalty = 100;
    if isnan(Life)
        cost = BigPenalty;
        return;
    end
    LifeNorm = Life / LifeScale;
    SecondaryNorm = SecondaryTemp / SecondaryScale;
    if Drop > MaxDP
        PressPenalty = 5 + 20 * ((Drop - MaxDP) / MaxDP)^2;
    else
        PressPenalty = 0;
    end
    cost = -LifeNorm + SecondaryTempWeight * SecondaryNorm + PressPenalty;
end

function cost = CostWrapperFcn(x9, NC_Fixed, Data, R_chamber, ODLimit, ThroatExcludeFrac, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight)
% fmincon-facing scalar objective for the polish phase (NC held fixed).
    [Life, Drop, ~, ~, SecondaryTemp] = EvaluatePhysicsRaw([x9, NC_Fixed], Data, R_chamber, ODLimit, ThroatExcludeFrac);
    cost = CostFromRaw(Life, Drop, SecondaryTemp, MaxDP, LifeScale, SecondaryScale, SecondaryTempWeight);
end

function OD = LinerODCalc(X, R_chamber)
    OD = 2 * (R_chamber + X(:, 1) + X(:, 7) .* X(:, 4));
end

function [Archives, Improved] = UpdateArchives(Geom, Cost, Life, Drop, Temp, Peak, Secondary, NumDims, Archives)
% Merged & fixed from BoundaryExplorer.m. Evaluates the design against
% all NumDims 1D archives (one bin per dimension slice) and records it
% wherever it beats the current best (lowest) cost for that slice.
%
% NOTE: the original BoundaryExplorer version took Archives as an input
% and mutated it locally without returning it - since MATLAB structs/
% cells are pass-by-value, none of those updates ever reached the
% caller's copy, so its archive plots would have stayed empty. Archives
% is now an explicit output so updates persist.
    Improved = false;

    for d = 1:NumDims
        binIdx = find(Geom(d) >= Archives{d}.Edges(1:end-1) & Geom(d) < Archives{d}.Edges(2:end), 1, 'first');
        if isempty(binIdx); binIdx = length(Archives{d}.MinCost); end % catch boundary edge cases

        if Cost < Archives{d}.MinCost(binIdx)
            Archives{d}.MinCost(binIdx) = Cost;
            Archives{d}.Geom(binIdx, :) = Geom;
            Archives{d}.Life(binIdx) = Life;
            Archives{d}.Drop(binIdx) = Drop;
            Archives{d}.Temp(binIdx) = Temp;
            Archives{d}.Peak(binIdx) = Peak;
            Archives{d}.Secondary(binIdx) = Secondary;
            Improved = true;
        end
    end
end