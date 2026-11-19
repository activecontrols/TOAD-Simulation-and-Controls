%% Special script meant to test SKIPPER Sensitivity to Heatflux Fit to TADPOLE Data
%
% Requires passing an additional parameter to SKRegen2.m with the heatflux
% factor, and rerunning CEA previous to running this script. Using iPDR
% geometry as of 8/3/26.
%
% Pablo Plata

%% Preload Material data for passage into regen script
% Only load the material and calls actually being used rn
addpath("cea\");
addpath("IPA Data\");
addpath("Material Data\");
addpath("Contours\");
Data = LoadData();

%% Optimal Geometry
NC = 65;
WT = [0.05, 0.05, 0.05];
AR = [1.62 3.4 4.5];
CW = [0.05, 0.020, 0.02];
CH = [AR(1)*CW(1), AR(2:3) * CW(2)];

% Loop
Bounds = [0.3, 3];
NumIter = 100;
CoefVector = linspace(Bounds(1), Bounds(2), NumIter);
Lifespan = zeros(NumIter, 1);
PressDrop = zeros(NumIter, 1);
TempMax = zeros(NumIter, 1);

for i = 1:NumIter
    FitCoef = CoefVector(i);
    try 
        [Lifespan(i), PressDrop(i), TempArray] = SKRegen2_ElectricBoogalo(Data, NC, WT, AR, CW, 0, FitCoef);
        TempMax(i) = max(TempArray);
    catch ME
        Lifespan(i) = 0;
        PressDrop(i) = 0;
        TempMax(i) = 0;
    end
end

%% Plots & format handling 
% Create a figure for plotting
figure;
ValidIdx = Lifespan > 0;
dropOffIdx = find(ValidIdx, 1, 'last');

% Plot Lifespan vs Fit Coefficient
subplot(2, 1, 1);
hold on; grid on; % Moved 
plot(CoefVector, Lifespan, '-r', 'LineWidth', 2);
xlim([CoefVector(1), CoefVector(end)]);
xline(CoefVector(dropOffIdx), '--');
xlabel('Fit Coefficient');
ylabel('Lifespan');
title('Lifespan vs Fit Coefficient');
legend('Valid Solutions', 'Drop-off Point');

x_val = CoefVector(dropOffIdx);
y_val1 = Lifespan(dropOffIdx);
text_str1 = sprintf('  Drop-off (Fit Coef: %.2f, Lifespan: %.2f)', x_val, y_val1);

text(x_val, y_val1, text_str1, ...
    'HorizontalAlignment', 'left', ...   
    'VerticalAlignment', 'bottom', ...   
    'FontSize', 10, ...
    'FontWeight', 'bold');


% Plot Pressure Drop vs Fit Coefficient
subplot(2, 1, 2);
hold on; grid on;
plot(CoefVector, TempMax, '-r', 'LineWidth', 2);
xlim([CoefVector(1), CoefVector(end)]);
xline(CoefVector(dropOffIdx), '--');
xlabel('Fit Coefficient');
ylabel('Max Temperature [K]');
title('Max Temperature vs Fit Coefficient');
legend('Valid Solutions', 'Drop-off Point');

x_val = CoefVector(dropOffIdx);
y_val1 = TempMax(dropOffIdx);
text_str1 = sprintf('  Drop-off (Fit Coef: %.2f, Lifespan: %.2f)', x_val, y_val1);

text(x_val, y_val1, text_str1, ...
    'HorizontalAlignment', 'left', ...   
    'VerticalAlignment', 'bottom', ...   
    'FontSize', 10, ...
    'FontWeight', 'bold');


