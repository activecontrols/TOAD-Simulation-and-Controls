%% Data Table Extraction Script for SKIPPER Testing
% Extract Data for ALL Tests
addpath('DATA');
fileList = dir('DATA\*.xlsx');
Data = struct();
KGM3_TO_LBIN3 = 3.61273e-5;
PSI_TO_PA = 6894.7;
Rho_MIN = FastFluidProps("Oxygen", 'From_P_T', 550, 130).rho;

for i = 1:length(fileList)
    currentFile = fullfile(fileList(i).folder, fileList(i).name);
    flightData = readtable(currentFile);
    
    % Data struct w/ data for all tests
    % Test type 
    if contains(fileList(i).name, "HOTFIRE")
        Data(i).TestType = "HF";
    elseif contains(fileList(i).name, "CV") || contains(fileList(i).name, "CL")
        if contains(fileList(i).name, "LOX")
            Data(i).TestType = "LOX CF";
        else
            Data(i).TestType = "IPA CF";
        end
    else
        Data(i).TestType = "Undetermined";
    end

    Data(i).Time = flightData.time;
    Data(i).Phase = flightData.phase;

    % PT Data
    Data(i).Sensor.PT_TCA_Chamber = flightData.chamber_pressure;
    Data(i).Sensor.PT_OX_ValveUp = flightData.lox_valve_upstream_pressure;
    Data(i).Sensor.PT_OX_ValveDown = flightData.lox_valve_downstream_pressure;
    Data(i).Sensor.PT_OX_ValveDownCalc = flightData.lox_valve_downstream_pressure_calc;
    Data(i).Sensor.PT_OX_Venturi = flightData.lox_venturi_differential_pressure;
    Data(i).Sensor.PT_OX_Mdot = flightData.lox_mdot;
    Data(i).Sensor.PT_IPA_ValveUp = flightData.ipa_valve_upstream_pressure;
    Data(i).Sensor.PT_IPA_ValveDown = flightData.ipa_valve_downstream_pressure;
    Data(i).Sensor.PT_IPA_ValveDownCalc = flightData.ipa_valve_downstream_pressure_calc;
    Data(i).Sensor.PT_IPA_Venturi = flightData.ipa_venturi_differential_pressure;
    Data(i).Sensor.PT_IPA_Mdot = flightData.ipa_mdot;

    % TC Data
    Data(i).Sensor.TC_OX_Venturi = flightData.lox_venturi_temperature;
    Data(i).Sensor.TC_OX_Valve = flightData.lox_valve_temperature;

    % Valve Command / Control Data
    Data(i).Control.Cmd.Thrust = flightData.thrust_cmd;
    Data(i).Control.Cmd.Pos = [flightData.lox_pos_cmd, flightData.ipa_pos_cmd];
    Data(i).Control.PI.Chamber = [flightData.chamber_pressure_controller_p_component, flightData.chamber_pressure_controller_i_component];
    Data(i).Control.PI.LOX = [flightData.lox_angle_controller_p_component, flightData.lox_angle_controller_i_component];
    Data(i).Control.PI.IPA = [flightData.ipa_angle_controller_p_component, flightData.ipa_angle_controller_i_component];
    Data(i).Control.OL.Mdot = [flightData.ol_lox_mdot, flightData.ol_ipa_mdot];
    Data(i).Control.OL.Angle = [flightData.ol_lox_angle, flightData.ol_ipa_angle];

    % Valve Encoder Data
    Data(i).Encoder.Pos = [flightData.lox_pos, flightData.ipa_pos];
    Data(i).Encoder.Vel = [flightData.lox_vel, flightData.ipa_vel];
end

%% Symmetric Filtering 
DataF = Data;
cutoffFreq = 3; % Hz

for i = 1:length(Data)
    dt = mean(diff(Data(i).Time));
    fs = 1 / dt;
    [b, a] = butter(2, cutoffFreq / (fs / 2), 'low');

    DataF(i).Time = Data(i).Time;
    DataF(i).Phase = Data(i).Phase;

    % PT Data
    DataF(i).Sensor.PT_TCA_Chamber = filtfilt(b, a, Data(i).Sensor.PT_TCA_Chamber);
    DataF(i).Sensor.PT_OX_ValveUp = filtfilt(b, a, Data(i).Sensor.PT_OX_ValveUp);
    DataF(i).Sensor.PT_OX_ValveDown = filtfilt(b, a, Data(i).Sensor.PT_OX_ValveDown);
    DataF(i).Sensor.PT_OX_ValveDownCalc = filtfilt(b, a, Data(i).Sensor.PT_OX_ValveDownCalc);
    DataF(i).Sensor.PT_OX_Venturi = filtfilt(b, a, Data(i).Sensor.PT_OX_Venturi);
    DataF(i).Sensor.PT_OX_Mdot = filtfilt(b, a, Data(i).Sensor.PT_OX_Mdot);
    DataF(i).Sensor.PT_IPA_ValveUp = filtfilt(b, a, Data(i).Sensor.PT_IPA_ValveUp);
    DataF(i).Sensor.PT_IPA_ValveDown = filtfilt(b, a, Data(i).Sensor.PT_IPA_ValveDown);
    DataF(i).Sensor.PT_IPA_ValveDownCalc = filtfilt(b, a, Data(i).Sensor.PT_IPA_ValveDownCalc);
    DataF(i).Sensor.PT_IPA_Venturi = filtfilt(b, a, Data(i).Sensor.PT_IPA_Venturi);
    DataF(i).Sensor.PT_IPA_Mdot = filtfilt(b, a, Data(i).Sensor.PT_IPA_Mdot);

    % TC Data
    DataF(i).Sensor.TC_OX_Venturi = filtfilt(b, a, Data(i).Sensor.TC_OX_Venturi);
    DataF(i).Sensor.TC_OX_Valve = filtfilt(b, a, Data(i).Sensor.TC_OX_Valve);

    % Valve Command / Control Data
    DataF(i).Control.Cmd.Thrust = filtfilt(b, a, Data(i).Control.Cmd.Thrust);
    DataF(i).Control.Cmd.Pos = filtfilt(b, a, Data(i).Control.Cmd.Pos);
    DataF(i).Control.PI.Chamber = filtfilt(b, a, Data(i).Control.PI.Chamber);
    DataF(i).Control.PI.LOX = filtfilt(b, a, Data(i).Control.PI.LOX);
    DataF(i).Control.PI.IPA = filtfilt(b, a, Data(i).Control.PI.IPA);
    DataF(i).Control.OL.Mdot = filtfilt(b, a, Data(i).Control.OL.Mdot);
    DataF(i).Control.OL.Angle = filtfilt(b, a, Data(i).Control.OL.Angle);

    % Valve Encoder Data
    DataF(i).Encoder.Pos = filtfilt(b, a, Data(i).Encoder.Pos);
    DataF(i).Encoder.Vel = filtfilt(b, a, Data(i).Encoder.Vel);
end


%% Analysis 
% Extract measured Cv values for each valve 
Params = struct();
for i = 1:length(fileList)

    % Extract OX Cv
    if Data(i).TestType ~= "IPA CF"
        mdot = DataF(i).Sensor.PT_OX_Mdot;
        P_tank = DataF(i).Sensor.PT_OX_ValveUp;
        P_valve = DataF(i).Sensor.PT_OX_ValveDown;
        T_inj = DataF(i).Sensor.TC_OX_Venturi;
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_tank * PSI_TO_PA, T_inj); 
        rho = OX_Props.rho * KGM3_TO_LBIN3;
        rho_water = 0.0361; 
        Params(i).CvOX = abs(mdot ./(rho .* rho_water .* (P_tank - P_valve)).^0.5 * 60 / 231);
    end

    % Extract FU Cv
    if Data(i).TestType ~= "LOX CF"
        mdot = DataF(i).Sensor.PT_IPA_Mdot;
        P_tank = DataF(i).Sensor.PT_IPA_ValveUp;
        P_valve = DataF(i).Sensor.PT_IPA_ValveDown;
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_valve * PSI_TO_PA, 293.17); 
        rho = IPA_Props.rho * KGM3_TO_LBIN3;
        rho_water = 0.0361; 
        Params(i).CvIPA = abs(mdot ./(rho .* rho_water .* (P_tank - P_valve)).^0.5 * 60 / 231);
    end

    % Extract Cstar for hotfires
    if DataF(i).TestType == "HF"
        g0 = 32.17;
        mdot_tot = DataF(i).Sensor.PT_IPA_Mdot + DataF(i).Sensor.PT_OX_Mdot;
        Pc = DataF(i).Sensor.PT_TCA_Chamber;
        At = 1.679;
        Params(i).Cstar = Pc .* At ./ mdot_tot .* g0;

        % Cstar Eff
        Params(i).OF    = DataF(i).Sensor.PT_OX_Mdot ./ DataF(i).Sensor.PT_IPA_Mdot;
        Params(i).Cstar = Pc .* At ./ mdot_tot .* g0;

        % Interpolate theoretical C* and get C* efficiency
        [Params(i).EtaCstar, Params(i).Cstar_theo] = getCstarEff(Pc, Params(i).OF, Params(i).Cstar);
    end

    % Extract CdA OX
    if Data(i).TestType ~= "IPA CF"
        gc = 32.174 * 12;
        mdot = DataF(i).Sensor.PT_OX_Mdot;
        P_c = DataF(i).Sensor.PT_TCA_Chamber;
        P_inj = DataF(i).Sensor.PT_OX_ValveDown;
        T_inj = DataF(i).Sensor.TC_OX_Valve;
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_inj * PSI_TO_PA, T_inj); 
        rho = OX_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_OX = abs(mdot./((P_inj - P_c) .* 2 * gc .*rho).^0.5);
        Params(i).CdA_OX(rho < Rho_MIN * KGM3_TO_LBIN3) = NaN;
    end

    % Extract CdA FU
    if Data(i).TestType ~= "LOX CF"
        gc = 32.174 * 12;
        mdot = DataF(i).Sensor.PT_IPA_Mdot;
        P_c = DataF(i).Sensor.PT_TCA_Chamber;
        P_inj = DataF(i).Sensor.PT_IPA_ValveDown;
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_inj * PSI_TO_PA, 293.17); 
        rho = IPA_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_IPA = abs(mdot./((P_inj - P_c) .*2 * gc .*rho).^0.5);
    end

    % Extract Momentum Ratio for hotfires
    if DataF(i).TestType == "HF"
        mdot_ox = DataF(i).Sensor.PT_OX_Mdot;
        mdot_ipa = DataF(i).Sensor.PT_IPA_Mdot;

        % Compute densities at the injector manifold
        P_inj_ox = DataF(i).Sensor.PT_OX_ValveDown;
        T_inj_ox = DataF(i).Sensor.TC_OX_Valve;
        rho_ox = FastFluidProps("Oxygen", 'From_P_T', P_inj_ox * PSI_TO_PA, T_inj_ox).rho * KGM3_TO_LBIN3;
        rho_ox(rho_ox < Rho_MIN * KGM3_TO_LBIN3) = NaN;
        P_inj_ipa = DataF(i).Sensor.PT_IPA_ValveDown;
        rho_ipa = FastFluidProps("IPA", 'From_P_T', P_inj_ipa * PSI_TO_PA, 293.17).rho * KGM3_TO_LBIN3;

        % Flow vels
        Params(i).Vel_OX  = mdot_ox ./ (rho_ox .* Params(i).CdA_OX * 12);
        Params(i).Vel_IPA = mdot_ipa ./ (rho_ipa .* Params(i).CdA_IPA * 12);

        % Using effective flow areas from above
        Mom_OX = (mdot_ox.^2) ./ (rho_ox .* Params(i).CdA_OX);
        Mom_IPA = (mdot_ipa.^2) ./ (rho_ipa .* Params(i).CdA_IPA);

        % Momentum Ratio (OX / Fuel)
        Params(i).MomRatio = Mom_OX ./ Mom_IPA;
    end

    fprintf('Test %i / %i done!\n', i, length(fileList));
end

%% Plotting
Trigger = [2.7; 1]; % [Start offset; End offset]
TimeExclude = [8.5, 10.5; 
               15.5, 17.5];

% Generate consistent color palette and clean test labels for legends
colors = lines(length(fileList));
testNames = strings(length(fileList), 1);
for i = 1:length(fileList)
    [~, name, ~] = fileparts(fileList(i).name);
    testNames(i) = strrep(name, '_', ' ');
end

% OX & IPA Cv vs. Valve Angle
figure('Name', 'Cv vs Angle');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        subplot(2,1,1);
        plot(Data(i).Encoder.Pos(idx, 1) * 360, Params(i).CvOX(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        subplot(2,1,2);
        plot(Data(i).Encoder.Pos(idx, 2) * 360, Params(i).CvIPA(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
subplot(2,1,1);
xlabel('LOX Valve Angle (deg)');
ylabel('Cv_{OX}');
ylim([0, 3.2]);
xlim([20, 90]);
title('LOX Valve Cv vs. Position');
legend('show', 'Location', 'best');

subplot(2,1,2);
xlabel('IPA Valve Angle (deg)');
ylabel('Cv_{IPA}');
ylim([0, 3.2]);
xlim([20, 90]);
title('IPA Valve Cv vs. Position');
legend('show', 'Location', 'best');
sgtitle('Flow Coefficient vs. Valve Angle');

% C* vs. Chamber Pressure
figure('Name', 'Cstar vs Chamber Pressure');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if DataF(i).TestType == "HF"
        plot(DataF(i).Sensor.PT_TCA_Chamber(idx), Params(i).Cstar(idx), 'x', ...
            'Color', colors(i,:), 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Chamber Pressure (psi)');
xlim([80, 300]);
ylabel('C^* (ft/s)');
title('Hotfire C^* vs. Chamber Pressure');
legend('show', 'Location', 'best');

% Injector CdA vs. Time (Solid = OX, Dashed = IPA, Color = Test)
figure('Name', 'Injector CdA vs Time');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        plot(DataF(i).Time(idx), Params(i).CdA_OX(idx), '-', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i) + " (OX)"); 
        hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        plot(DataF(i).Time(idx), Params(i).CdA_IPA(idx), '--', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i) + " (IPA)"); 
        hold on; grid on;
    end
end
ylim([0, 0.05]);
xlabel('Time (s)');
ylabel('CdA (in^2)');
title('Injector CdA vs. Time (Solid = OX, Dashed = IPA)');
legend('show', 'Location', 'best');

% Momentum Ratio & C* vs. Time for hotfires
figure('Name', 'Momentum Ratio vs Time');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if DataF(i).TestType == "HF"
        subplot(1,2,1);
        plot(DataF(i).Time(idx), Params(i).MomRatio(idx), ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
        
        subplot(1,2,2);
        plot(DataF(i).Time(idx), Params(i).EtaCstar(idx), ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
subplot(1,2,1);
xlabel('Time (s)');
ylabel('Momentum Ratio (J_{OX} / J_{IPA})');
title('Hotfire Momentum Ratio vs. Time');
ylim([0, 8]);
legend('show', 'Location', 'best');

subplot(1,2,2);
xlabel('Time (s)');
ylabel('C^* (ft/s)');
title('Hotfire \eta C^* vs. Time');
legend('show', 'Location', 'best');

% TMR vs. C* (Steady Time Gate)
figure('Name', 'TMR vs. Cstar');
for i = 1:length(fileList)
    t = Data(i).Time;
    inExcludedWindow = any((t >= TimeExclude(:,1)') & (t <= TimeExclude(:,2)'), 2);
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2))) & ~inExcludedWindow;
      
    if DataF(i).TestType == "HF"
        plot(Params(i).MomRatio(idx), Params(i).EtaCstar(idx), 'x', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Momentum Ratio (J_{OX} / J_{IPA})');
ylabel('C^* (ft/s)');
title('Hotfire TMR vs. \eta C^* (Time-Gated)');
legend('show', 'Location', 'best');

% Combined Figure: Left = TMR vs. C*, Right = Split Velocity Scatter Colormaps
figure('Name', 'TMR and Injector Velocities vs. Cstar');

% --- LEFT COLUMN: TMR vs. C* (Spans rows 1 & 3 of a 2x2 grid) ---
subplot(2, 2, [1, 3]);
all_TMR = [];
all_Vel_IPA = [];
all_Vel_OX  = [];
all_CstarEff   = [];

for i = 1:length(fileList)
    t = Data(i).Time;
    inExcludedWindow = any((t >= TimeExclude(:,1)') & (t <= TimeExclude(:,2)'), 2);
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2))) & ~inExcludedWindow;
      
    if DataF(i).TestType == "HF"
        plot(Params(i).MomRatio(idx), Params(i).EtaCstar(idx), 'x', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
        
        % Pool data for the RHS scatter plots
        all_TMR = [all_TMR; Params(i).MomRatio(idx)];
        all_Vel_IPA = [all_Vel_IPA; Params(i).Vel_IPA(idx)];
        all_Vel_OX  = [all_Vel_OX;  Params(i).Vel_OX(idx)];
        all_CstarEff   = [all_CstarEff;   Params(i).EtaCstar(idx)];
    end
end
[p, S] = polyfit(all_TMR, all_CstarEff, 1);
MomRatioSweep = linspace(1, 3, 50);
plot(MomRatioSweep, polyval(p, MomRatioSweep), 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Curve Fit');
fprintf('TMR vs. CstarEff R^2: %.3f\n', S.rsquared);
xlabel('Momentum Ratio (J_{OX} / J_{IPA})');
ylabel('\eta C^*');
title('Hotfire TMR vs. C^* (Time-Gated)');
legend('show', 'Location', 'southeast');

% Clean NaNs and lock shared colorbar limits across both RHS subplots
valid = ~isnan(all_Vel_IPA) & ~isnan(all_Vel_OX) & ~isnan(all_CstarEff);
all_Vel_IPA = all_Vel_IPA(valid);
all_Vel_OX  = all_Vel_OX(valid);
all_CstarEff   = all_CstarEff(valid);

c_bounds = [min(all_CstarEff), max(all_CstarEff)];
vel_threshold = 110; % Split threshold (ft/s) between low and high regimes

% High vel
ax_high = subplot(2, 2, 2);
idx_high = all_Vel_IPA >= vel_threshold;
scatter(all_Vel_IPA(idx_high), all_Vel_OX(idx_high), 25, all_CstarEff(idx_high), 'filled');
grid on;
colormap(ax_high, turbo);
clim(ax_high, c_bounds);
cb1 = colorbar;
cb1.Label.String = '\eta C^* (ft/s)';
xlabel('Fuel (IPA) Injection Velocity (ft/s)');
ylabel('LOX Injection Velocity (ft/s)');
title('High Injection Velocity C^* Map');

% Low vel
ax_low = subplot(2, 2, 4);
idx_low = all_Vel_IPA < vel_threshold;
scatter(all_Vel_IPA(idx_low), all_Vel_OX(idx_low), 25, all_CstarEff(idx_low), 'filled');
grid on;
colormap(ax_low, turbo);
clim(ax_low, c_bounds);
cb2 = colorbar;
cb2.Label.String = '\eta C^*';
xlabel('Fuel (IPA) Injection Velocity (ft/s)');
ylabel('LOX Injection Velocity (ft/s)');
title('Low Injection Velocity C^* Map');

% Time vs. OF (Steady Time Gate)
figure('Name', 'OF Ratio Traces');
for i = 1:length(fileList)
    t = Data(i).Time;
    idx = (t > Trigger(1)) & (t < (t(end) - Trigger(2)));

    if DataF(i).TestType == "HF"
        plot(Data(i).Time(idx), Params(i).OF(idx), 'x', ...
            'Color', colors(i,:), 'LineWidth', 1.2, 'DisplayName', testNames(i)); 
        hold on; grid on;
    end
end
xlabel('Time [S]');
ylabel('OF Ratio');
title('OF Ratio vs Time Traces');
legend('show', 'Location', 'best');