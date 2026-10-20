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
cutoffFreq = 10; % Hz

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
        % plot(Data(i).Encoder.Pos(:, 1) * 360, Params(i).CvOX, 'bx'); hold on; grid on;
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
        % plot(DataF(i).Sensor.PT_TCA_Chamber, Params(i).Cstar, 'gx'); hold on; grid on;
    end

    % Extract CdA OX
    if Data(i).TestType ~= "IPA CF"
        gc = 32.174 * 12;
        mdot = DataF(i).Sensor.PT_OX_Mdot;
        P_c = DataF(i).Sensor.PT_TCA_Chamber;
        P_inj = DataF(i).Sensor.PT_OX_ValveDown;
        T_inj = DataF(i).Sensor.TC_OX_Valve;
        OX_Props = FastFluidProps("Oxygen", 'From_P_T', P_valve * PSI_TO_PA, T_inj); 
        rho = OX_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_OX = abs(mdot./((P_inj - P_c) .* 2 * gc .*rho).^0.5);
        Params(i).CdA_OX(rho < Rho_MIN) = NaN;
        % plot(DataF(i).Time, Params(i).CdA_OX, 'b'); hold on; grid on;
    end

    % Extract CdA FU
    if Data(i).TestType ~= "LOX CF"
        mdot = DataF(i).Sensor.PT_IPA_Mdot;
        P_c = DataF(i).Sensor.PT_TCA_Chamber;
        P_inj = DataF(i).Sensor.PT_IPA_ValveDown;
        IPA_Props = FastFluidProps("IPA", 'From_P_T', P_valve * PSI_TO_PA, 293.17); 
        rho = IPA_Props.rho * KGM3_TO_LBIN3;
        Params(i).CdA_IPA = abs(mdot./((P_inj - P_c) .*2 * gc .*rho).^0.5);
        % plot(DataF(i).Time, Params(i).CdA_IPA, 'r'); hold on; grid on;
    end
    fprintf('Test %i / %i done!\n', i, length(fileList));
end

%% Plotting
Trigger = [3; 2]; % [Start offset; End offset] 

% OX Cv vs. Valve Angle
figure('Name', 'Cv vs Angle');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        subplot(2,1,1);
        plot(Data(i).Encoder.Pos(idx, 1) * 360, Params(i).CvOX(idx), 'bx'); hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        subplot(2,1,2);
        plot(Data(i).Encoder.Pos(idx, 2) * 360, Params(i).CvIPA(idx), 'rx'); hold on; grid on;
    end
end
subplot(2,1,1);
xlabel('LOX Valve Angle (deg)');
ylabel('Cv_{OX}');
ylim([0, 3.2]);
xlim([20, 90]);
title('LOX Valve Cv vs. Position');
subplot(2,1,2);
xlabel('IPA Valve Angle (deg)');
ylabel('Cv_{IPA}');
ylim([0, 3.2]);
xlim([20, 90]);
title('IPA Valve Cv vs. Position');
sgtitle('Flow Coefficient vs. Valve Angle')

% C* vs. Chamber Pressure
figure('Name', 'Cstar vs Chamber Pressure');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if DataF(i).TestType == "HF"
        plot(DataF(i).Sensor.PT_TCA_Chamber(idx), Params(i).Cstar(idx), 'rx'); hold on; grid on;
    end
end
xlabel('Chamber Pressure (psi)');
xlim([80, 300]);
ylabel('C^* (ft/s)');
title('Hotfire C^* vs. Chamber Pressure');
% Injector CdA vs. Time
figure('Name', 'Injector CdA vs Time');
for i = 1:length(fileList)
    idx = (Data(i).Time > Trigger(1)) & (Data(i).Time < (Data(i).Time(end) - Trigger(2)));
    if Data(i).TestType ~= "IPA CF"
        plot(DataF(i).Time(idx), Params(i).CdA_OX(idx), 'b'); hold on; grid on;
    end
    if Data(i).TestType ~= "LOX CF"
        plot(DataF(i).Time(idx), Params(i).CdA_IPA(idx), 'r'); hold on; grid on;
    end
end
ylim([0, 0.05]);
xlabel('Time (s)');
ylabel('CdA (in^2)');
title('Injector CdA vs. Time (Blue = OX, Red = IPA)');